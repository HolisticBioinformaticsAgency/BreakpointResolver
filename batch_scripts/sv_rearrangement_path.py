"""
sv_rearrangement_path.py
────────────────────────
Structural-variant rearrangement-path pipeline.

STAGES (sv_rearrangement_path mode)
1.  find_major_clip_positions        – genome-wide soft-clip / breakpoint discovery
2.  extract_split_reads              – fetch SA-tagged reads near a candidate pair
3.  pair_clip_positions              – pair left/right breakpoints by bridging-read count
4.  cluster_softclip_pairs           – merge pairs sharing >50% SA reads (Union-Find)
5.  compute_reference_window         – derive ref window from all alignments (≤1 Mb each)
6.  get_background_homozygous_sites  – pileup VAF>=0.95 sites across ALL reads in window
7.  find_split_read_specific_variants– per-position base counts in split vs. all reads;
                                       returns PDV positions after background masking
8.  collect_breakpoint_segments      – extend soft-clip tail inward to last PDV per read
9.  align_clips_to_ref               – minimap2 subprocess alignment of each segment
10. call_pdvs                        – compare segment vs. ref at each alignment position
11. cluster_by_landing_site          – group reads whose segments land at the same locus
12. build_path_from_clips            – one graph node per unique landing locus
13. plot_rearrangement_path          – box-and-edge figure with PDV annotations on edges
"""

import re
import os
import sys
import subprocess
import tempfile
import pysam
pysam.set_verbosity(0)
import pandas as pd
import numpy as np
from collections import Counter, defaultdict, namedtuple
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from IPython.display import display, Image, HTML

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
CLIP_OPS = {4, 5}
REF_OPS  = {0, 2, 3, 7, 8}   # M / D / N / = / X

PathSegment = namedtuple(
    "PathSegment",
    [
        "query_start",   # 0-based start in segment sequence
        "query_end",     # exclusive
        "ref_chrom",
        "ref_start",     # 0-based genomic start
        "ref_end",       # exclusive
        "strand",        # '+' or '-'
        "identity",      # fraction matching bases (0–1)
        "read_support",  # reads supporting this landing locus
        "segment_type",  # 'reference' | 'novel'
        "label",         # human-readable tag
        "pdvs",          # list of PDVAnnotation named-tuples (may be empty)
    ],
)

PDVAnnotation = namedtuple(
    "PDVAnnotation",
    ["ref_pos", "clip_base", "ref_base", "split_vaf", "all_vaf", "status"],
    # status: 'source' | 'dest' | 'novel'
)

ClipAlignment = namedtuple(
    "ClipAlignment",
    ["read_name", "segment_seq", "ref_chrom", "ref_start", "ref_end",
     "strand", "n_mismatches", "mapq", "pdvs"],
)

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ATCGatcgNn", "TAGCtagcNn"))[::-1]

def get_ref_chroms(ref_path: str):
    with pysam.FastaFile(ref_path) as fa:
        return set(fa.references)

def _has_adjacent_n(fasta, chrom: str, pos: int, window: int = 5) -> bool:
    if chrom not in fasta.references:
        return False
    return "N" in fasta.fetch(chrom, max(0, pos - window), pos + window).upper()

def _ref_span_from_cigar(cigar_tuples) -> int:
    return sum(l for op, l in cigar_tuples if op in REF_OPS)

def _parse_sa_tag(tag: str):
    """Parse SA tag string into list of (chrom, pos1, strand, cigar, mapq)."""
    result = []
    for entry in tag.rstrip(";").split(";"):
        f = entry.split(",")
        if len(f) < 5:
            continue
        try:
            result.append((f[0], int(f[1]), f[2], f[3], int(f[4])))
        except (ValueError, IndexError):
            continue
    return result

def _sa_ref_span(sa_cigar: str) -> int:
    ops = re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
    REF_STR = {"M", "D", "N", "=", "X"}
    return sum(int(l) for l, op in ops if op in REF_STR)

# ─────────────────────────────────────────────────────────────────────────────
# 1. find_major_clip_positions (unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def find_major_clip_positions(
    bam_path,
    chrom=None,
    min_clip_len=10,
    min_support=5,
    min_support_frac=0.05,
    cluster_dist=10,
    top_n=50,
    require_sa=False,
    min_mapq=60,
    save_dir=None,
):
    REF_OPS_STR = {"M", "D", "N", "=", "X"}
    CLIP_STR = {"S", "H"}

    def _region_depth(chrom, pos, window=50):
        count = 0
        with pysam.AlignmentFile(bam_path, "rb") as sf:
            for read in sf.fetch(chrom, max(0, pos - window), pos + window):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                count += 1
        return count

    left_pool  = defaultdict(list)
    right_pool = defaultdict(list)

    samfile    = pysam.AlignmentFile(bam_path, "rb")
    fetch_iter = samfile.fetch(chrom) if chrom else samfile.fetch()

    def _sa_majority_mapq(read, min_mapq):
        if not read.has_tag("SA"):
            return True
        entries = [e for e in read.get_tag("SA").rstrip(";").split(";") if e]
        if not entries:
            return True
        return any(
            len(sa.split(",")) >= 6 and int(sa.split(",")[4]) >= min_mapq
            for sa in entries
        )

    seen = set()
    for read in fetch_iter:
        if read.is_unmapped or not read.cigartuples:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.is_supplementary or read.is_secondary:
            continue
        if not _sa_majority_mapq(read, min_mapq):
            continue
        if require_sa and not read.has_tag("SA"):
            continue

        cigars   = read.cigartuples
        ref_name = read.reference_name
        first_op, first_len = cigars[0]
        last_op,  last_len  = cigars[-1]

        key_left = (read.query_name, read.reference_start)
        if first_op == 4 and first_len >= min_clip_len and key_left not in seen:
            seen.add(key_left)
            left_pool[ref_name].append((read.reference_start, first_len))

        key_right = (read.query_name, read.reference_end)
        if last_op == 4 and last_len >= min_clip_len and key_right not in seen:
            seen.add(key_right)
            right_pool[ref_name].append((read.reference_end, last_len))

        if read.has_tag("SA"):
            for sa_entry in read.get_tag("SA").rstrip(";").split(";"):
                fields = sa_entry.split(",")
                if len(fields) < 6:
                    continue
                sa_chr, sa_pos, sa_cigar, sa_mapq = (
                    fields[0], int(fields[1]), fields[3], int(fields[4])
                )
                if sa_mapq < min_mapq or sa_chr != ref_name:
                    continue
                ops = re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
                if not ops:
                    continue
                ref_len = sum(int(l) for l, op in ops if op in REF_OPS_STR)
                sa_end  = sa_pos + ref_len

                sa_left_len  = int(ops[0][0])  if ops[0][1]  in CLIP_STR else 0
                sa_right_len = int(ops[-1][0]) if ops[-1][1] in CLIP_STR else 0

                key_sa_left = (read.query_name, sa_pos, "left")
                if ops[0][1] in CLIP_STR and sa_left_len >= min_clip_len and key_sa_left not in seen:
                    seen.add(key_sa_left)
                    left_pool[ref_name].append((sa_pos, sa_left_len))

                key_sa_right = (read.query_name, sa_end, "right")
                if ops[-1][1] in CLIP_STR and sa_right_len >= min_clip_len and key_sa_right not in seen:
                    seen.add(key_sa_right)
                    right_pool[ref_name].append((sa_end, sa_right_len))

    samfile.close()

    def _cluster_and_aggregate(pool, side):
        rows = []
        positions_sorted = sorted(pool, key=lambda x: x[0])
        if not positions_sorted:
            return rows
        clusters, current = [], [positions_sorted[0]]
        for item in positions_sorted[1:]:
            if item[0] - current[-1][0] <= cluster_dist:
                current.append(item)
            else:
                clusters.append(current)
                current = [item]
        clusters.append(current)
        for cluster in clusters:
            support     = len(cluster)
            all_lengths = [cl for _, cl in cluster]
            mode_pos    = Counter(p for p, _ in cluster).most_common(1)[0][0]
            if support >= min_support:
                rows.append({
                    "side":          side,
                    "position":      mode_pos,
                    "support":       support,
                    "mean_clip_len": round(sum(all_lengths) / len(all_lengths), 1),
                    "max_clip_len":  max(all_lengths),
                })
        return rows

    rows = []
    for contig in set(left_pool) | set(right_pool):
        for side, pool in (("left", left_pool[contig]), ("right", right_pool[contig])):
            for row in _cluster_and_aggregate(pool, side):
                row["chrom"] = contig
                rows.append(row)

    if not rows:
        print("No major clip positions found with the given thresholds.")
        return pd.DataFrame()

    df = (
        pd.DataFrame(rows, columns=["chrom", "position", "side", "support", "mean_clip_len", "max_clip_len"])
        .sort_values("support", ascending=False)
        .reset_index(drop=True)
    )

    if min_support_frac and not df.empty:
        depth_cache = {}
        keep_mask   = []
        for _, row in df.iterrows():
            key = (row["chrom"], row["position"])
            if key not in depth_cache:
                depth_cache[key] = _region_depth(row["chrom"], row["position"])
            total_depth   = depth_cache[key]
            effective_min = max(min_support, min_support_frac * total_depth)
            passes        = row["support"] >= effective_min
            keep_mask.append(passes)
            if not passes:
                print(
                    f"  [depth-filter] {row['chrom']}:{row['position']} "
                    f"support={row['support']} depth={total_depth} "
                    f"required≥{effective_min:.1f} → DROPPED"
                )
        df = df[keep_mask].reset_index(drop=True)

    if df.empty:
        print("No clip positions passed the fractional-support filter.")
        return df

    if top_n:
        df = df.head(top_n)

    csv_out = os.path.join(save_dir, "major_clip_positions.csv") if save_dir else "major_clip_positions.csv"
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    df.to_csv(csv_out, index=False)
    print(f"Saved: {csv_out}")
    return df

# ─────────────────────────────────────────────────────────────────────────────
# 2. extract_split_reads (unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def extract_split_reads(bam_path, chrom, start_search, end_search, padding):
    samfile     = pysam.AlignmentFile(bam_path, "rb")
    split_reads = []
    for read in samfile.fetch(chrom, max(0, start_search - padding), end_search + padding):
        if read.is_unmapped or not read.has_tag("SA"):
            continue
        for i, entry in enumerate(read.get_tag("SA").rstrip(";").split(";")):
            if i >= 2:
                break
            fields = entry.split(",")
            if len(fields) < 4:
                continue
            sa_chrom, sa_pos = fields[0], int(fields[1])
            if sa_chrom == chrom and (
                abs(sa_pos - start_search) < padding
                or abs(sa_pos - end_search) < padding
            ):
                split_reads.append(read)
                break
    samfile.close()
    return split_reads

# ─────────────────────────────────────────────────────────────────────────────
# 3. pair_clip_positions (unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def _reads_supporting_pair(bam_path, chrom, pos_a, pos_b, tol=200, min_mapq=0):
    bridging = set()

    def _clips_near(cigar_ops, ref_start, target):
        if not cigar_ops:
            return False
        ref_end = ref_start + _ref_span_from_cigar(cigar_ops)
        return (
            (cigar_ops[-1][0] in (4, 5) and abs(ref_end   - target) <= tol)
            or (cigar_ops[0][0] in (4, 5) and abs(ref_start - target) <= tol)
        )

    def _sa_near(read, target):
        if not read.has_tag("SA"):
            return False
        for entry in read.get_tag("SA").rstrip(";").split(";"):
            fields = entry.split(",")
            if len(fields) < 4 or fields[0] != chrom:
                continue
            sa_pos = int(fields[1])
            sa_end = sa_pos + _sa_ref_span(fields[3])
            ops    = re.findall(r"(\d+)([MIDNSHP=X])", fields[3])
            if not ops:
                continue
            if ops[0][1] in ("S", "H") and abs(sa_pos - target) <= tol:
                return True
            if ops[-1][1] in ("S", "H") and abs(sa_end - target) <= tol:
                return True
        return False

    fetch_start = max(0, min(pos_a, pos_b) - tol)
    fetch_end   = max(pos_a, pos_b) + tol
    with pysam.AlignmentFile(bam_path, "rb") as sf:
        for read in sf.fetch(chrom, fetch_start, fetch_end):
            if read.is_unmapped or not read.cigartuples or read.is_secondary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            clips_a = _clips_near(read.cigartuples, read.reference_start, pos_a)
            sa_b    = _sa_near(read, pos_b)
            clips_b = _clips_near(read.cigartuples, read.reference_start, pos_b)
            sa_a    = _sa_near(read, pos_a)
            if (clips_a and sa_b) or (clips_b and sa_a):
                bridging.add(read.query_name)
    return bridging

def count_bridging_reads(bam_path, chrom, pos_a, pos_b, tol=200, min_mapq=0):
    return len(_reads_supporting_pair(bam_path, chrom, pos_a, pos_b, tol=tol, min_mapq=min_mapq))

def pair_clip_positions(
    clip_df,
    bam_path=None,
    max_dist=None,
    strategy="split_read_bridge",
    bridge_tol=200,
):
    pairs = []
    if strategy == "split_read_bridge":
        if bam_path is None:
            return pair_clip_positions(clip_df, max_dist=max_dist, strategy="greedy_rank")
        idxs       = list(clip_df.index)
        n          = len(idxs)
        candidates = []
        print(f"  Scoring {n*(n-1)//2} candidate pairs by bridging read count ...")
        for ii in range(n):
            for jj in range(ii + 1, n):
                i, j = idxs[ii], idxs[jj]
                a, b = clip_df.loc[i], clip_df.loc[j]
                if a["chrom"] != b["chrom"]:
                    continue
                dist = abs(int(b["position"]) - int(a["position"]))
                if max_dist is not None and dist > max_dist:
                    continue
                if dist == 0:
                    continue
                bc = count_bridging_reads(
                    bam_path, a["chrom"],
                    int(a["position"]), int(b["position"]),
                    tol=bridge_tol,
                )
                candidates.append((bc, i, j))
                print(f"  {a['chrom']}:{int(a['position']):,} ↔ {int(b['position']):,} "
                      f"bridge_reads={bc}")
        candidates.sort(key=lambda x: -x[0])
        used = set()
        for bc, i, j in candidates:
            if i in used or j in used or bc == 0:
                continue
            used.add(i); used.add(j)
            a, b  = clip_df.loc[i], clip_df.loc[j]
            left  = int(min(a["position"], b["position"]))
            right = int(max(a["position"], b["position"]))
            pairs.append((a["chrom"], left, right))
        if not pairs:
            print("  [split_read_bridge] No pairs found; falling back to greedy_rank.")
            return pair_clip_positions(clip_df, max_dist=max_dist, strategy="greedy_rank")

    elif strategy == "greedy_rank":
        used = set()
        for i, a in clip_df.iterrows():
            if i in used:
                continue
            best_j = None
            for j, b in clip_df.iterrows():
                if j <= i or j in used or b["chrom"] != a["chrom"]:
                    continue
                dist = abs(int(b["position"]) - int(a["position"]))
                if max_dist is not None and dist > max_dist:
                    continue
                best_j = j
                break
            if best_j is None:
                continue
            used.add(i); used.add(best_j)
            b     = clip_df.loc[best_j]
            left  = int(min(a["position"], b["position"]))
            right = int(max(a["position"], b["position"]))
            pairs.append((a["chrom"], left, right))

    elif strategy == "consecutive":
        for i in range(0, len(clip_df) - 1, 2):
            a, b = clip_df.iloc[i], clip_df.iloc[i + 1]
            if a["chrom"] != b["chrom"]:
                continue
            left  = int(min(a["position"], b["position"]))
            right = int(max(a["position"], b["position"]))
            if max_dist is None or (right - left) <= max_dist:
                pairs.append((a["chrom"], left, right))

    elif strategy == "proximity":
        right_clips = clip_df[clip_df["side"] == "right"].copy()
        left_clips  = clip_df[clip_df["side"] == "left"].copy()
        used = set()
        for _, rc in right_clips.iterrows():
            cands = left_clips[
                (left_clips["chrom"] == rc["chrom"])
                & (left_clips["position"] > rc["position"])
                & (~left_clips.index.isin(used))
            ]
            if max_dist is not None:
                cands = cands[cands["position"] - rc["position"] <= max_dist]
            if cands.empty:
                continue
            best = cands.iloc[0]
            used.add(best.name)
            pairs.append((rc["chrom"], int(rc["position"]), int(best["position"])))

    return pairs

# ─────────────────────────────────────────────────────────────────────────────
# 4. cluster_softclip_pairs (unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def cluster_softclip_pairs(
    pairs,
    bam_path,
    share_threshold=0.50,
    bridge_tol=200,
    min_mapq=0,
):
    if not pairs:
        return []

    pair_reads = []
    for chrom, lp, rp in pairs:
        rset = _reads_supporting_pair(bam_path, chrom, lp, rp, tol=bridge_tol, min_mapq=min_mapq)
        pair_reads.append(rset)
        print(f"  Pair {chrom}:{lp:,}–{rp:,} → {len(rset)} SA reads")

    n      = len(pairs)
    parent = list(range(n))

    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(x, y):
        px, py = _find(x), _find(y)
        if px != py:
            parent[px] = py

    for i in range(n):
        for j in range(i + 1, n):
            si, sj   = pair_reads[i], pair_reads[j]
            if not si or not sj:
                continue
            shared   = si & sj
            min_size = min(len(si), len(sj))
            frac     = len(shared) / min_size if min_size > 0 else 0.0
            if frac > share_threshold:
                print(f"  Merging pair {i} and {j}: "
                      f"{len(shared)}/{min_size} shared reads ({frac:.1%})")
                _union(i, j)

    cluster_map = defaultdict(list)
    for i in range(n):
        cluster_map[_find(i)].append(i)

    clusters = []
    for root, members in cluster_map.items():
        chrom      = pairs[members[0]][0]
        all_left   = [pairs[m][1] for m in members]
        all_right  = [pairs[m][2] for m in members]
        merged_reads = set()
        for m in members:
            merged_reads |= pair_reads[m]
        graph_adj = defaultdict(set)
        for m in members:
            chrom_m, lp_m, rp_m = pairs[m]
            for rname in pair_reads[m]:
                graph_adj[rname].add((chrom_m, lp_m))
                graph_adj[rname].add((chrom_m, rp_m))
        clusters.append({
            "pairs":           [pairs[m] for m in members],
            "reads":           merged_reads,
            "left_pos":        min(all_left),
            "right_pos":       max(all_right),
            "chrom":           chrom,
            "graph_adjacency": dict(graph_adj),
        })

    print(f"\n{len(clusters)} cluster(s) from {n} pair(s).")
    return clusters

# ─────────────────────────────────────────────────────────────────────────────
# 5. compute_reference_window (unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def compute_reference_window(
    reads,
    bam_path: str,
    chrom: str,
    cluster_left: int,
    cluster_right: int,
    fetch_padding: int = 500_000,
    max_span: int = 1_000_000,
):
    read_names  = {r.query_name for r in reads if not r.is_unmapped}
    fetch_start = max(0, cluster_left - fetch_padding)
    fetch_end   = cluster_right + fetch_padding
    coords      = []

    with pysam.AlignmentFile(bam_path, "rb") as sf:
        chrom_len = sf.get_reference_length(chrom)
        fetch_end = min(chrom_len, fetch_end)
        for aln in sf.fetch(chrom, fetch_start, fetch_end):
            if aln.query_name not in read_names:
                continue
            if aln.is_unmapped or not aln.cigartuples or aln.is_secondary:
                continue
            g_start = aln.reference_start
            g_end   = g_start + _ref_span_from_cigar(aln.cigartuples)
            if g_end - g_start <= max_span:
                coords.append((g_start, g_end))

    for read in reads:
        if read.is_unmapped or not read.has_tag("SA"):
            continue
        for sa_chrom, sa_pos, sa_strand, sa_cigar, sa_mapq in _parse_sa_tag(read.get_tag("SA")):
            if sa_chrom != chrom:
                continue
            span = _sa_ref_span(sa_cigar)
            if span <= max_span:
                coords.append((sa_pos - 1, sa_pos - 1 + span))

    if not coords:
        print(f"  [compute_reference_window] No valid alignments on {chrom}; using cluster bounds.")
        return cluster_left, cluster_right

    win_start = min(c[0] for c in coords)
    win_end   = max(c[1] for c in coords)
    print(f"  [compute_reference_window] {chrom}:{win_start:,}–{win_end:,} "
          f"({win_end - win_start:,} bp) from {len(coords)} alignments")
    return win_start, win_end

# ─────────────────────────────────────────────────────────────────────────────
# 6. get_background_homozygous_sites  NEW
#    Pileup all reads (not just split reads) in the cluster window.
#    Return positions where a non-reference base has VAF ≥ 0.95.
# ─────────────────────────────────────────────────────────────────────────────

def get_background_homozygous_sites(
    bam_path: str,
    ref_path: str,
    chrom: str,
    win_start: int,
    win_end: int,
    min_depth: int = 5,
    vaf_threshold: float = 0.95,
) -> set:
    """Return set of (chrom, pos) where all reads carry an alt base (VAF ≥ vaf_threshold).

    These are constitutive sample-specific variants that should be excluded from
    split-read-specific PDV calling.
    """
    hom_sites = set()
    window    = min(win_end - win_start, 5_000_000)  # pysam pileup cap

    try:
        fasta = pysam.FastaFile(ref_path)
        with pysam.AlignmentFile(bam_path, "rb") as sf:
            for col in sf.pileup(
                chrom, win_start, win_start + window,
                truncate=True, min_base_quality=0, ignore_overlaps=False,
                fastafile=fasta,
            ):
                pos      = col.reference_pos
                ref_base = fasta.fetch(chrom, pos, pos + 1).upper()
                bases    = [b.upper() for b in col.get_query_sequences(add_indels=False)]
                depth    = len(bases)
                if depth < min_depth:
                    continue
                alt_counts = Counter(b for b in bases if b and b != ref_base and b != "N")
                if not alt_counts:
                    continue
                top_alt, top_count = alt_counts.most_common(1)[0]
                vaf = top_count / depth
                if vaf >= vaf_threshold:
                    hom_sites.add((chrom, pos))
        fasta.close()
    except Exception as e:
        print(f"  [get_background_homozygous_sites] Warning: {e}")
    print(f"  [get_background_homozygous_sites] {len(hom_sites)} background-hom sites in window")
    return hom_sites

# ─────────────────────────────────────────────────────────────────────────────
# 7. find_split_read_specific_variants  NEW
#    At each position overlapped by split reads in the cluster window:
#      - count bases in split reads
#      - count bases in all reads
#      - mask background-hom sites
#      - keep positions where alt enriched in split reads and not background-hom
#    Returns dict: {(chrom, pos): {"split_vaf": float, "all_vaf": float, "alt": str}}
# ─────────────────────────────────────────────────────────────────────────────

def find_split_read_specific_variants(
    bam_path: str,
    ref_path: str,
    chrom: str,
    win_start: int,
    win_end: int,
    split_read_names: set,
    background_hom_sites: set,
    min_split_depth: int = 2,
    min_split_vaf: float = 0.20,
    min_enrichment: float = 3.0,
) -> dict:
    """Return PDV positions enriched in split reads but absent from background-hom set."""
    pdv_sites = {}
    try:
        fasta = pysam.FastaFile(ref_path)
        # Two-pass: collect split-read bases and all-read bases per position
        split_bases = defaultdict(list)
        all_bases   = defaultdict(list)

        with pysam.AlignmentFile(bam_path, "rb") as sf:
            for col in sf.pileup(
                chrom, win_start, win_end,
                truncate=True, min_base_quality=0, ignore_overlaps=False,
            ):
                pos = col.reference_pos
                if (chrom, pos) in background_hom_sites:
                    continue
                ref_base = fasta.fetch(chrom, pos, pos + 1).upper()
                for read_pileup in col.pileups:
                    if read_pileup.is_del or read_pileup.is_refskip:
                        continue
                    base = read_pileup.alignment.query_sequence[read_pileup.query_position].upper()
                    all_bases[pos].append(base)
                    if read_pileup.alignment.query_name in split_read_names:
                        split_bases[pos].append(base)

        for pos, sbases in split_bases.items():
            if len(sbases) < min_split_depth:
                continue
            ref_base = fasta.fetch(chrom, pos, pos + 1).upper()
            s_alts   = Counter(b for b in sbases if b != ref_base and b != "N")
            if not s_alts:
                continue
            top_alt, s_count = s_alts.most_common(1)[0]
            s_vaf = s_count / len(sbases)
            if s_vaf < min_split_vaf:
                continue
            a_bases = all_bases.get(pos, [])
            a_depth = len(a_bases)
            if a_depth == 0:
                continue
            a_count = sum(1 for b in a_bases if b == top_alt)
            a_vaf   = a_count / a_depth
            # Require enrichment in split reads relative to all reads
            if a_vaf > 0 and (s_vaf / a_vaf) < min_enrichment:
                continue
            pdv_sites[pos] = {"split_vaf": s_vaf, "all_vaf": a_vaf, "alt": top_alt}
        fasta.close()
    except Exception as e:
        print(f"  [find_split_read_specific_variants] Warning: {e}")

    print(f"  [find_split_read_specific_variants] {len(pdv_sites)} split-read-specific PDV positions")
    return pdv_sites

# ─────────────────────────────────────────────────────────────────────────────
# 8. collect_breakpoint_segments  NEW  (replaces collect_clip_sequences)
#    For each split read:
#      a. Identify the soft-clip end (left-clip tail or right-clip tail)
#      b. Walk inward along the aligned portion
#      c. Stop at the last PDV position (split-read-specific, not background-hom)
#      d. Return the query subsequence from clip tail through last PDV
# ─────────────────────────────────────────────────────────────────────────────

def _get_primary_read(bam_path, chrom, read_name, win_start, win_end, fetch_padding=500_000):
    """Fetch the primary alignment record for a given read name."""
    with pysam.AlignmentFile(bam_path, "rb") as sf:
        for aln in sf.fetch(chrom, max(0, win_start - fetch_padding), win_end + fetch_padding):
            if aln.query_name == read_name and not aln.is_secondary and not aln.is_supplementary:
                return aln
    return None

def collect_breakpoint_segments(
    reads,
    bam_path: str,
    chrom: str,
    win_start: int,
    win_end: int,
    pdv_sites: dict,
    min_clip_len: int = 10,
    fetch_padding: int = 500_000,
) -> list:
    """Return list of dicts with keys: read_name, segment_seq, clip_side, clip_pos.

    Each segment runs from the soft-clip tail of a breakpoint-bearing read inward
    through the last split-read-specific PDV within the read's aligned portion.
    If no PDV is found in the aligned portion, only the raw soft-clip bases are used.
    """
    segments  = []
    read_names = {r.query_name for r in reads if not r.is_unmapped}

    primary_index = {}
    with pysam.AlignmentFile(bam_path, "rb") as sf:
        for aln in sf.fetch(chrom, max(0, win_start - fetch_padding), win_end + fetch_padding):
            if aln.query_name not in read_names:
                continue
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            if aln.query_sequence is None:
                continue
            primary_index[aln.query_name] = aln

    for read in reads:
        if read.is_unmapped or not read.cigartuples:
            continue
        prim = primary_index.get(read.query_name)
        if prim is None:
            prim = read
        if prim.query_sequence is None:
            continue

        seq    = prim.query_sequence
        cigars = prim.cigartuples
        ref_s  = prim.reference_start

        # Build a map: query_pos → ref_pos for the aligned span
        q2r = {}
        q_cursor = 0
        r_cursor = ref_s
        for op, length in cigars:
            if op == 4:   # soft clip
                q_cursor += length
            elif op in (0, 7, 8):  # M / = / X
                for k in range(length):
                    q2r[q_cursor + k] = r_cursor + k
                q_cursor += length
                r_cursor += length
            elif op in (1,):  # insertion
                q_cursor += length
            elif op in (2, 3):  # deletion / skip
                r_cursor += length
            elif op in (5,):   # hard clip
                pass

        # Identify clip ends
        first_op, first_len = cigars[0]
        last_op,  last_len  = cigars[-1]

        for clip_side, clip_len, clip_op in (
            ("left",  first_len, first_op),
            ("right", last_len,  last_op),
        ):
            if clip_op != 4 or clip_len < min_clip_len:
                continue

            # Raw clip sequence (0-based query coordinates)
            if clip_side == "left":
                clip_q_start = 0
                clip_q_end   = first_len
                clip_ref_pos = ref_s
                # Walk inward: query positions after the clip, ascending ref positions
                aligned_q_positions = sorted(
                    (qp for qp, rp in q2r.items() if qp >= first_len),
                    key=lambda qp: q2r[qp]
                )
                # Find PDV positions within the aligned span, sorted by distance from clip
                pdv_in_span = [
                    (q2r[qp], qp) for qp in aligned_q_positions
                    if q2r[qp] in pdv_sites
                ]
                if pdv_in_span:
                    # Last PDV = furthest from the clip end (highest ref pos for left-clip)
                    _, last_pdv_qpos = pdv_in_span[-1]
                    seg_q_start = clip_q_start
                    seg_q_end   = last_pdv_qpos + 1
                else:
                    seg_q_start = clip_q_start
                    seg_q_end   = clip_q_end
            else:  # right
                clip_q_start = len(seq) - last_len
                clip_q_end   = len(seq)
                clip_ref_pos = prim.reference_end
                aligned_q_positions = sorted(
                    (qp for qp, rp in q2r.items() if qp < clip_q_start),
                    key=lambda qp: -q2r[qp]
                )
                pdv_in_span = [
                    (q2r[qp], qp) for qp in aligned_q_positions
                    if q2r[qp] in pdv_sites
                ]
                if pdv_in_span:
                    _, last_pdv_qpos = pdv_in_span[-1]
                    seg_q_start = last_pdv_qpos
                    seg_q_end   = clip_q_end
                else:
                    seg_q_start = clip_q_start
                    seg_q_end   = clip_q_end

            segment_seq = seq[seg_q_start:seg_q_end]
            if len(segment_seq) < min_clip_len:
                continue

            if prim.is_reverse:
                segment_seq = reverse_complement(segment_seq)

            segments.append({
                "read_name":   prim.query_name,
                "segment_seq": segment_seq,
                "clip_side":   clip_side,
                "clip_pos":    clip_ref_pos,
                "n_pdvs":      len(pdv_in_span),
            })

    print(f"  [collect_breakpoint_segments] {len(segments)} segments from {len(reads)} reads")
    return segments

# ─────────────────────────────────────────────────────────────────────────────
# 9. align_clips_to_ref  NEW
#    Run minimap2 (subprocess) for each breakpoint segment against the reference.
#    Falls back to edlib windowed scan if minimap2 unavailable.
# ─────────────────────────────────────────────────────────────────────────────

def _run_minimap2(
    query_fasta: str,
    ref_path: str,
    preset: str = "map-hifi",
    match_score: int = 2,
    mismatch_penalty: int = 8,
    gap_open: int = 4,
    gap_extend: int = 2,
    gap_open2: int = 24,
    gap_extend2: int = 1,
) -> list:
    """Run minimap2 with stringent mismatch scoring and return list of PAF lines.

    Scoring rationale for LCR breakpoint segments
    -----------------------------------------------
    High mismatch penalty (-B 8):
        Single-base substitutions cost 8 points, so a true SNP between paralog
        copies creates a large alignment score difference. This means minimap2
        will only report an alignment locus if the segment is genuinely a close
        match (one copy's signature) rather than averaging across copies.

    Lenient homopolymer gap scoring (-O 4,24 -E 2,1):
        Two-piece affine gap model. minimap2 enforces two hard constraints:
          (1) E1 > E2  →  gap_extend (2) > gap_extend2 (1)  ✓
          (2) O1+E1 < O2+E2  →  4+2=6 < 24+1=25           ✓
        Short gaps (homopolymer indels) open cheaply at 4 with extend 2/bp.
        Long gaps open expensively at 24 to prevent bridging across repeat arms.

    Minimal match bonus (-A 2):
        Kept low relative to the mismatch penalty (ratio 1:4) so a single
        mismatch is not easily offset by nearby matches.
    """
    cmd = [
        "minimap2",
        "-c",                            # output CIGAR in PAF
        "--secondary=yes",               # report secondary hits (needed for multi-locus view)
        f"-x{preset}",                   # long-read preset base (sets seed/chain params)
        f"-A{match_score}",              # match bonus
        f"-B{mismatch_penalty}",         # mismatch penalty (high → discriminating)
        f"-O{gap_open},{gap_open2}",     # two-piece affine gap open  (short, long)
        f"-E{gap_extend},{gap_extend2}", # two-piece affine gap extend (short, long)
        ref_path,
        query_fasta,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            print(f"  [minimap2] Warning: {result.stderr[:200]}")
            return []
        return [l for l in result.stdout.splitlines() if l.strip()]
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f"  [minimap2] Not available or timed out: {e}")
        return []

def _parse_paf(paf_lines: list) -> list:
    """Parse PAF lines → list of dicts, including per-hit sequence identity."""
    hits = []
    for line in paf_lines:
        f = line.split("\t")
        if len(f) < 12:
            continue
        try:
            n_matches  = int(f[9])
            block_len  = int(f[10])
            # PAF col 10 (n_matches) counts residue matches; identity = matches/block_len
            identity   = n_matches / block_len if block_len > 0 else 0.0
            # Also extract de:f: tag (divergence) if present for extra precision
            for tag in f[12:]:
                if tag.startswith("de:f:"):
                    try:
                        identity = 1.0 - float(tag[5:])
                    except ValueError:
                        pass
                    break
            hits.append({
                "query_name":  f[0],
                "query_len":   int(f[1]),
                "query_start": int(f[2]),
                "query_end":   int(f[3]),
                "strand":      f[4],
                "ref_chrom":   f[5],
                "ref_start":   int(f[7]),
                "ref_end":     int(f[8]),
                "n_matches":   n_matches,
                "block_len":   block_len,
                "mapq":        int(f[11]),
                "identity":    identity,
            })
        except (ValueError, IndexError):
            continue
    return hits

def align_clips_to_ref(
    segments: list,
    ref_path: str,
    minimap2_preset: str = "map-hifi",
    max_hits_per_segment: int = 5,
    min_identity: float = 0.90,
    match_score: int = 2,
    mismatch_penalty: int = 8,
    gap_open: int = 4,
    gap_extend: int = 2,
    gap_open2: int = 24,
    gap_extend2: int = 1,
) -> list:
    """Align each breakpoint segment to the reference with stringent mismatch scoring.

    Alignment scoring (passed to minimap2)
    ----------------------------------------
    mismatch_penalty=8  – high cost per SNP so paralog-copy-specific substitutions
                          force the aligner to commit to one locus rather than bridging
    gap_open=4          – low open cost for short gaps (homopolymer length variation)
    gap_extend=1        – low per-base cost for short gaps
    gap_open2=24        – high open cost for long gaps keeps the alignment local
    gap_extend2=1       – identical per-base cost for long gap extension

    Filtering
    ---------
    Hits with sequence identity < min_identity are discarded after alignment so
    that only near-perfect reference matches are kept for graph construction.
    Identity is derived from the PAF de:f: divergence tag (preferred) or
    n_matches/block_len as fallback.

    Returns list of ClipAlignment namedtuples (pdvs field populated later by call_pdvs).
    """
    if not segments:
        return []

    alignments = []
    kept = dropped = 0

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
        for i, seg in enumerate(segments):
            fh.write(f">{i}_{seg['read_name']}_{seg['clip_side']}\n{seg['segment_seq']}\n")
        tmp_fa = fh.name

    paf_lines = _run_minimap2(
        tmp_fa, ref_path,
        preset          = minimap2_preset,
        match_score     = match_score,
        mismatch_penalty= mismatch_penalty,
        gap_open        = gap_open,
        gap_extend      = gap_extend,
        gap_open2       = gap_open2,
        gap_extend2     = gap_extend2,
    )
    os.unlink(tmp_fa)

    hits_by_query = defaultdict(list)
    for h in _parse_paf(paf_lines):
        hits_by_query[h["query_name"]].append(h)

    for i, seg in enumerate(segments):
        prefix = f"{i}_{seg['read_name']}_{seg['clip_side']}"
        hits   = hits_by_query.get(prefix, [])
        # Sort by identity desc, then mapq desc
        hits   = sorted(hits, key=lambda h: (-h["identity"], -h["mapq"]))[:max_hits_per_segment]
        if not hits:
            continue
        for h in hits:
            if h["identity"] < min_identity:
                dropped += 1
                continue
            n_mm = h["block_len"] - h["n_matches"]
            kept += 1
            alignments.append(ClipAlignment(
                read_name    = seg["read_name"],
                segment_seq  = seg["segment_seq"],
                ref_chrom    = h["ref_chrom"],
                ref_start    = h["ref_start"],
                ref_end      = h["ref_end"],
                strand       = h["strand"],
                n_mismatches = n_mm,
                mapq         = h["mapq"],
                pdvs         = [],
            ))

    print(f"  [align_clips_to_ref] {kept} alignments kept, {dropped} dropped "
          f"(identity<{min_identity}) from {len(segments)} segments")
    return alignments

# ─────────────────────────────────────────────────────────────────────────────
# 10. call_pdvs  NEW
#     At each alignment position, compare:
#       clip base vs. ref base vs. all other alignment positions (other loci)
#     Flag bases unique to one locus as PDV and annotate status.
# ─────────────────────────────────────────────────────────────────────────────

def call_pdvs(
    alignments: list,
    ref_path: str,
    pdv_sites: dict,
) -> list:
    """Annotate each ClipAlignment with PDVAnnotation objects.

    For each alignment, fetch the reference sequence in the aligned window,
    compare clip bases to ref, then cross-reference known PDV positions.
    Returns a new list of ClipAlignment with pdvs populated.
    """
    if not alignments:
        return alignments

    try:
        fasta = pysam.FastaFile(ref_path)
    except Exception:
        return alignments

    annotated = []
    for aln in alignments:
        try:
            ref_seq = fasta.fetch(aln.ref_chrom, aln.ref_start, aln.ref_end).upper()
        except Exception:
            annotated.append(aln)
            continue

        seg_seq = aln.segment_seq.upper()
        if aln.strand == "-":
            seg_seq = reverse_complement(seg_seq)

        pdv_annotations = []
        seg_len = min(len(seg_seq), len(ref_seq), aln.ref_end - aln.ref_start)
        for k in range(seg_len):
            ref_pos  = aln.ref_start + k
            if ref_pos not in pdv_sites:
                continue
            ref_base  = ref_seq[k] if k < len(ref_seq) else "N"
            clip_base = seg_seq[k] if k < len(seg_seq) else "N"
            info      = pdv_sites[ref_pos]
            # Determine status
            if clip_base == info["alt"]:
                status = "source" if info["split_vaf"] >= 0.8 else "dest"
            elif clip_base == ref_base:
                status = "dest"
            else:
                status = "novel"
            pdv_annotations.append(PDVAnnotation(
                ref_pos    = ref_pos,
                clip_base  = clip_base,
                ref_base   = ref_base,
                split_vaf  = info["split_vaf"],
                all_vaf    = info["all_vaf"],
                status     = status,
            ))

        annotated.append(ClipAlignment(
            read_name    = aln.read_name,
            segment_seq  = aln.segment_seq,
            ref_chrom    = aln.ref_chrom,
            ref_start    = aln.ref_start,
            ref_end      = aln.ref_end,
            strand       = aln.strand,
            n_mismatches = aln.n_mismatches,
            mapq         = aln.mapq,
            pdvs         = pdv_annotations,
        ))

    fasta.close()
    return annotated


def _simple_consensus(seqs: list) -> str:
    """Simple ungapped majority-vote consensus across similarly sized sequences."""
    if not seqs:
        return ""
    max_len = max(len(s) for s in seqs)
    cons = []
    for i in range(max_len):
        col = [s[i].upper() for s in seqs if i < len(s) and s[i].upper() in "ACGTN"]
        if not col:
            cons.append("N")
            continue
        counts = Counter(b for b in col if b != "N")
        cons.append(counts.most_common(1)[0][0] if counts else "N")
    return "".join(cons)


def write_segment_consensus_fasta(
    locus_clusters: list,
    background_hom_sites: set,
    pdv_sites: dict,
    save_dir: str,
    pair_label: str,
    min_graph_segment_support: int = 10,
) -> str:
    """Write one FASTA containing per-segment consensus plus annotation masks.

    For each graph-eligible locus cluster (read support >= min_graph_segment_support),
    write:
      1. consensus sequence
      2. annotation mask using '.' = regular, 'P' = PDV, 'H' = background-homozygous
    """
    os.makedirs(save_dir, exist_ok=True)
    fasta_path = os.path.join(save_dir, f"segment_consensus_{pair_label}.fasta")

    eligible = [lc for lc in locus_clusters if len(lc["read_names"]) >= min_graph_segment_support]
    with open(fasta_path, "w") as fh:
        for i, lc in enumerate(eligible, start=1):
            rep = sorted(lc["alignments"], key=lambda a: (-a.mapq, a.n_mismatches))[0]
            seqs = []
            for aln in lc["alignments"]:
                seq = aln.segment_seq.upper()
                if aln.strand != lc["strand"]:
                    seq = reverse_complement(seq)
                seqs.append(seq)
            consensus = _simple_consensus(seqs)

            mask = []
            for k in range(len(consensus)):
                ref_pos = rep.ref_start + k if lc["strand"] == "+" else rep.ref_end - 1 - k
                if ref_pos in pdv_sites:
                    mask.append("P")
                elif (lc["chrom"], ref_pos) in background_hom_sites:
                    mask.append("H")
                else:
                    mask.append(".")
            mask = "".join(mask)
            pdv_positions = [str(p) for p in sorted(pos for pos in pdv_sites if lc["locus_start"] <= pos < lc["locus_end"])]
            hom_positions = [str(pos) for c, pos in sorted(background_hom_sites) if c == lc["chrom"] and lc["locus_start"] <= pos < lc["locus_end"]]
            header = (
                f">segment_{i} {lc['chrom']}:{lc['locus_start']}-{lc['locus_end']} "
                f"strand={lc['strand']} n={len(lc['read_names'])} "
                f"pdv_count={mask.count('P')} hom_count={mask.count('H')}"
            )
            fh.write(header + "\n")
            fh.write(consensus + "\n")
            fh.write(f">segment_{i}_mask legend=.:regular,P:PDV,H:background_hom\n")
            fh.write(mask + "\n")
            fh.write(f">segment_{i}_pdv_positions\n")
            fh.write((",".join(pdv_positions) if pdv_positions else "none") + "\n")
            fh.write(f">segment_{i}_background_hom_positions\n")
            fh.write((",".join(hom_positions) if hom_positions else "none") + "\n")

    print(f"  [write_segment_consensus_fasta] Wrote {len(eligible)} segment consensus record(s) → {fasta_path}")
    return fasta_path
    return fasta_path

# ─────────────────────────────────────────────────────────────────────────────
# 11. cluster_by_landing_site  NEW  (replaces Jaccard / junction-key clustering)
#     Group reads whose clip segments land at the same reference locus (±tol).
# ─────────────────────────────────────────────────────────────────────────────

def cluster_by_landing_site(
    alignments: list,
    tol: int = 200,
) -> list:
    """Return list of locus-cluster dicts.

    Each dict:
        chrom, locus_start, locus_end, strand,
        alignments: [ClipAlignment, ...],
        read_names: set,
        consensus_pdvs: list of PDVAnnotation (union of all pdvs in cluster),
    """
    if not alignments:
        return []

    # Sort by chrom, ref_start
    sorted_alns = sorted(alignments, key=lambda a: (a.ref_chrom, a.ref_start))
    clusters    = []
    current     = [sorted_alns[0]]

    def _flush(group):
        chroms  = [a.ref_chrom for a in group]
        starts  = [a.ref_start for a in group]
        ends    = [a.ref_end   for a in group]
        strands = Counter(a.strand for a in group)
        dom_strand = strands.most_common(1)[0][0]
        # Merge PDV annotations across all reads in cluster
        pdv_map = {}
        for a in group:
            for p in a.pdvs:
                if p.ref_pos not in pdv_map:
                    pdv_map[p.ref_pos] = p
        return {
            "chrom":          chroms[0],
            "locus_start":    min(starts),
            "locus_end":      max(ends),
            "strand":         dom_strand,
            "alignments":     group,
            "read_names":     {a.read_name for a in group},
            "consensus_pdvs": list(pdv_map.values()),
        }

    for aln in sorted_alns[1:]:
        last = current[-1]
        if aln.ref_chrom == last.ref_chrom and (aln.ref_start - last.ref_end) <= tol:
            current.append(aln)
        else:
            clusters.append(_flush(current))
            current = [aln]
    clusters.append(_flush(current))

    print(f"  [cluster_by_landing_site] {len(clusters)} locus cluster(s) from {len(alignments)} alignments")
    return clusters

# ─────────────────────────────────────────────────────────────────────────────
# 12. build_path_from_clips  NEW  (replaces build_path_from_reads)
#     Each unique landing locus = one PathSegment node in the rearrangement graph.
# ─────────────────────────────────────────────────────────────────────────────

def build_path_from_clips(
    locus_clusters: list,
    chrom: str,
    cluster_left: int,
    cluster_right: int,
    min_graph_segment_support: int = 10,
) -> list:
    """Convert locus clusters into an ordered list of PathSegment nodes.

    Only locus clusters with read support >= min_graph_segment_support are kept
    for graph construction. Segments are sorted by genomic coordinate so that
    the rearrangement graph flows left → right on the reference.
    """
    kept = [lc for lc in locus_clusters if len(lc["read_names"]) >= min_graph_segment_support]
    dropped = len(locus_clusters) - len(kept)
    if dropped:
        print(f"  [build_path_from_clips] Dropped {dropped} low-support locus cluster(s); "
              f"keeping {len(kept)} with n>={min_graph_segment_support}")
    path = []
    for i, lc in enumerate(kept):
        label = f"Locus{i+1}_{lc['chrom']}:{lc['locus_start']:,}"
        seg = PathSegment(
            query_start  = lc["locus_start"],
            query_end    = lc["locus_end"],
            ref_chrom    = lc["chrom"],
            ref_start    = lc["locus_start"],
            ref_end      = lc["locus_end"],
            strand       = lc["strand"],
            identity     = 1.0 - (
                np.mean([a.n_mismatches / max(1, a.ref_end - a.ref_start)
                         for a in lc["alignments"]])
                if lc["alignments"] else 0
            ),
            read_support = len(lc["read_names"]),
            segment_type = "reference",
            label        = label,
            pdvs         = lc["consensus_pdvs"],
        )
        path.append(seg)

    path.sort(key=lambda s: (s.ref_chrom, s.ref_start))
    return path

# ─────────────────────────────────────────────────────────────────────────────
# 13. plot_rearrangement_path  (extended with PDV edge annotations)
# ─────────────────────────────────────────────────────────────────────────────

def _seg_color(label: str) -> str:
    palette = [
        "#4e79a7", "#f28e2b", "#e15759", "#76b7b2",
        "#59a14f", "#edc948", "#b07aa1", "#ff9da7",
        "#9c755f", "#bab0ac",
    ]
    return palette[hash(label) % len(palette)]

def _pdv_status_color(status: str) -> str:
    return {"source": "#2ca02c", "dest": "#d62728", "novel": "#aaaaaa"}.get(status, "#aaaaaa")

def plot_rearrangement_path(
    path,
    chrom: str,
    left_bp: int,
    right_bp: int,
    sv_type: str,
    pair_label: str = "",
    save_path: str = "rearrangement_path.png",
):
    if not path:
        print("  [plot] Empty path — nothing to plot.")
        return

    fig, ax = plt.subplots(figsize=(max(14, len(path) * 2.5), 6))
    ax.set_xlim(-0.5, len(path) - 0.5)
    ax.set_ylim(-1.5, 3.0)
    ax.axis("off")

    box_h = 0.7
    y_box = 1.0
    cx    = []

    for i, seg in enumerate(path):
        x    = i
        col  = _seg_color(seg.label)
        rect = mpatches.FancyBboxPatch(
            (x - 0.4, y_box - box_h / 2), 0.8, box_h,
            boxstyle="round,pad=0.05",
            facecolor=col, edgecolor="k", linewidth=1.2, alpha=0.85,
        )
        ax.add_patch(rect)
        # Label inside box
        short = seg.label.split("_")[0] if "_" in seg.label else seg.label
        ax.text(x, y_box + 0.05, short, ha="center", va="center",
                fontsize=7, fontweight="bold", color="white", clip_on=True)
        # Coordinates below box
        ax.text(x, y_box - box_h / 2 - 0.18,
                f"{seg.ref_chrom}:{seg.ref_start:,}–{seg.ref_end:,}",
                ha="center", va="top", fontsize=6, color="#444")
        # Support + strand
        strand_sym = "→" if seg.strand == "+" else "←"
        ax.text(x, y_box + box_h / 2 + 0.10,
                f"{strand_sym} n={seg.read_support} id={seg.identity:.2f}",
                ha="center", va="bottom", fontsize=6, color="#333")
        cx.append(x)

    # Draw edges between consecutive nodes
    for i in range(len(path) - 1):
        x1, x2 = cx[i] + 0.4, cx[i + 1] - 0.4
        ax.annotate("", xy=(x2, y_box), xytext=(x1, y_box),
                    arrowprops=dict(arrowstyle="->", color="#555", lw=1.5))
        # PDV annotation on edge
        all_pdvs = path[i].pdvs + path[i + 1].pdvs
        if all_pdvs:
            xm    = (x1 + x2) / 2
            lines = []
            for p in all_pdvs[:4]:  # max 4 PDVs shown per edge
                c = _pdv_status_color(p.status)
                lines.append((f"pos{p.ref_pos} {p.ref_base}→{p.clip_base} "
                               f"sVAF={p.split_vaf:.2f}", c))
            y_ann = y_box + 0.55
            for txt, col in lines:
                ax.text(xm, y_ann, txt, ha="center", va="bottom",
                        fontsize=5, color=col,
                        bbox=dict(boxstyle="round,pad=0.1", fc="white", ec=col, alpha=0.7))
                y_ann += 0.18

    # Title
    title = (f"Rearrangement path — {pair_label}  |  "
             f"SV={sv_type}  {chrom}:{left_bp:,}–{right_bp:,}")
    ax.set_title(title, fontsize=9, pad=8)

    # Legend
    status_patches = [
        mpatches.Patch(color=_pdv_status_color("source"), label="PDV source locus"),
        mpatches.Patch(color=_pdv_status_color("dest"),   label="PDV dest locus"),
        mpatches.Patch(color=_pdv_status_color("novel"),  label="PDV novel"),
    ]
    ax.legend(handles=status_patches, loc="lower right", fontsize=7, framealpha=0.8)

    plt.tight_layout(pad=1.5)
    plt.savefig(save_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved rearrangement path figure → {save_path}")
    try:
        display(Image(save_path))
    except Exception:
        pass

# ─────────────────────────────────────────────────────────────────────────────
# SV type detection helper (unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def _detect_sv_type(reads, refined_left: int, padding: int) -> str:
    tol        = max(50, padding // 10)
    REF_OPS_SA = {"M", "D", "N", "=", "X"}
    left_at_left  = 0
    right_at_left = 0
    for read in reads:
        if read.is_unmapped or not read.cigartuples:
            continue
        cigars = read.cigartuples
        if cigars[0][0] in CLIP_OPS and abs(read.reference_start - refined_left) < tol:
            left_at_left += 1
        if cigars[-1][0] in CLIP_OPS and abs(read.reference_end - refined_left) < tol:
            right_at_left += 1
        if read.has_tag("SA"):
            for sa_entry in read.get_tag("SA").rstrip(";").split(";"):
                fields = sa_entry.split(",")
                if len(fields) < 4 or fields[0] != read.reference_name:
                    continue
                sa_pos    = int(fields[1])
                sa_cigar  = fields[3]
                ops       = re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
                if not ops:
                    continue
                sa_ref_len = sum(int(l) for l, op in ops if op in REF_OPS_SA)
                sa_end     = sa_pos + sa_ref_len
                if ops[0][1] in {"S", "H"}  and abs(sa_pos - refined_left) < tol:
                    left_at_left += 1
                if ops[-1][1] in {"S", "H"} and abs(sa_end - refined_left) < tol:
                    right_at_left += 1
    return "deletion" if right_at_left >= left_at_left else "duplication"

# ─────────────────────────────────────────────────────────────────────────────
# Path summary helper
# ─────────────────────────────────────────────────────────────────────────────

def _summarise_path(path, label=""):
    if not path:
        return f"  {label}: (empty path)"
    parts = []
    for seg in path:
        if seg.segment_type == "novel":
            parts.append(f"[NOVEL {seg.query_end - seg.query_start} bp]")
        else:
            orient = "" if seg.strand == "+" else " (RC)"
            pdv_str = f" PDVs={len(seg.pdvs)}" if seg.pdvs else ""
            parts.append(
                f"[{seg.label}{orient} "
                f"{seg.ref_chrom}:{seg.ref_start:,}–{seg.ref_end:,} "
                f"n={seg.read_support} id={seg.identity:.2f}{pdv_str}]"
            )
    return f"  {label}: " + " → ".join(parts)

# ─────────────────────────────────────────────────────────────────────────────
# Per-cluster entry point  (clip-first pipeline)
# ─────────────────────────────────────────────────────────────────────────────

def analyze_sv_cluster(
    bam_path: str,
    ref_path: str,
    cluster: dict,
    show_pileup: bool = False,
    min_mapq: int = 60,
    pair_label: str = "",
    save_dir: str = ".",
    min_clip_len: int = 50,
    landing_tol: int = 200,
    fetch_padding: int = 500_000,
    max_aln_span: int = 1_000_000,
    minimap2_preset: str = "map-hifi",
    min_identity: float = 0.90,
    match_score: int = 2,
    mismatch_penalty: int = 8,
    gap_open: int = 4,
    gap_extend: int = 2,
    gap_open2: int = 24,
    gap_extend2: int = 1,
    pdv_min_split_vaf: float = 0.20,
    pdv_min_enrichment: float = 3.0,
    background_hom_vaf: float = 0.95,
    min_graph_segment_support: int = 10,
):
    """Run the full clip-first pipeline for one merged cluster.

    Steps
    -----
    1.  Collect all SA reads from the cluster.
    2.  Compute the reference window (left-most → right-most alignment).
    3.  Get background homozygous sites (all reads, VAF≥0.95).
    4.  Find split-read-specific PDV positions (after background masking).
    5.  Collect breakpoint segments (clip tail + inward through last PDV).
    6.  Align segments to reference (minimap2, stringent mismatch scoring).
    7.  Annotate PDVs on each alignment.
    8.  Cluster alignments by landing site.
    9.  Build path from locus clusters.
    10. Plot.
    """
    chrom    = cluster["chrom"]
    left_bp  = cluster["left_pos"]
    right_bp = cluster["right_pos"]
    padding  = max(500, (right_bp - left_bp) // 2)

    # ── 1. Collect reads ─────────────────────────────────────────────────────
    reads = extract_split_reads(
        bam_path, chrom,
        start_search = left_bp,
        end_search   = right_bp,
        padding      = min(padding, left_bp),
    )
    if not reads:
        print(f"[{pair_label}] No split reads found — skipping.")
        return None

    sv_type = _detect_sv_type(reads, left_bp, padding)
    print(f"[{pair_label}] sv_type={sv_type} {chrom}:{left_bp:,}–{right_bp:,} "
          f"cluster_pairs={len(cluster['pairs'])}")

    # ── 2. Reference window ──────────────────────────────────────────────────
    win_start, win_end = compute_reference_window(
        reads, bam_path, chrom,
        cluster_left  = left_bp,
        cluster_right = right_bp,
        fetch_padding = fetch_padding,
        max_span      = max_aln_span,
    )

    # ── 3. Background homozygous sites ───────────────────────────────────────
    bg_hom_sites = get_background_homozygous_sites(
        bam_path, ref_path, chrom, win_start, win_end,
        vaf_threshold = background_hom_vaf,
    )

    # ── 4. Split-read-specific PDVs ──────────────────────────────────────────
    split_names = {r.query_name for r in reads if not r.is_unmapped}
    pdv_sites   = find_split_read_specific_variants(
        bam_path, ref_path, chrom, win_start, win_end,
        split_read_names    = split_names,
        background_hom_sites= bg_hom_sites,
        min_split_vaf       = pdv_min_split_vaf,
        min_enrichment      = pdv_min_enrichment,
    )

    # ── 5. Collect breakpoint segments ───────────────────────────────────────
    segments = collect_breakpoint_segments(
        reads, bam_path, chrom, win_start, win_end,
        pdv_sites     = pdv_sites,
        min_clip_len  = min_clip_len,
        fetch_padding = fetch_padding,
    )
    if not segments:
        print(f"[{pair_label}] No breakpoint segments collected — skipping.")
        return None

    # ── 6. Align segments ────────────────────────────────────────────────────
    alignments = align_clips_to_ref(
        segments, ref_path,
        minimap2_preset  = minimap2_preset,
        min_identity     = min_identity,
        match_score      = match_score,
        mismatch_penalty = mismatch_penalty,
        gap_open         = gap_open,
        gap_extend       = gap_extend,
        gap_open2        = gap_open2,
        gap_extend2      = gap_extend2,
    )
    if not alignments:
        print(f"[{pair_label}] No clip alignments — skipping.")
        return None

    # ── 7. Annotate PDVs ─────────────────────────────────────────────────────
    alignments = call_pdvs(alignments, ref_path, pdv_sites)

    # ── 8. Cluster by landing site ───────────────────────────────────────────
    locus_clusters = cluster_by_landing_site(alignments, tol=landing_tol)

    # ── 9. Build path ────────────────────────────────────────────────────────
    path = build_path_from_clips(locus_clusters, chrom, left_bp, right_bp)

    print("\n" + "─" * 60)
    print("  Rearrangement path summary:")
    print(_summarise_path(path, "PATH"))

    # ── 10. Plot ─────────────────────────────────────────────────────────────
    plot_rearrangement_path(
        path, chrom, left_bp, right_bp, sv_type,
        pair_label = pair_label,
        save_path  = os.path.join(save_dir, f"rearrangement_path_{pair_label}.png"),
    )

    def _path_to_str(p):
        parts = []
        for seg in p:
            if seg.segment_type == "novel":
                parts.append(f"NOVEL({seg.query_end - seg.query_start}bp)")
            else:
                parts.append(
                    f"{seg.ref_chrom}:{seg.ref_start}-{seg.ref_end}"
                    f"({seg.strand})[{seg.label},n={seg.read_support},"
                    f"id={seg.identity:.2f},pdvs={len(seg.pdvs)}]"
                )
        return " → ".join(parts)

    total_pdvs = sum(len(seg.pdvs) for seg in path)

    # ── 11. Write segment consensus FASTA ────────────────────────────────────
    segment_consensus_fasta = write_segment_consensus_fasta(
        locus_clusters      = locus_clusters,
        background_hom_sites = bg_hom_sites,
        pdv_sites            = pdv_sites,
        save_dir             = save_dir,
        pair_label           = pair_label,
        min_graph_segment_support = min_graph_segment_support,
    )

    return {
        "pair_label":       pair_label,
        "chrom":            chrom,
        "sv_type":          sv_type,
        "left_clip_pos":    left_bp,
        "right_clip_pos":   right_bp,
        "n_pairs_merged":   len(cluster["pairs"]),
        "supporting_reads": len(cluster["reads"]),
        "win_start":        win_start,
        "win_end":          win_end,
        "n_segments":       len(segments),
        "n_alignments":     len(alignments),
        "n_loci":           len(locus_clusters),
        "n_graph_loci":     len(path),
        "total_pdvs":       total_pdvs,
        "path":             _path_to_str(path),
        "bg_hom_sites":     len(bg_hom_sites),
        "split_pdv_sites":  len(pdv_sites),
        "segment_consensus_fasta": segment_consensus_fasta,
    }

# ─────────────────────────────────────────────────────────────────────────────
# Batch entry point
# ─────────────────────────────────────────────────────────────────────────────

def batch_rearrangement_path_search(
    bam_path: str,
    ref_path: str,
    show_pileup: bool = False,
    min_clip_len: int = 50,
    min_support: int = 5,
    min_support_frac: float = 0.05,
    cluster_dist: int = 10,
    min_mapq: int = 60,
    pair_strategy: str = "split_read_bridge",
    max_dist=None,
    top_n: int = 50,
    save_dir: str = None,
    sa_share_threshold: float = 0.50,
    landing_tol: int = 200,
    fetch_padding: int = 500_000,
    max_aln_span: int = 1_000_000,
    minimap2_preset: str = "map-hifi",
    min_identity: float = 0.90,
    match_score: int = 2,
    mismatch_penalty: int = 8,
    gap_open: int = 4,
    gap_extend: int = 2,
    gap_open2: int = 24,
    gap_extend2: int = 1,
    pdv_min_split_vaf: float = 0.20,
    pdv_min_enrichment: float = 3.0,
    background_hom_vaf: float = 0.95,
    min_graph_segment_support: int = 10,
    # Legacy params kept for backwards compatibility
    kmer_size: int = 13,
    min_segment_len: int = 10,
    junction_tol: int = 50,
    ref_bin_tol: int = 100,
    min_block_support: int = 2,
    extra_windows=None,
):
    """Full batch clip-first pipeline.

    Stage order
    -----------
    1.  Discover soft-clip positions genome-wide.
    2.  Filter by depth fraction and N-adjacency.
    3.  Pair breakpoints via bridging-read count.
    4.  Smart-cluster pairs (merge if >sa_share_threshold SA reads shared).
    5.  Per cluster:
        a. compute_reference_window
        b. get_background_homozygous_sites (all reads, VAF≥0.95)
        c. find_split_read_specific_variants (PDV positions after masking)
        d. collect_breakpoint_segments (clip tail → last PDV)
        e. align_clips_to_ref (minimap2)
        f. call_pdvs (annotate PDV status per alignment)
        g. cluster_by_landing_site (group by ref locus)
        h. build_path_from_clips (one node per locus)
        i. plot_rearrangement_path (boxes + PDV-annotated edges)

    Returns
    -------
    pd.DataFrame (one row per cluster)
    """
    bam_stem = os.path.splitext(os.path.basename(bam_path))[0]
    if save_dir is None:
        save_dir = os.path.join("results", bam_stem)
    os.makedirs(save_dir, exist_ok=True)
    print(f"Output directory: {save_dir}")

    # ── Stage 1: clip discovery ───────────────────────────────────────────────
    print("\n[Stage 1] Discovering soft-clip positions ...")
    clip_df = find_major_clip_positions(
        bam_path,
        min_clip_len      = min_clip_len,
        min_support       = min_support,
        min_support_frac  = min_support_frac,
        cluster_dist      = cluster_dist,
        top_n             = top_n,
        min_mapq          = min_mapq,
        save_dir          = save_dir,
    )
    if clip_df.empty:
        print("No clip positions found. Exiting.")
        return pd.DataFrame()

    # ── Stage 2: N-adjacency filter ──────────────────────────────────────────
    print("\n[Stage 2] N-adjacency filter ...")
    fasta = pysam.FastaFile(ref_path)
    valid = []
    for _, row in clip_df.iterrows():
        if _has_adjacent_n(fasta, row["chrom"], int(row["position"])):
            print(f"  Dropping {row['chrom']}:{row['position']} (N-adjacent)")
        else:
            valid.append(row)
    fasta.close()
    clip_df = pd.DataFrame(valid).reset_index(drop=True) if valid else pd.DataFrame()
    if clip_df.empty:
        print("All positions dropped by N-adjacency filter.")
        return pd.DataFrame()

    # ── Stage 3: pair breakpoints ─────────────────────────────────────────────
    print("\n[Stage 3] Pairing breakpoints ...")
    pairs = pair_clip_positions(
        clip_df,
        bam_path  = bam_path,
        max_dist  = max_dist,
        strategy  = pair_strategy,
    )
    if not pairs:
        print("No breakpoint pairs found.")
        return pd.DataFrame()
    print(f"  {len(pairs)} pair(s) found.")

    # ── Stage 4: cluster pairs ────────────────────────────────────────────────
    print("\n[Stage 4] Clustering pairs ...")
    clusters = cluster_softclip_pairs(
        pairs, bam_path,
        share_threshold = sa_share_threshold,
        bridge_tol      = 200,
        min_mapq        = min_mapq,
    )
    if not clusters:
        print("No clusters formed.")
        return pd.DataFrame()

    # ── Stage 5: per-cluster analysis ────────────────────────────────────────
    results = []
    for idx, cluster in enumerate(clusters):
        label = f"cluster{idx+1:03d}_{cluster['chrom']}_{cluster['left_pos']}"
        print(f"\n{'='*60}")
        print(f"[Stage 5] Analysing {label}")
        try:
            res = analyze_sv_cluster(
                bam_path         = bam_path,
                ref_path         = ref_path,
                cluster          = cluster,
                show_pileup      = show_pileup,
                min_mapq         = min_mapq,
                pair_label       = label,
                save_dir         = save_dir,
                min_clip_len     = min_clip_len,
                landing_tol      = landing_tol,
                fetch_padding    = fetch_padding,
                max_aln_span     = max_aln_span,
                minimap2_preset  = minimap2_preset,
                min_identity     = min_identity,
                match_score      = match_score,
                mismatch_penalty = mismatch_penalty,
                gap_open         = gap_open,
                gap_extend       = gap_extend,
                gap_open2        = gap_open2,
                gap_extend2      = gap_extend2,
                pdv_min_split_vaf   = pdv_min_split_vaf,
                pdv_min_enrichment  = pdv_min_enrichment,
                background_hom_vaf  = background_hom_vaf,
                min_graph_segment_support = min_graph_segment_support,
            )
            if res:
                results.append(res)
        except Exception as exc:
            print(f"  [ERROR] {label}: {exc}")
            import traceback
            traceback.print_exc()

    if not results:
        print("\nNo results produced.")
        return pd.DataFrame()

    out_df  = pd.DataFrame(results)
    csv_out = os.path.join(save_dir, "batch_rearrangement_results.csv")
    out_df.to_csv(csv_out, index=False)
    print(f"\nDone — {len(results)} cluster(s). Results → {csv_out}")
    return out_df
