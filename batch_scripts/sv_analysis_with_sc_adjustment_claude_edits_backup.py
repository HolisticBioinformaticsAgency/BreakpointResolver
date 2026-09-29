import re
import pysam
pysam.set_verbosity(0)
import pandas as pd
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from IPython.display import display, Image, HTML
from collections import defaultdict
from statistics import median
import os

CLIP_OPS = {4, 5}

# ─────────────────────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────────────────────

def get_ref_chroms(ref_path):
    """Return the set of contig names present in the FASTA reference."""
    with pysam.FastaFile(ref_path) as fa:
        return set(fa.references)



def _safe_fetch(fasta, chrom, start, end):
    """Boundary-safe wrapper around pysam.FastaFile.fetch.

    Clamps start and end to [0, contig_length] so that fetches near the
    beginning or end of a contig (e.g. chrM position 1) never raise
    ValueError: start out of range.  Returns an empty string when the
    clamped window is zero-length.
    """
    try:
        contig_len = fasta.get_reference_length(chrom)
    except KeyError:
        return ""
    start = max(0, int(start))
    end   = min(contig_len, int(end))
    if end <= start:
        return ""
    return fasta.fetch(chrom, start, end).upper()

def _has_adjacent_n(fasta, chrom, pos, window=5):
    """Return True if any base within *window* bp of *pos* in the reference is N."""
    if chrom not in fasta.references:
        return False
    region = _safe_fetch(fasta, chrom, max(0, pos - window), pos + window)
    return "N" in region


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence string."""
    return seq.translate(str.maketrans("ATCGatcg", "TAGCtagc"))[::-1]


# ─────────────────────────────────────────────────────────────────────────────
# 1. Extract split reads
# ─────────────────────────────────────────────────────────────────────────────

def extract_split_reads(bam_path, chrom, start_search, end_search, padding):
    """Extract split reads (reads with an SA supplementary alignment tag) from a BAM file within a padded window around [start_search, end_search] on *chrom*. Only reads whose SA partner maps near start_search or end_search are returned."""
    samfile = pysam.AlignmentFile(bam_path, "rb")
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
# 2. Genome-wide major soft-clip / breakpoint finder
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
    """Scan a BAM file genome-wide (or a single *chrom*) and return a DataFrame of high-confidence soft-clip positions, 
    each representing a candidate structural variant breakpoint. Reads with low mapping quality, secondary/supplementary alignments, 
    or SA tags below *min_mapq* are excluded. Clip positions are clustered within *cluster_dist* bp and filtered by *min_support* / *min_support_frac* of local depth. 
    Results are written to *save_dir*/major_clip_positions.csv."""
    REF_OPS = {"M", "D", "N", "=", "X"}
    CLIP_STR = {"S", "H"}

    def _region_depth(chrom, pos, window=50):
        """Count primary, non-secondary, non-supplementary reads overlapping a *window*-bp region around *pos* on *chrom*."""
        count = 0
        with pysam.AlignmentFile(bam_path, "rb") as sf:
            for read in sf.fetch(chrom, max(0, pos - window), pos + window):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                count += 1
        return count

    left_pool = defaultdict(list)
    right_pool = defaultdict(list)
    samfile = pysam.AlignmentFile(bam_path, "rb")
    fetch_iter = samfile.fetch(chrom) if chrom else samfile.fetch()

    def _sa_majority_mapq(read, min_mapq):
        """Return True if the read has no SA tag or at least one SA entry has mapping quality >= *min_mapq*."""
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
        cigars = read.cigartuples
        ref_name = read.reference_name
        first_op, first_len = cigars[0]
        last_op, last_len = cigars[-1]

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
                sa_chr = fields[0]
                sa_pos = int(fields[1])
                sa_cigar = fields[3]
                sa_mapq = int(fields[4])
                if sa_mapq < min_mapq:
                    continue
                if sa_chr != ref_name:
                    continue
                ops = re.findall(r'(\d+)([MIDNSHP=X])', sa_cigar)
                if not ops:
                    continue
                ref_len = sum(int(l) for l, op in ops if op in REF_OPS)
                sa_end = sa_pos + ref_len
                sa_left_len = int(ops[0][0]) if ops[0][1] in CLIP_STR else 0
                key_sa_left = (read.query_name, sa_pos, "left")
                if ops[0][1] in CLIP_STR and sa_left_len >= min_clip_len and key_sa_left not in seen:
                    seen.add(key_sa_left)
                    left_pool[ref_name].append((sa_pos, sa_left_len))
                sa_right_len = int(ops[-1][0]) if ops[-1][1] in CLIP_STR else 0
                key_sa_right = (read.query_name, sa_end, "right")
                if ops[-1][1] in CLIP_STR and sa_right_len >= min_clip_len and key_sa_right not in seen:
                    seen.add(key_sa_right)
                    right_pool[ref_name].append((sa_end, sa_right_len))
    samfile.close()

    def _cluster_and_aggregate(pool, side):
        """Cluster a list of (position, clip_length) tuples within *cluster_dist* bp of each other, then aggregate each cluster into a summary row containing support count, mean/max clip length, and modal position."""
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
            support = len(cluster)
            all_lengths = [cl for _, cl in cluster]
            mode_pos = Counter(p for p, _ in cluster).most_common(1)[0][0]
            if support >= min_support:
                rows.append({
                    "side": side,
                    "position": mode_pos,
                    "support": support,
                    "mean_clip_len": round(sum(all_lengths) / len(all_lengths), 1),
                    "max_clip_len": max(all_lengths),
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
        keep_mask = []
        for _, row in df.iterrows():
            key = (row["chrom"], row["position"])
            if key not in depth_cache:
                depth_cache[key] = _region_depth(row["chrom"], row["position"])
            total_depth = depth_cache[key]
            frac_thresh = min_support_frac * total_depth
            effective_min = max(min_support, frac_thresh)
            passes = row["support"] >= effective_min
            keep_mask.append(passes)
            if not passes:
                print(
                    f" [depth-filter] {row['chrom']}:{row['position']} "
                    f"support={row['support']} depth={total_depth} "
                    f"required≥{effective_min:.1f} ({min_support_frac*100:.0f}% of depth) "
                    f"→ DROPPED"
                )
        df = df[keep_mask].reset_index(drop=True)
        if df.empty:
            print("No clip positions passed the fractional-support filter.")
            return df

    if top_n:
        df = df.head(top_n)

    _csv_out = os.path.join(save_dir, "major_clip_positions.csv") if save_dir else "major_clip_positions.csv"
    os.makedirs(save_dir, exist_ok=True) if save_dir else None
    df.to_csv(_csv_out, index=False)
    print(f"Saved: {_csv_out}")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 3. Pair breakpoint positions
# ─────────────────────────────────────────────────────────────────────────────

def count_bridging_reads(bam_path, chrom, pos_a, pos_b, tol=200, min_mapq=0):
    """Count distinct read names that bridge both pos_a and pos_b on *chrom*: the primary alignment must clip near one position and the SA tag must map near the other, within *tol* bp."""
    REF_OPS = {"M", "D", "N", "=", "X"}
    bridging = set()

    def _clips_near(cigar_ops, ref_start, target):
        if not cigar_ops:
            return False
        ref_end = ref_start + sum(l for op, l in cigar_ops if op in (0, 2, 3, 7, 8))
        right_clip = cigar_ops[-1][0] in (4, 5) and abs(ref_end - target) <= tol
        left_clip = cigar_ops[0][0] in (4, 5) and abs(ref_start - target) <= tol
        return right_clip or left_clip

    def _sa_near(read, target):
        if not read.has_tag("SA"):
            return False
        for entry in read.get_tag("SA").rstrip(";").split(";"):
            fields = entry.split(",")
            if len(fields) < 4 or fields[0] != chrom:
                continue
            sa_pos = int(fields[1])
            sa_cigar = fields[3]
            ops = re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
            if not ops:
                continue
            ref_len = sum(int(l) for l, op in ops if op in REF_OPS)
            sa_end = sa_pos + ref_len
            if ops[0][1] in ("S", "H") and abs(sa_pos - target) <= tol:
                return True
            if ops[-1][1] in ("S", "H") and abs(sa_end - target) <= tol:
                return True
        return False

    with pysam.AlignmentFile(bam_path, "rb") as sf:
        fetch_start = max(0, min(pos_a, pos_b) - tol)
        fetch_end = max(pos_a, pos_b) + tol
        for read in sf.fetch(chrom, fetch_start, fetch_end):
            if read.is_unmapped or not read.cigartuples or read.is_secondary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            cigars = read.cigartuples
            ref_start = read.reference_start
            if _clips_near(cigars, ref_start, pos_a) and _sa_near(read, pos_b):
                bridging.add(read.query_name)
            elif _clips_near(cigars, ref_start, pos_b) and _sa_near(read, pos_a):
                bridging.add(read.query_name)
    return len(bridging)


def _reads_supporting_pair(bam_path, chrom, pos_a, pos_b, tol=200, min_mapq=0):
    """Return the *set* of read names bridging pos_a ↔ pos_b (needed for Union-Find)."""
    REF_OPS = {"M", "D", "N", "=", "X"}
    bridging = set()

    def _clips_near(cigar_ops, ref_start, target):
        if not cigar_ops:
            return False
        ref_end = ref_start + sum(l for op, l in cigar_ops if op in (0, 2, 3, 7, 8))
        return (
            (cigar_ops[-1][0] in (4, 5) and abs(ref_end - target) <= tol)
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
            sa_cigar = fields[3]
            ops = re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
            if not ops:
                continue
            ref_len = sum(int(l) for l, op in ops if op in REF_OPS)
            sa_end = sa_pos + ref_len
            if ops[0][1] in ("S", "H") and abs(sa_pos - target) <= tol:
                return True
            if ops[-1][1] in ("S", "H") and abs(sa_end - target) <= tol:
                return True
        return False

    with pysam.AlignmentFile(bam_path, "rb") as sf:
        fetch_start = max(0, min(pos_a, pos_b) - tol)
        fetch_end = max(pos_a, pos_b) + tol
        for read in sf.fetch(chrom, fetch_start, fetch_end):
            if read.is_unmapped or not read.cigartuples or read.is_secondary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            cigars = read.cigartuples
            ref_start = read.reference_start
            if _clips_near(cigars, ref_start, pos_a) and _sa_near(read, pos_b):
                bridging.add(read.query_name)
            elif _clips_near(cigars, ref_start, pos_b) and _sa_near(read, pos_a):
                bridging.add(read.query_name)
    return bridging


def pair_clip_positions(clip_df, bam_path=None, max_dist=None, strategy="split_read_bridge", bridge_tol=200):
    """Build (chrom, left_pos, right_pos) pairs from find_major_clip_positions() output."""
    pairs = []
    if strategy == "split_read_bridge":
        if bam_path is None:
            print(" [pair_clip_positions] bam_path not supplied; falling back to greedy_rank.")
            return pair_clip_positions(clip_df, bam_path=None, max_dist=max_dist, strategy="greedy_rank")
        idxs = list(clip_df.index)
        n = len(idxs)
        candidates = []
        print(f" Scoring {n*(n-1)//2} candidate pairs by bridging read count ...")
        for ii in range(n):
            for jj in range(ii + 1, n):
                i, j = idxs[ii], idxs[jj]
                row_a = clip_df.loc[i]
                row_b = clip_df.loc[j]
                if row_a["chrom"] != row_b["chrom"]:
                    continue
                dist = abs(int(row_b["position"]) - int(row_a["position"]))
                if max_dist is not None and dist > max_dist:
                    continue
                if dist == 0:
                    continue
                bc = count_bridging_reads(
                    bam_path, row_a["chrom"],
                    int(row_a["position"]), int(row_b["position"]),
                    tol=bridge_tol,
                )
                candidates.append((bc, i, j))
                print(f" {row_a['chrom']}:{int(row_a['position']):,} ↔ "
                      f"{int(row_b['position']):,} bridge_reads={bc}")
        candidates.sort(key=lambda x: -x[0])
        used = set()
        for bc, i, j in candidates:
            if i in used or j in used:
                continue
            if bc == 0:
                continue
            used.add(i); used.add(j)
            row_a = clip_df.loc[i]
            row_b = clip_df.loc[j]
            left = int(min(row_a["position"], row_b["position"]))
            right = int(max(row_a["position"], row_b["position"]))
            pairs.append((row_a["chrom"], left, right))
        if not pairs:
            print(" [split_read_bridge] No pairs with bridging reads found; falling back to greedy_rank.")
            return pair_clip_positions(clip_df, bam_path=None, max_dist=max_dist, strategy="greedy_rank")

    elif strategy == "greedy_rank":
        used = set()
        for i, row_a in clip_df.iterrows():
            if i in used:
                continue
            best_j = None
            for j, row_b in clip_df.iterrows():
                if j <= i or j in used:
                    continue
                if row_b["chrom"] != row_a["chrom"]:
                    continue
                dist = abs(int(row_b["position"]) - int(row_a["position"]))
                if max_dist is not None and dist > max_dist:
                    continue
                best_j = j
                break
            if best_j is None:
                continue
            used.add(i); used.add(best_j)
            row_b = clip_df.loc[best_j]
            left = int(min(row_a["position"], row_b["position"]))
            right = int(max(row_a["position"], row_b["position"]))
            pairs.append((row_a["chrom"], left, right))

    elif strategy == "consecutive":
        for i in range(0, len(clip_df) - 1, 2):
            a, b = clip_df.iloc[i], clip_df.iloc[i + 1]
            if a["chrom"] != b["chrom"]:
                continue
            left = int(min(a["position"], b["position"]))
            right = int(max(a["position"], b["position"]))
            if max_dist is None or (right - left) <= max_dist:
                pairs.append((a["chrom"], left, right))

    elif strategy == "proximity":
        right_clips = clip_df[clip_df["side"] == "right"].copy()
        left_clips = clip_df[clip_df["side"] == "left"].copy()
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
# 4. Cluster soft-clip pairs (Union-Find on shared SA reads)
# ─────────────────────────────────────────────────────────────────────────────

def cluster_softclip_pairs(
    pairs,
    bam_path,
    share_threshold=0.50,
    bridge_tol=200,
    min_mapq=0,
):
    """
    Merge breakpoint pairs that share > share_threshold of their SA-supporting
    reads (Jaccard over the smaller set).  Uses Union-Find.

    For each merged cluster the coordinate span is taken from the pair with the
    WIDEST span (min left, max right across all members).  Narrower pairs have
    their soft-clip tails effectively re-anchored to this wider span when reads
    are later extracted and consensus is built.

    Returns
    -------
    list of dicts, each with keys:
        pairs        – original (chrom, lp, rp) tuples that were merged
        reads        – union of bridging read-name sets
        chrom        – chromosome
        refined_left – widest left coordinate (min of all member lefts)
        refined_right– widest right coordinate (max of all member rights)
    """
    if not pairs:
        return []

    # Collect per-pair bridging read sets
    pair_reads = []
    for chrom, lp, rp in pairs:
        rset = _reads_supporting_pair(bam_path, chrom, lp, rp, tol=bridge_tol, min_mapq=min_mapq)
        pair_reads.append(rset)
        print(f" Pair {chrom}:{lp:,}–{rp:,}  {len(rset)} SA reads")

    n = len(pairs)
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
            si, sj = pair_reads[i], pair_reads[j]
            if not si or not sj:
                continue
            shared = si & sj
            min_size = min(len(si), len(sj))
            frac = len(shared) / min_size if min_size > 0 else 0.0
            if frac > share_threshold:
                print(
                    f" Merging pair {i} and {j}: "
                    f"{len(shared)}/{min_size} shared reads ({frac:.1%})"
                )
                _union(i, j)

    cluster_map = defaultdict(list)
    for i in range(n):
        cluster_map[_find(i)].append(i)

    clusters = []
    for root, members in cluster_map.items():
        chrom = pairs[members[0]][0]
        all_left  = [pairs[m][1] for m in members]
        all_right = [pairs[m][2] for m in members]
        merged_reads = set()
        for m in members:
            merged_reads |= pair_reads[m]

        # Widest span = outer envelope: min of all lefts, max of all rights.
        # Using the single pair with the largest (right-left) is WRONG because
        # a pair with a lower left anchor can win the span sort while another
        # pair reaches further right.  Direct min/max is always correct.
        refined_left  = min(all_left)
        refined_right = max(all_right)

        clusters.append({
            "pairs":          [pairs[m] for m in members],
            "reads":          merged_reads,
            "chrom":          chrom,
            "refined_left":   refined_left,
            "refined_right":  refined_right,
        })

    print(f"\n{len(clusters)} cluster(s) from {n} pair(s).")
    return clusters


# ─────────────────────────────────────────────────────────────────────────────
# 5. Build reference-anchored consensus
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# 3b. Single-pass clip-pair discovery  (replaces Stages 1-4)
# ─────────────────────────────────────────────────────────────────────────────

def find_major_clip_pairs(
    bam_path,
    ref_path,
    chrom=None,
    min_clip_len=50,
    min_support=5,
    min_support_frac=0.05,
    cluster_dist=10,
    top_n=50,
    min_mapq=60,
    max_dist=None,
    save_dir=None,
):
    """Single-pass discovery of soft-clip *pairs* directly from SA reads.

    For every read that carries an SA tag, the primary alignment and the
    SA entry with the highest MAPQ are selected (top-2 by MAPQ).  Each
    alignment must have a soft/hard clip of >= *min_clip_len* bases.  The
    clip position (left-end or right-end of the alignment, depending on
    which end is clipped) is recorded per read.

    After the full BAM pass, clip positions on each chromosome are binned
    using the same sliding-window _cluster_and_aggregate logic as
    find_major_clip_positions, producing a modal position per bin.  Reads
    are then counted per (chrom, modal_left, modal_right) pair key.

    Filters applied (in order):
        1. min_clip_len       - both clips must be >= this length
        2. min_mapq           - primary mapping_quality >= min_mapq
        3. same-chrom only    - inter-chromosomal pairs discarded
        4. max_dist           - pairs spanning > max_dist bp discarded
        5. min_support        - absolute read-count floor  (applied FIRST)
        6. top_n              - keep highest-support pairs (applied SECOND)
        7. min_support_frac   - pair support >= frac x min(depth_left,
                                depth_right); only run on top_n survivors

    Saves major_clip_pairs.csv to *save_dir* and returns a DataFrame with
    columns:
        chrom, pos_left, pos_right, support,
        mean_clip_len_left, mean_clip_len_right,
        mapq_left, mapq_right
    """
    import os, bisect, sys, time
    from collections import defaultdict, Counter

    # ── tqdm: use notebook version inside Jupyter, plain otherwise ───────────
    try:
        from tqdm.auto import tqdm as _tqdm
    except ImportError:
        _tqdm = None

    def _stage(msg):
        ts = time.strftime("%H:%M:%S")
        print(f"\n[{ts}] {msg}", flush=True)

    def _clip_pos_and_len(read):
        cigars = read.cigartuples
        if not cigars:
            return None
        llen = cigars[0][1]  if cigars[0][0]  in CLIP_OPS else 0
        rlen = cigars[-1][1] if cigars[-1][0] in CLIP_OPS else 0
        if llen == 0 and rlen == 0:
            return None
        if rlen >= llen:
            return (read.reference_end, rlen)
        else:
            return (read.reference_start, llen)

    def _sa_clip_pos_and_len(sa_entry_str):
        fields = sa_entry_str.split(",")
        if len(fields) < 5:
            return None
        sa_chrom = fields[0]
        sa_pos   = int(fields[1])
        sa_cigar = fields[3]
        try:
            sa_mapq = int(fields[4])
        except ValueError:
            sa_mapq = 0
        ops = re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
        if not ops:
            return None
        REF_OPS = {"M", "D", "N", "=", "X"}
        ref_len = sum(int(l) for l, op in ops if op in REF_OPS)
        sa_end  = sa_pos + ref_len
        llen = int(ops[0][0])  if ops[0][1]  in ("S", "H") else 0
        rlen = int(ops[-1][0]) if ops[-1][1] in ("S", "H") else 0
        if llen == 0 and rlen == 0:
            return None
        if rlen >= llen:
            return (sa_chrom, sa_end, rlen, sa_mapq)
        else:
            return (sa_chrom, sa_pos, llen, sa_mapq)

    # ── Stage A: BAM read scan ────────────────────────────────────────────────
    _stage("Stage A  |  Scanning BAM for SA-tagged reads ...")
    print(f"          BAM : {bam_path}", flush=True)

    pos_lists  = defaultdict(list)
    raw_pairs  = []

    n_primary  = 0
    n_sa_pass  = 0
    n_clip_pass = 0
    _t0        = time.time()

    with pysam.AlignmentFile(bam_path, "rb") as sf:
        # Get total mapped reads for progress bar (from index stats)
        try:
            total_mapped = sum(s.mapped for s in sf.get_index_statistics())
        except Exception:
            total_mapped = None

        fetch_iter = sf.fetch(chrom) if chrom else sf.fetch()

        if _tqdm is not None:
            pbar = _tqdm(
                fetch_iter,
                total=total_mapped,
                desc="  Scanning reads",
                unit=" reads",
                unit_scale=True,
                dynamic_ncols=True,
                miniters=100_000,
            )
        else:
            pbar = fetch_iter

        for read in pbar:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            n_primary += 1

            if not read.has_tag("SA"):
                continue
            if read.mapping_quality < min_mapq:
                continue

            n_sa_pass += 1

            prim = _clip_pos_and_len(read)
            if prim is None:
                continue
            p_pos, p_clip_len = prim
            if p_clip_len < min_clip_len:
                continue

            p_chrom = read.reference_name
            p_mapq  = read.mapping_quality

            best_sa      = None
            best_sa_mapq = -1
            for entry in read.get_tag("SA").rstrip(";").split(";"):
                if not entry:
                    continue
                parsed = _sa_clip_pos_and_len(entry)
                if parsed is None:
                    continue
                sa_chrom, sa_pos, sa_clip_len, sa_mapq = parsed
                if sa_clip_len < min_clip_len:
                    continue
                if sa_mapq > best_sa_mapq:
                    best_sa_mapq = sa_mapq
                    best_sa = (sa_chrom, sa_pos, sa_clip_len, sa_mapq)

            if best_sa is None:
                continue

            sa_chrom, sa_pos, sa_clip_len, sa_mapq = best_sa
            n_clip_pass += 1

            pos_lists[p_chrom].append(p_pos)
            pos_lists[sa_chrom].append(sa_pos)
            raw_pairs.append((
                p_chrom,  p_pos,  p_clip_len,  p_mapq,
                sa_chrom, sa_pos, sa_clip_len, sa_mapq,
            ))

        if _tqdm is not None:
            pbar.close()

    elapsed = time.time() - _t0
    print(
        f"          Done in {elapsed:.1f}s  |  "
        f"primary reads={n_primary:,}  "
        f"SA+MAPQ pass={n_sa_pass:,}  "
        f"clip pairs={n_clip_pass:,}",
        flush=True,
    )

    if not raw_pairs:
        print("  No SA-bridging reads found.")
        return pd.DataFrame()

    # ── Stage B: Position binning ─────────────────────────────────────────────
    _stage("Stage B  |  Binning clip positions (cluster_dist={}) ...".format(cluster_dist))

    def _build_bins(positions_unsorted):
        if not positions_unsorted:
            return []
        positions_sorted = sorted(positions_unsorted)
        bins = []
        cluster = [positions_sorted[0]]
        for p in positions_sorted[1:]:
            if p - cluster[-1] <= cluster_dist:
                cluster.append(p)
            else:
                modal = Counter(cluster).most_common(1)[0][0]
                bins.append((cluster[0], cluster[-1], modal))
                cluster = [p]
        modal = Counter(cluster).most_common(1)[0][0]
        bins.append((cluster[0], cluster[-1], modal))
        return bins

    modal_bins = {}
    for c, positions in pos_lists.items():
        modal_bins[c] = _build_bins(positions)
        print(f"          {c:<16}  raw positions={len(positions):>8,}  bins={len(modal_bins[c]):>6,}", flush=True)

    def _snap(chrom, pos):
        bins = modal_bins.get(chrom)
        if not bins:
            return pos
        starts = [b[0] for b in bins]
        idx = bisect.bisect_right(starts, pos) - 1
        if idx >= 0:
            b_start, b_end, modal = bins[idx]
            if b_start <= pos <= b_end + cluster_dist:
                return modal
        return pos

    # ── Stage C: Tally pair support ───────────────────────────────────────────
    _stage("Stage C  |  Tallying pair support ...")

    pair_support   = Counter()
    pair_clip_lens = defaultdict(lambda: {
        "left": [], "right": [], "mapq_left": [], "mapq_right": []
    })

    for (p_chrom, p_pos, p_clip_len, p_mapq,
         sa_chrom, sa_pos, sa_clip_len, sa_mapq) in raw_pairs:

        if p_chrom != sa_chrom:
            continue

        mp  = _snap(p_chrom,  p_pos)
        msa = _snap(sa_chrom, sa_pos)

        if mp == msa:
            continue
        left_pos  = min(mp, msa)
        right_pos = max(mp, msa)

        if max_dist is not None and (right_pos - left_pos) > max_dist:
            continue

        key = (p_chrom, left_pos, right_pos)
        pair_support[key] += 1
        d = pair_clip_lens[key]
        if mp <= msa:
            d["left"].append(p_clip_len);  d["mapq_left"].append(p_mapq)
            d["right"].append(sa_clip_len); d["mapq_right"].append(sa_mapq)
        else:
            d["left"].append(sa_clip_len);  d["mapq_left"].append(sa_mapq)
            d["right"].append(p_clip_len);  d["mapq_right"].append(p_mapq)

    print(f"          {len(pair_support):,} unique pair(s) tallied before filters.", flush=True)

    if not pair_support:
        print("  No same-chrom pairs within max_dist found.")
        return pd.DataFrame()

    # ── Stage D: Filtering ────────────────────────────────────────────────────
    _stage("Stage D  |  Filtering pairs ...")

    rows = []
    for (chrom_k, left_pos, right_pos), support in pair_support.items():
        d = pair_clip_lens[(chrom_k, left_pos, right_pos)]
        rows.append({
            "chrom":               chrom_k,
            "pos_left":            left_pos,
            "pos_right":           right_pos,
            "support":             support,
            "mean_clip_len_left":  round(sum(d["left"])  / len(d["left"]),  1) if d["left"]  else 0,
            "mean_clip_len_right": round(sum(d["right"]) / len(d["right"]), 1) if d["right"] else 0,
            "mapq_left":           round(sum(d["mapq_left"])  / len(d["mapq_left"]),  1) if d["mapq_left"]  else 0,
            "mapq_right":          round(sum(d["mapq_right"]) / len(d["mapq_right"]), 1) if d["mapq_right"] else 0,
        })

    df = (
        pd.DataFrame(rows)
        .sort_values("support", ascending=False)
        .reset_index(drop=True)
    )

    n_before = len(df)
    if min_support:
        df = df[df["support"] >= min_support].reset_index(drop=True)
    print(f"          min_support >= {min_support}: {n_before:,} → {len(df):,} pairs", flush=True)

    n_before = len(df)
    if top_n:
        df = df.head(top_n)
    print(f"          top_n={top_n}:  {n_before:,} → {len(df):,} pairs", flush=True)

    if min_support_frac and not df.empty:
        def _region_depth(chrom, pos, window=50):
            count = 0
            with pysam.AlignmentFile(bam_path, "rb") as sf:
                for read in sf.fetch(chrom, max(0, pos - window), pos + window):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    count += 1
            return count

        print(
            f"          min_support_frac={min_support_frac}: "
            f"checking depth at {len(df)*2} breakpoint(s) ...",
            flush=True,
        )
        depth_cache = {}
        keep_mask   = []
        for _, row in df.iterrows():
            for pos_key in [(row["chrom"], row["pos_left"]), (row["chrom"], row["pos_right"])]:
                if pos_key not in depth_cache:
                    depth_cache[pos_key] = _region_depth(pos_key[0], pos_key[1])
            d_left        = depth_cache[(row["chrom"], row["pos_left"])]
            d_right       = depth_cache[(row["chrom"], row["pos_right"])]
            min_depth     = min(d_left, d_right)
            effective_min = max(min_support, min_support_frac * min_depth)
            passes        = row["support"] >= effective_min
            keep_mask.append(passes)
            if not passes:
                print(
                    f"          [DROP] {row['chrom']}:{row['pos_left']:,}–{row['pos_right']:,} "
                    f"support={row['support']} min_depth={min_depth} "
                    f"required>={effective_min:.1f}",
                    flush=True,
                )
        df = df[keep_mask].reset_index(drop=True)
        if df.empty:
            print("  No pairs passed the fractional-support filter.")
            return df

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        bam_name = os.path.splitext(os.path.basename(bam_path))[0]
        out = os.path.join(save_dir, f"major_clip_pairs_{bam_name}.csv")
        df.to_csv(out, index=False)
        print(f"          Saved: {out}", flush=True)

    _stage(f"find_major_clip_pairs complete  |  {len(df)} pair(s) retained.")
    return df



def cluster_clip_pairs_by_reads(
    pairs_df,
    bam_path,
    cluster_dist=10,
    min_mapq=0,
    bridge_tol=200,
):
    """Read-centric Union-Find clustering of clip pairs.

    For every read with an SA tag that was already counted in *pairs_df*,
    we know it links two positions.  Pairs whose positions fall within
    *cluster_dist* bp of each other are treated as the same node.

    Union-Find is run directly on binned positions:
        Read bridges (pos_a, pos_b)  → union(bin(pos_a), bin(pos_b))
        Read bridges (pos_a, pos_c)  → union(bin(pos_a), bin(pos_c))
        Result: {pos_a, pos_b, pos_c} one cluster,
                refined_left  = min(pos_a, pos_b, pos_c)
                refined_right = max(pos_a, pos_b, pos_c)

    No additional BAM fetch is needed — the pair positions from pairs_df
    are used directly.  Each row in pairs_df already represents a
    pre-tallied (chrom, pos_left, pos_right) with support count.

    Returns the same cluster list format as cluster_softclip_pairs:
        [{ chrom, refined_left, refined_right, pairs, reads }]
    """
    if pairs_df.empty:
        return []

    # Collect all unique (chrom, pos) nodes
    # pos_key → canonical bin representative (itself initially)
    parent = {}  # pos_key → pos_key

    def _find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(x, y):
        px, py = _find(x), _find(y)
        if px != py:
            parent[px] = py

    # Initialise Union-Find nodes from pairs_df
    for _, row in pairs_df.iterrows():
        key_l = (row["chrom"], int(row["pos_left"]))
        key_r = (row["chrom"], int(row["pos_right"]))
        if key_l not in parent:
            parent[key_l] = key_l
        if key_r not in parent:
            parent[key_r] = key_r
        _union(key_l, key_r)

    # Optional: merge nodes within cluster_dist on the same chrom
    # (handles the rare case where two pair rows share nearly-identical breakpoints
    #  that were binned to slightly different modal positions)
    nodes_by_chrom = defaultdict(list)
    for key in parent:
        nodes_by_chrom[key[0]].append(key[1])
    for chrom_k, positions in nodes_by_chrom.items():
        sorted_pos = sorted(set(positions))
        for i in range(len(sorted_pos) - 1):
            if sorted_pos[i+1] - sorted_pos[i] <= cluster_dist:
                _union((chrom_k, sorted_pos[i]), (chrom_k, sorted_pos[i+1]))

    # Group by root
    cluster_map = defaultdict(list)
    for key in parent:
        cluster_map[_find(key)].append(key)

    # Build cluster dicts
    pairs_list = [
        (row["chrom"], int(row["pos_left"]), int(row["pos_right"]), int(row["support"]))
        for _, row in pairs_df.iterrows()
    ]

    clusters = []
    for root, members in cluster_map.items():
        chrom_k = members[0][0]
        positions = [m[1] for m in members]
        refined_left  = min(positions)
        refined_right = max(positions)

        # Collect original pair tuples that belong to this cluster
        member_pos_set = set(positions)
        member_pairs = [
            (c, lp, rp)
            for c, lp, rp, _ in pairs_list
            if c == chrom_k and lp in member_pos_set and rp in member_pos_set
        ]
        # Total support = sum of support across member pairs
        total_support = sum(
            s for c, lp, rp, s in pairs_list
            if c == chrom_k and lp in member_pos_set and rp in member_pos_set
        )

        clusters.append({
            "chrom":          chrom_k,
            "refined_left":   refined_left,
            "refined_right":  refined_right,
            "pairs":          member_pairs,
            "reads":          set(),   # populated below if needed
            "support":        total_support,
        })

    # Sort clusters by total support descending
    clusters.sort(key=lambda c: -c["support"])

    print(f"cluster_clip_pairs_by_reads: {len(clusters)} cluster(s) from {len(pairs_df)} pair(s).")
    for i, cl in enumerate(clusters, 1):
        print(
            f"  {i:>3}. {cl['chrom']} "
            f"{cl['refined_left']:,}–{cl['refined_right']:,} "
            f"support={cl['support']} ({len(cl['pairs'])} pair(s))"
        )
    return clusters


def get_homozygous_variants(bam_path, chrom, win_start, win_end, min_af=0.95, min_depth=5):
    """Return a dict {genomic_pos: base} for positions where a non-ref base has
    allele frequency >= min_af across ALL reads."""
    samfile = pysam.AlignmentFile(bam_path, "rb")
    position_counts = {}
    for read in samfile.fetch(chrom, win_start, win_end):
        if read.is_unmapped or read.query_sequence is None:
            continue
        ref_cursor = read.reference_start
        query_cursor = 0
        if not read.cigartuples:
            continue
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                for d in range(length):
                    rp = ref_cursor + d
                    if win_start <= rp < win_end:
                        base = read.query_sequence[query_cursor + d].upper()
                        if rp not in position_counts:
                            position_counts[rp] = Counter()
                        position_counts[rp][base] += 1
                ref_cursor += length
                query_cursor += length
            elif op in (2, 3):
                ref_cursor += length
            elif op == 1:
                query_cursor += length
            elif op in (4, 5):
                if op == 4:
                    query_cursor += length
    samfile.close()
    homozygous = {}
    for pos, counts in position_counts.items():
        total = sum(counts.values())
        if total < min_depth:
            continue
        top_base, top_count = counts.most_common(1)[0]
        if top_count / total >= min_af:
            homozygous[pos] = top_base
    return homozygous


def build_split_consensus(
    reads, refined_left, refined_right, flank, fasta, chrom, sv_type,
    bam_path, show_pileup=True, show_homozygous_variants=False,
    pileup_save_dir=None, pair_label=""
):
    """Build a reference-anchored consensus sequence for the left and right flanks of an SV breakpoint from a collection of split reads. For each flank, read columns are piled up over a *flank*-bp window, background homozygous variants are detected and masked, and a per-column majority-vote consensus is produced. Optionally renders and saves HTML / PNG pileup visualisations."""
    if chrom not in fasta.references:
        print(f"WARNING: {chrom} not in reference — cannot build consensus, skipping.")
        return None, None, None, None, refined_right, "", "", {}

    is_deletion = sv_type == "deletion"
    right_bp = refined_right - 1

    if is_deletion:
        ref_left  = _safe_fetch(fasta, chrom, refined_left - flank, refined_left)
        ref_right = _safe_fetch(fasta, chrom, right_bp, right_bp + flank)
    else:
        ref_left  = _safe_fetch(fasta, chrom, refined_left, refined_left + flank)
        ref_right = _safe_fetch(fasta, chrom, refined_right - flank, refined_right)

    def _cigar_walk_to_columns(read, win_start, win_end):
        """Walk a read's CIGAR string and return a dict mapping column index (relative to win_start) to the called base for every reference-consuming operation that falls within [win_start, win_end)."""
        bases = {}
        seq = read.query_sequence
        if seq is None:
            return bases
        ref_cursor = read.reference_start
        query_cursor = 0
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                for d in range(length):
                    rp = ref_cursor + d
                    if win_start <= rp < win_end:
                        bases[rp - win_start] = seq[query_cursor + d].upper()
                ref_cursor += length
                query_cursor += length
            elif op in (2, 3):
                for d in range(length):
                    rp = ref_cursor + d
                    if win_start <= rp < win_end:
                        bases[rp - win_start] = "-"
                ref_cursor += length
            elif op == 1:
                query_cursor += length
            elif op == 4:
                query_cursor += length
            elif op == 5:
                pass
        return bases

    def _build_consensus_from_columns(read_column_list, ref_anchor, flank, min_col_depth=3):
        """Build a consensus string of length *flank* from a list of column dicts. Columns with depth below *min_col_depth* fall back to the *ref_anchor* base."""
        columns = [Counter() for _ in range(flank)]
        for col_dict in read_column_list:
            for col, base in col_dict.items():
                if 0 <= col < flank and base != "N":
                    columns[col][base] += 1
        result = []
        low_cov_pos = []
        for i, c in enumerate(columns):
            depth = sum(c.values())
            if depth < min_col_depth:
                result.append(ref_anchor[i])
                if depth > 0:
                    low_cov_pos.append(i)
                continue
            top_base, top_count = c.most_common(1)[0]
            total = depth
            if top_base == "-" and top_count / total > 0.5:
                result.append("-")
            else:
                no_del = Counter({b: n for b, n in c.items() if b != "-"})
                if no_del:
                    result.append(no_del.most_common(1)[0][0])
                else:
                    result.append(ref_anchor[i])
        return "".join(result), low_cov_pos

    def _columns_to_padded_seq(col_dict, flank):
        """Convert a column dict to a fixed-length string of length *flank*, filling missing columns with 'N'."""
        return "".join(col_dict.get(i, "N") for i in range(flank))

    def _read_label(read):
        """Return a human-readable label string for a read, including its name, CIGAR string, and SA tag."""
        cigar_str = read.cigarstring or "none"
        sa_str = read.get_tag("SA") if read.has_tag("SA") else "none"
        return f"{read.query_name} CIGAR:{cigar_str} SA:{sa_str}"

    def _build_sa_index(bam_path, read_names, chrom):
        """Build two lookup dicts for a set of read names on *chrom*: supp_index maps (name, ref_start+1) → supplementary alignment; primary_index maps name → primary alignment."""
        supp_index = {}
        primary_index = {}
        with pysam.AlignmentFile(bam_path, "rb") as sf:
            for read in sf.fetch(chrom):
                if read.query_name not in read_names:
                    continue
                if read.is_unmapped or not read.cigartuples:
                    continue
                if read.is_supplementary:
                    supp_index[(read.query_name, read.reference_start + 1)] = read
                elif not read.is_secondary:
                    primary_index[read.query_name] = read
        return supp_index, primary_index

    def _cigar_walk_to_columns_from_aln(aln, win_start, win_end):
        """Thin wrapper around _cigar_walk_to_columns for a pre-fetched alignment object."""
        return _cigar_walk_to_columns(aln, win_start, win_end)

    def _cigar_walk_to_columns_sa_primary(primary_read, sa_pos, sa_cigar, sa_strand, win_start, win_end):
        """Extract column bases for an SA supplementary alignment using the query sequence from the primary read. Handles strand orientation via reverse complement and correctly slices the query sequence to the SA portion."""
        seq = primary_read.query_sequence
        if seq is None or not primary_read.cigartuples:
            return {}
        primary_consumes = sum(
            length for op, length in primary_read.cigartuples
            if op in (0, 1, 4, 7, 8)
        )
        primary_cigars = primary_read.cigartuples
        if primary_cigars[-1][0] in (4, 5):
            non_clip_consumes = sum(
                length for op, length in primary_cigars if op in (0, 1, 7, 8)
            )
            sa_query_start = non_clip_consumes
        else:
            sa_query_start = 0
        sa_ops = re.findall(r'(\d+)([MIDNSHP=X])', sa_cigar)
        if not sa_ops:
            return {}
        sa_query_len = sum(int(l) for l, op in sa_ops if op in ('M', 'I', 'S', '=', 'X'))
        sa_query_end = sa_query_start + sa_query_len
        sa_seq_raw = seq[sa_query_start:sa_query_end]
        if not sa_seq_raw:
            return {}
        if sa_strand == '-':
            sa_seq_raw = reverse_complement(sa_seq_raw)
        bases = {}
        ref_cursor = sa_pos
        q_cursor = 0
        op_start = 0
        if sa_ops[0][1] in ('S', 'H'):
            if sa_ops[0][1] == 'S':
                q_cursor += int(sa_ops[0][0])
            op_start = 1
        for l_str, op in sa_ops[op_start:]:
            length = int(l_str)
            if op in ('M', '=', 'X'):
                for d in range(length):
                    rp = ref_cursor + d
                    if win_start <= rp < win_end:
                        qi = q_cursor + d
                        if qi < len(sa_seq_raw):
                            bases[rp - win_start] = sa_seq_raw[qi].upper()
                ref_cursor += length
                q_cursor += length
            elif op in ('D', 'N'):
                for d in range(length):
                    rp = ref_cursor + d
                    if win_start <= rp < win_end:
                        bases[rp - win_start] = '-'
                ref_cursor += length
            elif op == 'I':
                q_cursor += length
            elif op == 'S':
                q_cursor += length
            elif op == 'H':
                pass
        return bases

    def _sa_clips_near(read, target_pos, clip_side, tol=50):
        """Return a list of (sa_pos, sa_cigar, sa_strand) tuples for all SA entries of *read* whose clip end (left or right, as specified by *clip_side*) falls within *tol* bp of *target_pos*."""
        if not read.has_tag("SA"):
            return []
        REF_OPS_SA = {'M', 'D', 'N', '=', 'X'}
        hits = []
        for sa_entry in read.get_tag("SA").rstrip(";").split(";"):
            fields = sa_entry.split(",")
            if len(fields) < 4 or fields[0] != read.reference_name:
                continue
            sa_pos = int(fields[1])
            sa_strand = fields[2] if len(fields) > 2 else '+'
            sa_cigar = fields[3]
            ops = re.findall(r'(\d+)([MIDNSHP=X])', sa_cigar)
            if not ops:
                continue
            ref_len = sum(int(l) for l, op in ops if op in REF_OPS_SA)
            sa_end = sa_pos + ref_len
            if clip_side == 'right' and ops[-1][1] in ('S', 'H'):
                if abs(sa_end - target_pos) < tol:
                    hits.append((sa_pos, sa_cigar, sa_strand))
            if clip_side == 'left' and ops[0][1] in ('S', 'H'):
                if abs(sa_pos - target_pos) < tol:
                    hits.append((sa_pos, sa_cigar, sa_strand))
        return hits

    def _primary_clips_near(read, target_pos, clip_side, tol=50):
        """Return True if the primary alignment has a soft/hard clip on *clip_side* within *tol* bp of *target_pos*."""
        cigars = read.cigartuples
        if not cigars:
            return False
        if clip_side == 'right' and cigars[-1][0] in CLIP_OPS:
            if abs(read.reference_end - target_pos) < tol:
                return True
        if clip_side == 'left' and cigars[0][0] in CLIP_OPS:
            if abs(read.reference_start - target_pos) < tol:
                return True
        return False

    def _get_background_homozygous_sites(bam_path, chrom, win_start, win_end, ref_seq, min_af=0.95, min_depth=8):
        """Use pysam pileup to identify positions within [win_start, win_end) where a non-reference base reaches allele frequency >= *min_af* and depth >= *min_depth*. Returns a dict mapping relative column index to (genomic_pos, alt_base, alt_count, depth, af)."""
        out = {}
        with pysam.AlignmentFile(bam_path, "rb") as samfile:
            for pileup_col in samfile.pileup(
                chrom, win_start, win_end,
                truncate=True, stepper="all", min_base_quality=0,
                ignore_overlaps=False, ignore_orphans=False,
            ):
                gpos = pileup_col.reference_pos
                if gpos < win_start or gpos >= win_end:
                    continue
                counts = Counter()
                for pr in pileup_col.pileups:
                    if pr.is_refskip or pr.is_del:
                        continue
                    read = pr.alignment
                    if read.is_unmapped or read.query_sequence is None:
                        continue
                    qpos = pr.query_position
                    if qpos is None:
                        continue
                    base = read.query_sequence[qpos].upper()
                    if base != "N":
                        counts[base] += 1
                depth = sum(counts.values())
                if depth < min_depth:
                    continue
                ref_base = ref_seq[gpos - win_start].upper()
                non_ref = Counter({b: n for b, n in counts.items() if b != ref_base})
                if not non_ref:
                    continue
                alt_base, alt_count = non_ref.most_common(1)[0]
                af = alt_count / depth
                if af >= min_af:
                    out[gpos - win_start] = (gpos, alt_base, alt_count, depth, af)
        return out

    if is_deletion:
        left_win_start,  left_win_end  = refined_left - flank, refined_left
        right_win_start, right_win_end = right_bp, right_bp + flank
    else:
        left_win_start,  left_win_end  = refined_left, refined_left + flank
        right_win_start, right_win_end = refined_right - flank, refined_right

    left_background_hom  = _get_background_homozygous_sites(bam_path, chrom, left_win_start,  left_win_end,  ref_left)
    right_background_hom = _get_background_homozygous_sites(bam_path, chrom, right_win_start, right_win_end, ref_right)

    if show_homozygous_variants:
        print(f"Left background homozygous sites: {len(left_background_hom)}")
        print(f"Right background homozygous sites: {len(right_background_hom)}")
        if left_background_hom:
            print("Left homozygous background variants:")
            for col in sorted(left_background_hom):
                gpos, alt, alt_count, depth, af = left_background_hom[col]
                print(f"  col={col:>4} genomic={gpos + 1:>12} ref={ref_left[col]} alt={alt} {alt_count}/{depth} AF={af:.3f}")
        if right_background_hom:
            print("Right homozygous background variants:")
            for col in sorted(right_background_hom):
                gpos, alt, alt_count, depth, af = right_background_hom[col]
                print(f"  col={col:>4} genomic={gpos + 1:>12} ref={ref_right[col]} alt={alt} {alt_count}/{depth} AF={af:.3f}")

    left_seqs,  right_seqs  = [], []
    left_cols,  right_cols  = [], []
    read_names = {r.query_name for r in reads if not r.is_unmapped}
    supp_index, primary_index = _build_sa_index(bam_path, read_names, chrom)

    for read in reads:
        if read.is_unmapped or not read.cigartuples:
            continue
        label = _read_label(read)

        if is_deletion:
            if _primary_clips_near(read, refined_left, 'right'):
                cd = _cigar_walk_to_columns(read, left_win_start, left_win_end)
                if cd:
                    left_cols.append(cd)
                    left_seqs.append((label, _columns_to_padded_seq(cd, flank)))
            else:
                for sa_pos, sa_cigar, sa_strand in _sa_clips_near(read, refined_left, 'right'):
                    supp = supp_index.get((read.query_name, sa_pos))
                    if supp is not None:
                        cd = _cigar_walk_to_columns_from_aln(supp, left_win_start, left_win_end)
                    else:
                        primary = primary_index.get(read.query_name)
                        cd = _cigar_walk_to_columns_sa_primary(primary, sa_pos, sa_cigar, sa_strand, left_win_start, left_win_end) if primary is not None else {}
                    if cd:
                        left_cols.append(cd)
                        left_seqs.append((label, _columns_to_padded_seq(cd, flank)))
                        break
            if _primary_clips_near(read, right_bp, 'left'):
                cd = _cigar_walk_to_columns(read, right_win_start, right_win_end)
                if cd:
                    right_cols.append(cd)
                    right_seqs.append((label, _columns_to_padded_seq(cd, flank)))
            else:
                for sa_pos, sa_cigar, sa_strand in _sa_clips_near(read, right_bp, 'left'):
                    supp = supp_index.get((read.query_name, sa_pos))
                    if supp is not None:
                        cd = _cigar_walk_to_columns_from_aln(supp, right_win_start, right_win_end)
                    else:
                        primary = primary_index.get(read.query_name)
                        cd = _cigar_walk_to_columns_sa_primary(primary, sa_pos, sa_cigar, sa_strand, right_win_start, right_win_end) if primary is not None else {}
                    if cd:
                        right_cols.append(cd)
                        right_seqs.append((label, _columns_to_padded_seq(cd, flank)))
                        break
        else:
            if _primary_clips_near(read, refined_left, 'left'):
                cd = _cigar_walk_to_columns(read, left_win_start, left_win_end)
                if cd:
                    left_cols.append(cd)
                    left_seqs.append((label, _columns_to_padded_seq(cd, flank)))
            else:
                for sa_pos, sa_cigar, sa_strand in _sa_clips_near(read, refined_left, 'left'):
                    supp = supp_index.get((read.query_name, sa_pos))
                    if supp is not None:
                        cd = _cigar_walk_to_columns_from_aln(supp, left_win_start, left_win_end)
                    else:
                        primary = primary_index.get(read.query_name)
                        cd = _cigar_walk_to_columns_sa_primary(primary, sa_pos, sa_cigar, sa_strand, left_win_start, left_win_end) if primary is not None else {}
                    if cd:
                        left_cols.append(cd)
                        left_seqs.append((label, _columns_to_padded_seq(cd, flank)))
                        break
            if _primary_clips_near(read, refined_right, 'right'):
                cd = _cigar_walk_to_columns(read, right_win_start, right_win_end)
                if cd:
                    right_cols.append(cd)
                    right_seqs.append((label, _columns_to_padded_seq(cd, flank)))
            else:
                for sa_pos, sa_cigar, sa_strand in _sa_clips_near(read, refined_right, 'right'):
                    supp = supp_index.get((read.query_name, sa_pos))
                    if supp is not None:
                        cd = _cigar_walk_to_columns_from_aln(supp, right_win_start, right_win_end)
                    else:
                        primary = primary_index.get(read.query_name)
                        cd = _cigar_walk_to_columns_sa_primary(primary, sa_pos, sa_cigar, sa_strand, right_win_start, right_win_end) if primary is not None else {}
                    if cd:
                        right_cols.append(cd)
                        right_seqs.append((label, _columns_to_padded_seq(cd, flank)))
                        break

    if not left_cols:
        print("Warning: no reads found spanning refined_left.")
    if not right_cols:
        print("Warning: no reads found spanning refined_right.")

    left_consensus,  left_low_cov  = _build_consensus_from_columns(left_cols,  ref_left,  flank, min_col_depth=5)
    right_consensus, right_low_cov = _build_consensus_from_columns(right_cols, ref_right, flank, min_col_depth=5)

    print(f"Left consensus  ({len(left_seqs):3d} reads): {left_consensus}")
    print(f"Right consensus ({len(right_seqs):3d} reads): {right_consensus}")
    print(f"Left reference:  {ref_left}")
    print(f"Right reference: {ref_right}")

    left_variant_mask  = [lc != lr for lc, lr in zip(left_consensus,  ref_left)]
    right_variant_mask = [rc != rr for rc, rr in zip(right_consensus, ref_right)]

    for col in left_background_hom:
        if 0 <= col < len(left_variant_mask):
            left_variant_mask[col] = False
    for col in right_background_hom:
        if 0 <= col < len(right_variant_mask):
            right_variant_mask[col] = False
    for col in left_low_cov:
        left_variant_mask[col] = False
    for col in right_low_cov:
        right_variant_mask[col] = False

    if left_low_cov:
        print(f"Left low-coverage positions (forced False): {left_low_cov}")
    if right_low_cov:
        print(f"Right low-coverage positions (forced False): {right_low_cov}")

    if show_pileup:
        _print_pileup("Left breakpoint pileup",  left_seqs,  left_consensus,  ref_left,  flank, after_bp=(not is_deletion))
        _print_pileup("Right breakpoint pileup", right_seqs, right_consensus, ref_right, flank, after_bp=is_deletion)

    if pileup_save_dir:
        os.makedirs(pileup_save_dir, exist_ok=True)
        _lbl = f"_{pair_label}" if pair_label else ""
        _save_pileup_png("Left breakpoint pileup",  left_seqs,  left_consensus,  ref_left,  flank, after_bp=(not is_deletion), save_path=os.path.join(pileup_save_dir, f"pileup_left{_lbl}.png"))
        _save_pileup_png("Right breakpoint pileup", right_seqs, right_consensus, ref_right, flank, after_bp=is_deletion,       save_path=os.path.join(pileup_save_dir, f"pileup_right{_lbl}.png"))

    return (left_consensus, right_consensus, left_variant_mask, right_variant_mask, right_bp, ref_left, ref_right)


# ─────────────────────────────────────────────────────────────────────────────
# 6. Pileup display helpers  (identical to sv_analysis.py)
# ─────────────────────────────────────────────────────────────────────────────

def _print_pileup(title, labeled_seqs, consensus, ref_anchor, flank, after_bp, max_reads=100):
    """Render a breakpoint pileup as an HTML block and display it in a Jupyter notebook. Shows each read's bases aligned to the reference, highlights mismatches, and prints per-column support statistics."""
    if not labeled_seqs:
        display(HTML(f"<pre>=== {title} ===\nNo sequences.</pre>"))
        return
    width = flank
    LABEL_WIDTH = 38

    def fmt_label(label):
        name = label.split(" ")[0]
        return name[:LABEL_WIDTH].ljust(LABEL_WIDTH)

    def align_to_ref(seq, ref):
        w = len(ref)
        result = []
        for i in range(w):
            if i < len(seq) and seq[i] == "-":
                result.append("-")
            elif i < len(seq) and seq[i] not in ("N",):
                result.append(seq[i].upper())
            else:
                result.append(" ")
        return result, {}

    seqs_up = [(label, s.upper()[:width]) for label, s in labeled_seqs if s is not None]
    lines = []
    lines.append(f"=== {title} === ({len(seqs_up)} reads)")
    lines.append(f"{'ref'.ljust(LABEL_WIDTH)} {ref_anchor[:width]}")
    lines.append(f"{chr(32) * LABEL_WIDTH} {chr(39).join(str((i // 10) % 10) for i in range(width))}")
    lines.append(f"{chr(32) * LABEL_WIDTH} {chr(39).join(str(i % 10) for i in range(width))}")
    lines.append(f"{chr(32) * LABEL_WIDTH} {chr(8212) * width}")

    for label, s in seqs_up[:max_reads]:
        result, inserts = align_to_ref(s, ref_anchor[:width])
        annotated_parts = []
        for j, base in enumerate(result):
            ins = inserts.get(j)
            if ins:
                annotated_parts.append(f'<span title="insertion: {ins}" style="color:blue;cursor:help">&</span>')
            if base == "-":
                annotated_parts.append("-")
            elif base == " ":
                annotated_parts.append(" ")
            elif base != ref_anchor[j]:
                annotated_parts.append(f'<span style="color:red">{base.lower()}</span>')
            else:
                annotated_parts.append(base)
        lines.append(f"{fmt_label(label)} {chr(39).join(annotated_parts)}")

    if len(seqs_up) > max_reads:
        lines.append(f"  … ({len(seqs_up) - max_reads} more reads not shown)")
    lines.append(f"{chr(32) * LABEL_WIDTH} {chr(8212) * width}")
    lines.append(f"{'consensus'.ljust(LABEL_WIDTH)} {consensus[:width]}")
    lines.append("\nColumn support:")
    all_results = [align_to_ref(s, ref_anchor[:width])[0] for _, s in seqs_up]
    for i in range(width):
        counts = Counter(r[i] for r in all_results if i < len(r) and r[i] not in (" ", "N"))
        total = sum(counts.values())
        major = counts.most_common(1)[0][0] if total else ref_anchor[i]
        frac  = counts.most_common(1)[0][1] / total if total else 0.0
        col_label = f"+{i:03d}" if after_bp else f"-{flank - 1 - i:03d}"
        lines.append(
            f"  {col_label} ref={ref_anchor[i]} cons={consensus[i] if i < len(consensus) else '?'}"
            f" top={major}:{frac:.2f} {dict(counts)}"
        )

    html_lines = []
    for line in lines:
        matched = False
        for label, _ in seqs_up[:max_reads]:
            truncated = fmt_label(label)
            if line.startswith(truncated.rstrip()):
                full_label = label.replace('"', '&quot;')
                rest = line[LABEL_WIDTH:]
                html_lines.append(f'<span title="{full_label}" style="cursor:help; background:#fffbe6">{truncated}</span>{rest}')
                matched = True
                break
        if not matched:
            html_lines.append(line)

    html = '<pre style="font-family: monospace; font-size: 13px; line-height: 1.4; '
    html += 'background: #f8f8f8; padding: 10px; border: 1px solid #ddd; '
    html += 'overflow-x: auto; max-height: 600px; overflow-y: auto">'
    html += "\n".join(html_lines) + "</pre>"
    display(HTML(html))


def _save_pileup_png(title, labeled_seqs, consensus, ref_anchor, flank, after_bp, save_path, max_reads=100):
    """Render a breakpoint pileup as a PNG file using matplotlib. Mismatching bases are colour-coded per nucleotide; the consensus is shown at the bottom."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    BASE_COLOR = {"A": "#27ae60", "C": "#2980b9", "G": "#f39c12", "T": "#c0392b", "-": "#7f8c8d"}
    width = flank
    LABEL_WIDTH = 38
    seqs_up = [(label, s.upper()[:width]) for label, s in labeled_seqs if s is not None]
    n_reads = min(len(seqs_up), max_reads)
    n_rows  = n_reads + 6
    row_h   = 0.18
    fig_w   = max(14, flank * 0.085 + 4)
    fig_h   = max(3.0, n_rows * row_h + 1.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.axis("off")
    FONT = {"fontfamily": "DejaVu Sans Mono", "fontsize": 7}
    FONT_TITLE = {"fontfamily": "DejaVu Sans Mono", "fontsize": 8}
    line_ys = list(range(n_rows, 0, -1))

    def fmt_label(label):
        return label.split(" ")[0][:LABEL_WIDTH].ljust(LABEL_WIDTH)

    y_idx = 0
    def next_y():
        nonlocal y_idx
        y = line_ys[y_idx]; y_idx += 1; return y

    ax.set_xlim(0, LABEL_WIDTH + width + 2)
    ax.set_ylim(0, n_rows + 2)
    ax.text(0, n_rows + 1.5, f"=== {title} === ({len(seqs_up)} reads)", **FONT_TITLE, fontweight="bold", va="top")
    y = next_y()
    ax.text(0, y, "ref".ljust(LABEL_WIDTH) + " " + ref_anchor[:width], **FONT, va="center")
    y = next_y()
    ax.text(0, y, " " * LABEL_WIDTH + " " + "".join(str((i // 10) % 10) for i in range(width)), **FONT, va="center", color="#888888")
    y = next_y()
    ax.text(0, y, " " * LABEL_WIDTH + " " + "".join(str(i % 10) for i in range(width)), **FONT, va="center", color="#888888")
    y = next_y()
    ax.text(0, y, " " * LABEL_WIDTH + " " + "─" * width, **FONT, va="center", color="#cccccc")
    for label, s in seqs_up[:max_reads]:
        y = next_y()
        ax.text(0, y, fmt_label(label) + " ", **FONT, va="center", color="#555555")
        x_off = LABEL_WIDTH + 2
        for col_i in range(width):
            if col_i >= len(s): break
            base = s[col_i].upper()
            ref_base = ref_anchor[col_i] if col_i < len(ref_anchor) else " "
            is_mismatch = base not in (" ", "N") and base != ref_base
            color = BASE_COLOR.get(base, "black") if is_mismatch else "#222222"
            ax.text(x_off + col_i, y, base, **FONT, va="center", ha="left", color=color)
    y = next_y()
    ax.text(0, y, " " * LABEL_WIDTH + " " + "─" * width, **FONT, va="center", color="#cccccc")
    y = next_y()
    ax.text(0, y, "consensus".ljust(LABEL_WIDTH) + " " + consensus[:width], **FONT, va="center", fontweight="bold")
    plt.tight_layout(pad=0.2)
    plt.savefig(save_path, dpi=120, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved pileup PNG: {save_path}")


# ─────────────────────────────────────────────────────────────────────────────
# 7. MMEJ microhomology detection  (identical to sv_analysis.py)
# ─────────────────────────────────────────────────────────────────────────────

def find_mmej_microhomology(
    left_consensus, left_variant_mask,
    right_consensus, right_variant_mask,
    ref_left, ref_right,
    refined_left, refined_right,
    flank, sv_type, min_mh=2,
    max_edge_suffix_prefix_mismatches=1,
    max_edge_shift=None,
    min_block_identity=0.8,
    edge_check=10,
    edge_buffer=10,
    fasta=None, chrom=None,
    save_dir=None, pair_label=None,
):
    """Detect MMEJ microhomology at an SV breakpoint by performing a suffix-prefix overlap
    search between the left and right consensus sequences. Variants within the matching
    block are reclassified against the opposite reference to distinguish junction-specific
    bases from haplotype-specific SNVs.

    Two code paths exist:
      1. The main scored suffix-prefix search (tolerant of mismatches via +2/-1 scoring),
         used when the matched block still carries reclassified variants.
      2. An "all-reference" edge search, used when every position in the matched block is
         confirmed non-variant on both sides. This path compares ref_left_aln vs
         ref_right_aln directly (not consensus) since microhomology is fundamentally a
         reference-genome property, and tolerates up to *max_edge_suffix_prefix_mismatches*
         mismatches so genuine single-base differences between the two flanks (e.g. TCTC vs
         TCTT) aren't truncated to a shorter exact match.

         Within this edge search, the match is first attempted anchored exactly at the
         soft-clip point (shift=0). If that anchored match doesn't reach *min_mh* matched
         bases, the search retries with the match window shifted progressively further away
         from the soft-clip point on BOTH flanks simultaneously (shift=1, 2, ... up to
         *max_edge_shift*). This accounts for cases where the true homology run extends
         beyond, or sits slightly offset from, the annotated soft-clip position (the
         soft-clip can land anywhere within a repeat/homology run, not necessarily at its
         outer edge). Applied identically for deletions (left_tail/right_head) and
         duplications (right_tail/left_head).

    Returns
    -------
    (mh_len, mh_seq_left, mh_seq_right, n_mismatches,
     left_mh_start, left_mh_end, right_mh_start, right_mh_end)

    mh_seq_left / mh_seq_right are reported separately (rather than a single mh_seq) so any
    mismatch between the two homology arms is visible, e.g. left="TCTC", right="TCTT".
    """
    is_deletion = sv_type == "deletion"
    left_win_start  = (refined_left  - flank) if is_deletion else refined_left
    right_win_start = refined_right if is_deletion else (refined_right - flank)

    def _strip_gaps(seq, mask, ref=None):
        """Remove gap characters ('-') from *seq* and *mask*, optionally applying the same
        filter to a *ref* string. Returns stripped copies of all supplied sequences/masks,
        plus *orig_idx*: orig_idx[stripped_i] = the index that position held in the
        pre-strip (raw window-offset) sequence, so genomic coordinates computed from
        stripped-array indices can be mapped back to true window offsets instead of
        silently drifting by one position per gap removed."""
        new_seq, new_mask, new_ref, orig_idx = [], [], [], []
        for i, (b, m) in enumerate(zip(seq, mask)):
            if b != "-":
                new_seq.append(b)
                new_mask.append(m)
                orig_idx.append(i)
                if ref is not None:
                    new_ref.append(ref[i] if i < len(ref) else 'N')
        if ref is not None:
            return "".join(new_seq), new_mask, "".join(new_ref), orig_idx
        return "".join(new_seq), new_mask, orig_idx

    left_consensus_aln,  left_variant_mask_aln,  ref_left_aln,  left_orig_idx  = _strip_gaps(left_consensus,  left_variant_mask,  ref_left)
    right_consensus_aln, right_variant_mask_aln, ref_right_aln, right_orig_idx = _strip_gaps(right_consensus, right_variant_mask, ref_right)

    def _normalise(seq, mask):
        """Replace variant positions (mask[i] == True) with 'X' so they are treated as
        wildcards during the microhomology overlap search."""
        return "".join('X' if (i < len(mask) and mask[i]) else b for i, b in enumerate(seq))

    def _refine_block_by_direct_match(l_arr, r_arr, l_anchor, r_anchor, length, max_shift, min_len):
        """Local, bounded re-anchoring: if l_arr/r_arr don't agree well over
        [l_anchor:l_anchor+length] vs [r_anchor:r_anchor+length], search small
        independent shifts (+/- max_shift) of the two start positions for the
        longest contiguous exact-match run near the original anchor. Bounded to
        a small neighbourhood (not the whole flank) so it can only correct a few
        bases of breakpoint-calling slop, not wander into an unrelated match
        elsewhere in a long repeat. Returns (l_start, r_start, run_len) for the
        best run found if it beats the original block's match count, else None."""
        orig_mismatches = sum(
            1 for a, b in zip(l_arr[l_anchor:l_anchor + length], r_arr[r_anchor:r_anchor + length])
            if a != b
        )
        orig_matches = length - orig_mismatches

        best = None  # (run_len, l_start, r_start)
        for dl in range(-max_shift, max_shift + 1):
            l0 = l_anchor + dl
            if l0 < 0 or l0 >= len(l_arr):
                continue
            for dr in range(-max_shift, max_shift + 1):
                r0 = r_anchor + dr
                if r0 < 0 or r0 >= len(r_arr):
                    continue
                max_run = min(len(l_arr) - l0, len(r_arr) - r0)
                n = 0
                while n < max_run and l_arr[l0 + n] == r_arr[r0 + n]:
                    n += 1
                if n >= min_len and (best is None or n > best[0]):
                    best = (n, l0, r0)

        if best is not None and best[0] > orig_matches:
            return best
        return None

    norm_left  = _normalise(left_consensus_aln,  left_variant_mask_aln)
    norm_right = _normalise(right_consensus_aln, right_variant_mask_aln)
    EDGE_SEARCH = min(len(norm_left), len(norm_right), flank)

    if is_deletion:
        query = norm_left[-EDGE_SEARCH:]; ref_seq = norm_right[:EDGE_SEARCH]
        left_offset = len(norm_left) - EDGE_SEARCH; right_offset = 0
    else:
        query = norm_right[-EDGE_SEARCH:]; ref_seq = norm_left[:EDGE_SEARCH]
        left_offset = 0; right_offset = len(norm_right) - EDGE_SEARCH

    best_score, best_overlap = -1, 0
    best_q_slice, best_r_slice, best_mismatches = "", "", 0
    for overlap in range(min_mh, min(len(query), len(ref_seq)) + 1):
        q_slice = query[-overlap:]; r_slice = ref_seq[:overlap]
        n_mismatch = sum(1 for q, r in zip(q_slice, r_slice) if q != r)
        score = 2 * (overlap - n_mismatch) - n_mismatch
        if score > best_score:
            best_score = score; best_overlap = overlap
            best_q_slice, best_r_slice, best_mismatches = q_slice, r_slice, n_mismatch

    best_identity = (best_overlap - best_mismatches) / best_overlap if best_overlap else 0.0
    print(f"\n--- Suffix-prefix overlap search (EDGE_SEARCH={EDGE_SEARCH}) ---")
    print(f"  best_overlap={best_overlap}  best_score={best_score}  "
          f"mismatches={best_mismatches}  identity={best_identity:.2f}")
    query_side, ref_side = ("left", "right") if is_deletion else ("right", "left")
    print(f"  query ({query_side} flank, breakpoint-anchored suffix): {best_q_slice!r}")
    print(f"  ref   ({ref_side} flank, breakpoint-anchored prefix): {best_r_slice!r}")
    if best_overlap < min_mh or best_score <= 0:
        print("No matching block found.")
        return 0, "", "", 0, None, None, None, None

    if is_deletion:
        best_l_start = left_offset + (EDGE_SEARCH - best_overlap); best_l_end = left_offset + EDGE_SEARCH
        best_r_start = right_offset; best_r_end = right_offset + best_overlap
    else:
        best_l_start = left_offset; best_l_end = left_offset + best_overlap
        best_r_start = right_offset + (EDGE_SEARCH - best_overlap); best_r_end = right_offset + EDGE_SEARCH

    best_k = best_overlap
    print(f"\n--- Best matching block (len={best_k}) ---")
    print(f"  Left  [{best_l_start}:{best_l_end}]  genomic {left_win_start + left_orig_idx[best_l_start]}–{left_win_start + left_orig_idx[best_l_end - 1]}")
    print(f"  Right [{best_r_start}:{best_r_end}]  genomic {right_win_start + right_orig_idx[best_r_start]}–{right_win_start + right_orig_idx[best_r_end - 1]}")

    left_variant_mask_aln  = list(left_variant_mask_aln)
    right_variant_mask_aln = list(right_variant_mask_aln)
    block_table = []
    print("\n--- Reclassifying variants within block against opposite reference ---")
    for k in range(best_k):
        l_pos = best_l_start + k; r_pos = best_r_start + k
        l_base = left_consensus_aln[l_pos]; r_base = right_consensus_aln[r_pos]
        rl_base = ref_left_aln[l_pos]  if l_pos < len(ref_left_aln)  else 'N'
        rr_base = ref_right_aln[r_pos] if r_pos < len(ref_right_aln) else 'N'
        l_var_pre = left_variant_mask_aln[l_pos]; r_var_pre = right_variant_mask_aln[r_pos]
        note_parts = []
        if l_var_pre:
            if l_base == rr_base: note_parts.append("L:junction")
            else: left_variant_mask_aln[l_pos] = False; note_parts.append("L:hap-specific→False")
        if r_var_pre:
            if r_base == rl_base: note_parts.append("R:junction")
            else: right_variant_mask_aln[r_pos] = False; note_parts.append("R:hap-specific→False")
        if not note_parts: note_parts.append("ref-match")
        block_table.append({
            "k": k, "l_pos": l_pos, "r_pos": r_pos,
            "l_genomic": left_win_start + left_orig_idx[l_pos] + 1, "r_genomic": right_win_start + right_orig_idx[r_pos] + 1,
            "l_cons": l_base, "r_cons": r_base, "ref_left": rl_base, "ref_right": rr_base,
            "l_var_final": left_variant_mask_aln[l_pos], "r_var_final": right_variant_mask_aln[r_pos],
            "note": " ".join(note_parts),
        })

    block_df = pd.DataFrame(block_table)
    _block_csv = (
        os.path.join(save_dir, f"block_table_{pair_label}.csv") if save_dir and pair_label
        else os.path.join(save_dir, "block_table.csv") if save_dir
        else "block_table.csv"
    )
    block_df.to_csv(_block_csv, index=False)
    print(f"Saved block.csv ({len(block_df)} rows)")

    left_vars_in_block  = left_variant_mask_aln[best_l_start : best_l_start + best_k]
    right_vars_in_block = right_variant_mask_aln[best_r_start : best_r_start + best_k]
    print("left variant mask in block: ",  left_vars_in_block)
    print("right variant mask in block: ", right_vars_in_block)

    if not any(left_vars_in_block) and not any(right_vars_in_block):
        print("\n--- All-reference case (post-reclassification): edge suffix-prefix search. ---")

        def _fetch_buffer(start, end):
            """Fetch edge_buffer bases reaching INTO the deleted region (deletions) or just
            outside the duplicated span (duplications). Returns "" if unavailable."""
            if fasta is None or chrom is None or edge_buffer <= 0:
                return ""
            return _safe_fetch(fasta, chrom, start, end)

        def _best_match_scan(a, a_suffix, b, b_suffix, max_mm, min_len):
            """Find the longest k where a/b (each sliced as a[-k:] if a_suffix else a[:k],
            and likewise for b) has >= min_len true matches and <= max_mm mismatches.

            a_suffix/b_suffix must match each string's own natural boundary-adjacency:
            a "tail"-style string (boundary-adjacent base at its END, e.g. left_tail, or a
            buffer fetched immediately BEFORE a clip point) uses suffix slicing (a[-k:]).
            A "head"-style string (boundary-adjacent base at its START, e.g. right_head, or
            a buffer fetched immediately AFTER a clip point) uses prefix slicing (a[:k]).
            Mixing conventions incorrectly (e.g. reversing a suffix-style buffer and then
            prefix-slicing it) silently compares the wrong bases -- this function requires
            the caller to specify the correct convention explicitly instead.

            Each k is tested independently (anchored exactly at the boundary; no shift).
            Returns (k, n_mismatches)."""
            max_len = min(len(a), len(b))
            for k in range(max_len, min_len - 1, -1):
                a_k = a[-k:] if a_suffix else a[:k]
                b_k = b[-k:] if b_suffix else b[:k]
                n_mm = sum(1 for x, y in zip(a_k, b_k) if x != y)
                n_matches = k - n_mm
                if n_mm <= max_mm and n_matches >= min_len:
                    return k, n_mm
            return 0, 0

        if is_deletion:
            left_tail  = ref_left_aln[-edge_check:]
            right_head = ref_right_aln[:edge_check]

            # Attempt 1: forward / forward -- left_tail is suffix-style, right_head is prefix-style
            edge_mh_len, edge_mismatches = _best_match_scan(
                left_tail, True, right_head, False, max_edge_suffix_prefix_mismatches, min_mh
            )
            edge_mode = "forward/forward"

            # Attempt 2 (only if 1 failed): keep left_tail as-is, test right_buf_raw --
            # a buffer fetched immediately BEFORE right_bp is suffix-style (boundary-adjacent
            # base at its END), same convention as left_tail. No reversal needed.
            if edge_mh_len == 0:
                right_buf_raw = _fetch_buffer(refined_right - edge_buffer, refined_right)
                if right_buf_raw:
                    edge_mh_len, edge_mismatches = _best_match_scan(
                        left_tail, True, right_buf_raw, True, max_edge_suffix_prefix_mismatches, min_mh
                    )
                    if edge_mh_len:
                        edge_mode = "forward-left/reverse-right"

            # Attempt 3 (only if 2 also failed): keep right_head as-is, test left_buf_raw --
            # a buffer fetched immediately AFTER refined_left is prefix-style (boundary-adjacent
            # base at its START), same convention as right_head. No reversal needed.
            if edge_mh_len == 0:
                left_buf_raw = _fetch_buffer(refined_left, refined_left + edge_buffer)
                if left_buf_raw:
                    edge_mh_len, edge_mismatches = _best_match_scan(
                        left_buf_raw, False, right_head, False, max_edge_suffix_prefix_mismatches, min_mh
                    )
                    if edge_mh_len:
                        edge_mode = "reverse-left/forward-right"

            if edge_mode == "forward-left/reverse-right":
                mh_seq_left  = left_tail[-edge_mh_len:] if edge_mh_len else ""
                mh_seq_right = right_buf_raw[-edge_mh_len:] if edge_mh_len else ""
                left_mh_genomic_start  = refined_left - edge_mh_len
                left_mh_genomic_end    = refined_left - 1
                right_mh_genomic_start = refined_right - edge_mh_len
                right_mh_genomic_end   = refined_right - 1
            elif edge_mode == "reverse-left/forward-right":
                mh_seq_left  = left_buf_raw[:edge_mh_len] if edge_mh_len else ""
                mh_seq_right = right_head[:edge_mh_len] if edge_mh_len else ""
                left_mh_genomic_start  = refined_left
                left_mh_genomic_end    = refined_left + edge_mh_len - 1
                right_mh_genomic_start = refined_right
                right_mh_genomic_end   = refined_right + edge_mh_len - 1
            else:  # forward/forward (whether it succeeded or all three attempts failed)
                mh_seq_left  = left_tail[-edge_mh_len:] if edge_mh_len else ""
                mh_seq_right = right_head[:edge_mh_len] if edge_mh_len else ""
                left_mh_genomic_start  = refined_left - edge_mh_len
                left_mh_genomic_end    = refined_left - 1
                right_mh_genomic_start = refined_right
                right_mh_genomic_end   = refined_right + edge_mh_len - 1

        else: # if it's a duplication...
            right_tail = ref_right_aln[-edge_check:]
            left_head  = ref_left_aln[:edge_check]

            # Attempt 1: forward / forward -- right_tail is suffix-style, left_head is prefix-style
            edge_mh_len, edge_mismatches = _best_match_scan(
                right_tail, True, left_head, False, max_edge_suffix_prefix_mismatches, min_mh
            )
            edge_mode = "forward/forward"

            # Attempt 2 (only if 1 failed): a buffer fetched immediately BEFORE refined_left
            # is suffix-style, same convention as right_tail. No reversal needed.
            if edge_mh_len == 0:
                left_buf_raw = _fetch_buffer(refined_left - edge_buffer, refined_left)
                if left_buf_raw:
                    edge_mh_len, edge_mismatches = _best_match_scan(
                        right_tail, True, left_buf_raw, True, max_edge_suffix_prefix_mismatches, min_mh
                    )
                    if edge_mh_len:
                        edge_mode = "forward-right/reverse-left"

            # Attempt 3 (only if 2 also failed): a buffer fetched immediately AFTER refined_right
            # is prefix-style, same convention as left_head. No reversal needed.
            if edge_mh_len == 0:
                right_buf_raw = _fetch_buffer(refined_right, refined_right + edge_buffer)
                if right_buf_raw:
                    edge_mh_len, edge_mismatches = _best_match_scan(
                        right_buf_raw, False, left_head, False, max_edge_suffix_prefix_mismatches, min_mh
                    )
                    if edge_mh_len:
                        edge_mode = "reverse-right/forward-left"

            if edge_mode == "forward-right/reverse-left":
                mh_seq_right = right_tail[-edge_mh_len:] if edge_mh_len else ""
                mh_seq_left  = left_buf_raw[-edge_mh_len:] if edge_mh_len else ""
                right_mh_genomic_start = refined_right - edge_mh_len
                right_mh_genomic_end   = refined_right - 1
                left_mh_genomic_start  = refined_left - edge_mh_len
                left_mh_genomic_end    = refined_left - 1
            elif edge_mode == "reverse-right/forward-left":
                mh_seq_right = right_buf_raw[:edge_mh_len] if edge_mh_len else ""
                mh_seq_left  = left_head[:edge_mh_len] if edge_mh_len else ""
                right_mh_genomic_start = refined_right
                right_mh_genomic_end   = refined_right + edge_mh_len - 1
                left_mh_genomic_start  = refined_left
                left_mh_genomic_end    = refined_left + edge_mh_len - 1
            else:  # forward/forward
                mh_seq_right = right_tail[-edge_mh_len:] if edge_mh_len else ""
                mh_seq_left  = left_head[:edge_mh_len] if edge_mh_len else ""
                right_mh_genomic_start = refined_right - edge_mh_len
                right_mh_genomic_end   = refined_right - 1
                left_mh_genomic_start  = refined_left
                left_mh_genomic_end    = refined_left + edge_mh_len - 1

        if edge_mode != "forward/forward" and edge_mh_len:
            print(f"    [edge-search] forward search insufficient; succeeded via {edge_mode} → "
                  f"{edge_mh_len} bp ({edge_mh_len - edge_mismatches} matches, {edge_mismatches} mismatches)")

        print(f"Edge MH detected: {edge_mh_len} bp (mismatches={edge_mismatches}, mode={edge_mode})  "
              f"left={mh_seq_left!r}  right={mh_seq_right!r}")
        print(f"  → Left MH genomic:  {left_mh_genomic_start}–{left_mh_genomic_end}")
        print(f"  → Right MH genomic: {right_mh_genomic_start}–{right_mh_genomic_end}")
        return (edge_mh_len, mh_seq_left, mh_seq_right, edge_mismatches,
                left_mh_genomic_start, left_mh_genomic_end,
                right_mh_genomic_start, right_mh_genomic_end)

    print(f"\n--- Variant walk ({'deletion' if is_deletion else 'duplication'} mode) ---")
    if is_deletion:
        first_left_var = next((i for i, v in enumerate(left_vars_in_block)  if v), -1)
        last_right_var = max((i for i, v in enumerate(right_vars_in_block) if v), default=-1)
        print(f"  Initial: first_left_var={first_left_var}  last_right_var={last_right_var}")
        if last_right_var >= first_left_var >= 0:
            last_right_var = max((i for i in range(first_left_var - 1, -1, -1) if right_vars_in_block[i]), default=-1)
            print(f"  Adjusted: last_right_var={last_right_var}")
        if first_left_var <= last_right_var >= 0:
            first_left_var = next((i for i in range(last_right_var + 1, best_k) if left_vars_in_block[i]), -1)
            print(f"  Adjusted: first_left_var={first_left_var}")
        mh_block_start = last_right_var + 1 if last_right_var >= 0 else 0
        mh_block_end   = first_left_var - 1 if first_left_var >= 0 else (best_k - 1)
    else:
        first_right_var = next((i for i, v in enumerate(right_vars_in_block) if v), -1)
        last_left_var   = max((i for i, v in enumerate(left_vars_in_block)  if v), default=-1)
        print(f"  Initial: last_left_var={last_left_var}  first_right_var={first_right_var}")
        if last_left_var >= first_right_var >= 0:
            last_left_var = max((i for i in range(first_right_var - 1, -1, -1) if left_vars_in_block[i]), default=-1)
            print(f"  Adjusted: last_left_var={last_left_var}")
        if first_right_var <= last_left_var >= 0:
            first_right_var = next((i for i in range(last_left_var + 1, best_k) if right_vars_in_block[i]), -1)
            print(f"  Adjusted: first_right_var={first_right_var}")
        mh_block_start = last_left_var + 1
        mh_block_end   = (first_right_var - 1) if first_right_var >= 0 else (best_k - 1)

    mh_len = mh_block_end - mh_block_start + 1
    if mh_len < min_mh:
        print(f"  → mh_len={mh_len} below min_mh={min_mh}, no MH called.")
        return 0, "", "", 0, None, None, None, None

    left_mh_idx  = best_l_start + mh_block_start
    right_mh_idx = best_r_start + mh_block_start
    mh_seq_left  = left_consensus_aln[left_mh_idx  : left_mh_idx  + mh_len]
    mh_seq_right = right_consensus_aln[right_mh_idx : right_mh_idx + mh_len]
    n_mismatches = sum(1 for a, b in zip(mh_seq_left, mh_seq_right) if a != b)

    block_identity = (mh_len - n_mismatches) / mh_len
    if block_identity < min_block_identity:
        shift = max_edge_shift if max_edge_shift is not None else 15
        print(f"  [shift-refine] block identity={block_identity:.2f} < {min_block_identity} "
              f"-> searching +/-{shift}bp local re-anchor ...")
        refined = _refine_block_by_direct_match(
            left_consensus_aln, right_consensus_aln,
            left_mh_idx, right_mh_idx, mh_len,
            max_shift=shift, min_len=min_mh,
        )
        if refined is not None:
            run_len, l0, r0 = refined
            new_left  = left_consensus_aln[l0:l0 + run_len]
            new_right = right_consensus_aln[r0:r0 + run_len]
            new_mismatches = sum(1 for a, b in zip(new_left, new_right) if a != b)
            print(f"  [shift-refine] found better local match: {run_len}bp (mismatches={new_mismatches}) "
                  f"vs original {mh_len}bp (mismatches={n_mismatches})")
            mh_len, mh_seq_left, mh_seq_right, n_mismatches = run_len, new_left, new_right, new_mismatches
            left_mh_idx, right_mh_idx = l0, r0
        else:
            print("  [shift-refine] no better local match found within the shift window; keeping original block.")

    left_mh_genomic_start  = left_win_start  + left_orig_idx[left_mh_idx]
    left_mh_genomic_end    = left_win_start  + left_orig_idx[left_mh_idx  + mh_len - 1]
    right_mh_genomic_start = right_win_start + right_orig_idx[right_mh_idx]
    right_mh_genomic_end   = right_win_start + right_orig_idx[right_mh_idx + mh_len - 1]
    print(f"\n → MH block: {mh_len} bp, mismatches={n_mismatches}")
    print(f" → Left MH genomic:  {left_mh_genomic_start} – {left_mh_genomic_end}")
    print(f" → Right MH genomic: {right_mh_genomic_start} – {right_mh_genomic_end}")
    print(f" → Left  MH seq: {mh_seq_left!r}")
    print(f" → Right MH seq: {mh_seq_right!r}")
    return (mh_len, mh_seq_left, mh_seq_right, n_mismatches,
            left_mh_genomic_start, left_mh_genomic_end,
            right_mh_genomic_start, right_mh_genomic_end)



# ─────────────────────────────────────────────────────────────────────────────
# 8. Visualisation  (identical to sv_analysis.py)
# ─────────────────────────────────────────────────────────────────────────────

def plot_microhomology(
    mh_seq_left, mh_seq_right,
    left_mh_genomic_start, right_mh_genomic_start, right_mh_genomic_end,
    left_junc_left_mh, left_junc_right_mh, right_junc_left_mh, right_junc_right_mh,
    sv_type, sample_id="", gene="", chrom="", save_path="microhomology_viz.png",
):
    """Generate and save a publication-style PNG visualising the proposed SV structure with
    the microhomology sequence highlighted. Shows the WT genomic context on both sides, pipe
    characters flanking the MH block, a brace indicating the deleted/duplicated segment, and
    a title with sample/gene annotation.

    mh_seq_left and mh_seq_right are rendered separately at their respective junction (rather
    than a single shared mh_seq), and any base where they differ is highlighted in red at
    both occurrences so a single-base mismatch between the two homology arms (e.g. left=TCTC,
    right=TCTT) is visible directly in the figure.
    """
    MONO = "DejaVu Sans Mono"
    SEQ_FS = 10.5
    MISMATCH_COLOR = "crimson"
    MATCH_COLOR    = "black"

    mh_seq_left  = (mh_seq_left  or "").lower()
    mh_seq_right = (mh_seq_right or "").lower()
    mh_len = max(len(mh_seq_left), len(mh_seq_right))

    def _mh_char_pieces(seq, other):
        """Return a list of (char, bold=True, color) tuples for *seq*, colouring any base
        that differs from the aligned base in *other* as a mismatch."""
        parts = []
        for i, ch in enumerate(seq):
            other_ch = other[i] if i < len(other) else None
            color = MISMATCH_COLOR if (other_ch is not None and ch != other_ch) else MATCH_COLOR
            parts.append((ch, True, color))
        return parts

    left_mh_pieces  = _mh_char_pieces(mh_seq_left,  mh_seq_right)
    right_mh_pieces = _mh_char_pieces(mh_seq_right, mh_seq_left)

    left_mh_genomic_end = left_mh_genomic_start + mh_len - 1   # FIXED: true end, was off by one
    interior_bp = (
        (right_mh_genomic_start - len(right_junc_left_mh))    # FIXED: now excludes right-side context too
        - (left_mh_genomic_end + len(left_junc_right_mh))
        - 1
    )
    gap_str = f" –({interior_bp:,}bp)– "
    chrom_d  = chrom or "chr?"
    sv_label = "duplication" if sv_type != "deletion" else "deletion"
    dup_s    = left_mh_genomic_start
    dup_e    = right_mh_genomic_end
    seq_s    = left_mh_genomic_start - len(left_junc_left_mh) - 1
    seq_e    = right_mh_genomic_end  + len(right_junc_right_mh) + 1
    wt_pre   = f"WT {chrom_d}:{seq_s:,} – "
    wt_suf   = f" – {chrom_d}:{seq_e:,}"

    pieces = [
        (wt_pre, False, "black"),
        (left_junc_left_mh.lower(), False, "black"), ("|", False, "black"),
        *left_mh_pieces, ("|", False, "black"),
        (left_junc_right_mh.lower(), False, "black"),
        (gap_str, False, "black"),
        (right_junc_left_mh.lower(), False, "black"), ("|", False, "black"),
        *right_mh_pieces, ("|", False, "black"),
        (right_junc_right_mh.lower(), False, "black"),
        (wt_suf, False, "black"),
    ]

    total_chars = sum(len(t) for t, _, _ in pieces)
    CHAR_W_IN = SEQ_FS * 0.6 / 72.0
    fig_w = max(14.0, total_chars * CHAR_W_IN + 1.2)
    fig, ax = plt.subplots(figsize=(fig_w, 3.6))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white"); ax.axis("off")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_bb    = ax.get_window_extent(renderer=renderer)

    def char_frac(txt, bold):
        t = ax.text(0, -5, txt, fontsize=SEQ_FS, fontfamily=MONO,
                    fontweight="bold" if bold else "normal",
                    transform=ax.transAxes, ha="left", va="center")
        bb = t.get_window_extent(renderer=renderer)
        t.remove()
        return bb.width / ax_bb.width

    widths = [char_frac(txt, bold) for txt, bold, _ in pieces]
    total_w = sum(widths)
    SEQ_Y   = 0.22
    x_start = (1.0 - total_w) / 2.0
    x = x_start
    pipe_x = []
    for i, (txt, bold, _) in enumerate(pieces):
        if txt == "|": pipe_x.append(x + widths[i] / 2)
        x += widths[i]
    brace_l   = pipe_x[0]; brace_r = pipe_x[3]; brace_mid = (brace_l + brace_r) / 2
    gene_str  = f" $\it{{{gene}}}$" if gene else ""
    title_str = (
        f"Proposed {sv_label} structure"
        + (f" for {sample_id}" if sample_id else "")
        + gene_str
        + (f" tandem {sv_label}" if sv_type != "deletion" else "")
    )
    ax.text(0.5, 0.97, title_str, ha="center", va="top", fontsize=13, fontweight="bold",
            fontfamily="DejaVu Sans", color="black", transform=ax.transAxes)
    seg_lbl = (
        f"{'Duplicated' if sv_type != 'deletion' else 'Deleted'} segment"
        f" – {chrom_d}:{dup_s:,} –{dup_e:,}"
    )
    SEG_Y = 0.72
    ax.text(brace_mid, SEG_Y, seg_lbl, ha="center", va="top", fontsize=10,
            fontfamily=MONO, color="black", transform=ax.transAxes)
    BRACE_Y = 0.52; TICK_H = 0.30
    ax.plot([brace_l, brace_r], [BRACE_Y, BRACE_Y], color="black", lw=1.5, transform=ax.transAxes, clip_on=False)
    ax.plot([brace_l, brace_l], [BRACE_Y, BRACE_Y - TICK_H], color="black", lw=1.5, transform=ax.transAxes, clip_on=False)
    ax.plot([brace_r, brace_r], [BRACE_Y, BRACE_Y - TICK_H], color="black", lw=1.5, transform=ax.transAxes, clip_on=False)
    ax.annotate("", xy=(brace_mid, BRACE_Y), xytext=(brace_mid, SEG_Y - 0.02),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-", color="black", lw=0.8))
    x = x_start
    for (txt, bold, color), w in zip(pieces, widths):
        if txt:
            ax.text(x, SEQ_Y, txt, fontsize=SEQ_FS, fontfamily=MONO,
                    fontweight="bold" if bold else "normal", color=color,
                    ha="left", va="center", transform=ax.transAxes)
        x += w
    plt.tight_layout(pad=0.4)
    plt.savefig(save_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved: {save_path}")
    display(Image(save_path))


# ─────────────────────────────────────────────────────────────────────────────
# 9. Per-cluster analysis
# ─────────────────────────────────────────────────────────────────────────────


def estimate_flank_from_primary_alignments(
    reads,
    min_flank=600,
    flank_fraction=0.5,
):
    """
    Estimate a single flank size from the median reference-consuming length
    of primary alignments from split reads.

    The flank is set to flank_fraction x median read length, subject to a
    minimum of min_flank.

    Returns
    -------
    flank : int
        Estimated flank size (bp), at least min_flank.
    """
    REF_OPS = {0, 2, 3, 7, 8}  # M, D, N, =, X
    lengths = []

    for read in reads:
        if read.is_unmapped or not read.cigartuples:
            continue
        if read.is_supplementary or read.is_secondary:
            continue

        ref_len = sum(length for op, length in read.cigartuples if op in REF_OPS)
        if ref_len > 0:
            lengths.append(ref_len)

    if not lengths:
        flank = min_flank
        med_len = None
    else:
        med_len = median(lengths)
        flank = max(min_flank, int(flank_fraction * med_len))

    print(
        f"  [flank] median_split_read_len={med_len} "
        f"flank_fraction={flank_fraction} "
        f"n_primary_alns={len(lengths)} "
        f"min_flank={min_flank}"
    )
    print(f"  [flank] flank={flank}")

    return flank

def analyze_sv_cluster(
    bam_path, ref_path, cluster,
    flank=None, min_mh=2, show_pileup=False, min_mapq=60,
    pair_label="", save_dir=".", show_homozygous_variants=False,
):
    """Run the full MH pipeline for one merged soft-clip cluster.

    The cluster dict carries refined_left / refined_right from the widest
    member pair.  All reads are extracted over the envelope span so that
    reads from narrower-span pairs are still collected, but the consensus
    windows are anchored to the wider coordinates.
    """
    chrom         = cluster["chrom"]
    refined_left  = cluster["refined_left"]
    refined_right = cluster["refined_right"]
    env_left      = cluster["refined_left"]
    env_right     = cluster["refined_right"]

    # Use the envelope for read extraction so we don't miss reads from
    # narrower member pairs whose clips fall outside the widest span.
    padding = max(500, (env_right - env_left) // 2)
    reads = extract_split_reads(
        bam_path, chrom,
        start_search=env_left, end_search=env_right,
        padding=min(padding, env_left),
    )
    if not reads:
        print(f"[{pair_label}] No split reads found — skipping.")
        return None

    # Auto-estimate flank from primary alignment lengths if not explicitly provided
    if flank is None:
        flank = estimate_flank_from_primary_alignments(reads)
    else:
        print(f"  [flank] Using user-supplied flank={flank}")

    # SV type detection
    tol = max(50, padding // 10)
    left_clips_at_left_bp, right_clips_at_left_bp = 0, 0
    REF_OPS_SA = {"M", "D", "N", "=", "X"}
    for read in reads:
        if read.is_unmapped or not read.cigartuples:
            continue
        cigars = read.cigartuples
        if cigars[0][0] in CLIP_OPS and abs(read.reference_start - refined_left) < tol:
            left_clips_at_left_bp  += 1
        if cigars[-1][0] in CLIP_OPS and abs(read.reference_end - refined_left) < tol:
            right_clips_at_left_bp += 1
        if read.has_tag("SA"):
            for sa_entry in read.get_tag("SA").rstrip(";").split(";"):
                fields = sa_entry.split(",")
                if len(fields) < 4 or fields[0] != chrom:
                    continue
                sa_pos = int(fields[1])
                sa_cigar = fields[3]
                ops = re.findall(r'(\d+)([MIDNSHP=X])', sa_cigar)
                if not ops:
                    continue
                sa_ref_len = sum(int(l) for l, op in ops if op in REF_OPS_SA)
                sa_end = sa_pos + sa_ref_len
                if ops[0][1] in {"S", "H"} and abs(sa_pos - refined_left) < tol:
                    left_clips_at_left_bp  += 1
                if ops[-1][1] in {"S", "H"} and abs(sa_end - refined_left) < tol:
                    right_clips_at_left_bp += 1

    sv_type = "deletion" if right_clips_at_left_bp >= left_clips_at_left_bp else "duplication"
    print(f"[{pair_label}] sv_type={sv_type}  {chrom}:{refined_left}–{refined_right}")
    print(f"[{pair_label}] merged pairs: {cluster['pairs']}")
    print(f"[{pair_label}] envelope: {env_left}–{env_right}  widest span used for consensus")

    fasta = pysam.FastaFile(ref_path)
    consensus_result = build_split_consensus(
        reads, refined_left, refined_right, flank, fasta, chrom,
        sv_type=sv_type, bam_path=bam_path,
        show_pileup=show_pileup,
        show_homozygous_variants=show_homozygous_variants,
        pileup_save_dir=save_dir,
        pair_label=pair_label,
    )

    if consensus_result[0] is None:
        print(f"[{pair_label}] skipping — {chrom} not found in reference.")
        return None

    left_consensus, right_consensus, left_var_mask, right_var_mask, right_bp, ref_left, ref_right = consensus_result

    (mh_len, mh_seq_left, mh_seq_right, mh_mismatches, left_mh_start, left_mh_end, right_mh_start, right_mh_end) = find_mmej_microhomology(
        left_consensus, left_var_mask,
        right_consensus, right_var_mask,
        ref_left, ref_right,
        refined_left, right_bp,
        flank=flank, sv_type=sv_type, min_mh=min_mh,
        fasta=fasta, chrom=chrom,
        pair_label=pair_label, save_dir=save_dir,
    )

    if mh_len and (mh_seq_left or mh_seq_right):
        FLANK = 10
        fasta = pysam.FastaFile(ref_path)
        left_junc_left_mh   = _safe_fetch(fasta, chrom, left_mh_start - FLANK,  left_mh_start)
        left_junc_right_mh  = _safe_fetch(fasta, chrom, left_mh_end + 1,        left_mh_end + 1 + FLANK)
        right_junc_left_mh  = _safe_fetch(fasta, chrom, right_mh_start - FLANK, right_mh_start)
        right_junc_right_mh = _safe_fetch(fasta, chrom, right_mh_end + 1,       right_mh_end + 1 + FLANK)
        fasta.close()
        plot_microhomology(
            mh_seq_left=mh_seq_left,
            mh_seq_right=mh_seq_right,
            left_mh_genomic_start=left_mh_start,
            right_mh_genomic_start=right_mh_start,
            right_mh_genomic_end=right_mh_end,
            left_junc_left_mh=left_junc_left_mh,
            left_junc_right_mh=left_junc_right_mh,
            right_junc_left_mh=right_junc_left_mh,
            right_junc_right_mh=right_junc_right_mh,
            sv_type=sv_type,
            sample_id=pair_label,
            chrom=chrom,
            save_path=os.path.join(save_dir, f"mh_{pair_label}.png"),
        )

    return {
        "pair_label":             pair_label,
        "chrom":                  chrom,
        "sv_type":                sv_type,
        "refined_left":           refined_left,
        "refined_right":          right_bp,
        "refined_left_envelope":  env_left,
        "refined_right_envelope": env_right,
        "n_pairs_merged":         len(cluster["pairs"]),
        "left_mh_genomic_start":  left_mh_start,
        "left_mh_genomic_end":    left_mh_end,
        "right_mh_genomic_start": right_mh_start,
        "right_mh_genomic_end":   right_mh_end,
        "microhomology_length":         mh_len,
        "microhomology_sequence_left":  mh_seq_left,
        "microhomology_sequence_right": mh_seq_right,
        "microhomology_mismatches":     mh_mismatches,
        "supporting_reads":             len(reads),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 10. Batch entry point
# ─────────────────────────────────────────────────────────────────────────────

def batch_microhomology_search(
    bam_path,
    ref_path,
    flank=600,
    min_mh=2,
    show_pileup=False,
    min_clip_len=50,
    min_support=5,
    min_support_frac=0.05,
    cluster_dist=10,
    min_mapq=60,
    pair_strategy="split_read_bridge",
    max_dist=None,
    top_n=50,
    save_dir=None,
    show_homozygous_variants=False,
    # ── sc_adjustment version ─────────────────────────────────────────────
    sa_share_threshold=0.50,    # kept for legacy; not used by fast path
    bridge_tol=200,             # SA-bridge search tolerance (bp)
):
    """
    Batch MMEJ microhomology search with soft-clip pair clustering.

    Pipeline stages
    ---------------
    A+B+C+D  find_major_clip_pairs    - single-pass SA-read discovery,
                                        position binning, pair tallying,
                                        support filtering.
    E        cluster_clip_pairs_by_reads - read-centric Union-Find, no
                                           extra BAM fetch.
    F        analyze_sv_cluster (x N)  - consensus + MH search per cluster.
    """
    import os, sys, time

    class _Tee:
        """Writes every print() to both the original stdout and a log file."""
        def __init__(self, *streams):
            self.streams = streams
        def write(self, data):
            for s in self.streams:
                s.write(data)
        def flush(self):
            for s in self.streams:
                s.flush()

    def _banner(stage, msg):
        ts  = time.strftime("%H:%M:%S")
        sep = "=" * 64
        print(f"\n{sep}", flush=True)
        print(f"[{ts}]  {stage}  |  {msg}", flush=True)
        print(sep, flush=True)

    bam_name = os.path.splitext(os.path.basename(bam_path))[0]
    if save_dir is None:
        save_dir = os.path.join("results", bam_name)
    os.makedirs(save_dir, exist_ok=True)

    log_path = os.path.join(save_dir, f"log_{bam_name}.txt")
    log_file = open(log_path, "w")
    _orig_stdout = sys.stdout
    sys.stdout = _Tee(_orig_stdout, log_file)

    try:
        _banner("PIPELINE START", f"sv_analysis_with_sc_adjustment")
        print(f"  BAM      : {bam_path}", flush=True)
        print(f"  Ref      : {ref_path}", flush=True)
        print(f"  Out dir  : {save_dir}", flush=True)
        print(f"  Log file : {log_path}", flush=True)

        # ── Stages A-D: single-pass pair discovery ────────────────────────────────
        _banner("STAGES A-D", "Single-pass clip-pair discovery")
        t0 = time.time()
        pairs_df = find_major_clip_pairs(
            bam_path, ref_path,
            chrom            = None,
            min_clip_len     = min_clip_len,
            min_support      = min_support,
            min_support_frac = min_support_frac,
            cluster_dist     = cluster_dist,
            top_n            = top_n,
            min_mapq         = min_mapq,
            max_dist         = max_dist,
            save_dir         = save_dir,
        )
        print(f"\n  Stages A-D elapsed: {time.time()-t0:.1f}s", flush=True)

        if pairs_df.empty:
            print("  No clip pairs found. Exiting.")
            return pd.DataFrame()

        # ── Stage E: Union-Find clustering ───────────────────────────────────────
        _banner("STAGE E", "Read-centric Union-Find clustering")
        t0 = time.time()
        clusters = cluster_clip_pairs_by_reads(
            pairs_df,
            bam_path,
            cluster_dist = cluster_dist,
            min_mapq     = min_mapq,
            bridge_tol   = bridge_tol,
        )
        print(f"  Stage E elapsed: {time.time()-t0:.1f}s", flush=True)

        if not clusters:
            print("  No clusters found. Exiting.")
            return pd.DataFrame()

        print(f"\n  {len(clusters)} cluster(s) queued for MH analysis:", flush=True)
        for i, cl in enumerate(clusters, 1):
            print(
                f"    {i:>3}. {cl['chrom']}  "
                f"{cl['refined_left']:,} – {cl['refined_right']:,}  "
                f"(span {cl['refined_right']-cl['refined_left']:,} bp  "
                f"support={cl['support']}  {len(cl['pairs'])} pair(s))",
                flush=True,
            )

        # ── Stage F: per-cluster MH analysis ────────────────────────────────────
        _banner("STAGE F", f"MH analysis for {len(clusters)} cluster(s)")
        results = []
        t_stage_f = time.time()

        for i, cluster in enumerate(clusters, 1):
            label = f"{bam_name}_c{i:02d}_{cluster['chrom']}_{cluster['refined_left']}"
            t0 = time.time()
            print(f"\n  [{time.strftime('%H:%M:%S')}]  Cluster {i}/{len(clusters)}: "
                  f"{cluster['chrom']} {cluster['refined_left']:,} – {cluster['refined_right']:,}",
                  flush=True)
            print(f"    pairs   : {cluster['pairs']}", flush=True)
            print(f"    support : {cluster['support']}", flush=True)
            row = analyze_sv_cluster(
                bam_path, ref_path, cluster,
                flank                    = flank,
                min_mh                   = min_mh,
                show_pileup              = show_pileup,
                min_mapq                 = min_mapq,
                pair_label               = label,
                save_dir                 = save_dir,
                show_homozygous_variants = show_homozygous_variants,
            )
            print(f"    elapsed : {time.time()-t0:.1f}s", flush=True)
            if row:
                results.append(row)
            break ### temporary to let top 1 MH cluster through only

        print(f"\n  Stage F total elapsed: {time.time()-t_stage_f:.1f}s", flush=True)

        summary_df = pd.DataFrame(results)
        if not summary_df.empty:
            out = os.path.join(save_dir, f"batch_mh_results_{bam_name}.csv")
            summary_df.to_csv(out, index=False)
            _banner("PIPELINE COMPLETE", f"{len(results)} result(s) → {out}")
            display(summary_df)
        else:
            _banner("PIPELINE COMPLETE", "No MH results found.")
        return summary_df
    finally:
        sys.stdout = _orig_stdout
        log_file.close()

