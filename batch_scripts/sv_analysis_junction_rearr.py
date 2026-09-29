"""
sv_analysis_junction_rearr.py
─────────────────────────────
SNP-anchored SV rearrangement graph builder.

Functions (in call order):
  Utilities:       _safe_fetch, _has_adjacent_n, reverse_complement, get_ref_chroms
  Breakpoint discovery:
    extract_split_reads           – pull SA-tagged reads near breakpoints
    find_major_clip_positions     – genome-wide soft-clip hotspot finder
    count_bridging_reads          – split_read_bridge strategy only
    _reads_supporting_pair        – internal helper for cluster_softclip_pairs
    pair_clip_positions           – pair hotspots via split_read_bridge
    cluster_softclip_pairs        – Union-Find merge of overlapping pairs
  SNP-node graph:
    decompose_alignment_to_snp_nodes  – emit (perfect_run → mismatch) nodes
    is_germline_mismatch              – filter SNP/noise nodes
    realign_node_to_reference         – confirm reference anchor precision
    build_snp_segment_graph           – build DiGraph across all reads
  Visualisation:
    plot_rearrangement_arc_diagram    – interactive HTML arc diagram
"""

import re
import os
import pysam
pysam.set_verbosity(0)
import pandas as pd
from collections import Counter, defaultdict
import networkx as nx


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

    Clamps start and end to [0, contig_length]. Returns an empty string
    when the clamped window is zero-length or the contig is unknown.
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
    """Extract split reads (reads with an SA supplementary alignment tag) from
    a BAM file within a padded window around [start_search, end_search] on
    *chrom*.  Only reads whose SA partner maps near start_search or end_search
    are returned.
    """
    samfile = pysam.AlignmentFile(bam_path, "rb")
    split_reads = []
    for read in samfile.fetch(chrom, max(0, start_search - padding),
                              end_search + padding):
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
    """Scan a BAM file genome-wide (or a single *chrom*) and return a DataFrame
    of high-confidence soft-clip positions, each representing a candidate SV
    breakpoint.

    Reads with low mapping quality, secondary/supplementary alignments, or SA
    tags below *min_mapq* are excluded.  Clip positions are clustered within
    *cluster_dist* bp and filtered by *min_support* / *min_support_frac* of
    local depth.  Results are written to *save_dir*/major_clip_positions.csv.
    """
    REF_OPS  = {"M", "D", "N", "=", "X"}
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
                sa_chr   = fields[0]
                sa_pos   = int(fields[1])
                sa_cigar = fields[3]
                sa_mapq  = int(fields[4])
                if sa_mapq < min_mapq:
                    continue
                if sa_chr != ref_name:
                    continue
                ops = re.findall(r'(\d+)([MIDNSHP=X])', sa_cigar)
                if not ops:
                    continue
                ref_len = sum(int(l) for l, op in ops if op in REF_OPS)
                sa_end  = sa_pos + ref_len
                sa_left_len = int(ops[0][0]) if ops[0][1] in CLIP_STR else 0
                key_sa_left = (read.query_name, sa_pos, "left")
                if (ops[0][1] in CLIP_STR and sa_left_len >= min_clip_len
                        and key_sa_left not in seen):
                    seen.add(key_sa_left)
                    left_pool[ref_name].append((sa_pos, sa_left_len))
                sa_right_len = int(ops[-1][0]) if ops[-1][1] in CLIP_STR else 0
                key_sa_right = (read.query_name, sa_end, "right")
                if (ops[-1][1] in CLIP_STR and sa_right_len >= min_clip_len
                        and key_sa_right not in seen):
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
        for side, pool in (("left", left_pool[contig]),
                           ("right", right_pool[contig])):
            for row in _cluster_and_aggregate(pool, side):
                row["chrom"] = contig
                rows.append(row)

    if not rows:
        print("No major clip positions found with the given thresholds.")
        return pd.DataFrame()

    df = (
        pd.DataFrame(rows,
                     columns=["chrom", "position", "side", "support",
                               "mean_clip_len", "max_clip_len"])
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
            frac_thresh   = min_support_frac * total_depth
            effective_min = max(min_support, frac_thresh)
            passes        = row["support"] >= effective_min
            keep_mask.append(passes)
            if not passes:
                print(
                    f" [depth-filter] {row['chrom']}:{row['position']} "
                    f"support={row['support']} depth={total_depth} "
                    f"required≥{effective_min:.1f} → DROPPED"
                )
        df = df[keep_mask].reset_index(drop=True)
        if df.empty:
            print("No clip positions passed the fractional-support filter.")
            return df

    if top_n:
        df = df.head(top_n)

    _csv_out = (os.path.join(save_dir, "major_clip_positions.csv")
                if save_dir else "major_clip_positions.csv")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    df.to_csv(_csv_out, index=False)
    print(f"Saved: {_csv_out}")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 3. Bridging read counter  (split_read_bridge strategy)
# ─────────────────────────────────────────────────────────────────────────────

def count_bridging_reads(bam_path, chrom, pos_a, pos_b, tol=200, min_mapq=0):
    """Count distinct read names that bridge both pos_a and pos_b on *chrom*
    using the split_read_bridge strategy:
      • the primary alignment must clip near one position, AND
      • the SA tag must map (with a clip) near the other position,
    within *tol* bp.
    """
    REF_OPS  = {"M", "D", "N", "=", "X"}
    bridging = set()

    def _clips_near(cigar_ops, ref_start, target):
        if not cigar_ops:
            return False
        ref_end     = ref_start + sum(l for op, l in cigar_ops
                                      if op in (0, 2, 3, 7, 8))
        right_clip  = cigar_ops[-1][0] in (4, 5) and abs(ref_end   - target) <= tol
        left_clip   = cigar_ops[0][0]  in (4, 5) and abs(ref_start - target) <= tol
        return right_clip or left_clip

    def _sa_near(read, target):
        if not read.has_tag("SA"):
            return False
        for entry in read.get_tag("SA").rstrip(";").split(";"):
            fields = entry.split(",")
            if len(fields) < 4 or fields[0] != chrom:
                continue
            sa_pos   = int(fields[1])
            sa_cigar = fields[3]
            ops      = re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
            if not ops:
                continue
            ref_len = sum(int(l) for l, op in ops if op in REF_OPS)
            sa_end  = sa_pos + ref_len
            if ops[0][1]  in ("S", "H") and abs(sa_pos - target) <= tol:
                return True
            if ops[-1][1] in ("S", "H") and abs(sa_end - target) <= tol:
                return True
        return False

    with pysam.AlignmentFile(bam_path, "rb") as sf:
        fetch_start = max(0, min(pos_a, pos_b) - tol)
        fetch_end   = max(pos_a, pos_b) + tol
        for read in sf.fetch(chrom, fetch_start, fetch_end):
            if read.is_unmapped or not read.cigartuples or read.is_secondary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            cigars    = read.cigartuples
            ref_start = read.reference_start
            if _clips_near(cigars, ref_start, pos_a) and _sa_near(read, pos_b):
                bridging.add(read.query_name)
            elif _clips_near(cigars, ref_start, pos_b) and _sa_near(read, pos_a):
                bridging.add(read.query_name)
    return len(bridging)


def _reads_supporting_pair(bam_path, chrom, pos_a, pos_b, tol=200, min_mapq=0):
    """Return the *set* of read names bridging pos_a ↔ pos_b (used by
    cluster_softclip_pairs for Union-Find overlap scoring).
    Identical logic to count_bridging_reads but returns the set, not the count.
    """
    REF_OPS  = {"M", "D", "N", "=", "X"}
    bridging = set()

    def _clips_near(cigar_ops, ref_start, target):
        if not cigar_ops:
            return False
        ref_end = ref_start + sum(l for op, l in cigar_ops
                                  if op in (0, 2, 3, 7, 8))
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
            sa_pos   = int(fields[1])
            sa_cigar = fields[3]
            ops      = re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
            if not ops:
                continue
            ref_len = sum(int(l) for l, op in ops if op in REF_OPS)
            sa_end  = sa_pos + ref_len
            if ops[0][1]  in ("S", "H") and abs(sa_pos - target) <= tol:
                return True
            if ops[-1][1] in ("S", "H") and abs(sa_end - target) <= tol:
                return True
        return False

    with pysam.AlignmentFile(bam_path, "rb") as sf:
        fetch_start = max(0, min(pos_a, pos_b) - tol)
        fetch_end   = max(pos_a, pos_b) + tol
        for read in sf.fetch(chrom, fetch_start, fetch_end):
            if read.is_unmapped or not read.cigartuples or read.is_secondary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            cigars    = read.cigartuples
            ref_start = read.reference_start
            if _clips_near(cigars, ref_start, pos_a) and _sa_near(read, pos_b):
                bridging.add(read.query_name)
            elif _clips_near(cigars, ref_start, pos_b) and _sa_near(read, pos_a):
                bridging.add(read.query_name)
    return bridging


def pair_clip_positions(clip_df, bam_path=None, max_dist=None, bridge_tol=200):
    """Build (chrom, left_pos, right_pos) pairs from find_major_clip_positions()
    output using the split_read_bridge strategy only.

    Each candidate pair is scored by the number of reads that have a clip near
    one position AND an SA alignment near the other.  Pairs are greedily
    selected in descending bridge-read count order; each clip position is
    used at most once.

    Parameters
    ----------
    clip_df    : DataFrame from find_major_clip_positions()
    bam_path   : path to the BAM file (required)
    max_dist   : optional maximum span (bp) between paired positions
    bridge_tol : tolerance (bp) for clip / SA proximity check

    Returns
    -------
    list of (chrom, left_pos, right_pos) tuples
    """
    if bam_path is None:
        raise ValueError("pair_clip_positions requires bam_path for "
                         "split_read_bridge strategy.")

    idxs = list(clip_df.index)
    n    = len(idxs)
    candidates = []
    print(f" Scoring {n*(n-1)//2} candidate pairs by bridging read count ...")

    for ii in range(n):
        for jj in range(ii + 1, n):
            i, j     = idxs[ii], idxs[jj]
            row_a    = clip_df.loc[i]
            row_b    = clip_df.loc[j]
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
    used  = set()
    pairs = []
    for bc, i, j in candidates:
        if i in used or j in used:
            continue
        if bc == 0:
            continue
        used.add(i)
        used.add(j)
        row_a  = clip_df.loc[i]
        row_b  = clip_df.loc[j]
        left   = int(min(row_a["position"], row_b["position"]))
        right  = int(max(row_a["position"], row_b["position"]))
        pairs.append((row_a["chrom"], left, right))

    if not pairs:
        print(" [split_read_bridge] No pairs with bridging reads found.")
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
    """Merge breakpoint pairs that share > share_threshold of their SA-supporting
    reads (Jaccard over the smaller set).  Uses Union-Find.

    For each merged cluster the coordinate span is taken as the outer envelope:
    min(all lefts), max(all rights).

    Parameters
    ----------
    pairs           : list of (chrom, left_pos, right_pos) from pair_clip_positions()
    bam_path        : BAM file path
    share_threshold : minimum read-sharing fraction to merge two pairs
    bridge_tol      : tolerance (bp) used when fetching bridging reads
    min_mapq        : minimum MAPQ for bridging read fetch

    Returns
    -------
    list of dicts with keys:
        pairs          – original (chrom, lp, rp) tuples merged
        reads          – union of bridging read-name sets
        chrom          – chromosome
        refined_left   – min of all member left positions
        refined_right  – max of all member right positions
    """
    if not pairs:
        return []

    pair_reads = []
    for chrom, lp, rp in pairs:
        rset = _reads_supporting_pair(bam_path, chrom, lp, rp,
                                      tol=bridge_tol, min_mapq=min_mapq)
        pair_reads.append(rset)
        print(f" Pair {chrom}:{lp:,}–{rp:,}  {len(rset)} SA reads")

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
            si, sj = pair_reads[i], pair_reads[j]
            if not si or not sj:
                continue
            shared   = si & sj
            min_size = min(len(si), len(sj))
            frac     = len(shared) / min_size if min_size > 0 else 0.0
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
        chrom_k       = pairs[members[0]][0]
        all_left      = [pairs[m][1] for m in members]
        all_right     = [pairs[m][2] for m in members]
        merged_reads  = set()
        for m in members:
            merged_reads |= pair_reads[m]

        clusters.append({
            "pairs":         [pairs[m] for m in members],
            "reads":         merged_reads,
            "chrom":         chrom_k,
            "refined_left":  min(all_left),
            "refined_right": max(all_right),
        })

    print(f"\n{len(clusters)} cluster(s) from {n} pair(s).")
    return clusters


# ─────────────────────────────────────────────────────────────────────────────
# 5. SNP-node decomposition
# ─────────────────────────────────────────────────────────────────────────────

def decompose_alignment_to_snp_nodes(read, fasta, min_perfect_run=10):
    """Walk a single alignment base-by-base using CIGAR + reference comparison.

    Emits one node per (perfect_run → mismatch) unit.  A node is only emitted
    when the perfectly-aligned run preceding the mismatch is >= min_perfect_run
    bases long.

    Deletions reset the run origin because they break the contiguous matching
    prefix required for a valid node.

    Each node dict:
        chrom          chromosome
        ref_start      0-based start of the perfect run
        ref_end        0-based position of the mismatch (exclusive run end)
        mismatch_ref   reference base at the mismatch position
        mismatch_read  read base at the mismatch position
        query_start    0-based query offset of the run start
        query_end      0-based query offset of the mismatch
        perfect_run    the matching reference sequence of the run (for realign)
    """
    nodes            = []
    ref_pos          = read.reference_start
    query_pos        = 0
    run_ref_start    = ref_pos
    run_query_start  = 0

    if not read.cigartuples or read.query_sequence is None:
        return nodes

    for op, length in read.cigartuples:
        if op in (0, 7):        # M or = (consuming both ref and query)
            seq     = read.query_sequence[query_pos : query_pos + length]
            ref_seq = _safe_fetch(fasta, read.reference_name,
                                  ref_pos, ref_pos + length)
            if len(ref_seq) < length:
                ref_pos   += length
                query_pos += length
                continue
            for i in range(length):
                if seq[i].upper() != ref_seq[i].upper():   # mismatch
                    perfect_len = ref_pos + i - run_ref_start
                    if perfect_len >= min_perfect_run:
                        run_ref_seq = _safe_fetch(
                            fasta, read.reference_name,
                            run_ref_start, ref_pos + i
                        )
                        nodes.append({
                            "chrom":         read.reference_name,
                            "ref_start":     run_ref_start,
                            "ref_end":       ref_pos + i,   # exclusive
                            "mismatch_ref":  ref_seq[i].upper(),
                            "mismatch_read": seq[i].upper(),
                            "query_start":   run_query_start,
                            "query_end":     query_pos + i,
                            "perfect_run":   run_ref_seq,
                        })
                    # Reset run immediately after the mismatch
                    run_ref_start   = ref_pos + i + 1
                    run_query_start = query_pos + i + 1
            ref_pos   += length
            query_pos += length

        elif op == 1:   # I – insertion: skip query, keep ref cursor
            query_pos += length

        elif op == 2:   # D – deletion: advance ref, reset run (gap breaks prefix)
            ref_pos       += length
            run_ref_start  = ref_pos
            run_query_start = query_pos

        elif op == 4:   # S – soft clip
            query_pos += length

        elif op == 5:   # H – hard clip: no query bases consumed
            pass

    return nodes


# ─────────────────────────────────────────────────────────────────────────────
# 6. Germline mismatch filter
# ─────────────────────────────────────────────────────────────────────────────

def is_germline_mismatch(
    bam_path,
    fasta,
    chrom,
    ref_pos,
    mismatch_read_base,
    window=5,
    min_non_split_support=3,
    min_non_split_frac=0.1,
):
    """Check whether a mismatch at (chrom, ref_pos, mismatch_read_base) also
    appears in reads that are NOT split reads (i.e. have no SA tag).

    Both thresholds must be met for a mismatch to be called germline:
        non_split_alt_count >= min_non_split_support   AND
        non_split_frac      >= min_non_split_frac

    Returns a dict:
        is_germline          bool
        non_split_alt_count  reads without SA showing this alt base
        non_split_total      total non-split reads covering the position
        non_split_frac       float
    """
    ref_base = _safe_fetch(fasta, chrom, ref_pos, ref_pos + 1)
    if not ref_base:
        return {"is_germline": False, "non_split_alt_count": 0,
                "non_split_total": 0, "non_split_frac": 0.0}

    alt_count = 0
    total     = 0

    with pysam.AlignmentFile(bam_path, "rb") as sf:
        for read in sf.fetch(chrom,
                             max(0, ref_pos - window),
                             ref_pos + window + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.has_tag("SA"):          # skip split reads
                continue
            if not read.cigartuples:
                continue
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos == ref_pos:
                    total += 1
                    if (read.query_sequence[qpos].upper()
                            == mismatch_read_base.upper()):
                        alt_count += 1
                    break

    frac       = alt_count / total if total > 0 else 0.0
    is_germline = (alt_count >= min_non_split_support) and (frac >= min_non_split_frac)

    return {
        "is_germline":         is_germline,
        "non_split_alt_count": alt_count,
        "non_split_total":     total,
        "non_split_frac":      round(frac, 4),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 7. Node reference-anchor verification
# ─────────────────────────────────────────────────────────────────────────────

def realign_node_to_reference(
    fasta,
    chrom,
    node_query_seq,
    expected_ref_start,
    search_radius=500,
    min_identity=0.95,
    require_exact_prefix=True,
):
    """Slide node_query_seq across the reference in a window centred on
    expected_ref_start ± search_radius, scoring each position by the
    fraction of matching bases.

    When require_exact_prefix=True (default), only candidate positions where
    all bases *except* the last (the mismatch) match the reference are scored.
    This enforces the key invariant: the perfect run must anchor cleanly.

    Returns a dict:
        best_ref_start   0-based best matching position (or None)
        best_ref_end     best_ref_start + len(node_query_seq) (or None)
        best_identity    fraction of bases identical
        confirmed        identity >= min_identity
        original_start   the expected_ref_start passed in
        offset           best_ref_start - expected_ref_start (or None)
    """
    query    = node_query_seq.upper()
    qlen     = len(query)
    win_start = max(0, expected_ref_start - search_radius)
    win_end   = expected_ref_start + search_radius + qlen

    ref_window = _safe_fetch(fasta, chrom, win_start, win_end)
    if not ref_window or len(ref_window) < qlen:
        return {
            "best_ref_start": None, "best_ref_end": None,
            "best_identity":  0.0,  "confirmed":    False,
            "original_start": expected_ref_start,  "offset": None,
        }

    best_score  = -1
    best_offset = 0

    for i in range(len(ref_window) - qlen + 1):
        ref_slice  = ref_window[i : i + qlen]

        if require_exact_prefix:
            # The last base is the mismatch; the prefix must match perfectly
            prefix_len = qlen - 1
            if ref_slice[:prefix_len] != query[:prefix_len]:
                continue

        matches  = sum(a == b for a, b in zip(query, ref_slice))
        identity = matches / qlen
        if identity > best_score:
            best_score  = identity
            best_offset = i

    best_ref_start = win_start + best_offset
    confirmed      = best_score >= min_identity

    return {
        "best_ref_start": best_ref_start,
        "best_ref_end":   best_ref_start + qlen,
        "best_identity":  round(best_score, 4),
        "confirmed":      confirmed,
        "original_start": expected_ref_start,
        "offset":         best_ref_start - expected_ref_start,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 8. Build SNP-segment rearrangement graph
# ─────────────────────────────────────────────────────────────────────────────

def build_snp_segment_graph(
    reads,
    fasta,
    bam_path,
    min_perfect_run=10,
    search_radius=500,
    min_identity=0.95,
    min_non_split_support=3,
    min_non_split_frac=0.1,
):
    """Build a directed NetworkX graph where:
      • each node  = unique (chrom, ref_end, mismatch_ref, mismatch_read)
      • each edge  = two consecutive nodes in the same read, in query order
      • edge weight = number of reads producing that consecutive pair

    For every new node key encountered:
        1. realign_node_to_reference  → node["confirmed"]
        2. is_germline_mismatch       → node["germline"]
        3. node["sv_informative"]     = confirmed AND NOT germline

    node["support"] counts how many reads carry this node (regardless of
    sv_informative status).

    Parameters
    ----------
    reads              : iterable of pysam.AlignedSegment (primary + supplementary)
    fasta              : open pysam.FastaFile
    bam_path           : BAM file path (for is_germline_mismatch fetches)
    min_perfect_run    : minimum perfect-match run length to emit a node
    search_radius      : search window for realign_node_to_reference
    min_identity       : minimum identity for a confirmed anchor
    min_non_split_support / min_non_split_frac : germline thresholds

    Returns
    -------
    nx.DiGraph  with node attributes:
        chrom, ref_start, ref_end, mismatch_ref, mismatch_read,
        query_start, query_end, perfect_run,
        support, confirmed, ref_offset, germline, germline_frac,
        sv_informative
    and edge attributes:
        weight (int)
    """
    G            = nx.DiGraph()
    node_checked = {}   # nkey → True once realign + germline have run

    for read in reads:
        if read.is_unmapped or read.query_sequence is None:
            continue

        nodes = decompose_alignment_to_snp_nodes(read, fasta, min_perfect_run)

        for i, node in enumerate(nodes):
            nkey = (node["chrom"], node["ref_end"],
                    node["mismatch_ref"], node["mismatch_read"])

            if nkey not in G:
                G.add_node(nkey, **node, support=0,
                           confirmed=False, ref_offset=None,
                           germline=None, germline_frac=None,
                           sv_informative=False)

            # Run per-node checks exactly once
            if nkey not in node_checked:
                node_checked[nkey] = True

                # 1. Verify reference anchor
                anchor = realign_node_to_reference(
                    fasta,
                    node["chrom"],
                    node_query_seq=node["perfect_run"] + node["mismatch_read"],
                    expected_ref_start=node["ref_start"],
                    search_radius=search_radius,
                    min_identity=min_identity,
                )
                G.nodes[nkey]["confirmed"]   = anchor["confirmed"]
                G.nodes[nkey]["ref_offset"]  = anchor["offset"]

                # 2. Check for germline / noise origin
                germ = is_germline_mismatch(
                    bam_path, fasta,
                    node["chrom"], node["ref_end"],
                    node["mismatch_read"],
                    min_non_split_support=min_non_split_support,
                    min_non_split_frac=min_non_split_frac,
                )
                G.nodes[nkey]["germline"]      = germ["is_germline"]
                G.nodes[nkey]["germline_frac"] = germ["non_split_frac"]

                # 3. SV-informative flag
                G.nodes[nkey]["sv_informative"] = (
                    anchor["confirmed"] and not germ["is_germline"]
                )

            G.nodes[nkey]["support"] += 1

            # Add / increment edge from previous node
            if i > 0:
                prev_node = nodes[i - 1]
                pkey = (prev_node["chrom"], prev_node["ref_end"],
                        prev_node["mismatch_ref"], prev_node["mismatch_read"])
                if G.has_edge(pkey, nkey):
                    G[pkey][nkey]["weight"] += 1
                else:
                    G.add_edge(pkey, nkey, weight=1)

    return G


# ─────────────────────────────────────────────────────────────────────────────
# 9. Interactive arc diagram
# ─────────────────────────────────────────────────────────────────────────────

def plot_rearrangement_arc_diagram(G, save_path="rearrangement_arc_diagram.html"):
    """Render an interactive SVG arc diagram of the SNP-segment rearrangement
    graph and save it to *save_path*.

    Visual encoding
    ───────────────
    Nodes   : circles scaled by support; teal glow = sv_informative = True
    Arcs    : forward edges (left→right in ref coords) = teal;
              reverse edges (right→left) = red (potential inversion signal).
              Arc height  ∝ genomic span between the two nodes.
              Arc opacity ∝ edge weight.
    Filters : minimum edge weight slider, sv-informative-only checkbox,
              show/hide germline nodes checkbox.
    Tooltips on hover for both nodes and arcs.
    """
    import json

    nodes_data = []
    for nkey, attr in G.nodes(data=True):
        nodes_data.append({
            "id":            str(nkey),
            "chrom":         attr.get("chrom", ""),
            "ref_end":       attr.get("ref_end", 0),
            "mismatch_ref":  attr.get("mismatch_ref", ""),
            "mismatch_read": attr.get("mismatch_read", ""),
            "support":       attr.get("support", 0),
            "confirmed":     attr.get("confirmed", False),
            "germline":      attr.get("germline", None),
            "germline_frac": attr.get("germline_frac", None),
            "sv_informative":attr.get("sv_informative", False),
            "ref_offset":    attr.get("ref_offset", None),
        })
    nodes_data.sort(key=lambda n: (n["chrom"], n["ref_end"]))

    edges_data = []
    for u, v, attr in G.edges(data=True):
        edges_data.append({
            "source": str(u),
            "target": str(v),
            "weight": attr.get("weight", 1),
        })

    nodes_json = json.dumps(nodes_data)
    edges_json = json.dumps(edges_data)

    # NOTE: plain string — no f-string — to avoid escaping every JS { }
    html_template = (
        '<!DOCTYPE html>\n'
        '<html lang="en">\n'
        '<head>\n'
        '<meta charset="UTF-8">\n'
        '<title>SV Rearrangement Arc Diagram</title>\n'
        '<style>\n'
        ':root {\n'
        '  --bg:#0f1117;--surface:#1a1d27;--border:#2a2d3a;\n'
        '  --text:#cdd6f4;--muted:#6c7086;\n'
        '  --teal:#4f98a3;--red:#dd6974;--gold:#e8af34;\n'
        '}\n'
        '*{box-sizing:border-box;margin:0;padding:0;}\n'
        'body{background:var(--bg);color:var(--text);\n'
        '     font-family:"JetBrains Mono","Fira Code",monospace;\n'
        '     font-size:13px;overflow-x:hidden;}\n'
        'header{padding:16px 24px;background:var(--surface);\n'
        '       border-bottom:1px solid var(--border);\n'
        '       display:flex;align-items:center;gap:16px;flex-wrap:wrap;}\n'
        'header h1{font-size:15px;font-weight:600;letter-spacing:.05em;color:var(--teal);}\n'
        '.controls{display:flex;gap:20px;flex-wrap:wrap;align-items:center;\n'
        '          padding:10px 24px;background:var(--surface);\n'
        '          border-bottom:1px solid var(--border);}\n'
        '.ctrl-group{display:flex;align-items:center;gap:8px;}\n'
        'label{color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.06em;}\n'
        'input[type=range]{accent-color:var(--teal);width:120px;}\n'
        'input[type=checkbox]{accent-color:var(--teal);width:14px;height:14px;}\n'
        '#stats{color:var(--muted);font-size:11px;margin-left:auto;}\n'
        '#diagram-wrap{overflow-x:auto;padding:0;}\n'
        'svg{display:block;background:var(--bg);}\n'
        '.node-circle{cursor:pointer;transition:r .15s;}\n'
        '.node-circle:hover{stroke-width:2.5!important;}\n'
        '.node-label{pointer-events:none;font-size:9px;fill:var(--muted);font-family:monospace;}\n'
        '.arc-path{cursor:pointer;fill:none;transition:opacity .15s;}\n'
        '.arc-path:hover{opacity:1!important;stroke-width:3!important;}\n'
        '#tooltip{\n'
        '  position:fixed;display:none;pointer-events:none;\n'
        '  background:#1e2030;border:1px solid var(--border);\n'
        '  border-radius:8px;padding:10px 14px;font-size:11px;\n'
        '  color:var(--text);max-width:320px;z-index:999;\n'
        '  line-height:1.7;box-shadow:0 8px 32px rgba(0,0,0,.4);\n'
        '}\n'
        '#tooltip .tt-key{color:var(--muted);}\n'
        '#tooltip .tt-sv{color:var(--teal);font-weight:700;}\n'
        '#tooltip .tt-germ{color:var(--gold);}\n'
        '#tooltip .tt-rev{color:var(--red);}\n'
        '.legend{display:flex;gap:18px;align-items:center;padding:8px 24px;\n'
        '        background:var(--surface);border-bottom:1px solid var(--border);\n'
        '        flex-wrap:wrap;}\n'
        '.leg-item{display:flex;align-items:center;gap:6px;font-size:11px;color:var(--muted);}\n'
        '.leg-swatch{width:24px;height:3px;border-radius:2px;}\n'
        '.leg-circle{width:10px;height:10px;border-radius:50%;}\n'
        '</style>\n'
        '</head>\n'
        '<body>\n'
        '<header>\n'
        '  <svg width="28" height="28" viewBox="0 0 28 28" fill="none">\n'
        '    <circle cx="14" cy="14" r="13" stroke="#4f98a3" stroke-width="1.5"/>\n'
        '    <path d="M7 14 Q14 4 21 14" stroke="#4f98a3" stroke-width="1.5" fill="none"/>\n'
        '    <circle cx="7" cy="14" r="2.5" fill="#4f98a3"/>\n'
        '    <circle cx="21" cy="14" r="2.5" fill="#dd6974"/>\n'
        '  </svg>\n'
        '  <h1>SV Rearrangement Arc Diagram</h1>\n'
        '  <span id="stats"></span>\n'
        '</header>\n'
        '<div class="legend">\n'
        '  <div class="leg-item"><div class="leg-circle" style="background:#4f98a3;box-shadow:0 0 6px #4f98a3aa"></div>SV-informative node</div>\n'
        '  <div class="leg-item"><div class="leg-circle" style="background:#e8af34"></div>Germline / noise</div>\n'
        '  <div class="leg-item"><div class="leg-circle" style="background:#6c7086"></div>Unconfirmed anchor</div>\n'
        '  <div class="leg-item"><div class="leg-swatch" style="background:#4f98a3"></div>Forward edge</div>\n'
        '  <div class="leg-item"><div class="leg-swatch" style="background:#dd6974"></div>Reverse edge (inversion signal)</div>\n'
        '</div>\n'
        '<div class="controls">\n'
        '  <div class="ctrl-group">\n'
        '    <label for="wt-slider">Min edge weight</label>\n'
        '    <input type="range" id="wt-slider" min="1" max="50" value="1" step="1">\n'
        '    <span id="wt-val" style="color:var(--text);width:28px">1</span>\n'
        '  </div>\n'
        '  <div class="ctrl-group">\n'
        '    <input type="checkbox" id="sv-only">\n'
        '    <label for="sv-only">SV-informative nodes only</label>\n'
        '  </div>\n'
        '  <div class="ctrl-group">\n'
        '    <input type="checkbox" id="hide-germ">\n'
        '    <label for="hide-germ">Hide germline nodes</label>\n'
        '  </div>\n'
        '</div>\n'
        '<div id="diagram-wrap"><svg id="svg" xmlns="http://www.w3.org/2000/svg"></svg></div>\n'
        '<div id="tooltip"></div>\n'
        '<script>\n'
        'const NODES_RAW = ##NODES_JSON##;\n'
        'const EDGES_RAW = ##EDGES_JSON##;\n'
        '\n'
        'const PAD_LEFT=80,PAD_RIGHT=80,PAD_TOP=40,ARC_AREA=260,LABEL_H=60;\n'
        'const MIN_R=4,MAX_R=18,MIN_ARC_W=1,MAX_ARC_W=5;\n'
        'let minWt=1,svOnly=false,hideGerm=false;\n'
        'const svg=document.getElementById("svg");\n'
        'const tt=document.getElementById("tooltip");\n'
        '\n'
        'function clamp(v,lo,hi){return Math.max(lo,Math.min(hi,v));}\n'
        'function getNodeColor(n){\n'
        '  if(n.germline) return "#e8af34";\n'
        '  if(!n.confirmed) return "#6c7086";\n'
        '  if(n.sv_informative) return "#4f98a3";\n'
        '  return "#6c7086";\n'
        '}\n'
        '\n'
        'function render(){\n'
        '  svg.innerHTML="";\n'
        '  let nodes=NODES_RAW.filter(n=>{\n'
        '    if(svOnly && !n.sv_informative) return false;\n'
        '    if(hideGerm && n.germline) return false;\n'
        '    return true;\n'
        '  });\n'
        '  const nodeSet=new Set(nodes.map(n=>n.id));\n'
        '  let edges=EDGES_RAW.filter(e=>e.weight>=minWt&&nodeSet.has(e.source)&&nodeSet.has(e.target));\n'
        '  document.getElementById("stats").textContent=`${nodes.length} nodes \xb7 ${edges.length} edges`;\n'
        '  if(nodes.length===0){\n'
        '    svg.setAttribute("width",600);svg.setAttribute("height",200);\n'
        '    const t=document.createElementNS("http://www.w3.org/2000/svg","text");\n'
        '    t.setAttribute("x","50%");t.setAttribute("y","50%");\n'
        '    t.setAttribute("text-anchor","middle");t.setAttribute("fill","#6c7086");\n'
        '    t.textContent="No nodes match current filters";svg.appendChild(t);return;\n'
        '  }\n'
        '  const totalW=Math.max(900,nodes.length*60+PAD_LEFT+PAD_RIGHT);\n'
        '  const nodeY=PAD_TOP+ARC_AREA;\n'
        '  const totalH=nodeY+LABEL_H+PAD_TOP;\n'
        '  svg.setAttribute("width",totalW);svg.setAttribute("height",totalH);\n'
        '  const xStep=(totalW-PAD_LEFT-PAD_RIGHT)/Math.max(nodes.length-1,1);\n'
        '  const nodeX={};\n'
        '  nodes.forEach((n,i)=>{nodeX[n.id]=PAD_LEFT+i*xStep;});\n'
        '  const supps=nodes.map(n=>n.support);\n'
        '  const minS=Math.min(...supps),maxS=Math.max(...supps);\n'
        '  function nodeR(s){if(maxS===minS)return(MIN_R+MAX_R)/2;return MIN_R+(MAX_R-MIN_R)*(s-minS)/(maxS-minS);}\n'
        '  const wts=edges.map(e=>e.weight);\n'
        '  const minW=wts.length?Math.min(...wts):1,maxW=wts.length?Math.max(...wts):1;\n'
        '  function arcW(w){if(maxW===minW)return(MIN_ARC_W+MAX_ARC_W)/2;return MIN_ARC_W+(MAX_ARC_W-MIN_ARC_W)*(w-minW)/(maxW-minW);}\n'
        '  function arcOp(w){if(maxW===minW)return 0.6;return 0.25+0.65*(w-minW)/(maxW-minW);}\n'
        '\n'
        '  // Draw arcs first (nodes sit on top)\n'
        '  const arcLayer=document.createElementNS("http://www.w3.org/2000/svg","g");\n'
        '  svg.appendChild(arcLayer);\n'
        '  edges.forEach(e=>{\n'
        '    const x1=nodeX[e.source],x2=nodeX[e.target];\n'
        '    if(x1===undefined||x2===undefined) return;\n'
        '    const isForward=x2>=x1;\n'
        '    const span=Math.abs(x2-x1);\n'
        '    const arcH=clamp(span*0.45,20,ARC_AREA-20);\n'
        '    const cx=(x1+x2)/2,cy=nodeY-arcH;\n'
        '    const color=isForward?"#4f98a3":"#dd6974";\n'
        '    const sw=arcW(e.weight),op=arcOp(e.weight);\n'
        '    const path=document.createElementNS("http://www.w3.org/2000/svg","path");\n'
        '    path.setAttribute("d",`M ${x1},${nodeY} Q ${cx},${cy} ${x2},${nodeY}`);\n'
        '    path.setAttribute("stroke",color);\n'
        '    path.setAttribute("stroke-width",sw);\n'
        '    path.setAttribute("opacity",op);\n'
        '    path.classList.add("arc-path");\n'
        '    const srcN=NODES_RAW.find(n=>n.id===e.source)||{};\n'
        '    const tgtN=NODES_RAW.find(n=>n.id===e.target)||{};\n'
        '    path.addEventListener("mousemove",ev=>{\n'
        '      tt.style.display="block";\n'
        '      tt.style.left=(ev.clientX+14)+"px";\n'
        '      tt.style.top=(ev.clientY-10)+"px";\n'
        '      tt.innerHTML=`<div><span class="tt-key">type</span> ${isForward\n'
        '        ?\'<span class="tt-sv">forward</span>\'\n'
        '        :\'<span class="tt-rev">reverse (inversion signal)</span>\'}</div>\n'
        '        <div><span class="tt-key">weight</span> ${e.weight}</div>\n'
        '        <div><span class="tt-key">source</span> ${srcN.chrom}:${srcN.ref_end} (${srcN.mismatch_ref}\u2192${srcN.mismatch_read})</div>\n'
        '        <div><span class="tt-key">target</span> ${tgtN.chrom}:${tgtN.ref_end} (${tgtN.mismatch_ref}\u2192${tgtN.mismatch_read})</div>`;\n'
        '    });\n'
        '    path.addEventListener("mouseleave",()=>{tt.style.display="none";});\n'
        '    arcLayer.appendChild(path);\n'
        '  });\n'
        '\n'
        '  // Draw nodes\n'
        '  const nodeLayer=document.createElementNS("http://www.w3.org/2000/svg","g");\n'
        '  svg.appendChild(nodeLayer);\n'
        '  nodes.forEach(n=>{\n'
        '    const x=nodeX[n.id],r=nodeR(n.support),col=getNodeColor(n);\n'
        '    const g=document.createElementNS("http://www.w3.org/2000/svg","g");\n'
        '    if(n.sv_informative){\n'
        '      const glow=document.createElementNS("http://www.w3.org/2000/svg","circle");\n'
        '      glow.setAttribute("cx",x);glow.setAttribute("cy",nodeY);\n'
        '      glow.setAttribute("r",r+5);glow.setAttribute("fill","none");\n'
        '      glow.setAttribute("stroke",col);glow.setAttribute("stroke-width","1.5");\n'
        '      glow.setAttribute("opacity","0.35");g.appendChild(glow);\n'
        '    }\n'
        '    const c=document.createElementNS("http://www.w3.org/2000/svg","circle");\n'
        '    c.setAttribute("cx",x);c.setAttribute("cy",nodeY);\n'
        '    c.setAttribute("r",r);c.setAttribute("fill",col);\n'
        '    c.setAttribute("stroke","#0f1117");c.setAttribute("stroke-width","1");\n'
        '    c.classList.add("node-circle");\n'
        '    c.addEventListener("mousemove",ev=>{\n'
        '      tt.style.display="block";\n'
        '      tt.style.left=(ev.clientX+14)+"px";\n'
        '      tt.style.top=(ev.clientY-10)+"px";\n'
        '      const germStr=n.germline===null?"unknown"\n'
        '                   :n.germline?`<span class="tt-germ">YES (frac=${n.germline_frac})</span>`\n'
        '                   :"no";\n'
        '      const svStr=n.sv_informative?\'<span class="tt-sv">\u2713 SV-informative</span>\':"&mdash;";\n'
        '      tt.innerHTML=`<div><span class="tt-key">locus</span> ${n.chrom}:${n.ref_end}</div>\n'
        '        <div><span class="tt-key">mismatch</span> ref=${n.mismatch_ref} \u2192 read=${n.mismatch_read}</div>\n'
        '        <div><span class="tt-key">support</span> ${n.support} read(s)</div>\n'
        '        <div><span class="tt-key">confirmed</span> ${n.confirmed}</div>\n'
        '        <div><span class="tt-key">germline</span> ${germStr}</div>\n'
        '        <div><span class="tt-key">sv_informative</span> ${svStr}</div>\n'
        '        <div><span class="tt-key">ref_offset</span> ${n.ref_offset}</div>`;\n'
        '    });\n'
        '    c.addEventListener("mouseleave",()=>{tt.style.display="none";});\n'
        '    g.appendChild(c);\n'
        '    const lb=document.createElementNS("http://www.w3.org/2000/svg","text");\n'
        '    lb.setAttribute("x",x);lb.setAttribute("y",nodeY+r+14);\n'
        '    lb.setAttribute("text-anchor","middle");\n'
        '    lb.classList.add("node-label");\n'
        '    lb.textContent=n.ref_end.toLocaleString();\n'
        '    g.appendChild(lb);\n'
        '    const lb2=document.createElementNS("http://www.w3.org/2000/svg","text");\n'
        '    lb2.setAttribute("x",x);lb2.setAttribute("y",nodeY+r+26);\n'
        '    lb2.setAttribute("text-anchor","middle");\n'
        '    lb2.classList.add("node-label");\n'
        '    lb2.textContent=`${n.mismatch_ref}\u2192${n.mismatch_read}`;\n'
        '    g.appendChild(lb2);\n'
        '    nodeLayer.appendChild(g);\n'
        '  });\n'
        '}\n'
        '\n'
        'const wtSlider=document.getElementById("wt-slider");\n'
        'const wtVal=document.getElementById("wt-val");\n'
        'wtSlider.addEventListener("input",()=>{\n'
        '  minWt=+wtSlider.value;wtVal.textContent=minWt;render();\n'
        '});\n'
        'document.getElementById("sv-only").addEventListener("change",e=>{\n'
        '  svOnly=e.target.checked;render();\n'
        '});\n'
        'document.getElementById("hide-germ").addEventListener("change",e=>{\n'
        '  hideGerm=e.target.checked;render();\n'
        '});\n'
        'const maxEdgeW=EDGES_RAW.length?Math.max(...EDGES_RAW.map(e=>e.weight)):1;\n'
        'wtSlider.max=maxEdgeW;\n'
        'render();\n'
        '</script>\n'
        '</body>\n'
        '</html>\n'
    )

    html = html_template.replace("##NODES_JSON##", nodes_json).replace("##EDGES_JSON##", edges_json)

    os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True) if os.path.dirname(save_path) else None
    with open(save_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    print(f"Saved arc diagram: {save_path}  ({len(G.nodes())} nodes, {len(G.edges())} edges)")

