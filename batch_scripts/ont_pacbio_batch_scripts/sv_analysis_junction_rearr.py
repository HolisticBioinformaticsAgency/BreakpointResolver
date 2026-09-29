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
    scan_context_copies               – paralogous-copy scan of node context
    build_snp_segment_graph           – build DiGraph across all reads
  Visualisation:
    plot_rearrangement_arc_diagram    – interactive HTML arc diagram
  Pipeline driver:
    run_pipeline                      – one-call four-stage pipeline
    main                              – argparse CLI entry point

Usage (notebook)
────────────────
    from sv_analysis_junction_rearr import run_pipeline
    results = run_pipeline(bam_path=BAM, ref_path=REF,
                           out_dir="rearrangement_arc", max_clusters=1)

Usage (CLI)
───────────
    python sv_analysis_junction_rearr.py --bam sample.bam --ref ref.fna \
        --out-dir rearrangement_arc --max-clusters 1

Patch notes (2026-08-19b)
──────────────────────────
  • cluster_softclip_pairs: added an explicit comment/docstring explanation
    of what a cluster's read-support actually represents.  A cluster merges
    one or more (chrom, left_pos, right_pos) breakpoint PAIRS (from
    pair_clip_positions) whose bridging-read sets overlap by more than
    share_threshold.  The cluster's "reads" field is the UNION of every
    merged pair's bridging-read set — that union count (len(cluster["reads"]))
    is what run_pipeline reports as "bridging read(s)" per cluster, i.e. the
    cluster-level read support.  This is distinct from a single pair's own
    support (len of _reads_supporting_pair for that ONE pair before any
    merging) — merging pairs can only ever grow or hold steady the reported
    support, never shrink it, because it's a set union, not a re-count.
    No functional change — documentation only.

Patch notes (2026-08-19)
─────────────────────────
  • Inversion signal now uses ACTUAL strand, not coordinate order.  Every
    segment (primary + each supplementary) recorded in build_snp_segment_
    graph's ordered_segments now carries its own read.is_reverse /
    supp_rec.is_reverse flag.  At each junction edge (the edge connecting
    the last node of one segment to the first node of the next), the graph
    now tracks a new edge attribute inversion_weight: the number of
    supporting reads for which the two segments' strands actually differed.
    This is separate from junction_weight, which only counts segment
    boundary crossings regardless of strand.
  • plot_rearrangement_arc_diagram colors an edge as the "inversion signal"
    (red) ONLY when inversion_weight > 0 for that edge — i.e. a genuine
    strand flip confirmed from the alignment records — instead of the old
    heuristic (comparing sorted reference-position rank of the two nodes,
    which never looked at is_reverse and could misclassify translocations
    or same-strand backward jumps as "inversions").  Tooltips show the
    inversion_weight / weight ratio; the legend text is updated to
    "Inversion junction (strand flip confirmed)" vs "Collinear edge".
  • run_pipeline per-cluster summary reports a distinct "Inversion-signal
    edges" count (edges with inversion_weight > 0), separate from the
    existing "Junction-crossing edges" count.

Patch notes (2026-08-18c)
──────────────────────────
  • realign_node_to_reference REMOVED entirely, along with _REALIGN_CACHE,
    the require_100 parameter (build_snp_segment_graph / run_pipeline / CLI
    --soft-anchor), and every node attribute it produced: confirmed,
    ref_offset, from_cache, ambiguous_anchor, n_anchors, n_verbatim_copies,
    verbatim_starts, realigned_start.  scan_context_copies (paralogous
    context-copy scan) is now the only per-node reference-check.
    sv_informative is redefined as (not context_paralogous) and
    (not germline) — a node is "informative" when its reference context is
    NOT found elsewhere in the search window (on either strand) and its
    mismatch is not explained by non-split-read (germline) evidence.
  • Single shared BAM handle: _reads_supporting_pair, count_bridging_reads,
    _fetch_supplementary_records and is_germline_mismatch all take an
    optional bam= parameter (an already-open pysam.AlignmentFile).  When
    supplied they read from it instead of opening (and re-reading the BAM
    index for) a fresh handle on every call.  pair_clip_positions,
    cluster_softclip_pairs and build_snp_segment_graph now open ONE handle
    per stage and pass it down through their O(n) / O(n²) inner loops.
  • _BRIDGE_CACHE: bridging-read sets from _reads_supporting_pair are now
    cached by (bam_path, chrom, min(pos_a,pos_b), max(pos_a,pos_b), tol,
    min_mapq).  pair_clip_positions and cluster_softclip_pairs both need the
    bridging set for the same breakpoint pairs; previously each recomputed
    it from scratch.
  • decompose_alignment_to_snp_nodes: the MD-tag fast path no longer walks
    aligned (M/=/X) blocks one base at a time.  It now consumes whole
    MD MATCH runs in a single step and only drops to per-base handling at
    an actual MISMATCH, so a long perfect-match run costs O(1) instead of
    O(run length).  Output is unchanged — pure performance change.  The
    non-MD fallback path is unchanged (still per-base).

Patch notes (2026-08-18)
────────────────────────
  • scan_context_copies: per-node paralogous-copy scan.  The node's reference
    context (ref_start − flank … ref_end + flank, default flank=50 → ~110+ bp)
    is searched across the split-read window on BOTH strands using a
    seed-and-verify strategy (exact seed hits via str.find, then full-context
    identity check against context_min_identity — set < 1.0 to tolerate
    diverged copies).  The self hit at the node's own locus is excluded;
    remaining hits are stored as ctx_n_copies / ctx_copy_sites and exposed in
    tooltips and per-cluster summaries.  The reference window string is
    cached module-level (_REF_WINDOW_CACHE) and shared across nodes in the
    same cluster, so the scan costs one FASTA fetch per cluster.
  • run_pipeline() incorporated — the four-stage Jupyter driver block now
    lives in the module, so notebooks pick up code changes on re-import
    without editing the cell.  main() adds an argparse CLI for HPC runs.
  • build_snp_segment_graph: min_node_support prunes low-support nodes.

Patch notes (2026-08-17)
────────────────────────
  • extract_split_reads: two-SA-entry cap removed (multi-split reads kept);
    secondary / supplementary / duplicate / QC-fail records excluded.
  • find_major_clip_positions: _sa_majority_mapq renamed to _sa_any_mapq to
    match its actual any-pass behaviour and made robust to malformed SA
    fields; local-depth estimation reuses the already-open BAM handle instead
    of reopening the file per candidate; duplicate / QC-fail reads excluded
    from both the scan and the depth estimate.
  • count_bridging_reads / _reads_supporting_pair: deduplicated into a single
    implementation (count delegates to the set-returning helper); SA
    positions converted from 1-based (SAM convention) to 0-based;
    supplementary / duplicate / QC-fail records excluded.
  • is_germline_mismatch: duplicate / QC-fail reads excluded.
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

    Patched: every SA entry is examined (the previous two-entry cap silently
    dropped the 3rd+ segments of multi-split reads), and secondary /
    supplementary / duplicate / QC-fail records are excluded.
    """
    split_reads = []
    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        for read in samfile.fetch(chrom, max(0, start_search - padding),
                                  end_search + padding):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate or read.is_qcfail
                    or not read.has_tag("SA")):
                continue
            for entry in read.get_tag("SA").rstrip(";").split(";"):
                fields = entry.split(",")
                if len(fields) < 4:
                    continue
                try:
                    sa_chrom, sa_pos = fields[0], int(fields[1])
                except ValueError:
                    continue
                if sa_chrom == chrom and (
                    abs(sa_pos - start_search) < padding
                    or abs(sa_pos - end_search) < padding
                ):
                    split_reads.append(read)
                    break
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

    Reads with low mapping quality, secondary/supplementary alignments,
    duplicates, QC-fail reads, or SA tags below *min_mapq* are excluded.
    Clip positions are clustered within *cluster_dist* bp and filtered by
    *min_support* / *min_support_frac* of local depth.  Results are written
    to *save_dir*/major_clip_positions.csv.

    Patched: the local-depth estimate reuses the single open BAM handle
    (previously the file was reopened for every candidate position), and the
    SA MAPQ gate is implemented by _sa_any_mapq (any-pass semantics).
    """
    REF_OPS  = {"M", "D", "N", "=", "X"}
    CLIP_STR = {"S", "H"}

    left_pool  = defaultdict(list)
    right_pool = defaultdict(list)
    samfile    = pysam.AlignmentFile(bam_path, "rb")
    fetch_iter = samfile.fetch(chrom) if chrom else samfile.fetch()

    def _region_depth(chrom, pos, window=50):
        """Primary-read depth around *pos*, reusing the open BAM handle."""
        count = 0
        for read in samfile.fetch(chrom, max(0, pos - window), pos + window):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate or read.is_qcfail):
                continue
            count += 1
        return count

    def _sa_any_mapq(read, min_mapq):
        """True when the read has no SA tag or at least one SA entry with
        MAPQ >= *min_mapq*.  (Formerly _sa_majority_mapq — the name implied a
        majority vote, but the logic has always been any-pass.)"""
        if not read.has_tag("SA"):
            return True
        entries = [e for e in read.get_tag("SA").rstrip(";").split(";") if e]
        if not entries:
            return True
        for sa in entries:
            fields = sa.split(",")
            if len(fields) < 6:
                continue
            try:
                if int(fields[4]) >= min_mapq:
                    return True
            except ValueError:
                continue
        return False

    seen = set()
    for read in fetch_iter:
        if read.is_unmapped or not read.cigartuples:
            continue
        if read.is_duplicate or read.is_qcfail:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.is_supplementary or read.is_secondary:
            continue
        if not _sa_any_mapq(read, min_mapq):
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
    # NOTE: samfile deliberately left open — _region_depth reuses it below.

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
        samfile.close()
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
            samfile.close()
            print("No clip positions passed the fractional-support filter.")
            return df

    samfile.close()

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

# Module-level cache: (bam_path, chrom, min_pos, max_pos, tol, min_mapq) →
# set of bridging read names.  pair_clip_positions and cluster_softclip_pairs
# both need the bridging-read set for the same breakpoint pairs; caching here
# means the second caller gets a free hit instead of re-scanning the BAM.
_BRIDGE_CACHE = {}


def count_bridging_reads(bam_path, chrom, pos_a, pos_b, tol=200, min_mapq=0, bam=None):
    """Count distinct read names that bridge both pos_a and pos_b on *chrom*
    using the split_read_bridge strategy:
      • the primary alignment must clip near one position, AND
      • the SA tag must map (with a clip) near the other position,
    within *tol* bp.

    Thin wrapper around _reads_supporting_pair so the bridging logic lives in
    exactly one place.

    Parameters
    ----------
    bam : optional already-open pysam.AlignmentFile.  When provided, reused
          instead of opening (and re-reading the index of) a fresh handle —
          pass one in when calling this in a loop (see pair_clip_positions).
    """
    return len(_reads_supporting_pair(
        bam_path, chrom, pos_a, pos_b, tol=tol, min_mapq=min_mapq, bam=bam
    ))


def _reads_supporting_pair(bam_path, chrom, pos_a, pos_b, tol=200, min_mapq=0, bam=None):
    """Return the *set* of read names bridging pos_a ↔ pos_b (used by
    cluster_softclip_pairs for Union-Find overlap scoring, and by
    count_bridging_reads / pair_clip_positions).

    This is the SINGLE-PAIR support set: the reads bridging exactly this one
    (pos_a, pos_b) breakpoint pair, before any Union-Find merging into a
    cluster.  See cluster_softclip_pairs for how multiple pairs' sets get
    combined into a cluster's overall read support.

    Patched: supplementary, duplicate and QC-fail records are excluded, and
    SA positions are converted from 1-based (SAM convention) to 0-based.
    Results are cached in _BRIDGE_CACHE (keyed by normalized position order)
    so pair_clip_positions and cluster_softclip_pairs never rescan the same
    pair twice.  An already-open *bam* handle can be supplied to avoid
    reopening the file on every call.
    """
    cache_key = (bam_path, chrom, min(pos_a, pos_b), max(pos_a, pos_b),
                 tol, min_mapq)
    if cache_key in _BRIDGE_CACHE:
        return set(_BRIDGE_CACHE[cache_key])

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
            try:
                sa_pos = int(fields[1]) - 1   # SA pos is 1-based
            except ValueError:
                continue
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

    def _scan(sf):
        fetch_start = max(0, min(pos_a, pos_b) - tol)
        fetch_end   = max(pos_a, pos_b) + tol
        for read in sf.fetch(chrom, fetch_start, fetch_end):
            if (read.is_unmapped or not read.cigartuples or read.is_secondary
                    or read.is_supplementary or read.is_duplicate
                    or read.is_qcfail):
                continue
            if read.mapping_quality < min_mapq:
                continue
            cigars    = read.cigartuples
            ref_start = read.reference_start
            if _clips_near(cigars, ref_start, pos_a) and _sa_near(read, pos_b):
                bridging.add(read.query_name)
            elif _clips_near(cigars, ref_start, pos_b) and _sa_near(read, pos_a):
                bridging.add(read.query_name)

    if bam is not None:
        _scan(bam)
    else:
        with pysam.AlignmentFile(bam_path, "rb") as sf:
            _scan(sf)

    _BRIDGE_CACHE[cache_key] = set(bridging)
    return bridging


def pair_clip_positions(clip_df, bam_path=None, max_dist=None, bridge_tol=200):
    """Build (chrom, left_pos, right_pos) pairs from find_major_clip_positions()
    output using the split_read_bridge strategy only.

    Each candidate pair is scored by the number of reads that have a clip near
    one position AND an SA alignment near the other.  Pairs are greedily
    selected in descending bridge-read count order; each clip position is
    used at most once.

    Patched: a single BAM handle is opened once and reused for every
    candidate pair in the O(n²) scoring loop below, instead of
    _reads_supporting_pair reopening the file per pair.

    Parameters
    ----------
    clip_df    : DataFrame from find_major_clip_positions()
    bam_path   : path to the BAM file (required)
    max_dist   : optional maximum span (bp) between paired positions
    bridge_tol : tolerance (bp) for clip / SA proximity check

    Returns
    -------
    list of (chrom, left_pos, right_pos) tuples.  Each tuple's own read
    support (before any cluster merging) is len(_reads_supporting_pair(...))
    for that pair — printed here as "bridge_reads=" during scoring, and
    recomputed for free from _BRIDGE_CACHE by cluster_softclip_pairs.
    """
    if bam_path is None:
        raise ValueError("pair_clip_positions requires bam_path for "
                         "split_read_bridge strategy.")

    idxs = list(clip_df.index)
    n    = len(idxs)
    candidates = []
    print(f" Scoring {n*(n-1)//2} candidate pairs by bridging read count ...")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
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
                    tol=bridge_tol, bam=bam,
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

    What a cluster's read-support actually means
    ──────────────────────────────────────────────
    Each input *pair* already has its own single-pair bridging-read set from
    _reads_supporting_pair (cached in _BRIDGE_CACHE, keyed by that pair's own
    positions).  When two or more pairs are merged into one cluster (because
    their bridging-read sets overlap by more than *share_threshold*), the
    cluster's "reads" field becomes the UNION of every merged pair's
    bridging-read set — NOT a re-scan, and never smaller than any single
    member pair's own support.  That union size, len(cluster["reads"]), is
    what run_pipeline reports and logs as the cluster's "bridging read(s)"
    count (n_bridge) — i.e. "cluster pair support" = the combined read
    evidence across every pair that got merged into this one cluster.  A
    cluster with pairs=[(chr1,100,900)] and one with
    pairs=[(chr1,100,900),(chr1,105,950)] merged together will report a
    reads count >= either pair's individual bridge_reads score printed by
    pair_clip_positions, since merging only ever adds read names to the set.

    Patched: a single BAM handle is opened once and reused across all pairs;
    _reads_supporting_pair also transparently hits _BRIDGE_CACHE for any pair
    already scored by pair_clip_positions, so this loop is typically free of
    BAM I/O entirely when called right after pair_clip_positions.

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
        pairs          – original (chrom, lp, rp) tuples merged into this
                         cluster (may be a single pair, or several that were
                         Union-Find-merged for sharing bridging reads)
        reads          – union of bridging read-name sets across all of
                         "pairs" above — this IS the cluster's read support
                         (see "What a cluster's read-support actually means"
                         above); downstream code uses len(reads) directly
                         (e.g. run_pipeline's n_bridge, printed as
                         "bridging read(s)")
        chrom          – chromosome
        refined_left   – min of all member left positions
        refined_right  – max of all member right positions
    """
    if not pairs:
        return []

    pair_reads = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, lp, rp in pairs:
            rset = _reads_supporting_pair(bam_path, chrom, lp, rp,
                                          tol=bridge_tol, min_mapq=min_mapq,
                                          bam=bam)
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
            merged_reads |= pair_reads[m]   # union, not re-scan — see docstring above

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

def _parse_md_tag(md_string):
    """Parse an MD tag string into a list of operations.

    Each element is one of:
        ("MATCH",    n)       – n consecutive reference-matching bases
        ("MISMATCH", ref_base) – single mismatch; ref_base is the reference allele
        ("DEL",      ref_bases) – deletion from reference (^ATG…)

    Examples
    --------
    "10A5"      → [("MATCH",10), ("MISMATCH","A"), ("MATCH",5)]
    "3^AC4"     → [("MATCH",3),  ("DEL","AC"),     ("MATCH",4)]
    "0C0T10"    → [("MISMATCH","C"), ("MISMATCH","T"), ("MATCH",10)]
    """
    ops  = []
    i    = 0
    s    = md_string
    while i < len(s):
        # ── match run (digits) ──────────────────────────────────────────
        m = re.match(r'(\d+)', s[i:])
        if m:
            n = int(m.group(1))
            if n > 0:
                ops.append(("MATCH", n))
            i += m.end()
            continue
        # ── deletion (^[ACGTN]+) ───────────────────────────────────────
        if s[i] == '^':
            i += 1
            mb = re.match(r'([ACGTNacgtn]+)', s[i:])
            if mb:
                ops.append(("DEL", mb.group(1).upper()))
                i += mb.end()
            continue
        # ── single-base mismatch ───────────────────────────────────────
        if s[i] in 'ACGTNacgtn':
            ops.append(("MISMATCH", s[i].upper()))
            i += 1
            continue
        i += 1   # skip unexpected character
    return ops


def decompose_alignment_to_snp_nodes(read, fasta, min_perfect_run=10):
    """Walk a single alignment and emit one node per (perfect_run → mismatch).

    Mismatch detection strategy
    ───────────────────────────
    Fast path  (preferred): parse the MD tag with _parse_md_tag.
        The MD tag encodes every mismatch and deletion explicitly, so no
        reference fetch is needed during the per-base walk.  The reference
        base at each mismatch comes directly from the MD operation, and the
        perfect-run sequence is reconstructed from the read query sequence
        (which is 100 % reference-identical over a confirmed match run).

        The MD fast path consumes whole MD MATCH runs in one step and only
        does per-base handling at an actual MISMATCH, so a long perfect-match
        run costs O(1) instead of O(run length).

    Fallback: if the MD tag is absent, fall back to fetching the reference
        via _safe_fetch and comparing seq[i] != ref_seq[i] base by base.
        This path is slower but produces identical results.

    CIGAR interpretation
    ────────────────────
    M / = : aligned columns (match or mismatch)
    X     : explicit mismatch (rare; treated same as M mismatch)
    I     : insertion into read — query_pos advances, ref_pos does not,
            run origin is NOT reset (insertions don't break the prefix anchor)
    D     : deletion from reference — ref_pos advances, run resets because
            a gap breaks the contiguous matching prefix
    S     : soft clip — query_pos advances, not part of any run
    H     : hard clip — no cursor movement

    Node dict fields
    ────────────────
    chrom          chromosome
    ref_start      0-based start of the perfect run on the reference
    ref_end        0-based position of the mismatch (= exclusive run end)
    mismatch_ref   reference base at the mismatch position (from MD or fasta)
    mismatch_read  read base at the mismatch position
    query_start    0-based query offset of the run start
    query_end      0-based query offset of the mismatch base
    perfect_run    the read sequence over the perfect run (identical to ref)
    """
    nodes = []
    if not read.cigartuples or read.query_sequence is None:
        return nodes

    seq       = read.query_sequence          # full query string
    ref_pos   = read.reference_start        # current reference cursor
    query_pos = 0                           # current query cursor

    run_ref_start   = ref_pos
    run_query_start = 0

    # ── Choose fast path (MD) vs fallback (fasta fetch) ──────────────────────
    use_md   = read.has_tag("MD")
    md_queue = []   # list of (op, value) remaining from MD parse

    if use_md:
        md_queue = _parse_md_tag(read.get_tag("MD"))
    md_idx     = 0    # index into md_queue
    md_op_used = 0    # bases consumed from the current MD op

    def _consume_md_del():
        """Advance past one DEL op in the MD queue (called on CIGAR D)."""
        nonlocal md_idx, md_op_used
        if md_idx < len(md_queue) and md_queue[md_idx][0] == "DEL":
            md_idx    += 1
            md_op_used = 0

    # ── Walk CIGAR ────────────────────────────────────────────────────────────
    for op, length in read.cigartuples:

        # ── Aligned columns: M, =, X ─────────────────────────────────────────
        if op in (0, 7, 8):
            if use_md:
                remaining   = length
                local_ref   = ref_pos
                local_query = query_pos
                while remaining > 0:
                    if md_idx >= len(md_queue):
                        # MD exhausted mid-block (shouldn't normally happen):
                        # treat the rest as an unresolved match run, same as
                        # the original per-base fallback behaviour.
                        local_ref   += remaining
                        local_query += remaining
                        remaining    = 0
                        break

                    md_op, md_val = md_queue[md_idx]

                    if md_op == "DEL":
                        # DEL ops belong to the CIGAR D handler, not here.
                        md_idx    += 1
                        md_op_used = 0
                        continue

                    if md_op == "MATCH":
                        avail = md_val - md_op_used
                        take  = min(avail, remaining)
                        local_ref   += take
                        local_query += take
                        remaining   -= take
                        md_op_used  += take
                        if md_op_used >= md_val:
                            md_idx    += 1
                            md_op_used = 0
                        continue

                    # md_op == "MISMATCH": exactly one base
                    read_base = seq[local_query].upper()
                    ref_base  = md_val
                    perfect_len = local_ref - run_ref_start
                    if perfect_len >= min_perfect_run:
                        perfect_run_seq = seq[run_query_start:local_query]
                        nodes.append({
                            "chrom":         read.reference_name,
                            "ref_start":     run_ref_start,
                            "ref_end":       local_ref,   # 0-based mismatch pos
                            "mismatch_ref":  ref_base,
                            "mismatch_read": read_base,
                            "query_start":   run_query_start,
                            "query_end":     local_query,
                            "perfect_run":   perfect_run_seq,
                        })
                    # Reset run to the base AFTER this mismatch
                    run_ref_start   = local_ref + 1
                    run_query_start = local_query + 1
                    local_ref   += 1
                    local_query += 1
                    remaining   -= 1
                    md_idx    += 1
                    md_op_used = 0

                ref_pos   = local_ref
                query_pos = local_query

            else:
                # Fallback: fetch reference lazily (one base at a time is
                # slow; batch the whole block and index in)
                ref_block = _safe_fetch(
                    fasta, read.reference_name, ref_pos, ref_pos + length
                )
                if len(ref_block) < length:
                    # Unresolvable fetch error — skip entire block
                    query_pos += length
                    ref_pos   += length
                    run_ref_start   = ref_pos
                    run_query_start = query_pos
                else:
                    for i in range(length):
                        read_base   = seq[query_pos + i].upper()
                        ref_base    = ref_block[i].upper()
                        is_mismatch = (read_base != ref_base)
                        if is_mismatch:
                            perfect_len = (ref_pos + i) - run_ref_start
                            if perfect_len >= min_perfect_run:
                                perfect_run_seq = seq[run_query_start : query_pos + i]
                                nodes.append({
                                    "chrom":         read.reference_name,
                                    "ref_start":     run_ref_start,
                                    "ref_end":       ref_pos + i,
                                    "mismatch_ref":  ref_base,
                                    "mismatch_read": read_base,
                                    "query_start":   run_query_start,
                                    "query_end":     query_pos + i,
                                    "perfect_run":   perfect_run_seq,
                                })
                            run_ref_start   = ref_pos + i + 1
                            run_query_start = query_pos + i + 1
                    ref_pos   += length
                    query_pos += length

        # ── Insertion: query advances, ref does not, run continues ────────────
        elif op == 1:
            # Insertions don't break the perfect-prefix anchor — the reference
            # bases on either side of the insertion still form a contiguous
            # match.  We simply skip the inserted query bases.
            query_pos += length

        # ── Deletion: ref advances, run resets ───────────────────────────────
        elif op == 2:
            if use_md:
                _consume_md_del()
            ref_pos         += length
            run_ref_start    = ref_pos
            run_query_start  = query_pos

        # ── Soft clip: query advances, not part of any alignment ─────────────
        elif op == 4:
            query_pos += length
            # Only reset the run if this is a TRAILING soft clip (ref_pos has
            # already advanced past the aligned bases).  A LEADING soft clip
            # (ref_pos == read.reference_start, i.e. no aligned bases yet)
            # must NOT reset the run — the clipped bases are the unmapped
            # overhang from the far side of a breakpoint, and the first
            # aligned base after the clip is a valid run origin at ref_pos.
            if ref_pos > read.reference_start:
                # trailing clip — alignment is done, reset is harmless
                run_ref_start   = ref_pos
                run_query_start = query_pos

        # ── Hard clip: no cursor movement ────────────────────────────────────
        elif op == 5:
            pass

        # ── Reference skip (N — intron): treat like deletion ─────────────────
        elif op == 3:
            ref_pos         += length
            run_ref_start    = ref_pos
            run_query_start  = query_pos

    return nodes


def is_germline_mismatch(
    bam_path,
    fasta,
    chrom,
    ref_pos,
    mismatch_read_base,
    window=5,
    min_non_split_support=3,
    min_non_split_frac=0.1,
    bam=None,
):
    """Check whether a mismatch at (chrom, ref_pos, mismatch_read_base) also
    appears in reads that are NOT split reads (i.e. have no SA tag).

    Both thresholds must be met for a mismatch to be called germline:
        non_split_alt_count >= min_non_split_support   AND
        non_split_frac      >= min_non_split_frac

    Patched: duplicate and QC-fail reads are excluded from the evidence count.
    An already-open *bam* handle can be supplied (see build_snp_segment_graph)
    to avoid reopening the BAM once per unique node.

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

    def _scan(sf):
        nonlocal alt_count, total
        for read in sf.fetch(chrom,
                             max(0, ref_pos - window),
                             ref_pos + window + 1):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate or read.is_qcfail):
                continue
            if read.has_tag("SA"):          # skip split reads
                continue
            if not read.cigartuples or read.query_sequence is None:
                continue
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if rpos == ref_pos:
                    total += 1
                    if (read.query_sequence[qpos].upper()
                            == mismatch_read_base.upper()):
                        alt_count += 1
                    break

    if bam is not None:
        _scan(bam)
    else:
        with pysam.AlignmentFile(bam_path, "rb") as sf:
            _scan(sf)

    frac       = alt_count / total if total > 0 else 0.0
    is_germline = (alt_count >= min_non_split_support) and (frac >= min_non_split_frac)

    return {
        "is_germline":         is_germline,
        "non_split_alt_count": alt_count,
        "non_split_total":     total,
        "non_split_frac":      round(frac, 4),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 7. Node context-copy scan (paralog / segmental-duplication check)
# ─────────────────────────────────────────────────────────────────────────────

# Module-level cache: (chrom, win_start, win_end) → reference window string.
# scan_context_copies searches this window per node; caching it means one
# FASTA fetch per cluster instead of one per node.
_REF_WINDOW_CACHE = {}


def _get_ref_window(fasta, chrom, win_start, win_end):
    """Fetch (and cache) the reference string for a search window."""
    key = (chrom, int(win_start), int(win_end))
    if key not in _REF_WINDOW_CACHE:
        # Windows are per-cluster; keep memory bounded by holding one.
        _REF_WINDOW_CACHE.clear()
        _REF_WINDOW_CACHE[key] = _safe_fetch(fasta, chrom, win_start, win_end)
    return _REF_WINDOW_CACHE[key]


def build_split_read_window(reads):
    """Compute the reference search window from a collection of split-read
    alignments (primary + supplementary segments of the same read set).

    The window is defined as:
        chrom       – chromosome of the primary alignment (most common)
        win_start   – min(reference_start) across all segments on that chrom
        win_end     – max(reference_end)   across all segments on that chrom

    This gives the tightest biologically-justified search space: the node's
    query sequence *must* anchor somewhere inside the span already covered by
    the split reads that produced it.

    Parameters
    ----------
    reads : iterable of pysam.AlignedSegment

    Returns
    -------
    dict  {chrom: str, win_start: int, win_end: int}
    or None if no mapped reads are found.
    """
    chrom_spans = defaultdict(lambda: [float("inf"), 0])

    for read in reads:
        if read.is_unmapped or read.reference_name is None:
            continue
        chrom = read.reference_name
        chrom_spans[chrom][0] = min(chrom_spans[chrom][0], read.reference_start)
        chrom_spans[chrom][1] = max(chrom_spans[chrom][1], read.reference_end or read.reference_start)

        # Also consume SA-tag segments so supplementary mates extend the window
        if read.has_tag("SA"):
            for entry in read.get_tag("SA").rstrip(";").split(";"):
                fields = entry.split(",")
                if len(fields) < 4:
                    continue
                sa_chrom = fields[0]
                try:
                    sa_pos = int(fields[1])
                except ValueError:
                    continue
                sa_cigar = fields[3]
                ref_ops  = {"M", "D", "N", "=", "X"}
                ops      = re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
                ref_len  = sum(int(l) for l, op in ops if op in ref_ops)
                chrom_spans[sa_chrom][0] = min(chrom_spans[sa_chrom][0], sa_pos)
                chrom_spans[sa_chrom][1] = max(chrom_spans[sa_chrom][1], sa_pos + ref_len)

    if not chrom_spans:
        return None

    # Pick the chromosome with the widest span (handles trans-chr edge cases
    # gracefully by still returning a sensible single-chrom window)
    best_chrom = max(chrom_spans, key=lambda c: chrom_spans[c][1] - chrom_spans[c][0])
    return {
        "chrom":     best_chrom,
        "win_start": int(chrom_spans[best_chrom][0]),
        "win_end":   int(chrom_spans[best_chrom][1]),
    }


def scan_context_copies(
    fasta,
    chrom,
    ref_start,
    ref_end,
    win_chrom,
    win_start,
    win_end,
    flank=50,
    seed_len=20,
    n_seeds=5,
    min_identity=1.0,
    search_revcomp=True,
    max_sites=50,
):
    """Search for copies of a node's reference context within the split-read
    window — the paralogous-region scan.

    The context is reference sequence from (ref_start − flank) to
    (ref_end + flank) on *chrom* (default ~110+ bp around the node).  Every
    occurrence of that context inside [win_start, win_end] on *win_chrom* is
    located with a seed-and-verify strategy:

      1.  *n_seeds* evenly spaced *seed_len*-mers from the context are
          located exactly via str.find (C speed) — on the forward strand and,
          if *search_revcomp*, on the reverse complement.
      2.  Each seed hit defines a candidate context start; the full context
          is compared base-by-base and kept when identity ≥ *min_identity*.
          With the default min_identity=1.0 this is an exact-copy scan;
          lower it (e.g. 0.95) to catch diverged segmental-duplication
          copies — seeds still prune the search space first.

    The hit at the node's own locus (forward strand, start == context_start)
    is classified as self and excluded from the copy count.  If the context
    cannot be found at its own locus (self_found=False) the node's reference
    coordinates are inconsistent with the FASTA — a red flag worth logging.

    Parameters
    ----------
    fasta          : open pysam.FastaFile
    chrom          : chromosome of the node (context is fetched from here)
    ref_start      : 0-based start of the node's perfect run
    ref_end        : 0-based position of the node's mismatch
    win_chrom      : chromosome of the search window (may differ for supp
                     segments on another contig)
    win_start, win_end : search window coordinates (build_split_read_window)
    flank          : bp of reference flank on each side (default 50)
    seed_len       : exact seed length for candidate generation
    n_seeds        : number of evenly spaced seeds per strand
    min_identity   : minimum full-context identity to accept a copy
    search_revcomp : also search the reverse complement (inverted copies)
    max_sites      : cap on recorded copy sites (count is always complete)

    Returns a dict:
        context_len      length of the searched context
        context_start    0-based reference start of the context on *chrom*
        context_end      0-based exclusive end
        self_found       True when the context matched its own locus
        n_copies         copies found EXCLUDING the self hit
        copy_sites       [{"start", "end", "strand", "identity"}, …]
        searched_revcomp bool
    """
    result = {
        "context_len":       0,
        "context_start":     None,
        "context_end":       None,
        "self_found":        False,
        "n_copies":          0,
        "copy_sites":        [],
        "searched_revcomp":  bool(search_revcomp),
    }

    # ref_end is the 0-based mismatch position — inclusive, hence +1.
    ctx_start = max(0, int(ref_start) - flank)
    ctx_end   = int(ref_end) + flank + 1
    context   = _safe_fetch(fasta, chrom, ctx_start, ctx_end)
    if not context:
        return result

    L = len(context)
    result.update(context_len=L, context_start=ctx_start, context_end=ctx_end)

    ref_window = _get_ref_window(fasta, win_chrom, win_start, win_end)
    if not ref_window or len(ref_window) < L:
        return result

    def _hits_for_strand(query_seq, strand):
        """Seed-and-verify scan of query_seq against ref_window."""
        if L >= seed_len:
            if n_seeds > 1:
                seed_offs = sorted({
                    round(i * (L - seed_len) / (n_seeds - 1))
                    for i in range(n_seeds)
                })
            else:
                seed_offs = [0]
        else:
            seed_offs = [0]

        candidates = set()
        for s_off in seed_offs:
            seed = query_seq[s_off : s_off + seed_len]
            pos  = 0
            while True:
                found = ref_window.find(seed, pos)
                if found < 0:
                    break
                cand = found - s_off
                if 0 <= cand <= len(ref_window) - L:
                    candidates.add(cand)
                pos = found + 1

        hits = []
        for cand in sorted(candidates):
            ref_slice = ref_window[cand : cand + L]
            identity  = sum(a == b for a, b in zip(query_seq, ref_slice)) / L
            if identity >= min_identity:
                hits.append({
                    "start":    win_start + cand,
                    "end":      win_start + cand + L,
                    "strand":   strand,
                    "identity": round(identity, 4),
                })
        return hits

    all_hits = _hits_for_strand(context, "+")
    if search_revcomp:
        all_hits += _hits_for_strand(reverse_complement(context), "-")

    copies = []
    for hit in all_hits:
        if hit["strand"] == "+" and hit["start"] == ctx_start:
            result["self_found"] = True
        else:
            copies.append(hit)

    result["n_copies"]   = len(copies)
    result["copy_sites"] = copies[:max_sites]
    return result


# ─────────────────────────────────────────────────────────────────────────────
# 8. Build SNP-segment rearrangement graph
# ─────────────────────────────────────────────────────────────────────────────

def _fetch_supplementary_records(read, bam_path, mapq_min=0, bam=None):
    """Re-fetch all supplementary alignments for a split read from the BAM.

    Parses the SA tag of *read* and fetches each supplementary record by
    querying the BAM at the reported SA coordinates.  Only records that match
    by query_name AND are flagged supplementary are returned.

    Parameters
    ----------
    read     : pysam.AlignedSegment  (the primary record carrying an SA tag)
    bam_path : str                   (path to the indexed BAM file)
    mapq_min : int                   (minimum MAPQ to accept a supp record)
    bam      : optional already-open pysam.AlignmentFile.  When provided,
               reused instead of opening (and re-reading the index of) a
               fresh handle per read — pass one in from build_snp_segment_graph.

    Returns
    -------
    list of pysam.AlignedSegment   (may be empty if none found / low MAPQ).
    Each returned record's .is_reverse flag reflects its OWN alignment
    strand — used by build_snp_segment_graph to detect real strand flips
    at junction edges (inversion signal).
    """
    if not read.has_tag("SA"):
        return []

    def _scan(handle):
        records = []
        for sa_entry in read.get_tag("SA").rstrip(";").split(";"):
            if not sa_entry:
                continue
            try:
                sa_chrom, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm = \
                    sa_entry.split(",")
                sa_pos  = int(sa_pos) - 1   # SA pos is 1-based
                sa_mapq = int(sa_mapq)
            except ValueError:
                continue

            if sa_mapq < mapq_min:
                continue

            # Fetch a small window around the reported SA start position
            # and match by query name + supplementary flag
            try:
                for candidate in handle.fetch(sa_chrom, max(0, sa_pos - 5),
                                              sa_pos + 5):
                    if (candidate.query_name == read.query_name
                            and candidate.is_supplementary
                            and candidate.reference_start == sa_pos
                            and candidate.query_sequence is not None):
                        records.append(candidate)
                        break
            except (ValueError, KeyError):
                continue   # contig not in BAM index
        return records

    if bam is not None:
        return _scan(bam)
    with pysam.AlignmentFile(bam_path, "rb") as handle:
        return _scan(handle)


def _add_nodes_from_segment(segment_nodes, all_nodes, read_label,
                             G, node_checked, fasta,
                             win_chrom, win_start, win_end,
                             bam_path, bam,
                             min_non_split_support, min_non_split_frac,
                             context_flank, context_min_identity,
                             context_search_revcomp):
    """Register nodes from one alignment segment into the graph.

    Shared by both the primary and every supplementary segment so the logic
    lives in one place.  Mutates G and node_checked in-place.

    Per unique node the function runs two validations:
      1. scan_context_copies      — paralogous copies of the ±flank reference
         context across the window, both strands (ctx_n_copies / ctx_copy_sites)
      2. is_germline_mismatch     — non-split-read alt support

    sv_informative = (not context_paralogous) and (not germline).

    Parameters
    ----------
    segment_nodes : list of node dicts from decompose_alignment_to_snp_nodes
    all_nodes     : the running ordered list for THIS READ (primary + supps
                    concatenated in query order) — used to build cross-segment
                    edges after all segments are processed
    read_label    : str  (for debug messages, e.g. "primary" / "supp:chr7:…")
    G, node_checked, fasta, win_chrom, win_start, win_end, bam_path,
    min_non_split_support, min_non_split_frac : as in caller
    bam           : already-open pysam.AlignmentFile, forwarded to
                    is_germline_mismatch to avoid reopening the BAM per node
    context_flank, context_min_identity, context_search_revcomp :
                    forwarded to scan_context_copies
    """
    for node in segment_nodes:
        nkey = (node["chrom"], node["ref_end"],
                node["mismatch_ref"], node["mismatch_read"])

        if nkey not in G:
            G.add_node(nkey, **node, support=0,
                       ctx_self_found=None, ctx_n_copies=0, ctx_copy_sites=None,
                       context_paralogous=False,
                       germline=None, germline_frac=None,
                       sv_informative=False, segment=read_label)

        if nkey not in node_checked:
            node_checked[nkey] = True

            ctx = scan_context_copies(
                fasta,
                node["chrom"], node["ref_start"], node["ref_end"],
                win_chrom, win_start, win_end,
                flank=context_flank,
                min_identity=context_min_identity,
                search_revcomp=context_search_revcomp,
            )
            G.nodes[nkey]["ctx_self_found"]     = ctx["self_found"]
            G.nodes[nkey]["ctx_n_copies"]       = ctx["n_copies"]
            G.nodes[nkey]["ctx_copy_sites"]     = ctx["copy_sites"]
            G.nodes[nkey]["context_paralogous"] = ctx["n_copies"] > 0

            germ = is_germline_mismatch(
                bam_path, fasta,
                node["chrom"], node["ref_end"],
                node["mismatch_read"],
                min_non_split_support=min_non_split_support,
                min_non_split_frac=min_non_split_frac,
                bam=bam,
            )
            G.nodes[nkey]["germline"]      = germ["is_germline"]
            G.nodes[nkey]["germline_frac"] = germ["non_split_frac"]

            G.nodes[nkey]["sv_informative"] = (
                not ctx["n_copies"] > 0
                and not germ["is_germline"]
            )

        G.nodes[nkey]["support"] += 1
        all_nodes.append((nkey, node))


def build_snp_segment_graph(
    reads,
    fasta,
    bam_path,
    min_perfect_run=10,
    min_non_split_support=3,
    min_non_split_frac=0.1,
    supp_mapq_min=0,
    min_node_support=1,
    context_flank=50,
    context_min_identity=1.0,
    context_search_revcomp=True,
):
    """Build a directed NetworkX graph where:
      • each node  = unique (chrom, ref_end, mismatch_ref, mismatch_read)
      • each edge  = two consecutive nodes within the same read (primary OR
                     supplementary), in query order, including edges that
                     cross the breakpoint junction from a primary segment
                     into a supplementary segment
      • edge weight = number of reads producing that consecutive node pair

    Supplementary alignments
    ────────────────────────
    For every read that carries an SA tag, all supplementary records are
    re-fetched from the BAM via _fetch_supplementary_records.  Their nodes
    are decomposed with decompose_alignment_to_snp_nodes and appended to the
    same read's ordered node list in query-coordinate order (sorted by
    query_start).  This means edges naturally cross the breakpoint:

        [ primary nodes … last_primary_node ]
              ↓  junction edge  (weight += 1 per read)
        [ first_supp_node … supp nodes … ]

    Inversion signal (real strand, not coordinate order)
    ──────────────────────────────────────────────────────
    Each segment (primary and every supplementary) carries its own strand
    (read.is_reverse / supp_rec.is_reverse).  At every junction edge (an
    edge whose two endpoint nodes come from different segments), the graph
    records inversion_weight: how many of the edge's supporting reads had
    DIFFERENT strands between the two segments they crossed.  This is a
    direct read-out of the alignment records' strand flags — it does not
    infer inversion from reference-coordinate order.  junction_weight still
    counts ALL segment-boundary crossings regardless of strand (so a
    junction with junction_weight=5, inversion_weight=0 is a real breakpoint
    that is NOT an inversion — e.g. a same-strand rearrangement).

    Patched (2026-08-19): inversion_weight added (see above).
    Patched (2026-08-18c): realign_node_to_reference removed; a single BAM
    handle is opened once for the whole read loop and passed down to
    _fetch_supplementary_records / is_germline_mismatch.

    Parameters
    ----------
    reads                 : iterable of pysam.AlignedSegment
    fasta                 : open pysam.FastaFile
    bam_path              : BAM file path
    min_perfect_run       : minimum perfect-match run length to emit a node
    min_non_split_support : germline filter absolute count threshold
    min_non_split_frac    : germline filter allele-fraction threshold
    supp_mapq_min         : minimum MAPQ to accept a supplementary record
    min_node_support      : drop nodes seen in fewer reads than this
                            (1 = keep everything; pruning removes their
                            incident edges too)
    context_flank         : bp of reference flank per side for the
                            paralogous-context scan (default 50)
    context_min_identity  : minimum identity for a context copy
                            (1.0 = exact copies only; lower to catch
                            diverged segmental duplications)
    context_search_revcomp: also scan the reverse complement (inverted copies)

    Returns
    -------
    nx.DiGraph  with node attributes:
        chrom, ref_start, ref_end, mismatch_ref, mismatch_read,
        query_start, query_end, perfect_run, segment,
        support, ctx_self_found, ctx_n_copies, ctx_copy_sites,
        context_paralogous, germline, germline_frac, sv_informative
    and edge attributes:
        weight (int), junction_weight (int), inversion_weight (int)
    """
    reads = list(reads)   # materialise so we can iterate twice

    # ── Build reference search window from primary + SA spans ─────────────
    window = build_split_read_window(reads)
    if window is None:
        print("[build_snp_segment_graph] No mapped reads — returning empty graph.")
        return nx.DiGraph()

    win_chrom = window["chrom"]
    win_start = window["win_start"]
    win_end   = window["win_end"]
    print(f"[build_snp_segment_graph] Search window: "
          f"{win_chrom}:{win_start:,}–{win_end:,} "
          f"({win_end - win_start:,} bp)  reads={len(reads)}")

    G            = nx.DiGraph()
    node_checked = {}

    n_supp_fetched   = 0
    n_junction_edges = 0
    n_inversion_edges = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in reads:
            if read.is_unmapped or read.query_sequence is None:
                continue

            # ── Decompose primary alignment ────────────────────────────────
            primary_nodes = decompose_alignment_to_snp_nodes(
                read, fasta, min_perfect_run
            )

            # ── Decompose each supplementary alignment ──────────────────────
            supp_segments = []   # (label, [nodes], query_start_of_segment, is_reverse)
            if read.has_tag("SA"):
                for supp_rec in _fetch_supplementary_records(
                        read, bam_path, mapq_min=supp_mapq_min, bam=bam):
                    supp_nodes = decompose_alignment_to_snp_nodes(
                        supp_rec, fasta, min_perfect_run
                    )
                    if supp_nodes:
                        # Use the query_start of the first node as sort key so
                        # we can order segments in query-coordinate order later
                        label = (f"supp:{supp_rec.reference_name}:"
                                 f"{supp_rec.reference_start}")
                        supp_segments.append((label, supp_nodes,
                                              supp_nodes[0]["query_start"],
                                              bool(supp_rec.is_reverse)))
                        n_supp_fetched += 1

            # ── Order all segments by query coordinate ──────────────────────
            # Primary query_start is the query_start of its first node (or 0
            # if it has no nodes — still place it first by reference_start).
            primary_qstart = (primary_nodes[0]["query_start"]
                              if primary_nodes else 0)
            ordered_segments = [
                ("primary", primary_nodes, primary_qstart, bool(read.is_reverse))
            ]
            ordered_segments += supp_segments
            ordered_segments.sort(key=lambda x: x[2])

            # ── Register all nodes and build the flat per-read node list ────
            all_read_nodes = []   # [(nkey, node_dict), …] in query order
            for label, seg_nodes, _, _ in ordered_segments:
                _add_nodes_from_segment(
                    seg_nodes, all_read_nodes, label,
                    G, node_checked, fasta,
                    win_chrom, win_start, win_end,
                    bam_path, bam,
                    min_non_split_support, min_non_split_frac,
                    context_flank, context_min_identity, context_search_revcomp,
                )

            # ── Add / increment edges along the flattened node list ────────
            for i in range(1, len(all_read_nodes)):
                pkey = all_read_nodes[i - 1][0]
                nkey = all_read_nodes[i][0]
                # Track whether this edge crosses a junction (segment
                # boundary) and, if so, whether the two segments' actual
                # alignment strands differ (real inversion signal).
                prev_idx = _segment_idx_for_node(i - 1, ordered_segments)
                curr_idx = _segment_idx_for_node(i,     ordered_segments)
                prev_label, _, _, prev_strand = ordered_segments[prev_idx]
                curr_label, _, _, curr_strand = ordered_segments[curr_idx]
                is_junction  = (prev_label != curr_label)
                is_inversion = is_junction and (prev_strand != curr_strand)

                if G.has_edge(pkey, nkey):
                    G[pkey][nkey]["weight"] += 1
                    if is_junction:
                        G[pkey][nkey]["junction_weight"] = \
                            G[pkey][nkey].get("junction_weight", 0) + 1
                    if is_inversion:
                        G[pkey][nkey]["inversion_weight"] = \
                            G[pkey][nkey].get("inversion_weight", 0) + 1
                else:
                    G.add_edge(pkey, nkey, weight=1,
                               junction_weight=1 if is_junction else 0,
                               inversion_weight=1 if is_inversion else 0)
                if is_junction:
                    n_junction_edges += 1
                if is_inversion:
                    n_inversion_edges += 1

    # ── Prune low-support nodes (removes their incident edges too) ───────
    n_pruned = 0
    if min_node_support > 1:
        low_support = [n for n, d in G.nodes(data=True)
                       if d.get("support", 0) < min_node_support]
        n_pruned = len(low_support)
        G.remove_nodes_from(low_support)

    print(f"[build_snp_segment_graph] Graph: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges  "
          f"(supp segments processed: {n_supp_fetched}, "
          f"junction edges: {n_junction_edges}, "
          f"inversion edges: {n_inversion_edges}, "
          f"pruned support<{min_node_support}: {n_pruned})")
    return G


def _segment_idx_for_node(flat_idx, ordered_segments):
    """Return which segment index a flat node index belongs to."""
    pos = 0
    for si, (_, seg_nodes, _, _) in enumerate(ordered_segments):
        pos += len(seg_nodes)
        if flat_idx < pos:
            return si
    return len(ordered_segments) - 1


# ─────────────────────────────────────────────────────────────────────────────
# 9. Interactive arc diagram
# ─────────────────────────────────────────────────────────────────────────────

def plot_rearrangement_arc_diagram(G, save_path="rearrangement_arc_diagram.html"):
    """Render an interactive SVG arc diagram of the SNP-segment rearrangement
    graph and save it to *save_path*.

    Visual encoding
    ───────────────
    Nodes   : circles scaled by support; teal glow = sv_informative = True
    Arcs    : an edge is drawn RED ("inversion junction") only when its
              inversion_weight > 0 — i.e. at least one supporting read
              crossed that junction between two segments with DIFFERENT
              actual alignment strands (read.is_reverse / supp_rec.is_reverse
              recorded in build_snp_segment_graph).  All other edges
              (including non-junction edges, and junction edges where the
              strands agree) are teal ("collinear").
              Arc height  ∝ genomic span between the two nodes.
              Arc opacity ∝ edge weight.
    Filters : minimum edge weight slider, sv-informative-only checkbox,
              show/hide germline nodes checkbox.
    Tooltips on hover for both nodes and arcs.  Edge tooltips show the
    inversion_weight / weight ratio; node tooltips include the paralogous
    context-copy count.
    """
    import json

    nodes_data = []
    for nkey, attr in G.nodes(data=True):
        nodes_data.append({
            "id":              str(nkey),
            "chrom":           attr.get("chrom", ""),
            "ref_end":         attr.get("ref_end", 0),
            "mismatch_ref":    attr.get("mismatch_ref", ""),
            "mismatch_read":   attr.get("mismatch_read", ""),
            "support":         attr.get("support", 0),
            "germline":        attr.get("germline", None),
            "germline_frac":   attr.get("germline_frac", None),
            "sv_informative":  attr.get("sv_informative", False),
            "ctx_copies":      attr.get("ctx_n_copies", 0),
            "ctx_sites":       (attr.get("ctx_copy_sites") or [])[:5],
        })
    nodes_data.sort(key=lambda n: (n["chrom"], n["ref_end"]))

    edges_data = []
    for u, v, attr in G.edges(data=True):
        edges_data.append({
            "source":           str(u),
            "target":           str(v),
            "weight":           attr.get("weight", 1),
            "junction_weight":  attr.get("junction_weight", 0),
            "inversion_weight": attr.get("inversion_weight", 0),
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
        '  <div class="leg-item"><div class="leg-circle" style="background:#6c7086"></div>Paralogous context / not informative</div>\n'
        '  <div class="leg-item"><div class="leg-swatch" style="background:#4f98a3"></div>Collinear edge</div>\n'
        '  <div class="leg-item"><div class="leg-swatch" style="background:#dd6974"></div>Inversion junction (strand flip confirmed)</div>\n'
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
        '    const isInversion=e.inversion_weight>0;\n'
        '    const span=Math.abs(x2-x1);\n'
        '    const arcH=clamp(span*0.45,20,ARC_AREA-20);\n'
        '    const cx=(x1+x2)/2,cy=nodeY-arcH;\n'
        '    const color=isInversion?"#dd6974":"#4f98a3";\n'
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
        '      tt.innerHTML=`<div><span class="tt-key">type</span> ${isInversion\n'
        '        ?\'<span class="tt-rev">inversion junction (strand flip confirmed)</span>\'\n'
        '        :\'<span class="tt-sv">collinear</span>\'}</div>\n'
        '        <div><span class="tt-key">weight</span> ${e.weight}</div>\n'
        '        <div><span class="tt-key">junction crossings</span> ${e.junction_weight}</div>\n'
        '        <div><span class="tt-key">strand-flip reads</span> ${e.inversion_weight} / ${e.weight}</div>\n'
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
        '        <div><span class="tt-key">germline</span> ${germStr}</div>\n'
        '        <div><span class="tt-key">sv_informative</span> ${svStr}</div>\n'
        '        <div><span class="tt-key">ctx copies</span> ${n.ctx_copies>0?n.ctx_copies+" @ "+n.ctx_sites.map(s=>s.start.toLocaleString()+"("+s.strand+")").join(", "):"0"}</div>`;\n'
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


# ─────────────────────────────────────────────────────────────────────────────
# 10. Pipeline driver  (the former Jupyter notebook block)
# ─────────────────────────────────────────────────────────────────────────────

def run_pipeline(
    bam_path,
    ref_path,
    sample_name,
    out_dir="rearrangement_arc",
    min_mapq=60,
    min_support=5,
    bridge_tol=200,
    share_threshold=0.50,
    padding=500,
    min_perfect_run=10,
    min_non_split_support=3,
    min_non_split_frac=0.1,
    supp_mapq_min=0,
    min_node_support=2,
    context_flank=50,
    context_min_identity=1.0,
    context_search_revcomp=True,
    max_clusters=None,
    show_df=True,
):
    """Run the full four-stage SV rearrangement pipeline end to end.

    Stages
    ──────
    1. find_major_clip_positions  – genome-wide soft-clip hotspot scan
    2. pair_clip_positions        – pair hotspots by split-read bridging
    3. cluster_softclip_pairs     – merge pairs sharing supporting reads
    4. build_snp_segment_graph + plot_rearrangement_arc_diagram per cluster

    Parameters
    ──────────
    bam_path, ref_path : input BAM (indexed) and reference FASTA (indexed)
    out_dir            : output directory for the hotspot CSV and arc diagrams
    min_mapq           : Stage 1 minimum mapping quality
    min_support        : Stage 1 minimum clipped-read support per hotspot
    bridge_tol         : Stage 2/3 clip↔SA proximity tolerance (bp)
    share_threshold    : Stage 3 minimum shared-read fraction for merging
                         (see cluster_softclip_pairs docstring for exactly
                         what "cluster read support" means once pairs merge)
    padding            : Stage 4 fetch padding around each cluster (bp)
    min_perfect_run    : minimum perfect anchor run for a SNP node.  Keep
                         ≥ ~8: shorter prefixes are noisier and less specific
    min_non_split_support / min_non_split_frac : germline filter thresholds
    supp_mapq_min      : minimum MAPQ for supplementary records
    min_node_support   : prune graph nodes seen in fewer reads than this
    context_flank      : bp of reference flank per side for the per-node
                         paralogous-context scan (default 50 → ~110 bp query)
    context_min_identity : minimum identity for a context copy
                         (1.0 = exact copies only; e.g. 0.95 catches diverged
                         segmental-duplication copies)
    context_search_revcomp : also scan the reverse complement for inverted copies
    max_clusters       : process only the first N clusters (None = all;
                         use 1 for debugging)
    show_df            : display() the hotspot DataFrame when in Jupyter,
                         plain print otherwise

    Returns
    ───────
    dict with keys:
        clip_df  – Stage 1 hotspot DataFrame
        pairs    – Stage 2 breakpoint pairs
        clusters – Stage 3 clusters (after max_clusters restriction).  Each
                   cluster's "reads" field is its combined pair support —
                   see cluster_softclip_pairs docstring.
        graphs   – list of {"cluster", "graph", "html"} per processed cluster
    """
    os.makedirs(out_dir, exist_ok=True)

    # ── Stage 1: find breakpoint hotspots ───────────────────────────────────
    print("=" * 60)
    print("Stage 1 — Major soft-clip positions")
    print("=" * 60)
    clip_df = find_major_clip_positions(
        bam_path,
        min_mapq=min_mapq,
        min_support=min_support,
        save_dir=out_dir,
    )
    print(f"  {len(clip_df)} clip position(s) found")
    if show_df:
        try:
            from IPython.display import display
            display(clip_df)
        except ImportError:
            print(clip_df)

    # ── Stage 2: pair hotspots via split-read bridging ──────────────────────
    print()
    print("=" * 60)
    print("Stage 2 — Pair clip positions (split_read_bridge)")
    print("=" * 60)
    pairs = pair_clip_positions(
        clip_df,
        bam_path=bam_path,
        bridge_tol=bridge_tol,
    )
    print(f"  {len(pairs)} pair(s) selected")
    for chrom, lp, rp in pairs:
        print(f"    {chrom}:{lp:,} ↔ {rp:,}  (span {rp - lp:,} bp)")

    # ── Stage 3: merge pairs sharing SA reads ───────────────────────────────
    print()
    print("=" * 60)
    print("Stage 3 — Cluster overlapping pairs")
    print("=" * 60)
    clusters = cluster_softclip_pairs(
        pairs,
        bam_path,
        share_threshold=share_threshold,
        bridge_tol=bridge_tol,
    )
    print(f"  {len(clusters)} cluster(s) after merging")

    # ── Optional restriction for debugging ──────────────────────────────────
    if max_clusters is not None:
        clusters = clusters[:max_clusters]
        if clusters:
            cl0 = clusters[0]
            print(f"  → Restricted to first {len(clusters)} cluster(s); first: "
                  f"{cl0['chrom']}:{cl0['refined_left']:,}–{cl0['refined_right']:,}  "
                  f"({len(cl0['reads'])} bridging read(s))")

    graphs = []

    if not clusters:
        print("  No clusters to process — check thresholds.")
    else:
        # ── Stage 4: SNP-node graph per cluster ────────────────────────────
        print()
        print("=" * 60)
        print("Stage 4 — Build SNP-segment rearrangement graphs")
        print("=" * 60)
        # Hygiene: clear per-run caches (keys are window/chrom-scoped so
        # cross-locus collisions cannot occur; clearing just bounds memory).
        _BRIDGE_CACHE.clear()
        _REF_WINDOW_CACHE.clear()

        fasta = pysam.FastaFile(ref_path)
        try:
            for ci, cl in enumerate(clusters):
                chrom    = cl["chrom"]
                lpos     = cl["refined_left"]
                rpos     = cl["refined_right"]
                # n_bridge = the cluster's combined read support: the union
                # of bridging reads across every pair merged into this
                # cluster by cluster_softclip_pairs (see that function's
                # docstring for exactly how pair-level support becomes
                # cluster-level support). This is what's printed below as
                # "bridging read(s)".
                n_bridge = len(cl["reads"])

                print(f"\n  Cluster {ci + 1}/{len(clusters)}: "
                      f"{chrom}:{lpos:,}–{rpos:,}  "
                      f"(span {rpos - lpos:,} bp, {n_bridge} bridging read(s))")

                reads = extract_split_reads(
                    bam_path, chrom, lpos, rpos, padding=padding
                )
                print(f"    extract_split_reads  → {len(reads)} split read(s)")

                if not reads:
                    print("    ⚠  No split reads extracted — skipping.")
                    continue

                n_with_sa = sum(1 for r in reads if r.has_tag("SA"))
                n_with_md = sum(1 for r in reads if r.has_tag("MD"))
                print(f"    SA tag coverage       → {n_with_sa}/{len(reads)}")
                print(f"    MD tag coverage       → {n_with_md}/{len(reads)} "
                      f"({'fast path' if n_with_md == len(reads) else 'mixed — some fallback fetches'})")

                G = build_snp_segment_graph(
                    reads, fasta, bam_path,
                    min_perfect_run=min_perfect_run,
                    min_non_split_support=min_non_split_support,
                    min_non_split_frac=min_non_split_frac,
                    supp_mapq_min=supp_mapq_min,
                    min_node_support=min_node_support,
                    context_flank=context_flank,
                    context_min_identity=context_min_identity,
                    context_search_revcomp=context_search_revcomp,
                )

                if G.number_of_nodes() == 0:
                    print("    ⚠  Empty graph — no (perfect_run → mismatch) nodes emitted.")
                    print("       Lower min_perfect_run cautiously (floor ~6–8).")
                    continue

                total_nodes    = G.number_of_nodes()
                sv_nodes       = [n for n, d in G.nodes(data=True) if d.get("sv_informative")]
                germ_nodes     = [n for n, d in G.nodes(data=True) if d.get("germline")]
                paralog_nodes  = [n for n, d in G.nodes(data=True)
                                  if d.get("context_paralogous")]
                no_self        = [n for n, d in G.nodes(data=True)
                                  if d.get("ctx_self_found") is False]
                supp_nodes     = [n for n, d in G.nodes(data=True)
                                  if d.get("segment", "").startswith("supp")]
                junction_edges  = [(u, v) for u, v, d in G.edges(data=True)
                                   if d.get("junction_weight", 0) > 0]
                inversion_edges = [(u, v) for u, v, d in G.edges(data=True)
                                   if d.get("inversion_weight", 0) > 0]

                print(f"    Graph                   → {total_nodes} node(s), {G.number_of_edges()} edge(s)")
                print(f"    SV-informative          → {len(sv_nodes)}")
                print(f"    Germline / noise        → {len(germ_nodes)}")
                print(f"    Paralogous context      → {len(paralog_nodes)}  "
                      f"(context copies in window, both strands)")
                if no_self:
                    print(f"    ⚠ Context self-miss      → {len(no_self)}  "
                          f"(reference/coordinate mismatch — check FASTA)")
                print(f"    From supp segments      → {len(supp_nodes)}")
                print(f"    Junction-crossing edges → {len(junction_edges)}")
                print(f"    Inversion-signal edges  → {len(inversion_edges)}  "
                      f"(confirmed via read.is_reverse strand flip)")

                out_html = os.path.join(out_dir, f"{sample_name}_arc_{chrom}_{lpos}_{rpos}.html")
                plot_rearrangement_arc_diagram(G, save_path=out_html)
                print(f"    Saved → {out_html}")

                graphs.append({"cluster": cl, "graph": G, "html": out_html})

        finally:
            fasta.close()

    print()
    print("=" * 60)
    print("Done.")
    print("=" * 60)

    return {
        "clip_df":  clip_df,
        "pairs":    pairs,
        "clusters": clusters,
        "graphs":   graphs,
    }


def main(argv=None):
    """argparse CLI entry point — runs run_pipeline() from the command line."""
    import argparse

    p = argparse.ArgumentParser(
        description="SNP-anchored SV rearrangement graph pipeline"
    )
    p.add_argument("--bam", required=True, help="Input BAM (indexed)")
    p.add_argument("--ref", required=True, help="Reference FASTA (indexed)")
    p.add_argument("--sample-name", required=True, help="Sample Name")
    p.add_argument("--out-dir", default="rearrangement_arc")
    p.add_argument("--min-mapq", type=int, default=60)
    p.add_argument("--min-support", type=int, default=5)
    p.add_argument("--bridge-tol", type=int, default=200)
    p.add_argument("--share-threshold", type=float, default=0.50)
    p.add_argument("--padding", type=int, default=500)
    p.add_argument("--min-perfect-run", type=int, default=10)
    p.add_argument("--min-node-support", type=int, default=1)
    p.add_argument("--supp-mapq-min", type=int, default=0)
    p.add_argument("--context-flank", type=int, default=50,
                   help="Reference flank per side for the context copy scan")
    p.add_argument("--context-min-identity", type=float, default=1.0,
                   help="Minimum identity for a context copy (1.0 = exact)")
    p.add_argument("--no-context-revcomp", action="store_true",
                   help="Do not scan the reverse complement for context copies")
    p.add_argument("--max-clusters", type=int, default=None,
                   help="Process only the first N clusters (default: all)")
    args = p.parse_args(argv)

    run_pipeline(
        bam_path=args.bam,
        ref_path=args.ref,
        sample_name=args.sample_name,
        out_dir=args.out_dir,
        min_mapq=args.min_mapq,
        min_support=args.min_support,
        bridge_tol=args.bridge_tol,
        share_threshold=args.share_threshold,
        padding=args.padding,
        min_perfect_run=args.min_perfect_run,
        supp_mapq_min=args.supp_mapq_min,
        min_node_support=args.min_node_support,
        context_flank=args.context_flank,
        context_min_identity=args.context_min_identity,
        context_search_revcomp=not args.no_context_revcomp,
        max_clusters=args.max_clusters,
        show_df=False,
    )


if __name__ == "__main__":
    main()
