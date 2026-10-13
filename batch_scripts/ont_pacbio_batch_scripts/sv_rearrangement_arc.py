"""
sv_rearrangement_arc.py
───────────────────────
SNP-node rearrangement graph + interactive HTML arc diagram for ONE
pre-defined breakpoint cluster.

Extracted from sv_analysis_junction_rearr.py (graph + visualisation half
only; function bodies unchanged).  Breakpoint discovery / clustering is NOT
done here — sv_analysis_with_sc_adjustment.py owns that, and
sv_analysis_sc_adjustment_arc.py feeds its clusters into these functions.

Functions (in call order):
  Utilities:       _safe_fetch, reverse_complement
  Split reads:     extract_split_reads           – SA-tagged reads near a cluster
  SNP-node graph:
    decompose_alignment_to_snp_nodes  – emit (perfect_run → mismatch) nodes
    is_germline_mismatch              – filter SNP/noise nodes
    scan_context_copies               – paralogous-copy scan of node context
    build_snp_segment_graph           – build DiGraph across all reads
  Visualisation:
    plot_rearrangement_arc_diagram    – interactive HTML arc diagram

Coordinates: all computation is 0-based (pysam), as in the original module.
"""

import re
import os
import pysam
pysam.set_verbosity(0)
from collections import defaultdict
import networkx as nx


# ─────────────────────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────────────────────

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


