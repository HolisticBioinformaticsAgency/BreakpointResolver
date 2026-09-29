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
import os



CLIP_OPS = {4, 5}

def get_ref_chroms(ref_path):
    """Return the set of contig names present in the FASTA reference."""
    with pysam.FastaFile(ref_path) as fa:
        return set(fa.references)

def _has_adjacent_n(fasta, chrom, pos, window=5):
    """Return True if any base within *window* bp of *pos* in the reference is N."""
    if chrom not in fasta.references:
        return False
    region = fasta.fetch(chrom, max(0, pos - window), pos + window).upper()
    return "N" in region


# ─────────────────────────────────────────────────────────────────────────────
# 1. Extract split reads
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
                abs(sa_pos - start_search) < padding or
                abs(sa_pos - end_search)   < padding
            ):
                split_reads.append(read)
                break
    samfile.close()
    return split_reads



def reverse_complement(seq):
    return seq.translate(str.maketrans("ATCGatcg", "TAGCtagc"))[::-1]

# ─────────────────────────────────────────────────────────────────────────────
# Genome-wide major soft-clip / breakpoint finder
# ─────────────────────────────────────────────────────────────────────────────
def find_major_clip_positions(
    bam_path,
    chrom=None,
    min_clip_len=10,
    min_support=5,
    min_support_frac=0.05,    # ← new: drop clip sites with support < frac × local depth
    cluster_dist=10,
    top_n=50,
    require_sa=False,
    min_mapq=60,
    save_dir=None,
):
    REF_OPS  = {"M", "D", "N", "=", "X"}
    CLIP_STR = {"S", "H"}

    def _region_depth(chrom, pos, window=50):
        """Count primary mapped reads overlapping *pos* within *window* bp."""
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
        # Keep the read if ANY SA entry has a confident placement —
        # low-MAPQ SA entries (e.g. mapping into a repeat) are expected
        # for split reads spanning repeat-mediated SVs and are handled
        # by the per-entry MAPQ filter inside the SA pool loop.
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
                ops     = re.findall(r'(\d+)([MIDNSHP=X])', sa_cigar)
                if not ops:
                    continue
                ref_len = sum(int(l) for l, op in ops if op in REF_OPS)
                sa_end  = sa_pos + ref_len

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
        pd.DataFrame(rows, columns=["chrom", "position", "side", "support",
                                    "mean_clip_len", "max_clip_len"])
        .sort_values("support", ascending=False)
        .reset_index(drop=True)
    )
    # ── Fractional-support filter ─────────────────────────────────────────────
    # Require support >= max(min_support, min_support_frac × local_depth).
    # This prevents low-coverage artefacts from passing just because they clear
    # the absolute min_support floor.
    if min_support_frac and not df.empty:
        depth_cache = {}   # (chrom, pos) → total read depth
        keep_mask   = []

        for _, row in df.iterrows():
            key = (row["chrom"], row["position"])
            if key not in depth_cache:
                depth_cache[key] = _region_depth(row["chrom"], row["position"])

            total_depth  = depth_cache[key]
            frac_thresh  = min_support_frac * total_depth
            effective_min = max(min_support, frac_thresh)
            passes = row["support"] >= effective_min
            keep_mask.append(passes)

            if not passes:
                print(
                    f"  [depth-filter] {row['chrom']}:{row['position']}  "
                    f"support={row['support']}  depth={total_depth}  "
                    f"required≥{effective_min:.1f} ({min_support_frac*100:.0f}% of depth)  "
                    f"→ DROPPED"
                )

        df = df[keep_mask].reset_index(drop=True)

        if df.empty:
            print("No clip positions passed the fractional-support filter.")
            return df
    # ─────────────────────────────────────────────────────────────────────────
    if top_n:
        df = df.head(top_n)
        _csv_out = os.path.join(save_dir, "major_clip_positions.csv") if save_dir else "major_clip_positions.csv"
        os.makedirs(save_dir, exist_ok=True) if save_dir else None
        df.to_csv(_csv_out, index=False)
        print(f"Saved: {_csv_out}")
    return df


def get_homozygous_variants(bam_path, chrom, win_start, win_end, min_af=0.95, min_depth=5):
    """
    Return a dict of {genomic_pos: base} for positions where a non-ref base
    has allele frequency >= min_af across ALL reads (not just split reads).
    These are homozygous background variants and should NOT be treated as
    junction-derived or haplotype-specific variants.
    """
    samfile = pysam.AlignmentFile(bam_path, "rb")
    position_counts = {}  # {pos: Counter({base: count})}

    for read in samfile.fetch(chrom, win_start, win_end):
        if read.is_unmapped or read.query_sequence is None:
            continue
        ref_cursor   = read.reference_start
        query_cursor = 0
        if not read.cigartuples:
            continue
        for op, length in read.cigartuples:
            if op in (0, 7, 8):  # M/=/X
                for d in range(length):
                    rp = ref_cursor + d
                    if win_start <= rp < win_end:
                        base = read.query_sequence[query_cursor + d].upper()
                        if rp not in position_counts:
                            position_counts[rp] = Counter()
                        position_counts[rp][base] += 1
                ref_cursor   += length
                query_cursor += length
            elif op in (2, 3):  # D/N
                ref_cursor += length
            elif op == 1:       # I
                query_cursor += length
            elif op in (4, 5):  # S/H
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


# ─────────────────────────────────────────────────────────────────────────────
# 4. Build reference-anchored consensus
# ─────────────────────────────────────────────────────────────────────────────
def build_split_consensus(
    reads, refined_left, refined_right, flank,
    fasta, chrom, sv_type, bam_path, show_pileup=True,
    show_homozygous_variants=False, pileup_save_dir=None, pair_label=""
):
    if chrom not in fasta.references:
        print(f"WARNING: {chrom} not in reference — cannot build consensus, skipping.")
        return None, None, None, None, refined_right, "", "", {}
    
    is_deletion = sv_type == "deletion"

    # ── Coordinate adjustment ─────────────────────────────────────────────────
    # refined_right from find_major_clip_positions() is the pysam 0-based
    # reference_start of the SA alignment, which is one base to the right of
    # the true breakpoint base visible in IGV (1-based). Subtract 1 so that
    # the right window opens at the correct base.
    right_bp = refined_right  # ← ADJUSTED: shift window 1 base left

    if is_deletion:
        ref_left  = fasta.fetch(chrom, refined_left - flank, refined_left).upper()
        ref_right = fasta.fetch(chrom, right_bp,             right_bp + flank).upper()
    else:
        ref_left  = fasta.fetch(chrom, refined_left,          refined_left  + flank).upper()
        ref_right = fasta.fetch(chrom, refined_right - flank, refined_right        ).upper()

    def _cigar_walk_to_columns(read, win_start, win_end):
        bases = {}
        seq = read.query_sequence
        if seq is None:
            return bases
        ref_cursor   = read.reference_start
        query_cursor = 0
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                for d in range(length):
                    rp = ref_cursor + d
                    if win_start <= rp < win_end:
                        bases[rp - win_start] = seq[query_cursor + d].upper()
                ref_cursor   += length
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
        columns = [Counter() for _ in range(flank)]
        for col_dict in read_column_list:
            for col, base in col_dict.items():
                if 0 <= col < flank and base != "N":
                    columns[col][base] += 1

        result      = []
        low_cov_pos = []  # local column indices with insufficient depth

        for i, c in enumerate(columns):
            depth = sum(c.values())

            if depth < min_col_depth:
                result.append(ref_anchor[i])   # insufficient depth → fall back to ref
                if depth > 0:
                    low_cov_pos.append(i)      # partially covered (>0 but below threshold)
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
        out = []
        for i in range(flank):
            out.append(col_dict.get(i, "N"))
        return "".join(out)

    def _read_label(read):
        cigar_str = read.cigarstring or "none"
        sa_str    = read.get_tag("SA") if read.has_tag("SA") else "none"
        return f"{read.query_name}  CIGAR:{cigar_str}  SA:{sa_str}"

    def _build_sa_index(bam_path, read_names, chrom):
        """Build two lookup dicts for SA-based fallback column extraction.

        Returns
        -------
        supp_index : dict  (query_name, ref_start_1based) -> AlignedSegment
            Actual supplementary BAM records (flag 0x800). Present in ONT/
            minimap2 data. These have their own query_sequence and cigartuples
            and are used directly.
        primary_index : dict  query_name -> AlignedSegment
            Primary alignments for reads in read_names. Used as a fallback for
            SA-tag-only data (e.g. PacBio CCS / pbmm2) where no supplementary
            record is written. In this case we reconstruct the SA segment from
            the primary's query_sequence using the SA CIGAR + strand field.
        """
        supp_index    = {}
        primary_index = {}
        with pysam.AlignmentFile(bam_path, "rb") as sf:
            for read in sf.fetch(chrom):
                if read.query_name not in read_names:
                    continue
                if read.is_unmapped or not read.cigartuples:
                    continue
                if read.is_supplementary:
                    supp_index[(read.query_name, read.reference_start + 1)] = read  # +1: SA pos is 1-based
                elif not read.is_secondary:
                    primary_index[read.query_name] = read
        return supp_index, primary_index

    def _cigar_walk_to_columns_from_aln(aln, win_start, win_end):
        """Same logic as _cigar_walk_to_columns but accepts any AlignedSegment."""
        return _cigar_walk_to_columns(aln, win_start, win_end)

    def _cigar_walk_to_columns_sa_primary(primary_read, sa_pos, sa_cigar, sa_strand,
                                          win_start, win_end):
        """Reconstruct columns for a SA segment using the primary read's sequence.

        Used when no supplementary BAM record exists (SA-tag-only aligners such
        as pbmm2 for PacBio CCS). The SA tag encodes:
          sa_pos    – 0-based leftmost reference position of the SA alignment
          sa_cigar  – CIGAR string for the SA alignment
          sa_strand – '+' or '-'

        The primary read's query_sequence contains ALL bases (hard-clip is not
        applied on the primary). The SA segment's bases occupy the portion of
        query_sequence that is NOT consumed by the primary alignment — i.e. the
        soft/hard-clipped tail (or head) of the primary CIGAR.

        If sa_strand == '-', the SA segment is reverse-complemented relative to
        the primary sequence, so we reverse-complement the relevant slice before
        walking the SA CIGAR.
        """
        seq = primary_read.query_sequence
        if seq is None or not primary_read.cigartuples:
            return {}

        # ── Determine how many query bases the primary alignment consumes ──────
        primary_query_len = sum(
            length for op, length in primary_read.cigartuples
            if op in (0, 1, 4, 7, 8)   # M/I/S/=/X consume query
        )
        # Hard-clipped bases are absent from query_sequence entirely.
        # primary_read.query_length gives the length of query_sequence as stored.
        total_query_stored = len(seq)

        # Count query bases consumed by the primary CIGAR (M/I/S/=/X, not H)
        primary_consumes = sum(
            length for op, length in primary_read.cigartuples
            if op in (0, 1, 4, 7, 8)
        )

        # The SA segment occupies the bases NOT consumed by the primary.
        # For a deletion: primary aligns LEFT side → SA aligns RIGHT side.
        #   Primary CIGAR ends with S/H  → SA bases are at the END of query_sequence.
        # For a duplication: primary aligns RIGHT side → SA aligns LEFT side.
        #   Primary CIGAR starts with S/H → SA bases are at the START of query_sequence.
        primary_cigars = primary_read.cigartuples
        if primary_cigars[-1][0] in (4, 5):   # trailing clip → SA is at end
            sa_query_start = primary_consumes - (
                primary_cigars[-1][1] if primary_cigars[-1][0] == 4 else 0
            )
            # For hard clip on primary end, SA starts right after the aligned portion
            # Count only non-soft-clip consuming ops for the offset
            non_clip_consumes = sum(
                length for op, length in primary_cigars
                if op in (0, 1, 7, 8)   # M/I/=/X  (not S)
            )
            sa_query_start = non_clip_consumes
        else:                                  # leading clip → SA is at start
            sa_query_start = 0

        sa_ops = re.findall(r'(\d+)([MIDNSHP=X])', sa_cigar)
        if not sa_ops:
            return {}

        # SA query length (bases consumed by SA CIGAR incl. clips)
        sa_query_len = sum(
            int(l) for l, op in sa_ops if op in ('M', 'I', 'S', '=', 'X')
        )
        sa_query_end = sa_query_start + sa_query_len

        sa_seq_raw = seq[sa_query_start:sa_query_end]
        if not sa_seq_raw:
            return {}

        # Reverse-complement if SA is on the minus strand
        if sa_strand == '-':
            sa_seq_raw = reverse_complement(sa_seq_raw)

        # ── Walk SA CIGAR against sa_seq_raw ──────────────────────────────────
        bases       = {}
        ref_cursor  = sa_pos
        q_cursor    = 0
        op_start    = 0
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
                q_cursor   += length
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
        """Return list of (sa_pos, sa_cigar, sa_strand) for SA entries whose
        clip side lands within *tol* bases of *target_pos* on the same chrom."""
        if not read.has_tag("SA"):
            return []
        REF_OPS_SA = {'M', 'D', 'N', '=', 'X'}
        hits = []
        for sa_entry in read.get_tag("SA").rstrip(";").split(";"):
            fields = sa_entry.split(",")
            if len(fields) < 4 or fields[0] != read.reference_name:
                continue
            sa_pos    = int(fields[1])
            sa_strand = fields[2] if len(fields) > 2 else '+'
            sa_cigar  = fields[3]
            ops       = re.findall(r'(\d+)([MIDNSHP=X])', sa_cigar)
            if not ops:
                continue
            ref_len = sum(int(l) for l, op in ops if op in REF_OPS_SA)
            sa_end  = sa_pos + ref_len
            if clip_side == 'right' and ops[-1][1] in ('S', 'H'):
                if abs(sa_end - target_pos) < tol:
                    hits.append((sa_pos, sa_cigar, sa_strand))
            if clip_side == 'left' and ops[0][1] in ('S', 'H'):
                if abs(sa_pos - target_pos) < tol:
                    hits.append((sa_pos, sa_cigar, sa_strand))
        return hits


    def _primary_clips_near(read, target_pos, clip_side, tol=50):
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

    def _get_background_homozygous_sites(
        bam_path, chrom, win_start, win_end, ref_seq, min_af=0.95, min_depth=8
    ):
        """
        Return {local_col: (genomic_pos, alt_base, alt_count, depth, af)}
        for sites where ALL reads show a non-reference allele at high AF.
        These should be masked FALSE in the split-read variant mask.
        """
        out = {}
        with pysam.AlignmentFile(bam_path, "rb") as samfile:
            for pileup_col in samfile.pileup(
                chrom,
                win_start,
                win_end,
                truncate=True,
                stepper="all",
                min_base_quality=0,
                ignore_overlaps=False,
                ignore_orphans=False,
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
        right_win_start, right_win_end = right_bp,             right_bp + flank
    else:
        left_win_start,  left_win_end  = refined_left,          refined_left  + flank
        right_win_start, right_win_end = refined_right - flank, refined_right

    # ── Background homozygous variants from ALL reads ─────────────────────────
    left_background_hom = _get_background_homozygous_sites(
        bam_path, chrom, left_win_start, left_win_end, ref_left, min_af=0.95
    )
    right_background_hom = _get_background_homozygous_sites(
        bam_path, chrom, right_win_start, right_win_end, ref_right, min_af=0.95
    )

    if show_homozygous_variants:
        print(f"Left background homozygous sites:  {len(left_background_hom)}")
        print(f"Right background homozygous sites: {len(right_background_hom)}")

        if left_background_hom:
            print("Left homozygous background variants:")
            for col in sorted(left_background_hom):
                gpos, alt, alt_count, depth, af = left_background_hom[col]
                print(
                    f"  col={col:>4}  genomic={gpos + 1:>12}  "
                    f"ref={ref_left[col]}  alt={alt}  {alt_count}/{depth}  AF={af:.3f}"
                )

        if right_background_hom:
            print("Right homozygous background variants:")
            for col in sorted(right_background_hom):
                gpos, alt, alt_count, depth, af = right_background_hom[col]
                print(
                    f"  col={col:>4}  genomic={gpos + 1:>12}  "
                    f"ref={ref_right[col]}  alt={alt}  {alt_count}/{depth}  AF={af:.3f}"
                )

    left_seqs, right_seqs = [], []
    left_cols, right_cols = [], []

    # Build an index of supplementary alignments fetched directly from the BAM.
    # This is necessary because reconstructing SA bases from the primary read's
    # query_sequence using the SA CIGAR offset is unreliable for hard-clipped
    # aligners (e.g. minimap2 on ONT data) and strand-reversed SA segments.
    # Fetching the supplementary record directly gives us its own query_sequence
    # and cigartuples, which are always correct.
    read_names = {r.query_name for r in reads if not r.is_unmapped}
    supp_index, primary_index = _build_sa_index(bam_path, read_names, chrom)

    for read in reads:
        if read.is_unmapped or not read.cigartuples:
            continue
        label = _read_label(read)

        if is_deletion:
            # ── Left breakpoint: primary ends (right-clip) near refined_left ──
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
                        # No supplementary record (SA-tag-only aligner, e.g. pbmm2);
                        # reconstruct from primary query_sequence with strand correction.
                        primary = primary_index.get(read.query_name)
                        cd = _cigar_walk_to_columns_sa_primary(
                            primary, sa_pos, sa_cigar, sa_strand, left_win_start, left_win_end,
                        ) if primary is not None else {}
                    if cd:
                        left_cols.append(cd)
                        left_seqs.append((label, _columns_to_padded_seq(cd, flank)))
                        break

            # ── Right breakpoint: primary starts (left-clip) near right_bp ───
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
                        # No supplementary record (SA-tag-only aligner, e.g. pbmm2);
                        # reconstruct from primary query_sequence with strand correction.
                        primary = primary_index.get(read.query_name)
                        cd = _cigar_walk_to_columns_sa_primary(
                            primary, sa_pos, sa_cigar, sa_strand, right_win_start, right_win_end,
                        ) if primary is not None else {}
                    if cd:
                        right_cols.append(cd)
                        right_seqs.append((label, _columns_to_padded_seq(cd, flank)))
                        break

        else:
            # ── Left breakpoint: primary starts (left-clip) near refined_left ─
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
                        # No supplementary record (SA-tag-only aligner, e.g. pbmm2);
                        # reconstruct from primary query_sequence with strand correction.
                        primary = primary_index.get(read.query_name)
                        cd = _cigar_walk_to_columns_sa_primary(
                            primary, sa_pos, sa_cigar, sa_strand, left_win_start, left_win_end,
                        ) if primary is not None else {}
                    if cd:
                        left_cols.append(cd)
                        left_seqs.append((label, _columns_to_padded_seq(cd, flank)))
                        break

            # ── Right breakpoint: primary ends (right-clip) near refined_right ─
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
                        # No supplementary record (SA-tag-only aligner, e.g. pbmm2);
                        # reconstruct from primary query_sequence with strand correction.
                        primary = primary_index.get(read.query_name)
                        cd = _cigar_walk_to_columns_sa_primary(
                            primary, sa_pos, sa_cigar, sa_strand, right_win_start, right_win_end,
                        ) if primary is not None else {}
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

    print(f"Left  consensus ({len(left_seqs):3d} reads): {left_consensus}")
    print(f"Right consensus ({len(right_seqs):3d} reads): {right_consensus}")
    print(f"Left  reference: {ref_left}")
    print(f"Right reference: {ref_right}")

    left_variant_mask  = [lc != lr for lc, lr in zip(left_consensus,  ref_left)]
    right_variant_mask = [rc != rr for rc, rr in zip(right_consensus, ref_right)]

    # ── Mask homozygous background variants to FALSE ──────────────────────────
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
        print(f"Left  low-coverage positions (forced False): {left_low_cov}")
    if right_low_cov:
        print(f"Right low-coverage positions (forced False): {right_low_cov}")

    # print("left_variant_mask: ",  left_variant_mask)
    # print("right_variant_mask: ", right_variant_mask)

    if show_pileup:
        _print_pileup(
            "Left breakpoint pileup",  left_seqs,  left_consensus,  ref_left,  flank,
            after_bp=(not is_deletion),
        )
        _print_pileup(
            "Right breakpoint pileup", right_seqs, right_consensus, ref_right, flank,
            after_bp=is_deletion,
        )
    if pileup_save_dir:
        os.makedirs(pileup_save_dir, exist_ok=True)
        _lbl = f"_{pair_label}" if pair_label else ""
        _save_pileup_png(
            "Left breakpoint pileup",  left_seqs,  left_consensus,  ref_left,  flank,
            after_bp=(not is_deletion),
            save_path=os.path.join(pileup_save_dir, f"pileup_left{_lbl}.png"),
        )
        _save_pileup_png(
            "Right breakpoint pileup", right_seqs, right_consensus, ref_right, flank,
            after_bp=is_deletion,
            save_path=os.path.join(pileup_save_dir, f"pileup_right{_lbl}.png"),
        )

    return (
        left_consensus,
        right_consensus,
        left_variant_mask,
        right_variant_mask,
        right_bp,
        ref_left,
        ref_right,
    )





# ─────────────────────────────────────────────────────────────────────────────
# 5. Pileup display  (unchanged)
# ─────────────────────────────────────────────────────────────────────────────
def _print_pileup(title, labeled_seqs, consensus, ref_anchor, flank, after_bp, max_reads=100):
    if not labeled_seqs:
        display(HTML(f"<pre>=== {title} ===\nNo sequences.</pre>"))
        return

    width       = flank
    LABEL_WIDTH = 38

    def fmt_label(label):
        name = label.split("  ")[0]
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
    lines.append(f"=== {title} ===  ({len(seqs_up)} reads)")
    lines.append(f"{'ref'.ljust(LABEL_WIDTH)}  {ref_anchor[:width]}")
    lines.append(f"{' ' * LABEL_WIDTH}  {''.join(str((i // 10) % 10) for i in range(width))}")
    lines.append(f"{' ' * LABEL_WIDTH}  {''.join(str(i % 10)        for i in range(width))}")
    lines.append(f"{' ' * LABEL_WIDTH}  {'─' * width}")

    for label, s in seqs_up[:max_reads]:
        result, inserts = align_to_ref(s, ref_anchor[:width])
        annotated_parts = []
        for j, base in enumerate(result):
            ins = inserts.get(j)
            if ins:
                annotated_parts.append(
                    f'<span title="insertion: {ins}" style="color:blue;cursor:help">^</span>'
                )
            if base == "-":
                annotated_parts.append('<span style="color:grey">-</span>')
            elif base == " ":
                annotated_parts.append(" ")
            elif base != ref_anchor[j]:
                annotated_parts.append(f'<span style="color:red">{base.lower()}</span>')
            else:
                annotated_parts.append(base)
        lines.append(f"{fmt_label(label)}  {''.join(annotated_parts)}")

    if len(seqs_up) > max_reads:
        lines.append(f"  … ({len(seqs_up) - max_reads} more reads not shown)")

    lines.append(f"{' ' * LABEL_WIDTH}  {'─' * width}")
    lines.append(f"{'consensus'.ljust(LABEL_WIDTH)}  {consensus[:width]}")
    lines.append("\nColumn support:")

    all_results = [align_to_ref(s, ref_anchor[:width])[0] for _, s in seqs_up]
    for i in range(width):
        counts = Counter(
            r[i] for r in all_results
            if i < len(r) and r[i] not in (" ", "N")
        )
        total = sum(counts.values())
        major = counts.most_common(1)[0][0] if total else ref_anchor[i]
        frac  = counts.most_common(1)[0][1] / total if total else 0.0
        col_label = f"+{i:03d}" if after_bp else f"-{flank - 1 - i:03d}"
        lines.append(
            f"  {col_label}  ref={ref_anchor[i]}  cons={consensus[i] if i < len(consensus) else '?'}"
            f"  top={major}:{frac:.2f}  {dict(counts)}"
        )

    html_lines = []
    for line in lines:
        matched = False
        for label, _ in seqs_up[:max_reads]:
            truncated = fmt_label(label)
            if line.startswith(truncated.rstrip()):
                full_label = label.replace('"', '&quot;')
                rest = line[LABEL_WIDTH:]
                html_lines.append(
                    f'<span title="{full_label}" style="cursor:help; background:#fffbe6">'
                    f'{truncated}</span>{rest}'
                )
                matched = True
                break
        if not matched:
            html_lines.append(line)

    html  = '<pre style="font-family: monospace; font-size: 13px; line-height: 1.4; '
    html += 'background: #f8f8f8; padding: 10px; border: 1px solid #ddd; '
    html += 'overflow-x: auto; max-height: 600px; overflow-y: auto;">'
    html += "\n".join(html_lines) + "</pre>"
    display(HTML(html))





# ─────────────────────────────────────────────────────────────────────────────
# 5b. Pileup PNG renderer
# ─────────────────────────────────────────────────────────────────────────────
def _save_pileup_png(title, labeled_seqs, consensus, ref_anchor, flank, after_bp,
                     save_path, max_reads=100):
    """Render a pileup as a PNG file (mirrors _print_pileup layout)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colormaps

    BASE_COLOR = {"A": "#27ae60", "C": "#2980b9", "G": "#f39c12", "T": "#c0392b",
                  "-": "#7f8c8d"}

    width       = flank
    LABEL_WIDTH = 38
    seqs_up     = [(label, s.upper()[:width]) for label, s in labeled_seqs if s is not None]
    n_reads     = min(len(seqs_up), max_reads)

    # Extra rows: ref, ruler×2, divider, consensus, divider = 6
    n_rows = n_reads + 6
    row_h  = 0.18   # inches per row
    fig_w  = max(14, flank * 0.085 + 4)
    fig_h  = max(3.0, n_rows * row_h + 1.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.axis("off")

    # Use a fixed-width font
    FONT       = {"fontfamily": "DejaVu Sans Mono", "fontsize": 7}
    FONT_TITLE = {"fontfamily": "DejaVu Sans Mono", "fontsize": 8}
    line_ys = list(range(n_rows, 0, -1))   # top-down y positions

    def fmt_label(label):
        name = label.split("  ")[0]
        return name[:LABEL_WIDTH].ljust(LABEL_WIDTH)

    y_idx = 0

    def next_y():
        nonlocal y_idx
        y = line_ys[y_idx]
        y_idx += 1
        return y

    # Title
    ax.set_xlim(0, LABEL_WIDTH + width + 2)
    ax.set_ylim(0, n_rows + 2)
    ax.text(0, n_rows + 1.5, f"=== {title} ===  ({len(seqs_up)} reads)",
            **FONT_TITLE, fontweight="bold", va="top")

    y = next_y()
    ax.text(0, y, "ref".ljust(LABEL_WIDTH) + "  " + ref_anchor[:width], **FONT, va="center")

    # Ruler rows
    y = next_y()
    ruler1 = "".join(str((i // 10) % 10) for i in range(width))
    ax.text(0, y, " " * LABEL_WIDTH + "  " + ruler1, **FONT, va="center", color="#888888")
    y = next_y()
    ruler2 = "".join(str(i % 10) for i in range(width))
    ax.text(0, y, " " * LABEL_WIDTH + "  " + ruler2, **FONT, va="center", color="#888888")

    # Divider
    y = next_y()
    ax.text(0, y, " " * LABEL_WIDTH + "  " + "─" * width, **FONT, va="center", color="#cccccc")

    # Read rows
    for label, s in seqs_up[:max_reads]:
        y = next_y()
        lbl_str = fmt_label(label)
        ax.text(0, y, lbl_str + "  ", **FONT, va="center", color="#555555")
        x_off = LABEL_WIDTH + 2
        for col_i in range(width):
            if col_i >= len(s):
                break
            base = s[col_i].upper() if col_i < len(s) else " "
            ref_base = ref_anchor[col_i] if col_i < len(ref_anchor) else " "
            is_mismatch = base not in (" ", "N") and base != ref_base
            color = BASE_COLOR.get(base, "black") if is_mismatch else "#222222"
            ax.text(x_off + col_i, y, base, **FONT, va="center",
                    ha="left", color=color)

    # Divider
    y = next_y()
    ax.text(0, y, " " * LABEL_WIDTH + "  " + "─" * width, **FONT, va="center", color="#cccccc")

    # Consensus row
    y = next_y()
    ax.text(0, y, "consensus".ljust(LABEL_WIDTH) + "  " + consensus[:width], **FONT, va="center",
            fontweight="bold")

    plt.tight_layout(pad=0.2)
    plt.savefig(save_path, dpi=120, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved pileup PNG: {save_path}")

# ─────────────────────────────────────────────────────────────────────────────
# 6. MMEJ microhomology detection
# ─────────────────────────────────────────────────────────────────────────────
def find_mmej_microhomology(
    left_consensus,   left_variant_mask,
    right_consensus,  right_variant_mask,
    ref_left,         ref_right,
    refined_left,     refined_right,
    flank,
    sv_type,
    min_mh=2,
    edge_check=10,
    save_dir=None, pair_label=None,
):
    is_deletion = sv_type == "deletion"

    left_win_start  = (refined_left  - flank) if is_deletion else refined_left
    right_win_start = refined_right            if is_deletion else (refined_right - flank)

    def _strip_gaps(seq, mask, ref=None):
        new_seq, new_mask, new_ref = [], [], []
        for i, (b, m) in enumerate(zip(seq, mask)):
            if b != "-":
                new_seq.append(b)
                new_mask.append(m)
                if ref is not None:
                    new_ref.append(ref[i] if i < len(ref) else 'N')
        if ref is not None:
            return "".join(new_seq), new_mask, "".join(new_ref)
        return "".join(new_seq), new_mask

    left_consensus_aln,  left_variant_mask_aln,  ref_left_aln  = \
        _strip_gaps(left_consensus,  left_variant_mask,  ref_left)
    right_consensus_aln, right_variant_mask_aln, ref_right_aln = \
        _strip_gaps(right_consensus, right_variant_mask, ref_right)

    # ── Step 1: suffix-prefix block search on normalised sequences ────────────
    # Must run BEFORE the all-reference check so reclassification can happen first.
    def _normalise(seq, mask):
        return "".join(
            'X' if (i < len(mask) and mask[i]) else b
            for i, b in enumerate(seq)
        )

    norm_left  = _normalise(left_consensus_aln,  left_variant_mask_aln)
    norm_right = _normalise(right_consensus_aln, right_variant_mask_aln)

    EDGE_SEARCH = min(len(norm_left), len(norm_right), flank)

    if is_deletion:
        query        = norm_left[-EDGE_SEARCH:]
        ref          = norm_right[:EDGE_SEARCH]
        left_offset  = len(norm_left)  - EDGE_SEARCH
        right_offset = 0
    else:
        query        = norm_right[-EDGE_SEARCH:]
        ref          = norm_left[:EDGE_SEARCH]
        left_offset  = 0
        right_offset = len(norm_right) - EDGE_SEARCH

    best_score, best_overlap = -1, 0
    for overlap in range(min_mh, min(len(query), len(ref)) + 1):
        q_slice = query[-overlap:]
        r_slice = ref[:overlap]
        score   = sum(2 if q == r else -1 for q, r in zip(q_slice, r_slice))
        if score > best_score:
            best_score   = score
            best_overlap = overlap

    print(f"\n--- Suffix-prefix overlap search (EDGE_SEARCH={EDGE_SEARCH}) ---")
    print(f"  best_overlap={best_overlap}  best_score={best_score}")

    if best_overlap < min_mh or best_score <= 0:
        print("No matching block found.")
        return 0, "", None, None, None, None

    if is_deletion:
        best_l_start = left_offset  + (EDGE_SEARCH - best_overlap)
        best_l_end   = left_offset  + EDGE_SEARCH
        best_r_start = right_offset + 0
        best_r_end   = right_offset + best_overlap
    else:
        best_l_start = left_offset  + 0
        best_l_end   = left_offset  + best_overlap
        best_r_start = right_offset + (EDGE_SEARCH - best_overlap)
        best_r_end   = right_offset + EDGE_SEARCH

    best_k = best_overlap

    print(f"\n--- Best matching block (len={best_k}) ---")
    print(f"  Left  [{best_l_start}:{best_l_end}]  "
          f"genomic {left_win_start + best_l_start}–{left_win_start + best_l_end - 1}")
    print(f"  seq:  {left_consensus_aln[best_l_start:best_l_end]}")
    print(f"  Right [{best_r_start}:{best_r_end}]  "
          f"genomic {right_win_start + best_r_start}–{right_win_start + best_r_end - 1}")
    print(f"  seq:  {right_consensus_aln[best_r_start:best_r_end]}")

    # ── Step 2: reclassify variants within block against opposite reference ───
    # Do this BEFORE deciding whether all variants are False.
    left_variant_mask_aln  = list(left_variant_mask_aln)
    right_variant_mask_aln = list(right_variant_mask_aln)

    block_table = []

    print(f"\n--- Reclassifying variants within block against opposite reference ---")
    # print(f"  {'k':>4}  {'l_pos':>6}  {'r_pos':>6}  {'l_genomic':>12}  {'r_genomic':>12}  "
    #       f"{'l_cons':>6}  {'r_cons':>6}  {'ref_l':>6}  {'ref_r':>6}  "
    #       f"{'l_var':>6}  {'r_var':>6}  note")
    # print(f"  {'─'*4}  {'─'*6}  {'─'*6}  {'─'*12}  {'─'*12}  "
    #       f"{'─'*6}  {'─'*6}  {'─'*6}  {'─'*6}  {'─'*6}  {'─'*6}  {'─'*30}")

    for k in range(best_k):
        l_pos = best_l_start + k
        r_pos = best_r_start + k

        l_base    = left_consensus_aln[l_pos]
        r_base    = right_consensus_aln[r_pos]
        rl_base   = ref_left_aln[l_pos]  if l_pos < len(ref_left_aln)  else 'N'
        rr_base   = ref_right_aln[r_pos] if r_pos < len(ref_right_aln) else 'N'
        l_var_pre = left_variant_mask_aln[l_pos]
        r_var_pre = right_variant_mask_aln[r_pos]

        note_parts = []
        if l_var_pre:
            if l_base == rr_base:
                note_parts.append("L:junction")
            else:
                left_variant_mask_aln[l_pos] = False
                note_parts.append("L:hap-specific→False")
        if r_var_pre:
            if r_base == rl_base:
                note_parts.append("R:junction")
            else:
                right_variant_mask_aln[r_pos] = False
                note_parts.append("R:hap-specific→False")
        if not note_parts:
            note_parts.append("ref-match")

        block_table.append({
            "k":           k,
            "l_pos":       l_pos,
            "r_pos":       r_pos,
            "l_genomic":   left_win_start  + l_pos + 1,
            "r_genomic":   right_win_start + r_pos + 1,
            "l_cons":      l_base,
            "r_cons":      r_base,
            "ref_left":    rl_base,
            "ref_right":   rr_base,
            "l_var_final": left_variant_mask_aln[l_pos],
            "r_var_final": right_variant_mask_aln[r_pos],
            "note":        "  ".join(note_parts),
        })

        # print(f"  {k:>4}  {l_pos:>6}  {r_pos:>6}  "
        #       f"{left_win_start + l_pos + 1:>12}  {right_win_start + r_pos + 1:>12}  "
        #       f"{l_base:>6}  {r_base:>6}  {rl_base:>6}  {rr_base:>6}  "
        #       f"{str(left_variant_mask_aln[l_pos]):>6}  {str(right_variant_mask_aln[r_pos]):>6}  "
        #       f"{'  '.join(note_parts)}")

    block_df = pd.DataFrame(block_table)
    _block_csv = (
        os.path.join(save_dir, f"block_table_{pair_label}.csv")
        if save_dir and pair_label else
        os.path.join(save_dir, "block_table.csv")
        if save_dir else
        "block_table.csv"
    )
    block_df.to_csv(_block_csv, index=False)
    print(f"Saved block.csv ({len(block_df)} rows)")

    # ── Step 3: NOW check if all variants are False after reclassification ────
    left_vars_in_block  = left_variant_mask_aln[best_l_start : best_l_start + best_k]
    right_vars_in_block = right_variant_mask_aln[best_r_start : best_r_start + best_k]

    print("left variant mask in block: ",  left_vars_in_block)
    print("right variant mask in block: ", right_vars_in_block)

    if not any(left_vars_in_block) and not any(right_vars_in_block):
        print("\n--- All-reference case (post-reclassification): edge suffix-prefix search. Please assess the soft-clip edges to properly assess the breakpoint! ---")
        if is_deletion:
            left_tail  = left_consensus_aln[-edge_check:]
            right_head = right_consensus_aln[:edge_check]
            edge_mh_len = 0
            for k in range(min(len(left_tail), len(right_head)), 0, -1):
                if left_tail[-k:] == right_head[:k]:
                    edge_mh_len = k
                    break
            mh_seq = left_tail[-edge_mh_len:]
            left_mh_genomic_start  = refined_left  - edge_mh_len
            left_mh_genomic_end    = refined_left  - 1
            right_mh_genomic_start = refined_right
            right_mh_genomic_end   = refined_right + edge_mh_len - 1
        else:
            right_tail = right_consensus_aln[-edge_check:]
            left_head  = left_consensus_aln[:edge_check]
            edge_mh_len = 0
            for k in range(min(len(right_tail), len(left_head)), 0, -1):
                if right_tail[-k:] == left_head[:k]:
                    edge_mh_len = k
                    break
            mh_seq = right_tail[-edge_mh_len:]
            right_mh_genomic_start = refined_right - edge_mh_len
            right_mh_genomic_end   = refined_right - 1
            left_mh_genomic_start  = refined_left
            left_mh_genomic_end    = refined_left  + edge_mh_len - 1
        print(f"  → Edge MH detected ({edge_mh_len} bp): {mh_seq!r}")
        print(f"  → Left  MH genomic: {left_mh_genomic_start}–{left_mh_genomic_end}")
        print(f"  → Right MH genomic: {right_mh_genomic_start}–{right_mh_genomic_end}")
        return len(mh_seq), mh_seq, left_mh_genomic_start, left_mh_genomic_end, \
               right_mh_genomic_start, right_mh_genomic_end

    # ── Step 4: Variant walking to identify the MH region (or approximate region where breakpoints occur (i.e. left junction transitions to right....)) ─────
    print(f"\n--- Variant walk ({'deletion' if is_deletion else 'duplication'} mode) ---")
    print(f"  left_vars_in_block:  {left_vars_in_block}")
    print(f"  right_vars_in_block: {right_vars_in_block}")

    if is_deletion:
        first_left_var = -1
        for i, v in enumerate(left_vars_in_block):
            if v:
                first_left_var = i
                break

        last_right_var = -1
        for i, v in enumerate(right_vars_in_block):
            if v:
                last_right_var = i

        print(f"  Initial: first_left_var={first_left_var}  last_right_var={last_right_var}")

        if last_right_var >= first_left_var and first_left_var >= 0:
            last_right_var = -1
            for i in range(first_left_var - 1, -1, -1):
                if right_vars_in_block[i]:
                    last_right_var = i
                    break
            print(f"  Adjusted: last_right_var={last_right_var}")

        if first_left_var <= last_right_var and last_right_var >= 0:
            first_left_var = -1
            for i in range(last_right_var + 1, best_k):
                if left_vars_in_block[i]:
                    first_left_var = i
                    break
            print(f"  Adjusted: first_left_var={first_left_var}")

        mh_block_start = last_right_var  + 1  if last_right_var  >= 0 else 0
        mh_block_end   = first_left_var  - 1  if first_left_var  >= 0 else (best_k - 1)

        print(f"  first_left_var genomic:  "
              f"{left_win_start  + best_l_start + first_left_var  if first_left_var  >= 0 else 'none'}")
        print(f"  last_right_var genomic:  "
              f"{right_win_start + best_r_start + last_right_var  if last_right_var  >= 0 else 'none'}")

    else:
        first_right_var = -1
        for i, v in enumerate(right_vars_in_block):
            if v:
                first_right_var = i
                break

        last_left_var = -1
        for i, v in enumerate(left_vars_in_block):
            if v:
                last_left_var = i

        print(f"  Initial: last_left_var={last_left_var}  first_right_var={first_right_var}")

        if last_left_var >= first_right_var and first_right_var >= 0:
            last_left_var = -1
            for i in range(first_right_var - 1, -1, -1):
                if left_vars_in_block[i]:
                    last_left_var = i
                    break
            print(f"  Adjusted: last_left_var={last_left_var}")

        if first_right_var <= last_left_var and last_left_var >= 0:
            first_right_var = -1
            for i in range(last_left_var + 1, best_k):
                if right_vars_in_block[i]:
                    first_right_var = i
                    break
            print(f"  Adjusted: first_right_var={first_right_var}")

        mh_block_start = last_left_var   + 1
        mh_block_end   = (first_right_var - 1) if first_right_var >= 0 else (best_k - 1)

        print(f"  last_left_var genomic:   "
              f"{left_win_start  + best_l_start + last_left_var   if last_left_var   >= 0 else 'none'}")
        print(f"  first_right_var genomic: "
              f"{right_win_start + best_r_start + first_right_var if first_right_var >= 0 else 'none'}")

    mh_len = mh_block_end - mh_block_start + 1
    if mh_len < min_mh:
        print(f"  → mh_len={mh_len} below min_mh={min_mh}, no MH called.")
        return 0, "", None, None, None, None

    mh_seq = left_consensus_aln[
        best_l_start + mh_block_start :
        best_l_start + mh_block_start + mh_len
    ]

    left_mh_genomic_start  = left_win_start  + best_l_start + mh_block_start
    left_mh_genomic_end    = left_win_start  + best_l_start + mh_block_end
    right_mh_genomic_start = right_win_start + best_r_start + mh_block_start
    right_mh_genomic_end   = right_win_start + best_r_start + mh_block_end

    print(f"\n  → MH block indices: [{mh_block_start} – {mh_block_end}]  ({mh_len} bp)")
    print(f"  → Left  MH genomic: {left_mh_genomic_start} – {left_mh_genomic_end}")
    print(f"  → Right MH genomic: {right_mh_genomic_start} – {right_mh_genomic_end}")
    print(f"  → mh_seq={mh_seq!r}")

    return mh_len, mh_seq, left_mh_genomic_start, left_mh_genomic_end, \
           right_mh_genomic_start, right_mh_genomic_end




# ─────────────────────────────────────────────────────────────────────────────
# 7. Visualisation
# ─────────────────────────────────────────────────────────────────────────────
def plot_microhomology(
    mh_seq,
    left_mh_genomic_start,
    right_mh_genomic_start,
    right_mh_genomic_end,
    left_junc_left_mh,    # flanking bases LEFT  of left  MH copy
    left_junc_right_mh,   # flanking bases RIGHT of left  MH copy (into interior)
    right_junc_left_mh,   # flanking bases LEFT  of right MH copy (out of interior)
    right_junc_right_mh,  # flanking bases RIGHT of right MH copy
    sv_type,
    sample_id="",
    gene="",
    chrom="",
    save_path="microhomology_viz.png",
):
    MONO   = "DejaVu Sans Mono"
    SEQ_FS = 10.5

    mh = (mh_seq or "").lower()

    left_mh_genomic_end  = left_mh_genomic_start  + len(mh)
    # right_mh_genomic_end passed in directly

    interior_bp = (
        right_mh_genomic_start
        - left_mh_genomic_start
        - len(mh)
        - len(left_junc_right_mh)
    )
    gap_str = f" \u2013({interior_bp:,}bp)\u2013 "

    chrom_d  = chrom or "chr?"
    sv_label = "duplication" if sv_type != "deletion" else "deletion"

    dup_s = left_mh_genomic_start
    dup_e = right_mh_genomic_end
    seq_s = left_mh_genomic_start  - len(left_junc_left_mh)  - 1
    seq_e = right_mh_genomic_end   + len(right_junc_right_mh) + 1

    wt_pre = f"WT  {chrom_d}:{seq_s:,} \u2013 "
    wt_suf = f" \u2013 {chrom_d}:{seq_e:,}"

    pieces = [
        (wt_pre,                      False, "black"),
        (left_junc_left_mh.lower(),   False, "black"),
        ("|",                         False, "black"),
        (mh,                          True,  "black"),
        ("|",                         False, "black"),
        (left_junc_right_mh.lower(),  False, "black"),
        (gap_str,                     False, "black"),
        (right_junc_left_mh.lower(),  False, "black"),
        ("|",                         False, "black"),
        (mh,                          True,  "black"),
        ("|",                         False, "black"),
        (right_junc_right_mh.lower(), False, "black"),
        (wt_suf,                      False, "black"),
    ]

    total_chars = sum(len(t) for t, _, _ in pieces)
    CHAR_W_IN   = SEQ_FS * 0.6 / 72.0
    fig_w = max(14.0, total_chars * CHAR_W_IN + 1.2)
    fig_h = 3.6

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.axis("off")

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

    widths  = [char_frac(txt, bold) for txt, bold, _ in pieces]
    total_w = sum(widths)

    SEQ_Y   = 0.22
    x_start = (1.0 - total_w) / 2.0

    x = x_start
    pipe_x = []
    for i, (txt, bold, _) in enumerate(pieces):
        if txt == "|":
            pipe_x.append(x + widths[i] / 2)
        x += widths[i]

    # pipe_x = [left_outer, left_inner, right_inner, right_outer]
    brace_l   = pipe_x[0]
    brace_r   = pipe_x[3]
    brace_mid = (brace_l + brace_r) / 2

    gene_str  = f" $\\it{{{gene}}}$" if gene else ""
    title_str = (
        f"Proposed {sv_label} structure"
        + (f" for {sample_id}" if sample_id else "")
        + gene_str
        + (f" tandem {sv_label}" if sv_type != "deletion" else "")
    )
    ax.text(0.5, 0.97, title_str,
            ha="center", va="top", fontsize=13, fontweight="bold",
            fontfamily="DejaVu Sans", color="black",
            transform=ax.transAxes)

    seg_lbl = (
        f"{'Duplicated' if sv_type != 'deletion' else 'Deleted'} segment"
        f"  \u2013  {chrom_d}:{dup_s:,} \u2013{dup_e:,}"
    )
    SEG_Y = 0.72
    ax.text(brace_mid, SEG_Y, seg_lbl,
            ha="center", va="top", fontsize=10,
            fontfamily=MONO, color="black",
            transform=ax.transAxes)

    BRACE_Y = 0.52
    TICK_H  = 0.30

    ax.plot([brace_l, brace_r], [BRACE_Y, BRACE_Y],
            color="black", lw=1.5, transform=ax.transAxes, clip_on=False)
    ax.plot([brace_l, brace_l], [BRACE_Y, BRACE_Y - TICK_H],
            color="black", lw=1.5, transform=ax.transAxes, clip_on=False)
    ax.plot([brace_r, brace_r], [BRACE_Y, BRACE_Y - TICK_H],
            color="black", lw=1.5, transform=ax.transAxes, clip_on=False)
    ax.annotate("", xy=(brace_mid, BRACE_Y),
                xytext=(brace_mid, SEG_Y - 0.02),
                xycoords="axes fraction", textcoords="axes fraction",
                arrowprops=dict(arrowstyle="-", color="black", lw=0.8))

    x = x_start
    for (txt, bold, color), w in zip(pieces, widths):
        if txt:
            ax.text(x, SEQ_Y, txt,
                    fontsize=SEQ_FS, fontfamily=MONO,
                    fontweight="bold" if bold else "normal",
                    color=color, ha="left", va="center",
                    transform=ax.transAxes)
        x += w

    plt.tight_layout(pad=0.4)
    plt.savefig(save_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved: {save_path}")
    display(Image(save_path))



# ─────────────────────────────────────────────────────────────────────────────
# 8. Main entry point
# ─────────────────────────────────────────────────────────────────────────────
def analyze_sv_pair(
    bam_path, ref_path,
    chrom, refined_left, refined_right,
    flank=600, min_mh=2, show_pileup=False,
    min_mapq=60, pair_label="", save_dir=".",
    show_homozygous_variants=False
):
    """Run the full MH pipeline for a single breakpoint pair. Returns a result dict."""
    padding = max(500, (refined_right - refined_left) // 2)
    reads = extract_split_reads(bam_path, chrom,
                            start_search=refined_left,
                            end_search=refined_right,
                            padding=min(padding, refined_left))
    if not reads:
        print(f"[{pair_label}] No split reads found — skipping.")
        return None

    # ── SV type detection (same logic as current auto function) ──────────────
    tol = max(50, padding // 10)
    left_clips_at_left_bp, right_clips_at_left_bp = 0, 0
    REF_OPS_SA = {"M", "D", "N", "=", "X"}
    for read in reads:
        if read.is_unmapped or not read.cigartuples:
            continue
        cigars = read.cigartuples
        if cigars[0][0] in CLIP_OPS and abs(read.reference_start - refined_left) < tol:
            left_clips_at_left_bp += 1
        if cigars[-1][0] in CLIP_OPS and abs(read.reference_end - refined_left) < tol:
            right_clips_at_left_bp += 1
        if read.has_tag("SA"):
            for sa_entry in read.get_tag("SA").rstrip(";").split(";"):
                fields = sa_entry.split(",")
                if len(fields) < 4 or fields[0] != chrom:
                    continue
                sa_pos   = int(fields[1])
                sa_cigar = fields[3]
                ops = re.findall(r'(\d+)([MIDNSHP=X])', sa_cigar)
                if not ops:
                    continue
                sa_ref_len = sum(int(l) for l, op in ops if op in REF_OPS_SA)
                sa_end = sa_pos + sa_ref_len
                if ops[0][1] in {"S", "H"} and abs(sa_pos - refined_left) < tol:
                    left_clips_at_left_bp += 1
                if ops[-1][1] in {"S", "H"} and abs(sa_end - refined_left) < tol:
                    right_clips_at_left_bp += 1

    sv_type = "deletion" if right_clips_at_left_bp >= left_clips_at_left_bp else "duplication"
    print(f"[{pair_label}] sv_type={sv_type}  {chrom}:{refined_left}–{refined_right}")

    fasta = pysam.FastaFile(ref_path)
    consensus_result = build_split_consensus(
        reads, refined_left, refined_right, flank, fasta, chrom, sv_type=sv_type,
        bam_path=bam_path, show_pileup=show_pileup,
        show_homozygous_variants=show_homozygous_variants,
        pileup_save_dir=save_dir, pair_label=pair_label,
    )
    fasta.close()

    if consensus_result[0] is None:   # chrom not in reference — build_split_consensus bailed out early
        print(f"{pair_label}: skipping — {chrom} not found in reference.")
        return None

    left_consensus, right_consensus, left_var_mask, right_var_mask, right_bp, ref_left, ref_right = consensus_result

    (mh_len, mh_seq, left_mh_start, left_mh_end,
     right_mh_start, right_mh_end) = find_mmej_microhomology(
        left_consensus,  left_var_mask,
        right_consensus, right_var_mask,
        ref_left, ref_right,
        refined_left, right_bp,
        flank=flank, sv_type=sv_type, min_mh=min_mh,
        pair_label=pair_label,
        save_dir=save_dir,
    )

    if mh_len and mh_seq:
        FLANK = 10
        fasta = pysam.FastaFile(ref_path)
        left_junc_left_mh   = fasta.fetch(chrom, left_mh_start  - FLANK, left_mh_start).upper()
        left_junc_right_mh  = fasta.fetch(chrom, left_mh_end    + 1,     left_mh_end   + 1 + FLANK).upper()
        right_junc_left_mh  = fasta.fetch(chrom, right_mh_start - FLANK, right_mh_start).upper()
        right_junc_right_mh = fasta.fetch(chrom, right_mh_end   + 1,     right_mh_end  + 1 + FLANK).upper()
        fasta.close()
        plot_microhomology(
            mh_seq=mh_seq,
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
        "pair_label":              pair_label,
        "chrom":                   chrom,
        "sv_type":                 sv_type,
        "left_clip_pos":           refined_left,
        "right_clip_pos":          right_bp,
        "left_mh_genomic_start":   left_mh_start,
        "left_mh_genomic_end":     left_mh_end,
        "right_mh_genomic_start":  right_mh_start,
        "right_mh_genomic_end":    right_mh_end,
        "microhomology_length":    mh_len,
        "microhomology_sequence":  mh_seq,
        "supporting_reads":        len(reads),
    }


def count_bridging_reads(bam_path, chrom, pos_a, pos_b, tol=200, min_mapq=0):
    """
    Count reads whose primary alignment clips near pos_a (or pos_b) AND whose
    SA tag places a supplementary alignment near the other position.
    This is the direct evidence that two clip sites belong to the same SV.
    """
    REF_OPS = {"M", "D", "N", "=", "X"}
    bridging = set()

    def _clips_near(cigar_ops, ref_start, target):
        if not cigar_ops:
            return False
        ref_end = ref_start + sum(l for op, l in cigar_ops if op in (0,2,3,7,8))
        right_clip = cigar_ops[-1][0] in (4, 5) and abs(ref_end   - target) <= tol
        left_clip  = cigar_ops[0][0]  in (4, 5) and abs(ref_start - target) <= tol
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
            ops = re.findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
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
        fetch_end   =        max(pos_a, pos_b) + tol
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

def pair_clip_positions(clip_df, bam_path=None, max_dist=None,
                        strategy="split_read_bridge", bridge_tol=200):
    """
    Build (chrom, left_pos, right_pos) pairs from find_major_clip_positions() output.

    strategy="split_read_bridge" : (default) score every candidate pair by the number
                             of reads whose SA tag bridges BOTH positions. Pairs are
                             ranked by bridge_count descending and committed greedily.
                             Requires bam_path. Falls back to greedy_rank if not given.
    strategy="greedy_rank" : iterate clips in support-rank order; pair each unmatched
                             clip with the highest-ranked partner on the same chrom
                             within max_dist. No BAM access needed.
    strategy="proximity"   : match each right clip to the nearest downstream left clip
                             on the same chrom (requires correctly labelled sides).
    strategy="consecutive" : pair rank-1+2, rank-3+4, etc.
    """
    pairs = []

    if strategy == "split_read_bridge":
        if bam_path is None:
            print("  [pair_clip_positions] bam_path not supplied; falling back to greedy_rank.")
            return pair_clip_positions(clip_df, bam_path=None, max_dist=max_dist,
                                       strategy="greedy_rank")

        idxs = list(clip_df.index)
        n    = len(idxs)
        candidates = []

        print(f"  Scoring {n*(n-1)//2} candidate pairs by bridging read count ...")
        for ii in range(n):
            for jj in range(ii + 1, n):
                i, j  = idxs[ii], idxs[jj]
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
                    bam_path,
                    row_a["chrom"],
                    int(row_a["position"]),
                    int(row_b["position"]),
                    tol=bridge_tol,
                )
                candidates.append((bc, i, j))
                print(f"    {row_a['chrom']}:{int(row_a['position']):,} ↔ "
                      f"{int(row_b['position']):,}  bridge_reads={bc}")

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
            left  = int(min(row_a["position"], row_b["position"]))
            right = int(max(row_a["position"], row_b["position"]))
            pairs.append((row_a["chrom"], left, right))

        if not pairs:
            print("  [split_read_bridge] No pairs with bridging reads found; "
                  "falling back to greedy_rank.")
            return pair_clip_positions(clip_df, bam_path=None, max_dist=max_dist,
                                       strategy="greedy_rank")

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
            left  = int(min(row_a["position"], row_b["position"]))
            right = int(max(row_a["position"], row_b["position"]))
            pairs.append((row_a["chrom"], left, right))

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
        used        = set()
        for _, rc in right_clips.iterrows():
            cands = left_clips[
                (left_clips["chrom"] == rc["chrom"]) &
                (left_clips["position"] > rc["position"]) &
                (~left_clips.index.isin(used))
            ]
            if max_dist is not None:
                cands = cands[cands["position"] - rc["position"] <= max_dist]
            if cands.empty:
                continue
            best = cands.iloc[0]
            used.add(best.name)
            pairs.append((rc["chrom"], int(rc["position"]), int(best["position"])))

    return pairs


def batch_microhomology_search(
    bam_path, ref_path,
    flank=600, min_mh=2, show_pileup=False,
    min_clip_len=50, min_support=5, min_support_frac=0.05,  # ← new
    cluster_dist=10, min_mapq=60,
    pair_strategy="split_read_bridge", max_dist=None, top_n=50,
    save_dir=None, show_homozygous_variants=False
):
    import os
    bam_stem = os.path.splitext(os.path.basename(bam_path))[0]
    if save_dir is None:
        save_dir = os.path.join("results", bam_stem)
    os.makedirs(save_dir, exist_ok=True)
    print(f"Output directory: {save_dir}")

    # ── 1. Discover all major clip positions ──────────────────────────────────
    clip_df = find_major_clip_positions(
        bam_path, chrom=None,
        min_clip_len=min_clip_len, min_support=min_support,
        min_support_frac=min_support_frac,
        cluster_dist=cluster_dist, top_n=top_n,
        require_sa=False, min_mapq=min_mapq,
        save_dir=save_dir,
    )
    if clip_df.empty:
        print("No clip positions found.")
        return pd.DataFrame()
    
    # ── Drop clip positions whose chrom is absent from the reference ──────────
    ref_chroms = get_ref_chroms(ref_path)
    missing = set(clip_df["chrom"].unique()) - ref_chroms
    if missing:
        print(f"WARNING: {len(missing)} contig(s) in BAM not found in reference — skipping: {missing}")
        clip_df = clip_df[clip_df["chrom"].isin(ref_chroms)].reset_index(drop=True)
    if clip_df.empty:
        print("No clip positions remain after filtering to reference contigs.")
        return pd.DataFrame()

    # ── 2. Drop N-adjacent positions ──────────────────────────────────────────
    fasta = pysam.FastaFile(ref_path)
    n_mask = clip_df.apply(
        lambda r: _has_adjacent_n(fasta, r["chrom"], r["position"]), axis=1
    )
    fasta.close()
    clip_df = clip_df[~n_mask].reset_index(drop=True)

    # ── 3. Build pairs ────────────────────────────────────────────────────────
    pairs = pair_clip_positions(clip_df, bam_path=bam_path, max_dist=max_dist,
                              strategy=pair_strategy)
    print(f"\n{len(pairs)} pairs to analyse (strategy='{pair_strategy}'):")
    for i, (chrom, lp, rp) in enumerate(pairs, 1):
        print(f"  {i:>3}.  {chrom}  {lp:>12,} – {rp:>12,}  (span {rp-lp:,} bp)")

    # ── 4. Run per-pair MH search ─────────────────────────────────────────────
    results = []
    for i, (chrom, left_bp, right_bp) in enumerate(pairs, 1):
        label = f"p{i:02d}_{chrom}_{left_bp}"
        print(f"\n{'='*64}")
        print(f"  Pair {i}/{len(pairs)}:  {chrom}  {left_bp:,} – {right_bp:,}")
        print(f"{'='*64}")
        row = analyze_sv_pair(
            bam_path, ref_path,
            chrom=chrom, refined_left=left_bp, refined_right=right_bp,
            flank=flank, min_mh=min_mh, show_pileup=show_pileup,
            min_mapq=min_mapq, pair_label=label, save_dir=save_dir,
            show_homozygous_variants=show_homozygous_variants
        )
        if row:
            results.append(row)

    # ── 5. Save summary ───────────────────────────────────────────────────────
    summary_df = pd.DataFrame(results)
    if not summary_df.empty:
        out = os.path.join(save_dir, "batch_mh_results.csv")
        summary_df.to_csv(out, index=False)
        print(f"\nBatch complete. Results saved → {out}")
        display(summary_df)
    return summary_df


