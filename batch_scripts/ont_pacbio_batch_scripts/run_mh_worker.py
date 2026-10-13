#!/usr/bin/env python3
"""
run_mh_worker.py
Runs the configured pipeline mode for ONE barcode + platform.
Called by the Slurm array job with:
    python3 run_mh_worker.py <barcode> <ONT|PacBio> [params_batch.yaml]


All tunable parameters are read from params_batch.yaml (default: same
directory as this script).  The 'mode' key selects which analysis module
and entry-point function to use:


    sv_analysis_sc_adjustment_arc  → sv_analysis_sc_adjustment_arc.py  (DEFAULT)
                                     SC-adjustment MH + interactive arc diagram
                                     for every SC-adjustment cluster
    sv_analysis_with_sc_adjustment → sv_analysis_with_sc_adjustment.py
                                     SC-adjustment MH only
"""
import os, sys, importlib.util, traceback


# ── Argument parsing ──────────────────────────────────────────────────────────
if len(sys.argv) < 3:
    print("Usage: run_mh_worker.py <barcode> <ONT|PacBio> [params_batch.yaml]")
    sys.exit(1)


BARCODE  = sys.argv[1]
PLATFORM = sys.argv[2]
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


PARAMS_FILE = sys.argv[3] if len(sys.argv) >= 4 else os.path.join(SCRIPT_DIR, "params_batch.yaml")
if not os.path.isfile(PARAMS_FILE):
    print(f"ERROR: params file not found — {PARAMS_FILE}")
    sys.exit(1)


# ── Load YAML ─────────────────────────────────────────────────────────────────
try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not installed. Run: pip install pyyaml")
    sys.exit(1)


with open(PARAMS_FILE) as fh:
    P = yaml.safe_load(fh)


# ── BAM path templates ────────────────────────────────────────────────────────
# ONT barcodes 1113 and 1114 live in a different demultiplex batch
# (dorado_demultiplex_round2_full_ont) than the rest of the ONT barcodes.
ONT_DEMUX_FOLDERS = {
    "default": "dorado_demultiplex_142pods_barcode_both_ends",
}
ONT_ROUND2_BARCODES = {
    "1113", "1114", "1115", "1116", "1117", "1118", "1119", "1120",
    "1121", "1122", "1123", "1124", "1125", "1126", "1127", "1128"
}

ONT_BAM_TMPL_DEFAULT = (
    "/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/"
    "dorado_demultiplex_142pods_barcode_both_ends/bams_dedup_cleaned/"
    "barcode{bc}_q10_mq20_dedup.bam"
)
ONT_BAM_TMPL_ROUND2 = (
    "/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/"
    "dorado_demultiplex_round2_full_ont/bams_dedup_cleaned/"
    "barcode{bc}_q10_mq20_dedup.bam"
)
PACBIO_BAM_TMPL = (
    "/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/"
    "demultiplex_pb/workflow_HiFiTargetEnrichment/batches/target_hifit/"
    "bc{bc}/realigned_minimap2/bc{bc}.minimap2.GRCh38.bam"
)


def resolve_ont_bam_tmpl(barcode: str) -> str:
    """Return the correct ONT BAM template for the given barcode."""
    if barcode in ONT_ROUND2_BARCODES:
        return ONT_BAM_TMPL_ROUND2
    return ONT_BAM_TMPL_DEFAULT


if PLATFORM not in ("ONT", "PacBio"):
    print(f"ERROR: platform must be ONT or PacBio, got: {PLATFORM!r}")
    sys.exit(1)


if PLATFORM == "ONT":
    bam_path = resolve_ont_bam_tmpl(BARCODE).format(bc=BARCODE)
else:
    bam_path = PACBIO_BAM_TMPL.format(bc=BARCODE)


if not os.path.isfile(bam_path):
    print(f"WARNING: BAM not found — {bam_path}")
    print("Nothing to do. Exiting cleanly.")
    sys.exit(0)


# ── Resolve output directory ───────────────────────────────────────────────────
mode     = P.get("mode", "sv_analysis_sc_adjustment_arc")
ref_path = P["ref_path"]


raw_save = P.get("save_dir", "")
if raw_save:
    save_dir = raw_save
else:
    bam_stem = os.path.splitext(os.path.basename(bam_path))[0]
    save_dir = os.path.join(SCRIPT_DIR, "results", f"{BARCODE}_{PLATFORM}", bam_stem, mode)


os.makedirs(save_dir, exist_ok=True)


print(f"Barcode : {BARCODE}")
print(f"Platform: {PLATFORM}")
print(f"Mode    : {mode}")
print(f"BAM     : {bam_path}")
print(f"Ref     : {ref_path}")
print(f"Out dir : {save_dir}")
print(f"Params  : {PARAMS_FILE}")


# ── Module map ────────────────────────────────────────────────────────────────
MODULE_FILES = {
    "sv_analysis_sc_adjustment_arc":  "sv_analysis_sc_adjustment_arc.py",
    "sv_analysis_with_sc_adjustment": "sv_analysis_with_sc_adjustment.py",
}
if mode not in MODULE_FILES:
    print(f"ERROR: unknown mode {mode!r}. "
          f"Choose one of: {list(MODULE_FILES)}")
    sys.exit(1)


module_path = os.path.join(SCRIPT_DIR, MODULE_FILES[mode])
if not os.path.isfile(module_path):
    print(f"ERROR: module file not found — {module_path}")
    sys.exit(1)


spec = importlib.util.spec_from_file_location(mode, module_path)
sv   = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sv)


# ── Helper: pull a param, with a fallback default ─────────────────────────────
def p(key, default=None):
    val = P.get(key, default)
    return None if val in ("", "null", "None", None) else val


# ── Dispatch by mode ──────────────────────────────────────────────────────────
# Parameters shared by both modes (SC-adjustment discovery, clustering, MH).
kwargs = dict(
    # core
    bam_path                 = bam_path,
    ref_path                 = ref_path,
    save_dir                 = save_dir,
    # Stages A-D — clip discovery / pairing
    top_n                    = p("top_n",                20),
    min_clip_len             = p("min_clip_len",         50),
    min_support              = p("min_support",          5),
    min_support_frac         = p("min_support_frac",     0.05),
    cluster_dist             = p("cluster_dist",         10),
    min_mapq                 = p("min_mapq",             60),
    max_dist                 = p("max_dist",             500_000),
    # in-CIGAR deletion reads
    min_del_len              = p("min_del_len",          50),
    del_match_tol            = p("del_match_tol",        20),
    # contig filtering
    edge_buffer              = p("edge_buffer",          None),
    allowed_chroms           = p("allowed_chroms",       None),
    # Stage E — clustering
    bridge_tol               = p("bridge_tol",           200),
    # Stage F — MH analysis
    flank                    = p("flank",                None),
    min_mh                   = p("min_mh",               2),
    # display / output
    show_pileup              = p("show_pileup",              False),
    show_homozygous_variants = p("show_homozygous_variants", False),
)

try:
    if mode == "sv_analysis_sc_adjustment_arc":
        # MH for every SC-adjustment cluster, plus an arc diagram per cluster
        # via batch_microhomology_search(per_cluster_hook=...).
        kwargs.update(
            arc_padding                = p("arc_padding",                  500),
            arc_min_perfect_run        = p("arc_min_perfect_run",          10),
            arc_min_non_split_support  = p("arc_min_non_split_support",    3),
            arc_min_non_split_frac     = p("arc_min_non_split_frac",       0.1),
            arc_supp_mapq_min          = p("arc_supp_mapq_min",            0),
            arc_min_node_support       = p("arc_min_node_support",         2),
            arc_context_flank          = p("arc_context_flank",            50),
            arc_context_min_identity   = p("arc_context_min_identity",     1.0),
            arc_context_search_revcomp = p("arc_context_search_revcomp",   True),
            arc_max_clusters           = p("arc_max_clusters",             None),
        )
        results = sv.batch_microhomology_arc_search(**kwargs)

    elif mode == "sv_analysis_with_sc_adjustment":
        results = sv.batch_microhomology_search(**kwargs)

    # batch_microhomology_search writes batch_mh_results_<bam>.csv itself.
    if results is not None and not results.empty:
        print(f"Done — {len(results)} MH result row(s) in {save_dir}")
    else:
        print("Done — no MH pairs detected.")


except Exception:
    traceback.print_exc()
    sys.exit(1)
