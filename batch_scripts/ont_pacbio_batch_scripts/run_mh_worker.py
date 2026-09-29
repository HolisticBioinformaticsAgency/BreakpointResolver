#!/usr/bin/env python3
"""
run_mh_worker.py
Runs the configured pipeline mode for ONE barcode + platform.
Called by the Slurm array job with:
    python3 run_mh_worker.py <barcode> <ONT|PacBio> [params_batch.yaml]


All tunable parameters are read from params_batch.yaml (default: same
directory as this script).  The 'mode' key selects which analysis module
and entry-point function to use:


    sv_analysis                    → sv_analysis.py
    sv_analysis_with_sc_adjustment → sv_analysis_with_sc_adjustment.py
    sv_rearrangement_path          → sv_rearrangement_path.py
    sv_analysis_junction_rearr     → sv_analysis_junction_rearr.py
    sv_analysis_sc_adjustment_arc  → sv_analysis_sc_adjustment_arc.py
                                     (INTEGRATED SC-adjustment MH + arc output
                                     for every exact SC-adjustment cluster)
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
mode     = P.get("mode", "sv_analysis_with_sc_adjustment")
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
    "sv_analysis":                    "sv_analysis.py",
    "sv_analysis_with_sc_adjustment": "sv_analysis_with_sc_adjustment.py",
    "sv_rearrangement_path":          "sv_rearrangement_path.py",
    "sv_analysis_junction_rearr":     "sv_analysis_junction_rearr.py",
    "sv_analysis_sc_adjustment_arc":  "sv_analysis_sc_adjustment_arc.py",
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
try:
    if mode in ("sv_analysis", "sv_analysis_with_sc_adjustment"):
        kwargs = dict(
            # core
            bam_path             = bam_path,
            ref_path             = ref_path,
            save_dir             = save_dir,
            # Stage 1 — clip discovery
            top_n                = p("top_n",                20),
            min_clip_len         = p("min_clip_len",         50),
            min_support          = p("min_support",          5),
            min_support_frac     = p("min_support_frac",     0.05),
            cluster_dist         = p("cluster_dist",         10),
            min_mapq             = p("min_mapq",             60),
            # Stage 3 — pairing
            max_dist             = p("max_dist",             500_000),
            # Stage 5 — analysis
            flank                = p("flank",                None),
            min_mh               = p("min_mh",               2),
            # display / output
            show_pileup              = p("show_pileup",              False),
            show_homozygous_variants = p("show_homozygous_variants", False),
        )
        # sc_adjustment-only params
        if mode == "sv_analysis":
            kwargs["pair_strategy"] = p("pair_strategy",        "split_read_bridge"),
        
        elif mode == "sv_analysis_with_sc_adjustment":
            kwargs["sa_share_threshold"] = p("sa_share_threshold", 0.50)
            kwargs["bridge_tol"]         = p("bridge_tol",         200)


        results = sv.batch_microhomology_search(**kwargs)


        if results is not None and not results.empty:
            out_csv = os.path.join(save_dir, "batch_mh_results.csv")
            results.to_csv(out_csv, index=False)
            print(f"Done — {len(results)} pair(s) found. Saved to {out_csv}")
        else:
            print("Done — no MH pairs detected.")


    elif mode == "sv_rearrangement_path":
        kwargs = dict(
            # core
            bam_path                  = bam_path,
            ref_path                  = ref_path,
            save_dir                  = save_dir,
            # Stage 1 — clip discovery (shared)
            top_n                     = p("top_n",                    20),
            min_clip_len              = p("min_clip_len",              50),
            min_support               = p("min_support",               5),
            min_support_frac          = p("min_support_frac",          0.05),
            cluster_dist              = p("cluster_dist",              10),
            min_mapq                  = p("min_mapq",                  60),
            # Stage 3 — pairing (shared)
            pair_strategy             = p("pair_strategy",             "split_read_bridge"),
            max_dist                  = p("max_dist",                  500_000),
            # rearrangement-path-specific
            min_segment_len           = p("min_segment_len",           10),
            junction_tol              = p("junction_tol",              50),
            ref_bin_tol               = p("ref_bin_tol",               100),
            min_block_support         = p("min_block_support",         2),
            fetch_padding             = p("fetch_padding",             500_000),
            landing_tol               = p("landing_tol",               200),
            max_aln_span              = p("max_aln_span",              1_000_000),
            minimap2_preset           = p("minimap2_preset",           "map-hifi"),
            min_identity              = p("min_identity",              0.90),
            match_score               = p("match_score",               2),
            mismatch_penalty          = p("mismatch_penalty",          8),
            gap_open                  = p("gap_open",                  4),
            gap_extend                = p("gap_extend",                2),
            gap_open2                 = p("gap_open2",                 24),
            gap_extend2               = p("gap_extend2",               1),
            pdv_min_split_vaf         = p("pdv_min_split_vaf",         0.20),
            pdv_min_enrichment        = p("pdv_min_enrichment",        3.0),
            background_hom_vaf        = p("background_hom_vaf",        0.95),
            min_graph_segment_support = p("min_graph_segment_support", 10),
        )
        results = sv.batch_rearrangement_path_search(**kwargs)


        if results is not None and not results.empty:
            out_csv = os.path.join(save_dir, "batch_rearrangement_results.csv")
            results.to_csv(out_csv, index=False)
            print(f"Done — {len(results)} pair(s) found. Saved to {out_csv}")
        else:
            print("Done — no rearrangement paths detected.")


    elif mode == "sv_analysis_junction_rearr":
        # Standalone SNP-anchored rearrangement arc diagram mode. It runs its
        # OWN independent clip discovery / pairing / clustering pipeline.
        # Use sv_analysis_sc_adjustment_arc below when MH and arcs must be
        # generated for the same SC-adjustment cluster set.
        kwargs = dict(
            bam_path               = bam_path,
            ref_path               = ref_path,
            sample_name            = f"{BARCODE}_{PLATFORM}",
            out_dir                = save_dir,
            min_mapq               = p("min_mapq",                     60),
            min_support            = p("min_support",                  5),
            bridge_tol             = p("bridge_tol",                   200),
            share_threshold        = p("sa_share_threshold",           0.50),
            padding                = p("arc_padding",                  500),
            min_perfect_run        = p("arc_min_perfect_run",          10),
            min_non_split_support  = p("arc_min_non_split_support",    3),
            min_non_split_frac     = p("arc_min_non_split_frac",       0.1),
            supp_mapq_min          = p("arc_supp_mapq_min",            0),
            min_node_support       = p("arc_min_node_support",         2),
            context_flank          = p("arc_context_flank",            50),
            context_min_identity   = p("arc_context_min_identity",     1.0),
            context_search_revcomp = p("arc_context_search_revcomp",   True),
            max_clusters           = p("arc_max_clusters",             None),
            show_df                = False,
        )
        result = sv.run_pipeline(**kwargs)


        graphs = (result or {}).get("graphs", [])
        if graphs:
            out_csv = os.path.join(save_dir, "arc_clusters_summary.csv")
            import pandas as pd
            pd.DataFrame([
                {
                    "chrom":         g["cluster"]["chrom"],
                    "refined_left":  g["cluster"]["refined_left"],
                    "refined_right": g["cluster"]["refined_right"],
                    "n_nodes":       g["graph"].number_of_nodes(),
                    "n_edges":       g["graph"].number_of_edges(),
                    "html":          g["html"],
                }
                for g in graphs
            ]).to_csv(out_csv, index=False)
            print(f"Done — {len(graphs)} cluster arc diagram(s) generated. "
                  f"Summary saved to {out_csv}")
        else:
            print("Done — no clusters produced a non-empty graph.")


    elif mode == "sv_analysis_sc_adjustment_arc":
        # INTEGRATED mode: SC-adjustment's find_major_clip_pairs() and
        # cluster_clip_pairs_by_reads() define the one shared cluster set.
        # For EVERY retained cluster, the integration module runs both
        # analyze_sv_cluster() (MH/consensus output) and the arc graph/HTML.
        # This is intentionally different from sv_analysis_junction_rearr:
        # there is no second, conflicting rediscovery/reclustering step.
        kwargs = dict(
            # Shared SC-adjustment discovery + cluster parameters
            bam_path                     = bam_path,
            ref_path                     = ref_path,
            save_dir                     = save_dir,
            top_n                        = p("top_n",                20),
            min_clip_len                 = p("min_clip_len",         50),
            min_support                  = p("min_support",          5),
            min_support_frac             = p("min_support_frac",     0.05),
            cluster_dist                 = p("cluster_dist",         10),
            min_mapq                     = p("min_mapq",             60),
            max_dist                     = p("max_dist",             500_000),
            bridge_tol                   = p("bridge_tol",           200),
            # SC-adjustment MH parameters
            flank                        = p("flank",                None),
            min_mh                       = p("min_mh",               2),
            show_pileup                  = p("show_pileup",          False),
            show_homozygous_variants     = p("show_homozygous_variants", False),
            # Arc-only graph/context parameters
            arc_padding                  = p("arc_padding",                  500),
            arc_min_perfect_run          = p("arc_min_perfect_run",          10),
            arc_min_non_split_support    = p("arc_min_non_split_support",    3),
            arc_min_non_split_frac       = p("arc_min_non_split_frac",       0.1),
            arc_supp_mapq_min            = p("arc_supp_mapq_min",            0),
            arc_min_node_support         = p("arc_min_node_support",         2),
            arc_context_flank            = p("arc_context_flank",            50),
            arc_context_min_identity     = p("arc_context_min_identity",     1.0),
            arc_context_search_revcomp   = p("arc_context_search_revcomp",   True),
            arc_max_clusters             = p("arc_max_clusters",             None),
        )
        results = sv.batch_microhomology_arc_search(**kwargs)


        if results is not None and not results.empty:
            out_csv = os.path.join(save_dir, "batch_mh_arc_results.csv")
            results.to_csv(out_csv, index=False)
            print(f"Done — integrated MH + arc output for {len(results)} cluster(s). "
                  f"Saved to {out_csv}")
        else:
            print("Done — no integrated MH result rows detected. Check "
                  "arc_clusters_summary.csv for arc-only status/errors.")


except Exception:
    traceback.print_exc()
    sys.exit(1)
