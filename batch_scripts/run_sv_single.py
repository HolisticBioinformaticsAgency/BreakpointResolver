#!/usr/bin/env python3
"""
run_sv_single.py
────────────────
Unified single-BAM runner for all three SV analysis pipelines.

Usage
-----
python3 run_sv_single.py                       # uses params_single.yaml in same directory
python3 run_sv_single.py /path/to/my_params_single.yaml

Modes (set mode: in params_single.yaml)
-----------------------------------
sv_analysis
    → batch_microhomology_search() in sv_analysis.py
    Pileup-based left+right consensus, MMEJ microhomology search, pileup PNGs.

sv_analysis_with_sc_adjustment
    → batch_microhomology_search() in sv_analysis_with_sc_adjustment.py
    Same as sv_analysis but breakpoint pairs that share >sa_share_threshold of
    their SA-supporting reads are merged (Union-Find).  The widest span in the
    cluster sets refined_left / refined_right; narrower pairs are re-anchored
    to that wider span before consensus / MH search.

sv_rearrangement_path
    → batch_rearrangement_path_search() in sv_rearrangement_path.py
    Clip-first, PDV-aware, minimap2-based k-mer path decomposition.
"""
import os
import sys
import importlib.util
import traceback

try:
    import yaml
except ImportError:
    sys.exit(
        "ERROR: PyYAML not installed.\n"
        "  Fix: pip install pyyaml (or: conda install pyyaml)"
    )

# ─────────────────────────────────────────────────────────────────────────────
# Locate and load params
# ─────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
params_file = (
    sys.argv[1]
    if len(sys.argv) > 1
    else os.path.join(SCRIPT_DIR, "params_single.yaml")
)

if not os.path.isfile(params_file):
    sys.exit(f"ERROR: params file not found: {params_file}")

with open(params_file) as fh:
    P = yaml.safe_load(fh)

# ─────────────────────────────────────────────────────────────────────────────
# Required fields
# ─────────────────────────────────────────────────────────────────────────────
bam_path = P.get("bam_path", "").strip()
ref_path = P.get("ref_path", "").strip()

if not bam_path:
    sys.exit("ERROR: bam_path is not set in params_single.yaml")
if not ref_path:
    sys.exit("ERROR: ref_path is not set in params_single.yaml")
if not os.path.isfile(bam_path):
    sys.exit(f"ERROR: BAM not found: {bam_path}")
if not os.path.isfile(ref_path):
    sys.exit(f"ERROR: Reference not found: {ref_path}")

# ─────────────────────────────────────────────────────────────────────────────
# Mode
# ─────────────────────────────────────────────────────────────────────────────
VALID_MODES = (
    "sv_analysis",
    "sv_analysis_with_sc_adjustment",
    "sv_rearrangement_path",
)
mode = P.get("mode", "sv_analysis").strip().lower()

if mode not in VALID_MODES:
    sys.exit(
        f"ERROR: Unknown mode '{mode}'\n"
        f"  Valid choices: {VALID_MODES}"
    )

# ─────────────────────────────────────────────────────────────────────────────
# Output directory
# ─────────────────────────────────────────────────────────────────────────────
save_dir = P.get("save_dir", "").strip() or None
if save_dir is None:
    bam_stem = os.path.splitext(os.path.basename(bam_path))[0]
    save_dir = os.path.join(SCRIPT_DIR, "results", bam_stem, mode)
os.makedirs(save_dir, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# Shared parameter helpers
# ─────────────────────────────────────────────────────────────────────────────
def _int(key, default):
    v = P.get(key, default)
    return None if v in (None, "", "null", "None") else int(v)

def _float(key, default):
    v = P.get(key, default)
    return None if v in (None, "", "null", "None") else float(v)

def _bool(key, default):
    v = P.get(key, default)
    if isinstance(v, bool):
        return v
    return str(v).lower() in ("true", "1", "yes")

def _str(key, default):
    return str(P.get(key, default)).strip()

# max_dist: null / "" / "None"  →  Python None (unlimited)
max_dist = _int("max_dist", 500_000)

# ─────────────────────────────────────────────────────────────────────────────
# Dynamically load the module that matches the chosen mode.
# sv_analysis and sv_analysis_with_sc_adjustment both expose
# batch_microhomology_search(); sv_rearrangement_path exposes
# batch_rearrangement_path_search().
# ─────────────────────────────────────────────────────────────────────────────
module_file = os.path.join(SCRIPT_DIR, f"{mode}.py")
if not os.path.isfile(module_file):
    sys.exit(
        f"ERROR: {mode}.py not found in {SCRIPT_DIR}\n"
        f"  Expected path: {module_file}"
    )

spec = importlib.util.spec_from_file_location(mode, module_file)
sv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sv)

# Disk cache for BAM scans (only sv_analysis_with_sc_adjustment has one)
if hasattr(sv, "USE_SCAN_CACHE"):
    sv.USE_SCAN_CACHE = _bool("use_scan_cache", True)

# ─────────────────────────────────────────────────────────────────────────────
# Print run banner
# ─────────────────────────────────────────────────────────────────────────────
W = 70
print("=" * W)
print(f" SV analysis runner — mode: {mode}")
print("=" * W)
print(f"  BAM         : {bam_path}")
print(f"  Ref         : {ref_path}")
print(f"  Save dir    : {save_dir}")
print(f"  Params file : {params_file}")
if hasattr(sv, "USE_SCAN_CACHE"):
    print(f"  Scan cache  : {sv.SCAN_CACHE_DIR if sv.USE_SCAN_CACHE else 'off'}")
print("-" * W)

# Shared
print(f"  top_n              : {_int('top_n', 20)}")
print(f"  max_dist           : {max_dist if max_dist is not None else 'unlimited'}")
print(f"  min_clip_len       : {_int('min_clip_len', 50)}")
print(f"  min_support        : {_int('min_support', 5)}")
print(f"  min_support_frac   : {_float('min_support_frac', 0.05)}")
print(f"  cluster_dist       : {_int('cluster_dist', 10)}")
print(f"  min_mapq           : {_int('min_mapq', 60)}")

if mode == "sv_analysis":
    print(f"  flank              : {_int('flank', 600)}")
    print(f"  min_mh             : {_int('min_mh', 2)}")
    print(f"  show_pileup        : {_bool('show_pileup', False)}")
    print(f"  show_hom_vars      : {_bool('show_homozygous_variants', False)}")

elif mode == "sv_analysis_with_sc_adjustment":
    print(f"  flank         : {_int('flank', 600)}")
    print(f"  min_mh             : {_int('min_mh', 2)}")
    print(f"  show_pileup        : {_bool('show_pileup', False)}")
    print(f"  show_hom_vars      : {_bool('show_homozygous_variants', False)}")
    print(f"  sa_share_threshold : {_float('sa_share_threshold', 0.50)}")
    print(f"  bridge_tol         : {_int('bridge_tol', 200)}")

elif mode == "sv_rearrangement_path":
    print(f"  sa_share_threshold : {_float('sa_share_threshold', 0.50)}")
    print(f"  landing_tol        : {_int('landing_tol', 200)}")
    print(f"  fetch_padding      : {_int('fetch_padding', 500_000):,}")
    print(f"  max_aln_span       : {_int('max_aln_span', 1_000_000):,}")
    print(f"  minimap2_preset    : {_str('minimap2_preset', 'map-hifi')}")
    print(f"  min_identity       : {_float('min_identity', 0.90)}")
    print(f"  mismatch_penalty   : {_int('mismatch_penalty', 8)}")
    print(f"  gap_open           : {_int('gap_open', 4)},{_int('gap_open2', 24)}")
    print(f"  gap_extend         : {_int('gap_extend', 2)},{_int('gap_extend2', 1)}")
    print(f"  pdv_min_split_vaf  : {_float('pdv_min_split_vaf', 0.20)}")
    print(f"  pdv_min_enrich     : {_float('pdv_min_enrichment', 3.0)}")
    print(f"  background_hom_vaf : {_float('background_hom_vaf', 0.95)}")
    print(f"  min_graph_seg_sup  : {_int('min_graph_segment_support', 10)}")

print("=" * W)

# ─────────────────────────────────────────────────────────────────────────────
# Dispatch
# ─────────────────────────────────────────────────────────────────────────────
try:
    # ── sv_analysis ──────────────────────────────────────────────────────────
    if mode == "sv_analysis":
        results = sv.batch_microhomology_search(
            bam_path  = bam_path,
            ref_path  = ref_path,
            save_dir  = save_dir,
            # clip discovery
            min_clip_len       = _int("min_clip_len", 50),
            min_support        = _int("min_support", 5),
            min_support_frac   = _float("min_support_frac", 0.05),
            cluster_dist       = _int("cluster_dist", 10),
            top_n              = _int("top_n", 20),
            min_mapq           = _int("min_mapq", 60),
            # pairing
            max_dist           = max_dist,
            pair_strategy      = _str("pair_strategy", "split_read_bridge"),
            # sv_analysis-specific
            flank                   = _int("flank", 600),
            min_mh                  = _int("min_mh", 2),
            show_pileup             = _bool("show_pileup", False),
            show_homozygous_variants= _bool("show_homozygous_variants", False),
        )
        if results is not None and not results.empty:
            out_csv = os.path.join(save_dir, "batch_mh_results.csv")
            results.to_csv(out_csv, index=False)
            print(f"\nDone — {len(results)} pair(s) found.")
            print(f"Results → {out_csv}")
        else:
            print("\nDone — no MH pairs detected.")

    # ── sv_analysis_with_sc_adjustment ───────────────────────────────────────
    elif mode == "sv_analysis_with_sc_adjustment":
        results = sv.batch_microhomology_search(
            bam_path  = bam_path,
            ref_path  = ref_path,
            save_dir  = save_dir,
            # clip discovery
            min_clip_len       = _int("min_clip_len", 50),
            min_support        = _int("min_support", 5),
            min_support_frac   = _float("min_support_frac", 0.05),
            cluster_dist       = _int("cluster_dist", 10),
            top_n              = _int("top_n", 20),
            min_mapq           = _int("min_mapq", 60),
            # pairing
            max_dist           = max_dist,
            # sv_analysis_with_sc_adjustment-specific
            flank                   = _int("flank", 600),
            min_mh                  = _int("min_mh", 2),
            show_pileup             = _bool("show_pileup", False),
            show_homozygous_variants= _bool("show_homozygous_variants", False),
            sa_share_threshold      = _float("sa_share_threshold", 0.50),
            bridge_tol              = _int("bridge_tol", 200),
        )
        if results is not None and not results.empty:
            out_csv = os.path.join(save_dir, "batch_mh_results.csv")
            results.to_csv(out_csv, index=False)
            print(f"\nDone — {len(results)} cluster(s) found.")
            print(f"Results → {out_csv}")
        else:
            print("\nDone — no MH clusters detected.")

    # ── sv_rearrangement_path ─────────────────────────────────────────────────
    elif mode == "sv_rearrangement_path":
        results = sv.batch_rearrangement_path_search(
            bam_path  = bam_path,
            ref_path  = ref_path,
            save_dir  = save_dir,
            # clip discovery
            min_clip_len       = _int("min_clip_len", 50),
            min_support        = _int("min_support", 5),
            min_support_frac   = _float("min_support_frac", 0.05),
            cluster_dist       = _int("cluster_dist", 10),
            top_n              = _int("top_n", 20),
            min_mapq           = _int("min_mapq", 60),
            # pairing
            max_dist           = max_dist,
            pair_strategy      = _str("pair_strategy", "split_read_bridge"),
            # sv_rearrangement_path-specific
            sa_share_threshold      = _float("sa_share_threshold", 0.50),
            landing_tol             = _int("landing_tol", 200),
            fetch_padding           = _int("fetch_padding", 500_000),
            max_aln_span            = _int("max_aln_span", 1_000_000),
            minimap2_preset         = _str("minimap2_preset", "map-hifi"),
            min_identity            = _float("min_identity", 0.90),
            match_score             = _int("match_score", 2),
            mismatch_penalty        = _int("mismatch_penalty", 8),
            gap_open                = _int("gap_open", 4),
            gap_extend              = _int("gap_extend", 2),
            gap_open2               = _int("gap_open2", 24),
            gap_extend2             = _int("gap_extend2", 1),
            pdv_min_split_vaf       = _float("pdv_min_split_vaf", 0.20),
            pdv_min_enrichment      = _float("pdv_min_enrichment", 3.0),
            background_hom_vaf      = _float("background_hom_vaf", 0.95),
            min_graph_segment_support = _int("min_graph_segment_support", 10),
        )
        if results is not None and not results.empty:
            out_csv = os.path.join(save_dir, "batch_rearrangement_results.csv")
            results.to_csv(out_csv, index=False)
            print(f"\nDone — {len(results)} cluster(s) found.")
            print(f"Results → {out_csv}")
        else:
            print("\nDone — no SV clusters detected.")

except Exception:
    traceback.print_exc()
    sys.exit(1)
