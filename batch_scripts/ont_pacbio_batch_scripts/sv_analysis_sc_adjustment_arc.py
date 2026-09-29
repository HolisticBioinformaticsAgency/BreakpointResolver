"""
sv_analysis_sc_adjustment_arc.py
────────────────────────────────
Integrated soft-clip-adjustment microhomology (MH) and rearrangement-arc
analysis for long-read SV breakpoint clusters.

This module deliberately leaves sv_analysis_with_sc_adjustment.py untouched.
It reuses that module's SINGLE-PASS discovery and read-centric clustering as
the one source of truth for candidate clusters:

    find_major_clip_pairs()
        → cluster_clip_pairs_by_reads()
            → for EACH resulting cluster:
                1. analyze_sv_cluster()
                   (flank consensus, optional pileup PNGs, MMEJ MH call,
                    microhomology structure PNG)
                2. build_snp_segment_graph() + plot_rearrangement_arc_diagram()
                   (SNP-node rearrangement graph + interactive arc HTML)

The central design choice is that MH and arc outputs are generated for the
EXACT SAME cluster coordinates.  We do NOT call run_pipeline() from
sv_analysis_junction_rearr.py because it independently rediscovers and
reclusters breakpoint pairs, which can yield a different candidate set.

Usage
─────
Called by run_mh_worker.py when params_batch.yaml contains:

    mode: "sv_analysis_sc_adjustment_arc"

Direct Python usage:

    import sv_analysis_sc_adjustment_arc as integrated
    results = integrated.batch_microhomology_arc_search(
        bam_path=BAM,
        ref_path=REF,
        save_dir="results/sample",
    )

Project constraint
──────────────────
Do not touch sv_analysis mode.  This module adds a separate integration mode
and does not modify sv_analysis.py or sv_analysis_with_sc_adjustment.py.
"""

import os
import sys
import importlib.util
from datetime import datetime

import pandas as pd
import pysam


_MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


def _load_local_module(module_name, filename):
    """Load a sibling Python module by file path without relying on PYTHONPATH."""
    path = os.path.join(_MODULE_DIR, filename)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Required module file not found: {path}")
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _banner(stage, message):
    sep = "=" * 76
    print(f"\n{sep}", flush=True)
    print(f"[{datetime.now().strftime('%H:%M:%S')}]  {stage}  |  {message}", flush=True)
    print(sep, flush=True)


def _arc_summary_row(cluster, graph, html_path, mh_row, arc_error=None):
    """Build one stable summary row joining the MH result to its arc result."""
    row = {
        "pair_label":                  mh_row.get("pair_label") if mh_row else None,
        "chrom":                       cluster["chrom"],
        "refined_left":                cluster["refined_left"],
        "refined_right":               cluster["refined_right"],
        "cluster_support":             cluster.get("support"),
        "n_cluster_pairs":             len(cluster.get("pairs", [])),
        "arc_html":                    html_path,
        "arc_error":                   arc_error,
        "arc_n_nodes":                 None,
        "arc_n_edges":                 None,
        "arc_n_sv_informative_nodes":  None,
        "arc_n_paralogous_nodes":      None,
        "arc_n_germline_nodes":        None,
        "arc_n_junction_edges":        None,
        "arc_n_inversion_edges":       None,
    }
    if graph is not None:
        nodes = list(graph.nodes(data=True))
        edges = list(graph.edges(data=True))
        row.update({
            "arc_n_nodes":                graph.number_of_nodes(),
            "arc_n_edges":                graph.number_of_edges(),
            "arc_n_sv_informative_nodes": sum(
                1 for _, d in nodes if d.get("sv_informative")
            ),
            "arc_n_paralogous_nodes": sum(
                1 for _, d in nodes if d.get("context_paralogous")
            ),
            "arc_n_germline_nodes": sum(
                1 for _, d in nodes if d.get("germline")
            ),
            "arc_n_junction_edges": sum(
                1 for _, _, d in edges if d.get("junction_weight", 0) > 0
            ),
            "arc_n_inversion_edges": sum(
                1 for _, _, d in edges if d.get("inversion_weight", 0) > 0
            ),
        })
    return row


def run_arc_for_sc_cluster(
    bam_path,
    fasta,
    cluster,
    pair_label,
    save_dir,
    arc_padding=500,
    arc_min_perfect_run=10,
    arc_min_non_split_support=3,
    arc_min_non_split_frac=0.1,
    arc_supp_mapq_min=0,
    arc_min_node_support=2,
    arc_context_flank=50,
    arc_context_min_identity=1.0,
    arc_context_search_revcomp=True,
):
    """Build an arc diagram for one pre-existing SC-adjustment cluster.

    The cluster is produced by cluster_clip_pairs_by_reads in
    sv_analysis_with_sc_adjustment.py.  It is therefore exactly the same
    (chrom, refined_left, refined_right) cluster that receives MH consensus
    analysis — no duplicate discovery/reclustering occurs here.

    Returns
    -------
    dict with keys: graph, html, n_split_reads, error.
    """
    arc = _load_local_module(
        "_sv_analysis_junction_rearr_integrated",
        "sv_analysis_junction_rearr.py",
    )

    chrom = cluster["chrom"]
    left  = int(cluster["refined_left"])
    right = int(cluster["refined_right"])
    arc_dir = os.path.join(save_dir, "arc_diagrams")
    os.makedirs(arc_dir, exist_ok=True)

    print(f"    [arc] Extracting split reads: {chrom}:{left:,}–{right:,} "
          f"(padding={arc_padding})", flush=True)
    reads = arc.extract_split_reads(
        bam_path,
        chrom,
        left,
        right,
        padding=arc_padding,
    )
    print(f"    [arc] extract_split_reads → {len(reads)} read(s)", flush=True)

    if not reads:
        return {
            "graph": None,
            "html": None,
            "n_split_reads": 0,
            "error": "No split reads extracted for SC-adjustment cluster",
        }

    graph = arc.build_snp_segment_graph(
        reads,
        fasta,
        bam_path,
        min_perfect_run=arc_min_perfect_run,
        min_non_split_support=arc_min_non_split_support,
        min_non_split_frac=arc_min_non_split_frac,
        supp_mapq_min=arc_supp_mapq_min,
        min_node_support=arc_min_node_support,
        context_flank=arc_context_flank,
        context_min_identity=arc_context_min_identity,
        context_search_revcomp=arc_context_search_revcomp,
    )

    if graph.number_of_nodes() == 0:
        return {
            "graph": graph,
            "html": None,
            "n_split_reads": len(reads),
            "error": "Arc graph contains no nodes after decomposition/pruning",
        }

    html_path = os.path.join(
        arc_dir,
        f"{pair_label}_arc_{chrom}_{left}_{right}.html",
    )
    arc.plot_rearrangement_arc_diagram(graph, save_path=html_path)
    return {
        "graph": graph,
        "html": html_path,
        "n_split_reads": len(reads),
        "error": None,
    }


def batch_microhomology_arc_search(
    bam_path,
    ref_path,
    flank=None,
    min_mh=2,
    show_pileup=False,
    min_clip_len=50,
    min_support=5,
    min_support_frac=0.05,
    cluster_dist=10,
    min_mapq=60,
    max_dist=None,
    top_n=50,
    save_dir=None,
    show_homozygous_variants=False,
    bridge_tol=200,
    arc_padding=500,
    arc_min_perfect_run=10,
    arc_min_non_split_support=3,
    arc_min_non_split_frac=0.1,
    arc_supp_mapq_min=0,
    arc_min_node_support=2,
    arc_context_flank=50,
    arc_context_min_identity=1.0,
    arc_context_search_revcomp=True,
    arc_max_clusters=None,
):
    """Run MH plus an interactive arc diagram for every SC-adjustment cluster.

    Pipeline
    --------
    A-D : sc.find_major_clip_pairs()
          Single-pass SA read scan, clip-pair binning, pair support tallying,
          and support filtering.
    E   : sc.cluster_clip_pairs_by_reads()
          Read-centric Union-Find clustering.
    F   : For every retained SC-adjustment cluster:
          (a) sc.analyze_sv_cluster() → consensus, pileup PNG, MH result/PNG.
          (b) run_arc_for_sc_cluster() → SNP-node graph + arc HTML.

    `arc_max_clusters=None` means ALL clusters receive an arc diagram.  It is
    present solely for debugging; use an integer during short test jobs.

    Returns
    -------
    pandas.DataFrame
        One row per successfully MH-analyzed cluster, augmented with arc graph
        summary metrics and the path to the corresponding HTML (if generated).
        A CSV is written to save_dir/batch_mh_arc_results.csv.
    """
    sc = _load_local_module(
        "_sv_analysis_with_sc_adjustment_integrated",
        "sv_analysis_with_sc_adjustment.py",
    )

    if save_dir is None:
        bam_stem = os.path.splitext(os.path.basename(bam_path))[0]
        save_dir = os.path.join("results", bam_stem, "sv_analysis_sc_adjustment_arc")
    os.makedirs(save_dir, exist_ok=True)

    _banner("PIPELINE START", "sv_analysis_sc_adjustment_arc")
    print(f"  BAM      : {bam_path}", flush=True)
    print(f"  Ref      : {ref_path}", flush=True)
    print(f"  Out dir  : {save_dir}", flush=True)
    print("  Strategy : SC-adjustment clusters are the shared source of truth; "
          "each cluster receives both MH and arc analysis.", flush=True)

    _banner("STAGES A-D", "SC-adjustment single-pass clip-pair discovery")
    pairs_df = sc.find_major_clip_pairs(
        bam_path,
        ref_path,
        chrom=None,
        min_clip_len=min_clip_len,
        min_support=min_support,
        min_support_frac=min_support_frac,
        cluster_dist=cluster_dist,
        top_n=top_n,
        min_mapq=min_mapq,
        max_dist=max_dist,
        save_dir=save_dir,
    )
    if pairs_df.empty:
        print("  No clip pairs found. Exiting.", flush=True)
        return pd.DataFrame()

    _banner("STAGE E", "SC-adjustment read-centric Union-Find clustering")
    clusters = sc.cluster_clip_pairs_by_reads(
        pairs_df,
        bam_path,
        cluster_dist=cluster_dist,
        min_mapq=min_mapq,
        bridge_tol=bridge_tol,
    )
    if not clusters:
        print("  No clusters found. Exiting.", flush=True)
        return pd.DataFrame()

    all_cluster_count = len(clusters)
    if arc_max_clusters is not None:
        clusters = clusters[:int(arc_max_clusters)]
        print(f"  DEBUG CAP: processing {len(clusters)}/{all_cluster_count} cluster(s) "
              f"because arc_max_clusters={arc_max_clusters}", flush=True)
    else:
        print(f"  Processing ALL {all_cluster_count} SC-adjustment cluster(s).", flush=True)

    _banner("STAGE F", "Integrated MH + rearrangement arc analysis")
    mh_rows = []
    arc_rows = []

    # The arc module's context cache is module-global.  Opening the FASTA once
    # here avoids per-cluster FASTA setup; its window cache itself remains
    # bounded by the implementation in sv_analysis_junction_rearr.py.
    fasta = pysam.FastaFile(ref_path)
    try:
        for i, cluster in enumerate(clusters, 1):
            label = f"c{i:02d}_{cluster['chrom']}_{cluster['refined_left']}"
            chrom = cluster["chrom"]
            left = cluster["refined_left"]
            right = cluster["refined_right"]

            print(f"\n  Cluster {i}/{len(clusters)}: {chrom}:{left:,}–{right:,} "
                  f"(support={cluster.get('support')}, "
                  f"pairs={len(cluster.get('pairs', []))})", flush=True)

            # MH output.  analyze_sv_cluster opens its own FASTA internally,
            # preserving the existing SC-adjustment implementation unchanged.
            print("    [MH] consensus + microhomology analysis", flush=True)
            mh_row = sc.analyze_sv_cluster(
                bam_path,
                ref_path,
                cluster,
                flank=flank,
                min_mh=min_mh,
                show_pileup=show_pileup,
                min_mapq=min_mapq,
                pair_label=label,
                save_dir=save_dir,
                show_homozygous_variants=show_homozygous_variants,
            )

            # Arc output is intentionally attempted even if the MH path could
            # not produce a result: that is often diagnostically useful for a
            # complex or misclassified event.
            print("    [ARC] SNP-node graph + interactive HTML", flush=True)
            try:
                arc_result = run_arc_for_sc_cluster(
                    bam_path=bam_path,
                    fasta=fasta,
                    cluster=cluster,
                    pair_label=label,
                    save_dir=save_dir,
                    arc_padding=arc_padding,
                    arc_min_perfect_run=arc_min_perfect_run,
                    arc_min_non_split_support=arc_min_non_split_support,
                    arc_min_non_split_frac=arc_min_non_split_frac,
                    arc_supp_mapq_min=arc_supp_mapq_min,
                    arc_min_node_support=arc_min_node_support,
                    arc_context_flank=arc_context_flank,
                    arc_context_min_identity=arc_context_min_identity,
                    arc_context_search_revcomp=arc_context_search_revcomp,
                )
            except Exception as exc:
                # Continue processing remaining clusters: a single malformed
                # alignment/edge case should not discard every MH result.
                arc_result = {
                    "graph": None,
                    "html": None,
                    "n_split_reads": None,
                    "error": f"{type(exc).__name__}: {exc}",
                }
                print(f"    [ARC] WARNING: {arc_result['error']}", flush=True)

            joined_row = _arc_summary_row(
                cluster,
                arc_result["graph"],
                arc_result["html"],
                mh_row,
                arc_error=arc_result["error"],
            )
            joined_row["arc_n_split_reads"] = arc_result["n_split_reads"]
            arc_rows.append(joined_row)

            if mh_row is not None:
                merged = dict(mh_row)
                merged.update({
                    key: value for key, value in joined_row.items()
                    if key not in merged or key.startswith("arc_")
                })
                mh_rows.append(merged)
            else:
                print("    [MH] No MH result row; arc result retained in "
                      "arc_clusters_summary.csv.", flush=True)

    finally:
        fasta.close()

    arc_summary = pd.DataFrame(arc_rows)
    arc_summary_path = os.path.join(save_dir, "arc_clusters_summary.csv")
    arc_summary.to_csv(arc_summary_path, index=False)

    results = pd.DataFrame(mh_rows)
    result_path = os.path.join(save_dir, "batch_mh_arc_results.csv")
    if not results.empty:
        results.to_csv(result_path, index=False)
        _banner("PIPELINE COMPLETE",
                f"MH + arc results for {len(results)} cluster(s) → {result_path}")
        try:
            from IPython.display import display
            display(results)
        except ImportError:
            print(results.to_string(index=False))
    else:
        _banner("PIPELINE COMPLETE",
                "No MH rows returned; arc-only status retained in "
                f"{arc_summary_path}")

    print(f"  Arc cluster index: {arc_summary_path}", flush=True)
    return results
