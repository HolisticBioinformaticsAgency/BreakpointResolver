"""
sv_analysis_sc_adjustment_arc.py
────────────────────────────────
Integrated soft-clip-adjustment microhomology (MH) and rearrangement-arc
analysis for long-read SV breakpoint clusters.

This module runs sv_analysis_with_sc_adjustment.batch_microhomology_search()
UNCHANGED (scan cache, contig filter, CIGAR-D deletion reads, clustering,
consensus, MMEJ MH, logging) and plugs an arc-diagram step into it through
its per_cluster_hook argument:

    batch_microhomology_search()
        Stages A-D  cached_find_major_clip_pairs → _filter_pairs_by_contig
        Stage  E    cluster_clip_pairs_by_reads
        Stage  F    for EACH cluster:
                      1. analyze_sv_cluster()        (MH outputs)
                      2. per_cluster_hook → run_arc_for_sc_cluster()
                         (SNP-node graph + interactive arc HTML, from
                          sv_rearrangement_arc.py)

MH and arc outputs are therefore generated for the EXACT SAME clusters and
share the same pair_label (mh_<label>.png  ↔  arc_diagrams/<label>_arc.html).

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

Outputs (in save_dir)
─────────────────────
    batch_mh_results_<bam>.csv   MH rows + arc_* columns (written by the sc module)
    arc_clusters_summary.csv     one row per cluster, incl. clusters with no MH row
    arc_diagrams/<label>_arc.html
"""

import os
import sys
import importlib.util

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


def _arc_summary_row(cluster, pair_label, arc_result):
    """One arc summary row for a cluster.  Coordinates are 1-based inclusive,
    matching the MH CSV (refined_left + 1 .. refined_right)."""
    graph = arc_result["graph"]
    row = {
        "pair_label":                  pair_label,
        "chrom":                       cluster["chrom"],
        "refined_left":                cluster["refined_left"] + 1,
        "refined_right":               cluster["refined_right"],
        "cluster_support":             cluster.get("support"),
        "n_cluster_pairs":             len(cluster.get("pairs", [])),
        "arc_html":                    arc_result["html"],
        "arc_error":                   arc_result["error"],
        "arc_n_split_reads":           arc_result["n_split_reads"],
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

    The cluster comes from cluster_clip_pairs_by_reads in
    sv_analysis_with_sc_adjustment.py (0-based half-open refined_left/right),
    so no duplicate discovery/reclustering occurs here.

    Returns
    -------
    dict with keys: graph, html, n_split_reads, error.
    """
    arc = _load_local_module("_sv_rearrangement_arc", "sv_rearrangement_arc.py")

    chrom = cluster["chrom"]
    left  = int(cluster["refined_left"])
    right = int(cluster["refined_right"])
    arc_dir = os.path.join(save_dir, "arc_diagrams")
    os.makedirs(arc_dir, exist_ok=True)

    print(f"    [arc] Extracting split reads: {chrom}:{left + 1:,}–{right:,} "
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

    html_path = os.path.join(arc_dir, f"{pair_label}_arc.html")
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
    min_del_len=50,
    del_match_tol=20,
    edge_buffer=None,
    allowed_chroms=None,
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
    """Run the SC-adjustment MH pipeline plus an arc diagram per cluster.

    Every cluster receives MH analysis.  `arc_max_clusters=None` means every
    cluster also receives an arc diagram; an integer N limits arcs to the
    top-N clusters (by support) for short debugging runs.

    Returns
    -------
    pandas.DataFrame
        The batch_microhomology_search result (one row per MH-analysed
        cluster) with arc_* summary columns.  arc_clusters_summary.csv is
        also written to save_dir, including clusters with no MH row.
    """
    sc = _load_local_module(
        "_sv_analysis_with_sc_adjustment_integrated",
        "sv_analysis_with_sc_adjustment.py",
    )

    if save_dir is None:
        bam_stem = os.path.splitext(os.path.basename(bam_path))[0]
        save_dir = os.path.join("results", bam_stem, "sv_analysis_sc_adjustment_arc")
    os.makedirs(save_dir, exist_ok=True)

    arc_rows = []
    fasta = pysam.FastaFile(ref_path)

    def _arc_hook(index, cluster, pair_label, mh_row):
        if arc_max_clusters is not None and index > int(arc_max_clusters):
            return None
        # Arc output is attempted even if the MH path produced no row: that
        # is often diagnostically useful for a complex or misclassified event.
        print("    [ARC] SNP-node graph + interactive HTML", flush=True)
        try:
            arc_result = run_arc_for_sc_cluster(
                bam_path=bam_path,
                fasta=fasta,
                cluster=cluster,
                pair_label=pair_label,
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

        row = _arc_summary_row(cluster, pair_label, arc_result)
        row["has_mh_row"] = mh_row is not None
        arc_rows.append(row)
        return {k: v for k, v in row.items() if k.startswith("arc_")}

    try:
        results = sc.batch_microhomology_search(
            bam_path,
            ref_path,
            flank                    = flank,
            min_mh                   = min_mh,
            show_pileup              = show_pileup,
            min_clip_len             = min_clip_len,
            min_support              = min_support,
            min_support_frac         = min_support_frac,
            cluster_dist             = cluster_dist,
            min_mapq                 = min_mapq,
            max_dist                 = max_dist,
            top_n                    = top_n,
            save_dir                 = save_dir,
            show_homozygous_variants = show_homozygous_variants,
            bridge_tol               = bridge_tol,
            min_del_len              = min_del_len,
            del_match_tol            = del_match_tol,
            edge_buffer              = edge_buffer,
            allowed_chroms           = allowed_chroms,
            per_cluster_hook         = _arc_hook,
        )
    finally:
        fasta.close()

    arc_summary_path = os.path.join(save_dir, "arc_clusters_summary.csv")
    pd.DataFrame(arc_rows).to_csv(arc_summary_path, index=False)
    print(f"  Arc cluster index ({len(arc_rows)} cluster(s)): {arc_summary_path}",
          flush=True)
    return results
