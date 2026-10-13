#!/usr/bin/env python3
"""
collate_mh_results.py
Collates the per-sample outputs of run_mh_worker.py into single tables, and
gathers the top-cluster arc diagrams into one folder.

Usage:
    python3 collate_mh_results.py [--mode MODE] [--results-dir results]
                                  [--match-tol 50] [--top-n 1]
                                  [--top-n-override 1113_ONT=2 ...]

Inputs (one set per barcode × platform):
    results/<bc>_<platform>/<bam_stem>/<mode>/batch_mh_results_<bam_stem>.csv
    results/<bc>_<platform>/<bam_stem>/<mode>/arc_clusters_summary.csv
    results/<bc>_<platform>/<bam_stem>/<mode>/arc_diagrams/<label>_arc.html

Outputs (in results/):
    collated_batch_mh_results_<mode>.csv    top-N MH rows per BAM (default 1;
                                            1113_ONT: 2)
    collated_batch_mh_results_all_clusters_<mode>.csv
                                            every MH row, all samples
    collated_arc_clusters_summary_<mode>.csv every arc summary row, all samples
    collated_sample_status_<mode>.csv       one row per task in mh_tasks.txt
    collated_arc_html_<mode>/               top-N arc HTMLs per BAM (same N),
                                            index.csv and index.html

"Top" = lowest cluster index (c01 = highest-support cluster), i.e. the first
rows of each per-sample batch_mh_results CSV.

Both collated tables get barcode / platform / bam_stem / source_csv columns
and an ONT-vs-PacBio event match:
    event_id            rows of the same barcode on the same chrom whose
                        refined_left AND refined_right are both within
                        --match-tol bp are linked (transitively) into one event
    platforms_detected  "ONT+PacBio", "ONT" or "PacBio" for that event
Coordinates in both tables are 1-based, so they are compared directly.
"""
import os, re, sys, glob, html, argparse

import pandas as pd


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Samples that get more than --top-n MH rows / arc diagrams collated.
DEFAULT_TOP_N_OVERRIDES = {"1113_ONT": 2}

ID_COLS = ["barcode", "sample_id", "platform", "bam_stem", "event_id", "platforms_detected"]


def load_sample_ids(path):
    """{barcode: sample_id} from a two-column TSV (barcode, sample_id); {} if absent."""
    if not path or not os.path.isfile(path):
        print(f"WARNING: sample ID table not found — {path}; sample_id left blank")
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(df["barcode"].str.strip(), df["sample_id"].str.strip()))


def add_sample_ids(df, sample_ids):
    """Insert a sample_id column right after barcode."""
    if df.empty:
        return df
    df = df.copy()
    if "sample_id" in df.columns:
        df = df.drop(columns="sample_id")
    df.insert(df.columns.get_loc("barcode") + 1, "sample_id",
              df["barcode"].astype(str).map(sample_ids))
    return df


def default_mode():
    """Mode from params_batch.yaml, falling back to the worker's default."""
    try:
        import yaml
        with open(os.path.join(SCRIPT_DIR, "params_batch.yaml")) as fh:
            return yaml.safe_load(fh).get("mode") or "sv_analysis_sc_adjustment_arc"
    except Exception:
        return "sv_analysis_sc_adjustment_arc"


def cluster_index(pair_label):
    """Cluster rank from a pair label like '<bam>_c03_chr17_43045093' → 3."""
    m = re.search(r"_c(\d+)_", str(pair_label))
    return int(m.group(1)) if m else None


def short_label(pair_label):
    """Label without the BAM stem: '<bam>_c03_chr17_43045093' → 'c03_chr17_43045093'."""
    m = re.search(r"_(c\d+_.+)$", str(pair_label))
    return m.group(1) if m else str(pair_label)


def find_sample_dirs(results_dir, mode):
    """{(barcode, platform): (bam_stem, mode_dir)} for every sample output folder."""
    out = {}
    for mode_dir in sorted(glob.glob(os.path.join(results_dir, "*", "*", mode))):
        sample = os.path.basename(os.path.dirname(os.path.dirname(mode_dir)))
        if "_" not in sample:
            continue
        barcode, platform = sample.split("_", 1)
        bam_stem = os.path.basename(os.path.dirname(mode_dir))
        out[(barcode, platform)] = (bam_stem, mode_dir)
    return out


def read_csv_or_empty(path):
    if not os.path.isfile(path):
        return None
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def collate(sample_dirs, filename_fn):
    """Stack one CSV per sample, adding barcode / platform / bam_stem / source_csv."""
    frames = []
    for (barcode, platform), (bam_stem, mode_dir) in sample_dirs.items():
        path = os.path.join(mode_dir, filename_fn(bam_stem))
        df = read_csv_or_empty(path)
        if df is None or df.empty:
            continue
        df.insert(0, "barcode", barcode)
        df.insert(1, "platform", platform)
        df.insert(2, "bam_stem", bam_stem)
        df["source_csv"] = path
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    df["cluster_index"] = df["pair_label"].map(cluster_index)
    return df


def assign_events(df, tol):
    """Add event_id / platforms_detected by linking rows of the same barcode and
    chrom whose refined_left and refined_right are both within *tol* bp."""
    if df.empty:
        return df
    df = df.copy()
    parent = list(range(len(df)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for _, grp in df.groupby(["barcode", "chrom"]):
        idx = grp.index.tolist()
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                i, j = idx[a], idx[b]
                if (abs(df.at[i, "refined_left"] - df.at[j, "refined_left"]) <= tol and
                        abs(df.at[i, "refined_right"] - df.at[j, "refined_right"]) <= tol):
                    parent[find(i)] = find(j)

    df["_root"] = [find(i) for i in range(len(df))]
    # Number events per barcode in genomic order.
    first = (df.groupby("_root")
               .agg(barcode=("barcode", "first"), chrom=("chrom", "first"),
                    left=("refined_left", "min"))
               .sort_values(["barcode", "chrom", "left"]))
    first["n"] = first.groupby("barcode").cumcount() + 1
    event_ids = {root: f"{r.barcode}_E{r.n:02d}" for root, r in first.iterrows()}
    df["event_id"] = df["_root"].map(event_ids)

    plats = df.groupby("event_id")["platform"].agg(lambda s: set(s))
    label = {e: ("ONT+PacBio" if {"ONT", "PacBio"} <= p else "+".join(sorted(p)))
             for e, p in plats.items()}
    df["platforms_detected"] = df["event_id"].map(label)

    df = df.drop(columns="_root")
    df = df.sort_values(["barcode", "platform", "cluster_index"]).reset_index(drop=True)
    front = [c for c in ID_COLS if c in df.columns]
    return df[front + [c for c in df.columns if c not in front]]


def read_tasks(results_dir):
    """(barcode, platform) list from mh_tasks.txt, or None if it is missing."""
    path = os.path.join(SCRIPT_DIR, "mh_tasks.txt")
    if not os.path.isfile(path):
        return None
    tasks = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                tasks.append((parts[0], parts[1]))
    return tasks


def sample_status(tasks, sample_dirs, mh):
    """One row per task: ok / no_mh_rows / incomplete / missing."""
    n_rows = mh.groupby(["barcode", "platform"]).size().to_dict() if not mh.empty else {}
    rows = []
    for barcode, platform in tasks:
        key = (barcode, platform)
        n = n_rows.get(key, 0)
        if key not in sample_dirs:
            status = "missing"
        elif n > 0:
            status = "ok"
        else:
            bam_stem, mode_dir = sample_dirs[key]
            log = os.path.join(mode_dir, f"log_{bam_stem}.txt")
            text = open(log).read() if os.path.isfile(log) else ""
            # The sc pipeline prints "... Exiting." when it stops early with
            # no pairs/clusters, and "PIPELINE COMPLETE" when it finishes.
            finished = "PIPELINE COMPLETE" in text or "Exiting." in text
            status = "no_mh_rows" if finished else "incomplete"
        rows.append({"barcode": barcode, "platform": platform,
                     "status": status, "n_mh_rows": n})
    return pd.DataFrame(rows)


def collect_arc_html(arc, sample_dirs, out_dir, top_n, overrides, sample_ids):
    """Copy the top-N arc HTMLs per sample into out_dir; return the index table."""
    os.makedirs(out_dir, exist_ok=True)
    # Clear HTMLs from a previous collation so removed/renamed ones don't linger.
    for old in glob.glob(os.path.join(out_dir, "*_arc.html")):
        os.remove(old)

    rows = []
    for (barcode, platform) in sorted(sample_dirs):
        sample = f"{barcode}_{platform}"
        n = overrides.get(sample, top_n)
        sub = (arc[(arc["barcode"] == barcode) & (arc["platform"] == platform)]
               if not arc.empty else arc)
        for rank in range(1, n + 1):
            hit = sub[sub["cluster_index"] == rank] if not sub.empty else sub
            row = {"sample": sample, "barcode": barcode,
                   "sample_id": sample_ids.get(barcode), "platform": platform,
                   "cluster_index": rank, "pair_label": None, "event_id": None,
                   "platforms_detected": None, "collated_html": None, "note": None}
            if hit.empty:
                row["note"] = "no cluster at this rank"
                rows.append(row)
                continue
            r = hit.iloc[0]
            row.update(pair_label=r["pair_label"], event_id=r["event_id"],
                       platforms_detected=r["platforms_detected"])
            src = r.get("arc_html")
            if not isinstance(src, str) or not os.path.isfile(src):
                row["note"] = (r.get("arc_error") if isinstance(r.get("arc_error"), str)
                               else "arc HTML not found")
            else:
                sid = sample_ids.get(barcode) or "unknown"
                lbl = short_label(r["pair_label"])
                name = f"{barcode}_{platform}_{re.sub(r'[^A-Za-z0-9]+', '_', sid)}_{lbl}_arc.html"
                copy_annotated_html(src, os.path.join(out_dir, name),
                                    f"{sid} · {barcode} {platform} · {lbl}")
                row["collated_html"] = name
            rows.append(row)

    index = pd.DataFrame(rows)
    index.to_csv(os.path.join(out_dir, "index.csv"), index=False)
    write_index_html(index, os.path.join(out_dir, "index.html"))
    return index


def copy_annotated_html(src, dst, caption):
    """Copy an arc HTML, adding *caption* (sample ID, barcode, cluster) to its
    <title> and <h1> so each diagram identifies its sample when opened."""
    with open(src, encoding="utf-8") as fh:
        text = fh.read()
    cap = html.escape(caption)
    text = re.sub(r"<title>(.*?)</title>", lambda m: f"<title>{cap} · {m.group(1)}</title>",
                  text, count=1)
    text = re.sub(r"<h1>(.*?)</h1>", lambda m: f"<h1>{m.group(1)} — {cap}</h1>",
                  text, count=1)
    with open(dst, "w", encoding="utf-8") as fh:
        fh.write(text)


def write_index_html(index, path):
    """Plain table linking every collated HTML, for browsing the folder."""
    body = []
    for r in index.itertuples():
        link = (f'<a href="{html.escape(r.collated_html)}">{html.escape(r.collated_html)}</a>'
                if isinstance(r.collated_html, str) else "")
        body.append(
            "<tr>" + "".join(f"<td>{c}</td>" for c in [
                html.escape(r.sample), html.escape(str(r.sample_id or "")), r.cluster_index,
                html.escape(str(r.event_id or "")),
                html.escape(str(r.platforms_detected or "")),
                link, html.escape(str(r.note or "")),
            ]) + "</tr>"
        )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            "<!doctype html><meta charset='utf-8'><title>Arc diagrams</title>"
            "<style>body{font-family:sans-serif;margin:16px}"
            "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:4px 8px;"
            "text-align:left}th{background:#eee}</style>"
            "<h1>Top-cluster arc diagrams</h1><table><tr><th>Sample</th><th>Sample ID</th><th>Cluster</th>"
            "<th>Event</th><th>Platforms</th><th>Arc diagram</th><th>Note</th></tr>"
            + "\n".join(body) + "</table>"
        )


def top_rows(df, top_n, overrides):
    """First top_n rows (by cluster index) per BAM, with per-sample overrides."""
    if df.empty:
        return df
    n = (df["barcode"].astype(str) + "_" + df["platform"]).map(
        lambda sample: overrides.get(sample, top_n))
    rank = df.groupby(["barcode", "platform"])["cluster_index"].rank(method="first")
    return df[rank <= n].reset_index(drop=True)


def parse_overrides(items):
    out = dict(DEFAULT_TOP_N_OVERRIDES)
    for item in items or []:
        sample, n = item.split("=")
        out[sample] = int(n)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", default=default_mode())
    ap.add_argument("--results-dir", default=os.path.join(SCRIPT_DIR, "results"))
    ap.add_argument("--sample-ids", default=os.path.join(SCRIPT_DIR, "barcode_sample_ids.tsv"),
                    help="TSV with columns barcode, sample_id")
    ap.add_argument("--match-tol", type=int, default=50,
                    help="bp tolerance on both breakpoints for ONT/PacBio event matching")
    ap.add_argument("--top-n", type=int, default=1,
                    help="MH rows and arc HTMLs to keep per BAM (top clusters by support)")
    ap.add_argument("--top-n-override", nargs="*", metavar="SAMPLE=N",
                    help=f"per-sample top-N, e.g. 1113_ONT=2 "
                         f"(default overrides: {DEFAULT_TOP_N_OVERRIDES})")
    args = ap.parse_args()

    mode, rdir = args.mode, args.results_dir
    sample_dirs = find_sample_dirs(rdir, mode)
    print(f"Mode        : {mode}")
    print(f"Results dir : {rdir}")
    print(f"Samples     : {len(sample_dirs)} output folder(s)")
    if not sample_dirs:
        sys.exit("No sample output folders found.")

    sample_ids = load_sample_ids(args.sample_ids)
    mh = add_sample_ids(
        assign_events(collate(sample_dirs, lambda s: f"batch_mh_results_{s}.csv"),
                      args.match_tol), sample_ids)
    arc = add_sample_ids(
        assign_events(collate(sample_dirs, lambda s: "arc_clusters_summary.csv"),
                      args.match_tol), sample_ids)

    overrides = parse_overrides(args.top_n_override)
    mh_top = top_rows(mh, args.top_n, overrides)

    mh_path = os.path.join(rdir, f"collated_batch_mh_results_{mode}.csv")
    mh_all_path = os.path.join(rdir, f"collated_batch_mh_results_all_clusters_{mode}.csv")
    arc_path = os.path.join(rdir, f"collated_arc_clusters_summary_{mode}.csv")
    mh_top.to_csv(mh_path, index=False)
    mh.to_csv(mh_all_path, index=False)
    arc.to_csv(arc_path, index=False)
    print(f"MH top rows : {len(mh_top):>4} → {mh_path}")
    print(f"MH all rows : {len(mh):>4} → {mh_all_path}")
    print(f"Arc rows    : {len(arc):>4} → {arc_path}")

    tasks = read_tasks(rdir) or sorted(sample_dirs)
    status = add_sample_ids(sample_status(tasks, sample_dirs, mh), sample_ids)
    status_path = os.path.join(rdir, f"collated_sample_status_{mode}.csv")
    status.to_csv(status_path, index=False)
    counts = status["status"].value_counts().to_dict()
    print(f"Status      : {len(status)} task(s): "
          + ", ".join(f"{v} {k}" for k, v in counts.items()) + f" → {status_path}")
    for r in status[status["status"] != "ok"].itertuples():
        print(f"              {r.barcode}_{r.platform}: {r.status}")

    if not mh.empty:
        ev = mh.drop_duplicates("event_id")["platforms_detected"].value_counts().to_dict()
        print(f"MH events   : {mh['event_id'].nunique()} (tol {args.match_tol} bp): "
              + ", ".join(f"{v} {k}" for k, v in ev.items()))

    html_dir = os.path.join(rdir, f"collated_arc_html_{mode}")
    index = collect_arc_html(arc, sample_dirs, html_dir, args.top_n, overrides,
                             sample_ids)
    n_ok = index["collated_html"].notna().sum()
    print(f"Arc HTML    : {n_ok}/{len(index)} copied → {html_dir}")
    for r in index[index["collated_html"].isna()].itertuples():
        print(f"              {r.sample} c{r.cluster_index:02d}: {r.note}")


if __name__ == "__main__":
    main()
