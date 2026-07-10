"""Genome-wide mRNA survival screen in NSCLC + pathway over-representation.

Answers two questions end to end:
  (1) which genes' expression has a *real* impact on overall survival, and
  (2) which biological pathways those genes concentrate in.

Pipeline
--------
1. Pull overall-survival clinical data for the NSCLC cohort (TCGA LUAD + LUSC).
2. Stream mRNA expression for every protein-coding gene from cBioPortal in
   chunks (log2-transformed).
3. Screen every gene with a fast, subtype-stratified Cox score test
   (:mod:`mirna_tcga.screen`) and FDR-correct.
4. Refit a full Cox model (lifelines) on the strongest hits for hazard ratios.
5. Test the FDR-significant genes for pathway over-representation
   (:mod:`mirna_tcga.enrich`) against MSigDB / Enrichr gene sets, split by
   direction (risk vs protective).

Requires the survival extra:  pip install 'mirna-tcga[survival]'
Run (full):   python scripts/09_survival_screen.py
Run (quick):  python scripts/09_survival_screen.py --max-genes 2000
"""

from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import urlopen

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, nsclc_clinical
from mirna_tcga.enrich import over_representation, parse_gmt
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.screen import cox_score_screen
from mirna_tcga.survival import coerce_clinical, cox_model

# Default self-contained gene-set libraries (reachable GitHub raw mirrors).
DEFAULT_GENE_SETS = {
    "KEGG_2016": "https://raw.githubusercontent.com/zqfang/GSEApy/master/tests/extdata/enrichr.KEGG_2016.gmt",
    "Hallmark": "https://raw.githubusercontent.com/zqfang/GSEApy/master/tests/extdata/h.all.v7.0.symbols.gmt",
}


def protein_coding_genes(client: CBioPortalClient) -> pd.DataFrame:
    """All protein-coding genes as a DataFrame[entrezGeneId, hugoGeneSymbol]."""
    genes = pd.DataFrame(client._get("genes", {"projection": "SUMMARY"}))
    genes = genes[genes["type"] == "protein-coding"]
    return genes[["entrezGeneId", "hugoGeneSymbol"]].dropna()


def fetch_expression_wide(client, profile, sample_list, entrez, id2sym, chunk=2000):
    """Stream a genes(HUGO) x samples log2 expression matrix in entrez chunks."""
    parts = []
    for i in range(0, len(entrez), chunk):
        sub = entrez[i : i + chunk]
        long = client.fetch_molecular_data(profile, sub, sample_list_id=sample_list)
        if long.empty:
            continue
        long["gene"] = long["entrezGeneId"].astype(int).map(id2sym)
        wide = long.pivot_table(index="gene", columns="sampleId", values="value", aggfunc="first")
        parts.append(wide.astype("float32"))
        del long, wide
    if not parts:
        return pd.DataFrame()
    mat = pd.concat(parts)
    return np.log2(mat.clip(lower=0) + 1.0)


def load_gene_sets(specs: dict[str, str], cache_dir: Path) -> dict[str, set[str]]:
    """Load (and disk-cache) GMT gene-set libraries from URLs or local paths."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    merged: dict[str, set[str]] = {}
    for label, src in specs.items():
        if src.startswith("http"):
            cache = cache_dir / f"{label}.gmt"
            if not cache.exists():
                cache.write_bytes(urlopen(src, timeout=60).read())  # noqa: S310
            text = cache.read_text()
        else:
            text = Path(src).read_text()
        sets = parse_gmt(text)
        # KEGG_2016 names carry a " Homo sapiens hsaNNNNN" suffix; trim it.
        for name, genes in sets.items():
            clean = name.split(" Homo sapiens ")[0].strip()
            merged[f"{label}: {clean}"] = genes
    return merged


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=None)
    parser.add_argument("--fdr", type=float, default=0.05, help="BH q threshold for 'real impact'")
    parser.add_argument("--top", type=int, default=25, help="hits to refit with full Cox")
    parser.add_argument("--max-genes", type=int, default=None, help="cap universe (for a quick run)")
    parser.add_argument("--save-dir", default=None, help="directory to write result CSVs")
    args = parser.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)

    # 1. survival outcomes for the NSCLC cohort
    clinical = nsclc_clinical(client, cfg, patient_level=True)
    surv = coerce_clinical(clinical)  # -> duration, event
    surv = surv.set_index("patientId")[["duration", "event", "cohort"]]
    print(f"NSCLC patients with overall-survival data: {len(surv)}")

    # 2. stream all protein-coding expression per study, log2 space
    gene_tbl = protein_coding_genes(client)
    if args.max_genes:
        gene_tbl = gene_tbl.head(args.max_genes)
    id2sym = dict(zip(gene_tbl["entrezGeneId"].astype(int), gene_tbl["hugoGeneSymbol"]))
    entrez = list(id2sym)
    print(f"Screening {len(entrez)} protein-coding genes across LUAD + LUSC ...")

    mats = []
    for key in cohort_study_keys(cfg, "nsclc"):
        mat = fetch_expression_wide(
            client, cfg.mrna_profile(key), cfg.all_samples_list(key), entrez, id2sym
        )
        print(f"  {key.upper()}: {mat.shape[0]} genes x {mat.shape[1]} samples")
        mats.append(mat)
    expr = pd.concat(mats, axis=1)  # genes x samples (both studies)

    # 3. align samples -> patients -> survival, then score-test screen
    X = expr.T  # samples x genes
    X["patientId"] = sample_to_patient(X.index)
    X = X[X["patientId"].isin(surv.index)]
    meta = surv.loc[X["patientId"]]
    strata = meta["cohort"].to_numpy()
    dur = pd.Series(meta["duration"].to_numpy(), index=X.index)
    evt = pd.Series(meta["event"].to_numpy(), index=X.index)
    Xg = X.drop(columns="patientId").dropna(axis=1)
    print(f"Analysis matrix: {Xg.shape[0]} samples x {Xg.shape[1]} genes "
          f"({int(evt.sum())} deaths)\n")

    res = cox_score_screen(Xg, dur, evt, pd.Series(strata, index=X.index))
    hits = res[res["q"] < args.fdr]
    risk = hits[hits["z"] > 0]
    prot = hits[hits["z"] < 0]
    print(f"Genes with real survival impact (BH q < {args.fdr}): {len(hits)}")
    print(f"  higher expression -> WORSE survival: {len(risk)}")
    print(f"  higher expression -> BETTER survival: {len(prot)}\n")
    print("Top 15 by significance:")
    print(res.head(15).assign(dir=np.where(res.head(15)["z"] > 0, "risk", "prot")).round(4), "\n")

    # 4. reportable hazard ratios via full Cox on the strongest hits
    print(f"Full Cox hazard ratios (per SD, subtype-adjusted) for top {args.top} hits:")
    top = hits.head(args.top)
    rows = []
    for gene in top.index:
        d = pd.DataFrame({
            "expr_z": (Xg[gene] - Xg[gene].mean()) / Xg[gene].std(),
            "is_lusc": (strata == "LUSC").astype(int),
            "duration": dur.to_numpy(),
            "event": evt.to_numpy(),
        }).dropna()
        try:
            cph = cox_model(d, ["expr_z", "is_lusc"])
            s = cph.summary.loc["expr_z"]
            rows.append((gene, s["exp(coef)"], s["exp(coef) lower 95%"],
                         s["exp(coef) upper 95%"], s["p"]))
        except Exception:  # pragma: no cover - numerical edge cases
            continue
    hr = pd.DataFrame(rows, columns=["gene", "HR", "HR_low", "HR_high", "cox_p"]).set_index("gene")
    print(hr.round(3), "\n")

    # 5. pathway over-representation of the significant genes
    cache = Path(args.save_dir or ".") / "genesets_cache"
    gene_sets = load_gene_sets(cfg.raw.get("gene_sets") or DEFAULT_GENE_SETS, cache)
    universe = list(Xg.columns)
    print(f"Pathway over-representation ({len(gene_sets)} sets, universe={len(universe)} genes):\n")
    ora_tables = {}
    for label, genes in (("all significant", hits.index),
                         ("risk (worse survival)", risk.index),
                         ("protective (better survival)", prot.index)):
        ora = over_representation(genes, universe, gene_sets)
        ora_tables[label] = ora
        sig = ora[ora["q"] < 0.05] if not ora.empty else ora
        print(f"--- {label}: {len(genes)} genes, {len(sig)} pathways at q<0.05 ---")
        if not ora.empty:
            print(ora.head(8)[["gene_set", "set_size", "overlap", "fold_enrichment", "q"]]
                  .to_string(index=False), "\n")

    # 6. persist
    if args.save_dir:
        out = Path(args.save_dir)
        out.mkdir(parents=True, exist_ok=True)
        res.to_csv(out / "survival_screen_genes.csv")
        hr.to_csv(out / "top_hits_hazard_ratios.csv")
        for label, ora in ora_tables.items():
            ora.to_csv(out / f"pathways_{label.split()[0]}.csv", index=False)
        print(f"Wrote results -> {out}/")


if __name__ == "__main__":
    main()
