# OmicOne

> 🌐 **Languages:** **English** · [中文版 README](README.zh-CN.md)

**Interactive multi-omics analysis on your own machine.** Launch a browser
interface with one command, pick the kind of data you have, and work through
that pipeline step by step — using your own computer's CPU and RAM. No cloud
server, nothing uploaded anywhere.

Every step is designed for beginners *and* experts:

- 💡 a plain-language **"What is this step?"** explainer (with a worked example)
- 🔧 a **method choice** (each step offers alternatives — e.g. UMAP *or* t-SNE)
- 🎚️ **adjustable thresholds** with sensible defaults
- 📊 a **result summary** and a large **preview plot**
- 🧾 a **reproducibility log** you can export as an R script or a narrated report
- 🎞️ a short **animation** in every step's empty canvas showing what the step does
  (preview them all in `docs/explain-gallery.html`)
- 🌗 light / dark themes (following your OS), English / 中文 — switchable from the top bar

> **Renamed:** this project was called **OMICstudio** up to v0.5.0. The name
> clashed with [OmicStudio](https://www.omicstudio.cn), an established cloud
> platform — an unhelpful collision for a tool whose whole point is that nothing
> leaves your machine. Same code, new name; see [NEWS.md](NEWS.md).

> **Status (v0.6.0):** the **single-cell RNA-seq** (21 steps) and **WES /
> somatic mutation** (12 steps) pipelines are complete and statically validated.
> The other three omics show their planned steps in the interface and are not
> implemented yet. Neither pipeline has been run end-to-end against a live
> install of its engine; see [Caveats](#caveats).

---

## The five pipelines

Choosing a card on the start screen routes you into that pipeline.

| Omics | Engine | Status |
|---|---|---|
| **Single-cell RNA-seq** | [scop](https://github.com/mengxu98/scop) + Seurat | ✅ complete (21 steps) |
| **WES / somatic mutations** | [maftools](https://bioconductor.org/packages/maftools) | ✅ complete (12 steps) |
| **Bulk RNA-seq** | TOmicsVis + DESeq2 | 🚧 roadmap shown in-app |
| **Spatial transcriptomics** | Seurat + SpatialExperiment | 🚧 roadmap shown in-app |
| **Multi-omics integration** | MOFA2 / iClusterPlus / SNFtool | 🚧 roadmap shown in-app |

Clinical follow-up is **shared across pipelines**: a cohort loaded in the
single-cell *Clinical & survival* step is offered to the WES *Mutation vs
survival* step as a source (the full cohort, never an analysis subset) — same
curves, same log-rank test, same code path.

maftools installs as an ordinary Bioconductor binary, so the WES pipeline needs
no source build and runs on a machine that blocks compilation.

---

## Why "localhost-first"?

Real omics analysis (Seurat/Bioconductor) needs native compute and real RAM.
Browser-only (WASM) apps can't run it. So OmicOne runs a **local server +
browser UI**: the interface is a web page, but all computation happens on *your*
machine. This is the same model as `cellxgene launch`.

---

## Three ways to run it

Pick the tier that matches how much you want to install.

### A. You already have R (lightest)
```r
# install.packages("remotes")
remotes::install_github("Tianqi-Ma/OmicOne")
OmicOne::run_app()   # opens your browser automatically
```
Needs R ≥ 4.1 with **shiny ≥ 1.7.4** and **bslib ≥ 0.7.0** (shiny ≥ 1.8.1 is
recommended), plus the heavy analysis packages for the steps you actually run.
Smallest download, one command.

### B. No R, no dependencies → Docker (recommended for most users)
Everything (R + Seurat + Bioconductor + scop + the app) is baked into one image.
You only need [Docker](https://www.docker.com/products/docker-desktop/).
```bash
docker build -t omicone .           # build once from this repo
docker run --rm -p 3838:3838 -m 16g omicone
# then open http://localhost:3838 in your browser
```
Upload your data through the browser (no volume mounting needed). Give Docker
enough memory (`-m 16g`, and raise Docker Desktop's memory limit for large data).

### C. Non-technical, double-click (planned)
A desktop installer (Tauri/Electron/electricShine) that bundles R and the app so
users just double-click — no terminal, no Docker. Planned for a later release.

---

## Try it instantly (no data needed)

On the **Import** step, keep the source on **Demo data** and click *Load demo
data*. Three demos are offered:

| Demo | What it is | Needs |
|---|---|---|
| **PBMC 3k** (default) | the classic 2,700-cell 10x PBMC dataset, **bundled in the package** | nothing — instant, offline |
| **Pancreas** | mouse pancreas with lineage structure and spliced/unspliced layers, for trajectory / velocity | `scop` installed |
| **Tiny example** | 500 × 300 synthetic matrix for a quick UI check | nothing |

You can also paste a URL to any `.rds` / `.h5` file.

## What you can upload

- **RDS** — a saved `Seurat` or `SingleCellExperiment` object (old Seurat objects
  are updated automatically)
- **10x `.h5`** — Cell Ranger HDF5
- **Counts table** — CSV/TSV, genes in rows, cells in columns
- **Clinical table** — CSV/TSV, one row per patient (id, follow-up time, outcome)

Browser uploads are capped high (5 GB by default; change with
`run_app(max_upload_mb = ...)`).

---

## The single-cell steps

Eight phases, twenty-three steps. Bold = the default method.

| Phase | # | Step | Methods offered |
|---|---|------|-----------------|
| **Data & QC** | 1 | Import | RDS / 10x .h5 / table / URL / demo |
| | 2 | Quality control | **MAD (adaptive)** / manual; mito, ribo, hemoglobin, dissociation % |
| | 3 | Doublets | **scDblFinder** / DoubletFinder |
| **Preprocess** | 4 | Normalize | **LogNormalize** / SCTransform |
| | 5 | Features / PCA | HVG **vst** / mvp / dispersion |
| | 6 | Integrate | **none** / Harmony / CCA / RPCA (Seurat v5), with a batch-mixing score |
| **Structure** | 7 | Cluster | **Leiden** / Louvain |
| | 8 | Embed | **UMAP** / t-SNE |
| **Identity** | 9 | Markers | **wilcox** / roc / MAST |
| | 10 | Annotate | **manual** / SingleR (per cluster, pruned labels) / Azimuth |
| | 11 | Enrichment / GSEA | ORA + GSEA (clusterProfiler via scop) |
| **Compare conditions** | 12 | Pseudobulk DE | counts summed per sample × cell type; **edgeR quasi-likelihood** / limma-voom / DESeq2; optional pairing / batch covariate |
| | 13 | Differential abundance | cell-type proportions per sample, propeller (logit + limma robust eBayes, BH) |
| **Trajectory** | 14 | Trajectory | **Slingshot** / Monocle2 / Monocle3 / PAGA / Palantir |
| | 15 | RNA velocity | scVelo (steady-state / stochastic / dynamical) |
| | 16 | Dynamic features | scop dynamic features + heatmap |
| **Advanced** | 17 | Cell cycle & signatures | Seurat cell cycle, UCell / AddModuleScore |
| | 18 | Cell communication | LIANA / CellChat |
| | 19 | Malignant / CNV | CopyKAT, stemness |
| | 20 | **Clinical & survival** | Kaplan-Meier + log-rank + univariable Cox; stratify by a clinical column or by per-sample cell-type composition (median / tertile / optimal cutpoint) |
| **Output** | 21 | Visualize | UMAP / violin / dotplot / feature / heatmap |
| | 22 | Report | narrated HTML report from the reproducibility log |
| | 23 | Export | .rds / .h5ad (best-effort) / figures / R script |

A **⤓ Export** menu in the top bar is available on *every* step (object,
metadata, counts, embeddings) so you never have to reach the last step to get
your data out.

Methods and defaults follow current (2023–2025) best practice: MAD-based QC,
scDblFinder, Leiden, Harmony, scop/SCP plotting throughout.

---

## The WES steps

Five phases, fifteen steps, all over maftools.

| Phase | # | Step | What it does |
|---|---|------|-----------------|
| **Input** | 1 | Import MAF | MAF (+ optional clinical table), or the bundled TCGA LAML demo |
| | 2 | Variant filters | FILTER = PASS, tumour depth, alt reads, VAF, population AF (gnomAD / ExAC), normal alt reads — a funnel of what each removed; hypermutator flag (optional step) |
| | 3 | Cohort summary | variant classes, per-sample burden, top mutated genes |
| **Landscape** | 4 | Oncoplot | the waterfall view, with clinical annotation bars |
| | 5 | TiTv / VAF / rainfall | mutation spectrum, allele frequencies, kataegis |
| | 6 | TMB | mutations per Mb, per sample, over your capture size |
| | 7 | Lollipop / domains | one gene's mutations along the protein, over its domains |
| | 8 | Drivers & interactions | oncodrive positional clustering; co-occurrence / mutual exclusivity |
| **Signatures** | 9 | Mutational signatures | trinucleotide matrix → de-novo signatures → COSMIC match |
| **Clinical & prognosis** | 10 | Clinical / pathway / drug | gene enrichment by clinical group, oncogenic pathways, drug-gene interactions |
| | 11 | Cohort comparison | Fisher test between two clinical groups, forest plot |
| | 12 | Mutation vs survival | mutant vs WT Kaplan-Meier + log-rank, via the **shared** survival layer |
| | 13 | TMB vs outcome | Cox HR or logistic OR per doubling of TMB (no cut-off), ROC AUC for a response; median-split KM as illustration |
| | 14 | Heterogeneity | VAF clustering per sample, MATH score |
| **Output** | 15 | Report | HTML report of the WES steps + the runnable R script |

Step 1 offers maftools' bundled **TCGA LAML** cohort (193 samples), so the whole
pipeline is explorable offline with no data of your own.

---

## Caveats

- **What has and has not been run (2026-10 audit).** The core single-cell chain
  (import → QC → doublets → normalize → PCA → integrate → cluster → UMAP →
  markers) and the clinical/survival maths were run on the bundled data with
  Seurat 4.3, and their logged R code replays to identical results. The WES
  pipeline was run on maftools' TCGA-LAML example (maftools 2.14). The scop-backed
  steps (enrichment, trajectory, velocity, dynamic features), CellChat / LIANA,
  copykat, SingleR, Azimuth and scDblFinder were checked against their current
  source code but **not executed**; run them once on your own install before
  relying on them.
- **`scop` is the plotting and compute engine** for the single-cell pipeline and
  is a GitHub package; if your machine blocks source builds, use the Docker
  image, which bakes in scop and a pre-built conda environment for the
  Python-backed steps (scVelo, PAGA, Palantir).
- **Heavy packages are optional at install time** (they live in `Suggests`). The
  UI loads without them; each compute step checks for what it needs and tells you
  what to install if it's missing.
- **Some steps need internet** the first time (SingleR/Azimuth reference download).
- **`.h5ad` export** in pure R relies on SeuratDisk (best-effort); `.rds` is the
  primary export format.
- **The "optimal" survival cutpoint is exploratory.** It is chosen to maximise
  separation, so its p-value is optimistic — report the median split, or validate
  the cutpoint in an independent cohort.
- **Mutational signatures need a BSgenome package** (`BSgenome.Hsapiens.UCSC.hg19`
  or hg38) plus `NMF`. That is a large download; every other WES step works
  without it.
- **VAF-based steps need a VAF column.** Many callers do not report one, in which
  case the VAF, rainfall and heterogeneity views say so rather than guessing.
- **Memory** scales with cell count. <100k cells is comfortable on 16–32 GB;
  larger needs more. The app warns you and offers downsampling.

---

## Development

R package layout (golem-style):

```
R/app_ui.R            top-level shell (top bar, sidebar, splash)
R/app_landing.R       omics chooser
R/app_server.R        omics routing + shared hub (rv$obj, rv$clinical, rv$status)
R/steps.R             per-omics step registries — the single source of truth
R/mod_*.R             one module per step
R/fct_compute.R       scop/Seurat compute wrappers
R/fct_wes.R           maftools wrappers + the MAF -> survival bridge
R/fct_plots.R         scop plotting wrappers + shared theme/palette
R/fct_survival.R      omics-agnostic survival maths (survival + ggplot2)
R/utils_ui.R          the shared step layout, bilingual helpers
```

Run locally during development with:
```r
pkgload::load_all(); run_app()
```

Before committing, the checks that must pass:
```r
testthat::test_dir("tests/testthat")   # includes a guard that every module's
                                       # outputs evaluate, and that no function
                                       # calls a name that does not exist
```
```sh
tools/check_mirror.sh                  # shared files must stay identical to
                                       # scStudio apart from the package name
```

The single-cell pipeline is also maintained standalone as
[scStudio](https://github.com/Tianqi-Ma/scStudio); shared changes are mirrored
between the two repositories.

## License

MIT.
