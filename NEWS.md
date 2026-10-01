# OmicOne 0.6.1 (unreleased)

## Added
- **Global progress chip** in the topbar (done/total with a mini bar) and a
  floating **"next step" chip** that appears bottom-right once the current
  step completes — click it to jump to the next unfinished step.
- **Auto-generated insight bars** above every WES result: one bilingual
  sentence computed from the actual numbers (e.g. median TMB and where it
  sits against the immunotherapy-relevant line; how many genes pass the FDR
  cutoff; MATH score with its clonality reading; log-rank p with a
  power caveat when groups are small).
- **Figure captions**: the primary plot of every step (WES and single-cell)
  carries a one-line caption to the left of the export row.
- **Step subtitles and reading guides**: every WES step now shows an
  always-visible one-line subtitle under its title, the explainer gains a
  "How to read the result / 如何解读结果" facet (axes, colours, which
  patterns matter, rough cutoffs), and the empty preview canvas shows a
  centred guide — what the step will draw and how to run it — that steps
  aside automatically once the plot renders (pure CSS, no JS). All 20
  single-cell steps gained the same three layers.
- **Figure export on every WES plot**: each preview now carries a download row
  with a format selector (PNG / JPEG / vector PDF) and, for raster formats, a
  DPI field (72–1200). The export replays the exact draw closure used on
  screen, so the file matches the preview — for maftools' base-graphics plots
  and the ggplot Kaplan-Meier curve alike.
- **TMB vs TCGA**: the Tumour mutational burden step gains a tab that plots the
  cohort's median TMB against all 33 TCGA cohorts
  (`maftools::tcgaCompare()`), labelled with the imported cohort's name.

## Changed
- **Adaptive plot fonts**: text now scales with how much is drawn — oncoplot
  gene labels with the number of genes, the cohort summary with sample count,
  interaction and oncodrive plots with genes tested. The lollipop plot repels
  labels (berryFunctions) and rotates them vertical when the labelled
  positions are densely packed, so hotspots like TP53 stay readable.
- **Plot text is larger overall**: the cohort-summary dashboard and every
  maftools plot without its own font argument get a base-graphics text boost
  (`with_text_boost()`), so axis and title text stays legible on wide
  browser panels.
- **Workspace layout polish**: the "What is this step?" explainer is now a
  full-width strip above the workspace instead of being squeezed into the
  280 px control rail; the left navigator's phases collapse into expandable
  groups with a done-count badge (only the current phase stays open); the
  nested tab card inside the plot area is flattened (no card-in-card
  borders); summary pills and the export row are visually quieter (the
  export row fades in on hover).

## Fixed
- **Export row alignment**: the DPI label is now an inline unit next to its
  field, so the download button, format selector and DPI input sit on one
  baseline.
- **Mutational signatures run end-to-end** (BSgenome.Hsapiens.UCSC.hg19 + NMF).
  Two environment pitfalls are worked around in `wes_signatures()`: NMF 0.28's
  internal `getGeneric("seed")` returns NULL once the Bioconductor stack is
  attached (NMF is now attached explicitly), and NMF's parallel mode expects a
  foreach backend Shiny never registers (extraction now runs sequentially).
- **TMB**: a blank or non-numeric capture size shows a clear message instead
  of an error.
- The APOBEC tab now explains when a cohort simply has no APOBEC enrichment
  instead of showing maftools' raw error.

## Validation
- Full WES pipeline re-validated against maftools 2.26.0 (Bioc 3.22, R 4.5):
  24/24 wrapper checks, 12/12 modules via `shiny::testServer` (now including
  signatures), 15/15 figure downloads byte-verified across PNG/JPEG/PDF.

---

# OmicOne 0.6.0

Renamed from **OMICstudio** to **OmicOne**.

The old name collided with [OmicStudio](https://www.omicstudio.cn), an
established cloud analysis platform — confusing in general, and especially so
here, since this project is the opposite thing: everything runs on your own
machine and nothing is uploaded. No functionality changed in this release; it is
the rename and nothing else.

What this means for you:

- The package is now `OmicOne`: `OmicOne::run_app()`,
  `remotes::install_github("Tianqi-Ma/OmicOne")`.
- The repository moved to `Tianqi-Ma/OmicOne`. GitHub redirects the old URL, and
  an existing clone keeps working, but update your remote:
  `git remote set-url origin git@github.com:Tianqi-Ma/OmicOne.git`
- If you had the old package installed, remove it:
  `remove.packages("OMICstudio")`.
- Exported filenames now start with `omicone_` instead of `omicstudio_`.

Release notes below 0.6.0 were published under the old name.

---

# OmicOne 0.5.0

The WES / somatic mutation pipeline. Two of the five pipelines are now complete.

## Added
- **Twelve WES steps over maftools**, wired into the omics router the same way
  single-cell is: Import MAF, Cohort summary, Oncoplot, TiTv/VAF/rainfall, TMB,
  Lollipop/domains, Drivers & interactions, Mutational signatures,
  Clinical/pathway/drug, Cohort comparison, Mutation vs survival, Heterogeneity.
  `R/fct_wes.R` holds the wrappers; each module keeps the same shape as the
  single-cell ones (explainer, method choice, thresholds, run, summary, preview).
- **Instant offline WES demo**: maftools ships the TCGA LAML cohort (193
  samples), so the whole pipeline is explorable with no data of your own —
  the counterpart to the bundled pbmc3k.
- **The survival layer pays off across omics.** *Mutation vs survival* does not
  reimplement anything: it turns a MAF into per-sample mutation status and hands
  that to `fct_survival.R`. A cohort loaded in the single-cell *Clinical &
  survival* step is picked up automatically, so the curves, the log-rank test and
  the Cox screen are identical on both sides.
- **Per-omics object slots.** `rv$obj` (Seurat) and `rv$maf` (MAF) are separate,
  and the sidebar reports whichever pipeline is active. Previously, switching
  from single-cell to another omics left the sidebar showing the Seurat object's
  cell and gene counts.
- **Landing-page availability badges**, so which pipelines actually run is
  visible before clicking into one.
- `render_base_plot()` for maftools' base-graphics output, with the same
  readable-error-on-canvas behaviour as the ggplot renderer.
- `wes_missing_api()`: checks an installed maftools against every entry point
  the modules call, so version drift is one command to diagnose.

## Notes
- maftools is a plain Bioconductor binary — no source build — so this pipeline
  runs on machines where the single-cell engine (scop, a GitHub package) cannot
  be installed.
- Mutational signatures additionally need a BSgenome package and `NMF`; every
  other WES step works without them.
- `maftools::pathways()` was `OncogenicPathways()` before 2.12; the wrapper
  handles both.

## Known limitations
- The WES pipeline is statically validated (parse, UI build, output evaluation,
  tests, install) but has **not been run against a live maftools install with
  real data**; some argument names may still need adjusting.
- Bulk, spatial and integration remain roadmap only.

---

# OmicOne 0.4.0

Multi-omics suite. OmicOne grew out of
[scStudio](https://github.com/Tianqi-Ma/scStudio) (single-cell only), which is
still maintained standalone; shared changes are mirrored between the two.

## Added
- **Omics landing page.** Five cards — single-cell RNA-seq, bulk RNA-seq, WES,
  spatial transcriptomics, multi-omics integration — each with a themed canvas
  animation on hover. Picking one routes into that pipeline; a "← Omics" link in
  the top bar returns to the chooser.
- **Per-omics step registries** (`R/steps.R`) as the single source of truth for
  the navigator and the workspace. Bulk / WES / spatial / integration show their
  planned steps through a placeholder module, each card naming its own step.
- **Clinical & survival step** (single-cell step 18), backed by a new
  omics-agnostic `R/fct_survival.R`: Kaplan-Meier curves drawn in ggplot2, the
  log-rank test, median survival per arm, and univariable Cox regression.
  Stratify by a clinical column or by **per-sample cell-type composition**
  (median / tertile / log-rank-optimal cutpoint). The maths uses `survival`,
  which ships with R, so nothing needs to be compiled. The normalised cohort is
  written to the shared hub (`rv$clinical`) for the other pipelines to reuse.
- **Import overview.** The Import step now previews the raw data before any
  filtering: per-cell QC violins, top-20 expressed genes, a counts-vs-genes
  scatter, the cell metadata table and a slice of the counts matrix. Species is
  detected from gene-name casing so mouse data no longer reads 0% mitochondrial.
- **Shared export menu** in the top bar, available on every step: the object
  (.rds), cell metadata (.csv), counts matrix (.rds) and embeddings (.csv).
- **Bundled pbmc3k.** The real 2,700-cell 10x PBMC dataset ships inside the
  package, so the recommended demo loads instantly and offline — no download and
  no SeuratData.

## Fixed
- The Import overview called `scstudio_theme()`, a name that does not exist in
  this package after the rename, so the whole overview panel rendered as an
  error message. Two tests now guard this: every module's outputs are forced to
  evaluate, and every function is checked for calls to names that do not resolve.
- The bundled launchers (`inst/launch.bat`, `inst/launch.command`) still invoked
  `scStudio::run_app()` and could never have worked.
- `DESCRIPTION` declared `bslib (>= 0.5.0)` while the UI uses
  `input_dark_mode()` (bslib 0.6.0) and `input_task_button()` (bslib 0.7.0);
  installing against an older bslib produced an app that crashed on startup.
  Now `bslib (>= 0.7.0)`.
- Interactive plots asked `ggplotly()` for a `text` tooltip that no layer
  provides, so hovering showed an empty box; and a multi-gene FeaturePlot (a
  patchwork) was handed to `ggplotly()`, which cannot render it. Interactivity is
  now used only where it works, and shows the real values.
- Static assets 404'd without `addResourcePath()`; plot failures surfaced as an
  empty "Plot error:" or an opaque `[object Object]`.
- Visualize / Report / Export never marked themselves complete, so the step
  navigator could not reach the end of the pipeline.
- The Export step could error on its own filename when the format input had not
  been initialised.
- Missing styles for the sidebar's empty-state and dataset readout; the Google
  font is now optional rather than a hard requirement for an offline app.

## Known limitations
- **Not validated end-to-end.** The single-cell pipeline is fully wired and
  statically validated (parse, UI build, output evaluation, tests, install) but
  has not been run against a live Seurat/scop install with real data; some scop
  argument names still need checking.
- Bulk, WES, spatial and integration are roadmap only.
- The "optimal" survival cutpoint is exploratory: chosen to maximise separation,
  so its p-value is optimistic.
- `.h5ad` export is best-effort (SeuratDisk); `.rds` is primary.

---

# scStudio 0.3.0 (inherited history)

- **scop as the engine.** Plotting (`CellDimPlot`, `GroupHeatmap`,
  `DynamicHeatmap`, `VolcanoPlot`, ...) and compute (`Standard_SCP`,
  `Integration_SCP`, `RunDEtest`, `RunEnrichment`, `RunGSEA`, `RunSlingshot`,
  Monocle2/3, PAGA, Palantir, WOT, scVelo, dynamic features) behind wrappers that
  degrade to a readable message when scop is absent.
- Eight new steps: enrichment/GSEA, trajectory, RNA velocity, dynamic features,
  cell cycle & signatures, cell communication (LIANA/CellChat), malignancy/CNV
  (CopyKAT), and the narrated report. Twenty steps across seven phases.
- Mascarade cluster outlining on the embedding and visualization steps.

# scStudio 0.2.0

- Dark theme by default with a light toggle; grouped left step navigator with
  per-step status; `input_task_button()` so a running step is visibly running;
  plot-first layout (narrow control rail, large canvas, collapsed explainer).
- Single-language interface (English or 中文, never both) via `data-en`/`data-zh`
  attributes swapped client-side — no server round-trip.
- Startup animation: scattered cells coalescing into UMAP-like clusters.

# scStudio 0.1.0

- First scaffold: localhost-first Shiny app launched with `run_app()`, twelve
  analysis modules following a uniform pattern, MAD-based QC, scDblFinder,
  Harmony, Leiden, UMAP/t-SNE, SingleR/Azimuth, a reproducibility log exportable
  as an R script, a Docker image, and a bilingual README.
