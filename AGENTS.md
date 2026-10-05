# AGENTS.md — design rules for OmicOne and scStudio

This file is the contract for anyone changing these two repositories, human or
model. Models are switched mid-project here, so every rule below exists to keep
one style across hands. Read it before the first edit; when a rule and the code
disagree, the rule wins and the code is a bug, unless the rule is outdated, in
which case update this file in the same change.

The two repos share most of their code:

| Repo | What it is | Source of truth for |
|---|---|---|
| **OmicOne** | multi-omics shell: landing page, per-omics pipelines (single-cell, WES; bulk / spatial / integration are placeholders) | every shared file |
| **scStudio** | the single-cell pipeline on its own | only its shell (`app_ui.R`, `app_server.R`, `mod_report.R`) |

`AGENTS.md` and `CLAUDE.md` are identical in both repos.

---

## 1. Ground rules

1. **Shared code is edited in OmicOne only**, then mirrored:
   `tools/sync_mirror.sh` (copies + renames `OmicOne`→`scStudio`), then
   `tools/check_mirror.sh` must print `Mirror OK`. Never hand-edit a shared file
   in scStudio. The files allowed to differ are listed in `EXPECTED` inside both
   scripts; keep the two lists identical.
2. **Science first, then UI.** A figure, label or sentence that overstates what
   the statistics show is a bug of the highest severity (see §7).
3. **No silent behaviour.** If the app drops cells, samples, genes or patients,
   it says how many and why, in the summary strip or the insight bar.
4. **Every change ships with**: tests (§9), a `NEWS.md` entry under the
   unreleased version (Added / Changed / Fixed), and — if it changes a rule
   here — an edit to this file.
5. **Run before you hand over**: `testthat` suite, `tools/check_mirror.sh`.
   State plainly what you could and could not run.

---

## 2. Architecture

```
R/
  run_app.R        entry point (run_app), static asset path
  app_ui.R         shell UI: top bar, left navigator, main area      [per repo]
  app_server.R     shell server: routing, navigator, progress, wiring [per repo]
  app_landing.R    omics picker (OmicOne only)
  steps.R          non-single-cell registries + all_step_registries() (OmicOne only)
  steps_sc.R       single-cell registry steps_sc() + phases_sc()      [shared]
  fct_state.R      pipeline state: deps, staleness, epochs, checkpoints [shared]
  utils_ui.R       UI building blocks (step_container, explainer, i18n, pills)
  utils_server.R   server helpers (require_pkgs, progress, log_step, %||%)
  utils_inputs.R   control helpers (parse_genes, keep_selected, categorical_cols)
  fct_plots.R      plotting layer: renderers, figure download, palette, explainer canvas
  fct_compute.R    single-cell compute on Seurat (QC → markers, annotation)
  fct_scop.R       everything through scop / UCell / CellChat / LIANA / copykat
  fct_survival.R   clinical/survival maths shared by every omics
  fct_object.R     defensive Seurat accessors, reduction / grouping choices
  fct_export.R     top-bar export menu
  fct_wes.R        WES compute wrappers, WES tables and exports (OmicOne only)
  mod_<step>.R     one module per pipeline step
inst/app/www/
  custom.css       design tokens + all styling                     [shared]
  app.js           language switch, splash, landing cards          [per repo]
  explain.js       step-animation engine                           [shared]
  explain-sc.js    single-cell scenes                              [shared]
  explain-wes.js   WES scenes (OmicOne only)
docs/explain-gallery.html   every animation on one page (open in a browser)
tools/   sync_mirror.sh, check_mirror.sh, check_explain.js
```

Layering is strict — dependencies point downwards only:

```
mod_*  (UI + reactive glue, no science)
  └─ fct_*  (pure functions: data in, data out; testable without Shiny)
       └─ external packages, always called as pkg::fun()
```

* A module never contains an algorithm longer than a few lines; it calls an
  `fct_*` function. If you need a new computation, add it to the matching
  `fct_*.R` with a roxygen block and a unit test.
* `fct_*` functions never touch `input`, `output`, `session` or `rv`.
* External packages are always namespaced (`Seurat::RunPCA`). Heavy packages live
  in `Suggests`; a step checks them with `require_pkgs()` before computing.
  Never `library()` inside the package (the one sanctioned exception is
  `attachNamespace("NMF")` in `wes_signatures()`, documented there).

### 2.1 The step registry

Each pipeline is an ordered list of steps (`steps_sc()`, `steps_wes()`, ...):

```r
list(v = "normalize", n = 4, phase = "sc_prep", en = "Normalize", zh = "归一化",
     ui = mod_normalize_ui, deps = c("qc", "doublet"))
```

* `v` is the step key = module id = navigator id = status key. Globally unique.
* `deps` lists the steps whose **output this step consumes**. It drives
  staleness. Only earlier steps may be listed. Output-only steps (report,
  export) list nothing.
* Adding a step = registry entry + module file + server wiring line in
  `app_server.R` (both repos if single-cell) + the module list in
  `tests/testthat/test-outputs.R`.

### 2.2 Shared state (`rv`) — who may write what

| Field | Meaning | Written by |
|---|---|---|
| `rv$omics` | active pipeline, NULL = landing page | shell |
| `rv$obj` | single-cell Seurat object | single-cell modules |
| `rv$maf` | maftools MAF (as imported, or filtered) | `mod_wes_import`; `mod_wes_filter` (always from its pinned input, `step_input()`) |
| `rv$wes_sequenced` | every sample of the imported MAF, incl. silent-only ones: a sample the filters emptied is still sequenced (WT, TMB 0) | `mod_wes_import` only |
| `rv$clinical` | the full normalised clinical cohort, shared by all omics; never an analysis subset | `mod_clinical` only |
| `rv$wes_clin_raw` | clinical table as read with the MAF (id column renamed to `Tumor_Sample_Barcode`) | `mod_wes_import` only |
| `rv$markers` | the current marker table | `mod_markers` only |
| `rv$cellcomm` | the current cell-communication result | `mod_cellcomm` only |
| `rv$status[[step]]` | `NULL` / `TRUE` / `"stale"` | `mark_done()`, `start_epoch()` only |
| `rv$epoch_<omics>` | data-version counter | `start_epoch()` only |
| `rv$ckpt` | pinned inputs of filtering steps (plain env) | `step_input()`, `ckpt_drop()` |

`log_rv` holds the reproducibility log; write it only through `log_step()`.

Facts about the single-cell object that later steps must agree on live in
`obj@misc` (keys prefixed with the package name, rewritten by the sync script):

| Key | Meaning | Written by |
|---|---|---|
| `omicone_reduction` | reduction to build graphs/embeddings on (pca, harmony, integrated.dr) | reduce, integrate |
| `omicone_cluster_col` / `omicone_cluster_cols` | the active resolution's column / all resolutions run | cluster |
| `omicone_trajectory` | pseudotime / lineage columns of the last trajectory run | trajectory |
| `omicone_import_actions` | what `standardize_sc_object()` changed on import | import |

Read them through `default_graph_reduction()`, `graph_reduction_choices()` and
`obj_misc()`; never guess "harmony first" or hard-code `seurat_clusters`.

### 2.3 Pipeline-state rules (fct_state.R)

* **Import** calls `start_epoch(rv, "<omics>", log_rv)` *before* writing the new
  object. That clears the omics' statuses, checkpoints and log, and resets every
  module store.
* **Module result stores** are created with
  `step_results(rv, "<omics>", field = init, ...)` (or `step_result()` for a
  single value), never with bare `reactiveValues()` — so a new dataset can never
  show the previous dataset's results. Exception: the import modules.
* **Finishing a step** calls `mark_done(rv, "<step>")` after `rv$obj` is
  written. It marks finished dependents `"stale"` (amber dot + banner) and drops
  their checkpoints.
* **Filtering steps** (anything that removes cells/samples) read their input as
  `step_input(rv, "<step>", rv$obj)`, so a re-run with new thresholds starts
  from the same cells instead of compounding.
* A step must be safe to run twice in a row with the same parameters and give
  the same result.

---

## 3. Module contract

Every step module has exactly this shape (copy an existing module, e.g.
`mod_normalize.R`, rather than starting from scratch):

```r
#' Module N: <Title>            # roxygen: what the step does, in one paragraph
#' @name mod_<step>
NULL

mod_<step>_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(title = list(en=, zh=),
    what = list(en=, zh=), why = list(en=, zh=), how = list(en=, zh=),
    read = list(en=, zh=), example = list(en=, zh=))
  controls <- shiny::tagList(
    label_with_help(...), <input>, ...,
    run_button(ns("run"), "<Verb phrase>", "<动词短语>"))
  step_container(
    title = list(en=, zh=), subtitle = list(en=, zh=),
    explainer = explainer, controls = controls,
    summary = shiny::uiOutput(ns("summary")),
    preview = preview_plot_ui(ns("preview"), download = TRUE,
                              guide = list(en=, zh=), caption = list(en=, zh=)))
}

mod_<step>_server <- function(id, rv, log_rv) {
  shiny::moduleServer(id, function(input, output, session) {
    res <- step_results(rv, "sc", done = FALSE, ...)
    shiny::observeEvent(input$run, {
      shiny::req(rv$obj)
      if (!require_pkgs(c(...), "<Step>")) return(NULL)
      obj <- with_progress_notify(<fct_ call>, message = "...")
      if (is.null(obj)) return(NULL)
      rv$obj <- obj
      mark_done(rv, "<step>")
      log_step(log_rv, "<Step name>", params = list(...), code = "<exact R code>")
    })
    output$summary <- shiny::renderUI(...)      # stat_tile() pills
    render_step_plot(output, input, "preview", function() {...}, name = "<step>")
  })
}
```

Rules:

* **Inputs** are read once at the top of the observer into locals; never read
  `input$x` inside a `with_progress_notify()` block after a long computation.
* **Run button**: one per step, `run_button()`, verb phrase label
  ("Normalize", "Compute TMB"), never "Submit"/"OK".
* **Summary strip**: a `shiny::tagList()` of `stat_tile()` pills (2–5 pills),
  or a single `omicone-placeholder` div saying what to do before the first run.
  Do not wrap pills in `bslib::layout_columns()`.
* **Insight bar** (`insight_bar()`): one or two sentences computed from the
  run's own numbers, placed first in `preview`. It states the result *and* its
  main caveat (small n, many tests, unadjusted p). It never interprets beyond
  the statistics.
* **Selectors** offering metadata columns use `sc_group_cols()` (groupings,
  batches; never QC metrics or scores) or `categorical_cols()`. A selector
  rebuilt by `renderUI()` keeps the user's choice with
  `keep_selected(shiny::isolate(input$x), choices, default)`. Typed gene lists
  go through `parse_genes()`. Numeric inputs are coerced (`int_input()`) before
  use; a half-typed number must never crash the session.
* **Explainer animation**: every step has a scene with the step key in
  `explain-sc.js` / `explain-wes.js` (rules in §5.1). `preview_plot_ui(guide =)`
  shows it automatically; a step whose output is not a plot puts
  `explain_scene("<step>", en, zh)` in its empty preview.
* **Every figure is downloadable**: `preview_plot_ui(download = TRUE)` paired
  with `render_step_plot()` (ggplot / ComplexHeatmap objects) or
  `render_base_plot()` + `register_figure_download()` (base graphics,
  maftools). The file must be the same closure as the screen.
* **Errors** never crash the session: compute goes through
  `with_progress_notify()`, drawing through the renderers, which turn errors
  into a readable message on the canvas.
* **Tables** are capped (≤ 5,000 rows) for the browser, with the full table as
  a CSV download: `render_tbl_wrap()` for simple cases, `wes_table()` in the
  WES modules.

---

## 4. Text, language and tone

* **Every user-visible string is bilingual** via `i18n(en, zh)` or a
  `list(en =, zh =)` argument. No bare English strings in the UI, including
  `showNotification()` messages and pill labels.
* Only one language is ever on screen (client-side swap in `app.js`). The
  choice persists in `localStorage`; first visit follows the browser language;
  the server sees it as `input$lang`.
* The explainer has five facets, each 1–3 sentences: **what** (the operation),
  **why** (what goes wrong without it), **how** (which controls matter and their
  safe defaults), **read** (how to read *this step's* figure: axes, colours,
  what pattern is good/bad), **example** (one concrete number-bearing case).
  Do not repeat a sentence across facets.
* `guide` (empty canvas) says what will be drawn; `caption` (under the figure)
  says what one mark means. Both must describe the figure the code actually
  draws — if you change the plot, change all three texts.
* Chinese text uses full-width punctuation（，。：；）and keeps technical terms
  in English where that is the field's usage (UMAP, Leiden, TMB, MAF).
* Tone: plain, specific, no exclamation marks, no marketing words ("powerful",
  "seamless", "一键"…). Numbers carry units.

---

## 5. Visual design

* Design tokens live at the top of `custom.css` (`--sc-*`). Use tokens, never
  literal colours, in CSS. Light and dark values are both defined; the app
  follows the OS preference (`input_dark_mode()` without `mode`).
* Class prefix is the package name (`omicone-` / `scstudio-`, rewritten by the
  sync script). No inline `style=` attributes except one-off layout tweaks.
* Layout per step: header (title, summary pills, subtitle) → collapsed explainer
  strip → control rail (left, 280 px) + plot area (fills the rest).
* Plots: white background (what gets downloaded), `omicone_theme()` for ggplot,
  `sc_palette()` for categorical colours, colour-blind-safe pairs for binary
  contrasts (`#3b6ea5` kept / `#c1476b` flagged). Axis titles carry units and
  scale ("UMI count (log10)").
* Fonts: system stack from `app_font_stack()`. No web-font downloads — the app
  must render offline and where Google is blocked.
* Nothing floats over the workspace. Global actions live in the top bar.
* Motion is quiet: a step fades in when opened, a new insight rises in; every
  animation is switched off under `prefers-reduced-motion`.

### 5.1 Step animations (explain*.js)

Each step's empty preview plays a short looping scene of what the step does —
the entry point for beginners and the seed of a future interactive tutorial.

* **No libraries, nothing fetched.** Plain Canvas 2D; the app must run offline.
  Do not add three.js / GSAP / CDN scripts to the app.
* Register with `OmicOneExplain.register("<step key>", { period, stages, still,
  init, draw })`. Draw in the fixed **560 × 260** logical space; the engine
  scales, handles HiDPI, pauses off-screen, and draws `still` under reduced
  motion.
* **Three stages**, each with a short bilingual caption (≤ ~55 characters):
  the input, the operation, the reading. The scene shows the *idea*, with
  invented but plausible data; it never shows real numbers or p-values.
* **Science rules apply to pictures too**: if the analysis has a caveat users
  routinely miss (doublets per sample, one row per patient, clonal peak at
  purity / 2), the scene shows it as a footnote (`g.note()`).
* Randomness comes only from the seeded `R()` passed to `init()`, so a scene is
  identical on every load. Colours come from the theme (`g.th`) or `H.PAL`;
  mixed colours may be mixed again (`H.mix` accepts `rgb()`).
* Check every change: `node tools/check_explain.js` (no exceptions, no
  non-finite coordinates, every step has a scene) and look at
  `docs/explain-gallery.html` in light and dark, EN and 中.

---

## 6. Reproducibility log

* `log_step(log_rv, step, params, code)` after every successful run. `code`
  must be R that reproduces exactly what ran, with the parameters actually used
  (no placeholders such as `keep_cells` that the script never defines).
* Re-running a step replaces its entry; a new import clears the omics' log.
* The exported script must run top to bottom on the user's input file after
  they edit one path variable.

---

## 7. Statistical and analysis standards

These are requirements, not suggestions; the user is a translational scientist
and results may go into publications.

1. **Define the analysis set first, derive cut-offs second.** Medians,
   optimal cut-points, quantiles and z-scores are computed only after the
   subset, de-duplication and missing-value removal are final.
2. **One row per patient for outcome analysis.** Single-cell features are
   aggregated per patient before any survival model; cells are never the unit
   of a survival test. Report n patients and n events.
3. **Multiple testing**: any screen over genes, pathways, clusters or
   signatures reports BH-adjusted q-values, and sorting/highlighting uses q.
4. **Optimal cut-points are exploratory.** If `surv_cutpoint`-style
   dichotomisation is offered, label it as such and show the median split
   alongside.
5. **Small samples**: warn when a group has < 10 patients or < 10 events, or
   when a Cox model has < 10 events per covariate.
6. **The p-value on a figure is the test named next to it** (log-rank, Wald,
   Fisher...). Never mix a Cox p-value onto a KM curve without saying so.
7. **Pseudo-replication**: cell-level tests (e.g. FindAllMarkers) are labelled
   as marker discovery, not as between-condition inference. Condition
   comparisons across samples use pseudo-bulk or sample-level statistics:
   the Pseudobulk DE and Differential abundance steps (`fct_sample.R`). A
   condition must be constant within a sample (`sample_table()` refuses a
   cell-level column), the FDR of pseudobulk DE is computed within each cell
   type, and a cell type with < 2 samples per group is skipped, not tested.
8. **Multi-sample data**: per-sample operations stay per-sample (doublet
   detection, QC thresholds where batches differ).
9. **State the denominator.** Percentages and burdens say what they are
   divided by (cells kept of cells in; mutations per captured Mb, with the
   capture size used).
10. **Continuous first.** For a continuous biomarker the primary result is the
    Cox HR per SD (`cox_continuous()`); the high/low Kaplan-Meier split is an
    illustration. An optimal cut-point is reported with its selection-adjusted
    p-value (`split_numeric(..., "optimal")` → `attr(, "p_adjusted")`).
11. **Enrichment backgrounds.** ORA uses the genes actually tested as the
    universe; GSEA ranks all tested genes (never a pre-filtered list).
12. **Who is wild-type.** In WES, a sequenced sample with no mutation record is
    wild-type, not missing (`wes_mutation_status(universe =)`), and its TMB is
    0. That includes samples the Variant filters step left with no variant:
    pass `rv$wes_sequenced` (`sequenced =` / `samples =`).

---

## 8. Naming and code style

* Files: `mod_<step>.R`, `fct_<domain>.R`, `utils_<ui|server>.R`.
* Functions: `snake_case`; module functions `mod_<step>_ui/_server`; compute
  wrappers `<domain>_<verb>` (`wes_tmb`, `qc_mad_keep`).
* One assignment per line; no `;` chaining.
* Comments explain *why* (the trap avoided, the paper followed), not *what*.
* roxygen block on every function, `@keywords internal` unless exported.
* 2-space indent, `<-` for assignment, double quotes, lines ≤ 100 chars.

---

## 9. Tests

* `tests/testthat/test-logic.R` — pure functions in `fct_*`; must run without
  Seurat, maftools or network.
* `test-state.R` — pipeline-state rules (deps, staleness, epochs, checkpoints).
* `test-outputs.R` — every module's outputs evaluate with and without data;
  add every new module to its list.
* `test-survival.R`, `test-wes.R` — statistical helpers on known inputs;
  WES tests `skip_if_not_installed("maftools")`.
* `test-sc-compute.R` — the Seurat compute chain on the bundled data, including
  replaying the logged code to identical results (`skip_if_not_installed`).
* `test-scop.R` — every `pkg::fun(arg =)` in fct_scop.R exists in the installed
  package with that argument (skips when the package is absent).
* `node tools/check_explain.js` — the step animations (see §5.1).
* A bug fix adds the test that would have caught it.

---

## 10. Hand-over checklist

- [ ] Rules in §3–§7 hold for every step you touched
- [ ] `testthat` passes; `tools/check_mirror.sh` prints `Mirror OK`;
      `node tools/check_explain.js` exits 0
- [ ] `NEWS.md` updated; this file updated if a rule changed
- [ ] You said what you could not run (missing packages, no live data)
