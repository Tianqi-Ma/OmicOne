#' Per-omics pipeline step registries
#'
#' Single source of truth for each omics type's ordered steps (value key, phase,
#' bilingual label, and the module UI function). The navigator and the hidden
#' tabset are both built from these. Single-cell is fully implemented; the other
#' omics show their planned steps via a placeholder module until implemented.
#'
#' @name steps
#' @keywords internal
NULL

#' Steps for the currently selected omics
#' @param omics One of "sc","wes","bulk","spatial","integration".
#' @keywords internal
steps_for <- function(omics) {
  switch(omics,
         sc          = steps_sc(),
         wes         = steps_wes(),
         bulk        = steps_bulk(),
         spatial     = steps_spatial(),
         integration = steps_integration(),
         list())
}

# The single-cell registry, steps_sc(), lives in steps_sc.R: it is shared
# verbatim with scStudio.

#' Every registry, keyed by omics (used by the pipeline-state helpers)
#' @keywords internal
all_step_registries <- function() {
  list(sc = steps_sc(), wes = steps_wes(), bulk = steps_bulk(),
       spatial = steps_spatial(), integration = steps_integration())
}

# Helper to build a placeholder ("coming soon") step.
.ph <- function(v, n, phase, en, zh) {
  list(v = v, n = n, phase = phase, en = en, zh = zh, ui = mod_placeholder_ui,
       deps = character(0))
}

#' WES (maftools) steps — fully implemented
#'
#' Every step reads the imported MAF and nothing else from the pipeline, so
#' each depends on `wes_import` only (the clinical table is loaded inside the
#' steps that need it).
#' @keywords internal
steps_wes <- function() {
  list(
    list(v = "wes_import", n = 1,  phase = "wes_io",   en = "Import MAF",        zh = "导入 MAF",      ui = mod_wes_import_ui,
         deps = character(0)),
    list(v = "wes_summary",n = 2,  phase = "wes_io",   en = "Cohort summary",    zh = "队列概览",      ui = mod_wes_summary_ui,
         deps = "wes_import"),
    list(v = "wes_onco",   n = 3,  phase = "wes_land", en = "Oncoplot",          zh = "Oncoplot",      ui = mod_wes_onco_ui,
         deps = "wes_import"),
    list(v = "wes_titv",   n = 4,  phase = "wes_land", en = "TiTv / VAF / rainfall", zh = "TiTv/VAF/rainfall", ui = mod_wes_titv_ui,
         deps = "wes_import"),
    list(v = "wes_tmb",    n = 5,  phase = "wes_land", en = "TMB",               zh = "突变负荷 TMB",  ui = mod_wes_tmb_ui,
         deps = "wes_import"),
    list(v = "wes_lolli",  n = 6,  phase = "wes_land", en = "Lollipop / domains",zh = "Lollipop/结构域", ui = mod_wes_lolli_ui,
         deps = "wes_import"),
    list(v = "wes_driver", n = 7,  phase = "wes_land", en = "Drivers & interactions", zh = "驱动基因与互作", ui = mod_wes_driver_ui,
         deps = "wes_import"),
    list(v = "wes_sig",    n = 8,  phase = "wes_sig",  en = "Mutational signatures", zh = "突变特征",  ui = mod_wes_sig_ui,
         deps = "wes_import"),
    list(v = "wes_clin",   n = 9,  phase = "wes_prog", en = "Clinical / pathway / drug", zh = "临床/通路/药物", ui = mod_wes_clin_ui,
         deps = "wes_import"),
    list(v = "wes_compare",n = 10, phase = "wes_prog", en = "Cohort comparison", zh = "队列比较",      ui = mod_wes_compare_ui,
         deps = "wes_import"),
    list(v = "wes_surv",   n = 11, phase = "wes_prog", en = "Mutation vs survival", zh = "突变-预后",  ui = mod_wes_surv_ui,
         deps = "wes_import"),
    list(v = "wes_hetero", n = 12, phase = "wes_prog", en = "Heterogeneity",     zh = "异质性",        ui = mod_wes_hetero_ui,
         deps = "wes_import")
  )
}

#' Bulk RNA-seq (TOmicsVis + backends) steps — planned (placeholder)
#' @keywords internal
steps_bulk <- function() {
  list(
    .ph("bulk_import", 1, "bulk_io",  "Import & groups",   "导入与分组"),
    .ph("bulk_qc",     2, "bulk_io",  "QC / overview",     "QC/概览"),
    .ph("bulk_de",     3, "bulk_de",  "Differential expression", "差异表达"),
    .ph("bulk_degviz", 4, "bulk_de",  "DEG plots",         "差异基因图"),
    .ph("bulk_sets",   5, "bulk_de",  "Set operations",    "集合运算"),
    .ph("bulk_trend",  6, "bulk_sub", "Trend clustering",  "趋势聚类"),
    .ph("bulk_subtype",7, "bulk_sub", "Molecular subtyping", "分子分型"),
    .ph("bulk_prog",   8, "bulk_sub", "Prognosis (subtype vs OS)", "预后(亚型-OS)"),
    .ph("bulk_enrich", 9, "bulk_enr", "Enrichment / GSEA", "富集/GSEA"),
    .ph("bulk_wgcna", 10, "bulk_enr", "WGCNA",             "WGCNA")
  )
}

#' Spatial (Seurat + SpatialExperiment) steps — planned (placeholder)
#' @keywords internal
steps_spatial <- function() {
  list(
    .ph("sp_import",   1, "sp_io",     "Import (Visium/Xenium/...)", "导入(Visium/Xenium)"),
    .ph("sp_qc",       2, "sp_io",     "QC",                "质控"),
    .ph("sp_norm",     3, "sp_struct", "Normalize / cluster", "归一化/聚类"),
    .ph("sp_view",     4, "sp_spatial","Spatial visualization", "空间可视化"),
    .ph("sp_svg",      5, "sp_spatial","Spatially variable genes", "空间可变基因"),
    .ph("sp_domain",   6, "sp_spatial","Spatial domains",   "空间域"),
    .ph("sp_deconv",   7, "sp_integ",  "Deconvolution",     "去卷积"),
    .ph("sp_ref",      8, "sp_integ",  "scRNA reference mapping", "scRNA 参考映射"),
    .ph("sp_ccc",      9, "sp_integ",  "Spatial cell communication", "空间细胞通讯")
  )
}

#' Multi-omics integration steps — planned (placeholder)
#' @keywords internal
steps_integration <- function() {
  list(
    .ph("int_import",  1, "int_io",    "Import & align omics", "导入与对齐"),
    .ph("int_prep",    2, "int_io",    "Preprocess / feature select", "预处理/特征选择"),
    .ph("int_cluster", 3, "int_integ", "Integrative subtyping (MOFA2/iCluster/SNF)", "整合分型"),
    .ph("int_prog",    4, "int_prog",  "Prognostic stratification", "预后分层"),
    .ph("int_assoc",   5, "int_prog",  "Cross-omics association", "跨组学关联"),
    .ph("int_report",  6, "int_prog",  "Report",            "报告")
  )
}

#' Phase labels (en/zh) covering every phase key used above
#' @keywords internal
app_phases <- function() {
  c(phases_sc(), list(
    # WES
    wes_io = list(en = "Input", zh = "输入"),
    wes_land = list(en = "Landscape", zh = "突变全景"),
    wes_sig = list(en = "Signatures", zh = "突变特征"),
    wes_prog = list(en = "Clinical & prognosis", zh = "临床与预后"),
    # bulk
    bulk_io = list(en = "Input & QC", zh = "输入与质控"),
    bulk_de = list(en = "Differential expression", zh = "差异表达"),
    bulk_sub = list(en = "Subtyping & prognosis", zh = "分型与预后"),
    bulk_enr = list(en = "Enrichment & networks", zh = "富集与网络"),
    # spatial
    sp_io = list(en = "Input & QC", zh = "输入与质控"),
    sp_struct = list(en = "Structure", zh = "结构"),
    sp_spatial = list(en = "Spatial analysis", zh = "空间分析"),
    sp_integ = list(en = "Deconvolution & integration", zh = "去卷积与整合"),
    # integration
    int_io = list(en = "Input", zh = "输入"),
    int_integ = list(en = "Integration", zh = "整合"),
    int_prog = list(en = "Prognosis & output", zh = "预后与产出")
  ))
}

#' Bilingual label of a step, looked up by its key across every registry
#'
#' Lets a shared module (e.g. the "coming soon" placeholder) name the step it is
#' standing in for, rather than showing the same anonymous card 37 times.
#'
#' @param v A step key such as "wes_onco".
#' @return list(en=, zh=); the key itself if it is not in any registry.
#' @keywords internal
step_label <- function(v) {
  for (om in c("sc", "wes", "bulk", "spatial", "integration")) {
    for (s in steps_for(om)) {
      if (identical(s$v, v)) return(list(en = s$en, zh = s$zh))
    }
  }
  list(en = v, zh = v)
}

#' Ordered phase keys for one omics (first-appearance order in its steps)
#' @keywords internal
phase_order_for <- function(omics) {
  unique(vapply(steps_for(omics), function(s) s$phase, character(1)))
}

#' The omics catalogue for the landing page
#' @keywords internal
omics_catalogue <- function() {
  list(
    list(v = "sc",          icon = "circle-nodes", en = "Single-cell RNA-seq", zh = "单细胞 RNA-seq",
         desc_en = "QC, clustering, annotation, trajectory, velocity, survival (scop engine).",
         desc_zh = "质控、聚类、注释、轨迹、速率、生存分析（scop 引擎）。", anim = "sc",
         ready = TRUE),
    list(v = "bulk",        icon = "chart-column", en = "Bulk RNA-seq", zh = "Bulk RNA-seq",
         desc_en = "DE, subtyping, enrichment, WGCNA, prognosis (TOmicsVis).",
         desc_zh = "差异表达、分子分型、富集、WGCNA、预后（TOmicsVis）。", anim = "bulk"),
    list(v = "wes",         icon = "dna", en = "WES / somatic mutations", zh = "WES / 体细胞突变",
         desc_en = "Oncoplot, TMB, signatures, drivers, mutation-survival (maftools).",
         desc_zh = "Oncoplot、TMB、突变特征、驱动基因、突变-预后（maftools）。", anim = "wes",
         ready = TRUE),
    list(v = "spatial",     icon = "border-all", en = "Spatial transcriptomics", zh = "空间转录组",
         desc_en = "Tissue-overlay viz, spatial domains, SVGs, deconvolution (Seurat).",
         desc_zh = "组织叠加图、空间域、空间可变基因、去卷积（Seurat）。", anim = "spatial"),
    list(v = "integration", icon = "layer-group", en = "Multi-omics integration", zh = "多组学整合",
         desc_en = "Integrative subtyping & prognostic stratification (MOFA2/iCluster/SNF).",
         desc_zh = "整合分型与预后分层（MOFA2/iCluster/SNF）。", anim = "integration")
  )
}
