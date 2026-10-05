#' Single-cell pipeline registry (shared by OmicOne and scStudio)
#'
#' Single source of truth for the single-cell steps: value key (= module id),
#' order, phase, bilingual label, module UI function, and `deps` -- the steps
#' whose output this step consumes. `deps` drives staleness (see fct_state.R):
#' when a step re-runs, every finished step that depends on it, directly or
#' transitively, is flagged for a re-run. Output steps (report, export) depend
#' on nothing so they are never flagged.
#'
#' @name steps_sc
#' @keywords internal
NULL

#' Single-cell steps
#' @keywords internal
steps_sc <- function() {
  list(
    list(v = "import",    n = 1,  phase = "sc_data",   en = "Import",        zh = "导入",       ui = mod_import_ui,
         deps = character(0)),
    list(v = "qc",        n = 2,  phase = "sc_data",   en = "Quality control", zh = "质控",     ui = mod_qc_ui,
         deps = "import"),
    list(v = "doublet",   n = 3,  phase = "sc_data",   en = "Doublets",      zh = "去双细胞", ui = mod_doublet_ui,
         deps = c("import", "qc")),
    list(v = "normalize", n = 4,  phase = "sc_prep",   en = "Normalize",     zh = "归一化", ui = mod_normalize_ui,
         deps = c("qc", "doublet")),
    list(v = "reduce",    n = 5,  phase = "sc_prep",   en = "Features / PCA",zh = "特征/PCA",   ui = mod_reduce_ui,
         deps = "normalize"),
    list(v = "integrate", n = 6,  phase = "sc_prep",   en = "Integrate",     zh = "整合",       ui = mod_integrate_ui,
         deps = "reduce"),
    list(v = "cluster",   n = 7,  phase = "sc_struct", en = "Cluster",       zh = "聚类",       ui = mod_cluster_ui,
         deps = c("reduce", "integrate")),
    list(v = "embed",     n = 8,  phase = "sc_struct", en = "Embed",         zh = "降维图", ui = mod_embed_ui,
         deps = c("reduce", "integrate")),
    list(v = "markers",   n = 9,  phase = "sc_id",     en = "Markers",       zh = "标志基因", ui = mod_markers_ui,
         deps = "cluster"),
    list(v = "annotate",  n = 10, phase = "sc_id",     en = "Annotate",      zh = "注释",       ui = mod_annotate_ui,
         deps = c("cluster", "markers")),
    list(v = "enrichment",n = 11, phase = "sc_id",     en = "Enrichment/GSEA", zh = "富集/GSEA", ui = mod_enrichment_ui,
         deps = c("cluster", "annotate")),
    list(v = "pseudobulk",n = 12, phase = "sc_compare", en = "Pseudobulk DE", zh = "Pseudobulk 差异", ui = mod_pseudobulk_ui,
         deps = c("cluster", "annotate")),
    list(v = "abundance", n = 13, phase = "sc_compare", en = "Differential abundance", zh = "细胞组成差异", ui = mod_abundance_ui,
         deps = c("cluster", "annotate")),
    list(v = "trajectory",n = 14, phase = "sc_traj",   en = "Trajectory",    zh = "轨迹",       ui = mod_trajectory_ui,
         deps = c("embed", "cluster", "annotate")),
    list(v = "velocity",  n = 15, phase = "sc_traj",   en = "RNA velocity",  zh = "RNA 速率",   ui = mod_velocity_ui,
         deps = "embed"),
    list(v = "dynamic",   n = 16, phase = "sc_traj",   en = "Dynamic features", zh = "动态特征", ui = mod_dynamic_ui,
         deps = "trajectory"),
    list(v = "cellcycle", n = 17, phase = "sc_adv",    en = "Cell cycle & signatures", zh = "周期与信号", ui = mod_cellcycle_signatures_ui,
         deps = "normalize"),
    list(v = "cellcomm",  n = 18, phase = "sc_adv",    en = "Cell communication", zh = "细胞通讯", ui = mod_cellcomm_ui,
         deps = c("cluster", "annotate")),
    list(v = "malignancy",n = 19, phase = "sc_adv",    en = "Malignant / CNV", zh = "恶性/CNV", ui = mod_malignancy_ui,
         deps = c("normalize", "cluster")),
    list(v = "clinical",  n = 20, phase = "sc_adv",    en = "Clinical & survival", zh = "临床与生存", ui = mod_clinical_ui,
         deps = c("cluster", "annotate", "cellcycle", "malignancy")),
    list(v = "viz",       n = 21, phase = "sc_out",    en = "Visualize",     zh = "可视化", ui = mod_viz_ui,
         deps = c("embed", "cluster", "annotate")),
    list(v = "report",    n = 22, phase = "sc_out",    en = "Report",        zh = "报告",       ui = mod_report_ui,
         deps = character(0)),
    list(v = "export",    n = 23, phase = "sc_out",    en = "Export",        zh = "导出",       ui = mod_export_ui,
         deps = character(0))
  )
}

#' Phase labels (en/zh) for the single-cell registry
#' @keywords internal
phases_sc <- function() {
  list(
    sc_data   = list(en = "Data & QC",  zh = "数据与质控"),
    sc_prep   = list(en = "Preprocess", zh = "预处理"),
    sc_struct = list(en = "Structure",  zh = "结构"),
    sc_id     = list(en = "Identity",   zh = "身份"),
    sc_compare = list(en = "Compare conditions", zh = "条件比较"),
    sc_traj   = list(en = "Trajectory", zh = "轨迹与动态"),
    sc_adv    = list(en = "Advanced",   zh = "高级"),
    sc_out    = list(en = "Output",     zh = "产出")
  )
}
