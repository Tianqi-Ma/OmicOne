#' Module 1: Import & inspect
#'
#' Upload a count matrix and turn it into the working object. Accepts a Seurat or
#' SingleCellExperiment `.rds`, a 10x directory/`.h5`, or a plain counts table
#' (csv/tsv, genes x cells). Old Seurat objects are updated automatically.
#'
#' @param id Module id.
#' @name mod_import
NULL

#' @rdname mod_import
#' @keywords internal
mod_import_ui <- function(id) {
  ns <- shiny::NS(id)
  explainer <- explainer_card(
    title = list(en = "Import your data", zh = "导入数据"),
    what = list(
      en = "Load your single-cell count data into the app as the working object.",
      zh = "将您的单细胞计数数据加载到应用中，作为工作对象。"),
    why  = list(
      en = "Every later step operates on this object. Counts are the raw number of
            transcripts detected per gene per cell.",
      zh = "之后的每一步都在这个对象上进行。计数即每个基因在每个细胞中检测到的原始转录本数量。"),
    how  = list(
      en = "<b>Just exploring?</b> Choose <b>Demo data</b> and click load to try the
            whole pipeline in seconds. To use your own data, pick <b>Upload file</b>
            and the format that matches it (RDS if it's a saved Seurat object).",
      zh = "<b>只是想体验一下？</b>选择<b>演示数据</b>并点击加载，几秒内即可试用整个流程。要使用自己的数据，
            请选择<b>上传文件</b>并选中与之匹配的格式（如果是已保存的 Seurat 对象则选 RDS）。"),
    read = list(
      en = "The overview is drawn once, at import, and does not change when
            later steps filter cells: violins of UMIs, genes and mitochondrial
            % per cell, the 20 genes taking the largest share of all counts,
            and UMIs against genes per cell. Check the matrix loaded the right
            way round (genes in rows) and the numbers are in the expected
            ballpark.",
      zh = "总览在导入时绘制一次，之后的步骤过滤细胞也不会改变它：每细胞 UMI、基因数与线粒体比例的小提琴图，占总计数比例最高的 20 个基因，以及每细胞 UMI 与基因数的散点图。请确认矩阵方向正确（基因为行）、数量级符合预期。"),
    example = list(
      en = "The bundled demo loads instantly with no download. Or upload
               <code>pbmc.rds</code> (a Seurat object), a 10x <code>.h5</code>, or a
               counts table (genes in rows, cells in columns).",
      zh = "内置演示数据无需下载，即刻加载。或上传 <code>pbmc.rds</code>（Seurat 对象）、
               10x <code>.h5</code> 文件，或计数表格（基因为行、细胞为列）。")
  )
  demos <- demo_datasets()
  demo_choices <- stats::setNames(demos$id, paste0(demos$name, "  (", demos$cells, " cells)"))

  controls <- shiny::tagList(
    label_with_help("Data source",
                    "New here? Pick 'Demo data' to try the app instantly. Otherwise upload your own file, or fetch one from a URL.",
                    label_zh = "数据来源",
                    tip_zh = "初次使用？选择“演示数据”即可立即试用应用。否则上传自己的文件，或从网址获取。"),
    shiny::radioButtons(ns("source"), NULL,
                        c("Demo data" = "demo",
                          "Upload file" = "upload",
                          "From URL" = "url"),
                        selected = "demo"),

    # --- Demo data ---
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'demo'", ns("source")),
      label_with_help("Demo dataset",
                      "PBMC 3k (real 10x data) and the tiny example ship with the app and load offline. The pancreas demo needs the scop package.",
                      label_zh = "演示数据集",
                      tip_zh = "PBMC 3k（真实 10x 数据）和小型示例随应用附带，可离线加载。胰腺演示数据需要安装 scop 包。"),
      shiny::selectInput(ns("demo_id"), NULL, choices = demo_choices, selected = "pbmc3k"),
      shiny::helpText(shiny::textOutput(ns("demo_desc"), inline = TRUE)),
      run_button(ns("load_demo"), "Load demo data", "加载演示数据")
    ),

    # --- Upload ---
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'upload'", ns("source")),
      label_with_help("Input format",
                      "RDS = a saved Seurat/SingleCellExperiment object. 10x = Cell Ranger output. Table = a CSV/TSV of counts (genes in rows, cells in columns).",
                      label_zh = "输入格式",
                      tip_zh = "RDS = 已保存的 Seurat/SingleCellExperiment 对象。10x = Cell Ranger 输出。表格 = 计数的 CSV/TSV 文件（基因为行、细胞为列）。"),
      shiny::selectInput(ns("fmt"), NULL,
                         choices = c("RDS (Seurat/SCE)" = "rds",
                                     "10x HDF5 (.h5)"    = "h5",
                                     "Counts table (csv/tsv)" = "table"),
                         selected = "rds"),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'table'", ns("fmt")),
        label_with_help("Separator",
                        "How columns are separated. Auto reads .csv as comma-separated and .tsv/.txt as tab-separated.",
                        label_zh = "分隔符",
                        tip_zh = "表格中各列之间的分隔方式。自动：.csv 按逗号，.tsv/.txt 按制表符读取。"),
        shiny::selectInput(ns("sep"), NULL,
                           choices = c("Auto (from file name)" = "auto", "Tab" = "\t",
                                       "Comma" = ","),
                           selected = "auto")
      ),
      shiny::fileInput(ns("file"), i18n("Choose file", "选择文件"),
                       accept = c(".rds", ".h5", ".csv", ".tsv", ".txt", ".gz")),
      run_button(ns("load"), "Load data", "加载数据")
    ),

    # --- From URL ---
    shiny::conditionalPanel(
      sprintf("input['%s'] == 'url'", ns("source")),
      label_with_help("File URL",
                      "Direct link to a .rds or 10x .h5 file. It is downloaded to a temporary file on your machine.",
                      label_zh = "文件网址",
                      tip_zh = "指向 .rds 或 10x .h5 文件的直接链接。文件会下载到您本机的临时文件中。"),
      shiny::textInput(ns("url"), NULL, placeholder = "https://.../data.h5"),
      shiny::selectInput(ns("url_fmt"), i18n("Format", "格式"),
                         choices = c("RDS (Seurat/SCE)" = "rds", "10x HDF5 (.h5)" = "h5"),
                         selected = "h5"),
      run_button(ns("load_url"), "Fetch & load", "获取并加载")
    )
  )
  tbl_out <- function(id) {
    if (has_pkg("DT")) DT::dataTableOutput(id) else shiny::verbatimTextOutput(id)
  }
  step_container(id = id, 
    title     = list(en = "Import & inspect", zh = "导入与检查"),
    subtitle  = list(en = "Load counts or a Seurat object and inspect the overview.",
                     zh = "载入计数或 Seurat 对象，并查看数据总览。"),
    explainer = explainer,
    controls  = controls,
    summary   = shiny::uiOutput(ns("summary")),
    preview   = bslib::navset_card_tab(
      bslib::nav_panel(i18n("Overview (at import)", "总览（导入时）"),
                       preview_plot_ui(ns("ov_plot"), download = TRUE,
        guide = list(en = "An overview of the dataset as imported will be drawn here.",
                     zh = "运行后，这里将绘制导入时数据集的总览。"),
        caption = list(en = "As imported, before any filtering: per-cell QC violins, top genes by share of counts, UMIs vs genes (one point = one cell).",
                       zh = "导入时、任何过滤之前的数据：每细胞 QC 小提琴图、按计数占比排序的头部基因、UMI 与基因数散点（每个点为一个细胞）。"))),
      bslib::nav_panel(i18n("Cell metadata (current)", "细胞元数据（当前）"),
                       tbl_out(ns("meta_tbl"))),
      bslib::nav_panel(i18n("Counts preview (current)", "表达矩阵预览（当前）"),
                       tbl_out(ns("counts_tbl")))
    )
  )
}

#' @rdname mod_import
#' @param rv Shared reactiveValues hub (with `$obj`, `$source`).
#' @param log_rv reactiveVal reproducibility log.
#' @param parent Parent session (unused; reserved for nav control).
#' @keywords internal
mod_import_server <- function(id, rv, log_rv, parent = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    # Show the selected demo's description under the dropdown.
    output$demo_desc <- shiny::renderText({
      ds <- demo_datasets()
      row <- ds[ds$id == input$demo_id, , drop = FALSE]
      if (nrow(row)) row$description else ""
    })

    # Overview of the dataset as imported (the import module is the documented
    # exception to step_results(): its own epoch bump would wipe it).
    snap <- shiny::reactiveVal(NULL)

    # Put a freshly built object into the hub, start a new epoch, and log the
    # exact code that rebuilds it.
    commit <- function(obj, source_label, log_params, code) {
      start_epoch(rv, "sc", log_rv)
      rv$obj <- obj
      rv$source <- source_label
      snap(list(ov = data_overview(obj),
                plot = tryCatch(overview_plots(obj, when = "at import"),
                                error = function(e) NULL)))
      mark_done(rv, "import")
      log_step(log_rv, "Import", params = log_params, code = code)
      shiny::showNotification(i18n("Data loaded.", "数据已加载。"), type = "message")
    }

    # What was read, for the log: the object family decides the rebuild code.
    import_kind <- function(loaded, fmt) {
      if (fmt %in% c("h5", "table")) return(fmt)
      if (methods::is(loaded, "Seurat")) return("seurat")
      if (methods::is(loaded, "SingleCellExperiment")) return("sce")
      "matrix"
    }

    # Shared loader: read a file of a given format -> Seurat object -> hub.
    load_into_hub <- function(path, fmt, sep = "auto", file_name = basename(path),
                              source_label = fmt, log_params = list(),
                              demo_file = NULL) {
      if (!require_pkgs(c("Seurat", "SeuratObject"), "Import")) return(invisible(NULL))
      sep_used <- if (identical(sep, "auto")) guess_table_sep(file_name) else sep
      got <- with_progress_notify({
        loaded <- switch(
          fmt,
          rds   = readRDS(path),
          h5    = Seurat::Read10X_h5(path),
          table = read_counts_table(path, sep = sep_used)
        )
        kind <- import_kind(loaded, fmt)
        list(obj = as_seurat(loaded), kind = kind)
      }, message = "Loading and building object...")
      if (is.null(got)) return(invisible(NULL))
      code <- import_log_code(got$kind,
                              actions = obj_misc(got$obj, "omicone_import_actions") %||% character(0),
                              sep = sep_used, file = file_name, demo_file = demo_file)
      if (identical(fmt, "table")) log_params$sep <- sep_used
      commit(got$obj, source_label, log_params, code)
    }

    # (a) Upload
    shiny::observeEvent(input$load, {
      shiny::req(input$file)
      fmt <- input$fmt
      file <- input$file
      load_into_hub(file$datapath, fmt, sep = input$sep, file_name = file$name,
                    source_label = switch(fmt, rds = "RDS", h5 = "10x .h5", table = "Table"),
                    log_params = list(source = "upload", format = fmt, file = file$name))
    })

    # (b) Demo data
    shiny::observeEvent(input$load_demo, {
      id <- input$demo_id
      got <- tryCatch(fetch_demo(id), error = function(e) {
        shiny::showNotification(conditionMessage(e), type = "error", duration = 12)
        NULL
      })
      shiny::req(got)
      nm <- demo_datasets()
      label <- nm$name[nm$id == id]
      if (!is.null(got$obj)) {
        # Demo shipped as an object (pancreas_sub from scop): load directly.
        if (!require_pkgs(c("Seurat", "SeuratObject"), "Import")) return(NULL)
        obj <- with_progress_notify(as_seurat(got$obj), message = "Loading demo...")
        if (is.null(obj)) return(NULL)
        code <- import_log_code("scop_demo",
                                actions = obj_misc(obj, "omicone_import_actions") %||% character(0))
        commit(obj, paste0("Demo: ", label), list(source = "demo", demo = id), code)
      } else {
        load_into_hub(got$path, got$format, file_name = got$file,
                      source_label = paste0("Demo: ", label),
                      log_params = list(source = "demo", demo = id),
                      demo_file = got$file)
      }
    })

    # (c) From URL
    shiny::observeEvent(input$load_url, {
      url <- input$url
      fmt <- input$url_fmt
      shiny::req(nzchar(url))
      ext <- if (fmt == "h5") ".h5" else ".rds"
      dest <- tempfile(fileext = ext)
      ok <- with_progress_notify(
        tryCatch(utils::download.file(url, dest, mode = "wb", quiet = TRUE) == 0,
                 error = function(e) FALSE),
        message = "Downloading...")
      if (!isTRUE(ok) || !file.exists(dest) || file.size(dest) == 0) {
        shiny::showNotification(i18n("Download failed (check the URL and your connection).",
                                     "下载失败（请检查网址和网络连接）。"),
                                type = "error", duration = 12)
        return(NULL)
      }
      load_into_hub(dest, fmt, file_name = basename(url), source_label = "URL",
                    log_params = list(source = "url", url = url, format = fmt))
    })

    output$summary <- shiny::renderUI({
      obj <- rv$obj
      if (is.null(obj)) {
        return(shiny::div(class = "omicone-placeholder",
                          i18n("No data yet. Tip: pick <b>Demo data</b> and click <b>Load demo data</b> to try it instantly.",
                               "尚无数据。提示：选择<b>演示数据</b>并点击<b>加载演示数据</b>即可立即试用。")))
      }
      ov <- snap()$ov %||% data_overview(obj)
      fmt <- function(x) if (is.na(x)) "-" else format(round(x), big.mark = ",")
      shiny::tagList(
        stat_tile(i18n("Cells at import", "导入细胞数"), fmt(ov$cells)),
        stat_tile(i18n("Genes", "基因"), fmt(ov$genes)),
        stat_tile(i18n("Median genes/cell", "中位基因/细胞"), fmt(ov$median_genes)),
        stat_tile(i18n("Median UMIs/cell", "中位UMI/细胞"), fmt(ov$median_umi)),
        stat_tile(i18n("Metadata cols", "元数据列"), ov$meta_cols)
      )
    })

    # Overview tab: the survey drawn at import (falls back to the current
    # object, labelled as such, if the data did not come through this module).
    render_step_plot(output, input, "ov_plot", function() {
      shiny::req(rv$obj)
      p <- snap()$plot
      if (!is.null(p)) return(p)
      overview_plots(rv$obj, when = "current object")
    }, name = "import_ov_plot")

    # Cell metadata tab: the real per-cell metadata table.
    output$meta_tbl <- render_tbl_wrap(function() {
      shiny::req(rv$obj)
      md <- obj_meta(rv$obj)
      utils::head(md, 5000)   # cap rows for the browser
    })
    # Counts preview tab: a small dense slice so users see the real matrix.
    output$counts_tbl <- render_tbl_wrap(function() {
      shiny::req(rv$obj)
      counts_preview(rv$obj)
    })
  })
}
