#' Object helpers shared across modules
#'
#' Thin, defensive wrappers around Seurat/SeuratObject so modules can query the
#' working object without every module repeating null-checks and class handling.
#' Accessors that differ between SeuratObject 4 (slots) and 5 (layers) are
#' resolved here, so the rest of the code never branches on the version.
#'
#' @name fct_object
#' @keywords internal
NULL

#' Dimensions (cells x genes) of the working object
#' @param obj A Seurat object.
#' @return list(cells, genes); NA if unknown.
#' @keywords internal
obj_dims <- function(obj) {
  if (is.null(obj)) return(list(cells = NA_integer_, genes = NA_integer_))
  cells <- tryCatch(ncol(obj), error = function(e) NA_integer_)
  genes <- tryCatch(nrow(obj), error = function(e) NA_integer_)
  list(cells = cells, genes = genes)
}

#' Cell metadata as a data.frame
#' @param obj A Seurat object.
#' @return data.frame (0-col if unavailable).
#' @keywords internal
obj_meta <- function(obj) {
  if (is.null(obj)) return(data.frame())
  md <- tryCatch(obj[[]], error = function(e) NULL)
  if (!is.data.frame(md)) md <- tryCatch(obj@meta.data, error = function(e) data.frame())
  if (!is.data.frame(md)) md <- data.frame()
  md
}

#' Names of metadata columns
#' @param obj A Seurat object.
#' @keywords internal
obj_meta_cols <- function(obj) {
  colnames(obj_meta(obj))
}

#' Available dimensional reductions (e.g. pca, umap, tsne, harmony)
#' @param obj A Seurat object.
#' @keywords internal
obj_reductions <- function(obj) {
  if (is.null(obj)) return(character(0))
  tryCatch(names(obj@reductions), error = function(e) character(0))
}

#' Does the object have a given reduction?
#' @param obj Seurat object; @param red reduction name.
#' @keywords internal
has_reduction <- function(obj, red) {
  red %in% obj_reductions(obj)
}

#' Assay names of the working object (empty if not a Seurat object)
#' @param obj A Seurat object.
#' @keywords internal
obj_assays <- function(obj) {
  if (is.null(obj)) return(character(0))
  tryCatch(SeuratObject::Assays(obj), error = function(e) character(0))
}

#' The default assay, or NULL when it cannot be read
#' @param obj A Seurat object.
#' @keywords internal
obj_default_assay <- function(obj) {
  if (is.null(obj)) return(NULL)
  tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) NULL)
}

#' One entry of `obj@misc`, or NULL
#'
#' OmicOne keeps its own bookkeeping there under `omicone_*` keys: the
#' reduction downstream steps should use (`omicone_reduction`), the active
#' clustering column (`omicone_cluster_col`) and what import standardisation
#' did (`omicone_import_actions`).
#' @param obj A Seurat object. @param key Entry name.
#' @keywords internal
obj_misc <- function(obj, key) {
  if (is.null(obj)) return(NULL)
  tryCatch(obj@misc[[key]], error = function(e) NULL)
}

#' A matrix layer of one assay, for SeuratObject 4 and 5 alike
#'
#' SeuratObject 5 stores matrices as layers (`LayerData()`); version 4 as
#' slots (`GetAssayData(slot = )`). Both are tried, newest first.
#' @param obj A Seurat object.
#' @param layer "counts", "data" or "scale.data".
#' @param assay Assay name; NULL = the default assay.
#' @return The matrix, or NULL if it is not available.
#' @keywords internal
obj_layer <- function(obj, layer = "counts", assay = NULL) {
  assay <- assay %||% obj_default_assay(obj)
  if (is.null(assay)) return(NULL)
  if (exists("LayerData", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    m <- tryCatch(SeuratObject::LayerData(obj, assay = assay, layer = layer),
                  error = function(e) NULL)
    if (!is.null(m)) return(m)
  }
  tryCatch(SeuratObject::GetAssayData(obj, assay = assay, slot = layer),
           error = function(e) NULL)
}

#' Raw counts of an assay (genes x cells)
#' @param obj A Seurat object. @param assay Assay name (default "RNA").
#' @keywords internal
obj_counts <- function(obj, assay = "RNA") {
  m <- obj_layer(obj, "counts", assay)
  if (is.null(m) || !nrow(m)) {
    stop("No raw counts found in assay '", assay, "'.")
  }
  m
}

#' Reductions a neighbour graph or an embedding can be built on
#'
#' Every reduction except the 2-D display maps (UMAP / t-SNE / PaCMAP), in the
#' order PCA, Harmony, Seurat's `integrated.dr`, then anything else (e.g. an
#' Azimuth reference space). Integrated spaces must be offered here, or the
#' CCA / RPCA output could never be used downstream.
#' @param obj A Seurat object.
#' @return Character vector of reduction names (possibly empty).
#' @keywords internal
graph_reduction_choices <- function(obj) {
  reds <- obj_reductions(obj)
  if (!length(reds)) return(character(0))
  linear <- reds[!grepl("umap|tsne|pacmap", tolower(reds))]
  first <- intersect(c("pca", "harmony", "integrated.dr"), linear)
  c(first, setdiff(linear, first))
}

#' The reduction downstream steps should use by default
#'
#' The Integrate step records its output in `obj@misc$omicone_reduction`
#' ("pca" when it ran with no correction); Features/PCA resets it to "pca".
#' @param obj A Seurat object.
#' @return A reduction name, or NULL when there is none.
#' @keywords internal
default_graph_reduction <- function(obj) {
  ch <- graph_reduction_choices(obj)
  if (!length(ch)) return(NULL)
  cur <- obj_misc(obj, "omicone_reduction")
  if (length(cur) == 1 && cur %in% ch) return(cur)
  if ("pca" %in% ch) "pca" else ch[1]
}

#' Metadata columns usable as a sample / batch / grouping column
#'
#' [categorical_cols()], minus per-cell QC metrics and scores by name: on a
#' small dataset an integer metric such as `nFeature_RNA` can have few enough
#' distinct values to pass as categorical, and it must never be offered as a
#' batch.
#' @param md Cell metadata. @param min_levels Minimum distinct values.
#' @keywords internal
sc_group_cols <- function(md, min_levels = 2) {
  cols <- categorical_cols(md, min_levels = min_levels)
  metric <- grepl("^(nCount|nFeature)_|^percent[._]|[._]?[Ss]core$|^doublet_score$", cols)
  cols[!metric]
}

#' Which reduction a rebuilt selector should show (pure)
#'
#' Keeps the user's choice across re-renders, except when the default itself
#' moved (an integration ran, or PCA was recomputed): then the new default
#' wins, so a fresh Harmony space is not hidden behind a stale "pca" choice.
#' @param current Current input value. @param choices New choices.
#' @param default New default. @param prev_default Default at the last render.
#' @keywords internal
reduction_selection <- function(current, choices, default, prev_default = NULL) {
  if (!identical(default, prev_default)) return(default)
  keep_selected(current, choices, default)
}

#' Remove dimensional reductions from an object
#' @param obj A Seurat object. @param reds Reduction names to drop.
#' @keywords internal
drop_reductions <- function(obj, reds) {
  for (r in intersect(reds, obj_reductions(obj))) obj[[r]] <- NULL
  obj
}

#' Build an "RNA" assay from a counts matrix, in the session's assay version
#'
#' Seurat 5 integrates through split layers, which only exist on `Assay5`;
#' Seurat 4 has only the classic `Assay`.
#' @param counts Genes x cells counts matrix.
#' @keywords internal
new_rna_assay <- function(counts) {
  v5 <- exists("CreateAssay5Object", envir = asNamespace("SeuratObject"),
               inherits = FALSE)
  if (v5 && !identical(getOption("Seurat.object.assay.version"), "v3")) {
    return(SeuratObject::CreateAssay5Object(counts = counts))
  }
  SeuratObject::CreateAssayObject(counts = counts)
}

#' Number of counts layers in an assay (1 for SeuratObject 4 assays)
#' @param obj A Seurat object. @param assay Assay name.
#' @keywords internal
n_counts_layers <- function(obj, assay = "RNA") {
  if (!exists("Layers", envir = asNamespace("SeuratObject"), inherits = FALSE)) return(1L)
  a <- tryCatch(obj[[assay]], error = function(e) NULL)
  if (is.null(a) || !methods::is(a, "Assay5")) return(1L)
  length(tryCatch(SeuratObject::Layers(a, search = "counts"),
                  error = function(e) "counts"))
}

#' Bring any imported object to the shape every step assumes
#'
#' Downstream code reads `nCount_RNA` / `nFeature_RNA` and the RNA counts, so
#' an imported object must have: an assay called "RNA" holding raw counts, set
#' as the default assay (Seurat v4 integrated objects default to
#' "integrated"), with one counts layer (`GetAssayData()` aborts on a v5
#' object whose layers are split by sample), and the per-cell totals in the
#' metadata. A missing RNA assay (e.g. "originalexp" after an SCE conversion)
#' is rebuilt from counts rather than renamed, because renaming does not
#' create `nCount_RNA`. What was done is recorded in
#' `obj@misc$omicone_import_actions`, so the import log can replay it.
#'
#' @param obj A Seurat object.
#' @return The standardised Seurat object.
#' @keywords internal
standardize_sc_object <- function(obj) {
  actions <- character(0)
  assays <- obj_assays(obj)
  if (!"RNA" %in% assays) {
    def <- obj_default_assay(obj)
    cand <- unique(c(def, assays))
    src <- NULL
    for (a in cand) {
      if (n_counts_layers(obj, a) > 1) obj[[a]] <- SeuratObject::JoinLayers(obj[[a]])
      m <- obj_layer(obj, "counts", a)
      if (!is.null(m) && nrow(m) > 0 && ncol(m) > 0) {
        src <- a
        break
      }
    }
    if (is.null(src)) {
      stop("The object has no assay with raw counts; OmicOne needs raw UMI counts.")
    }
    obj[["RNA"]] <- new_rna_assay(obj_layer(obj, "counts", src))
    SeuratObject::DefaultAssay(obj) <- "RNA"
    obj <- tryCatch({
      obj[[src]] <- NULL
      obj
    }, error = function(e) obj)
    actions <- c(actions, paste0("rebuild_rna:", src))
  }
  def <- obj_default_assay(obj)
  if (!identical(def, "RNA")) {
    SeuratObject::DefaultAssay(obj) <- "RNA"
    actions <- c(actions, paste0("default_assay:", def))
  }
  if (n_counts_layers(obj, "RNA") > 1) {
    obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
    actions <- c(actions, "join_layers")
  }
  md <- obj_meta(obj)
  if (is.null(md$nCount_RNA) || is.null(md$nFeature_RNA) ||
      "rebuild_rna" %in% sub(":.*$", "", actions)) {
    counts <- obj_counts(obj, "RNA")
    obj$nCount_RNA <- Matrix::colSums(counts)
    obj$nFeature_RNA <- Matrix::colSums(counts > 0)
    actions <- c(actions, "add_ncount")
  }
  obj@misc$omicone_import_actions <- actions
  obj
}

#' Extract a 2D embedding + metadata as a tidy data.frame for plotting
#'
#' @param obj A Seurat object.
#' @param reduction Reduction name (e.g. "umap", "tsne", "pca").
#' @param dims Which two dims to take. Default c(1, 2).
#' @param color_by Optional metadata column to attach as `color`.
#' @return data.frame with columns dim1, dim2, cell, and (optionally) color.
#' @keywords internal
embedding_df <- function(obj, reduction, dims = c(1, 2), color_by = NULL) {
  emb <- tryCatch(
    SeuratObject::Embeddings(obj, reduction = reduction),
    error = function(e) NULL
  )
  if (is.null(emb)) stop("Reduction '", reduction, "' not found in object.")
  df <- data.frame(
    dim1 = emb[, dims[1]],
    dim2 = emb[, dims[2]],
    cell = rownames(emb),
    stringsAsFactors = FALSE
  )
  if (!is.null(color_by)) {
    md <- obj_meta(obj)
    if (color_by %in% colnames(md)) df$color <- md[df$cell, color_by]
  }
  df
}

#' Coerce an uploaded object/matrix into a Seurat object
#'
#' Accepts: a Seurat object (updated if from an old version), a
#' SingleCellExperiment (rebuilt from its `counts` assay and `colData`), a
#' matrix / dgCMatrix / data.frame of counts (genes x cells), or the list a
#' multi-modal 10x file gives (its "Gene Expression" element is used). The
#' result always goes through [standardize_sc_object()].
#'
#' @param x The loaded R object.
#' @param project Project name for a freshly created object.
#' @return A Seurat object.
#' @keywords internal
as_seurat <- function(x, project = "OmicOne") {
  if (!has_pkg("Seurat") || !has_pkg("SeuratObject")) {
    stop("Package 'Seurat' is required to build the working object.")
  }
  cls <- class(x)[1]

  if (is.list(x) && !is.data.frame(x) && !methods::is(x, "Seurat")) {
    x <- pick_gene_expression(x)
  }

  if (methods::is(x, "Seurat")) {
    obj <- tryCatch(SeuratObject::UpdateSeuratObject(x), error = function(e) x)
    return(standardize_sc_object(obj))
  }

  # Rebuilt from counts + colData: Seurat::as.Seurat() would name the assay
  # "originalexp" and leave nCount_RNA undefined.
  if (methods::is(x, "SingleCellExperiment")) {
    if (!has_pkg("SingleCellExperiment")) {
      stop("Package 'SingleCellExperiment' is required to convert this object.")
    }
    if (!"counts" %in% SummarizedExperiment::assayNames(x)) {
      stop("This SingleCellExperiment has no 'counts' assay; OmicOne needs raw counts.")
    }
    m <- SummarizedExperiment::assay(x, "counts")
    if (!methods::is(m, "dgCMatrix")) m <- methods::as(m, "dgCMatrix")
    md <- as.data.frame(SummarizedExperiment::colData(x))
    md <- md[, vapply(md, is.atomic, logical(1)), drop = FALSE]
    obj <- Seurat::CreateSeuratObject(counts = m, project = project,
                                      meta.data = if (ncol(md)) md else NULL)
    return(standardize_sc_object(obj))
  }

  if (methods::is(x, "dgCMatrix") || is.matrix(x) || is.data.frame(x)) {
    m <- x
    if (is.data.frame(m)) m <- as.matrix(m)
    obj <- Seurat::CreateSeuratObject(counts = m, project = project)
    return(standardize_sc_object(obj))
  }

  stop("Unsupported object of class '", cls,
       "'. Expected Seurat, SingleCellExperiment, or a counts matrix/table.")
}

#' The gene-expression matrix of a multi-modal 10x read
#'
#' `Read10X()` / `Read10X_h5()` return a list (one matrix per feature type)
#' when the file also holds antibody or CRISPR features.
#' @param x A list of matrices, or a matrix (returned unchanged).
#' @keywords internal
pick_gene_expression <- function(x) {
  if (!is.list(x) || is.data.frame(x)) return(x)
  if ("Gene Expression" %in% names(x)) return(x[["Gene Expression"]])
  if (length(x) == 1) return(x[[1]])
  stop("This file holds several feature types (", paste(names(x), collapse = ", "),
       ") but none called 'Gene Expression'.")
}

#' Read a 10x HDF5 file, keeping only the gene-expression matrix
#' @param path Path to the `.h5` file.
#' @keywords internal
read_10x_h5 <- function(path) {
  pick_gene_expression(Seurat::Read10X_h5(path))
}

#' Field separator of a counts table, from its file name
#'
#' ".csv" (also ".csv.gz") is comma-separated; ".tsv", ".txt" and anything
#' else is read as tab-separated.
#' @param name File name (the original upload name, not the temp path).
#' @return "," or "\\t".
#' @keywords internal
guess_table_sep <- function(name) {
  n <- tolower(sub("\\.gz$", "", basename(name %||% "")))
  if (grepl("\\.csv$", n)) "," else "\t"
}

#' Read a table (csv/tsv) of counts into a matrix (genes x cells)
#'
#' The first column must hold gene names and every other column one cell's
#' counts; a non-numeric column (e.g. a second annotation column) is an error
#' rather than being silently coerced to NA.
#' @param path File path. @param sep Field separator, or "auto" to infer it
#'   from `name`. @param name File name used for "auto" (defaults to `path`).
#' @keywords internal
read_counts_table <- function(path, sep = "auto", name = path) {
  if (identical(sep, "auto")) sep <- guess_table_sep(name)
  df <- utils::read.delim(path, sep = sep, row.names = 1, check.names = FALSE)
  if (!ncol(df)) stop("The table has no cell columns (check the separator).")
  num <- vapply(df, is.numeric, logical(1))
  if (!all(num)) {
    bad <- names(df)[!num]
    stop("Counts table has non-numeric columns: ",
         paste(utils::head(bad, 5), collapse = ", "),
         if (length(bad) > 5) ", ..." else "",
         ". Expected gene names in the first column and counts in all others.")
  }
  as.matrix(df)
}

#' Reproducible R code for the Import step
#'
#' Produces a Seurat object exactly as [as_seurat()] does, including the
#' standardisation steps that ran (`actions`, from
#' `obj@misc$omicone_import_actions`). The script header defines
#' `input_path`; bundled demos point it at the package file instead.
#'
#' @param kind What was loaded: "seurat", "sce", "matrix" (an .rds holding a
#'   matrix), "h5", "table", or "scop_demo".
#' @param actions Character vector of standardisation actions.
#' @param sep Separator used for a table.
#' @param file Original file name (for a comment).
#' @param demo_file For bundled demos: file name under `inst/extdata`.
#' @param project Project name passed to `CreateSeuratObject()`.
#' @return Character vector of R code lines.
#' @keywords internal
import_log_code <- function(kind, actions = character(0), sep = "\t", file = NULL,
                            demo_file = NULL, project = "OmicOne") {
  head <- character(0)
  if (!is.null(demo_file)) {
    head <- sprintf('input_path <- system.file("extdata", %s, package = "OmicOne")  # bundled demo',
                    r_lit(demo_file))
  } else if (!is.null(file) && !identical(kind, "scop_demo")) {
    head <- sprintf("# input file: %s (set input_path at the top)", file)
  }
  create <- sprintf("obj <- Seurat::CreateSeuratObject(counts = counts, project = %s)",
                    r_lit(project))
  body <- switch(
    kind,
    seurat = c("obj <- readRDS(input_path)",
               "obj <- SeuratObject::UpdateSeuratObject(obj)"),
    sce = c("sce <- readRDS(input_path)",
            'counts <- methods::as(SummarizedExperiment::assay(sce, "counts"), "dgCMatrix")',
            "md <- as.data.frame(SummarizedExperiment::colData(sce))",
            "md <- md[, vapply(md, is.atomic, logical(1)), drop = FALSE]",
            sprintf(paste0("obj <- Seurat::CreateSeuratObject(counts = counts, project = %s,",
                           " meta.data = if (ncol(md)) md else NULL)"), r_lit(project))),
    matrix = c("counts <- readRDS(input_path)", create),
    h5 = c("counts <- Seurat::Read10X_h5(input_path)",
           'if (is.list(counts)) counts <- counts[["Gene Expression"]]',
           create),
    table = c(sprintf(paste0("counts <- as.matrix(utils::read.delim(input_path, sep = %s,",
                             " row.names = 1, check.names = FALSE))"), r_lit(sep)),
              create),
    scop_demo = c('utils::data("pancreas_sub", package = "scop")',
                  "obj <- pancreas_sub"),
    stop("Unknown import kind: ", kind)
  )
  std <- character(0)
  for (a in actions) {
    what <- sub(":.*$", "", a)
    arg <- sub("^[^:]*:", "", a)
    std <- c(std, switch(
      what,
      rebuild_rna = c(
        sprintf('counts <- SeuratObject::LayerData(obj, assay = %s, layer = "counts")', r_lit(arg)),
        'obj[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = counts)',
        'SeuratObject::DefaultAssay(obj) <- "RNA"',
        sprintf("obj[[%s]] <- NULL", r_lit(arg))),
      default_assay = 'SeuratObject::DefaultAssay(obj) <- "RNA"',
      join_layers = 'obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])',
      add_ncount = c(
        'counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")',
        "obj$nCount_RNA <- Matrix::colSums(counts)",
        "obj$nFeature_RNA <- Matrix::colSums(counts > 0)"),
      NULL))
  }
  c(head, body, std)
}
