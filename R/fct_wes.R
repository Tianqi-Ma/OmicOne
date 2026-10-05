#' WES / somatic mutation layer — wrappers around maftools
#'
#' Every function here is gated by `require_pkgs("maftools")` at the call site
#' and wrapped in tryCatch, so a missing package or a signature mismatch shows a
#' friendly message instead of crashing the app. maftools installs as a plain
#' Bioconductor binary, so unlike the single-cell engine it needs no source
#' build and runs on a locked-down Windows box.
#'
#' NOTE: maftools renamed and extended a few entry points across versions. The
#' oncogenic-pathway summary is `OncogenicPathways()` up to 2.14 and
#' `pathways()` (with `plotType = "treemap"/"bar"`) after 2.14; `tmb()` gained
#' `captureRegions`/`plotType` and `compareSignatures()` the `"SBS_v34"`
#' database after 2.14. Where that happened the wrapper checks `formals()` or
#' the installed files first, so both generations work.
#'
#' The maftools entry points this pipeline depends on, for checking against your
#' installed version on a first real run:
#' `read.maf`, `getSampleSummary`, `getGeneSummary`, `getClinicalData`,
#' `getFields`, `subsetMaf`, `plotmafSummary`, `oncoplot`, `titv`, `plotTiTv`,
#' `plotVaf`, `rainfallPlot`, `tmb`, `tcgaCompare`, `lollipopPlot`, `oncodrive`,
#' `plotOncodrive`, `somaticInteractions`, `pathways`, `trinucleotideMatrix`,
#' `extractSignatures`, `estimateSignatures`, `plotCophenetic`,
#' `compareSignatures`, `plotSignatures`, `plotApobecDiff`,
#' `clinicalEnrichment`, `plotEnrichmentResults`, `drugInteractions`,
#' `mafCompare`, `coBarplot`, `inferHeterogeneity`, `plotClusters`,
#' `math.score`.
#'
#' @name fct_wes
#' @keywords internal
NULL

#' Which of the maftools entry points this pipeline needs are actually present
#'
#' Used by the tests and useful at the console: a fast way to see whether an
#' installed maftools has drifted away from what the modules call.
#'
#' @return Character vector of missing function names (empty when all present);
#'   `NA_character_` when maftools itself is not installed.
#' @keywords internal
wes_missing_api <- function() {
  if (!has_pkg("maftools")) return(NA_character_)
  needed <- c("read.maf", "getSampleSummary", "getGeneSummary", "getClinicalData",
              "getFields", "subsetMaf", "plotmafSummary", "oncoplot", "titv",
              "plotTiTv", "plotVaf", "rainfallPlot", "tmb", "tcgaCompare",
              "lollipopPlot", "oncodrive", "plotOncodrive", "somaticInteractions",
              "trinucleotideMatrix", "extractSignatures", "estimateSignatures",
              "plotCophenetic", "compareSignatures", "plotSignatures",
              "plotApobecDiff", "clinicalEnrichment", "plotEnrichmentResults",
              "drugInteractions", "mafCompare", "coBarplot", "inferHeterogeneity",
              "plotClusters", "math.score")
  ns <- asNamespace("maftools")
  missing <- needed[!vapply(needed, exists, logical(1), where = ns, inherits = FALSE)]
  # pathways() (after 2.14) vs OncogenicPathways() (2.14 and earlier): either is fine
  if (!exists("pathways", where = ns, inherits = FALSE) &&
      !exists("OncogenicPathways", where = ns, inherits = FALSE)) {
    missing <- c(missing, "pathways / OncogenicPathways")
  }
  missing
}

#' Does an installed maftools function take a given argument?
#'
#' The version guard behind every optional argument (`tmb(plotType =)`,
#' `somaticInteractions(plotPadj =)`, `subsetMaf(ranges =)` ...).
#' @param fn maftools function name. @param arg Argument name.
#' @keywords internal
wes_has_arg <- function(fn, arg) {
  if (!has_pkg("maftools")) return(FALSE)
  f <- tryCatch(get(fn, envir = asNamespace("maftools"), inherits = FALSE),
                error = function(e) NULL)
  is.function(f) && arg %in% names(formals(f))
}

# ---- shared UI / server fragments ---------------------------------------------

#' The "no MAF yet" placeholder every WES step shows before Import has run
#' @keywords internal
wes_no_maf <- function() {
  shiny::div(class = "omicone-placeholder",
             i18n("Load a MAF file on the <b>Import MAF</b> step first.",
                  "请先在<b>导入 MAF</b> 步骤加载 MAF 文件。"))
}

#' The "run this step" prompt shown before a WES step has produced anything
#' @param en,zh Prompt text.
#' @keywords internal
wes_prompt <- function(en, zh) {
  shiny::div(class = "omicone-placeholder", i18n(en, zh))
}

#' Bilingual notification (every WES toast goes through this)
#' @param en,zh Message. @param type,duration Passed to `showNotification()`.
#' @keywords internal
wes_notify <- function(en, zh = en, type = "error", duration = 10) {
  if (!shiny::isRunning()) return(invisible(NULL))
  shiny::showNotification(i18n(en, zh), type = type, duration = duration)
  invisible(NULL)
}

#' A labelled column / option selector (the recurring VAF / protein-change /
#' sample pickers)
#'
#' @param ns Module namespace function. @param id Local input id.
#' @param label,tip Bilingual `list(en =, zh =)`.
#' @param choices Values (optionally named). @param selected Initial value.
#' @param none Optional label for an empty-string first choice (e.g. "(none)").
#' @param multiple,selectize Passed through.
#' @keywords internal
wes_col_select <- function(ns, id, label, tip, choices, selected = NULL, none = NULL,
                           multiple = FALSE, selectize = FALSE) {
  if (!is.null(none)) choices <- c(stats::setNames("", none), choices)
  input <- if (selectize || multiple) {
    shiny::selectizeInput(ns(id), NULL, choices = choices, selected = selected,
                          multiple = multiple)
  } else {
    shiny::selectInput(ns(id), NULL, choices = choices, selected = selected)
  }
  shiny::tagList(
    label_with_help(label$en, tip$en, label_zh = label$zh, tip_zh = tip$zh),
    input)
}

#' Results table of a WES step: capped on screen, complete as a CSV download
#'
#' Returns the three render functions a module assigns explicitly
#' (`output$<id>_slot`, `output$<id>`, `output$<id>_dl`), so the module's
#' output ids stay visible to the output tests. The on-screen table goes through
#' [render_tbl_wrap()] and is capped at `cap` rows; the download has every row.
#'
#' @param ns Module namespace function. @param id Local id of the table output.
#' @param data_fn Zero-argument function returning the full data.frame.
#' @param name File stem of the CSV.
#' @param has_maf,ready Zero-argument predicates: is a MAF loaded / has the step
#'   produced this table yet.
#' @param cap Rows shown in the browser.
#' @param not_ready Bilingual `list(en =, zh =)` shown while `ready()` is FALSE.
#' @return list(slot =, table =, download =).
#' @keywords internal
wes_table <- function(ns, id, data_fn, name, has_maf, ready, cap = 5000,
                      not_ready = list(en = "Run this step to show the results table.",
                                       zh = "运行本步骤后显示结果表。")) {
  slot <- shiny::renderUI({
    if (!isTRUE(has_maf())) return(wes_no_maf())
    if (!isTRUE(ready())) return(wes_prompt(not_ready$en, not_ready$zh))
    n <- tryCatch(nrow(data_fn()), error = function(e) NA_integer_)
    note <- NULL
    if (isTRUE(n > cap)) {
      note <- shiny::div(class = "omicone-muted",
                         i18n(sprintf("Showing the first %s of %s rows; the CSV has all of them.",
                                      format(cap, big.mark = ","), format(n, big.mark = ",")),
                              sprintf("仅显示前 %s 行（共 %s 行）；CSV 文件包含全部行。",
                                      format(cap, big.mark = ","), format(n, big.mark = ","))))
    }
    shiny::tagList(
      shiny::div(class = "omicone-tbl-dl",
                 shiny::downloadButton(ns(paste0(id, "_dl")),
                                       i18n("Download table (.csv)", "下载表格（.csv）"),
                                       class = "btn-sm")),
      note,
      if (has_pkg("DT")) DT::dataTableOutput(ns(id)) else shiny::verbatimTextOutput(ns(id)))
  })
  table <- render_tbl_wrap(function() {
    df <- data_fn()
    if (is.data.frame(df) && nrow(df) > cap) df <- utils::head(df, cap)
    df
  })
  download <- shiny::downloadHandler(
    filename = function() sprintf("%s_%s.csv", name, format(Sys.time(), "%Y%m%d_%H%M%S")),
    content = function(file) utils::write.csv(data_fn(), file, row.names = FALSE))
  list(slot = slot, table = table, download = download)
}

#' Draw a figure once on a null device, so a step can tell a figure that will
#' fail from one that will draw before it marks itself done
#' @param draw Zero-argument draw function.
#' @return TRUE (invisibly); errors propagate.
#' @keywords internal
wes_dry_run <- function(draw) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
  invisible(TRUE)
}

#' Format a number for a pill / sentence ("-" for missing)
#' @param x Number. @param digits Rounding (NULL = integer with separators).
#' @keywords internal
wes_fmt <- function(x, digits = NULL) {
  if (is.null(x) || !length(x) || is.na(x[1])) return("-")
  if (is.null(digits)) return(format(round(x[1]), big.mark = ","))
  formatC(x[1], format = "f", digits = digits, big.mark = ",")
}

#' A positive whole number from a numeric input, or `default` when it is empty
#' @param x Input value. @param default Fallback. @param lo,hi Clamp.
#' @keywords internal
wes_int <- function(x, default, lo = -Inf, hi = Inf) {
  v <- suppressWarnings(as.integer(x %||% NA))
  if (length(v) != 1 || is.na(v)) v <- as.integer(default)
  as.integer(max(lo, min(hi, v)))
}

#' A probability-style threshold from a numeric input, or `default`
#' @param x Input value. @param default Fallback.
#' @keywords internal
wes_prob <- function(x, default) {
  v <- suppressWarnings(as.numeric(x %||% NA))
  if (length(v) != 1 || !is.finite(v) || v <= 0 || v >= 1) v <- default
  v
}

# ---- reproducibility code ---------------------------------------------------

#' Mark a string to be pasted into logged code verbatim (e.g. `wes_raw("NULL")`
#' for an argument that must be passed as an explicit NULL)
#' @param x Code text.
#' @keywords internal
wes_raw <- function(x) structure(x, class = "wes_raw")

#' One line of runnable R code from a function name and the arguments used
#'
#' Values are deparsed, so the logged code carries the parameters that actually
#' ran; pass symbols (`quote(maf)`) for objects defined by earlier log lines and
#' [wes_raw()] for literal code. `NULL` arguments are omitted (the default).
#' @param fn Function name, e.g. "maftools::tmb".
#' @param args Named list of argument values or symbols.
#' @param assign Optional variable name to assign the result to.
#' @return Character scalar.
#' @keywords internal
wes_code <- function(fn, args = list(), assign = NULL) {
  args <- args[!vapply(args, is.null, logical(1))]
  dep <- function(v) {
    if (inherits(v, "wes_raw")) return(unclass(v))
    paste(deparse(v, width.cutoff = 500L), collapse = " ")
  }
  parts <- vapply(names(args), function(k) paste(k, "=", dep(args[[k]])), character(1))
  call <- sprintf("%s(%s)", fn, paste(parts, collapse = ", "))
  if (is.null(assign)) call else paste(assign, "<-", call)
}

# ---- input ------------------------------------------------------------------

#' The default non-synonymous classes of maftools, optionally plus Splice_Region
#' @param splice_region Count `Splice_Region` as non-synonymous.
#' @keywords internal
wes_vc_nonsyn <- function(splice_region = FALSE) {
  vc <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
          "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
          "Missense_Mutation")
  if (isTRUE(splice_region)) c(vc, "Splice_Region") else vc
}

#' Field separator of a delimited file, from its (original) file name
#' @param name File name (the upload's `name`, not its temp path).
#' @keywords internal
wes_sep_for <- function(name) {
  ext <- tolower(sub("\\.gz$", "", basename(name %||% ""), ignore.case = TRUE))
  if (grepl("\\.csv$", ext)) "," else "\t"
}

#' Case-insensitive column lookup
#' @param cols Column names. @param target Name wanted.
#' @return The matching column name, or NULL.
#' @keywords internal
wes_find_col <- function(cols, target) {
  hit <- cols[tolower(cols) == tolower(target)]
  if (length(hit)) hit[1] else NULL
}

#' Read a clinical table without mangling sample ids
#'
#' Everything is read as text first, so ids such as "001" stay "001"; columns
#' other than the id are then type-converted, except those with leading zeros
#' (other ids). The separator follows the file extension (.csv = comma).
#' @param path File path. @param name Original file name (for the extension).
#' @param id_col Column to keep as text (any case); NULL to guess
#'   `Tumor_Sample_Barcode`.
#' @return data.frame; attribute "text_cols" lists the columns kept as text.
#' @keywords internal
wes_read_clinical <- function(path, name = path, id_col = NULL) {
  sep <- wes_sep_for(name)
  df <- utils::read.delim(path, sep = sep, check.names = FALSE, colClasses = "character",
                          na.strings = c("", "NA"), stringsAsFactors = FALSE)
  id <- if (!is.null(id_col)) wes_find_col(names(df), id_col) else NULL
  id <- id %||% wes_find_col(names(df), "Tumor_Sample_Barcode")
  lead0 <- vapply(df, function(v) any(grepl("^0[0-9]+$", v[!is.na(v)])), logical(1))
  text_cols <- union(id, names(df)[lead0])
  conv <- setdiff(names(df), text_cols)
  df[conv] <- lapply(df[conv], utils::type.convert, as.is = TRUE)
  attr(df, "text_cols") <- text_cols
  attr(df, "sep") <- sep
  df
}

#' Standard MAF column names, matched case-insensitively against a header
#' @param header Column names found in the file.
#' @return Named character vector `found_name = standard_name` of the columns
#'   that need renaming (empty when the header is already standard).
#' @keywords internal
wes_maf_renames <- function(header) {
  std <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
           "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
           "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode")
  out <- character(0)
  for (s in std) {
    if (s %in% header) next
    hit <- header[tolower(header) == tolower(s)]
    if (length(hit)) out[hit[1]] <- s
  }
  out
}

#' Where a MAF's header is and which of its columns need renaming
#' @param path MAF path (.gz allowed).
#' @return list(skip = comment lines before the header, renamed = named vector
#'   `found = standard`).
#' @keywords internal
wes_maf_header <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else file(path)
  head_lines <- tryCatch(readLines(con, n = 500, warn = FALSE), finally = close(con))
  comment <- grepl("^#", head_lines)
  skip <- if (any(!comment)) which(!comment)[1] - 1L else 0L
  header <- strsplit(head_lines[skip + 1L], "\t", fixed = TRUE)[[1]]
  list(skip = skip, renamed = wes_maf_renames(header))
}

#' Sample barcodes of a MAF file, read without loading the whole table
#' @param path MAF path. @param hdr [wes_maf_header()] result.
#' @param is_tcga Truncate to 12 characters, as `read.maf(isTCGA = TRUE)` does.
#' @keywords internal
wes_maf_barcodes <- function(path, hdr = wes_maf_header(path), is_tcga = FALSE) {
  col <- names(hdr$renamed)[hdr$renamed == "Tumor_Sample_Barcode"]
  if (!length(col)) col <- "Tumor_Sample_Barcode"
  tsb <- data.table::fread(path, skip = hdr$skip, sep = "\t", quote = "", fill = TRUE,
                           header = TRUE, select = col[1], colClasses = "character")[[1]]
  tsb <- unique(as.character(tsb))
  if (isTRUE(is_tcga)) tsb <- unique(substr(tsb, 1, 12))
  tsb
}

#' Align clinical ids to the MAF's spelling when they differ only in case or
#' surrounding spaces
#'
#' read.maf() joins clinical rows on exact barcodes, so "tcga-ab-2802 " would
#' silently lose its clinical data. Only unambiguous matches are rewritten.
#' @param ids Clinical ids. @param maf_ids MAF barcodes.
#' @return list(ids = aligned ids, n = number rewritten).
#' @keywords internal
wes_align_ids <- function(ids, maf_ids) {
  ids <- as.character(ids)
  key_m <- wes_norm_id(maf_ids)
  unique_key <- !key_m %in% key_m[duplicated(key_m)]
  hit <- match(wes_norm_id(ids), key_m[unique_key])
  fix <- !is.na(hit) & !ids %in% maf_ids
  out <- ids
  out[fix] <- maf_ids[unique_key][hit[fix]]
  list(ids = out, n = sum(fix))
}

#' Read a MAF file into a maftools MAF object
#'
#' MAF column names are matched case-insensitively (some pipelines write
#' `tumor_sample_barcode`); when the header is already standard, the path goes
#' to `read.maf()` unchanged.
#'
#' @param path Path to a `.maf` / `.maf.gz` / tab-delimited mutation table.
#' @param clinical Optional clinical data.frame; must contain
#'   `Tumor_Sample_Barcode`.
#' @param vc_nonSyn Optional character vector overriding which
#'   `Variant_Classification` values count as non-synonymous.
#' @param is_tcga,rm_flags Passed to `read.maf(isTCGA =, rmFlags =)`.
#' @return A maftools MAF object.
#' @keywords internal
wes_read_maf <- function(path, clinical = NULL, vc_nonSyn = NULL, is_tcga = FALSE,
                         rm_flags = FALSE) {
  hdr <- wes_maf_header(path)
  maf_in <- path
  if (length(hdr$renamed)) {
    # same reader settings as read.maf(), whose skip = "Hugo_Symbol" is
    # case-sensitive and therefore cannot find a lower-case header
    maf_in <- data.table::fread(path, skip = hdr$skip, sep = "\t", quote = "", fill = TRUE,
                                header = TRUE, stringsAsFactors = FALSE)
    data.table::setnames(maf_in, names(hdr$renamed), unname(hdr$renamed))
    if ("Tumor_Sample_Barcode" %in% names(maf_in)) {
      maf_in[["Tumor_Sample_Barcode"]] <- as.character(maf_in[["Tumor_Sample_Barcode"]])
    }
  }
  args <- list(maf = maf_in, isTCGA = isTRUE(is_tcga), rmFlags = isTRUE(rm_flags))
  if (!is.null(clinical)) args$clinicalData <- clinical
  if (!is.null(vc_nonSyn)) args$vc_nonSyn <- vc_nonSyn
  do.call(maftools::read.maf, args)
}

#' Runnable R code for the WES import that just ran
#'
#' @param maf_file,clin_file Displayed file names (upload) or NULL.
#' @param demo Was it the bundled demo?
#' @param id_from Original name of the clinical id column (renamed to
#'   `Tumor_Sample_Barcode`), or NULL.
#' @param text_cols Clinical columns kept as text.
#' @param renamed Named vector of MAF header renames (see [wes_maf_renames()]).
#' @param skip Comment lines before the MAF header.
#' @param sep Clinical separator. @param read_args Extra read.maf() arguments.
#' @param aligned Number of clinical ids aligned to the MAF spelling (> 0 adds
#'   the alignment lines).
#' @keywords internal
wes_import_code <- function(maf_file, clin_file, demo, id_from = NULL, text_cols = NULL,
                            renamed = character(0), skip = 0L, sep = "\t",
                            read_args = list(), aligned = 0L) {
  q <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  paths <- if (isTRUE(demo)) {
    c('maf_path <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")',
      if (!is.null(clin_file))
        'clin_path <- system.file("extdata", "tcga_laml_annot.tsv", package = "maftools")')
  } else {
    c(sprintf("maf_path <- %s   # <-- edit: path to your MAF", q(maf_file)),
      if (!is.null(clin_file))
        sprintf("clin_path <- %s   # <-- edit: path to your clinical table", q(clin_file)))
  }
  clin <- if (!is.null(clin_file)) {
    c(sprintf(paste0("clin_raw <- read.delim(clin_path, sep = %s, check.names = FALSE, ",
                     "colClasses = \"character\", na.strings = c(\"\", \"NA\"))"), q(sep)),
      sprintf("conv <- setdiff(names(clin_raw), %s)   # ids stay text", q(text_cols)),
      "clin_raw[conv] <- lapply(clin_raw[conv], utils::type.convert, as.is = TRUE)",
      if (!is.null(id_from) && !identical(id_from, "Tumor_Sample_Barcode")) {
        c(paste0('names(clin_raw)[names(clin_raw) == "Tumor_Sample_Barcode"] <- ',
                 '"Tumor_Sample_Barcode_orig"'),
          sprintf('names(clin_raw)[names(clin_raw) == %s] <- "Tumor_Sample_Barcode"', q(id_from)))
      })
  }
  maf_src <- quote(maf_path)
  reader <- NULL
  if (length(renamed)) {
    # read.maf() looks for a case-sensitive "Hugo_Symbol" header, so the
    # lower-case columns are renamed before it sees the table
    reader <- c(sprintf(paste0("maf_dt <- data.table::fread(maf_path, skip = %d, sep = \"\\t\", ",
                               "quote = \"\", fill = TRUE, header = TRUE)"), as.integer(skip)),
                sprintf("data.table::setnames(maf_dt, %s, %s)", q(unname(names(renamed))),
                        q(unname(renamed))),
                "maf_dt$Tumor_Sample_Barcode <- as.character(maf_dt$Tumor_Sample_Barcode)")
    maf_src <- quote(maf_dt)
  }
  clin_obj <- quote(clin_raw)
  align <- NULL
  if (!is.null(clin_file) && isTRUE(aligned > 0)) {
    tsb_col <- names(renamed)[renamed == "Tumor_Sample_Barcode"]
    tsb_col <- if (length(tsb_col)) tsb_col[1] else "Tumor_Sample_Barcode"
    align <- c(
      sprintf(paste0("maf_tsb <- unique(data.table::fread(maf_path, skip = %d, sep = \"\\t\", ",
                     "quote = \"\", select = %s, colClasses = \"character\")[[1]])"),
              as.integer(skip), q(tsb_col)),
      if (isTRUE(read_args$isTCGA)) "maf_tsb <- unique(substr(maf_tsb, 1, 12))",
      "key <- function(x) toupper(trimws(as.character(x)))",
      "hit <- match(key(clin_raw$Tumor_Sample_Barcode), key(maf_tsb))",
      "clin_for_maf <- clin_raw   # ids differing only in case / spaces take the MAF spelling",
      paste0("clin_for_maf$Tumor_Sample_Barcode <- ifelse(is.na(hit), ",
             "clin_raw$Tumor_Sample_Barcode, maf_tsb[hit])"))
    clin_obj <- quote(clin_for_maf)
  }
  args <- c(list(maf = maf_src), if (!is.null(clin_file)) list(clinicalData = clin_obj),
            read_args)
  c(paths, clin, align, reader, wes_code("maftools::read.maf", args, assign = "maf"))
}

#' Path to maftools' bundled TCGA LAML example (MAF + clinical annotation)
#'
#' Ships with maftools, so the WES pipeline has an instant offline demo in the
#' same spirit as the bundled pbmc3k on the single-cell side. The annotation
#' table (`tcga_laml_annot.tsv`, 200 patients, 7 of them without any MAF
#' record) is the WES demo clinical table: it exercises the survival step's
#' "no MAF record = sequenced wild-type" option without any extra file.
#'
#' @return list(maf=, clinical=); components are "" when maftools is absent.
#' @keywords internal
wes_demo_paths <- function() {
  list(
    maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
    clinical = system.file("extdata", "tcga_laml_annot.tsv", package = "maftools")
  )
}

#' Sample ids in a MAF (every sample with at least one variant, silent included)
#' @param maf A MAF object.
#' @keywords internal
wes_samples <- function(maf) {
  s <- tryCatch(as.character(maftools::getSampleSummary(maf)$Tumor_Sample_Barcode),
                error = function(e) character(0))
  s[!is.na(s)]
}

#' Gene names in a MAF, most-mutated first
#' @param maf A MAF object. @param n How many to return (Inf for all).
#' @keywords internal
wes_genes <- function(maf, n = Inf) {
  g <- tryCatch(as.character(maftools::getGeneSummary(maf)$Hugo_Symbol),
                error = function(e) character(0))
  if (is.finite(n)) utils::head(g, n) else g
}

#' Map typed gene symbols onto the MAF's own spelling (case-insensitive)
#' @param typed Character vector. @param known Gene names in the MAF.
#' @return Character vector of MAF gene names (unknown ones dropped).
#' @keywords internal
wes_match_genes <- function(typed, known) {
  hit <- known[match(toupper(typed), toupper(known))]
  unique(hit[!is.na(hit)])
}

#' Clinical columns carried inside the MAF object
#' @param maf A MAF object.
#' @keywords internal
wes_clinical_cols <- function(maf) {
  cd <- tryCatch(maftools::getClinicalData(maf), error = function(e) NULL)
  if (is.null(cd)) return(character(0))
  setdiff(colnames(cd), "Tumor_Sample_Barcode")
}

#' Columns of the MAF itself (used to find the VAF column)
#' @param maf A MAF object.
#' @keywords internal
wes_fields <- function(maf) {
  tryCatch(as.character(maftools::getFields(maf)), error = function(e) character(0))
}

#' Pick the tumour variant-allele-frequency column out of a set of MAF fields
#'
#' MAFs disagree on this completely: TCGA uses `i_TumorVAF_WU`, others use
#' `VAF`, `tumor_vaf`, `tumor_f`, or only the read counts. Population allele
#' frequencies (ExAC, gnomAD, 1000 Genomes, ESP) look like VAF columns by name
#' but are germline frequencies, so the fuzzy pass never returns them.
#'
#' @param fields Character vector of column names.
#' @return A column name, or NULL.
#' @keywords internal
wes_pick_vaf_col <- function(fields) {
  if (!length(fields)) return(NULL)
  exact <- c("t_vaf", "i_TumorVAF_WU", "VAF", "vaf", "tumor_vaf", "TumorVAF", "tumor_f")
  hit <- exact[exact %in% fields]
  if (length(hit)) return(hit[1])
  population <- grepl("exac|gnomad|1000g|1kg|thousand|esp6500|esp_|dbsnp|popfreq|cosmic|^af$|_af$|^af_",
                      fields, ignore.case = TRUE)
  hit <- fields[grepl("vaf|allele_freq|allele_frac", fields, ignore.case = TRUE) & !population]
  if (length(hit)) return(hit[1])
  NULL
}

#' Can maftools derive the VAF itself (a `t_vaf` column, or the read counts)?
#' @param fields Character vector of column names.
#' @keywords internal
wes_vaf_auto <- function(fields) {
  "t_vaf" %in% fields || all(c("t_ref_count", "t_alt_count") %in% fields)
}

#' Pick the protein-change column out of a set of MAF field names
#' @param fields Character vector of column names.
#' @keywords internal
wes_pick_aa_col <- function(fields) {
  if (!length(fields)) return(NULL)
  exact <- c("HGVSp_Short", "AAChange", "Protein_Change", "amino_acid_change",
             "AAChange.refGene")
  hit <- fields[fields %in% exact]
  if (length(hit)) return(hit[1])
  hit <- grep("aachange|hgvsp|protein_change", fields, ignore.case = TRUE, value = TRUE)
  if (length(hit)) return(hit[1])
  NULL
}

#' Guess the VAF column of a MAF object
#' @param maf A MAF object.
#' @keywords internal
wes_guess_vaf_col <- function(maf) wes_pick_vaf_col(wes_fields(maf))

#' Guess the protein-change column of a MAF object
#' @param maf A MAF object.
#' @keywords internal
wes_guess_aa_col <- function(maf) wes_pick_aa_col(wes_fields(maf))

#' Choices and default of a VAF selector
#'
#' When maftools can compute the VAF itself (`t_vaf`, or `t_ref_count` +
#' `t_alt_count`), the empty choice means "let maftools do it" and is the
#' default unless a named VAF column exists.
#' @param fields MAF field names.
#' @return list(none = label of the empty choice, selected =).
#' @keywords internal
wes_vaf_choice <- function(fields) {
  auto <- wes_vaf_auto(fields)
  guess <- wes_pick_vaf_col(fields)
  if (identical(guess, "t_vaf")) guess <- NULL
  list(none = if (auto) "(auto: t_vaf or t_alt / (t_ref + t_alt))" else "(none)",
       selected = guess %||% "",
       auto = auto)
}

# ---- cohort summaries -------------------------------------------------------

#' Headline numbers for the summary tiles
#' @param maf A MAF object.
#' @return list(samples, genes, variants (non-synonymous), median_per_sample,
#'   top_gene, top_pct).
#' @keywords internal
wes_overview <- function(maf) {
  ss <- tryCatch(maftools::getSampleSummary(maf), error = function(e) NULL)
  gs <- tryCatch(maftools::getGeneSummary(maf), error = function(e) NULL)
  n_samples <- if (!is.null(ss)) nrow(ss) else NA_integer_
  per_sample <- if (!is.null(ss) && "total" %in% colnames(ss)) ss$total else NULL
  list(
    samples  = n_samples,
    genes    = if (!is.null(gs)) nrow(gs) else NA_integer_,
    variants = if (!is.null(per_sample)) sum(per_sample, na.rm = TRUE) else NA_real_,
    median_per_sample = if (!is.null(per_sample)) stats::median(per_sample, na.rm = TRUE)
                        else NA_real_,
    top_gene = if (!is.null(gs) && nrow(gs)) as.character(gs$Hugo_Symbol[1]) else NA_character_,
    top_pct  = if (!is.null(gs) && nrow(gs) && !is.na(n_samples) && n_samples > 0)
                 100 * gs$MutatedSamples[1] / n_samples else NA_real_
  )
}

#' Per-gene mutation frequency table
#' @param maf A MAF object. @param n Number of genes.
#' @return data.frame: Hugo_Symbol, MutatedSamples, total, pct (of all samples).
#' @keywords internal
wes_gene_table <- function(maf, n = 50) {
  gs <- as.data.frame(maftools::getGeneSummary(maf))
  n_samples <- nrow(maftools::getSampleSummary(maf))
  gs <- utils::head(gs, n)
  keep <- intersect(c("Hugo_Symbol", "MutatedSamples", "total"), colnames(gs))
  out <- gs[, keep, drop = FALSE]
  if ("MutatedSamples" %in% keep && n_samples > 0) {
    out$pct <- round(100 * out$MutatedSamples / n_samples, 1)
  }
  out
}

#' Mutated-sample counts of a set of genes, most-mutated first
#' @param maf A MAF object. @param genes Gene names (NULL = all).
#' @return data.frame: gene, mutated, pct.
#' @keywords internal
wes_gene_freq <- function(maf, genes = NULL) {
  gs <- as.data.frame(maftools::getGeneSummary(maf))
  n <- nrow(maftools::getSampleSummary(maf))
  if (!is.null(genes)) gs <- gs[gs$Hugo_Symbol %in% genes, , drop = FALSE]
  ms <- if ("MutatedSamples" %in% colnames(gs)) gs$MutatedSamples else gs$total
  out <- data.frame(gene = as.character(gs$Hugo_Symbol), mutated = ms,
                    pct = if (n > 0) 100 * ms / n else NA_real_, stringsAsFactors = FALSE)
  out[order(-out$mutated), , drop = FALSE]
}

# ---- burden -----------------------------------------------------------------

#' Merge overlapping / abutting intervals of a BED-style table
#' @param bed data.frame with chr, start (0-based), end.
#' @return Merged data.frame (chr, start, end), sorted.
#' @keywords internal
wes_merge_bed <- function(bed) {
  if (!nrow(bed)) return(data.frame(chr = character(0), start = numeric(0), end = numeric(0)))
  pieces <- lapply(split(bed, bed$chr), function(d) {
    d <- d[order(d$start, d$end), , drop = FALSE]
    prev_end <- c(-Inf, cummax(d$end)[-nrow(d)])
    grp <- cumsum(d$start > prev_end)
    data.frame(chr = d$chr[1], start = as.numeric(tapply(d$start, grp, min)),
               end = as.numeric(tapply(d$end, grp, max)), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

#' Read a BED file of capture targets (first three columns), merging overlaps
#'
#' BED is 0-based half-open, so the territory is `sum(end - start)` after
#' merging; vendor BEDs often list overlapping probes, which would otherwise
#' inflate the denominator.
#' @param path File path.
#' @return Merged data.frame (chr, start, end); attribute "n_raw" = input rows.
#' @keywords internal
wes_read_bed <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !grepl("^(#|track|browser)", lines)]
  parts <- strsplit(lines, "[\t ]+")
  pick <- function(i) vapply(parts, function(p) if (length(p) >= i) p[i] else NA_character_,
                             character(1))
  bed <- data.frame(chr = pick(1), start = suppressWarnings(as.numeric(pick(2))),
                    end = suppressWarnings(as.numeric(pick(3))), stringsAsFactors = FALSE)
  bed <- bed[!is.na(bed$chr) & !is.na(bed$start) & !is.na(bed$end) & bed$end > bed$start, ,
             drop = FALSE]
  if (!nrow(bed)) stop("No usable intervals in the BED file (need chrom, start, end).")
  out <- wes_merge_bed(bed)
  attr(out, "n_raw") <- length(lines)
  out
}

#' Territory of a merged BED in megabases
#' @param bed Output of [wes_read_bed()] / [wes_merge_bed()].
#' @keywords internal
wes_bed_size_mb <- function(bed) sum(bed$end - bed$start) / 1e6

#' Do most chromosome names carry a "chr" prefix?
#' @param chrom Chromosome values.
#' @keywords internal
wes_chr_has_prefix <- function(chrom) {
  chrom <- as.character(chrom)
  chrom <- chrom[!is.na(chrom) & nzchar(chrom)]
  length(chrom) > 0 && mean(grepl("^chr", chrom, ignore.case = TRUE)) > 0.5
}

#' Put chromosome names in the MAF's style (with or without "chr")
#' @param chrom Chromosome values. @param prefix Should they carry "chr"?
#' @keywords internal
wes_harmonise_chr <- function(chrom, prefix) {
  bare <- sub("^chr", "", as.character(chrom), ignore.case = TRUE)
  if (isTRUE(prefix)) paste0("chr", bare) else bare
}

#' Tumour mutational burden per sample
#'
#' Non-synonymous mutations divided by the captured territory. With a BED,
#' variants are first restricted to the targets (`subsetMaf(ranges =)`, BED
#' starts shifted to 1-based) and the territory is the merged BED size. Every
#' sample of the MAF is returned: samples left with no variant in the targets
#' come back with TMB 0 instead of silently vanishing from the median.
#'
#' @param maf A MAF object.
#' @param capture_size Captured territory in Mb (ignored when `bed` is given).
#' @param bed Optional merged BED (see [wes_read_bed()]).
#' @param log_scale Passed to `tmb(logScale =)` (plotting only).
#' @return list(df = data.frame(Tumor_Sample_Barcode, total, total_perMB),
#'   maf = the MAF the burden was counted on, capture =, n_zero =).
#' @keywords internal
wes_tmb <- function(maf, capture_size = NULL, bed = NULL, log_scale = TRUE) {
  sub <- maf
  if (!is.null(bed)) {
    if (!wes_has_arg("subsetMaf", "ranges")) {
      stop("This maftools version cannot restrict variants to regions (subsetMaf(ranges =)).")
    }
    has_chr <- wes_chr_has_prefix(maf@data$Chromosome)
    rng <- data.frame(Chromosome = wes_harmonise_chr(bed$chr, has_chr),
                      Start_Position = bed$start + 1, End_Position = bed$end,
                      stringsAsFactors = FALSE)
    sub <- maftools::subsetMaf(maf = maf, ranges = rng)
    capture_size <- wes_bed_size_mb(bed)
  }
  cap <- suppressWarnings(as.numeric(capture_size %||% NA))
  if (length(cap) != 1 || !is.finite(cap) || cap <= 0) {
    stop("Capture size must be a positive number of megabases.")
  }
  args <- list(maf = sub, captureSize = cap, logScale = isTRUE(log_scale))
  if (wes_has_arg("tmb", "plotType")) args$plotType <- NA
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  t <- as.data.frame(do.call(maftools::tmb, args))
  df <- data.frame(Tumor_Sample_Barcode = as.character(t$Tumor_Sample_Barcode),
                   total = as.numeric(t$total), stringsAsFactors = FALSE)
  miss <- setdiff(wes_samples(maf), df$Tumor_Sample_Barcode)
  if (length(miss)) {
    df <- rbind(df, data.frame(Tumor_Sample_Barcode = miss, total = 0,
                               stringsAsFactors = FALSE))
  }
  df$total_perMB <- df$total / cap
  df <- df[order(-df$total_perMB), , drop = FALSE]
  rownames(df) <- NULL
  list(df = df, maf = sub, capture = cap, n_zero = sum(df$total == 0))
}

#' Headline TMB numbers, computed once for the pills and the insight
#' @param df `wes_tmb()$df`. @param ref Reference line (mut/Mb).
#' @keywords internal
wes_tmb_stats <- function(df, ref = 10) {
  v <- df$total_perMB
  list(n = length(v), median = stats::median(v), min = min(v), max = max(v),
       n_high = sum(v >= ref), n_zero = sum(df$total == 0), ref = ref)
}

#' Per-sample TMB, sorted, as a ggplot
#'
#' Zero-mutation samples cannot sit on a log axis; they are drawn as open
#' triangles at the axis floor and counted in the legend, never dropped.
#' @param df `wes_tmb()$df`. @param capture Mb used. @param log_scale Log y.
#' @param ref Reference line (mut/Mb).
#' @keywords internal
wes_tmb_plot <- function(df, capture, log_scale = TRUE, ref = 10) {
  d <- df[order(df$total_perMB), , drop = FALSE]
  d$rank <- seq_len(nrow(d))
  d$kind <- ifelse(d$total == 0, "0 mutations", ">= 1 mutation")
  nz <- d$total_perMB[d$total_perMB > 0]
  floor_y <- if (length(nz)) min(nz) / 2 else 0.01
  d$y <- if (isTRUE(log_scale)) ifelse(d$total == 0, floor_y, d$total_perMB) else d$total_perMB
  med <- stats::median(d$total_perMB)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$rank, y = .data$y, shape = .data$kind)) +
    ggplot2::geom_point(colour = "#3b6ea5", alpha = 0.65, size = 1.9) +
    ggplot2::scale_shape_manual(values = c(">= 1 mutation" = 16, "0 mutations" = 2),
                                name = NULL) +
    ggplot2::geom_hline(yintercept = ref, linetype = 3, colour = "#8b98a5") +
    ggplot2::annotate("text", x = 1, y = ref, hjust = 0, vjust = -0.5, size = 3.3,
                      colour = "#5f6b76",
                      label = "10 mut/Mb: FoundationOne CDx cut-off (KEYNOTE-158), not WES-calibrated") +
    ggplot2::labs(x = sprintf("Samples sorted by TMB (n = %d)", nrow(d)),
                  y = sprintf("Non-synonymous mutations per Mb (capture %s Mb)",
                              format(signif(capture, 4))),
                  title = sprintf("Tumour mutational burden — median %.2f mut/Mb", med)) +
    omicone_theme()
  if (med > 0 || !isTRUE(log_scale)) {
    p <- p + ggplot2::geom_hline(yintercept = med, linetype = 2, colour = "#c1476b")
  }
  if (isTRUE(log_scale)) p <- p + ggplot2::scale_y_log10()
  p
}

#' Runnable R code for the TMB that just ran (with or without a BED)
#' @param capture Mb used. @param bed_name BED file name or NULL.
#' @param maf_prefix Does the MAF use "chr" names? @param log_scale Log axis.
#' @param cohort tcgaCompare() cohort label.
#' @keywords internal
wes_tmb_code <- function(capture, bed_name = NULL, maf_prefix = FALSE, log_scale = TRUE,
                         cohort = "Input") {
  if (is.null(bed_name)) {
    head <- c(wes_code("maftools::tmb", list(maf = quote(maf), captureSize = capture,
                                             logScale = log_scale), assign = "tmb_res"),
              "tmb_maf <- maf",
              sprintf("capture_mb <- %s", format(capture, digits = 10)))
  } else {
    head <- c(
      sprintf("bed_lines <- readLines(%s)   # <-- edit: your kit's BED", deparse(bed_name)),
      "bed_lines <- bed_lines[nzchar(trimws(bed_lines)) & !grepl(\"^(#|track|browser)\", bed_lines)]",
      "bed <- utils::read.table(text = bed_lines, header = FALSE, fill = TRUE)[, 1:3]",
      "names(bed) <- c(\"chr\", \"start\", \"end\")",
      "bed <- bed[bed$end > bed$start, ]",
      "merge_bed <- function(d) {",
      "  d <- d[order(d$start), ]",
      "  prev_end <- c(-Inf, cummax(d$end)[-nrow(d)])",
      "  g <- cumsum(d$start > prev_end)",
      "  data.frame(chr = d$chr[1], start = tapply(d$start, g, min), end = tapply(d$end, g, max))",
      "}",
      "bed <- do.call(rbind, lapply(split(bed, bed$chr), merge_bed))   # overlaps merged",
      if (isTRUE(maf_prefix)) "bed$chr <- paste0(\"chr\", sub(\"^chr\", \"\", bed$chr))"
      else "bed$chr <- sub(\"^chr\", \"\", bed$chr)",
      "capture_mb <- sum(bed$end - bed$start) / 1e6   # BED is 0-based, half-open",
      paste0("tmb_maf <- maftools::subsetMaf(maf = maf, ranges = data.frame(Chromosome = bed$chr, ",
             "Start_Position = bed$start + 1, End_Position = bed$end))"),
      wes_code("maftools::tmb", list(maf = quote(tmb_maf), captureSize = quote(capture_mb),
                                     logScale = log_scale), assign = "tmb_res"))
  }
  c(head,
    "all_tsb <- as.character(maftools::getSampleSummary(maf)$Tumor_Sample_Barcode)",
    "n_absent <- length(setdiff(all_tsb, as.character(tmb_res$Tumor_Sample_Barcode)))",
    "tmb_all <- c(tmb_res$total_perMB, rep(0, n_absent))   # samples with no variant left = 0",
    "median(tmb_all)",
    "sum(tmb_all >= 10)   # FoundationOne CDx cut-off (KEYNOTE-158), not WES-calibrated",
    wes_code("maftools::tcgaCompare",
             list(maf = quote(tmb_maf), cohortName = cohort, capture_size = quote(capture_mb),
                  logscale = log_scale, rm_zero = TRUE)))
}

#' The 33 TCGA cohort codes tcgaCompare() draws
#' @keywords internal
wes_tcga_codes <- function() {
  f <- if (has_pkg("maftools")) system.file("extdata", "tcga_cohort.txt.gz", package = "maftools")
       else ""
  codes <- if (nzchar(f) && file.exists(f)) {
    tryCatch(unique(as.character(data.table::fread(f, select = "cohort")$cohort)),
             error = function(e) NULL)
  }
  codes %||% c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM",
               "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
               "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD",
               "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
}

#' Cohort label for tcgaCompare(), renamed when it collides with a TCGA code
#'
#' tcgaCompare() colours every cohort whose name equals `cohortName` as the
#' input, so an input called "ESCA" would paint TCGA-ESCA as yours too.
#' @param source Import label / file name. @param codes TCGA codes.
#' @keywords internal
wes_tcga_label <- function(source, codes = wes_tcga_codes()) {
  lab <- source %||% "This cohort"
  lab <- sub("\\.(maf|maf\\.gz|txt|tsv|csv|gz)$", "", basename(lab), ignore.case = TRUE)
  lab <- trimws(gsub("[()]", "", sub("(?i)\\bdemo\\b\\s*:?", "", lab, perl = TRUE)))
  if (!nzchar(lab)) lab <- "This cohort"
  if (nchar(lab) > 18) lab <- paste0(substr(lab, 1, 17), "~")
  if (toupper(lab) %in% toupper(codes)) lab <- paste0(lab, "_input")
  lab
}

#' Transition / transversion summary
#' @param maf A MAF object. @param use_syn Include synonymous variants.
#' @keywords internal
wes_titv <- function(maf, use_syn = TRUE) {
  maftools::titv(maf = maf, plot = FALSE, useSyn = use_syn)
}

#' Pooled share of the six base changes and of transitions
#'
#' Counts are summed over the cohort before dividing, so a sample with three
#' SNVs does not weigh as much as one with three hundred (the per-sample mean
#' of shares does exactly that).
#' @param tv A [wes_titv()] result.
#' @return Named numeric vector of percentages (six classes, "Ti", "Tv"), plus
#'   attribute "n_snv".
#' @keywords internal
wes_titv_pooled <- function(tv) {
  rc <- as.data.frame(tv$raw.counts)
  cls <- intersect(c("C>A", "C>G", "C>T", "T>C", "T>A", "T>G"), colnames(rc))
  tot <- colSums(rc[cls], na.rm = TRUE)
  n <- sum(tot)
  pct <- if (n > 0) 100 * tot / n else tot * NA
  ti <- sum(pct[intersect(c("C>T", "T>C"), names(pct))])
  out <- c(pct, Ti = ti, Tv = if (n > 0) 100 - ti else NA_real_)
  attr(out, "n_snv") <- n
  out
}

# ---- lollipop ---------------------------------------------------------------

#' Amino-acid positions from protein-change strings, parsed as lollipopPlot() does
#'
#' Handles HGVSp short / long ("p.R882H", "p.Arg882His", "p.C229Lfs*18",
#' "p.Asn1986GlnfsTer13") and ranges ("p.761_762del" -> 761).
#' @param x Protein-change strings.
#' @return Numeric positions (NA where none can be read).
#' @keywords internal
wes_aa_position <- function(x) {
  conv <- vapply(strsplit(as.character(x), ".", fixed = TRUE),
                 function(p) if (length(p)) p[length(p)] else NA_character_, character(1))
  pos <- gsub("Ter.*", "", conv)
  pos <- gsub("[[:alpha:]]", "", pos)
  pos <- gsub("\\*$", "", pos)
  pos <- gsub("^\\*", "", pos)
  pos <- gsub("\\*.*", "", pos)
  pos <- vapply(strsplit(pos, "_", fixed = TRUE),
                function(p) if (length(p)) p[1] else NA_character_, character(1))
  abs(suppressWarnings(as.numeric(pos)))
}

#' Protein positions of a gene hit at least `min_n` times (non-synonymous)
#' @param maf A MAF object. @param gene Gene symbol.
#' @param aa_col Protein-change column. @param min_n Minimum occurrences.
#' @return list(positions = recurrent positions, all = every parsed position).
#' @keywords internal
wes_recurrent_positions <- function(maf, gene, aa_col, min_n = 2) {
  d <- as.data.frame(maftools::subsetMaf(maf = maf, genes = gene, includeSyn = FALSE,
                                         mafObj = FALSE))
  if (!nrow(d) || is.null(aa_col) || !aa_col %in% names(d)) {
    return(list(positions = numeric(0), all = numeric(0)))
  }
  pos <- wes_aa_position(d[[aa_col]])
  pos <- pos[!is.na(pos)]
  tab <- table(pos)
  list(positions = sort(as.numeric(names(tab)[tab >= min_n])), all = pos)
}

#' Transcripts maftools has protein domains for, longest first
#' @param gene Gene symbol.
#' @return data.frame: refseq.ID, protein.ID, aa.length (possibly empty).
#' @keywords internal
wes_gene_transcripts <- function(gene) {
  f <- if (has_pkg("maftools")) system.file("extdata", "protein_domains.RDs", package = "maftools")
       else ""
  empty <- data.frame(refseq.ID = character(0), protein.ID = character(0),
                      aa.length = numeric(0), stringsAsFactors = FALSE)
  if (!nzchar(f) || !file.exists(f) || is.null(gene) || !nzchar(gene)) return(empty)
  pd <- as.data.frame(readRDS(f))
  pd <- pd[pd$HGNC %in% gene, c("refseq.ID", "protein.ID", "aa.length"), drop = FALSE]
  pd <- pd[!duplicated(pd$refseq.ID), , drop = FALSE]
  pd[order(-pd$aa.length), , drop = FALSE]
}

# ---- drivers ----------------------------------------------------------------

#' Positional clustering driver detection (oncodrive)
#'
#' maftools falls back to a preset background (mean 0.279, SD 0.13) when the
#' cohort has too few synonymous variants to build one; that changes every
#' p-value, so the message is captured and shown instead of lost in the console.
#' @param maf A MAF object. @param aa_col Protein-change column.
#' @param min_mut Minimum mutations per gene to test.
#' @return list(res = oncodrive table, bg_note = character).
#' @keywords internal
wes_oncodrive <- function(maf, aa_col = NULL, min_mut = 5) {
  notes <- character(0)
  res <- withCallingHandlers(
    maftools::oncodrive(maf = maf, AACol = aa_col, minMut = min_mut, pvalMethod = "zscore"),
    message = function(m) {
      txt <- trimws(conditionMessage(m))
      if (grepl("predefined values", txt, ignore.case = TRUE)) notes <<- c(notes, txt)
    })
  list(res = as.data.frame(res), bg_note = unique(notes))
}

#' Mutually exclusive / co-occurring gene pairs, with BH-adjusted p-values
#' @param maf A MAF object. @param top Number of top genes to test.
#' @return data.frame, one row per gene pair (pValue, pAdj, Event, counts).
#' @keywords internal
wes_interactions <- function(maf, top = 25) {
  args <- list(maf = maf, top = top, pvalue = c(0.05, 0.1))
  if (wes_has_arg("somaticInteractions", "plotPadj")) args$plotPadj <- TRUE
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  as.data.frame(do.call(maftools::somaticInteractions, args))
}

#' Oncogenic pathway summary, across the maftools rename
#'
#' `OncogenicPathways()` (2.14 and earlier) draws pathway / fraction of the
#' pathway's genes mutated / fraction of samples affected; after 2.14
#' `pathways(plotType = "bar")` draws the same three panels.
#' @param maf A MAF object.
#' @keywords internal
wes_pathways <- function(maf) {
  if (exists("pathways", where = asNamespace("maftools"), inherits = FALSE)) {
    return(maftools::pathways(maf = maf, plotType = "bar"))
  }
  fn <- utils::getFromNamespace("OncogenicPathways", "maftools")
  fn(maf = maf)
}

# ---- signatures -------------------------------------------------------------

#' The BSgenome package a reference build needs
#' @param build "hg19" or "hg38".
#' @keywords internal
wes_bsgenome_pkg <- function(build = c("hg19", "hg38")) {
  build <- match.arg(build)
  if (build == "hg19") "BSgenome.Hsapiens.UCSC.hg19" else "BSgenome.Hsapiens.UCSC.hg38"
}

#' Reference build from MAF `NCBI_Build` values
#' @param x Values of the NCBI_Build column.
#' @return "hg19", "hg38" or NA.
#' @keywords internal
wes_build_from_ncbi <- function(x) {
  v <- toupper(trimws(as.character(x)))
  v <- v[!is.na(v) & nzchar(v)]
  if (!length(v)) return(NA_character_)
  v <- names(sort(table(v), decreasing = TRUE))[1]
  if (v %in% c("37", "GRCH37", "HG19", "GRCH37-LITE", "B37")) return("hg19")
  if (v %in% c("38", "GRCH38", "HG38")) return("hg38")
  NA_character_
}

#' Reference build declared by a MAF
#' @param maf A MAF object.
#' @keywords internal
wes_guess_build <- function(maf) {
  x <- tryCatch(maf@data$NCBI_Build, error = function(e) NULL)
  if (is.null(x)) return(NA_character_)
  wes_build_from_ncbi(x)
}

#' Check MAF reference alleles against the chosen BSgenome
#'
#' A wrong build still produces a 96-channel matrix — just of the wrong
#' contexts. Up to `n` evenly spaced SNVs are looked up; skipped (NULL) when
#' BSgenome or the build's package is not installed.
#' @param maf A MAF object. @param build "hg19"/"hg38". @param n SNVs to check.
#' @return list(n, concordance, missing_chr) or NULL.
#' @keywords internal
wes_ref_check <- function(maf, build, n = 500) {
  pkg <- wes_bsgenome_pkg(build)
  if (!has_pkg("BSgenome") || !has_pkg(pkg)) return(NULL)
  d <- as.data.frame(maf@data)
  snv <- d[d$Variant_Type %in% "SNP" & d$Reference_Allele %in% c("A", "C", "G", "T"), ,
           drop = FALSE]
  if (!nrow(snv)) return(NULL)
  idx <- unique(round(seq(1, nrow(snv), length.out = min(n, nrow(snv)))))
  snv <- snv[idx, , drop = FALSE]
  genome <- BSgenome::getBSgenome(pkg)
  chr <- as.character(snv$Chromosome)
  if (!wes_chr_has_prefix(chr)) chr <- paste0("chr", chr)
  pos <- as.numeric(snv$Start_Position)
  ok <- chr %in% BSgenome::seqnames(genome) & is.finite(pos)
  if (!any(ok)) return(list(n = 0L, concordance = 0, missing_chr = sum(!ok)))
  ref <- BSgenome::getSeq(genome, names = chr[ok], start = pos[ok], end = pos[ok],
                          as.character = TRUE)
  list(n = sum(ok), concordance = mean(toupper(ref) == snv$Reference_Allele[ok]),
       missing_chr = sum(!ok))
}

#' Trinucleotide context matrix (the input to signature extraction)
#'
#' UCSC BSgenomes name chromosomes "chr1"; a "chr" prefix is added only when
#' the MAF does not already carry one (adding it twice drops every variant).
#' @param maf A MAF object. @param build "hg19"/"hg38".
#' @return The trinucleotideMatrix() result; attribute "prefix_added".
#' @keywords internal
wes_trinuc <- function(maf, build = "hg19") {
  add <- !wes_chr_has_prefix(maf@data$Chromosome)
  args <- list(maf = maf, ref_genome = wes_bsgenome_pkg(build))
  if (add) {
    args$prefix <- "chr"
    args$add <- TRUE
  }
  tnm <- do.call(maftools::trinucleotideMatrix, args)
  attr(tnm, "prefix_added") <- add
  tnm
}

#' The newest COSMIC SBS database the installed maftools ships
#' @keywords internal
wes_sig_db <- function() {
  v34 <- if (has_pkg("maftools")) {
    system.file("extdata", "SBS_v34_signatures.RDs", package = "maftools")
  } else ""
  if (nzchar(v34) && file.exists(v34)) "SBS_v34" else "SBS"
}

#' Run an NMF-based maftools call, adding pConstant only if it needs it
#'
#' maftools documents `pConstant` as a fix for the "non-conformable arrays"
#' error; NMF itself refuses a matrix with an all-zero channel ("null or
#' NA-filled row"), which low-burden exome cohorts often have. Only those two
#' errors trigger a retry with 1e-4; adding it unconditionally shifts every
#' count.
#' @param fn maftools function. @param args Argument list.
#' @return list(res =, pconstant = NULL or 1e-4).
#' @keywords internal
wes_nmf_retry <- function(fn, args) {
  if (!"package:NMF" %in% search()) suppressMessages(attachNamespace("NMF"))
  # drop arguments this maftools version does not have (formals guard)
  args <- args[names(args) %in% names(formals(fn))]
  first <- tryCatch(list(res = do.call(fn, args), pconstant = NULL),
                    error = function(e) e)
  if (!inherits(first, "error")) return(first)
  if (!grepl("non-conformable|null or NA-filled", conditionMessage(first), ignore.case = TRUE)) {
    stop(first)
  }
  args$pConstant <- 1e-4
  list(res = do.call(fn, args), pconstant = 1e-4)
}

#' Extract de-novo mutational signatures and match them to COSMIC
#' @param tnm A [wes_trinuc()] result. @param n Number of signatures.
#' @return list(sig =, cmp =, db =, pconstant =).
#' @keywords internal
wes_signatures <- function(tnm, n = 3) {
  # NMF (0.28, current CRAN) resolves its internal `seed` S4 generic with a
  # bare getGeneric("seed"), which returns NULL once the Bioconductor stack —
  # attached by trinucleotideMatrix's BSgenome load — is on the search path
  # while NMF itself is not. Attaching NMF (in wes_nmf_retry) restores it.
  # parallel = NULL: NMF's parallel mode (.opt = "P4") dies inside Shiny where
  # no foreach backend is registered; sequential NMF is seconds on this matrix.
  out <- wes_nmf_retry(maftools::extractSignatures,
                       list(mat = tnm, n = n, parallel = NULL))
  db <- wes_sig_db()
  cmp <- tryCatch(maftools::compareSignatures(nmfRes = out$res, sig_db = db, verbose = FALSE),
                  error = function(e) NULL)
  list(sig = out$res, cmp = cmp, db = db, pconstant = out$pconstant)
}

#' Estimate the number of signatures (cophenetic correlation over ranks)
#' @param tnm A [wes_trinuc()] result. @param n_try Largest rank tried.
#' @param nrun NMF runs per rank.
#' @return list(res =, pconstant =).
#' @keywords internal
wes_estimate_rank <- function(tnm, n_try = 6, nrun = 10) {
  ns <- asNamespace("maftools")
  if (!exists("estimateSignatures", where = ns, inherits = FALSE) ||
      !exists("plotCophenetic", where = ns, inherits = FALSE)) {
    stop("This maftools version has no estimateSignatures() / plotCophenetic().")
  }
  args <- list(mat = tnm, nMin = 2, nTry = n_try, nrun = nrun, parallel = NULL)
  if (wes_has_arg("estimateSignatures", "verbose")) args$verbose <- FALSE
  # estimateSignatures() draws its cophenetic plot as a side effect
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  wes_nmf_retry(maftools::estimateSignatures, args)
}

#' Best COSMIC match per extracted signature, with its cosine similarity
#' @param cmp compareSignatures() result.
#' @return data.frame: signature, best_match, cosine, aetiology.
#' @keywords internal
wes_sig_matches <- function(cmp) {
  cs <- cmp$cosine_similarities
  if (is.null(cs) || !length(cs)) return(NULL)
  cs <- as.matrix(cs)
  best <- colnames(cs)[apply(cs, 1, which.max)]
  ae <- cmp$aetiology_db
  aet <- if (!is.null(ae)) {
    ae <- as.data.frame(ae)
    vals <- ae[[ncol(ae)]]
    keys <- if (!is.null(rownames(ae)) && all(best %in% rownames(ae))) rownames(ae) else ae[[1]]
    as.character(vals[match(best, keys)])
  } else rep(NA_character_, length(best))
  data.frame(signature = rownames(cs), best_match = best,
             cosine = round(apply(cs, 1, max), 3), aetiology = aet,
             stringsAsFactors = FALSE, row.names = NULL)
}

# ---- clinical / comparison --------------------------------------------------

#' Clinical columns fit for an enrichment test
#'
#' 2–`max_levels` non-missing levels, each carried by at least `min_n` samples:
#' with fewer, the per-level Fisher test has nothing to say.
#' @param cd Clinical data.frame. @param max_levels,min_n Limits.
#' @keywords internal
wes_enrichment_features <- function(cd, max_levels = 10, min_n = 3) {
  cd <- as.data.frame(cd)
  cols <- setdiff(categorical_cols(cd, min_levels = 2, max_levels = max_levels,
                                   allow_na = TRUE), "Tumor_Sample_Barcode")
  cols[vapply(cols, function(cl) {
    v <- as.character(cd[[cl]])
    v <- v[!is.na(v) & nzchar(v)]
    tab <- table(v)
    length(tab) >= 2 && all(tab >= min_n)
  }, logical(1))]
}

#' Genes enriched in one level of a clinical feature
#'
#' Samples with a missing value are removed and the rest passed as
#' `annotationDat`: otherwise clinicalEnrichment() counts them in every
#' level's "rest" group.
#' @param maf A MAF object. @param feature Clinical column name.
#' @param min_mut Genes need more than this many mutated samples.
#' @return clinicalEnrichment() result; attributes "n_used", "n_missing".
#' @keywords internal
wes_clinical_enrichment <- function(maf, feature, min_mut = 5) {
  cd <- as.data.frame(maftools::getClinicalData(maf))
  if (!feature %in% names(cd)) stop("Clinical column not found: ", feature)
  v <- as.character(cd[[feature]])
  keep <- !is.na(v) & nzchar(trimws(v))
  anno <- data.frame(Tumor_Sample_Barcode = as.character(cd$Tumor_Sample_Barcode[keep]),
                     v[keep], stringsAsFactors = FALSE)
  names(anno)[2] <- feature
  res <- maftools::clinicalEnrichment(maf = maf, clinicalFeature = feature,
                                      annotationDat = anno, minMut = min_mut)
  attr(res, "n_used") <- sum(keep)
  attr(res, "n_missing") <- sum(!keep)
  res
}

#' Groupwise enrichment table with BH q-values
#' @param enr clinicalEnrichment() result.
#' @keywords internal
wes_enrichment_table <- function(enr) {
  d <- as.data.frame(enr$groupwise_comparision)
  if (!"fdr" %in% names(d)) d$fdr <- stats::p.adjust(d$p_value, method = "BH")
  d[order(d$fdr, d$p_value), , drop = FALSE]
}

#' Significant enrichment results, counted in genes rather than rows
#'
#' One gene can be enriched in several levels; the table has a row per
#' gene x level. Only OR > 1 rows are what plotEnrichmentResults() can draw.
#' @param enr clinicalEnrichment() result. @param q FDR cutoff.
#' @return list(tested_genes, sig_genes, sig_rows, drawn (data.frame of rows
#'   with q < cutoff, OR > 1 and p < 0.05)).
#' @keywords internal
wes_enrichment_sig <- function(enr, q = 0.05) {
  d <- wes_enrichment_table(enr)
  sig <- d[!is.na(d$fdr) & d$fdr < q, , drop = FALSE]
  drawn <- sig[!is.na(sig$OR) & sig$OR > 1 & sig$p_value < 0.05, , drop = FALSE]
  list(tested_genes = length(unique(d$Hugo_Symbol)),
       sig_genes = length(unique(sig$Hugo_Symbol)),
       sig_rows = nrow(sig), sig = sig, drawn = drawn,
       genes = unique(as.character(sig$Hugo_Symbol)))
}

#' Feature levels usable for a two-group comparison (>= `min_n` samples each)
#' @param cd Clinical data.frame. @param feature Column. @param min_n Minimum.
#' @keywords internal
wes_feature_levels <- function(cd, feature, min_n = 2) {
  v <- as.character(as.data.frame(cd)[[feature]])
  tab <- sort(table(v[!is.na(v) & nzchar(v)]), decreasing = TRUE)
  names(tab)[tab >= min_n]
}

#' Compare two cohorts defined by a clinical column's levels
#'
#' @param maf A MAF object.
#' @param feature Clinical column to split on.
#' @param level1,level2 The two levels to compare.
#' @param min_mut A gene is tested if it is mutated in at least this many
#'   samples of either group (mafCompare's own rule).
#' @return list(m1=, m2=, res=, n1=, n2=) — the two sub-MAFs, the test, and the
#'   sample sizes mafCompare actually used.
#' @keywords internal
wes_compare_cohorts <- function(maf, feature, level1, level2, min_mut = 5) {
  cd <- as.data.frame(maftools::getClinicalData(maf))
  v <- as.character(cd[[feature]])
  s1 <- as.character(cd$Tumor_Sample_Barcode[which(v == level1)])
  s2 <- as.character(cd$Tumor_Sample_Barcode[which(v == level2)])
  if (length(s1) < 2 || length(s2) < 2) {
    stop("Each group needs at least 2 samples (got ", length(s1), " and ",
         length(s2), ").")
  }
  m1 <- maftools::subsetMaf(maf = maf, tsb = s1, mafObj = TRUE)
  m2 <- maftools::subsetMaf(maf = maf, tsb = s2, mafObj = TRUE)
  res <- maftools::mafCompare(m1 = m1, m2 = m2, m1Name = level1, m2Name = level2,
                              minMut = min_mut)
  ss <- as.data.frame(res$SampleSummary)
  list(m1 = m1, m2 = m2, res = res,
       n1 = as.numeric(ss$SampleSize[1]), n2 = as.numeric(ss$SampleSize[2]))
}

#' Forest-plot data from a mafCompare() result, significance by BH q
#'
#' Fisher's p and q come from the uncorrected 2x2 table. A zero cell makes the
#' odds ratio 0 or Inf; for drawing only, those rows get a pseudo-count of 1 in
#' every cell (mafCompare's `pseudoCount` rule) and are flagged.
#' @param res mafCompare() result. @param q FDR cutoff. @param n1,n2 Group sizes.
#' @keywords internal
wes_forest_data <- function(res, q = 0.05, n1, n2) {
  d <- as.data.frame(res$results)
  if (!"adjPval" %in% names(d)) d$adjPval <- stats::p.adjust(d$pval, method = "BH")
  d <- d[!is.na(d$adjPval) & d$adjPval < q, , drop = FALSE]
  if (!nrow(d)) return(d)
  a <- d[[2]]
  b <- d[[3]]
  zero <- a == 0 | b == 0 | a == n1 | b == n2
  or_pc <- ((a + 1) * (n2 - b + 1)) / ((n1 - a + 1) * (b + 1))
  se_pc <- sqrt(1 / (a + 1) + 1 / (n1 - a + 1) + 1 / (b + 1) + 1 / (n2 - b + 1))
  d$or_plot <- ifelse(zero, or_pc, d$or)
  d$lo_plot <- ifelse(zero, exp(log(or_pc) - 1.96 * se_pc), d$ci.low)
  d$hi_plot <- ifelse(zero, exp(log(or_pc) + 1.96 * se_pc), d$ci.up)
  d$pseudo <- zero
  d[order(d$adjPval), , drop = FALSE]
}

#' Forest plot (log odds-ratio axis) of the genes passing the FDR cutoff
#' @param fd [wes_forest_data()] output. @param l1,l2 Group names.
#' @param n1,n2 Group sizes. @param q FDR cutoff.
#' @keywords internal
wes_forest_plot <- function(fd, l1, l2, n1, n2, q) {
  if (!nrow(fd)) {
    stop(sprintf("No gene passes FDR < %g. The Results tab lists every tested gene.", q))
  }
  fd$label <- sprintf("%s  (%d/%d vs %d/%d)", fd$Hugo_Symbol, as.integer(fd[[2]]),
                      as.integer(n1), as.integer(fd[[3]]), as.integer(n2))
  fd$label <- factor(fd$label, levels = rev(fd$label))
  fd$kind <- ifelse(fd$pseudo, "zero cell: +1 pseudo-count (display only)",
                    "Fisher estimate")
  ggplot2::ggplot(fd, ggplot2::aes(x = .data$or_plot, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "#8b98a5") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo_plot, xmax = .data$hi_plot),
                           width = 0.25, orientation = "y", colour = "#3b6ea5") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$kind), size = 2.6, colour = "#3b6ea5") +
    ggplot2::scale_shape_manual(values = c("Fisher estimate" = 16,
                                           "zero cell: +1 pseudo-count (display only)" = 1),
                                name = NULL) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = sprintf("Odds ratio, %s vs %s (log scale; 95%% CI)", l1, l2),
                  y = NULL,
                  title = sprintf("Genes at FDR < %g (Fisher's exact test, BH)", q)) +
    omicone_theme() +
    ggplot2::theme(legend.position = "bottom")
}

# ---- survival ---------------------------------------------------------------

#' Normalise sample / patient ids for joining
#'
#' Upper-case and trimmed, so "tcga-ab-2802 " joins "TCGA-AB-2802"; optionally
#' cut to the 12-character TCGA patient barcode (TCGA-XX-XXXX).
#' @param x Ids. @param tcga12 Truncate TCGA barcodes to 12 characters.
#' @keywords internal
wes_norm_id <- function(x, tcga12 = FALSE) {
  if (is.null(x)) return(character(0))
  out <- toupper(trimws(as.character(x)))
  if (isTRUE(tcga12)) {
    tc <- !is.na(out) & grepl("^TCGA-", out)
    out[tc] <- substr(out[tc], 1, 12)
  }
  out
}

#' Per-sample mutation status of a gene set
#'
#' The bridge from a MAF to the shared survival layer: which samples carry a
#' non-synonymous mutation in any of `genes`. Every sample in the MAF is
#' returned, and so is every id in `universe` (samples known to be sequenced
#' but absent from the MAF, e.g. clinical rows that read.maf dropped because
#' they had no variant at all): they are "WT", never missing.
#'
#' @param maf A MAF object. @param genes Character vector of gene symbols.
#' @param universe Optional extra sample ids to return as WT when absent.
#' @param tcga12 Join on 12-character TCGA patient barcodes.
#' @return data.frame: `.id` (normalised), `in_maf`, `mutated`, `status`.
#' @keywords internal
wes_mutation_status <- function(maf, genes, universe = NULL, tcga12 = FALSE) {
  all_s <- wes_samples(maf)
  if (!length(all_s)) stop("No samples found in the MAF.")
  dat <- as.data.frame(maftools::subsetMaf(maf = maf, genes = genes, includeSyn = FALSE,
                                           mafObj = FALSE))
  hit <- if (nrow(dat)) unique(as.character(dat$Tumor_Sample_Barcode)) else character(0)
  maf_ids <- unique(wes_norm_id(all_s, tcga12))
  hit_ids <- unique(wes_norm_id(hit, tcga12))
  extra <- wes_norm_id(universe, tcga12)
  ids <- unique(c(maf_ids, extra[!is.na(extra) & nzchar(extra)]))
  mutated <- ids %in% hit_ids
  data.frame(.id = ids, in_maf = ids %in% maf_ids, mutated = mutated,
             status = ifelse(mutated, "Mutant", "WT"), stringsAsFactors = FALSE)
}

#' Build the per-patient survival analysis set for a mutation split
#'
#' Order of operations (AGENTS.md §7.1-2): normalise ids -> attach mutation
#' status -> collapse to one row per patient (mutant if any sample is) ->
#' optionally keep clinical rows without a MAF record as sequenced WT -> drop
#' missing / non-finite follow-up. Each count along the way is returned in
#' attribute "flow" for the insight bar.
#'
#' @param clin Clinical data.frame (raw, one row per sample or patient).
#' @param status [wes_mutation_status()] result.
#' @param id_col,time_col,event_col Columns of `clin`.
#' @param time_unit Unit of `time_col` ("days", "months", "years").
#' @param patient_col Optional column grouping several samples per patient.
#' @param tcga12 Join on 12-character TCGA barcodes.
#' @param unmatched_wt Keep clinical rows with no MAF record as WT.
#' @return [normalise_clinical()] output with `.group` (WT/Mutant), `.in_maf`.
#' @keywords internal
wes_surv_data <- function(clin, status, id_col, time_col, event_col, time_unit = "days",
                          patient_col = NULL, tcga12 = FALSE, unmatched_wt = TRUE) {
  clin <- as.data.frame(clin, stringsAsFactors = FALSE)
  for (cl in c(id_col, time_col, event_col, patient_col)) {
    if (!cl %in% names(clin)) stop("Column not found in the clinical table: ", cl)
  }
  n_rows <- nrow(clin)
  sid <- wes_norm_id(clin[[id_col]], tcga12)
  has_id <- !is.na(sid) & nzchar(sid)
  clin <- clin[has_id, , drop = FALSE]
  sid <- sid[has_id]
  pid <- if (!is.null(patient_col)) wes_norm_id(clin[[patient_col]]) else sid
  pid[is.na(pid) | !nzchar(pid)] <- sid[is.na(pid) | !nzchar(pid)]
  in_maf <- sid %in% status$.id[status$in_maf]
  mutated <- sid %in% status$.id[status$mutated]

  # several rows per patient are fine only if they agree on follow-up
  key_out <- paste(clin[[time_col]], clin[[event_col]], sep = "\r")
  dup_ids <- unique(pid[duplicated(pid)])
  conflict <- dup_ids[vapply(dup_ids, function(k) length(unique(key_out[pid == k])) > 1,
                             logical(1))]
  if (length(conflict)) {
    stop(sprintf("%d id(s) appear on several rows with different follow-up (e.g. %s). Fix the clinical table, or pick the patient-id column that groups them.",
                 length(conflict), paste(utils::head(conflict, 3), collapse = ", ")))
  }
  first <- !duplicated(pid)
  one <- clin[first, , drop = FALSE]
  one$.pid <- pid[first]
  one$.in_maf <- as.logical(tapply(in_maf, pid, any)[one$.pid])
  one$.mutated <- as.logical(tapply(mutated, pid, any)[one$.pid])
  n_ids <- nrow(one)
  n_matched <- sum(one$.in_maf)
  if (!isTRUE(unmatched_wt)) one <- one[one$.in_maf, , drop = FALSE]

  # Inf / non-numeric follow-up is missing, not "survived forever"
  tt <- suppressWarnings(as.numeric(one[[time_col]]))
  tt[!is.finite(tt)] <- NA_real_
  one[[time_col]] <- tt
  norm <- normalise_clinical(one, ".pid", time_col, event_col, time_unit = time_unit)
  norm$.group <- factor(ifelse(norm$.mutated, "Mutant", "WT"), levels = c("WT", "Mutant"))
  maf_ids <- status$.id[status$in_maf]
  attr(norm, "flow") <- list(
    clin_rows = n_rows, no_id = sum(!has_id), dup_rows = sum(has_id) - n_ids,
    clin_ids = n_ids, maf_samples = length(maf_ids), matched = n_matched,
    unmatched = n_ids - n_matched, unmatched_wt = isTRUE(unmatched_wt),
    maf_only = sum(!maf_ids %in% unique(c(sid, pid))),
    dropped_na = nrow(one) - nrow(norm), n = nrow(norm), events = sum(norm$.event),
    n_mut = sum(norm$.group == "Mutant"), n_wt = sum(norm$.group == "WT"))
  norm
}

#' R code that encodes an event column the way [encode_event()] did
#' @param raw The raw event column. @param col Its name.
#' @keywords internal
wes_event_code <- function(raw, col) {
  q <- function(x) deparse(x)
  if (is.logical(raw)) {
    return(sprintf("clin$.event <- as.integer(clin[[%s]])", q(col)))
  }
  if (is.numeric(raw)) {
    u <- sort(unique(stats::na.omit(raw)))
    if (length(u) == 2 && all(u == c(1, 2))) {
      return(sprintf("clin$.event <- as.integer(clin[[%s]] == 2)", q(col)))
    }
    return(sprintf("clin$.event <- as.integer(clin[[%s]] > 0)", q(col)))
  }
  c(sprintf("ev <- tolower(trimws(as.character(clin[[%s]])))", q(col)),
    paste0('clin$.event <- ifelse(ev %in% c("1", "true", "yes", "y", "dead", "deceased", ',
           '"death", "event", "progressed", "progression", "recurrence", "relapse"), 1L, ',
           'ifelse(ev %in% c("0", "false", "no", "n", "alive", "living", "censored", ',
           '"censor", "no event", "disease free", "disease-free"), 0L, NA))'))
}

#' Runnable R code for the mutation-vs-survival analysis that just ran
#' @param p list(genes, source, id_col, time_col, event_col, time_unit,
#'   patient_col, tcga12, unmatched_wt, raw_event).
#' @keywords internal
wes_surv_code <- function(p) {
  q <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")
  scale <- switch(p$time_unit, days = " / 30.4375", years = " * 12", "")
  norm_fn <- if (isTRUE(p$tcga12)) {
    paste0('norm_id <- function(x) { x <- toupper(trimws(as.character(x))); ',
           'tc <- grepl("^TCGA-", x); x[tc] <- substr(x[tc], 1, 12); x }')
  } else {
    "norm_id <- function(x) toupper(trimws(as.character(x)))"
  }
  src <- if (identical(p$source, "shared")) {
    c('clin_shared <- read.csv("PATH/TO/clinical_cohort.csv", check.names = FALSE)  # <-- edit: the cohort used in Clinical & survival',
      "clin <- clin_shared")
  } else {
    "clin <- clin_raw   # the clinical table read in the WES import step"
  }
  pid_line <- if (!is.null(p$patient_col)) {
    sprintf("clin$.pid <- norm_id(clin[[%s]])", q(p$patient_col))
  } else {
    "clin$.pid <- clin$.sid"
  }
  event <- if (identical(p$source, "shared")) character(0) else wes_event_code(p$raw_event, p$event_col)
  time_line <- if (identical(p$source, "shared")) {
    "clin$.time[!is.finite(clin$.time)] <- NA   # .time is already in months"
  } else {
    sprintf("clin$.time <- suppressWarnings(as.numeric(clin[[%s]]))%s   # months",
            q(p$time_col), scale)
  }
  c(sprintf("genes <- %s", q(p$genes)),
    norm_fn,
    paste0("mut_tsb <- unique(as.character(maftools::subsetMaf(maf = maf, genes = genes, ",
           "includeSyn = FALSE, mafObj = FALSE)$Tumor_Sample_Barcode))"),
    "maf_ids <- norm_id(maftools::getSampleSummary(maf)$Tumor_Sample_Barcode)",
    src,
    sprintf("clin$.sid <- norm_id(clin[[%s]])", q(p$id_col)),
    "clin <- clin[!is.na(clin$.sid) & nzchar(clin$.sid), ]",
    pid_line,
    "clin$.mut <- ave(clin$.sid %in% norm_id(mut_tsb), clin$.pid, FUN = any)",
    "clin$.in_maf <- ave(clin$.sid %in% maf_ids, clin$.pid, FUN = any)",
    "clin <- clin[!duplicated(clin$.pid), ]   # one row per patient",
    if (!isTRUE(p$unmatched_wt)) "clin <- clin[clin$.in_maf, ]   # no MAF record = not analysed"
    else "# clinical rows without a MAF record stay in, as sequenced wild-type",
    time_line,
    event,
    "clin <- clin[!is.na(clin$.time) & clin$.time >= 0 & !is.na(clin$.event), ]",
    'clin$.group <- factor(ifelse(clin$.mut, "Mutant", "WT"), levels = c("WT", "Mutant"))',
    "table(clin$.group, clin$.event)",
    "fit <- survival::survfit(survival::Surv(.time, .event) ~ .group, data = clin)",
    "summary(fit)$table   # medians (NA = not reached)",
    "survival::survdiff(survival::Surv(.time, .event) ~ .group, data = clin)",
    "summary(survival::coxph(survival::Surv(.time, .event) ~ .group, data = clin))")
}

# ---- heterogeneity ----------------------------------------------------------

#' Infer clonal structure of one sample from its VAF distribution
#'
#' `vaf_col = NULL` lets maftools use `t_vaf` or compute it from
#' `t_ref_count` / `t_alt_count`.
#' @param maf A MAF object. @param sample Tumor_Sample_Barcode.
#' @param vaf_col VAF column name, or NULL.
#' @keywords internal
wes_heterogeneity <- function(maf, sample, vaf_col = NULL) {
  args <- list(maf = maf, tsb = sample)
  if (!is.null(vaf_col)) args$vafCol <- vaf_col
  het <- do.call(maftools::inferHeterogeneity, args)
  if (is.null(het) || is.null(het$clusterData) || !nrow(het$clusterData)) {
    stop("No clusters for ", sample, ": fewer than 3 variants with a usable VAF ",
         "(mclust needs at least 3), or none outside copy-number-altered regions.")
  }
  het
}

#' Number of VAF clusters, not counting outliers and copy-number-altered variants
#' @param cluster The `cluster` column of `inferHeterogeneity()$clusterData`.
#' @keywords internal
wes_clone_count <- function(cluster) {
  cl <- as.character(cluster)
  length(unique(cl[!is.na(cl) & !cl %in% c("outlier", "CN_altered")]))
}

#' MATH scores of the whole cohort (maftools::math.score, VAF >= 0.075,
#' samples with >= 5 such variants)
#' @param maf A MAF object. @param vaf_col VAF column or NULL.
#' @keywords internal
wes_math_scores <- function(maf, vaf_col = NULL) {
  args <- list(maf = maf)
  if (!is.null(vaf_col)) args$vafCol <- vaf_col
  as.data.frame(suppressMessages(do.call(maftools::math.score, args)))
}

#' Where one sample's MATH score falls in the cohort's distribution
#' @param scores [wes_math_scores()] output. @param sample Sample id.
#' @return list(math, percentile, n) — NA when the sample has no score.
#' @keywords internal
wes_math_percentile <- function(scores, sample) {
  v <- scores$MATH[!is.na(scores$MATH)]
  x <- scores$MATH[as.character(scores$Tumor_Sample_Barcode) == sample][1]
  if (!length(v) || is.null(x) || is.na(x)) {
    return(list(math = NA_real_, percentile = NA_real_, n = length(v)))
  }
  list(math = x, percentile = 100 * mean(v <= x), n = length(v))
}

# ---- exports ----------------------------------------------------------------

#' Top-bar export menu buttons for the WES pipeline
#' @keywords internal
wes_export_items <- function() {
  shiny::tagList(
    export_button("dl_maf_rds",     "MAF object (.rds)",     "MAF 对象 (.rds)"),
    export_button("dl_maf_tsv",     "Mutation table (.tsv)", "突变表 (.tsv)"),
    export_button("dl_maf_samples", "Sample summary (.csv)", "样本汇总 (.csv)"),
    export_button("dl_maf_genes",   "Gene summary (.csv)",   "基因汇总 (.csv)"))
}

#' Download handlers behind [wes_export_items()]
#' @param output Server output. @param rv Shared hub (uses `rv$maf`).
#' @keywords internal
register_wes_exports <- function(output, rv) {
  stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
  need_maf <- function() {
    if (is.null(rv$maf)) {
      wes_notify("Import a MAF first, then export.", "请先导入 MAF，再导出。",
                 type = "warning", duration = 6)
      return(FALSE)
    }
    TRUE
  }

  output$dl_maf_rds <- shiny::downloadHandler(
    filename = function() paste0("omicone_maf_", stamp(), ".rds"),
    content = function(file) {
      if (!need_maf()) return(saveRDS(NULL, file))
      saveRDS(rv$maf, file)
    }
  )

  # Non-synonymous and silent variants together, as one flat table
  output$dl_maf_tsv <- shiny::downloadHandler(
    filename = function() paste0("omicone_mutations_", stamp(), ".tsv"),
    content = function(file) {
      if (!need_maf()) return(utils::write.table(data.frame(), file))
      m <- rv$maf
      parts <- Filter(function(x) !is.null(x) && nrow(x),
                      list(as.data.frame(m@data), as.data.frame(m@maf.silent)))
      cols <- unique(unlist(lapply(parts, names)))
      parts <- lapply(parts, function(x) {
        x[setdiff(cols, names(x))] <- NA
        x[cols]
      })
      utils::write.table(do.call(rbind, parts), file, sep = "\t",
                         quote = FALSE, row.names = FALSE, na = "")
    }
  )

  output$dl_maf_samples <- shiny::downloadHandler(
    filename = function() paste0("omicone_sample_summary_", stamp(), ".csv"),
    content = function(file) {
      if (!need_maf()) return(utils::write.csv(data.frame(), file))
      utils::write.csv(as.data.frame(maftools::getSampleSummary(rv$maf)),
                       file, row.names = FALSE)
    }
  )

  output$dl_maf_genes <- shiny::downloadHandler(
    filename = function() paste0("omicone_gene_summary_", stamp(), ".csv"),
    content = function(file) {
      if (!need_maf()) return(utils::write.csv(data.frame(), file))
      utils::write.csv(as.data.frame(maftools::getGeneSummary(rv$maf)),
                       file, row.names = FALSE)
    }
  )
}
