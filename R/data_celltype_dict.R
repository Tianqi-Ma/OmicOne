#' Plain-language cell-type explanation dictionary
#'
#' A lookup used by the annotation module's cell-type table: given a predicted
#' cell-type label, show a one-line, beginner-friendly explanation of what that
#' cell type is. Matching is case-insensitive, on whole words, and tolerant of
#' the label styles of celldex (SingleR), Azimuth and CellTypist.
#'
#' This is intentionally a general starter set covering frequent immune,
#' stromal, epithelial, and developmental types. Extend `celltype_dictionary()`
#' as needed for your tissue of interest.
#'
#' @name data_celltype_dict
#' @keywords internal
NULL

#' The cell-type explanation table
#'
#' @return A data.frame with columns `key` (lowercase canonical token, words
#'   separated by single spaces) and `explanation` (one-line description).
#' @keywords internal
celltype_dictionary <- function() {
  d <- c(
    "t cell"            = "Immune cell that coordinates and executes adaptive immune responses.",
    "cd4 t cell"        = "Helper T cell; directs other immune cells (a coordinator).",
    "cd4 t"             = "Helper T cell; directs other immune cells (a coordinator).",
    "helper t cell"     = "Helper T cell; directs other immune cells (a coordinator).",
    "cd4 naive"         = "Naive helper T cell that has not yet met its antigen.",
    "cd4 tcm"           = "Central-memory helper T cell; long-lived, recirculates through lymph nodes.",
    "cd4 tem"           = "Effector-memory helper T cell; responds quickly in tissues.",
    "cd4 ctl"           = "Cytotoxic CD4 T cell; a helper T cell that can kill target cells.",
    "cd8 t cell"        = "Cytotoxic T cell; kills infected or abnormal cells.",
    "cd8 t"             = "Cytotoxic T cell; kills infected or abnormal cells.",
    "cytotoxic t cell"  = "Cytotoxic T cell; kills infected or abnormal cells.",
    "cd8 naive"         = "Naive cytotoxic T cell that has not yet met its antigen.",
    "cd8 tcm"           = "Central-memory cytotoxic T cell; long-lived, recirculates through lymph nodes.",
    "cd8 tem"           = "Effector-memory cytotoxic T cell; kills quickly on re-encounter.",
    "regulatory t cell" = "Treg; dampens immune responses to prevent over-reaction.",
    "treg"              = "Treg; dampens immune responses to prevent over-reaction.",
    "mait"              = "MAIT cell; innate-like T cell that senses bacterial metabolites.",
    "gdt"               = "Gamma-delta T cell; innate-like T cell common in barrier tissues.",
    "gamma delta t cell"= "Gamma-delta T cell; innate-like T cell common in barrier tissues.",
    "dnt"               = "Double-negative T cell (neither CD4 nor CD8).",
    "b cell"            = "Immune cell that produces antibodies.",
    "b naive"           = "Naive B cell that has not yet met its antigen.",
    "b memory"          = "Memory B cell; responds fast when the same antigen returns.",
    "b intermediate"    = "B cell between the naive and memory states.",
    "plasma cell"       = "Mature B cell specialised for mass antibody production.",
    "plasma"            = "Mature B cell specialised for mass antibody production.",
    "plasmablast"       = "Early antibody-secreting B cell, precursor of plasma cells.",
    "nk cell"           = "Natural killer cell; kills stressed/infected cells without prior priming.",
    "nk"                = "Natural killer cell; kills stressed/infected cells without prior priming.",
    "ilc"               = "Innate lymphoid cell; a lymphocyte without antigen receptors.",
    "monocyte"          = "Circulating immune cell that becomes a macrophage in tissue.",
    "mono"              = "Circulating immune cell that becomes a macrophage in tissue.",
    "cd14 mono"         = "Classical (CD14+) monocyte; the main circulating monocyte.",
    "cd16 mono"         = "Non-classical (CD16+) monocyte; patrols blood vessel walls.",
    "macrophage"        = "Tissue immune cell that engulfs debris and pathogens.",
    "dendritic cell"    = "Antigen-presenting cell that activates T cells.",
    "dc"                = "Dendritic cell; antigen-presenting cell that activates T cells.",
    "cdc1"              = "Conventional dendritic cell type 1; primes cytotoxic T cells.",
    "cdc2"              = "Conventional dendritic cell type 2; primes helper T cells.",
    "dc2"               = "Conventional dendritic cell type 2; primes helper T cells.",
    "pdc"               = "Plasmacytoid dendritic cell; produces type I interferon against viruses.",
    "asdc"              = "AXL+ SIGLEC6+ dendritic cell, a rare blood dendritic-cell subset.",
    "neutrophil"        = "Fast-responding immune cell against bacterial infection.",
    "eosinophil"        = "Granulocyte involved in parasite defence and allergy.",
    "basophil"          = "Rare granulocyte releasing histamine in allergic responses.",
    "mast cell"         = "Immune cell releasing histamine; involved in allergy.",
    "hspc"              = "Haematopoietic stem/progenitor cell; gives rise to all blood cells.",
    "hsc"               = "Haematopoietic stem cell; gives rise to all blood cells.",
    "erythrocyte"       = "Red blood cell; carries oxygen (often a contaminant in scRNA-seq).",
    "eryth"             = "Red blood cell lineage; carries oxygen (often a contaminant in scRNA-seq).",
    "erythroblast"      = "Red-blood-cell precursor.",
    "platelet"          = "Cell fragment involved in blood clotting.",
    "doublet"           = "Two cells captured together; not a real cell type.",
    "epithelial cell"   = "Cell forming the lining of organs, glands and surfaces.",
    "epithelial"        = "Cell forming the lining of organs, glands and surfaces.",
    "keratinocyte"      = "Main cell of the skin's outer layer (epidermis).",
    "endothelial cell"  = "Cell lining blood and lymphatic vessels.",
    "endothelial"       = "Cell lining blood and lymphatic vessels.",
    "fibroblast"        = "Structural cell producing extracellular matrix (connective tissue).",
    "smooth muscle cell"= "Involuntary muscle cell in vessels and organs.",
    "smooth muscle"     = "Involuntary muscle cell in vessels and organs.",
    "pericyte"          = "Cell wrapping capillaries to regulate blood flow.",
    "adipocyte"         = "Fat cell; stores energy as lipid.",
    "chondrocyte"       = "Cartilage cell.",
    "osteoblast"        = "Bone-forming cell.",
    "melanocyte"        = "Pigment-producing cell of skin and eye.",
    "myocyte"           = "Muscle cell.",
    "cardiomyocyte"     = "Heart muscle cell responsible for contraction.",
    "neuron"            = "Nerve cell that transmits electrical/chemical signals.",
    "astrocyte"         = "Support (glial) cell in the brain maintaining neurons.",
    "oligodendrocyte"   = "Glial cell that myelinates neurons in the CNS.",
    "microglia"         = "Resident immune cell of the brain.",
    "hepatocyte"        = "Main functional cell of the liver.",
    "stem cell"         = "Undifferentiated cell that can self-renew and specialise.",
    "progenitor cell"   = "Partially committed cell on its way to a mature type.",
    "progenitor"        = "Partially committed cell on its way to a mature type.",
    "proliferating cell"= "Actively dividing cell (high cell-cycle gene expression).",
    "proliferating"     = "Actively dividing cell (high cell-cycle gene expression)."
  )
  data.frame(
    key = names(d),
    explanation = unname(d),
    stringsAsFactors = FALSE
  )
}

#' Normalise a cell-type label for dictionary lookup (pure)
#'
#' Lower case; `+`, `_`, `-`, `/` and brackets become spaces; plural words lose
#' their final "s" ("cells" -> "cell", but "class" is kept).
#' @param label Character vector.
#' @keywords internal
normalize_celltype_label <- function(label) {
  x <- tolower(trimws(as.character(label)))
  x <- gsub("[+]", "", x)
  x <- gsub("[_/()-]+", " ", x)
  x <- gsub("\\b([a-z]{2,}[^s ])s\\b", "\\1", x, perl = TRUE)
  x <- gsub("\\s+", " ", trimws(x))
  x
}

#' Look up a plain-language explanation for a cell-type label
#'
#' Exact match first; otherwise the longest dictionary key found as whole
#' words, so "Mast cells" is a mast cell (not a T cell, as a substring match
#' on "t cell" would say) and "CD14 Mono" beats "Mono".
#' @param label Character vector of predicted cell-type labels.
#' @return Character vector of explanations (empty string if not found).
#' @keywords internal
explain_celltype <- function(label) {
  dict <- celltype_dictionary()
  norm <- normalize_celltype_label(label)
  pat <- paste0("(^| )", dict$key, "( |$)")
  out <- vapply(norm, function(x) {
    if (is.na(x) || !nzchar(x)) return("")
    hit <- which(dict$key == x)
    if (length(hit)) return(dict$explanation[hit[1]])
    part <- which(vapply(pat, function(p) grepl(p, x), logical(1)))
    if (!length(part)) return("")
    best <- part[which.max(nchar(dict$key[part]))]
    dict$explanation[best]
  }, character(1))
  unname(out)
}
