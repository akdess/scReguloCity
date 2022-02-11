


#' @details
#' The main functions you will need to use are CreateSCReguloCityObject() and applyReguloCity(scReguloCity_object).
#' For additional details on running the analysis step by step, please refer to the example vignette.
#' @aliases scReguloCity-package
"_PACKAGE"


#' The scReguloCity Class
#'
#'
#' The scReguloCity Class
#' The scReguloCity object is required for finding local RNA Velocity patterns. It stores all information
#' associated with the dataset.
#'
#' @name scReguloCity
#' @rdname scReguloCity
#' @aliases scReguloCity-class
#' @exportClass scReguloCity

scReguloCity <- methods::setClass("scReguloCity", slots = c(emat = "ANY", 
	smat = "ANY", nmat = "ANY", emb = "ANY",  species ="ANY", gradient="ANY", 
	gradGenes= "ANY", gradGenesMasked= "ANY", filtered="ANY", cell.types="ANY",
	motif_ref="ANY", vel="ANY",   pattern_clusters="ANY",  vel.corrected="ANY", tfCorNet="ANY", tfCor="ANY"))

setMethod(f = "show", signature = "scReguloCity", definition = function(object) {
       invisible(x = NULL)
})

