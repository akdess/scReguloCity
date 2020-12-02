


#' @details
#' The main functions you will need to use are CreateSCReguloCityObject() and runscReguloCity(scReguloCity_object).
#' For additional details on running the analysis step by step, please refer to the example vignette.
#' @aliases scReguloCity-package
"_PACKAGE"


#' The scReguloCity Class
#'
#'
#' The scReguloCity Class
#' The scReguloCity object is required for performing CNV analysis on single-cell and bulk RNA-Seq. It stores all information
#' associated with the dataset, including data, smoothed data, baf values, annotations, scale specific segments, scale specific large scale events etc.
#'
#' @name scReguloCity
#' @rdname scReguloCity
#' @aliases scReguloCity-class
#' @exportClass scReguloCity

scReguloCity <- methods::setClass("scReguloCity", slots = c(emat = "ANY", 
	smat = "ANY", nmat = "ANY", emb = "ANY", mappability= "ANY", species ="ANY", gradient="ANY", 
	gradGenes= "ANY", gradGenesMasked= "ANY", filtered="ANY",
	motif_ref="ANY", vel="ANY", vel.corrected="ANY", tfCorNet="ANY", tfCor="ANY"))

setMethod(f = "show", signature = "scReguloCity", definition = function(object) {
       invisible(x = NULL)
})

