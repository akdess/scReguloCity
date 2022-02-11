
#' @title CreateSCReguloCityObject
#'
#' @param raw.data  the matrix of genes (rows) vs. cells (columns) containing the raw counts
#' 
#' @description Creation of a casper object.
#'
#' @return casper
#'
#' @export
#'
CreateSCReguloCityObject <- function(emat, nmat, smat, emb, species, motif_ref, cell.types, 
    ...) {
    
    object <- new(Class = "scReguloCity", emat = emat, nmat = nmat, smat = smat,
     emb=emb, species=species, motif_ref=motif_ref, cell.types=cell.types)
   
    common <- intersect(rownames(object@emat), rownames(object@nmat))
   
    if(!is.null(object@smat))  {
      common <- intersect(common, rownames(object@smat))
      object@smat <- object@smat[match(common, rownames(object@smat)),]
    }
    object@emat <- object@emat[match(common, rownames(object@emat)),]
    object@nmat <- object@nmat[match(common, rownames(object@nmat)),]
    

    if(is.null(object@smat)) { object <- calculateVelocity(object) }
    
    if(!is.null(object@smat))  {object <- calculateVelocityWithOffset(object) }

 #  object <- correctGCContentAndMappability(object)
    
    return(object)
}


#' @title calculateVelocity
#'
#' @param raw.data  the matrix of genes (rows) vs. cells (columns) containing the raw counts
#' 
#' @description Creation of a casper object.
#'
#' @return casper
#'
#' @export
#'
calculateVelocity <- function(object)
{
    emat <- object@emat
    nmat <- object@nmat
    smat <- object@smat
    vel.mat <- NULL
    vel.mat <- data.frame(do.call(rbind,lapply(1:length(rownames(emat)),function(y) {
          df <- data.frame(n=(nmat[y,]),e=(emat[y,]))
          df$o <- 0
          d <- lm(n~e+offset(o), data=df)
          #d <- lm(n~e, data=df)
          return( resid(d))
        })))
    rownames(vel.mat) <- rownames(emat)
    colnames(vel.mat) <- colnames(emat)

    
    common <- intersect(colnames(vel.mat), rownames(object@emb))
    object@emb <- object@emb[match(common, rownames(object@emb)),]
    object@vel <- vel.mat[, match(common, colnames(vel.mat))]
    object@cell.types <- object@cell.types[match(common, rownames(object@emb))]
    return(object)
}

calculateVelocityWithOffset <- function(object)
{
 

    emat <- object@emat
    nmat <- object@nmat
    smat <- object@smat
    vel.mat <- NULL
      
    sfit <- data.frame(do.call(rbind,lapply(1:length(rownames(emat)),function(y) {
      df <- data.frame(n=(nmat[y,]),e=(emat[y,]),s=smat[y,])
      sd <- lm(n~s,data=df)
      r <- with(df[df$s>0,],cor(n,s,method='spearman'),3)
      return(c(o=pmax(0,as.numeric(sd$coef[1])),s=as.numeric(sd$coef[2]),r=r))
    })))
    rownames(sfit) <- rownames(emat)

    
    vel.mat <- data.frame(do.call(rbind,lapply(1:length(rownames(emat)),function(gn) {
        df <- data.frame(n=(nmat[gn,]),e=(emat[gn,]),o=sfit[gn,'o'])
        d <- lm(n~e+offset(o), data=df)
        return(resid(d))
    })))

    rownames(vel.mat) <- rownames(emat)
    colnames(vel.mat) <- colnames(emat)
    
    common <- intersect(colnames(vel.mat), rownames(object@emb))
    object@emb <- object@emb[match(common, rownames(object@emb)),]
    object@vel <- vel.mat[, match(common, colnames(vel.mat))]
    object@cell.types <- object@cell.types[match(common, rownames(object@emb))]
    return(object)

}

