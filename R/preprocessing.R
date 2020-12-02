
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
CreateSCReguloCityObject <- function(emat, nmat, smat, emb, species, motif_ref, mappability,
    ...) {
    
    object <- new(Class = "scReguloCity", emat = emat, nmat = nmat, smat = smat,
     emb=emb, species=species, motif_ref=motif_ref, mappability=mappability)
   
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
    return(object)

}


#    mappabilty <- read.table("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\mappability\\mm10_gencode.vM21_intrex_mmap_seq_stats.txt", header=T)


#' @title ProcessData()
#'
#' @description Processing expression signal. Step 1. Recursively iterative median filtering  Step 2. Center Normalization Step 3. Control Normalization
#'
#' @param object casper object
#' 
#' @return object
#'
#' @export
#'
#' 
correctGCContentAndMappability <- function(object) {

    map <- object@mappability
 
    map$exonMapp <- map$MM_EXON/map$EXON_LENGTH 
    map$intronMapp <- map$MM_INTRON/map$INTRON_LENGTH 

    exon_all <- apply(map[, c("EXON_A", "EXON_C", "EXON_G", "EXON_T")], 1, sum)
    exon_gc <- apply(map[, c( "EXON_C", "EXON_G")], 1, sum)

    map$exonGC <- exon_gc/exon_all

    intron_all <- apply(map[, c("INTRON_A", "INTRON_C", "INTRON_G", "INTRON_T")], 1, sum)
    intron_gc <- apply(map[, c( "INTRON_C", "INTRON_G")], 1, sum)

    map$intronGC <- intron_gc/intron_all
    map$gc <- map$intronGC/map$exonGC
    map$d <- map$intronMapp/map$exonMapp
    map <- na.omit(map)
    common <- intersect(rownames(object@vel), map$GENE_ID)
    map <- map[match(common, map$GENE_ID), ]
   

    vel <- object@vel[match(common, rownames(object@vel)), ]

    vel.corrected <- data.frame(do.call(cbind, lapply(1:ncol(vel),function(y) {

        df <- data.frame(n=vel[,y],m=map$d, gc=map$gc)
        vel.cor <- lm(df$n~df$gc+df$m) 
        if(is.null(resid(vel.cor))) { print(y) }
        return(resid(vel.cor))
    })))

    rownames(vel.corrected) <- rownames(vel)
    colnames(vel.corrected) <- colnames(vel)
    
    object@vel.corrected <- vel.corrected
    
    return(object)
}