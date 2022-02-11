#' @title readIntrexStats()
#'
#' @description Read the output of IntrexStats algorithm
#'
#' @param file.dir file directory
#' 
#' @param file.name file name
#' 
#' @return accuracy measures 
#'
#' @export
#'
#' 
readIntrexStats <- function(file.dir, file.name, emb) {
 
    intrex <- read.delim(file.path(file.dir, file.name), stringsAsFactor=F)
   	common_colnames <- intersect(rownames(emb), colnames(intrex)[-(1:2)])
		emb <- emb [match(common_colnames, rownames(emb)), ]

    emat <- intrex[intrex$EXON_INTRON_IDENTIFIER=="EXONLY", -(1:2) ]
    rownames(emat) <- intrex[intrex$EXON_INTRON_IDENTIFIER=="EXONLY", "GENE_NAME"]
   	
    nmat <- intrex[intrex$EXON_INTRON_IDENTIFIER=="INTRONLY", -(1:2) ]
    rownames(nmat) <- intrex[intrex$EXON_INTRON_IDENTIFIER=="INTRONLY", "GENE_NAME"]
   
    smat <- intrex[intrex$EXON_INTRON_IDENTIFIER=="INTREXIC", -(1:2) ]
    rownames(smat) <- intrex[intrex$EXON_INTRON_IDENTIFIER=="INTREXIC", "GENE_NAME"]
   	
		emat <- emat [, match(common_colnames, colnames(emat)) ]
		nmat <- nmat [ ,match(common_colnames, colnames(nmat)) ]
		smat <- smat [ ,match(common_colnames, colnames(smat)) ]

    return(list(emat=emat, nmat=nmat, smat=smat, emb=emb))
}


#' @title normalizeAndSmoothIntrexStats()
#'
#' @description Normalize and Smooth Intrex statistics
#'
#' @param emat spliced read count 
#' 
#' @param nmat ambigious read count 
#' 
#' @param smat unspliced read count 
#' 
#' @param k number of 
#' 
#' @return normalized and smoothed values
#'
#' @export
#'
#' 
normalizeAndSmoothIntrexStats <- function(emat, nmat, smat, k=15, emb=emb )
{
	emat.size <- Matrix::colSums(emat); 
	nmat.size <- Matrix::colSums(nmat); 
	
	emat.cs <- emat.size[colnames(emat)]/1e3;
	nmat.cs <- nmat.size[colnames(nmat)]/1e3;
	
	conv.emat.norm <- log(t(t(as.matrix(emat))/emat.cs)+1);
	conv.nmat.norm <- log(t(t(as.matrix(nmat))/nmat.cs)+1);
	conv.smat.norm <- NULL

	if(!is.null(smat)) 
	{
		smat.size <- Matrix::colSums(smat); 
		smat.cs <- smat.size[colnames(smat)]/1e3;
		conv.smat.norm <- log(t(t(as.matrix(smat))/smat.cs)+1);
		colnames(conv.smat.norm) <- colnames(smat)
	}
	
	conv.nmat.norm <- knn_smoothing(mat=conv.nmat.norm, k=k)
	conv.emat.norm <- knn_smoothing(mat=conv.emat.norm, k=k)
	

	if(!is.null(conv.smat.norm)) {
		  conv.smat.norm <- knn_smoothing(mat=conv.smat.norm, k=15)
			stats <- list(emat=conv.emat.norm, nmat=conv.nmat.norm, smat=conv.smat.norm)
	} else {
			stats <- list(emat=conv.emat.norm, nmat=conv.nmat.norm)
	}


	common_rownames <- Reduce(intersect, lapply(stats, rownames))
	stats  <- lapply(stats, function(x) { x[rownames(x) %in% common_rownames,] })

	common_colnames <- Reduce(intersect, lapply(stats, colnames))
	stats  <- lapply(stats, function(x) { x[, colnames(x) %in% common_colnames] })

	return(stats)

}


### code from https://github.com/yanailab/knn-smoothing/blob/master/knn_smooth.R
# K-nearest neighbor smoothing for high-throughput single-cell RNA-Seq data
# (R implementation.)

# Authors:
#   Yun Yan <yun.yan@nyumc.org>
#   Florian Wagner <florian.wagner@nyu.edu>
# Copyright (c) 2017, 2018 New York University

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(rsvd))

randomized_pca <- function(tmat, d, seed){
  # @param tmat A non-negative matrix with samples by features
  # @return A matrix with features by samples projected on PCA space
  set.seed(seed)
  #rsvd_obj <- rsvd(scale(tmat, center = TRUE, scale = FALSE), k=d)
  #rsvd_obj$u %*% diag(rsvd_obj$d)
  rpca_obj <- rpca(tmat, k=d, center=T, scale=F, retx=T, p=10, q=7)
  rpca_obj$x
}

normalization_median <- function(mat){
  # Median normalization
  # @param mat A non-negative matrix with genes by samples
  num_transcripts <- Matrix::colSums(mat)
  size_factor <- median(num_transcripts, na.rm = T) / num_transcripts
  t(t(mat) * size_factor)
}

freeman_tukey_transform <- function(mat){
  sqrt(mat) + sqrt(mat + 1)
}

pdist <- function(tmat){
  # @param tmat A non-negative matrix with samples by features
  # @reference http://r.789695.n4.nabble.com/dist-function-in-R-is-very-slow-td4738317.html
  mtm <- Matrix::tcrossprod(tmat)
  sq <- rowSums(tmat^2)
  out0 <- outer(sq, sq, "+") - 2 * mtm
  out0[out0 < 0] <- 0

  sqrt(out0)
}

smoother_aggregate_nearest_nb <- function(mat, D, k){
  # @param mat A matrix in a shape of #genes x #samples.
  # @param D A predefined distance matrix in a shape of #samples x #samples.
  # @param k An integer to choose \code{k} nearest samples (self-inclusive) to
  #  aggregate based on the distance matrix \code{D}. If \code{k} is greater than
  #  #samples, \code{k} is forced to be #samples to continue aggregation.
  sapply(seq_len(ncol(mat)), function(cid){
    nb_cid <- head(order(D[cid, ]), k)
    closest_mat <- mat[, nb_cid, drop=FALSE]
    return(Matrix::rowSums(closest_mat))
  })
}

#' 
knn_smoothing <- function(mat, k, d=10, seed=42){
  #' KNN-smoothing on UMI-filtered single-cell RNA-seq data
  #'
  #' @param mat A numeric matrix with gene names on rows and cell names on columns.
  #' @param k Number of nearest neighbours to aggregate.
  #' @param d Number of Principal components.
  #' @param seed Seed number. (default=42)
  #' @return A smoothed numeric matrix.
  #' @examples
  #' X <- matrix(abs(sin(seq(from=1, to=1000, length.out = 1000))),
  #' nrow = 25, byrow = T)
  #' y <- rep(1:4, each=10)
  #' dim(X)
  #' colnames(X) <- as.character(paste0("s", seq_len(ncol(X))))
  #' rownames(X) <- as.character(paste0("g", seq_len(nrow(X))))
  #' S <- knn_smoother(X, k=5)
  #' plot(X[1, ], X[3, ], col=factor(y), main="original")
  #' plot(S[1, ], S[3, ], col=factor(y), main="smoothed")
  #' @export

  cname <- colnames(mat)
  gname <- rownames(mat)

  num_steps <- ceiling(log2(k + 1))
  S <- mat
  for (p in seq(1, num_steps)){
    k_step <- min(2^p - 1, k)
    message(paste0('Step ', p, '/', num_steps, ': ',
                   'Smoothing using k=', k_step))
    Y <- freeman_tukey_transform(normalization_median(S))
    if (! is.null(d)) {
      Y <- t(randomized_pca(t(Y), d=d, seed=seed))
    }
    D <- pdist(t(Y))
    S <- smoother_aggregate_nearest_nb(mat, D, k_step + 1)
  }
  if (! is.null(cname)) colnames(S) <- cname
  if (! is.null(gname)) rownames(S) <- gname

  S
}

