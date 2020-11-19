
pat <- function(k){
 if(is.na(k)) return ("none");
 if(k<0) return("down");
 if(k>0) return("up")
}

#' @title loadIntrexStats()
#'
#' @description Calculates tpr and fpr values using genotyping array as gold standard
#'
#' @param chrMat large scale event matrix generated using CaSpER
#' 
#' @return accuracy measures 
#'
#' @export
#'
#' 
loadIntrexStats <- function(file.path, file.name) {
 
    intrex <- read.delim(file.path(file.dir, file.name), stringsAsFactor=F)
    emat <- intrex[intrex$EXON_INTRON_IDENTIFIER=="EXONLY", -(1:2) ]
    rownames(emat) <- intrex[intrex$EXON_INTRON_IDENTIFIER=="EXONLY", "GENE_NAME"]
   
    nmat <- intrex[intrex$EXON_INTRON_IDENTIFIER=="INTRONLY", -(1:2) ]
    rownames(nmat) <- intrex[intrex$EXON_INTRON_IDENTIFIER=="INTRONLY", "GENE_NAME"]
   
    smat <- intrex[intrex$EXON_INTRON_IDENTIFIER=="INTREXIC", -(1:2) ]
    rownames(smat) <- intrex[intrex$EXON_INTRON_IDENTIFIER=="INTREXIC", "GENE_NAME"]
   
    return(list(emat=emat, nmat=nmat, smat=smat))
}

filter.genes.by.cluster.expression <- function(emat,clusters,min.max.cluster.average=0.1) {
  if(!any(colnames(emat) %in% names(clusters))) stop("provided clusters do not cover any of the emat cells!")
  vc <- intersect(colnames(emat),names(clusters))
  cl.emax <- apply(do.call(cbind,tapply(vc,as.factor(clusters[vc]),function(ii) Matrix::rowMeans(emat[,ii]))),1,max)
  vi <- cl.emax>min.max.cluster.average;
  emat[vi,]
}

#' adjust colors, while keeping the vector names
#' 
#' @param x color vector
#' @param alpha transparenscy value (passed to adjustcolors as alpha.f)
#' @param ... parameters passsed to adjustcolor
#' @export
ac <- function(x, alpha=1, ...) { y <- adjustcolor(x, alpha.f=alpha, ...); names(y) <- names(x); return(y)}

# quick function to map value vector to colors
val2col <- function(x,gradientPalette=NULL,zlim=NULL,gradient.range.quantile=0.95) {
  if(all(sign(x)>=0)) {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c('gray90','red'), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- as.numeric(quantile(na.omit(x),p=c(1-gradient.range.quantile,gradient.range.quantile)))
      if(diff(zlim)==0) {
        zlim <- as.numeric(range(na.omit(x)))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])
    
  } else {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- c(-1,1)*as.numeric(quantile(na.omit(abs(x)),p=gradient.range.quantile))
      if(diff(zlim)==0) {
        zlim <- c(-1,1)*as.numeric(na.omit(max(abs(x))))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])
    
  }
  
  gp <- gradientPalette[x*(length(gradientPalette)-1)+1]
  if(!is.null(names(x))) { names(gp) <- names(x) }
  gp
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


normalizeAndSmoothIntrexStats <- function(emat, nmat, smat, cell.colors )
{


	emat.size <- Matrix::colSums(emat2); 
	nmat.size <- Matrix::colSums(nmat2); 
	
	emat.cs <- emat.size[colnames(emat2)]/1e3;
	nmat.cs <- nmat.size[colnames(nmat2)]/1e3;
	
	conv.emat.norm <- log(t(t(as.matrix(emat2))/emat.cs)+1);
	conv.nmat.norm <- log(t(t(as.matrix(nmat2))/nmat.cs)+1);
	conv.smat.norm <- NULL

	if(!is.null(smat)) 
	{
		smat.size <- Matrix::colSums(smat2); 
		smat.cs <- smat.size[colnames(smat2)]/1e3;
		conv.smat.norm <- log(t(t(as.matrix(smat2))/smat.cs)+1);
		colnames(conv.smat.norm) <- colnames(smat2)
	}

	colnames(conv.emat.norm) <- colnames(emat2)
	colnames(conv.nmat.norm) <- colnames(nmat2)
	common <-  intersect(rownames(conv.nmat.norm), rownames(conv.emat.norm))
	
	if(!is.null(conv.smat.norm)) {
		common <- intersect(common, rownames(conv.smat.norm))
		conv.smat.norm <- conv.smat.norm[match(common, rownames(conv.smat.norm)), ]
	}
	
	conv.nmat.norm <- conv.nmat.norm[match(common, rownames(conv.nmat.norm)), ]
	conv.emat.norm <- conv.emat.norm[match(common, rownames(conv.emat.norm)), ]	
	
	conv.nmat.norm <- knn_smoothing(mat=conv.nmat.norm, k=15)
	conv.emat.norm <- knn_smoothing(mat=conv.emat.norm, k=15)
	
	if(!is.null(conv.smat.norm)) {conv.smat.norm <- knn_smoothing(mat=conv.smat.norm, k=15)}
	
	return(list(emat=conv.emat.norm, nmat=conv.nmat.norm, smat=conv.smat.norm))

}

plotFunctionalEnrichment<- function(){


	numPatterns <- unique(all$finalPattern)
	genes2 <-unique(as.character(sub.p$geneId))
	universe.id <-as.vector(unique(na.omit(unlist(mget (genes2, org.Mm.egALIAS2EG,ifnotfound=NA)))))
	go.BP.plot<- list()
	go.MF.plot<- list()
	kegg.plot<- list()
	i <- 1
	for (k in 1:length(numPatterns)) {

		genes <-as.character(all[all$finalPattern==numPatterns[k],1])
		entrez.id <-as.vector(unique(na.omit(unlist(mget (genes,org.Mm.egALIAS2EG,ifnotfound=NA)))))


		go.BP<-NULL; go.MF<-NULL; kegg=NULL;


		go.BP<-go.enrichment.BP(genes=entrez.id , ontology="BP",universe=universe.id,  
		pvalue=0.05, annotation='org.Mm.eg.db', conditionalSearch=TRUE, genes2=genes)
		go.BP <-go.BP[order(go.BP$Pvalue),]

		write.table(go.BP, paste0(p[k], "goBP.txt"), quote=F, sep="\t")
		go.MF<-go.enrichment.MF(genes=entrez.id,
		ontology="MF", universe=universe.id,
		pvalue=0.05, annotation='org.Mm.eg.db', conditionalSearch=TRUE, genes2=genes)
		go.MF <-go.MF[order(go.MF$Pvalue),]

		write.table(go.MF, paste0(p[k], "goMF.txt"), quote=F, sep="\t")

		kegg<-kegg.enrichment(genes=entrez.id,universe=universe.id,
		pvalue=0.05, annotation='org.Mm.eg.db', genes2=genes)
		kegg <-kegg[order(kegg$Pvalue),]

		write.table(kegg, paste0(p[k], "kegg.txt"), quote=F, sep="\t")

		if(nrow(go.BP)>10) go.BP <- go.BP[1:10,]
		#go.BP$Term <- factor(go.BP$Term, levels=rev(go.BP$Term))

		# Diverging Barcharts
		go.BP.plot[[i]]<- ggplot(go.BP, aes(x=Term, y= Pvalue , label=FDR)) + 
		geom_bar(stat='identity', width=.5,position="dodge")  +
		labs(title= "GO BP") + 
		coord_flip() + theme_bw()

		if(nrow(go.MF)>10) go.MF <- go.MF[1:10,]

		#go.MF$Term <- factor(go.MF$Term, levels=rev(go.MF$Term))
		# 
		go.MF.plot[[i]]<- ggplot(go.MF, aes(x=Term, y= Pvalue , label=FDR)) + 
		geom_bar(stat='identity', width=.5,position="dodge")  +
		labs(title= "GO MF") + 
		coord_flip() + theme_bw()

		if(nrow(kegg)>10) kegg <- kegg[1:10,]

		#kegg$KEGGID <- factor(kegg$KEGGID, levels=rev(kegg$KEGGID))
		# 
		kegg.plot[[i]]<- ggplot(kegg, aes(x=Term, y= Pvalue , label=FDR)) + 
		geom_bar(stat='identity', width=.5,position="dodge")  +
		labs(title= "KEGG") + 
		coord_flip() + theme_bw()
		i <- i+ 1
	}

		
	for (i in 1:10)
	{
		pdf(paste0("p", i,"_enrichment.pdf"))
		par(mfrow=c(1, 3)) 
		print(go.BP.plot[[i]])
		print(go.MF.plot[[i]])
		print(kegg.plot[[i]])
		dev.off()
	}

}




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

  num_steps <- ceiling(log2(k + 1))
  S <- mat
  for (p in seq(1, num_steps)){
    k_step <- min(2^p - 1, k)
      message(paste0('Step ', p, '/', num_steps, ': ',
                     'Smoothing using k=', k_step))
    D <- 1-cor(S);
    S <- smoother_aggregate_nearest_nb(mat, D, k_step + 1)
  }

  S
}

