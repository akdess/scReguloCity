plottfCorNet <- function(object) {

  df_result2 <-   object@tfCorNet
  df_result2$weight <- df_result2$cor
  g <- graph.data.frame( df_result2[,1:2], directed=FALSE)
  # Convert edge weights to absolute values
  V(g)$shape <- "circle"


  # Change colour of graph vertices
  V(g)$color <- "yellow"
  V(g)$color[V(g)$name %in% TF_list ] <- "orange"

  #V(g)$vertex.frame.color <- "black"
  V(g)$shape[V(g)$name %in% TF_list] <- "rectangle"

  #E(g)$color <- as.factor(df_result2$pattern)
  #E(g)$weight <- (df_result$weight)

  # Change colour of vertex frames
  V(g)$vertex.frame.color <- "white"
  g <- igraph::simplify(g,remove.loops=TRUE, remove.multiple=F)

  e <- get.edgelist(g,names=F)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
        area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))

  plot(g, layout=l,vertex.size=6, label.color="black", asp=0,
        edge.arrow.size=0.2,   edge.width =  2)
}


plotPatterns <- function(object) {
    all <- object@gradGenes
    emb <- object@emb
    df_result2 <- object@tfCorNet
    all$finalPattern <- paste0(all$pattern, "_", all$classes)

    tfs <- unique(as.character(df_result2[,1]))
    plot(emb[,1:2],pch=21,col=as.factor(cell.types[rownames(emb)]), cex=0.8,xlab='',ylab=''); 

    plot(emb[,1:2],pch=21,col=as.factor(cell.types[rownames(emb)]), 
      cex=0.8,xlab='tsne1',ylab='tsne2'); 
    legend("topright" , pch=21,legend=as.character( unique(as.factor(cell.types[rownames(emb)])),
    fill= c("red", "black", "green", "blue"), col=c("red", "black", "green", "blue")))

    pdf("TF.Targets.final.chromaffin.pdf")
    par(mfrow=c(4,4)) 

    for (k in 1:length(tfs)) {
        genes <- (df_result2 %>% filter(TF %in% tfs[k]))$linked_gene
        genes <- genes [!(genes %in% tfs[k])]
        # plot(emb[,1:2],pch=21,col=as.factor(cell.types[rownames(emb)]), cex=0.8,xlab='',ylab='',axes=F); box();
        plot(emb[,1:2],pch=21,col=val2col(res[tfs[k],],gradientPalette=NULL),
             bg=val2col(res[tfs[k],],gradientPalette=NULL),
             cex=0.4,xlab='tsne1',ylab='tsne2', main=tfs[k]); 
        targets <- colMeans(res[genes, , drop = FALSE])
        plot(emb[,1:2],pch=21,col=val2col(targets,gradientPalette=NULL),
             bg=val2col(targets,gradientPalette=NULL),
             cex=0.4,xlab='tsne1',ylab='tsne2',main=paste0(tfs[k], 
                                                 " targets, ngenes:", length(genes))); 
        
    }
    dev.off()

}


plotClusteredPatterns <- function(object, minGenes=100)
{
  
  cent <- object@pattern_clusters 
  emb <- object@emb
  all <- object@gradGenes
  p <- names(which(table(all$finalPattern)>minGenes))
  numOfGenes <- table(all$finalPattern)[(which(table(all$finalPattern)>minGenes))]

  pdf("pattern.summary.5.pdf")
  par(mfrow=c(2,3)) 
  for (k in 1:length(p)) {

    plot(emb[,1:2],pch=21,col=val2col(cent[,k],gradientPalette=NULL),
    bg=val2col(cent[,k],gradientPalette=NULL),cex=0.4,xlab='tsne1',ylab='tsne2',main=paste0("Pattern:", k, " numOfGenes", numOfGenes[k]));

  }
  dev.off()

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

#' adjust colors, while keeping the vector names
#' 
#' @param x color vector
#' @param alpha transparenscy value (passed to adjustcolors as alpha.f)
#' @param ... parameters passsed to adjustcolor
#' @export
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

#' adjust colors, while keeping the vector names
#' 
#' @param x color vector
#' @param alpha transparenscy value (passed to adjustcolors as alpha.f)
#' @param ... parameters passsed to adjustcolor
#' @export
plotTFNetwork <- function(scRegulo.r.f.RN, TF_list ){

  library(igraph)
  library(qgraph)

  
  df_result2 <- scRegulo.r.f.RN@tfCorNet
  df_result2$weight <- df_result2$cor
  g <- graph.data.frame( df_result2[,1:2], directed=FALSE)

  # Convert edge weights to absolute values
  V(g)$shape <- "circle"
  # Change colour of graph vertices
  V(g)$color <- "yellow"
  V(g)$color[V(g)$name %in% TF_list ] <- "orange"

  #V(g)$vertex.frame.color <- "black"
  V(g)$shape[V(g)$name %in% TF_list] <- "rectangle"

  #E(g)$color <- as.factor(df_result2$pattern)
  #E(g)$weight <- (df_result$weight)

  # Change colour of vertex frames
  V(g)$vertex.frame.color <- "white"
  g <- igraph::simplify(g,remove.loops=TRUE, remove.multiple=F)

  e <- get.edgelist(g,names=F)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
    area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))

  pdf("TF.CorNet.pdf")
  plot(g, layout=l,vertex.size=6, label.color="black", asp=0,
        edge.arrow.size=0.2,   edge.width =  2)
  dev.off()
}

plotTFTargets <-  function (scRegulo.r.f.RN, tfs)
{
  tfs <- sort(unique(scRegulo.r.f.RN@tfCorNet[,1]))
  res <- scRegulo.r.f.RN@vel
  df_result2 <- scRegulo.r.f.RN@tfCorNet
  emb <- scRegulo.r.f.RN@emb
  pdf("TF.TargetsV4.pdf")
  par(mfrow=c(3,2)) 

  for (k in 1:length(tfs)) {
    #  genes <- (df_result2 %>% filter(TF %in% tfs[k]))$linked_gene
      genes <- df_result2[df_result2$TF %in% tfs[k], ]$linked_gene
      
      genes <- as.character(genes [!(genes %in% tfs[k])])
      # plot(emb[,1:2],pch=21,col=as.factor(cell.types[rownames(emb)]), cex=0.8,xlab='',ylab='',axes=F); box();
      plot(emb[,1:2],pch=21,col=val2col(t(res[tfs[k],]),gradientPalette=NULL),
           bg=val2col(t(res[tfs[k],]),gradientPalette=NULL),
           cex=0.5,xlab='UMAP_1',ylab='UMAP_2', main=tfs[k]); 



      targets <- colMeans(res[genes, , drop = FALSE])
      corV <- round(cor(t(res[tfs[k],]), targets), digits=2)
      plot(emb[,1:2],pch=21,col=val2col(targets,gradientPalette=NULL),
           bg=val2col(targets,gradientPalette=NULL),
           cex=0.5,xlab='UMAP_1',ylab='UMAP_2',main=paste0(tfs[k], " cor:", corV, 
                                               " targets, ngenes:", length(genes))); 
      
  }
  dev.off()

}

plotTFTargetsExpression <-  function (scRegulo.r.f.RN, tfs)
{
  res <- scRegulo.r.f.RN@emat
  df_result2 <- scRegulo.r.f.RN@tfCorNet
  
  pdf("TF.Targets4.pdf")
  par(mfrow=c(3,3)) 

  for (k in 1:length(tfs)) {
      genes <- (df_result2 %>% filter(TF %in% tfs[k]))$linked_gene
      genes <- as.character(genes [!(genes %in% tfs[k])])
      # plot(emb[,1:2],pch=21,col=as.factor(cell.types[rownames(emb)]), cex=0.8,xlab='',ylab='',axes=F); box();
      plot(emb[,1:2],pch=21,col=val2col(t(res[tfs[k],]),gradientPalette=NULL),
           bg=val2col(t(res[tfs[k],]),gradientPalette=NULL),
           cex=0.4,xlab='tsne1',ylab='tsne2', main=tfs[k]); 

      targets <- colMeans(res[genes, , drop = FALSE])
      corV <- round(cor((res[tfs[k],]), targets), digits=2)
      plot(emb[,1:2],pch=21,col=val2col(targets,gradientPalette=NULL),
           bg=val2col(targets,gradientPalette=NULL),
           cex=0.4,xlab='tsne1',ylab='tsne2',main=paste0(tfs[k], " cor:", corV, 
                                               " targets, ngenes:", length(genes))); 
      
  }
  dev.off()

}

plottfCorNet <- function(object) {

  df_result2 <-   object@tfCorNet
  df_result2$weight <- df_result2$cor
  g <- graph.data.frame( df_result2[,1:2], directed=FALSE)
  # Convert edge weights to absolute values
  V(g)$shape <- "circle"


  # Change colour of graph vertices
  V(g)$color <- "yellow"
  V(g)$color[V(g)$name %in% TF_list ] <- "orange"

  #V(g)$vertex.frame.color <- "black"
  V(g)$shape[V(g)$name %in% TF_list] <- "rectangle"

  #E(g)$color <- as.factor(df_result2$pattern)
  #E(g)$weight <- (df_result$weight)

  # Change colour of vertex frames
  V(g)$vertex.frame.color <- "white"
  g <- igraph::simplify(g,remove.loops=TRUE, remove.multiple=F)

  e <- get.edgelist(g,names=F)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
        area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))

  plot(g, layout=l,vertex.size=6, label.color="black", asp=0,
        edge.arrow.size=0.2,   edge.width =  2)
}



#' adjust colors, while keeping the vector names
#' 
#' @param x color vector
#' @param alpha transparenscy value (passed to adjustcolors as alpha.f)
#' @param ... parameters passsed to adjustcolor
#' @export
findPatterns <- function(object, minGenes)
{
  library(cluster)
  library(fpc)

  maskedMat2 <- object@gradGenesMasked
  all <- object@gradGenes
  maskedMat2 <- object@gradGenesMasked
  all$pam_class <- ""
  p <- names(which(table(all$pattern)>minGenes))
  for (k in 1:length(p)) {
      print(k)
      d = dist(maskedMat2[p[k]==(all$pattern),])
      pamk.best <- pamk(d)
      all$pam_class[match(names(pamk.best$pamobject$clustering), all$id2)] <- pamk.best$pamobject$clustering
      
  }
  all$finalPattern <- paste0(all$pattern, "_", all$pam_class)
  all$l <- paste0(all$geneId, "_", all$clusterId)

  object@gradGenes <- all 

  all$finalPattern <- paste0(all$pattern, "_", all$pam_class)

  p <- names(which(table(all$finalPattern)>minGenes))
  cent <- NULL
  for (k in 1:length(p)) {
      s <- all[all$finalPattern==p[k],]
      cent <- cbind(cent, colMeans(maskedMat2[rownames(maskedMat2) %in%  as.character(s$l), , drop = FALSE]))
  }
  object@pattern_clusters <- cent

  library(enrichR)
  enr <- list()
  cent <- object@pattern_clusters 
  emb <- object@emb
  all <- object@gradGenes
  p <- names(which(table(all$finalPattern)>minGenes))
  numOfGenes <- table(all$finalPattern)[(which(table(all$finalPattern)>minGenes))]

  pdf("pattern.summary.pdf")
  par(mfrow=c(2,3)) 
  for (k in 1:length(p)) {

    plot(emb[,1:2],pch=21,col=val2col(cent[,k],gradientPalette=NULL),
    bg=val2col(cent[,k],gradientPalette=NULL),cex=0.4,xlab='tsne1',ylab='tsne2',main=paste0("Pattern:", k, " numOfGenes", numOfGenes[k]));

  }
  dev.off()
  
  for (k in 1:length(p)) {
      genes <-as.character(all[all$finalPattern==p[k],1])
      #if(length(genes[genes %in% tfs])>0){
      print(p[k])
      print(genes[genes %in% tfs])
      print(length(genes))
      dbs <- c("GO_Molecular_Function_2018", "GO_Biological_Process_2018", "KEGG_2019_Human", "WikiPathways_2019_Human", "MSigDB_Hallmark_2020")
      enr[[k]] <- enrichr(as.character(genes), dbs)
      fl <- lapply(enr[[k]], function(x) x[x$Adjusted.P.value<0.05, ])
      fl <- fl[names(which(lapply(fl, function(x) dim(x)[1])>1))]
      if(length(fl)>0) {
        printEnrich(fl, prefix=paste0("pattern", k)) 
        summary <- data.frame(geneSymbol=genes, isTF="no")
        summary$isTF[genes %in% TF_list] <- "yes"
        write.table(summary, paste0("pattern",  k, ".txt"), quote=F, sep="\t", row.names=F)
      }
  }
}



plotClusteredPatterns <- function(object, minGenes=100)
{
  
  cent <- object@pattern_clusters 
  emb <- object@emb
  all <- object@gradGenes
  p <- names(which(table(all$finalPattern)>minGenes))
  numOfGenes <- table(all$finalPattern)[(which(table(all$finalPattern)>minGenes))]

  pdf("pattern.summary.5.pdf")
  par(mfrow=c(2,3)) 
  for (k in 1:length(p)) {

    plot(emb[,1:2],pch=21,col=val2col(cent[,k],gradientPalette=NULL),
    bg=val2col(cent[,k],gradientPalette=NULL),cex=0.4,xlab='tsne1',ylab='tsne2',main=paste0("Pattern:", k, " numOfGenes", numOfGenes[k]));

  }
  dev.off()

}

plotGeneVelocity  <- function(object, geneName="Dock10" )
{
    direction.l <- lapply(object@gradient, function(x) if(!is.null(x$direction.l.k)) x$direction.l.k)
    ids <- lapply(direction.l, function(x) names(x[1][1]))
  
    df <- data.frame(d1= direction.l[grep(paste0(geneName, "_"), ids)][[1]][[1]])

  pdf(paste0(geneName, ".plots.pdf"))
    #par(mfrow=c(1, 4))

    p2<- ggplot(df,aes(x=d1)) + 
    stat_density(
        geom = "tile", 
        aes(fill = ..density..)

    ) + 
    scale_fill_gradientn(colours=rev(rainbow(32))) + 
    coord_polar() + 
    scale_x_continuous(limits = c(0, 360))  + 
    theme_void()+ theme(legend.position = "none")

    p2
    plot(emb[,1:2],pch=21,col=val2col(t(scRegulo.r.f.RN@vel[rownames(scRegulo.r.f.RN@vel) %in% geneName,]),gradientPalette=NULL),
    bg=val2col(t(scRegulo.r.f.RN@vel[rownames(scRegulo.r.f.RN@vel) %in% geneName,]),gradientPalette=NULL),cex=0.4,xlab='tsne1',ylab='tsne2',main=geneName);

     plot(emb[,1:2],pch=21,col=val2col(t(scRegulo.r.f.RN@emat[rownames(scRegulo.r.f.RN@vel) %in% geneName,]),gradientPalette=NULL),
    bg=val2col(t(scRegulo.r.f.RN@emat[rownames(scRegulo.r.f.RN@emat) %in% geneName,]),gradientPalette=NULL),cex=0.4,xlab='tsne1',ylab='tsne2',main=geneName);

    plot(emb[,1:2],pch=21,col=val2col(t(object@gradGenesMasked[ids[grep(paste0(geneName, "_"), ids)][[1]],]),gradientPalette=NULL),
    bg=val2col(t(object@gradGenesMasked[ids[grep(paste0(geneName, "_"), ids)][[1]],]),gradientPalette=NULL),cex=0.4,xlab='tsne1',ylab='tsne2',main=geneName);

    plot(emb[,1:2],pch=21,col=val2col(t(object@gradGenesMasked[ids[grep(paste0(geneName, "_"), ids)][[1]],]),gradientPalette=NULL),
    bg=val2col(t(object@gradGenesMasked[ids[grep(paste0(geneName, "_"), ids)][[1]],]),gradientPalette=NULL),cex=0.4,xlab='tsne1',ylab='tsne2',main=geneName);


Dock10
    gsubs <- lapply(scRegulo.r.f.RN@gradient, function(x) if(!is.null(x$masked)) x$masked)
    g.sub <- gsubs[grep(paste0(geneName, "_"), ids)][[1]][[1]]

    plot(g.sub, layout = as.matrix(scRegulo.r.f.RN@emb[, 1:2]), vertex.size=1, edge.arrow.size=0.2,
    edge.width =  E(g.sub)$weight*0.1 , vertex.label=NA, asp=0)
  
  dev.off()

}

