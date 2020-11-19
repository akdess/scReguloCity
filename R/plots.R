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


plotClusteredPatterns <- function(object)
{
  all <- object@gradGenes

  filt <- names(which(sort(table(all$finalPattern))>1))

  sub.p <- all[order(all$rayleigh.pval, decreasing=F), ]
  sub.p2 <- sub.p[sub.p$finalPattern %in% filt, ]
  sub.p <- sub.p2[order(sub.p2$rayleigh.pval, decreasing=F), ]
  #sub.p <- sub.p[order(sub.p$classes, decreasing=F), ]
  sub.p$l <- paste0(sub.p$geneId, "_", sub.p$clusterId)
  p <- unique(sub.p$finalPattern)

  cent <- NULL
  for (k in 1:length(p)) {
    s <- sub.p[sub.p$finalPattern==p[k],]
    cent <- cbind(cent, colMeans(maskedMat2[rownames(maskedMat2) %in%  as.character(s$l), , drop = FALSE]))

  }

  masked.mat<- object@gradGenesMasked
  angles <- as.vector(unlist(lapply(split(sub.p, sub.p$finalPattern2), function(x) mean(x$meanAngle))))
  #ang <-lapply(split(sub.p, sub.p$finalPattern2), function(x) direction.l[[as.character(x$geneId)]][[x$clusterId]))

  df <- NULL
  df$d1 <- direction.l.k[[as.character(sub.p$geneId[k])]][[sub.p$clusterId[k]]]
  #  plot.circular( df$d1)

  df <- data.frame(df)
     
  m<- do.call(k, rbind)

  c <- circular(k[[1]], type = "direction", units="radians")
  df <-  NULL
  df$d1 <- k[[1]]
  df <- as.data.frame(df)

  p2<- ggplot(df,aes(x=d1)) + 
  stat_density(
      geom = "tile", 
      aes(fill = ..density..)

  ) + 
  scale_fill_gradientn(colours=rev(rainbow(32))) + 
  coord_polar() + 
  scale_x_continuous(limits = c(0, 360))  + 
  theme_void()+ theme(legend.position = "none")

  sub.p$spatial.adjust.log10 <- (-log10(sub.p[,"spatial.adjust"]))
  ggplot(sub.p, aes(finalPattern, spatial.adjust.log10)) + geom_boxplot()+ theme(axis.text.x = element_text(angle = 90))



  itr <- 1
  pAll <- list()
  pdf("pattern.summary.updated.all.pdf")
  par(mfrow=c(3,4)) 
  #for (k in c(2, 4, 5, 6, 7, 8, 9, 14, 15)) {
  for (k in 1:length(p)) {

    # plot(emb[,1:2],pch=21,col=as.factor(cell.types[rownames(emb)]), cex=0.8,xlab='',ylab='',axes=F); box();
    plot(emb[,1:2],pch=21,col=val2col(cent[,k],gradientPalette=NULL),
    bg=val2col(cent[,k],gradientPalette=NULL),cex=0.4,xlab='tsne1',ylab='tsne2',main=paste0("Pattern:", itr, ""));



  df <-  NULL
  df$d1 <- ang[[k]]
  df <- as.data.frame(df)

  p2<- ggplot(df,aes(x=d1)) + 
  stat_density(
    geom = "tile", 
    aes(fill = ..density..)

  ) + 
  scale_fill_gradientn(colours=rev(rainbow(32))) + 
  coord_polar() + 
  scale_x_continuous(limits = c(0, 360))  + 
  theme_void()+ theme(legend.position = "none")




  pAll[[itr]] <- plot_grid(base2grob(~       plot(emb[,1:2],pch=21,col=val2col(cent[,k],gradientPalette=NULL),
  bg=val2col(cent[,k],gradientPalette=NULL),
   cex=0.4,xlab='',ylab='',main=paste0("Pattern:", itr, ""))
  ), p2,   rel_widths  = c(4,1))

    itr <- itr + 1
  }
  dev.off()


  pdf("finalChromaffin_vfinal_v33.pdf")
  par(mfrow=c(2,1)) 
  for(k in 1:10) {
    print(pAll[[k]])
  }
  dev.off()


  pdf("pattern.summary.2.pdf")
  par(mfrow=c(4,1)) 
  for (k in 1:length(p)) {
   # plot(emb[,1:2],pch=21,col=as.factor(cell.types[rownames(emb)]), cex=0.8,xlab='',ylab='',axes=F); box();
    plot(emb[,1:2],pch=21,col=val2col(cent[,k],gradientPalette=NULL),
      bg=val2col(cent[,k],gradientPalette=NULL),
       cex=0.4,xlab='',ylab='',axes=F, main=p[k]); box();
  }
  dev.off()


}
