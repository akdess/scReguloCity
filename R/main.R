

filterGenes <- function(object, minPval=0.05)
{

  all <- do.call(rbind , lapply(object@gradient, function(x) if(!is.null(x[[1]])) do.call(rbind, x$all)))
  all$id2 <- paste0(all$geneId, "_", all$clusterId)
#  rownames(all) <- all$id2
  masked.mat <- do.call(rbind , lapply(object@gradient, function(x) if(!is.null(x$masked.mat)) do.call(rbind, x$masked.mat)))
#  rownames(masked.mat) <- all$id2 

  all$adjust <- p.adjust(all$rayleigh.pval, method="fdr")
  all$spatial.adjust <- p.adjust(all$spatcor.pval ,method="fdr")
  all$fisherTest.adjust <- p.adjust(all$fisherTest ,method="fdr")

  all <- all[which(all$adjust<minPval),]
  all <- all[which(all$spatial.adjust<minPval),]
  all <- all[which(all$fisherTest.adjust<minPval),]
  all <- all[order(all$fisherTest.adjust),]
  masked.final <- masked.mat[match(all$id2, rownames(masked.mat)), ]

  object@gradGenes <- all
  object@gradGenesMasked <- masked.final

  return(object)
}


pat <- function(k){
 if(is.na(k)) return ("none");
 if(k<0) return("down");
 if(k>0) return("up")
}

applyReguloCity <- function (object, neigh, minsize=10, minPval=0.05)
{

    cell.types <- object@cell.types
    emb <- object@emb
    res <- object@vel[, rownames(emb)]
    nb_15nn<- frNN(emb[, 1:2], ep=neigh)
    data_list <- nb_15nn$id
    data_list<- lapply(data_list, function(x) if(length(table(cell.types[x]))>1) return(x))
    #data_list <- as.list(data.frame(t(nb_15nn$id)))
    #cl <- parallel::makeCluster(4) 
    #registerDoParallel(cl)
   lapply.res <- foreach(i=1:nrow(res))  %do% {
    #  lapply.res <- lapply(1:nrow(res), function(i) {
            print(i)
            ge <- as.numeric(res[i,])
            data_list2<- lapply(1:length(data_list), function(x) data_list[[x]][which(ge[x]-ge[data_list[[x]]]>2)])
            names(data_list2) <- names(data_list)
            g2 <- graph_from_adj_list(data_list2, mode="out")
            e <- as_edgelist(g2)

                masked <- list()
                direction.l.k <- list()
                masked.mat <- list()
                all <- list()
                itr <- 1
            if(dim(e)[1]>0) {
           
                thr <- c(-1, 1)
                
                tail <- ge[e[,1]]
                head <- ge[e[,2]]

                filter <- (tail <=thr[1] | tail >=thr[2]) & (head <=thr[1] |  head >=thr[2]) & ((tail*head)<0)
                E(g2)$weight <- as.numeric(tail-head)
                
                weight.all <- as.numeric(E(g2)$weight)
                V(g2)$color <- val2col (as.numeric(res[i, ]),gradientPalette=NULL)
                V(g2)$weight <-as.numeric(res[i, ])
               
                E(g2)$color <- val2col (as.numeric(E(g2)$weight),
                    gradientPalette=colorRampPalette(c("gray90","black"),space = "Lab")(1024))
                
                g.sub <- delete_edges(g2, E(g2)[!filter])
                E(g.sub)$arrow[E(g.sub)$weight>0] <- 1
                E(g.sub)$arrow[E(g.sub)$weight<0] <- 2

                 # plot(g.sub, layout = as.matrix(emb[, 1:2]), vertex.size=2, edge.arrow.size=0.5, edge.width =  E(g.sub)$weight , vertex.label=NA)

                 if(length(E(g.sub))>1) {
                      comps <- components(g.sub)
               
                      E(g.sub)$color <- val2col (as.numeric(E(g.sub)$weight),
                        gradientPalette=colorRampPalette(c("gray90","black"), space = "Lab")(1024))

                      clusters <- which(comps$csize>minsize)
                     
                      if(length(clusters)>0) {
                            
                            e2 <- as_edgelist(g.sub)
                            e2.1 <- e2[,1]
                            e2.2 <- e2[,2]
                            x <- emb[e2.1,1:2]
                            y <- emb[e2.2,1:2]
                            x1 <- x[,1]
                            x2 <- x[,2]
                            y1 <- y[,1]
                            y2 <- y[,2]

                            direction <- rep(NA, nrow(x))
                            for (m in 1:nrow(x)) {
                              dx <- x1[m]-y1[m]
                              dy <- x2[m]-y2[m]
                              direction[m] <- atan2(dy, dx)
                            }

                            E(g.sub)$direction <- direction
                      
                    for (cl in 1:length(clusters)) {
                     #lapply.res2 <- foreach(cl=1:length(clusters)) %do% {
                         
                          #filt <- sort(unique(as.numeric(nb_15nn$id[ which(components(g.sub)$membership==clusters[cl]),])))
                          filt <- sort(unique(which(components(g.sub)$membership==clusters[cl])))
                          tab<- table( cell.types [filt])
                          tab <- tab[which(tab>0)]
                          
                          if (length(tab)==2 & sum(tab>10)==2) {
                                  class1_mean <- NA; class2_mean <- NA; 
                                  cells_class1 <- NA;  cells_class2 <- NA;  
                                  fisherTest <- NA; pattern <- "";  
                                  dat <- NULL;
                                  cell.type.means <- NA; 
                                  fDF <- NULL;
                                  
                                  cell.type.means <- sapply(names(tab), 
                                      function(x) 
                                          mean(as.numeric(ge[ intersect(filt, which(cell.types==x))])))

                                  mcls <-res[i, filt]
                                  names(mcls) <- rownames(emb)[filt]
                                 
                                  class1_mean <- cell.type.means[1]
                                  class2_mean <- cell.type.means[2]
                                 
                                  cells_class1 <-  names(mcls)[mcls>0]
                                  cells_class2 <-  names(mcls)[mcls<0]
                                 
                                  tab_class1<- table( cell.types [names(cell.types ) %in% cells_class1])
                                  mclust_class1 <- paste0( paste0(tab_class1[tab_class1>0], collapse=","), ";", paste0(names(tab_class1[tab_class1>0]), collapse=","), sep="")
                                  tab_class2<- table( cell.types[names(cell.types ) %in% cells_class2])
                                  mclust_class2 <- paste0( paste0(tab_class2[tab_class2>0], collapse=","), ";", paste0(names(tab_class2[tab_class2>0]), collapse=","), sep="")
      
                                 

                                  if(dim(tab_class2)==2 & dim(tab_class1)==2 &(class1_mean*class2_mean)<0 & abs(class1_mean)>1 & abs(class2_mean)>1 )
                                  {
                                         fDF <- cbind(tab_class1, tab_class2)
                                       fisherTest<- fisher.test(fDF)$p.val
                                        
                 
                                        d <- !(e2.1 %in% filt &  e2.2 %in% filt)
                                        g.sub2 <- delete_edges(g.sub, E(g.sub)[d])

                                        g <- induced_subgraph(g.sub2, sort(filt))

                                        # plot(g, layout = as.matrix(emb[components(g.sub)$membership==clusters[cl], 1:2]), vertex.size=3,    edge.arrow.size=1, edge.width =  E(g)$weight , vertex.label=NA)

                                        xval <- sum(cos(E(g)$direction) * as.numeric(E(g)$weight))
                                        yval <- sum(sin(E(g)$direction) * as.numeric(E(g)$weight))
                                        meanAngle <- ((atan2(yval, xval)*180/pi) + 360) %% 360 

                                        xv <- cos(E(g)$direction) * as.numeric(E(g)$weight)
                                        yv <- sin(E(g)$direction) * as.numeric(E(g)$weight)
                                        xpro <- atan2(yv, xv)
                                        df <-  NULL
                                        df$d1 <- ((atan2(yv, xv)*180/pi) + 360) %% 360
                                        df <- as.data.frame(df)

                                        # p2<- ggplot(df,aes(x=d1)) + 
                                        #   stat_density(
                                        #     geom = "tile", 
                                        #     aes(fill = ..density..)

                                        #   ) + 
                                        #   scale_fill_gradientn(colours=rev(rainbow(32))) + 
                                        #   coord_polar() + 
                                        #   scale_x_continuous(limits = c(0, 360))  + 
                                        #   theme_void()+ theme(legend.position = "none")

                                        c <- circular(xpro, type = "direction", units="radians")

                                        watson.stat <- watson.test(c)$stat
                                        kuiper.stat <- kuiper.test(c)$stat
                                        rayleigh <- rayleigh.test(c)

                                        sub <- emb[filt, 1:2]
                                        dists <- as.matrix(dist(sub))   
                                        dists.inv <- 1/dists
                                        diag(dists.inv) <- 0
                                        spatcor.pval <- 1
                                        spatcor.observed <- NA
                                        spatcor.expected <- NA

                                        try({mor <- Moran.I(as.numeric(res[i,filt]), dists.inv)  
                                        spatcor.pval <- mor$p.value
                                        spatcor.observed <- mor$observed
                                        spatcor.expected <- mor$expected
                                        }, silent=T)      
                                   
                    
                                      pattern<- paste0(names(tab), "_", sapply(names(tab), 
                                      function(x) 
                                          pat(k=mean(as.numeric(ge[ intersect(filt, which(cell.types==x))])))), collapse=",")


                                      dat <- data.frame(geneId=rownames(res)[i], 
                                      NumOfCells=length(filt),
                                      numbers= paste( tab, collapse=","), 
                                      clusters= paste( names(tab), collapse=","), 
                                      pattern=pattern,
                                      fisherTest=fisherTest,
                                      numOfTsneClusters=length(which(tab>0)), 
                                      spatcor.pval=spatcor.pval,
                                       spatcor.observed=spatcor.observed,
                                       spatcor.expected=spatcor.expected,
                                       #direction.classes=direction.classes,
                                       scale=neigh,
                                      cluster=cl,
                                      clusterId =itr,
                                      meanAngle=meanAngle,
                                      cellType1_mean=class1_mean, 
                                      cellType2_mean=class2_mean, 
                                      mclust_class1=mclust_class1,
                                      mclust_class2=mclust_class2,
                                      nclass1=tab[1],
                                      nclass2=tab[2],
                                      rayleigh.pval=rayleigh$p.value,
                                     
                                      #classes = classes, 
                                      #watson.pval=watson.pval, 
                                     # kuiper.pval=kuiper.pval,
                                      watson.stat=watson.stat, 
                                      kuiper.stat=kuiper.stat)

                                      all[[itr]] <- dat
                                      names(all)[itr] <- paste0(rownames(res)[i],"_", itr)
              
                                      mas <- rep (0, nrow(emb))
                                      mas [rownames(emb) %in% cells_class1] <- 1
                                      mas [rownames(emb) %in% cells_class2] <- (-1)
                                      masked.mat[[itr]]<- mas
                                      names(masked.mat)[itr] <- paste0(rownames(res)[i],"_", itr)
              
                                      masked[[itr]] <-g.sub2
                                      direction.l.k[[itr]] <-  df$d1
                                      names(masked)[itr] <- paste0(rownames(res)[i],"_", itr)
                                      names(direction.l.k)[itr] <- paste0(rownames(res)[i],"_", itr)
                                        
                                 #     l<- list(all=dat, masked.mat=mas, masked = g.sub2, direction.l.k=df$d1)        
                                      itr <- itr +1          
                                          
                                  }    
                             }
                             #return(l)
                          }
                   } 
                }
          }
            return(list(all=all, masked.mat=masked.mat, masked = masked, direction.l.k=direction.l.k))

    }

    #stopCluster(cl)

    object@gradient <- lapply.res
    object<- filterGenes (object, minPval=0.05)

    return(object)
}

# modified from RcisTarget package
extractRegulatoryNetwork <- function(object,  minNumGenesInPattern=5, motifAnnotations=motifAnnotations)
{

    motifRankings <- importRankings(object@motif_ref)
    gene_names = colnames(motifRankings@rankings)
    TF_list = unique(motifAnnotations$TF)
    Cor.table.filt.all <- NULL
    f <- object@vel[ as.character(unique(object@gradGenes$geneId)),]
    Norm.interest.corr <- corr.test( t(f), method="pearson", ci=F)  
    Norm.interest.corr$p[!(rownames(Norm.interest.corr$p) %in% TF_list), ]=NA
    Pval.adj <- as.data.frame(as.table(Norm.interest.corr$p))
    Norm.interest.corr$r [!(rownames(Norm.interest.corr$r) %in% TF_list)]=NA
    Correlation <- as.data.frame(as.table(Norm.interest.corr$r))
    Cor.table <- na.exclude(cbind( Correlation, Pval.adj))[,c(1,2,3,6)]

    colnames(Cor.table) <- c("TF","gene","cor","p.adj")
    Cor.table.filt.all <- Cor.table [((Cor.table[,3])>0.1 & Cor.table[,4] <0.05 ),]
   #Cor.table.filt.all <- Cor.table [((Cor.table[,3])>0.5 & Cor.table[,4] <0.01 ),]


    df_link <- Cor.table.filt.all
    df_link  = df_link %>% dplyr::filter(TF %in% gene_names)
    df_link  = df_link %>% dplyr::filter(gene %in% gene_names)

    object@tfCor <- df_link %>% mutate_if(is.factor, as.character) 
    gene_list <- split(df_link[, 2], df_link[,1])
    gene_list <- gene_list[which(lapply(gene_list, length)>10)]

    df_result = lapply(1:length(gene_list), function(x) {
            print(x)
            TF_id = names(gene_list)[x]
            gene_vector = as.character(gene_list[[x]])
            df_tmp <- NULL
            if(length(gene_vector)>0) {
            df_tmp = TF_link_gene_list(gene_vector, TF_id, motifRankings, motifAnnotations)
          }
            return(df_tmp)
    })

  
    df_result2 = do.call(rbind, df_result)
    df_result2 = df_result2 %>% dplyr::filter(linked_gene != "No_gene")
    colnames(df_link)[2]<- "linked_gene"
    df_result2 = left_join(df_result2, df_link)
    df_result2$TF_link = str_c(df_result2$TF, df_result2$linked_gene, sep = "_")
    df_result2$Conf = "Motif"
    object@tfCorNet <- df_result2
    return(object)
}


# modified from RcisTarget package

TF_link_gene_list <- function(gene_list, TF_name, motifRankings, motifAnnotations) {
    gene_list = as.character(gene_list)
    motifEnrichmentTable_wGenes <- RcisTarget::cisTarget(gene_list, motifRankings,
                                 motifAnnot=motifAnnotations)
    if(nrow(motifEnrichmentTable_wGenes)>0){
      df_tmp_high = motifEnrichmentTable_wGenes %>% filter(str_detect(TF_highConf, TF_name))
      linked_gene_high = (unique(as.character((str_split(df_tmp_high$enrichedGenes,pattern = ";", simplify = T)))))
      linked_gene_high = linked_gene_high[linked_gene_high != ""]
      if(length(linked_gene_high) > 0) {
          df_combine_high = data.frame("TF" = TF_name, "linked_gene" = linked_gene_high, "Conf" = "high")
      } else {
          df_combine_high = data.frame("TF" = TF_name, "linked_gene" = "No_gene", "Conf" = "high")
      }
      
      df_tmp_low = motifEnrichmentTable_wGenes %>% filter(str_detect(TF_lowConf, TF_name))
      linked_gene_low = (unique(as.character((str_split(df_tmp_low$enrichedGenes,pattern = ";", simplify = T)))))
      linked_gene_low = linked_gene_low[!(linked_gene_low %in% linked_gene_high)]
      linked_gene_low = linked_gene_low[linked_gene_low != ""]
      if(length(linked_gene_low) > 0) {
          df_combine_low = data.frame("TF" = TF_name, "linked_gene" = linked_gene_low, "Conf" = "low")
      } else {
          df_combine_low = data.frame("TF" = TF_name, "linked_gene" = "No_gene", "Conf" = "low")
      }
      
      df_combine = rbind(df_combine_high, df_combine_low)
    return(df_combine)
  }
  else {
    return(NULL)
  }
}



clusterGenes <- function(object, minGenes=10)
{
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
  
 
  p <- names(which(table(all$finalPattern)>minGenes))
  cent <- NULL
  for (k in 1:length(p)) {
    s <- all[all$finalPattern==p[k],]
    cent <- cbind(cent, colMeans(maskedMat2[rownames(maskedMat2) %in%  as.character(s$l), , drop = FALSE]))
  }

  object@pattern_clusters <- cent
  return(object)
}




