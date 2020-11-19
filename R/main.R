
applyReguloCity <- function (object, minNumCells, neigh, cell.types, minsize=5)
{
    library(ape)
    emb <- object@emb
    res <- object@vel.corrected
    nb_15nn<- frNN(emb[, 1:2], ep=neigh)
    data_list <- nb_15nn$id
    g2 <- graph_from_adj_list(data_list, mode="out")
    e <- as_edgelist(g2)
    cell.types <- as.character(cell.types)
    lapply.res <- foreach(i=1:nrow(res))  %do% {

            print(i)
            ge <- as.numeric(res[i,])
            thr <- c(-1, 1)
            
            tail <- ge[e[,1]]
            head <- ge[e[,2]]

            filter <- (tail <=thr[1] | tail >=thr[2]) & (head <=thr[1] |  head >=thr[2]) & ((tail*head)<0)
            E(g2)$weight <- as.numeric(tail-head)
           
            weight.all <- as.numeric(E(g2)$weight)
            V(g2)$color <- val2col (as.numeric(res[i, rownames(emb)]),gradientPalette=NULL)
            
            E(g2)$color <- val2col (as.numeric(E(g2)$weight),
              gradientPalette=colorRampPalette(c("gray90","black"), space = "Lab")(1024))
            
            g.sub.t <- delete_edges(g2, E(g2)[!filter])
            g.sub <- delete_edges(g.sub.t, E(g.sub.t)[which((as.numeric(E(g.sub.t)$weight))<2)])


            masked <- list()
            direction.l.k <- list()
            masked.mat <- list()
            all <- list()

            itr <- 1
            
            if(length(E(g.sub))>1) {
                  E(g.sub)$color <- val2col (as.numeric(E(g.sub)$weight),
                    gradientPalette=colorRampPalette(c("gray90","black"), space = "Lab")(1024))

                  comps <- components(g.sub)
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

                  clusters <- which(comps$csize>minsize)
          
                  if(length(clusters)>0) {
                  for (cl in 1:length(clusters)) {
                      #filt <- sort(unique(as.numeric(nb_15nn$id[ which(components(g.sub)$membership==clusters[cl]),])))
                      filt <- sort(unique(which(components(g.sub)$membership==clusters[cl])))
                
                      d <- !(e2.1 %in% filt &  e2.2 %in% filt)
                      g.sub2 <- delete_edges(g.sub, E(g.sub)[d])

                      g <- induced_subgraph(g.sub2, sort(filt))

                      xval <- sum(cos(E(g)$direction) * as.numeric(E(g)$weight))
                      yval <- sum(sin(E(g)$direction) * as.numeric(E(g)$weight))
                      meanAngle <- ((atan2(yval, xval)*180/pi) + 360) %% 360 

                      xv <- cos(E(g)$direction) * as.numeric(E(g)$weight)
                      yv <- sin(E(g)$direction) * as.numeric(E(g)$weight)
                      xpro <- atan2(yv, xv)
                      df <-  NULL
                      df$d1 <- ((atan2(yv, xv)*180/pi) + 360) %% 360
                      df <- as.data.frame(df)
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


                      class1_mean <- NA
                      class2_mean <- NA
                      mclust_class1 <- NA
                      mclust_class2 <- NA

                      mcls <-res[i, filt]
                      names(mcls) <- rownames(emb)[filt]
                      class1_mean <- mean( mcls[mcls>0])
                      class2_mean <- mean(mcls[mcls<0])
                      cells_class1 <-  names(mcls)[mcls>0]
                      cells_class2 <-  names(mcls)[mcls<0]
                      tab2 <- table( cell.types [names(cell.types ) %in% cells_class1])
                      mclust_class1 <- paste0( paste0(tab2[tab2>0], collapse=","), ";", paste0(names(tab2[tab2>0]), collapse=","), sep="")
                      tab3 <- table( cell.types[names(cell.types ) %in% cells_class2])
                      mclust_class2 <- paste0( paste0(tab3[tab3>0], collapse=","), ";", paste0(names(tab3[tab3>0]), collapse=","), sep="")

                      tab<- table( cell.types [filt])


                      if(!(is.na(class1_mean) | is.na(class2_mean)))
                      {
                          
                        velocityMean <- paste0(as.numeric(sapply(names(tab), 
                          function(x) mean(as.numeric(ge[ intersect(filt, which(cell.types==x))])))), collapse=",")

                        aboveOne <- as.vector(unlist((lapply(strsplit(as.character(velocityMean), split=","),
                          function(x) all(abs(as.numeric(x))>1)))))
                        f <- aboveOne | (abs(class1_mean)>1 &  abs(class2_mean)>1) 

                        pass <- f

                        if(length(which(tab>0))>1) { pass <- !(unlist(lapply(strsplit(as.character(velocityMean), split=","), 
                        function(x) all(x>0) | all(x<0)))) & f  }


                              if(pass){

                                if(length(which(tab>0))>1) {
                                      pattern<- paste0(names(tab), "_", sapply(names(tab), 
                                      function(x) 
                                          pat(k=mean(as.numeric(ge[ intersect(filt, which(cell.types==x))])))), collapse=",")
                                 }
                                if(length(which(tab>0))==1) {
                                    pattern <-names(tab)
                                 }     
                                    dat <- NULL
                                    dat <- data.frame(geneId=rownames(res)[i], 
                                    NumOfCells=length(filt),
                                    numbers= paste( tab, collapse=","), 
                                    clusters= paste( names(tab), collapse=","), 
                                    velocityMean=velocityMean, 
                                    pattern=pattern,
                                    numOfTsneClusters=length(which(tab>0)), 
                                    spatcor.pval=spatcor.pval,
                                     spatcor.observed=spatcor.observed,
                                     spatcor.expected=spatcor.expected,
                                     #direction.classes=direction.classes,
                                     scale=neigh,
                                    cluster=cl,
                                    clusterId =itr,
                                    meanAngle=meanAngle,
                                    class1_mean=class1_mean, 
                                    class2_mean=class2_mean, 
                                    class1=mclust_class1,
                                    class2=mclust_class2,
                                    nclass1=length(cells_class1),
                                    nclass2=length(cells_class2),
                                    rayleigh.pval=rayleigh$p.value,
                                   
                                    #classes = classes, 
                                    #watson.pval=watson.pval, 
                                   # kuiper.pval=kuiper.pval,
                                    watson.stat=watson.stat, 
                                    kuiper.stat=kuiper.stat)

                                    all[[itr]] <- dat
                                    
                                    mas <- rep (0, nrow(emb))
                                    mas [rownames(emb) %in% cells_class1] <- 1
                                    mas [rownames(emb) %in% cells_class2] <- (-1)
                                    masked.mat[[itr]]<- mas
                                    names(masked.mat)[itr] <- paste0(rownames(res)[i],"_", itr)
            
                                    masked[[itr]] <-g.sub2
                                    direction.l.k[[itr]] <-  df$d1
                                    names(masked)[itr] <- paste0(rownames(res)[i],"_", itr)
                                    names(direction.l.k)[itr] <- paste0(rownames(res)[i],"_", itr)
                                
                                    itr <- itr +1          
                              }
                      }    
                    }
                 }
               } 

            list(all=all, masked.mat=masked.mat, masked = masked, direction.l.k=direction.l.k)
    }

    object@gradient <- lapply.res
    return(object)
}

filterGenes <- function(object, minPval=0.05, minNumCells=5)
{
  all <- do.call(rbind , lapply(object@gradient, function(x) if(!is.null(x$all)) do.call(rbind, x$all)))
  all$id2 <- paste0(all$geneId, "_", all$clusterId)

  masked.mat <- do.call(rbind , lapply(object@gradient, function(x) if(!is.null(x$masked.mat)) do.call(rbind, x$masked.mat)))
  rownames(masked.mat) <- all$id2 

  all$adjust <- p.adjust(all$rayleigh.pval, method="fdr")
  all$spatial.adjust <- p.adjust(all$spatcor.pval ,method="fdr")

  all <- all[which(all$adjust<minPval),]
  all <- all[which(all$spatial.adjust<minPval),]
  all <- all[all$nclass1>minNumCells & all$nclass2>minNumCells,  ]

  masked.final <- masked.mat[match(all$id2, rownames(masked.mat)), ]

  object@gradGenes <- all
  object@gradGenesMasked <- masked.final

  return(object)
}

extractRegulatoryNetwork <- function(object,  minNumGenesInPattern=5)
{

    if(object@species=="mouse") 
    {
      data(motifAnnotations_mgi)
      motifAnnotations <- motifAnnotations_mgi 
    }
    if(object@species=="human")
    {
      data(motifAnnotations_hgnc)
      motifAnnotations <- motifAnnotations_hgnc
    }

    motifRankings <- importRankings(motif_ref)
    gene_names = colnames(motifRankings@rankings)
    TF_list = unique(motifAnnotations$TF)
    patterns <- names(which(table(object@gradGenes$pattern)>minNumGenesInPattern))
    Cor.table.filt.all <- NULL
    f <- object@vel.corrected[ as.character(unique(object@gradGenes$geneId)),]
    Norm.interest.corr <- corr.test( t(f), method="pearson", ci=F)  
    Norm.interest.corr$p[!(rownames(Norm.interest.corr$p) %in% TF_list), ]=NA
    Pval.adj <- as.data.frame(as.table(Norm.interest.corr$p))
    Norm.interest.corr$r [!(rownames(Norm.interest.corr$r) %in% TF_list)]=NA
    Correlation <- as.data.frame(as.table(Norm.interest.corr$r))
    Cor.table <- na.exclude(cbind( Correlation, Pval.adj))[,c(1,2,3,6)]

    colnames(Cor.table) <- c("TF","gene","cor","p.adj")
    Cor.table.filt.all <- Cor.table [((Cor.table[,3])>0.3 & Cor.table[,4] <0.01 ),]


    df_link <- Cor.table.filt.all
    df_link  = df_link %>% filter(TF %in% gene_names)
    df_link  = df_link %>% filter(gene %in% gene_names)

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
    df_result2 = df_result2 %>% filter(linked_gene != "No_gene")
    colnames(df_link)[2]<- "linked_gene"
    df_result2 = left_join(df_result2, df_link)
    df_result2$TF_link = str_c(df_result2$TF, df_result2$linked_gene, sep = "_")
    df_result2$Conf = "Motif"
    object@tfCorNet <- df_result2
    return(object)
}



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



clusterGenes <- function(object)
{
  all <- object@gradGenes
  for (k in 1:length(unique(all$pattern))) {

      d = dist(maskedMat2[unique(all$pattern)[k]==(all$pattern),])
      pamk.best <- pamk(d)
      all$pam_class[match(names(pamk.best$pamobject$clustering), all$id2)] <- pamk.best$pamobject$clustering

  }

  all$finalPattern <- paste0(all$pattern, "_", all$pam_class)
  object@gradGenes <- all 
  return(object)
}



