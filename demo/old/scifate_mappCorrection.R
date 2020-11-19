

pat <- function(k){
 if(k<0) return("down");
 if(k>0) return("up")
}

library(dbscan)
library(igraph)
library(circular)
library(arules)
library(reshape2)
library(wvtool)
library(ape)
library(spdep)
library(OpenImageR)
library(Rmagic)
library(ggplotify) 
library(gridExtra)
library(KEGG.db)
library(WriteXLS)
library(org.Mm.eg.db)
library(ggplot2)
library(foreach)
library(DESeq2)
require(BiocParallel)
require(RcisTarget)
library(dplyr)  
library(glmnet)
library(raster)

source("C:\\Users\\aharmanci\\Google Drive\\uthealth\\codebase_yale\\codebase_uthealth\\scPathVelo\\utilities.R")
setwd("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate")

setwd("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate")

load("scifate.rda")
load("190501_data.RData")


setwd("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate")
load("scRegulo.scifate.rda")

res <- object@vel
colnames(res) <- gsub("-", ".",colnames(res))

library(DESeq2)

colnames(gene_count_all) <- gsub("-", ".", colnames(gene_count_all))
colnames(gene_count_new) <- gsub("-", ".", colnames(gene_count_new))

umi_sum <- apply(as.matrix(gene_count_all[, colnames(res)]), 2, sum)
gene_count_all.f <- gene_count_all[, -which(umi_sum<10000)]

gene.filt1 <- apply(as.matrix(gene_count_all[, ]), 1, function(x) (length(which(x==0))/length(x))>0.9)
gene.filt2 <- apply(as.matrix(gene_count_new[, ]), 1, function(x) (length(which(x==0))/length(x))>0.9)
genes.filt<- intersect(names(gene.filt1)[gene.filt1], names(gene.filt2)[gene.filt2])

gene_count_all.f <- gene_count_all.f[!rownames(gene_count_all.f) %in% genes.filt, ]
gene_count_new <- gene_count_new[!rownames(gene_count_new) %in% genes.filt, ]

commonRow<- intersect(rownames(gene_count_all.f), rownames(gene_count_new))
commonCol<- intersect(colnames(gene_count_all.f), colnames(gene_count_new))

gene_count_all <- gene_count_all.f[rownames(gene_count_all.f) %in% commonRow, colnames(gene_count_all.f) %in% commonCol ]
gene_count_new <- gene_count_new[rownames(gene_count_new) %in% commonRow, colnames(gene_count_new) %in% commonCol ]

estAll <- estimateSizeFactorsForMatrix(gene_count_all)
gene_count_all.c <- t(scale(t(log2(gene_count_all/estAll+ 1))))
gene_count_new.c <- t(scale(t(log2(gene_count_new/estAll + 1)))

gene_count_all <- gene_count_all.c 
gene_count_new <-gene_count_new.c

colnames(res) <- gsub("-", ".", colnames(res))

commonC<- intersect(commonCol, colnames(res))
res <- res[, commonC]

df_cell$sample <- gsub("-", ".", df_cell$sample)
cell.colors <- paste0(df_cell$cell_cycle_state, "_", df_cell$Cluster_DEX_module, sep="")
cell.colors <- paste0(df_cell$cell_cycle_state)
names(cell.colors) <- df_cell$sample
cell.types <- cell.colors

#r <- readRDS("df_main_umap.RDS")
#r$sample <- gsub("-", ".", r$sample)

common<- intersect(colnames(res),  df_cell$sample)
#emb <- r[match(common, r$sample), ]
emb <- df_cell[match(common, df_cell$sample), c("umap_1_cell_cycle", "umap_2_cell_cycle")]
rownames(emb) <- common

res <- res[, match(common, colnames(res))]
sds <- apply(res, 1, sd)
res <- res[sds>0.1, ]

  
newVsVel <- rep(-1, dim(res)[1])

for (i in 1:dim(res)[1])
{
  print(i)
  #genes <- links[links$TF==TFS[i], "linked_gene"]
  gene <-  rownames(res)[i]
  gn <- intersect(rownames(gene_count_new), as.character(df_gene[df_gene$gene_short_name==as.character(gene),1]))

  if(length(gn)>0){
  if(length(gn)>1) {
  #expressed <- as.numeric(apply(as.matrix(gene_count_all[gn,  colnames(gene_count_all) %in% rownames(emb)]), 2, mean))
  new <- as.numeric(apply(as.matrix(gene_count_new[gn,  colnames(gene_count_new) %in% rownames(emb)]), 2, mean))
  }
  if(length(gn)==1) {
   # expressed <- as.numeric((as.matrix(gene_count_all[gn, colnames(gene_count_all) %in% rownames(emb)])))
    new <- as.numeric((as.matrix(gene_count_new[gn,colnames(gene_count_new) %in%  rownames(emb)])))
    }
  vel <- res[i, colnames(res) %in% rownames(emb)]
  if( sum(abs(as.numeric(vel)))>0){
    newVsVel[i] <- cor (as.numeric(vel), new)
   # expressedVsNew[i] <- cor (expressed, new)
  #  expressedVsVel[i] <- cor (t(vel), expressed)
  }
  }
}

names(newVsVel) <- rownames(res)
newVsVel <- newVsVel[-which(newVsVel==(-1))]

mappabilty <- read.table("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\mappability\\hg19_gencode.v19_intrex_mmap_seq_stats.txt", header=T )
map <- mappabilty[match(rownames(res), mappabilty$GENE_ID), ]
map$exonMapp <- map$MM_EXON/map$EXON_LENGTH 
map$intronMapp <- map$MM_INTRON/map$INTRON_LENGTH 

exon_all <- apply(map[, c("EXON_A", "EXON_C", "EXON_G", "EXON_T")], 1, sum)
exon_gc <- apply(map[, c( "EXON_C", "EXON_G")], 1, sum)

map$exonGC <- exon_gc/exon_all

intron_all <- apply(map[, c("INTRON_A", "INTRON_C", "INTRON_G", "INTRON_T")], 1, sum)
intron_gc <- apply(map[, c( "INTRON_C", "INTRON_G")], 1, sum)

map$intronGC <- intron_gc/intron_all

common <- intersect(rownames(res), intersect(map[,1], names(newVsVel)))
map <- map[match(common, map[,1]),]
newVsVel <- newVsVel[common]
cor(newVsVel, map$exonMapp)
cor(na.omit(data.frame(newVsVel, map$intronMapp)))
res <- res[common, ]


# set up cut-off values 

# bucketing values into bins
map$tag_exon <- cut(map$exonMapp,
                  breaks= seq(min(map$exonMapp), max(map$exonMapp), length=10), 
                  include.lowest=TRUE, 
                   right=FALSE)
map$tag_intron <- cut(map$intronMapp,
                  breaks= seq(min(na.omit(map$intronMapp)), max(na.omit(map$intronMapp)), length=10), 
                  include.lowest=TRUE, 
                   right=FALSE)


map$gc <- map$intronGC/map$exonGC
map$gc_tag <- cut(map$gc,
                  breaks= seq(min(na.omit(map$gc)), max(na.omit(map$gc)), length=10), 
                  include.lowest=TRUE, 
                   right=FALSE)


map$d <- map$intronMapp/map$exonMapp
map$d_tag <- cut(map$d,
                  breaks= seq(min(na.omit(map$d)), max(na.omit(map$d)), length=10), 
                  include.lowest=TRUE, 
                   right=FALSE)



dat <- data.frame(vel=res[,1], m=map$d, gc=map$gc)
rownames(dat) <- rownames(res)
dat <- na.omit(dat)

vel.cor <- matrix(0, nrow=dim(dat[1]), ncol=dim(res)[2])
rownames(vel.cor) <- rownames(dat)

for (i in 1:dim(res)[2])
{
  print(i)
  dat$vel <- res[rownames(dat) ,i]
  gcCount.loess <- lm((dat$vel)~dat$gc+dat$m) 
  vel.cor[,i] <-resid(gcCount.loess)
}
colnames(vel.cor) <- colnames(res)


newVsVel2 <- rep(-1, dim(vel.cor)[1])

for (i in 1:dim(vel.cor)[1])
{

  #genes <- links[links$TF==TFS[i], "linked_gene"]
  gene <-  rownames(vel.cor)[i]
  gn <- intersect(rownames(gene_count_new), as.character(df_gene[df_gene$gene_short_name==as.character(gene),1]))

  if(length(gn)>0){
  if(length(gn)>1) {
  #expressed <- as.numeric(apply(as.matrix(gene_count_all[gn,  colnames(gene_count_all) %in% rownames(emb)]), 2, mean))
  new <- as.numeric(apply(as.matrix(gene_count_new[gn,  colnames(gene_count_new) %in% rownames(emb)]), 2, mean))
  }
  if(length(gn)==1) {
   # expressed <- as.numeric((as.matrix(gene_count_all[gn, colnames(gene_count_all) %in% rownames(emb)])))
    new <- as.numeric((as.matrix(gene_count_new[gn,colnames(gene_count_new) %in%  rownames(emb)])))
    }
  vel <- vel.cor[i, colnames(vel.cor) %in% rownames(emb)]
  if( sum(abs(as.numeric(vel)))>0){
    newVsVel2[i] <- cor (as.numeric(vel), new)
   # expressedVsNew[i] <- cor (expressed, new)
  #  expressedVsVel[i] <- cor (t(vel), expressed)
  }
  }
}
names( newVsVel2) <- rownames(vel.cor)
map2 <- map[match(names( newVsVel2),map$GENE_ID) , ]
newVsVel <- newVsVel[match(names( newVsVel2), names(newVsVel))]
map2$newVsVel <- newVsVel
map2$newVsVel2 <- newVsVel2
#map2 <- map2[map2$GENE_ID %in% names(which(fs2>=0.1)), ]
df <- melt(map2, measure.vars = c("newVsVel", "newVsVel2"))
p1 <- ggplot(data = df, mapping = aes(value, colour=variable, fill=variable)) + 
    
geom_density(alpha=0.2) + 
labs(x='(intron gc content)/(exon gc content)', y='NewlyExpressedVsVelocity Correlation') +
guides(color=FALSE) + 
#ylim(c(0,1)) + 
theme_minimal() + ggtitle("Before Correction") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


df <- melt(map2, measure.vars = c("newVsVel", "newVsVel2"))
df$correction<- as.character(df$variable)
df$correction[df$variable=="newVsVel"] <- "before"
df$correction[df$variable=="newVsVel2"] <- "after"
    
p1 <- ggplot(data = df, mapping = aes(value, colour=correction, fill=correction)) + 
    
    geom_density(alpha=0.2) + 
    labs(x='Newly Expressed vs Velocity Correlation') +
    guides(color=FALSE) + 
    #ylim(c(0,1)) + 
    theme_minimal() + ggtitle("Mappability and GC Correction") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
