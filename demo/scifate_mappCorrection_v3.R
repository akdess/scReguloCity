library(scReguloCity)

load("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\scifate.rda")
load("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\190501_data.RData")

emb <- df_cell[, c("umap_1_cell_cycle" ,"umap_2_cell_cycle")]
rownames(emb) <- df_cell$sample

cell.colors <- paste0(df_cell$cell_cycle_state, df_cell$Cluster_DEX_module, sep="_")
names(cell.colors) <- df_cell[,1]
#motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\mm9-tss-centered-10kb-10species.mc9nr.feather"
motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\hg19-tss-centered-10kb-7species.mc9nr.feather"



 cds_TF_gene_linkage_analysis<- function(cds_all, cds_new, TF_list, output_folder,
                                        gene_filter_rate = 0.1, cell_filter_UMI = 10000,
                                        core_n_lasso, core_n_filtering,
                                        motif_ref, df_gene_TF_link_ENCODE) {
    # process the cds files to generate the TF expression and gene expression matrix for linkage analysis
    input_data = cds_processing_TF_link(cds_all, cds_new, TF_list, 
                                            gene_filter_rate, cell_filter_UMI)
    # Apply LASSO for TF-gene linkage analysis
    mm_1 = input_data[[1]]
    mm_2 = input_data[[2]]
    df_cell = input_data[[3]]
    df_gene = input_data[[4]]
    df_gene_TF = input_data[[5]]
    TF_matrix = input_data[[6]]
    new_size_factor = scale(df_cell$Size_Factor)
    
    tmp = as.data.frame(t(TF_matrix))
    tmp$new_size_factor = new_size_factor
    if(sd(df_cell$labeling_rate) != 0) {
        labeling_ratio = scale(df_cell$labeling_rate)
        tmp$labeling_ratio = labeling_ratio
    }
    tmp = as.matrix(t(tmp))
    
    # link TF and genes based on covariance
    link_result = link_TF_gene_analysis(tmp, mm_2, df_gene_TF, core_num = core_n_lasso)
    save(link_result, file = file.path(output_folder, "raw_links.RDS"))
    # filtering the links using TF-gene binding data and store the result in the target folder
    TF_gene_filter_links(link_result, df_gene, output_folder, core_n_filtering, motif_ref, df_gene_TF_link_ENCODE)
}


link_TF_gene_analysis <- function(TF_matrix, gene_matrix, df_gene_TF, core_num = 10, cor_thresh = 0.03, seed = 123456) {
    
    gene_list = (rownames(gene_matrix))
    
    link_result = bplapply(gene_list, function(x) {
        # if the gene name is in the TF_matrix, then remove this gene in the TF matrix
        if(x %in% df_gene_TF$gene_id) {
            TF_name = (df_gene_TF %>% filter(gene_id == x))$gene_short_name
            input_TF_matrix = TF_matrix[rownames(TF_matrix) != TF_name, ]
        } else {
            input_TF_matrix = TF_matrix
        }
        
        result = TF_gene_link(input_TF_matrix, gene_matrix[as.character(x),], cor_thresh = cor_thresh, seed = seed)
        return(result)
    }, BPPARAM = SnowParam(workers=core_num))
    return(list(gene_list, link_result))
}

cds_processing_TF_link<- function(cds_all, cds_new, TF_list, gene_filter_rate = 0.1, cell_filter_UMI = 10000) {
    gene_filter_1 = (Matrix::rowSums(exprs(cds_all) > 0)) > (gene_filter_rate * nrow(pData(cds_all)))
    gene_filter_2 = (Matrix::rowSums(exprs(cds_new) > 0)) > (gene_filter_rate * nrow(pData(cds_new)))
    cat("\nOriginal gene number: ", nrow(fData(cds_all)))
    cat("\nGene number after filtering: ", sum(gene_filter_1 * gene_filter_2))
    cds_all = cds_all[gene_filter_2 * gene_filter_1 > 0, ]
    cds_new = cds_new[gene_filter_2 * gene_filter_1 > 0, ]
    
    cell_filter = ((Matrix::colSums(exprs(cds_all))) > cell_filter_UMI)
    cat("\nOriginal cell number: ", nrow(pData(cds_all)))
    cds_all = cds_all[, cell_filter]
    cds_new = cds_new[, cell_filter]
    cat("\nCell number after filtering: ", nrow(pData(cds_all)))
    
    # generate the expression matrix for downstream analysis
    cds_all = estimateSizeFactors(cds_all)
    mm_1 = t(t(exprs(cds_all)) / cds_all$Size_Factor)
    mm_2 = t(t(exprs(cds_new)) / cds_all$Size_Factor)

    mm_1 = log(mm_1 + 0.1)
    mm_1 = t(scale(t(mm_1)))

    mm_2 = log(mm_2 + 0.1)
    mm_2 = t(scale(t(mm_2)))
    
    # compute the labeling reads rate in each cell
    df_cell = pData(cds_all)
    df_gene = fData(cds_all)
    df_cell$labeling_rate = (Matrix::colSums(exprs(cds_new))) / (Matrix::colSums(exprs(cds_all)))
    
    # extract the TF matrix
    df_gene_TF = df_gene %>% filter(gene_short_name %in% TF_list)
    cat("\nNumber of TFs found in the list: ", nrow(df_gene_TF))
    
    TF_matrix = mm_1[as.character(df_gene_TF$gene_id), ]
    rownames(TF_matrix) = as.character(df_gene_TF$gene_short_name)
    return(list(mm_1, mm_2, df_cell, df_gene, df_gene_TF, TF_matrix))
}

gene_count_all.c <- scale(log2(gene_count_all/estAll+ 1))
gene_count_new.c <- scale(log2(gene_count_new/estAll + 1))





setwd("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate")
load("scRegulo.scifate.rda")
#load("scifate.rda")
df_gene <- readRDS ("df_gene.RDS")
load("linkage_function.RData")

library(DESeq2)
require(BiocParallel)
require(RcisTarget)
library(dplyr)  
library(glmnet)

data(motifAnnotations_hgnc)
TF_list = unique(motifAnnotations_hgnc$TF)

# define the location of the motif reference downloaded from RcisTarget: https://resources.aertslab.org/cistarget/
# for human motif matrix, it can be downloaded from: https://shendure-web.gs.washington.edu/content/members/cao1025/public/nobackup/sci_fate/data/hg19-tss-centered-10kb-7species.mc9nr.feather
motif_ref = "./hg19-tss-centered-10kb-7species.mc9nr.feather"
# define the location of chip-seq peak matrix downloaded from https://amp.pharm.mssm.edu/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets
df_gene_TF_link_ENCODE = "./df_TF_gene_ENCODE.RData"

# define core number for linkage analysis (core_n_lasso)
core_n_lasso = 5
# define core number for Rcistarget analysis (core_n_filtering)
core_n_filtering = 1

# define output folder
output_folder = "./output/TF_gene_analysis"
dir.create(output_folder)

# load the monocle2 cds object for full gene expression and newly synthesized reads of each cell
# downloaded from: https://shendure-web.gs.washington.edu/content/members/cao1025/public/nobackup/sci_fate/data/cds_all_new.RDS
load("./cds_all_new.RDS")
res <- scRegulo@vel
colnames(res) <- gsub("-", ".",colnames(res))
# check if the cell names and gene names are the same
all(rownames(cds_all) == rownames(cds_new))
all(colnames(cds_all) == colnames(cds_new))

# screening TF ang gene regulatory links based on the covariance between TF expression
# and gene newly synthesis rate as well as the enrichment of TF binding-motif/chip-seq peaks near gene's
# promoter

# Main parameters:
# cds_all: monocle2 cds object for full gene expression of cells
# cds_new: monocle2 cds object for newly synthesized gene expression of cells
# TF_list: gene names of TFs for linkage analysis
# output_links_new_RNA: output folder
# gene_filter_rate: minimum percentage of expressed cells for gene filtering 
# cell_filter_UMI: minimum number of UMIs for cell filtering
# core_n_lasso: number of cores for lasso regression in linkage analysis
# core_n_filtering: number of cores for filtering TF-gene links
# motif_ref: TF binding motif data as described above
# df_gene_TF_link_ENCODE: TF chip-seq data as described above

output_links_new_RNA = file.path(output_folder, "Link_new")
dir.create(output_links_new_RNA)
input_data = cds_processing_TF_link(cds_all, cds_new, TF_list, 
                                     gene_filter_rate=0.1, cell_filter_UMI=10000)

# cds_TF_gene_linkage_analysis(cds_all, cds_new, TF_list, output_links_new_RNA,
#                                         gene_filter_rate = 0.1, cell_filter_UMI = 10000,
#                                         core_n_lasso, core_n_filtering,
#                                         motif_ref, df_gene_TF_link_ENCODE)

 
gene_count_all <- input_data[[1]]
gene_count_new <- input_data[[2]]

all(rownames(gene_count_all) == rownames(gene_count_new))
all(colnames(gene_count_all) == colnames(gene_count_new))

colnames(gene_count_all) <- gsub("-", ".", colnames(gene_count_all))
colnames(gene_count_new) <- gsub("-", ".", colnames(gene_count_new))
df_cell$sample<- gsub("-", ".", df_cell$sample)

# r <- readRDS ("df_main_umap.RDS")

common<- intersect(colnames(gene_count_all), intersect(colnames(res),  df_cell$sample))
#emb <- r[match(common, r$sample), ]
df_cell <- df_cell[match(common, df_cell$sample),]
res <- res[, match(common, colnames(res))]
gene_count_all <- gene_count_all[, match(common, colnames(gene_count_all))]
gene_count_new <- gene_count_new[, match(common, colnames(gene_count_new))]

cell.colors <- paste0(df_cell$cell_cycle_state, "_", df_cell$Cluster_DEX_module, sep="")
cell.colors <- paste0(df_cell$cell_cycle_state)
names(cell.colors) <- df_cell$sample
cell.types <- cell.colors

links <- read.delim("TF_gene_links.txt")
df_cell$sample <- gsub("-", ".", df_cell$sample)
TFS <- c("E2F1", "E2F2", "E2F7", "EZH2", "BRCA1", "E2F8", "FOXM1", "MYBL2" )
commonGenes<- intersect(rownames(gene_count_all),  as.character(df_gene$gene_id))
df_gene <- df_gene[df_gene$gene_id %in% commonGenes, ]

sds <- apply(res, 1, sd)
res <- res[sds>0.1, ]


newVsVel <- rep(-1, dim(res)[1])
expressedVsNew <- rep(-1, dim(res)[1])
expressedVsVel <- rep(-1, dim(res)[1])

for (i in 1:dim(res)[1])
{

  #genes <- links[links$TF==TFS[i], "linked_gene"]
  gene <-  rownames(res)[i]
  gn <- as.character(df_gene[df_gene$gene_short_name==as.character(gene),1])
  if(length(gn)>0) {
	  if(length(gn)>1) {
		  expressed <- as.numeric(apply(as.matrix(gene_count_all[gn,  ]), 2, mean))
		  new <- as.numeric(apply(as.matrix(gene_count_new[gn,  ]), 2, mean))
	  }
	  if(length(gn)==1) {
		  expressed <- as.numeric((as.matrix(gene_count_all[gn, ])))
		  new <- as.numeric((as.matrix(gene_count_new[gn, ])))
	  }
	  vel <- as.numeric(res[i, ])

	  if( sum(abs(as.numeric(vel)))>0){
		  newVsVel[i] <- cor (vel, new)
		  expressedVsNew[i] <- cor (expressed, new)
		  expressedVsVel[i] <- cor (vel, expressed)
	  }
	 }
}



par(mfrow=c(1,3)) 
boxplot(newVsVel[1:2800])
boxplot(expressedVsNew[1:2800])
boxplot(expressedVsVel[1:2800])


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
colnames(vel.cor) <- colname

vel.cor <- scRegulo@vel.corrected
colnames(vel.cor) <- gsub("-", ".",colnames(vel.cor))

vel.cor <- vel.cor[rownames(vel.cor) %in% rownames(res),
colnames(vel.cor) %in% colnames(res)]

vel.cor <- vel.cor[, match(colnames(res), colnames(vel.cor))]

newVsVel2 <- rep(-1, dim(vel.cor)[1])

for (i in 1:dim(vel.cor)[1])
{

  #genes <- links[links$TF==TFS[i], "linked_gene"]
  gene <-  rownames(vel.cor)[i]
  gn <- intersect(rownames(gene_count_new), as.character(df_gene[df_gene$gene_short_name==as.character(gene),1]))

  if(length(gn)>0){
  if(length(gn)>1) {
  #expressed <- as.numeric(apply(as.matrix(gene_count_all[gn,  colnames(gene_count_all) %in% rownames(emb)]), 2, mean))
  new <- as.numeric(apply(as.matrix(gene_count_new[gn,  ]), 2, mean))
  }
  if(length(gn)==1) {
   # expressed <- as.numeric((as.matrix(gene_count_all[gn, colnames(gene_count_all) %in% rownames(emb)])))
    new <- as.numeric((as.matrix(gene_count_new[gn, ])))
    }
  vel <- vel.cor[i, ]
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
