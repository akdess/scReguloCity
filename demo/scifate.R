library(scReguloCity)

load("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\scifate.rda")
load("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\190501_data.RData")

emb <- df_cell[, c("umap_1_cell_cycle" ,"umap_2_cell_cycle")]
rownames(emb) <- df_cell$sample

#cell.colors <- paste0(df_cell$cell_cycle_state, df_cell$Cluster_DEX_module, sep="_")
cell.colors <- paste0(df_cell$cell_cycle_state, sep="_")
names(cell.colors) <- df_cell[,1]
#motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\mm9-tss-centered-10kb-10species.mc9nr.feather"
motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\hg19-tss-centered-10kb-7species.mc9nr.feather"


colnames(smat) <- gsub(".bam", "", gsub("sci-A549:", "", colnames(smat)))
colnames(emat) <- gsub(".bam", "", gsub("sci-A549:", "", colnames(emat)))
colnames(nmat) <- gsub(".bam", "", gsub("sci-A549:", "", colnames(nmat)))

common1 <- intersect(df_cell$sample, colnames(smat))
common2 <- intersect(df_cell$sample, colnames(emat))
common3 <- intersect(df_cell$sample, colnames(nmat))

common <- intersect(common1, intersect(common2, common3))
emat2 <- emat[, match(common, colnames(emat))]
nmat2 <- nmat[, match(common, colnames(nmat))]
smat2 <- smat[, match(common, colnames(smat))]

emat2 <- as(as.matrix(emat2), "dgCMatrix")
nmat2 <- as(as.matrix(nmat2), "dgCMatrix")
smat2 <- as(as.matrix(smat2), "dgCMatrix")

emat2 <- filter.genes.by.cluster.expression(emat2,cell.colors,min.max.cluster.average = 0.5)
nmat2 <- filter.genes.by.cluster.expression(nmat2,cell.colors,min.max.cluster.average = 0.1)


stats <- normalizeAndSmoothIntrexStats (emat=emat2, nmat=nmat2, smat=NULL, cell.colors=cell.colors)

colnames(stats$emat) <- colnames(emat2)
colnames(stats$nmat) <- colnames(nmat2)
colnames(stats$smat) <- colnames(smat2)

mappability <- read.table("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\mappability\\hg19_gencode.v19_intrex_mmap_seq_stats.txt", header=T )
c <- intersect(colnames(stats$emat), colnames(stats$nmat) )

emb <- emb[match(c, rownames(emb)), ]
colnames(stats$emat) <- rownames(emb)
colnames(stats$nmat) <- rownames(emb)
c <- intersect(rownames(stats$emat), mappability$GENE_ID )

stats$emat <- stats$emat[match(c, rownames(stats$emat)), ]
stats$nmat <- stats$nmat[match(c, rownames(stats$nmat)), ]
if (!is.null(smat)) stats$smat <- stats$smat[match(c, rownames(stats$smat)), ]
mappability<- mappability[match(c, mappability$GENE_ID), ]

scRegulo <- CreateSCReguloCityObject(emat = stats$emat, nmat = stats$nmat, smat = NULL,
 emb=emb, mappability=mappability, species="human" , motif_ref=motif_ref)

save(scRegulo, file="scRegulo.scifate.rda")

#save(scRegulo, file="scRegulo.scifate.v2.withNmatcutoff.rda")

emb <- df_cell[, c("umap_1_cell_cycle" ,"umap_2_cell_cycle")]
rownames(emb) <- df_cell$sample

c <- intersect(colnames(scRegulo@vel), rownames(emb))

emb <- emb[match(c, rownames(emb)), ]
res <- scRegulo@vel["EZH2", match(c, colnames(scRegulo@vel))]
res <- t(res)
plot(emb[c,],pch=21,col=val2col(res[],gradientPalette=NULL),bg=val2col(res[],gradientPalette=NULL),
  cex=0.8,xlab='',ylab='',axes=F); box();


sds <- apply(scRegulo@vel.corrected, 1, sd)
res <- res[which(sds>0.1), ]

minNumCells <- 50

kNNdistplot(emb[, 1:2], k = 20)
abline(h=0.04)

neigh <- 0.04

cell.colors <- df_cell$cell_cycle_state
names(cell.colors) <- df_cell[,1]


scRegulo.r <- applyReguloCity  (scRegulo, minNumCells=50, cell.types=cell.colors, neigh=0.04, minsize=5)
scRegulo.r.f <- filterGenes (object=scRegulo.r, minPval=0.05, minNumCells=5)
scRegulo.r.f.RN <- extractRegulatoryNetwork (object=scRegulo.r.f,  minNumGenesInPattern=5)
