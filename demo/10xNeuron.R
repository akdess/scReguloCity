library(scReguloCity)


cell.colors <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/cell.colors.rds"))
emb <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/embedding.rds"))


cell.types <- rep("", length(cell.colors))
cell.types[cell.colors=="#FF0000FF"]<- "Differentiation"
cell.types[cell.colors=="#CCFF00FF"] <- "Differentiation"
cell.types[cell.colors=="#00FF66FF"] <- "Chromaffin"
cell.types[cell.colors=="#0066FFFF"] <- "SCP"
cell.types[cell.colors=="#CC00FFFF"] <- "Unknown"
names(cell.types) <- names(cell.colors)
cell.types <- cell.types[rownames(emb)]


motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\mm9-tss-centered-10kb-10species.mc9nr.feather"
#motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\hg19-tss-centered-10kb-7species.mc9nr.feather"

setwd("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\chromaffin_Intex\\Chromaffin_Data\\intrex_stats_ST.txt\\")
#setwd("/mnt/c/Users/aharmanci/Google\ Drive/uthealth/scPathVelo/chromaffin_Intex/Chromaffin_Data/intrex_stats_ST.txt")

intex <- read.delim("intrex_stats_ST.txt", stringsAsFactor=F)
emat2 <- intex[intex$EXON_INTRON_IDENTIFIER=="EXONLY", -(1:2) ]
rownames(emat2) <- intex[intex$EXON_INTRON_IDENTIFIER=="EXONLY", "GENE_NAME"]
emat2 <- emat2[-(1:2),]

nmat2 <- intex[intex$EXON_INTRON_IDENTIFIER=="INTRONLY", -(1:2) ]
rownames(nmat2) <- intex[intex$EXON_INTRON_IDENTIFIER=="INTRONLY", "GENE_NAME"]
nmat2 <- nmat2[-(1:2),]

smat2 <- intex[intex$EXON_INTRON_IDENTIFIER=="INTREXIC", -(1:2) ]
rownames(smat2) <- intex[intex$EXON_INTRON_IDENTIFIER=="INTREXIC", "GENE_NAME"]
smat2 <- smat2[-(1:2),]

emat2 <- as(as.matrix(emat2), "dgCMatrix")
nmat2 <- as(as.matrix(nmat2), "dgCMatrix")
smat2 <- as(as.matrix(smat2), "dgCMatrix")

emat2 <- filter.genes.by.cluster.expression(emat2,cell.colors,min.max.cluster.average = 5)
nmat2 <- filter.genes.by.cluster.expression(nmat2,cell.colors,min.max.cluster.average = 1)


stats <- normalizeAndSmoothIntrexStats (emat=emat2, nmat=nmat2, smat=smat2, cell.colors=cell.colors)
colnames(stats$emat) <- colnames(emat2)
colnames(stats$nmat) <- colnames(nmat2)
colnames(stats$smat) <- colnames(smat2)

mappability <- read.table("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\mappability\\mm10_gencode.vM21_intrex_mmap_seq_stats.txt", header=T)
c <- intersect(colnames(stats$emat), colnames(stats$nmat) )
c <- intersect(c, rownames(emb) )
stats$emat <- stats$emat[, match(c, colnames(stats$emat)) ]
stats$nmat <- stats$nmat[, match(c, colnames(stats$nmat)) ]

emb <- emb[match(c, rownames(emb)), ]
colnames(stats$emat) <- rownames(emb)
colnames(stats$nmat) <- rownames(emb)
c <- intersect(rownames(stats$emat), mappability$GENE_ID )

stats$emat <- stats$emat[match(c, rownames(stats$emat)), ]
stats$nmat <- stats$nmat[match(c, rownames(stats$nmat)), ]
if (!is.null(smat)) stats$smat <- stats$smat[match(c, rownames(stats$smat)), ]
mappability<- mappability[match(c, mappability$GENE_ID), ]

scRegulo <- CreateSCReguloCityObject(emat = stats$emat, nmat = stats$nmat, smat =NULL,
 emb=emb, mappability=mappability, species="mouse" , motif_ref=motif_ref)


plot(scRegulo@emb[,],pch=21,col=val2col( t(scRegulo@vel.corrected["SLC35A3",]) ,gradientPalette=NULL),bg=val2col(t(scRegulo@vel.corrected["SLC35A3",]),gradientPalette=NULL),
  cex=0.8,xlab='',ylab='',axes=F); box();

c <- intersect(colnames(scRegulo@vel), rownames(emb))

emb <- emb[match(c, rownames(emb)), ]
res <- scRegulo@vel.corrected["Serpine2", match(c, colnames(scRegulo@vel))]
res <- t(res)
plot(emb,pch=21,col=val2col(res[],gradientPalette=NULL),bg=val2col(res[],gradientPalette=NULL),
  cex=0.8,xlab='',ylab='',axes=F); box();

plot(emb[,],pch=21,col=val2col( t(scRegulo@vel.corrected["Serpine2",]) ,gradientPalette=NULL),bg=val2col(t(scRegulo@vel.corrected["Serpine2",]),gradientPalette=NULL),
  cex=0.8,xlab='',ylab='',axes=F); box();

plot(emb[,],pch=21,col=val2col(t(res["Serpine2",rownames(emb)]),gradientPalette=NULL),bg=val2col(t(res["Serpine2",rownames(emb)]),gradientPalette=NULL),
     cex=0.8,xlab='',ylab='',axes=F); box();



sds <- apply(scRegulo@vel.corrected, 1, sd)
res <- res[which(sds>0.1), ]

minNumCells <- 50
emb <- scRegulo@emb
kNNdistplot(emb[, 1:2], k = 10)
abline(h=6)
neigh <- 6 

cell.types <- cell.types[rownames(scRegulo@emb)]

scRegulo.r <- applyReguloCity  (scRegulo, minNumCells=5, cell.types=cell.types, neigh=6, minsize=5)
scRegulo.r.f <- filterGenes (object=scRegulo.r, minPval=0.05, minNumCells=5)
scRegulo.r.f.RN <- extractRegulatoryNetwork (object=scRegulo.r.f,  minNumGenesInPattern=5)

GPRASP2