source("C:\\Users\\aharmanci\\Google Drive\\uthealth\\codebase_yale\\codebase_uthealth\\scPathVelo\\utilities.R")
source("/mnt/c/Users/aharmanci/Google\ Drive/uthealth/codebase_yale/codebase_uthealth/scPathVelo/utilities.R")

library(Seurat) 
library(velocyto.R)
#setwd("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\chromaffin_Intex\\Chromaffin_Data\\intrex_stats_ST.txt\\")
setwd("/mnt/c/Users/aharmanci/Google\ Drive/uthealth/scPathVelo/DG/")

setwd("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\DG\\")
emat <- t(readMM("spliced.mtx"))
nmat <- t(readMM("unspliced.mtx"))
smat <- t(readMM("ambiguous.mtx"))
pheno <- read.csv("phenotypes.csv", stringsAsFactor=F)
genes <- read.csv("genes.csv", stringsAsFactor=F)[, 2]

colnames(emat) <- pheno[,1]
colnames(nmat) <- pheno[,1]
colnames(smat) <- pheno[,1]

rownames(emat) <- genes
rownames(nmat) <- genes
rownames(smat) <- genes

postal <- CreateSeuratObject(emat)
postal@meta.data$clusters <- pheno$clusters
emb <- read.table("tsne.txt")
rownames(emb) <- pheno$index

cell.types <-postal@meta.data$clusters
names(cell.types) <- rownames(postal@meta.data)

cell.colors <- as.factor(cell.types)
names(cell.colors) <- rownames(postal@meta.data)
plot(emb[,1:2],pch=21,col=as.factor(cell.types),
    cex=0.8,xlab='',ylab='',axes=F);

library(RColorBrewer)
n <- 14
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col_vector

ct <- as.factor(cell.types[rownames(emb)])
cv <-col_vector[as.numeric(ct)]

plot(emb[,1:2],pch=21, 
  cex=0.8,xlab='tsne1',ylab='tsne2', col=cv); 
legend("topright" , pch=21,legend=as.character(unique(ct)), 
	col=unique(col_vector), fill=unique(col_vector))

library("Matrix")
emat2 <- as(as.matrix(emat), "dgCMatrix")
nmat2 <- as(as.matrix(nmat), "dgCMatrix")
smat2 <- as(as.matrix(smat), "dgCMatrix")

emat2 <- filter.genes.by.cluster.expression(emat2,cell.colors,min.max.cluster.average = 0.5)
nmat2 <- filter.genes.by.cluster.expression(nmat2,cell.colors,min.max.cluster.average = 0.1)


stats <- normalizeAndSmoothIntrexStats (emat=emat2, nmat=nmat2, smat=NULL, cell.colors=cell.colors)

colnames(stats$emat) <- colnames(emat2)
colnames(stats$nmat) <- colnames(nmat2)
colnames(stats$smat) <- colnames(smat2)

mappability <- read.table("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\mappability\\mm10_gencode.vM21_intrex_mmap_seq_stats.txt", header=T)
c <- intersect(colnames(stats$emat), colnames(stats$nmat) )

motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\mm9-tss-centered-10kb-10species.mc9nr.feather"
#motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\hg19-tss-centered-10kb-7species.mc9nr.feather"

emb <- emb[match(c, rownames(emb)), ]
colnames(stats$emat) <- rownames(emb)
colnames(stats$nmat) <- rownames(emb)
c <- intersect(rownames(stats$emat), mappability$GENE_ID )

stats$emat <- stats$emat[match(c, rownames(stats$emat)), ]
stats$nmat <- stats$nmat[match(c, rownames(stats$nmat)), ]
if (!is.null(smat)) stats$smat <- stats$smat[match(c, rownames(stats$smat)), ]
mappability<- mappability[match(c, mappability$GENE_ID), ]

scRegulo <- CreateSCReguloCityObject(emat = stats$emat, nmat = stats$nmat, smat = NULL,
 emb=emb, mappability=mappability, species="mouse" , motif_ref=motif_ref)

save(scRegulo, file="scRegulo.DG.rda")
load("scRegulo.DG.rda")

minNumCells <- 50

kNNdistplot(emb[, 1:2], k = 20)
abline(h=0.8)

neigh <- 0.8

scRegulo.r <- applyReguloCity  (scRegulo, minNumCells=50, cell.types=cell.colors, neigh=0.8, minsize=5)
scRegulo.r.f <- filterGenes (object=scRegulo.r, minPval=0.05, minNumCells=5)
scRegulo.r.f.RN <- extractRegulatoryNetwork (object=scRegulo.r.f,  minNumGenesInPattern=5)
 gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)

tfs <- c("Sox10", "Egr1")
df_result2 <- scRegulo.r.f.RN@tfCorNet
res <- scRegulo.r.f.RN@vel.corrected
pdf("TFs.pdf")
par(mfrow=c(2,2)) 

for (k in 1:length(tfs)) {
    genes <- (df_result2 %>% filter(TF %in% tfs[k]))$linked_gene
    genes <- genes [!(genes %in% tfs[k])]
    # plot(emb[,1:2],pch=21,col=as.factor(cell.types[rownames(emb)]), cex=0.8,xlab='',ylab='',axes=F); box();
    plot(emb[,1:2],pch=21,col=val2col(t(res[tfs[k],rownames(emb)]),gradientPalette=NULL),
         bg=val2col(t(res[tfs[k],rownames(emb)]),gradientPalette=NULL),
         cex=0.4,xlab='tsne1',ylab='tsne2', main=tfs[k]); 
    targets <- colMeans(res[genes, rownames(emb), drop = FALSE])
    plot(emb[,1:2],pch=21,col=val2col(targets,gradientPalette=NULL),
         bg=val2col(targets,gradientPalette=NULL),
         cex=0.4,xlab='tsne1',ylab='tsne2',main=paste0(tfs[k], 
                                             " targets, ngenes:", length(genes))); 
    
}
dev.off()


library(cluster)
library(fpc)
patterns <- names(which(table(all$pattern)>10))

all <- scRegulo.r.f.RN@gradGenes
all$pam_class <- NA
maskedMat2 <-  scRegulo.r.f.RN@gradGenesMasked
for (k in 1:length((patterns))) {

	print(k)
    d = dist(maskedMat2[patterns[k]==(all$pattern),])
    pamk.best <- pamk(d)
    all$pam_class[match(names(pamk.best$pamobject$clustering), all$id2)] <- pamk.best$pamobject$clustering

}

all$finalPattern <- paste0(all$pattern, "_", all$pam_class)
all <- all[!is.na(all$pam_class ),]
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
   #cent <- cbind(cent, apply(maskedMat2[rownames(maskedMat2) %in%  as.character(s$l),, drop = FALSE],2, median ))


}

itr <- 1
pAll <- list()
pdf("pattern.summary.updated.all2..pdf")
par(mfrow=c(3,3)) 
#for (k in c(2, 4, 5, 6, 7, 8, 9, 14, 15)) {
for (k in 1:length(p)) {

# plot(emb[,1:2],pch=21,col=as.factor(cell.types[rownames(emb)]), cex=0.8,xlab='',ylab='',axes=F); box();
plot(emb[,1:2],pch=21,col=val2col(cent[,k],gradientPalette=NULL),
bg=val2col(cent[,k],gradientPalette=NULL),
 cex=0.4,xlab='tsne1',ylab='tsne2',main=paste0("Pattern:", itr, ""));


itr <- itr + 1
}

