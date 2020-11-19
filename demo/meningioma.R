library(scReguloCity)
setwd("C:\\Users\\aharmanci\\Google Drive\\uthealth\\akash\\scregulocity")

emb <- read.table ("030620.Scell.MN.combined.TSNE.coord.txt")

clusters <- read.table ("combined.meta.data.txt")
clusters <- data.frame(ID=rownames(clusters), clusters)

emb <- data.frame(ID=rownames(emb), emb)

m<- merge(clusters, emb, by.x="ID", by.y="ID")

#motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\mm9-tss-centered-10kb-10species.mc9nr.feather"
motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\hg19-tss-centered-10kb-7species.mc9nr.feather"
intex <- read.delim(".\\Patel_etal_Meningioma_SC\\Ant\\intrex_stats_ST.txt", stringsAsFactor=F)
emat2 <- intex[intex$EXON_INTRON_IDENTIFIER=="EXONLY", -(1:2) ]
rownames(emat2) <- intex[intex$EXON_INTRON_IDENTIFIER=="EXONLY", "GENE_NAME"]
emat2 <- emat2[-(1:2),]
colnames(emat2) <- paste0("frontal__", colnames(emat2))

nmat2 <- intex[intex$EXON_INTRON_IDENTIFIER=="INTRONLY", -(1:2) ]
rownames(nmat2) <- intex[intex$EXON_INTRON_IDENTIFIER=="INTRONLY", "GENE_NAME"]
nmat2 <- nmat2[-(1:2),]
colnames(nmat2) <- paste0("frontal__", colnames(nmat2))

smat2 <- intex[intex$EXON_INTRON_IDENTIFIER=="INTREXIC", -(1:2) ]
rownames(smat2) <- intex[intex$EXON_INTRON_IDENTIFIER=="INTREXIC", "GENE_NAME"]
smat2 <- smat2[-(1:2),]
colnames(smat2) <- paste0("frontal__", colnames(smat2))

intex <- read.delim(".\\Patel_etal_Meningioma_SC\\Post\\intrex_stats_ST.txt", stringsAsFactor=F)
emat <- intex[intex$EXON_INTRON_IDENTIFIER=="EXONLY", -(1:2) ]
rownames(emat) <- intex[intex$EXON_INTRON_IDENTIFIER=="EXONLY", "GENE_NAME"]
emat <- emat[-(1:2),]
colnames(emat) <- paste0("postal__", colnames(emat))

nmat<- intex[intex$EXON_INTRON_IDENTIFIER=="INTRONLY", -(1:2) ]
rownames(nmat) <- intex[intex$EXON_INTRON_IDENTIFIER=="INTRONLY", "GENE_NAME"]
nmat <- nmat[-(1:2),]
colnames(nmat) <- paste0("postal__", colnames(nmat))

smat <- intex[intex$EXON_INTRON_IDENTIFIER=="INTREXIC", -(1:2) ]
rownames(smat) <- intex[intex$EXON_INTRON_IDENTIFIER=="INTREXIC", "GENE_NAME"]
smat <- smat[-(1:2),]
colnames(smat) <- paste0("postal__", colnames(smat))

emat3 <- cbind(emat, emat2)
nmat3 <- cbind(nmat, nmat2)
smat3 <- cbind(smat, smat2)

emat2 <- as(as.matrix(emat3), "dgCMatrix")
nmat2 <- as(as.matrix(nmat3), "dgCMatrix")
smat2 <- as(as.matrix(smat3), "dgCMatrix")


colnames(emat2) <- gsub("\\.1", "", colnames(emat2))
colnames(nmat2) <- gsub("\\.1", "", colnames(nmat2))
colnames(smat2) <- gsub("\\.1", "", colnames(smat2))

emat2 <- emat2[, colnames(emat2) %in% m$ID]
nmat2 <- nmat2[, colnames(nmat2) %in% m$ID]
smat2 <- smat2[, colnames(smat2) %in% m$ID]


m <- m[match(colnames(emat2), m$ID),]
cell.colors <- m$seurat_clusters
names(cell.colors) <- m$ID
#emat <- filter.genes.by.cluster.expression(emat2,cell.colors,min.max.cluster.average = 1)
#nmat <- filter.genes.by.cluster.expression(nmat2,cell.colors,min.max.cluster.average = 0.5)


emat <- emat2
nmat <- nmat2

stats <- normalizeAndSmoothIntrexStats (emat=emat, nmat=nmat, smat=NULL, cell.colors=cell.colors)
colnames(stats$emat) <- colnames(emat2)
colnames(stats$nmat) <- colnames(nmat2)
#colnames(stats$smat) <- colnames(smat2)


mappability <- read.table("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\mappability\\hg19_gencode.v19_intrex_mmap_seq_stats.txt", header=T )
c <- intersect(colnames(stats$emat), colnames(stats$nmat) )
c <- intersect(c, rownames(emb) )
stats$emat <- stats$emat[, match(c, colnames(stats$emat)) ]
stats$nmat <- stats$nmat[, match(c, colnames(stats$nmat)) ]

emb <- read.table ("030620.Scell.MN.combined.TSNE.coord.txt")

emb <- emb[match(c, rownames(emb)),]
colnames(stats$emat) <- rownames(emb)
colnames(stats$nmat) <- rownames(emb)

ge <- data.frame(stats$emat)
ensembl <- unlist(lapply(strsplit(as.character(rownames(stats$emat)), split="\\."), function(x) x[[1]]))
ge$ensembl <- ensembl
symbols <- mapIds(org.Hs.eg.db, keys = ensembl, keytype = "ENSEMBL", column="SYMBOL")
ids <- data.frame(ID=rownames(data.frame(symbols)), data.frame(symbols))
ge2 <- merge(ge, ids, by.x="ensembl", by.y="ID", all.x=T)
ge3<- ge2[, -c(1:3,163)]
rownames(stats$emat)  <- as.vector(ge2$symbols)
rownames(stats$nmat)  <- as.vector(ge2$symbols)

c <- intersect(rownames(stats$emat), mappability$GENE_ID )

stats$emat <- stats$emat[match(c, rownames(stats$emat)), ]
stats$nmat <- stats$nmat[match(c, rownames(stats$nmat)), ]
if (!is.null(smat)) stats$smat <- stats$smat[match(c, rownames(stats$smat)), ]
mappability<- mappability[match(c, mappability$GENE_ID), ]

scRegulo <- CreateSCReguloCityObject(emat = stats$emat, nmat = stats$nmat, smat =NULL,
 emb=emb, mappability=mappability, species="human" , motif_ref=motif_ref)


plot(scRegulo@emb[,],pch=21,col=val2col( t(scRegulo@vel.corrected["Serpine2",]) ,gradientPalette=NULL),bg=val2col(t(scRegulo@vel.corrected["Serpine2",]),gradientPalette=NULL),
  cex=0.8,xlab='',ylab='',axes=F); box();

c <- intersect(colnames(scRegulo@vel), rownames(emb))

emb <- emb[match(c, rownames(emb)), ]
res <- scRegulo@vel.corrected["Serpine2", match(c, colnames(scRegulo@vel))]
res <- t(res)
plot(emb,pch=21,col=val2col(res[],gradientPalette=NULL),bg=val2col(res[],gradientPalette=NULL),
  cex=0.8,xlab='',ylab='',axes=F); box();

plot(emb[,],pch=21,col=val2col( t(scRegulo@vel.corrected["GPRASP2",rownames(emb)]) ,gradientPalette=NULL),bg=val2col(t(scRegulo@vel.corrected["GPRASP2",]),gradientPalette=NULL),
  cex=0.8,xlab='',ylab='',axes=F); box();

plot(emb[,],pch=21,col=val2col(t(res["Serpine2",rownames(emb)]),gradientPalette=NULL),bg=val2col(t(res["Serpine2",rownames(emb)]),gradientPalette=NULL),
     cex=0.8,xlab='',ylab='',axes=F); box();



sds <- apply(scRegulo@vel.corrected, 1, sd)
res <- res[which(sds>0.1), ]

minNumCells <- 50
emb <- scRegulo@emb
kNNdistplot(emb[, 1:2], k = 20)
abline(h=3)
neigh <-3 

cell.types <- cell.colors
cell.types <- cell.types[rownames(scRegulo@emb)]

scRegulo.r <- applyReguloCity  (scRegulo, minNumCells=5, cell.types=cell.types, neigh=3, minsize=10)
scRegulo.r.f <- filterGenes (object=scRegulo.r, minPval=0.05, minNumCells=10)
scRegulo.r.f.RN <- extractRegulatoryNetwork (object=scRegulo.r.f,  minNumGenesInPattern=5)
