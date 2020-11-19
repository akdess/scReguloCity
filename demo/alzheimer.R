source("C:\\Users\\aharmanci\\Google Drive\\uthealth\\codebase_yale\\codebase_uthealth\\scPathVelo\\utilities.R")
#source("/mnt/c/Users/aharmanci/Google\ Drive/uthealth/codebase_yale/codebase_uthealth/scPathVelo/utilities.R")

library(Seurat) 
library(velocyto.R)
library(RColorBrewer)

setwd("C:\\Users\\aharmanci\\Google Drive\\uthealth\\alzheimer_Scell\\velocyto")
load("velocyto.Alzheimer.rda")

load("rvel.cd.updated.rda")

setwd("C:\\Users\\aharmanci\\Google Drive\\uthealth\\alzheimer_Scell\\tsai_sCELL\\LargeScale\\")

motif_ref="C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\scifate\\hg19-tss-centered-10kb-7species.mc9nr.feather"

all_Large_Scale_events <- NULL
for (j in 1:42){
  
  load(file=paste0(as.character(unique(final$Subcluster)[j]), ".chrMat.rda"))
  all_Large_Scale_events <- rbind(all_Large_Scale_events, finalChrMat)
}


cell.colors <- as.factor(final$broad.cell.type)
names(cell.colors) <- final$IDVel
cols <-  data.frame(row.names = levels(final$broad.cell.type), color = brewer.pal(nlevels(final$broad.cell.type), name = 'Set1'))

cell.colors <- cols[match(final$broad.cell.type, rownames(cols)), 1]
names(cell.colors) <- final$IDVel

common <- intersect(final$TAG, rownames(all_Large_Scale_events) )
final <- final[match(common, final$TAG),]
all_Large_Scale_events <- all_Large_Scale_events[match(common, rownames(all_Large_Scale_events)), ]

final <- cbind(final, all_Large_Scale_events)

conv.nmat.norm <- rvel.cd$conv.nmat.norm
conv.emat.norm <- rvel.cd$conv.emat.norm

conv.nmat.norm <- conv.nmat.norm[, match(final$IDVel, colnames(conv.nmat.norm))]
conv.emat.norm <- conv.emat.norm[, match(final$IDVel, colnames(conv.emat.norm))]

emb <-final[, c("tsne1", "tsne2")]
rownames(emb) <- final$IDVel


mappability <- read.table("C:\\Users\\aharmanci\\Google Drive\\uthealth\\scPathVelo\\mappability\\hg19_gencode.v19_intrex_mmap_seq_stats.txt", header=T )

scRegulo <- CreateSCReguloCityObject(emat = conv.emat.norm, nmat = conv.nmat.norm, smat =NULL,
 emb=emb, mappability=mappability, species="human" , motif_ref=motif_ref)


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
kNNdistplot(emb[, 1:2], k = 20)
abline(h=6)
neigh <- 2 

cell.types <- cell.colors

scRegulo.r <- applyReguloCity  (scRegulo, minNumCells=5, cell.types=cell.types, neigh=2, minsize=100)
scRegulo.r.f <- filterGenes (object=scRegulo.r, minPval=0.05, minNumCells=5)
scRegulo.r.f.RN <- extractRegulatoryNetwork (object=scRegulo.r.f,  minNumGenesInPattern=5)
