setwd("/projects/ps-renlab/chz272/02.REACH-seq/05.salk_tissue/14.combined/Combined_matrix/filter_BlackList/20190317")

library(SnapATAC)
dna_matrix = t(readMM("../../DNA_filtered/matrix.mtx"))
bc<-read.csv("../../DNA_filtered/barcodes.tsv", head=F); rownames(dna_matrix)<-bc[,1]
peaks = read.table("../../DNA_filtered/peaks.bed"); peaks.gr = GRanges(peaks[,1], IRanges(peaks[,2], peaks[,3]))
mat=Matrix(as.matrix((dna_matrix)), sparse=TRUE);gc()
### creates snapATAC object, filter DNA bins
x.sp = createSnapFromPmat(mat,  barcode=rownames(mat),  peaks=peaks.gr)
cols<-colSums(mat); scols<-scale(log10(cols))
hist(scols, breaks=100)
x.sp = makeBinary(x.sp, "pmat"); 
x.sp = filterBins(x.sp, low.threshold=-2, high.threshold=2, mat="pmat")
###
dim(x.sp@pmat);dim(mat)
x.sp = calJaccard(x.sp, ncell.chunk=18808, mat="pmat", max.var=18808, seed.use=10, norm.method="normOVE", row.center=TRUE, row.scale=TRUE, low.threshold=-5, high.threshold=5, keep.jmat=FALSE,do.par = TRUE, num.cores = 8)
x.sp = runPCA(x.sp, pc.num=50, input.mat = "nmat", method="svd", weight.by.sd = FALSE, center=TRUE, scale=FALSE, seed.use=10)
plotPCA(x.sp, method="elbow")
plotPCA(x.sp, method="pairwise"); par(mfrow=c(1,1))
dna_jar<-x.sp@nmat
dna_pca<-x.sp@smat
write.table(dna_jar, col.names=F, row.names=F, file="dna_jar.xls", quote=F)
write.table(dna_pca, col.names=F, row.names=F, file="dna_pca.xls", quote=F)