library(SnapATAC)
rna_matrix = t(readMM("../../RNA/matrix.mtx"));bc<-read.csv("../../DNA/barcodes.tsv", head=F); rownames(rna_matrix)<-bc[,1]; gmat=Matrix(as.matrix((rna_matrix)), sparse=TRUE); genes<-read.csv("../../RNA/genes.tsv", sep="\t", head=F); colnames(gmat)<-genes[,1]
s<-colSums(gmat)
gmat<-as.matrix(gmat)
f<-scale(log10(s+1)); par(mfrow=c(1,1)); hist(f, breaks=100)
v<- f>=-1 & f<=2
fgmat<-gmat[,v]; fgmat<-Matrix(fgmat, sparse=TRUE); x.sp@gmat<-fgmat; x.sp = makeBinary(x.sp, "gmat"); dim(x.sp@gmat)
x.sp = calJaccard(x.sp, ncell.chunk=18808, mat="gmat", max.var=18808, seed.use=10, norm.method="normOVE", row.center=TRUE, row.scale=TRUE, low.threshold=-5, high.threshold=5, keep.jmat=FALSE,do.par = TRUE, num.cores = 8)
x.sp = runPCA(x.sp, pc.num=50, input.mat = "nmat", method="svd", weight.by.sd = FALSE, center=TRUE, scale=FALSE, seed.use=10)
plotPCA(x.sp, method="elbow")
plotPCA(x.sp, method="pairwise"); par(mfrow=c(1,1))
rna_jar<-x.sp@nmat
rna_pca<-x.sp@smat
write.table(dna_jar, col.names=F, row.names=F, file="rna_jar.xls", quote=F)
write.table(dna_pca, col.names=F, row.names=F, file="rna_pca.xls", quote=F)