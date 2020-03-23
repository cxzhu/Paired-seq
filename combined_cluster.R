library(matrixcalc)
library(uwot)
library(RColorBrewer)
library(igraph)
library(FNN)
library(Matrix)
col2<-brewer.pal(n=12, name="Paired")
col1<-brewer.pal(n=8, name="Dark2")
col<-c(col2,col1)
setwd("/projects/ps-renlab/chz272/02.REACH-seq/05.salk_tissue/14.combined/Combined_matrix/filter_BlackList/20190318")
dna_jar<-read.csv("dna_jar.xls", sep=" ", head=F)
rna_jar<-read.csv("rna_jar.xls", sep=" ", head=F)

dna_nonNeg<-dna_jar + min(dna_jar)
rna_nonNeg<-rna_jar + min(rna_jar)
had_jar<-hadamard.prod(dna_nonNeg, rna_nonNeg)
pc.had<-prcomp_irlba(had_jar, n=30)

par(mfrow=c(3,3))
for(i in 1:9){
  plot(pc.had$x[,(i*2-1):(i*2)], pch=16, cex=0.1, col="grey")
}
par(mfrow=c(1,1))
set.seed(131); umap.had=umap(pc.had$x[,1:9], n_neighbors=15 , min_dist=0.01, verbose=T, n_threads=8, ret_nn = T, nn_method="FNN")
plot(umap.had$embedding, col="grey", pch=16, cex=0.25)

clust2<-function(matrix1, matrix2){
  k=200
  knn.norm = FNN::get.knn(as.matrix(matrix1), k = k/2)
  knn.norm2 = FNN::get.knn(as.matrix(matrix2), k=k/2) 
  knn.norm$nn.index <- cbind(knn.norm$nn.index, knn.norm2$nn.index)
  knn.norm$nn.dist <- cbind(knn.norm$nn.dist, knn.norm2$nn.dist)
  knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),k), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
  nw.norm = graph_from_data_frame(knn.norm, directed = FALSE)
  nw.norm = simplify(nw.norm)
  lc.norm.rna = cluster_louvain(nw.norm)
  return(lc.norm.rna)
}

clust<-clust2(pc.had$x[,1:9])
