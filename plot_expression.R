{

  prom<-read.csv("filt_TSS_enrichment.xls", sep="\t", head=F)
  expr<-read.csv("filt_TTS_enrichment.xls", sep="\t", head=F)
  
  r.s2<-expr[order(prom[,2], decreasing=TRUE),]
  d.s2<-prom[order(prom[,2], decreasing=TRUE),]
  x<-r.s2[,2];y<-d.s2[,2]
  plot(x, pch=16)
  
  d_merge<-c(prom[,2], prom[,3], prom[,4], prom[,5], prom[,5], prom[,7], prom[,8], prom[,9], prom[,10])
  r_merge<-c(expr[,2], expr[,3], expr[,4], expr[,5], expr[,5], expr[,7], expr[,8], expr[,9], expr[,10])
  
  d_filt<-d_merge[d_merge*r_merge!=0]
  r_filt<-r_merge[d_merge*r_merge!=0]

  
  d_sort<-d_filt[order(d_filt, decreasing=FALSE)]
  r_sort<-r_filt[order(d_filt, decreasing=FALSE)]
  
  m<-matrix(r_sort[1:51575], ncol=5); 
  boxplot(m, pch=16, cex=0.2, ylim=c(0,350))
  
  
  d_mean<-rowMeans(prom[2:10])
  r_mean<-rowMeans(expr[2:10])
  
  dm_sort<-d_mean[order(d_mean, decreasing=FALSE)]
  rm_sort<-r_mean[order(d_mean, decreasing=FALSE)]
  
  dm_filt<-dm_sort[r_mean!=0]
  rm_filt<-rm_sort[r_mean!=0]
  
  mm<-matrix(c(rm_filt[1:12697]), ncol=10)
  colnames(mm)=c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")

  col_grey<-colorRampPalette(c("grey50", "grey85"))(10)
  pdf("All_cells_10_quantile.pdf")
  boxplot(mm, pch=16, cex=0, ylim=c(0,250), medianLwd=0.1, col=col_grey, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="All cells")#, labels=c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10"))
  dev.off()
  
  ### Asc
  d_asc<-prom[,2]
  r_asc<-expr[,2]
  
  rs_asc<-r_asc[order(d_asc, decreasing=FALSE)]
  ds_asc<-d_asc[order(d_asc, decreasing=FALSE)]
  rs_asc<-rs_asc[r_asc!=0]
  m_asc<-matrix(rs_asc, ncol=4);colnames(m_asc)<-c("Q1", "Q2", "Q3", "Q4")
  col<-colorRampPalette(c("firebrick3", "black"))(10); boxplot(m_asc, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="AS")
  
  ### Mic
  d_asc<-prom[,3]
  r_asc<-expr[,3]
  rs_asc<-r_asc[order(d_asc, decreasing=FALSE)]
  ds_asc<-d_asc[order(d_asc, decreasing=FALSE)]
  rs_asc<-rs_asc[r_asc!=0]
  m_mic<-matrix(rs_asc, ncol=4);colnames(m_mic)<-c("Q1", "Q2", "Q3", "Q4")
  col<-colorRampPalette(c("orange4", "black"))(10); boxplot(m_mic, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="MG")
  
  ### Oligodendrocyte
  d_asc<-prom[,4]
  r_asc<-expr[,4]
  rs_asc<-r_asc[order(d_asc, decreasing=FALSE)]
  ds_asc<-d_asc[order(d_asc, decreasing=FALSE)]
  rs_asc<-rs_asc[r_asc!=0]
  m_oli<-matrix(rs_asc, ncol=4);colnames(m_oli)<-c("Q1", "Q2", "Q3", "Q4")
  col<-colorRampPalette(c("gold", "black"))(10); boxplot(m_oli, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="OC")
  
  ### Ex1
  d_asc<-prom[,5]
  r_asc<-expr[,5]
  rs_asc<-r_asc[order(d_asc, decreasing=FALSE)]
  ds_asc<-d_asc[order(d_asc, decreasing=FALSE)]
  rs_asc<-rs_asc[r_asc!=0]
  m_ex1<-matrix(rs_asc, ncol=4);colnames(m_ex1)<-c("Q1", "Q2", "Q3", "Q4")
  col<-colorRampPalette(c("green3", "black"))(10); boxplot(m_ex1, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="Ex1")
  
  ### Ex2
  d_asc<prom[,6]
  r_asc<-expr[,6]
  rs_asc<-r_asc[order(d_asc, decreasing=FALSE)]
  ds_asc<-d_asc[order(d_asc, decreasing=FALSE)]
  rs_asc<-rs_asc[r_asc!=0]
  m_ex2<-matrix(rs_asc, ncol=4);colnames(m_ex2)<-c("Q1", "Q2", "Q3", "Q4")
  col<-colorRampPalette(c("green4", "black"))(10); boxplot(m_ex2, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="Ex2")
  
  ### Ex3
  d_asc<prom[,7]
  r_asc<-expr[,7]
  rs_asc<-r_asc[order(d_asc, decreasing=FALSE)]
  ds_asc<-d_asc[order(d_asc, decreasing=FALSE)]
  rs_asc<-rs_asc[r_asc!=0]
  m_ex3<-matrix(rs_asc, ncol=4);colnames(m_ex3)<-c("Q1", "Q2", "Q3", "Q4")
  col<-colorRampPalette(c("darkgreen", "black"))(10); boxplot(m_ex3, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="Ex3")
  
  ### In1
  d_asc<prom[,8]
  r_asc<-expr[,8]
  rs_asc<-r_asc[order(d_asc, decreasing=FALSE)]
  ds_asc<-d_asc[order(d_asc, decreasing=FALSE)]
  rs_asc<-rs_asc[r_asc!=0]
  m_in1<-matrix(rs_asc, ncol=4);colnames(m_in1)<-c("Q1", "Q2", "Q3", "Q4")
  col<-colorRampPalette(c("lightskyblue", "black"))(10); boxplot(m_in1, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="In1")
  
  ### In2
  d_asc<prom[,9]
  r_asc<-expr[,9]
  rs_asc<-r_asc[order(d_asc, decreasing=FALSE)]
  ds_asc<-d_asc[order(d_asc, decreasing=FALSE)]
  rs_asc<-rs_asc[r_asc!=0]
  m_in2<-matrix(rs_asc, ncol=4);colnames(m_in2)<-c("Q1", "Q2", "Q3", "Q4")
  col<-colorRampPalette(c("deepskyblue3", "black"))(10); boxplot(m_in2, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="In2")
  
  ### Ex3
  d_asc<prom[,10]
  r_asc<-expr[,10]
  rs_asc<-r_asc[order(d_asc, decreasing=FALSE)]
  ds_asc<-d_asc[order(d_asc, decreasing=FALSE)]
  rs_asc<-rs_asc[r_asc!=0]
  m_in3<-matrix(rs_asc, ncol=4);colnames(m_in3)<-c("Q1", "Q2", "Q3", "Q4")
  col<-colorRampPalette(c("deepskyblue4", "black"))(10); boxplot(m_in3, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="In3")
  
  pdf("clusters_quantile.pdf")
  par(mfrow=c(2,5))
  col<-colorRampPalette(c("firebrick3", "white"))(10); boxplot(m_asc, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="AS")
  col<-colorRampPalette(c("orange4", "white"))(10); boxplot(m_mic, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="MG")
  col<-colorRampPalette(c("gold", "white"))(10); boxplot(m_oli, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="OC")
  col<-colorRampPalette(c("green3", "white"))(10); boxplot(m_ex1, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="Ex1")
  col<-colorRampPalette(c("green4", "white"))(10); boxplot(m_ex2, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="Ex2")
  col<-colorRampPalette(c("darkgreen", "white"))(10); boxplot(m_ex3, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="Ex3")
  col<-colorRampPalette(c("lightskyblue", "white"))(10); boxplot(m_in1, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="In1")
  col<-colorRampPalette(c("deepskyblue3", "white"))(10); boxplot(m_in2, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="In2")
  col<-colorRampPalette(c("deepskyblue4", "white"))(10); boxplot(m_in3, pch=16, cex=0.0, ylim=c(0,300), col=col, xlab="Promoter accessibility (from low to high)", ylab="Gene expression (TPM)", main="In3")
  par(mfrow=c(1,1))
  dev.off()
  

}

