
  ## read TSS annotation
  # tss<-read.csv("mm10_tss_all.xls_dedup.xls", sep="\t", head=F); rownames(tss)<-tss[,1]; tss<-tss[,2:4]; colnames(tss)<-c("chr", "pss", "pse")
  
  sample<-"encode_sort.xls_pseudo_cell_E16.xls";
  #sample<-"shuffle_list.xls";
  promoter_table_file<-paste(sample, "TSS_count.xls", sep="_"); expression_table_file<-paste(sample, "TTS_count.xls", sep="_"); distal_table_file<-paste(sample, "REG_count.xls", sep="_")
  
  tss_position<-read.csv("mm10_tss_all.xls_dedup.xls", sep="\t", head=F); rownames(tss_position)<-tss_position[,1]; tss_position<-tss_position[,2:4]; colnames(tss_position)<-c("chr", "pss", "pse");
  promoter_table<-read.csv(promoter_table_file, sep="\t", head=F); rownames(promoter_table)<-promoter_table[,1]; promoter_table<-promoter_table[, 2:dim(promoter_table)[2]];promoter_table<-log(promoter_table+1)
  expression_table<-read.csv(expression_table_file, sep="\t", head=F); rownames(expression_table)<-expression_table[,1]; expression_table<-expression_table[, 2:dim(expression_table)[2]];expression_table<-log(expression_table+1)
  
  reg_table<-read.csv(distal_table_file, sep="\t", head=F); rownames(reg_table)<-reg_table[,1]; reg_table<-reg_table[,2:dim(reg_table)[2]]; reg_table<-reg_table[rowSums(reg_table)!=0,]; reg_table<-log(reg_table+1)
  reg_table<-reg_table[rowSums(reg_table)!=0,]
  reg_position <- data.frame(do.call(rbind, strsplit(rownames(reg_table), "[:-]"))); rownames(reg_position) <- rownames(reg_table); colnames(reg_position)<-c("chr", "pss", "pse")
  reg_position <- as.matrix(reg_position)
  
  ## cal each gene
  cal_gene<-function(gene){
    gene_pos <- tss_position[gene,]; gene_pos<-t(gene_pos); gchr<-gene_pos[1]; gpss<-as.numeric(gene_pos[2]); gpse<-as.numeric(gene_pos[3])
    reg_list<- reg_position[,1]==gchr & (as.numeric(reg_position[,2]) < (gpse + 500000)) & (as.numeric(reg_position[,3]) > (gpss - 500000))
    test_list<-reg_table[reg_list,]
    gene_vector<-expression_table[gene,]; gene_vector<-as.vector(t(gene_vector));
    cor<-array(cor(as.matrix(t(test_list)), gene_vector)); 
    cor_mat<-cbind(rep(gene, dim(test_list)[1]), rownames(test_list), cor)
    colnames(cor_mat)<-c("Gene", "Reg", "Cor")

    
    return(cor_mat)
  }
  first_gene<-rownames(expression_table)[1]
  gene_distal<-cal_gene(first_gene)
  t=0;
  for(i in 2:length(rownames(expression_table))){
    cur<-cal_gene(rownames(expression_table)[i])
    gene_distal<-rbind(gene_distal, cur)
    t=t+1;print(t)
  }
  write.table(gene_distal, col.names=F, row.names=F, quote=F, file="Gene_distal_E16_cells.xls")
  #write.table(gene_distal, col.names=F, row.names=F, quote=F, file="Gene_distal_shuffle_cells.xls")
  
  

  ## read TSS annotation
  # tss<-read.csv("mm10_tss_all.xls_dedup.xls", sep="\t", head=F); rownames(tss)<-tss[,1]; tss<-tss[,2:4]; colnames(tss)<-c("chr", "pss", "pse")
  
  #sample<-"Pseudo_cells.xls";
  sample<-"encode_sort.xls_pseudo_cell_E16.xls";
  promoter_table_file<-paste(sample, "TSS_count.xls", sep="_"); expression_table_file<-paste(sample, "TTS_count.xls", sep="_"); distal_table_file<-paste(sample, "REG_count.xls", sep="_")
  
  tss_position<-read.csv("mm10_tss_all.xls_dedup.xls", sep="\t", head=F); rownames(tss_position)<-tss_position[,1]; tss_position<-tss_position[,2:4]; colnames(tss_position)<-c("chr", "pss", "pse");
  promoter_table<-read.csv(promoter_table_file, sep="\t", head=F); rownames(promoter_table)<-promoter_table[,1]; promoter_table<-promoter_table[, 2:dim(promoter_table)[2]];promoter_table<-log(promoter_table+1)
  expression_table<-read.csv(expression_table_file, sep="\t", head=F); rownames(expression_table)<-expression_table[,1]; expression_table<-expression_table[, 2:dim(expression_table)[2]];expression_table<-log(expression_table+1)
  
  reg_table<-read.csv(distal_table_file, sep="\t", head=F); rownames(reg_table)<-reg_table[,1]; reg_table<-reg_table[,2:dim(reg_table)[2]]; reg_table<-reg_table[rowSums(reg_table)!=0,]; reg_table<-log(reg_table+1)
  reg_table<-reg_table[rowSums(reg_table)!=0,]
  reg_position <- data.frame(do.call(rbind, strsplit(rownames(reg_table), "[:-]"))); rownames(reg_position) <- rownames(reg_table); colnames(reg_position)<-c("chr", "pss", "pse")
  reg_position <- as.matrix(reg_position)
  cal_promoter<-function(gene){
    gene_pos <- tss_position[gene,]; gene_pos<-t(gene_pos); gchr<-gene_pos[1]; gpss<-as.numeric(gene_pos[2]); gpse<-as.numeric(gene_pos[3])
    reg_list<- reg_position[,1]==gchr & (as.numeric(reg_position[,2]) < (gpse + 500000)) & (as.numeric(reg_position[,3]) > (gpss - 500000))
    test_list<-reg_table[reg_list,]
    gene_vector<-promoter_table[gene,]; gene_vector<-as.vector(t(gene_vector));
    cor<-array(cor(as.matrix(t(test_list)), gene_vector)); 
    cor_mat<-cbind(rep(gene, dim(test_list)[1]), rownames(test_list), cor)
    colnames(cor_mat)<-c("Gene", "Reg", "Cor")

    return(cor_mat)
  }
  
  
  
  first_gene<-rownames(promoter_table)[1]
  promoter_distal<-cal_promoter(first_gene)
  t=0;
  for(i in 2:length(rownames(promoter_table))){
    cur<-cal_promoter(rownames(promoter_table)[i])
    promoter_distal<-rbind(promoter_distal, cur)
    t=t+1;print(t)
  }
  #write.table(gene_distal, col.names=F, row.names=F, quote=F, file="Promoter_distal_pseudo_cells.xls")
  write.table(promoter_distal, col.names=F, row.names=F, quote=F, file="Promoter_distal_E16_cells.xls")