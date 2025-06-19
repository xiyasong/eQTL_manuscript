install.packages("xlsx")
library('xlsx')
library(dplyr)
library(data.table)
library(biomaRt)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                     mirror = "useast")




KM_0.05$simple_gene_id <- sapply(strsplit(KM_0.05$gene,'\\.'),'[',1)

list <- getBM(attributes=c('gene_biotype','ensembl_gene_id','external_gene_name'), 
              filters = 'ensembl_gene_id',
              values = KM_0.05$simple_gene_id, mart = ensembl)


for (row in 1:nrow(KM_0.05)) {
    KM_0.05$gene_name[row] <- Japan_eqtl_0.05[Japan_eqtl_0.05$gene_id == KM_0.05$gene[row],]$gene_name
    if (KM_0.05$simple_gene_id[row] %in% list$ensembl_gene_id){
    KM_0.05$gene_biotype[row] <- list[list$ensembl_gene_id ==KM_0.05$simple_gene_id[row],]$gene_biotype
    }
}

KM_0.05 <- KM_0.05[order(KM_0.05$p),]
write.xlsx(KM_0.05,"/Users/mstermo/eQTL-continue/All_supple_materials/Univariate-significant_eGenes.xlsx", sheetName = "KM", 
           col.names = TRUE, row.names = FALSE, append = TRUE)



##########annotate cox 
colnames(cox_0.05)[1] <- "gene"
cox_0.05$simple_gene_id <- sapply(strsplit(cox_0.05$gene,'\\.'),'[',1)

list <- getBM(attributes=c('gene_biotype','ensembl_gene_id','external_gene_name'), 
              filters = 'ensembl_gene_id',
              values = cox_0.05$simple_gene_id, mart = ensembl)


for (row in 1:nrow(cox_0.05)) {
    cox_0.05$gene_name[row] <- Japan_eqtl_0.05[Japan_eqtl_0.05$gene_id ==cox_0.05$gene[row],]$gene_name
    if (cox_0.05$simple_gene_id[row] %in% list$ensembl_gene_id){
    cox_0.05$gene_biotype[row] <- list[list$ensembl_gene_id ==cox_0.05$simple_gene_id[row],]$gene_biotype
    
  }
  
}
cox_0.05 <- cox_0.05[order(cox_0.05$P),]
write.xlsx(cox_0.05,"/Users/mstermo/eQTL-continue/All_supple_materials/Univariate-significant_eGenes.xlsx", sheetName = "Cox", 
           col.names = TRUE, row.names = FALSE, append = TRUE)

#########annotate the common prognositic eGenes
KM_cox_0.05_genes <- as.data.frame(KM_cox_0.05_genes)

colnames(KM_cox_0.05_genes)[1] <- "gene"
KM_cox_0.05_genes$simple_gene_id <- sapply(strsplit(KM_cox_0.05_genes$gene,'\\.'),'[',1)

list <- getBM(attributes=c('gene_biotype','ensembl_gene_id','external_gene_name'), 
              filters = 'ensembl_gene_id',
              values = KM_cox_0.05_genes$simple_gene_id, mart = ensembl)


for (row in 1:nrow(KM_cox_0.05_genes)) {
  
    KM_cox_0.05_genes$gene_name[row] <- Japan_eqtl_0.05[Japan_eqtl_0.05$gene_id ==KM_cox_0.05_genes$gene[row],]$gene_name
    KM_cox_0.05_genes$KM_p_value[row] <- KM_0.05[KM_0.05$simple_gene_id ==KM_cox_0.05_genes$simple_gene_id[row],]$p
    KM_cox_0.05_genes$Cox_p_value[row] <- cox_0.05[cox_0.05$simple_gene_id==KM_cox_0.05_genes$simple_gene_id[row],]$P
    if (KM_cox_0.05_genes$simple_gene_id[row] %in% list$ensembl_gene_id){
      KM_cox_0.05_genes$gene_biotype[row] <- list[list$ensembl_gene_id ==KM_cox_0.05_genes$simple_gene_id[row],]$gene_biotype
    
  }
  
}

KM_cox_0.05_genes <- KM_cox_0.05_genes[order(KM_cox_0.05_genes$KM_p_value),]
write.xlsx(KM_cox_0.05_genes,"/Users/mstermo/eQTL-continue/All_supple_materials/Univariate-significant_eGenes.xlsx", sheetName = "KM_Cox_common", 
           col.names = TRUE, row.names = FALSE, append = TRUE)


#########annotate signif_pairs_with_rs_id file
#Supp_2 <- fread("/Users/xiyas/eQTL-continue/Manuscript_Systematically_identification/All_supple_materials/Supplementary_table2_signif_pairs.xlsx")
Supp_2 <- fread("/Users/mstermo/Degree_Project/data/eQTL-ccRCC-final.signifpairs.txt")

#colnames(cox_0.05)[1] <- "gene"
Supp_2$simple_gene_id <- sapply(strsplit(Supp_2$gene_id,'\\.'),'[',1)

list <- getBM(attributes=c('gene_biotype','ensembl_gene_id','external_gene_name'), 
              filters = 'ensembl_gene_id',
              values = Supp_2$simple_gene_id, mart = ensembl)


for (row in 1:nrow(Supp_2)) {
  Supp_2$gene_name[row] <- Japan_eqtl[Japan_eqtl$gene_id ==Supp_2$gene_id[row],]$gene_name
  if (Supp_2$simple_gene_id[row] %in% list$ensembl_gene_id){
    Supp_2$gene_biotype[row] <- list[list$ensembl_gene_id ==Supp_2$simple_gene_id[row],]$gene_biotype
    
  }
  
}

write.xlsx(Supp_2,"/Users/mstermo/eQTL-continue/All_supple_materials/Supplementary_table2_signif_pairs.xlsx", sheetName = "signif_pairs", 
           col.names = TRUE, row.names = FALSE, append = TRUE)



####write Supple_1

Japan_eqtl <- Japan_eqtl[order(Japan_eqtl$qval),]
Japan_eqtl_0.05 <- Japan_eqtl_0.05[order(Japan_eqtl_0.05$qval),]
library(openxlsx)
dataset_names <- list("All_tested_genes" = Japan_eqtl, "Significant_eGenes" = Japan_eqtl_0.05)
openxlsx::write.xlsx(dataset_names,"/Users/mstermo/eQTL-continue/All_supple_materials/Supplementary_table1_Japan_eGenes.xlsx",
           colNames = TRUE, rowNames = FALSE, append = TRUE)


#####Annotate variants for the 4558 all eQTLS survival analysis
colnames(KM_p_0.05)[1] <- "variant_id"
KM_p_0.05 <- merge(x=KM_p_0.05,y=Supp_2,by="variant_id")
KM_p_0.05 <- KM_p_0.05[order(KM_p_0.05$p),]

colnames(cox_p_0.05)[1] <- "variant_id"
cox_p_0.05 <- merge(x=cox_p_0.05,y=Supp_2,by="variant_id")
cox_p_0.05 <- cox_p_0.05[order(cox_p_0.05$P),]


test <- intersect(KM_p_0.05$genotype,cox_p_0.05$target)
KM_COX_SIG_eqtls <- test
test <- as.data.frame(test)
colnames(test) <- "variant_id"
#can not achieve by this:
#test$KM_p_value[row] <- KM_p_0.05[KM_p_0.05$genotype ==test$variant_id[row],]$p
#test$Cox_p_value[row] <- cox_p_0.05[cox_p_0.05$target==test$variant_id[row],]$P
#test<- merge(x=test,y=Supp_2,by="variant_id")


dataset_names <- list("KM_eQTLs" = KM_p_0.05, "Cox_eQTLs" = cox_p_0.05,"KM_Cox_common_eQTLs" = test)
openxlsx::write.xlsx(dataset_names,"/Users/mstermo/eQTL-continue/All_supple_materials/Supplementary_table4_survival_eqtls.xlsx",
                     colNames = TRUE, rowNames = FALSE, append = TRUE)



