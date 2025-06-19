# 1. load the library =====================
# Load all the libraries in one go
libraries <- c("knitr", "dplyr", "tidyr", "purrr", "survival", "ggplot2", 
               "tibble", "ggsurvfit", "stringr", "clusterProfiler", "enrichplot", 
               "org.Hs.eg.db", "data.table", "gtsummary", "readxl", "openxlsx", 
               "igraph", "survivalAnalysis", "survminer", "patchwork", 
               "tidyverse", "biomaRt", "grid", "forestploter")

# Use lapply to load them all
lapply(libraries, require, character.only = TRUE)

#install.packages("forestploter")

setwd('/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis')
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")

# 2.JP Survival analysis on clinical variables ================================
JP_data <- Cheng_readFile("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/ccRCC_tpm_target_gene_exp.txt")
JP_data$nT <- as.numeric(str_extract(JP_data$Stage, "T[0-9]") %>% str_replace("T", ""))
JP_data$N <- as.numeric(str_extract(JP_data$Stage, "N[0-9a-zA-Z]+(?=M)") %>% str_replace("N", ""))
JP_data$M <- as.numeric(str_extract(JP_data$Stage, "M[0-9a-zA-Z]") %>% str_replace("M", ""))

LivingDays <- as.numeric(JP_data$LivingDays)
SurvInput <- as.data.frame(LivingDays)
SurvInput$EXP <- as.numeric(JP_data[, gene])
SurvInput$Gender <- as.factor(JP_data[, 'Gender'])
SurvInput$nT <- JP_data[, 'nT']
SurvInput$N <- JP_data[, 'N']
SurvInput$M <- JP_data[, 'M']
SurvInput$DeadInd <- JP_data$Status %in% c('dead')
SurvInput$Age <- as.numeric(JP_data[, 'Age'])
SurvInput$SurvObj <- with(SurvInput, Surv(LivingDays, DeadInd))
##Plot Forest Plot -------------
covariate_names <- list(
  Gender = "Gender",
  Age = "Age",
  nT = "Tumor Stage (nT)",
  N = "Lymph Node Involvement (N)",
  M = "Metastasis (M)"
)
forest_plot <- map(vars(Gender, Age,nT, N,M), function(by)
{
  analyse_multivariate(SurvInput,
                       vars(LivingDays, DeadInd),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="OS"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 10))

# show the plot
forest_plot
dev.off()
ggsave("FigureS2.pdf", plot = forest_plot, device = "pdf",width = 8, height = 4, units = "in")

# 3.eGene survival analysis both adjusted and unadjusted =================
## KM and Cox analysis performed by separate R scripts.
## eGene_survival_analysis.R
dim(KM_0.05)

# Test all genes : 432 genes which combined KM,COX,COX_adjusted
#union_result <- union(KM_0.05$gene,cox_0.05$gene) %>% as.data.frame()
#colnames(union_result) <- "gene"
# new: delete KM analysis

# 3.1 Figure 4A data Venn plot----------------------
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
  #ggsave(my_venn, file="my_venn.png",dpi = 600, width = 14, height =8,device = "png")
}

x <- list(
  A = cox_0.05$gene,
  B = cox_0.05_adjusted$gene
)
y <- list(
  A = sig_cox_results$genotype,
  B = sig_cox_results_adjusted$genotype
)
#png(filename = "Venn_diagram_4_section.png",width = 8, height = 5, res = 600)
display_venn(
  x,
  #filename = NULL,
  category.names = c("Unadjusted Cox Gene","Adjusted Cox Gene"),
  # Circles
  lwd = 2,
  lty = 'blank',
  #"#56B4E9", "#E69F00"
  #fill = c("#ffcc66","grey","#53c68c","#99d6ff"),
  fill = c("#ffce80","#56B4E9"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.1,0.1)
)
display_venn(
  y,
  #filename = NULL,
  category.names = c("Unadjusted Cox Gene","Adjusted Cox Gene"),
  # Circles
  lwd = 2,
  lty = 'blank',
  #"#56B4E9", "#E69F00"
  #fill = c("#ffcc66","grey","#53c68c","#99d6ff"),
  fill = c("#ffce80","#56B4E9"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.1,0.1)
)
#ggsave(filename = "Venn_final.png",dpi = 600, width = 12, height =8,device = "png")


# 3.2 Union_resuts dataset
fdr_eGene <- cox_0.05_annotated[cox_0.05_annotated$FDR < 0.05,]
fdr_eGene_adjusted <- cox_0.05_adjusted[cox_0.05_adjusted$FDR < 0.05,]
union_result <- union(cox_0.05_adjusted$gene,cox_0.05$gene)%>% as.data.frame()
colnames(union_result) <- "gene"

length(cox_0.05_adjusted$gene)
### 82 
length(union(cox_0.05_adjusted$gene,cox_0.05$gene))
### 187
dim(union_result)
#add information from original table
union_result <- merge(union_result,cox_0.05[,c("COX_coef","P","gene")],by = "gene",all.x=TRUE)
colnames(union_result)[2:3]<- c("COX_coef_univariate","P_univariate")

union_result <- merge(union_result,cox_0.05_adjusted[,c("COX_coef","P","gene")],by = "gene",all.x=TRUE)
colnames(union_result)[4:5]<- c("COX_coef_multivariate","P_multivariate")
#union_result <- merge(union_result,KM_0.05[,c("coef","gene")],by = "gene",all.x=TRUE)
#colnames(union_result)[4]<- "KM_binary_coef"

### Identify the rows where COX_coef and coef have different signs: rows_with_different_sign =================
check_signs <- function(row) {
  # Extract coefficients from the row and convert to numeric
  # coefs <- c(as.numeric(row["COX_coef_univariate"]), as.numeric(row["COX_coef_multivariate"]), as.numeric(row["KM_binary_coef"]))
  coefs <- c(as.numeric(row["COX_coef_univariate"]), as.numeric(row["COX_coef_multivariate"]))
  # Remove NAs
  coefs <- na.omit(coefs)
  # Check if there are any differences in signs
  return(length(unique(sign(coefs))) > 1)
}

# Apply the function to each row and add the result as a new column
union_result$different_signs <- apply(union_result, 1, check_signs)
# Filter the genes with different signs
genes_with_different_signs <- union_result %>% filter(different_signs == TRUE)
# Display the results
print(genes_with_different_signs)
#1 ENSG00000138395.14 
#2  ENSG00000237674.1  
# ___________-- after remove KM results there is no genes with different signs
#annotate the gene name 

union_result$simple_gene_id <- sapply(strsplit(union_result$gene,'\\.'),'[',1)
#ensembl <- useEnsembl(biomart = "ensembl", 
#                      dataset = "hsapiens_gene_ensembl", 
#                      mirror = "useast")
#list <- getBM(attributes=c('gene_biotype','ensembl_gene_id','external_gene_name'), 
##              filters = 'ensembl_gene_id',
#              values = union_result$simple_gene_id, mart = ensembl)

# left 187 eGenes after KM removal ,,,
union_result <- union_result %>%
  filter(!gene %in% genes_with_different_signs$gene)
dim(union_result)
#[1] 187  7

### How many fav/unfavorable:114 TRUE, 74 FALSE =================
# select the cox1 model coef, as representative, and for cox only variants, select the cox2 coef model
#union_result$overall_coef <- ifelse(!is.na(union_result$KM_binary_coef), 
#                                             union_result$KM_binary_coef,
#                                             ifelse(!is.na(union_result$COX_coef_univariate), 
#                                                    union_result$COX_coef_univariate, union_result$COX_coef_multivariate))

union_result$overall_coef <- ifelse(!is.na(union_result$COX_coef_univariate), 
                            union_result$COX_coef_univariate, union_result$COX_coef_multivariate)

# Determine the sign of the overall coefficient
#union_result$overall_sign <- sign(union_result$overall_coef)
# mark fav/unfav
union_result$marker <- ifelse(union_result$overall_coef > 0, "Unfavorable", "Favorable")
#Favorable Unfavorable 
#302         128 
table(union_result$marker)

#> table(union_result$marker)
#Favorable Unfavorable 
#105          82 

union_result <- union_result[order(union_result$overall_coef), ]
for (row in 1:nrow(union_result)) {
  if (union_result$simple_gene_id[row] %in% list$ensembl_gene_id){
    union_result$gene_biotype[row] <- list[list$ensembl_gene_id ==union_result$simple_gene_id[row],]$gene_biotype
    union_result$gene_name[row] <- list[list$ensembl_gene_id ==union_result$simple_gene_id[row],]$external_gene_name
  }
}
### Supp_table 6 new --------------------
### wrote union_result
union_result$gene_name <- signif_egenes$gene_name[match(union_result$gene, signif_egenes$gene_id)]
dim(union_result)
#[1] 187  11
colnames(union_result)
#[1] "gene"                  "COX_coef_univariate"   "P_univariate"          "COX_coef_multivariate" "P_multivariate"       
#[6] "different_signs"       "simple_gene_id"        "overall_coef"          "marker"                "gene_biotype"         
#[11] "gene_name"     
openxlsx::write.xlsx(union_result,"/Users/xiyas/eQTL-continue/Manuscript_Systematically_identification/eQTL_manuscript_1119/PLOS_Revise1/Supplementary_Materials/Supplementary_table5_Cox_significant_eGenes.xlsx", sheetName = "Cox_adjusted&unadjusted", 
     colNames = TRUE, rowNames = FALSE)

### GO term analysis: no significant overlaps =======================
eGenes_db <- readxl::read_excel("/Users/xiyas/eQTL-continue/Manuscript_Systematically_identification/eQTL_manuscript/All_Supple_Materials/Supplementary_table2_Japan_eGenes.xlsx",sheet = 'Significant_eGenes')
eGenes_db$simple_id <- sapply(strsplit(eGenes_db$gene_id,'\\.'),'[',1)
list <- getBM(attributes=c('gene_biotype','ensembl_gene_id','external_gene_name'), 
              filters = 'ensembl_gene_id',
              values = eGenes_db$simple_id, mart = ensembl)
for (row in 1:nrow(eGenes_db)) {
  if (eGenes_db$simple_id [row] %in% list$ensembl_gene_id){
    eGenes_db$gene_biotype[row] <- list[list$ensembl_gene_id ==eGenes_db$simple_id[row],]$gene_biotype
  }
}

union_result$gene_symbol <- eGenes_db$gene_name[match(union_result$gene,eGenes_db$gene_id)]
union_result$Ensembl <- sapply(strsplit(union_result$gene,'\\.'),'[',1)
favorable_eGene <- union_result %>% filter(overall_coef < 0 )
unfavorable_eGene <- union_result %>% filter(overall_coef > 0 )
go_results <- enrichGO(gene = favorable_eGene$Ensembl, 
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "ALL",  # Biological Process ontology,
                       maxGSSize = 5000,
                       minGSSize = 5,
                       pvalueCutoff  = 0.2,
                       readable = TRUE)
#par(mfrow = c(1,2))
#plot(sort(go_results@result$pvalue))
dotplot(go_results, showCategory = 20)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

# Open a text file for writing
file_path <- "favorable_gene_symbols.txt"
file_conn <- file(file_path, "w")

# Write gene symbols separated by commas to the file
gene_symbols <- paste(favorable_eGene$gene_symbol, collapse = ",")
writeLines(gene_symbols, file_conn)
# Close the file connection
close(file_conn)
# Print confirmation message
cat("Gene symbols have been successfully written to", file_path, "\n")

# 4.JP eQTL survival analysis ====================
vcf_data<-fread("/Users/xiyas/Degree_Project/data/df1_for_eqtl.vcf.gz",fill = TRUE, header =TRUE)

##4.1 Clean vcf genotype matrix ====================
vcf_data_final <- vcf_data[,c(3,10:109)]
vcf_data_final_a = names(vcf_data_final) %>%
  map(
    function(x) 
      vcf_data_final %>% 
      select(x) %>% 
      separate(x, 
               into = paste0(x, c("_genotype", "aa")), 
               sep = ":")  %>% 
      select(paste0(x, "_genotype"))
  ) %>%
  bind_cols()

## 4.2 Get the significant eQTLs list from JP eqtl analysis  ====================
signif_pairs <- fread("/Users/xiyas/Degree_Project/data/eQTL-ccRCC-final.signifpairs.txt")
###this is for all signif pairs km analsis
###4558 significant variants 
sig_variant_genotype_japan_4558 = vcf_data_final_a[vcf_data_final_a$ID_genotype %in% signif_pairs$variant_id,]
sig_variant_genotype_japan_4558[sig_variant_genotype_japan_4558 == '0/0'] <- '0'
sig_variant_genotype_japan_4558[sig_variant_genotype_japan_4558 == '0/1' | sig_variant_genotype_japan_4558 == '0|1'] <- '1'
sig_variant_genotype_japan_4558[sig_variant_genotype_japan_4558 == '1/1' | sig_variant_genotype_japan_4558 == '1|1' ] <- '2'
sig_variant_genotype_japan_4558[sig_variant_genotype_japan_4558 == './.'] <- '0'
sig_variant_genotype_japan_4558[sig_variant_genotype_japan_4558 == '0|0'] <- '0'

## column, row T transfer
sig_variant_genotype_japan <- t(sig_variant_genotype_japan_4558)
colnames(sig_variant_genotype_japan) <- sig_variant_genotype_japan[1,]
## colnames(sig_variant_genotype_japan) <- new_name
sig_variant_genotype_japan <- sig_variant_genotype_japan[-1,]
##change row names 
a <-rownames(sig_variant_genotype_japan)
a <-strsplit(a,split = '_')

a = sapply(a,"[[",2)
a = gsub("-tumor","",a)
a = gsub("\\-","_",a)
rownames(sig_variant_genotype_japan) <- a

## 4.3 Attention: any variants with genotype less than 3 indivdual should be removed  =================================
genotype_counts_list <- apply(sig_variant_genotype_japan, 1, function(x) as.data.frame(table(factor(x, levels = c(0, 1, 2)))))
genotype_counts_df <- do.call(rbind, lapply(seq_along(genotype_counts_list), function(i) {
  df <- genotype_counts_list[[i]]
  df$SNP <- rownames(sig_variant_genotype_japan)[i]
  df
}))
genotype_counts_wide <- genotype_counts_df %>%
  pivot_wider(names_from = Var1, values_from = Freq, values_fill = list(Freq = 0))
colnames(genotype_counts_wide) <- c("SNP", "Genotype_0", "Genotype_1", "Genotype_2")
# View the result
print(genotype_counts_wide)

# 61 SNPs with few genotype : Here I mean if any of the heterozygous/homozygous is larger than 3, 
#then it is being kept; but should be changed to need
# all genotypes is larger than 3? 

snps_with_few_genotype <- genotype_counts_wide %>%
  filter(Genotype_1 < 3 & Genotype_2 < 3 ) %>%
  pull(SNP)

## 4.4 KM for variant (Categorical) =================================

### Make combined_surv_genotype table =================================
setwd("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis")
source('Cheng_toolbox_beta.R')
cancerType<-"ccRCC_genotype_4558"
#combined_surv_genotype <- merge(surv_data,sig_variant_genotype_japan_4558,by = "row.names",all.x = FALSE)
#write.table(combined_surv_genotype,file= "ccRCC_genotype_4558_tpm_target_gene_exp.txt",row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")

combined_surv_genotype_test<-as.matrix(read.table("ccRCC_genotype_4558_tpm_target_gene_exp.txt",header=T,sep="\t"))
#sapply(exp,class)
# all character value
genotype_list<-as.matrix(colnames(combined_surv_genotype_test)[8:dim(combined_surv_genotype_test)[2]])
gnedataDir=path_raw
path_out_pdf<-paste0(path_raw,"KM_genotype_pdf_4558/",sep = "")
path_out_pdf
#dir.create(path_out_pdf)
setwd(path_out_pdf)

### Modify the new generate KM plot function for 3 kinds of genoytpes ----------------

Cheng_generateKMplot_genoptype <- function(SurvInput,outFile,figureTitle) { 
  SurvInput$EXP[SurvInput$EXP==0]=0
  SurvInput$EXP[SurvInput$EXP==1]=1
  SurvInput$EXP[SurvInput$EXP==2]=2
  
  survOut = survfit(SurvObj ~ EXP,SurvInput)
  pdf(file = paste(outFile,"pdf",sep = "."),width=5.8,height=6)
  plot(survOut, col=c("#008B45B2","#3C5488B2","#DC0000B2"), mark.time=T, cex=1.4,xlab="Time (year)",xscale = 365,lty =1, ylab = "Survival Probability",las=1, cex.lab=1.4)
  group1legend = paste("AA genotype",paste("(n =",as.character(sum(SurvInput$EXP==0)),")",sep = ""), sep = " ")
  group2legend = paste("Aa genotype",paste("(n =",as.character(sum(SurvInput$EXP==1)),")",sep = ""), sep = " ")
  group3legend = paste("aa genotype",paste("(n =",as.character(sum(SurvInput$EXP==2)),")",sep = ""), sep = " ")
  if (!hasArg("figureTitle")) {
    #title("ASS1")
  } else {
    title(figureTitle)
  }
  legend(
    "bottomleft",
    legend=c(group1legend,group2legend,group3legend),
    col=c("#008B45B2","#3C5488B2","#DC0000B2"),
    horiz=FALSE,
    lty=1,
    bty='n',
    cex = 1.4)
  #  title("ASS1") 
  res = survdiff(SurvObj ~ EXP, SurvInput)
  logRankP = 1 - pchisq(res$chisq, length(res$n)-1)
  legend("topright",legend =c(paste0("P=",as.character(format(logRankP,scientific = TRUE,digits = 3)))),text.font=2,bty="n",cex = 1.4)
  dev.off()
  sum_result<-summary(coxph(SurvObj ~ EXP, SurvInput))
  coef<-sum_result$coefficients[1]
  result = as.data.frame(logRankP)
  #result$EXPcut = cutoff
  result$coef=coef
  return(result)
}


### Perform KM variant analysis ====================
### for only one genotype
#genotype_list[3,1]

#output<-Cheng_generateSurvInputfromTCGA(genotype_list[1,1],cancerType,dataDir)
#res = survdiff(output$SurvObj ~ output$EXP, output)
#result<-Cheng_generateKMplot_genoptype(output,0.8,outFile=paste0(cancerType,"_",genotype_list[1,1]))
cancerType <- "ccRCC_genotype_4558"
#### Batch  processing KM variant analysis : KM_p_0.05====================
p_list<-NULL
cutoff_list<-NULL
coef_list<-NULL
for (j in 1:length(genotype_list)){
  print(j)
  output<-Cheng_generateSurvInputfromTCGA(genotype_list[j,1],cancerType,dataDir)
  setwd(path_out_pdf)
  result<-Cheng_generateKMplot_genoptype(output,0.8,outFile=paste0(cancerType,"_",genotype_list[j,1]))
  log_rank_p<-as.matrix(result$logRankP)
  cut_off<-0.8
  coef<-as.matrix(result$coef)
  p_list<-rbind(p_list,log_rank_p)
  cutoff_list<-rbind(cutoff_list,cut_off)
  coef_list<-rbind(coef_list,coef)
  rm(log_rank_p,cut_off,coef)
}
fdr_list <-p.adjust(p_list,method="BH",n=length(p_list))
final_genotype_result<-cbind(as.matrix(genotype_list),coef_list,cutoff_list,p_list,fdr_list)
colnames(final_genotype_result)<-c("genotype","coef","cutoff","p","fdr")

setwd("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis")
write.table(final_genotype_result,file="genotype_coef_cutoff_list.txt",sep="\t",row.names=F,col.names=T,quote=F)

###select the variants with fdr <0.05
final_genotype_result <- as.data.frame(final_genotype_result)
KM_p_0.05 <- filter(final_genotype_result,p <= 0.05)
KM_p_0.05 <- KM_p_0.05[order(KM_p_0.05$p), ]
# 253 variants
KM_p_0.05 <- KM_p_0.05 %>% filter(KM_p_0.05$p !=0)
# After REMOVE : 242 variants 
KM_p_0.05 <- KM_p_0.05 %>% filter(! genotype %in% snps_with_few_genotype)

# 4.5. Cox variant analysis: both categorical and numeric ====================

setwd("/Users/mstermo/eQTL-continue/eQTL-continue-survival-analysis")
#survival table didn't remove the lost patients from ccRCC_1 to ccRCC_106
surv_data_test <-  surv_data[-c(33,47,57,63,74,101),]

survival<-surv_data_test[,c("Status","LivingDays")]
survival[which(survival[,1]=="alive"),1]=0
survival[which(survival[,1]=="dead"),1]=1
survival[,2]=gsub(" ","",survival[,2])
i <- c(1,2)
survival[ , i] <- apply(survival[ , i], 2,            # Specify own function within apply
                        function(x) as.numeric(as.character(x)))

survival$T <- as.numeric(str_extract(surv_data_test$Stage, "T[0-9]") %>% str_replace("T", ""))
survival$N <- as.numeric(str_extract(surv_data_test$Stage, "N[0-9a-zA-Z]+(?=M)") %>% str_replace("N", ""))
survival$M <- as.numeric(str_extract(surv_data_test$Stage, "M[0-9a-zA-Z]") %>% str_replace("M", ""))

#check the survival metadata is the numeric number
sapply(survival,class)
coef_list<-NULL
p_list<-NULL

#sig_variant_genotype_japan <- t(sig_variant_genotype_japan)

#for (i in 1:dim(sig_variant_genotype_japan)[1]){
#  cox_result=summary(coxph(Surv(survival[,2], survival[,1]) ~ as.matrix(sig_variant_genotype_japan[i,])))
#  coef=as.matrix(cox_result$coefficients)[1,"coef"]
#  p=as.matrix(cox_result$coefficients)[1,"Pr(>|z|)"]
# coef_list<-rbind(coef_list,coef)
#  p_list=rbind(p_list,p)
#}

###4.5.1 Single-Cox analysis test on single variants: if adjusted for TNM stage ====================
cox_result_adjusted=summary(coxph(Surv(survival[,2], survival[,1]) ~
                           factor(sig_variant_genotype_japan['chr2_185802394_T_C',])+ 
                           survival$T +
                           survival$N + 
                          survival$M))

cox_result=summary(coxph(Surv(survival[,2], survival[,1]) ~
                           factor(sig_variant_genotype_japan['chr10_80032436_A_G',])))
                         
cox_result=summary(coxph(Surv(survival[,2], survival[,1]) ~
                           as.numeric(sig_variant_genotype_japan['chr10_80032436_A_G',])))

#coef=as.matrix(cox_result$coefficients)[1,"coef"]
#p=as.matrix(cox_result$coefficients)[1,"Pr(>|z|)"]

### 4.5.2 176 is old eQTLs from previous results, not used ------------------
# if only categorical variables and first level of P value is in cox_p_0.05, and also this 176 variants result didn't remove the low genotypes counts
#setdiff(cox_p_0.05$genotype,sig_cox_results[sig_cox_results$P_categorical_het <0.05,]$genotype)
#[1] "chr1_145888034_T_G" "chr1_150260371_T_C" "chr1_157747046_T_G" "chr1_157767141_C_T" "chr1_157778492_C_T" "chr1_157797616_T_G" "chr1_157798395_T_C"
#[8] "chr1_157801238_T_C" "chr2_186505386_A_C" "chr7_139631031_T_C" "chr19_49399658_A_C" "chrX_135465460_T_G"

fdr_list<-p.adjust(p_list,method="BH",n=length(p_list))
cox_list<-cbind(coef_list,p_list,fdr_list)
colnames(cox_list)<-c("COX_coef","P","FDR")
target <- rownames(sig_variant_genotype_japan)
target<-cbind(target,cox_list)
#count the significant genes
#class(target)
target <- as.data.frame(target)
cox_p_0.05_new<- filter(target,P <= 0.05)

########cox_p_0.05
cox_0.05 <- filter(target,P <= 0.05)
cox_0.001 <- filter(target, P <= 0.001)
#write?
write.table(target,file="target_one_cox_in_TCGA.txt",sep="\t",row.names=F,col.names=T,quote=F)


## 4.5.3 NOT ADJUSTED Get sig_cox_results without adjust tumor grade: A new cox code which test both categorical and numeric value in SNP ========================

# Initialize lists to store results
coef_list_additive <- c()
p_list_additive <- c()
coef_list_categorical_het <- c()
coef_list_categorical_hom <- c()
p_list_categorical_het <- c()
p_list_categorical_hom <- c()
lrt_list <- c()
wald_list <- c()
logrank_list <- c()

# Loop over each variant
for (i in 1:dim(sig_variant_genotype_japan)[1]) {
  # Additive model (numeric)
  additive_model <- coxph(Surv(survival[, 2], survival[, 1]) ~ as.numeric(sig_variant_genotype_japan[i, ]))
  additive_result <- summary(additive_model)
  coef_additive <- as.matrix(additive_result$coefficients)[1, "coef"]
  p_additive <- as.matrix(additive_result$coefficients)[1, "Pr(>|z|)"]
  coef_list_additive <- rbind(coef_list_additive, coef_additive)
  p_list_additive <- rbind(p_list_additive, p_additive)
  # Categorical model (factor)
  categorical_model <- tryCatch({
    coxph(Surv(survival[, 2], survival[, 1]) ~ factor(sig_variant_genotype_japan[i, ]))
  }, error = function(e) {
    NULL  # Return NULL if only one level is present
  })
  if (!is.null(categorical_model)) {
    categorical_result <- summary(categorical_model)
    # Extract coefficients and p-values for heterozygous and homozygous minor genotypes
    if ("factor(sig_variant_genotype_japan[i, ])1" %in% rownames(categorical_result$coefficients)) {
      coef_categorical_het <- as.matrix(categorical_result$coefficients)["factor(sig_variant_genotype_japan[i, ])1", "coef"]
      p_categorical_het <- as.matrix(categorical_result$coefficients)["factor(sig_variant_genotype_japan[i, ])1", "Pr(>|z|)"]
    } else {
      coef_categorical_het <- NA
      p_categorical_het <- NA
    }
    if ("factor(sig_variant_genotype_japan[i, ])2" %in% rownames(categorical_result$coefficients)) {
      coef_categorical_hom <- as.matrix(categorical_result$coefficients)["factor(sig_variant_genotype_japan[i, ])2", "coef"]
      p_categorical_hom <- as.matrix(categorical_result$coefficients)["factor(sig_variant_genotype_japan[i, ])2", "Pr(>|z|)"]
    } else {
      coef_categorical_hom <- NA
      p_categorical_hom <- NA
    }
    # Extract LRT, Wald, and Logrank test p-values
    lrt_pvalue <- categorical_result$logtest["pvalue"]
    wald_pvalue <- categorical_result$waldtest["pvalue"]
    logrank_pvalue <- categorical_result$sctest["pvalue"]
    coef_list_categorical_het <- rbind(coef_list_categorical_het, coef_categorical_het)
    coef_list_categorical_hom <- rbind(coef_list_categorical_hom, coef_categorical_hom)
    p_list_categorical_het <- rbind(p_list_categorical_het, p_categorical_het)
    p_list_categorical_hom <- rbind(p_list_categorical_hom, p_categorical_hom)
    lrt_list <- rbind(lrt_list, lrt_pvalue)
    wald_list <- rbind(wald_list, wald_pvalue)
    logrank_list <- rbind(logrank_list, logrank_pvalue)
  } else {
    # If only one level is present, assign NA values for coefficients and p-values
    coef_list_categorical_het <- rbind(coef_list_categorical_het, NA)
    coef_list_categorical_hom <- rbind(coef_list_categorical_hom, NA)
    p_list_categorical_het <- rbind(p_list_categorical_het, NA)
    p_list_categorical_hom <- rbind(p_list_categorical_hom, NA)
    lrt_list <- rbind(lrt_list, NA)
    wald_list <- rbind(wald_list, NA)
    logrank_list <- rbind(logrank_list, NA)
  }
}
#### Combine results into cox_table =======================
cox_table <- data.frame(
  genotype= rownames(sig_variant_genotype_japan),
  COX_coef_additive = coef_list_additive,
  P_additive = p_list_additive,
  COX_coef_categorical_het = coef_list_categorical_het,
  P_categorical_het = p_list_categorical_het,
  COX_coef_categorical_hom = coef_list_categorical_hom,
  P_categorical_hom = p_list_categorical_hom,
  LRT_pvalue = lrt_list,
  Wald_pvalue = wald_list,
  Logrank_pvalue = logrank_list
)
colnames(cox_table)[c(8:10)] = c('LRT_pvalue','Wald_pvalue','Logrank_pvalue')

#### Remove few genotype less than 3 individulas  ========================
##### 4558 left 4497 snps
cox_table <- cox_table %>% filter(!genotype %in% snps_with_few_genotype)
dim(cox_table)
#[1] 4497   10
## make a plot by plot_3_p_value
#? so which p-value in the whole model test should used? > wald test

# View the result table
#cox_table$cox_table$P_additive < 0.05) 

# Plot using ggplot2
cox_table_long <- cox_table %>%
  pivot_longer(cols = c(LRT_pvalue, Wald_pvalue, Logrank_pvalue), 
               names_to = "Test", 
               values_to = "P_value")
ggplot(cox_table_long, aes(x = genotype, y = -log10(P_value), color = Test)) +
  geom_point() +
  labs(title = "Comparison of P-values from LRT, Wald, and Logrank Tests",
       x = "genotype",
       y = "-log10(P-value)",
       color = "Test") +
  theme_minimal() 
#### Filter significant results: if take wald test as whole test  ==================
sig_cox_results <- cox_table[(!is.na(cox_table$P_additive) & cox_table$P_additive < 0.05) | 
                                      (!is.na(cox_table$P_categorical_het) & cox_table$P_categorical_het < 0.05) | 
                                      (!is.na(cox_table$P_categorical_hom) & cox_table$P_categorical_hom < 0.05) |
                                     (!is.na(cox_table$Wald_pvalue) & cox_table$Wald_pvalue< 0.05),]

# 367 sig_cox_results
dim(sig_cox_results)
#sig_cox_results <- cox_table[(!is.na(cox_table$Wald_pvalue) & cox_table$Wald_pvalue< 0.05),]
#sig_cox_results <- cox_table[(!is.na(cox_table$P_additive) & cox_table$P_additive< 0.05) |
                                      # (!is.na(cox_table$Wald_pvalue) & cox_table$Wald_pvalue< 0.05),]
# Count significant SNPs in each model
num_significant_additive <- sum(!is.na(sig_cox_results$P_additive) & sig_cox_results$P_additive < 0.05)
num_significant_categorical_het <- sum(!is.na(sig_cox_results$P_categorical_het) & sig_cox_results$P_categorical_het < 0.05)
num_significant_categorical_hom <- sum(!is.na(sig_cox_results$P_categorical_hom) & sig_cox_results$P_categorical_hom < 0.05)

# 4.5.4 ADJUSTED Get adjusted sig_cox_results which is included tumor grade as adjusted variables ========================

# Initialize lists to store results
coef_list_additive <- c()
p_list_additive <- c()
coef_list_categorical_het <- c()
coef_list_categorical_hom <- c()
p_list_categorical_het <- c()
p_list_categorical_hom <- c()
lrt_list <- c()
wald_list <- c()
logrank_list <- c()

# Loop over each variant
for (i in 1:dim(sig_variant_genotype_japan)[1]) {
  # Additive model (numeric)
  additive_model <- coxph(Surv(survival[, 2], survival[, 1]) ~ as.numeric(sig_variant_genotype_japan[i, ]) + survival$T + survival$N + survival$M)
  additive_result <- summary(additive_model)
  coef_additive <- as.matrix(additive_result$coefficients)[1, "coef"]
  p_additive <- as.matrix(additive_result$coefficients)[1, "Pr(>|z|)"]
  coef_list_additive <- rbind(coef_list_additive, coef_additive)
  p_list_additive <- rbind(p_list_additive, p_additive)
  
  # Categorical model (factor)
  categorical_model <- tryCatch({
    coxph(Surv(survival[, 2], survival[, 1]) ~ factor(sig_variant_genotype_japan[i, ]) + survival$T + survival$N + survival$M)
  }, error = function(e) {
    NULL  # Return NULL if only one level is present
  })
  
  if (!is.null(categorical_model)) {
    categorical_result <- summary(categorical_model)
    # Extract coefficients and p-values for heterozygous and homozygous minor genotypes
    if ("factor(sig_variant_genotype_japan[i, ])1" %in% rownames(categorical_result$coefficients)) {
      coef_categorical_het <- as.matrix(categorical_result$coefficients)["factor(sig_variant_genotype_japan[i, ])1", "coef"]
      p_categorical_het <- as.matrix(categorical_result$coefficients)["factor(sig_variant_genotype_japan[i, ])1", "Pr(>|z|)"]
    } else {
      coef_categorical_het <- NA
      p_categorical_het <- NA
    }
    
    if ("factor(sig_variant_genotype_japan[i, ])2" %in% rownames(categorical_result$coefficients)) {
      coef_categorical_hom <- as.matrix(categorical_result$coefficients)["factor(sig_variant_genotype_japan[i, ])2", "coef"]
      p_categorical_hom <- as.matrix(categorical_result$coefficients)["factor(sig_variant_genotype_japan[i, ])2", "Pr(>|z|)"]
    } else {
      coef_categorical_hom <- NA
      p_categorical_hom <- NA
    }
    
    # Extract LRT, Wald, and Logrank test p-values
    lrt_pvalue <- categorical_result$logtest["pvalue"]
    wald_pvalue <- categorical_result$waldtest["pvalue"]
    logrank_pvalue <- categorical_result$sctest["pvalue"]
    
    coef_list_categorical_het <- rbind(coef_list_categorical_het, coef_categorical_het)
    coef_list_categorical_hom <- rbind(coef_list_categorical_hom, coef_categorical_hom)
    p_list_categorical_het <- rbind(p_list_categorical_het, p_categorical_het)
    p_list_categorical_hom <- rbind(p_list_categorical_hom, p_categorical_hom)
    lrt_list <- rbind(lrt_list, lrt_pvalue)
    wald_list <- rbind(wald_list, wald_pvalue)
    logrank_list <- rbind(logrank_list, logrank_pvalue)
  } else {
    # If only one level is present, assign NA values for coefficients and p-values
    coef_list_categorical_het <- rbind(coef_list_categorical_het, NA)
    coef_list_categorical_hom <- rbind(coef_list_categorical_hom, NA)
    p_list_categorical_het <- rbind(p_list_categorical_het, NA)
    p_list_categorical_hom <- rbind(p_list_categorical_hom, NA)
    lrt_list <- rbind(lrt_list, NA)
    wald_list <- rbind(wald_list, NA)
    logrank_list <- rbind(logrank_list, NA)
  }
}

#### Combine results into "cox_table_adjusted" -----------
cox_table_adjusted <- data.frame(
  genotype = rownames(sig_variant_genotype_japan),
  COX_coef_additive = coef_list_additive,
  P_additive = p_list_additive,
  COX_coef_categorical_het = coef_list_categorical_het,
  P_categorical_het = p_list_categorical_het,
  COX_coef_categorical_hom = coef_list_categorical_hom,
  P_categorical_hom = p_list_categorical_hom,
  LRT_pvalue = lrt_list,
  Wald_pvalue = wald_list,
  Logrank_pvalue = logrank_list
)
## need to rename 
colnames(cox_table_adjusted)[c(8:10)] = c('LRT_pvalue','Wald_pvalue','Logrank_pvalue')

##### 4558 left 4497 snps
cox_table_adjusted <- cox_table_adjusted %>% filter(!genotype %in% snps_with_few_genotype)

#### get "sig_cox_results_adjusted" -----------
### in this case, Wald_pvalue is showing all significant results, cause the tumor grade is extremely significant. 
sig_cox_results_adjusted <- cox_table_adjusted[(!is.na(cox_table_adjusted$P_additive) & cox_table_adjusted$P_additive < 0.05) | 
                               (!is.na(cox_table_adjusted$P_categorical_het) & cox_table_adjusted$P_categorical_het < 0.05) | 
                               (!is.na(cox_table_adjusted$P_categorical_hom) & cox_table_adjusted$P_categorical_hom < 0.05),]



#  531 sig_cox_results_adjusted 
#dim(sig_cox_results)
dim(sig_cox_results_adjusted)
# compare the results
num_significant_additive <- sum(!is.na(sig_cox_results$P_additive) & sig_cox_results$P_additive < 0.05)
num_significant_categorical_het <- sum(!is.na(sig_cox_results$P_categorical_het) & sig_cox_results$P_categorical_het < 0.05)
num_significant_categorical_hom <- sum(!is.na(sig_cox_results$P_categorical_hom) & sig_cox_results$P_categorical_hom < 0.05)

num_significant_additive <- sum(!is.na(sig_cox_results_adjusted$P_additive) & sig_cox_results_adjusted$P_additive < 0.05)
num_significant_categorical_het <- sum(!is.na(sig_cox_results_adjusted$P_categorical_het) & sig_cox_results_adjusted$P_categorical_het < 0.05)
num_significant_categorical_hom <- sum(!is.na(sig_cox_results_adjusted$P_categorical_hom) & sig_cox_results_adjusted$P_categorical_hom < 0.05)

# conclusion: seems like if adjusted, can get more significant results.
# see Unadjusted_vs_adjusted.R

# 5.NEW Get the variants either significant in KM or Cox all models ======================
# Updates: no KM analysis anymore
# union_variant_result <- union(KM_p_0.05$genotype,sig_cox_results$genotype) 
union_variant_result <- union(sig_cox_results$genotype,sig_cox_results_adjusted$genotype) %>% as.data.frame()
colnames(union_variant_result) <- "genotype"
#colnames(cox_p_0.05)[1] <- "genotype"
union_variant_result <- merge(union_variant_result,sig_cox_results,by = "genotype",all.x=TRUE)
# allow to have different signs currently
#variants_with_different_signs <- union_variant_result %>% filter(different_signs == TRUE)
# Merge with sig_cox_results_adjusted
union_variant_result <- merge(union_variant_result, sig_cox_results_adjusted, by = "genotype", all.x = TRUE)
colnames(union_variant_result)
dim(union_variant_result)
#> dim(union_variant_result)
#[1] 711  19
# Rename the columns from sig_cox_results_adjusted to indicate adjusted values
adjusted_columns <- colnames(sig_cox_results_adjusted)[-1] # exclude the genotype column
new_adjusted_columns <- paste0(adjusted_columns, "_adjusted")

# Rename columns in the union_variant_result
colnames(union_variant_result)[(ncol(union_variant_result) - length(adjusted_columns) + 1):ncol(union_variant_result)] <- new_adjusted_columns
colnames(union_variant_result) <- sub("\\.x$", "", colnames(union_variant_result))
colnames(union_variant_result)
## how many unfavorable:114 (genes) and 215 (variants)
#table(union_variant_result$COX_coef > 0 | union_variant_result$coef > 0)
#union_variant_result <- merge(union_variant_result,KM_p_0.05[,c("coef","genotype","p")],by = "genotype",all.x=TRUE)
print(colnames(union_variant_result))
#colnames(union_variant_result)[(ncol(union_variant_result) - 1):ncol(union_variant_result)] <- c('KM_coef', 'KM_P')

# 6.Check the direction of SNP allelic affect ================================
# Define a function to check if a variant has different signs in its coefficients
check_signs <- function(row) {
  # Extract coefficients from the row and convert to numeric
  coefs <- c(
    as.numeric(row["COX_coef_additive"]),
    as.numeric(row["COX_coef_categorical_het"]),
    as.numeric(row["COX_coef_categorical_hom"]),
    #as.numeric(row["KM_coef"]),
    as.numeric(row["COX_coef_additive_adjusted"]),
    as.numeric(row["COX_coef_categorical_het_adjusted"]),
    as.numeric(row["COX_coef_categorical_hom_adjusted"])
  )
  # Remove NAs
  coefs <- na.omit(coefs)
  # Check if there are any differences in signs
  return(length(unique(sign(coefs))) > 1)
}

# Apply the function to each row and add the result as a new column
union_variant_result$different_signs <- apply(union_variant_result, 1, check_signs)
table(union_variant_result$different_signs)
# Step 1: Filter "sure_ones" where `different_signs` is FALSE

sure_ones <- union_variant_result[union_variant_result$different_signs == FALSE, ]

# Step 2: Calculate the `overall_coef`
# Using the same logic you applied earlier for genes, choose between different coefficients
sure_ones <- sure_ones %>%
  mutate(overall_coef = coalesce(
    COX_coef_additive,
    COX_coef_categorical_het,
    COX_coef_categorical_hom,
    COX_coef_additive_adjusted,
    COX_coef_categorical_het_adjusted,
    COX_coef_categorical_hom_adjusted
  ))

# Step 3: Determine if it's "Favorable" or "Unfavorable"
sure_ones <- sure_ones %>%
  mutate(marker = ifelse(overall_coef > 0, "Unfavorable", "Favorable"))


# You can now check the results, e.g.:
table(sure_ones$marker)

# after removing KM variant analysis: 526 FALSE 185 TRUE
# > length(unique(union_variant_result$genotype))
# [1] 711
#> sum(!is.na(union_variant_result$COX_coef_categorical_het) & union_variant_result$P_categorical_het <0.05)
#[1] 164 These are previous cox_p_0.05 results after removing genotypes filters

## supplementary table 7 ----------------------------
## adding rsID
#ensembl_snp <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
SNP_lookup <- read.table("/Users/xiyas/Degree_Project/data/snp_lookup.txt",header = TRUE)

#test <- read.table("/Users/mstermo/Degree_Project/data/eQTL-ccRCC-final.genes.annotated.txt",header = TRUE)
##write.table(Japan_eqtl_final_0.05, file = "/Users/mstermo/eQTL-continue/All_supple_materials/Supple_table1", quote= FALSE)

for (row in 1:nrow(union_variant_result)) {
  K=strsplit(union_variant_result$genotype[row],split = "\\_")
  if (nchar(K[[1]][3]) == nchar(K[[1]][4])) {
    union_variant_result$variant_type[row] <- "SNP"
    }
  else if (nchar(K[[1]][3])> nchar(K[[1]][4])) {
    union_variant_result$variant_type[row] <- "Del"
  }
  else if (nchar(K[[1]][3])< nchar(K[[1]][4])) {
    union_variant_result$variant_type[row] <- "Ins"
  }
  else {
    union_variant_result$variant_type[row] <- "Other"
  }
}

for (row in 1:nrow(union_variant_result)) {
  if (union_variant_result$genotype[row] %in% SNP_lookup$variant_id){
    union_variant_result$rs_id[row] <- 
      SNP_lookup[SNP_lookup$variant_id ==union_variant_result$genotype[row] ,]$rs_id_dbsnp
  }
}

# modify column order
# Reorder columns by specifying the desired column order
union_variant_result <- union_variant_result[, c("genotype", "variant_type", "rs_id", 
                                                 "COX_coef_additive", "P_additive", 
                                                 "COX_coef_categorical_het", "P_categorical_het", 
                                                 "COX_coef_categorical_hom", "P_categorical_hom", 
                                                 "LRT_pvalue", "Wald_pvalue", "Logrank_pvalue", 
                                                 "COX_coef_additive_adjusted", "P_additive_adjusted", 
                                                 "COX_coef_categorical_het_adjusted", "P_categorical_het_adjusted", 
                                                 "COX_coef_categorical_hom_adjusted", "P_categorical_hom_adjusted", 
                                                 "LRT_pvalue_adjusted", "Wald_pvalue_adjusted", 
                                                 "Logrank_pvalue_adjusted", "different_signs")]

dim(union_variant_result)
#711 22
## supplementary table 7 
#penxlsx::write.xlsx(union_variant_result,"/Users/xiyas/eQTL-continue/Manuscript_Systematically_identification/eQTL_manuscript/All_Supple_Materials/Supplementary_table6_Cox_significant_eQTLs.xlsx", sheetName = "Cox_adjusted_unadjusted", 
#                   colNames = TRUE, rowNames = FALSE)


## get the Prognostic genes identified by KM/Cox models ------------------------
## fav.unfav genes and eQTLs
#KM_0.05_annotated <- read_excel("/Users/xiyas/eQTL-continue/Manuscript_Systematically_identification/eQTL_manuscript/All_supple_materials_new/Supplementary_table3_Univariate-significant_eGenes.xlsx", sheet = 'KM+Cox_model_1')
#cox_0.05_annotated <-read_excel("/Users/xiyas/eQTL-continue/Manuscript_Systematically_identification/eQTL_manuscript/All_supple_materials_new/Supplementary_table3_Univariate-significant_eGenes.xlsx", sheet = 'Cox_model_2')


# 7. JP Find final consensus eqtl-egene pairs ================================
# both significant on gene and eqtl level
# union_variant_results: 388 variants >>>>>>>>>>>>>>>>.???????????????????????????
# get the target_pairs, that is the variants that is the eqtl for the eGenes and meanwhile have direct prognostic effect on JP

target_pairs <- signif_pairs %>% filter(variant_id %in% union_variant_result$genotype) %>% 
  filter(gene_id %in% union_result$gene)

# merge
target_pairs <- merge(target_pairs,union_result,by.x = "gene_id", by.y = "gene",all.x = TRUE)
colnames(target_pairs)[2]<- "genotype"

#colnames(union_variant_result)[4] <- "genotype_overall_coef"
target_pairs <- merge(target_pairs,union_variant_result,by='genotype',all.x = TRUE)
dim(target_pairs)
#[1] 356  43
# ===356 pairs now,start to further slopes filtering
#union_variant_result$overall_coef <- ifelse(is.na(union_variant_result$genotype_coef), union_variant_result$genotype_COX_coef,
#                                            union_variant_result$genotype_coef)

#target_pairs$genotype_overall_coef <- union_variant_result$overall_coef[match(target_pairs$variant_id, union_variant_result$genotype)]
# overall_coef : this is the coef for eGenes on two models

## 7.1 Testing slopes ================================

# Ensure the coefficients are numeric
target_pairs$overall_coef <- as.numeric(target_pairs$overall_coef)
target_pairs$COX_coef_additive <- as.numeric(target_pairs$COX_coef_additive)
target_pairs$COX_coef_categorical_het <- as.numeric(target_pairs$COX_coef_categorical_het)
target_pairs$COX_coef_categorical_hom <- as.numeric(target_pairs$COX_coef_categorical_hom)
#target_pairs$KM_coef <- as.numeric(target_pairs$KM_coef)
#target_pairs$genotype_overall_coef <- as.numeric(target_pairs$genotype_overall_coef)
target_pairs$COX_coef_additive_adjusted <- as.numeric(target_pairs$COX_coef_additive_adjusted)
target_pairs$COX_coef_categorical_het_adjusted <- as.numeric(target_pairs$COX_coef_categorical_het_adjusted)
target_pairs$COX_coef_categorical_hom_adjusted <- as.numeric(target_pairs$COX_coef_categorical_hom_adjusted)

# Rename
colnames(target_pairs)[which(colnames(target_pairs) == "different_signs.y")] <- "genotype_different_signs"
colnames(target_pairs)[which(colnames(target_pairs) == "different_signs.x")] <- "gene_different_signs"

# 7.2 Define the filtering condition ==========================
filtered_target_pairs <- target_pairs %>%
  filter(
    (slope > 0 & ((sign(overall_coef) == sign(COX_coef_additive)) |
                    (sign(overall_coef) == sign(COX_coef_categorical_het)) |
                    (sign(overall_coef) == sign(COX_coef_categorical_hom)) |
                    #(sign(overall_coef) == sign(KM_coef)) |
                    (sign(overall_coef) == sign(COX_coef_additive_adjusted)) |
                    (sign(overall_coef) == sign(COX_coef_categorical_het_adjusted)) |
                    (sign(overall_coef) == sign(COX_coef_categorical_hom_adjusted))
    )) |
      (slope < 0 & ((sign(overall_coef) != sign(COX_coef_additive)) |
                      (sign(overall_coef) != sign(COX_coef_categorical_het)) |
                      (sign(overall_coef) != sign(COX_coef_categorical_hom)) |
                      #(sign(overall_coef) != sign(KM_coef)) |
                      (sign(overall_coef) != sign(COX_coef_additive_adjusted)) |
                      (sign(overall_coef) != sign(COX_coef_categorical_het_adjusted)) |
                      (sign(overall_coef) != sign(COX_coef_categorical_hom_adjusted))
      ))
  )


# View the filtered results
print(filtered_target_pairs)
dim(filtered_target_pairs)
#329 43

length(unique(filtered_target_pairs$genotype))
length(unique(filtered_target_pairs$gene_id))
length(unique(filtered_target_pairs$gene_name))
#223 and 54
# 534 pairs 2024.6 and 223 target snps,94 target eGenes
# after removing KM results: 54 eggenes with 223 eqtls within 329 pairs

### supp_table 8 newest =================
filtered_target_pairs$gene_name <- eGenes_db$gene_name[match(filtered_target_pairs$gene_id,eGenes_db$gene_id)]
openxlsx::write.xlsx(filtered_target_pairs,"/Users/xiyas/eQTL-continue/Manuscript_Systematically_identification/eQTL_manuscript/All_Supple_Materials/Supplementary_table8_filtered_target_pairs.xlsx", sheetName = "target_pairs", 
               colNames = TRUE, rowNames = FALSE, append = TRUE)

## GO term analysis in eGenes in filtered_target_pairs ==================

target_gene <- as.data.frame(unique(sapply(strsplit(filtered_target_pairs$gene_id,'\\.'),'[',1)))
target_gene$somatically <- target_gene[,1]%in% gene_count$Gene_ID
#target_gene<- unique(filtered_target_pairs$gene_name)
#file_path <- "/Users/xiyas/eQTL-continue/Manuscript_Systematically_identification/eQTL_manuscript/All_supple_materials_new/gene_names.txt"
# Write the gene names to the text file
#writeLines(eGenes_db$gene_name, file_path)

#favorable_eGene <- union_result %>% filter(COX_coef < 0 | coef < 0)
#unfavorable_eGene <- union_result %>% filter(COX_coef > 0 | coef > 0)

#ordered <- filtered_target_pairs[order(filtered_target_pairs$overall_coef), ]
#filtered_target_pairs$gene_name <- factor(filtered_target_pairs$gene_name, levels = unique(filtered_target_pairs$gene_name[order(filtered_target_pairs$overall_coef)]))
filtered_target_pairs <- filtered_target_pairs[order(filtered_target_pairs$overall_coef), ]

union_result_new <- union_result %>% filter(gene %in% unique(filtered_target_pairs$gene_id))
union_result_new$overall_coef <- as.numeric(union_result_new$overall_coef )
union_result_new <- union_result_new[order(union_result_new$overall_coef), ]
union_result_new$gene_name <- union_result_new$gene_symbol
#union_result_new$gene_name <- factor(union_result_new$gene_name,levels = union_result_new$gene_name[order(union_result_new$overall_coef)])
union_result_new$gene_name <- factor(union_result_new$gene_name, levels = unique(union_result_new$gene_name[order(union_result_new$overall_coef)]))
# Ensure that the dataframe is ordered by overall_coef
union_result_new <- union_result_new %>% arrange(overall_coef)
# Reorder the gene_name factor based on the ordered overall_coef
#union_result$gene_symbol <- eGenes_db$gene_name[match(union_result$gene,eGenes_db$gene_id)]
union_result_new$gene_symbol <- factor(union_result_new$gene_name, levels = union_result_new$gene_name)

## 7.2 Figure 4C Plot final prognostic eGene's Coefficients:  =================================
ggplot(union_result_new, aes(x = gene_symbol, y = overall_coef, fill = overall_coef > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#00468b", "#ad002a"), labels = c("Favorable", "Unfavorable")) +
  labs(title = "Coefficients of Genes",
       x = "Gene Name",
       y = "Coefficient Value",
       fill = "eGene prognosis function") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

biotypes <- unique(union_result_new$gene_biotype)

# Generate a palette with enough colors
biotype_palette <- pal_npg()(9)
library(ggplot2)
library(ggnewscale)

# Bar plot with the first fill scale for favorable/unfavorable
plot<- ggplot(union_result_new, aes(x = gene_symbol, y = overall_coef, fill = marker)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#00468b", "#ad002a"), labels = c("Favorable", "Unfavorable")) +
  # Add the second fill scale for gene_biotype using ggnewscale
  ggnewscale::new_scale_fill() +
  # Tile plot with the second fill scale for gene_biotype
  geom_tile(aes(x = gene_symbol, y = -8, fill = gene_biotype), width = 0.8, height = 1) +
  scale_fill_manual(values = biotype_palette, name = "Gene Biotype") +
  scale_y_continuous(breaks = seq(-7, 5, by = 1)) +
  labs(title = "Coefficients of Genes", x = "Gene Name", y = "Coefficient Value") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
        axis.text.y = element_text(size = 12),                       
        axis.title.x = element_text(size = 16),                    
        axis.title.y = element_text(size = 16),                    
        plot.title = element_text(size = 18, face = "bold"),     
        legend.text = element_text(size = 12),                    
        legend.title = element_text(size = 14),                       
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") 
plot

ggsave("Figure_4C-fav_unfav_target_pair_genes.pdf", plot = plot, width = 14.5, height = 8,dpi = 600,units = "in")

# Assuming union_result_new contains a column 'gene_biotype'

# Create a frequency table for the gene_biotype column
library(ggsci)  # For pal_npg()
library(ggplot2)

# Assuming union_result_new contains a column 'gene_biotype'

# Create a frequency table for the gene_biotype column
biotype_freq <- table(union_result_new$gene_biotype)

# Define the biotype palette using pal_npg (9 colors)
biotype_palette <- pal_npg()(length(unique(union_result_new$gene_biotype)))

# Create the pie chart
pie(biotype_freq, 
    col = biotype_palette,  # Use the same palette as geom_tile
    main = "Distribution of Gene Biotypes", 
    labels = paste(names(biotype_freq), round(100 * biotype_freq / sum(biotype_freq), 1), "%"))


### Create the pie chart Figure 4C --------------
pie(biotype_freq, 
    col = biotype_palette,  # Use the same palette as geom_tile
    main = "Distribution of Gene Biotypes", labels = NA)

# 8. Finally, the survival analysis on TCGA ========================

## 8.1.Capture the target variants from joint-calling VCF,by python file "joint-calling-target-capture-and-new-id.ipynb" =================
#  writeLines(unique(filtered_target_pairs$genotype), '/Users/xiyas/eQTL-continue/joint-calling/target_snp_202406.txt')
# 386 eQTLs,???

## 8.2 Read target vcf data and clean/fomat the data ========================

vcf_data<-fread("/Users/xiyas/eQTL-continue/joint-calling/capture-targets.vcf.gz",fill = TRUE, header =TRUE,sep = "\t",skip = '#CHROM')
vcf_data_test <- vcf_data[,c(3,10:296)]
vcf_data_test <- as.data.frame(vcf_data_test)
#vcf_data_test %>% select("V17")

# Separate each column, keeping only those ending with "_genotype"
library(purrr)
vcf_data_final_test = names(vcf_data_test) %>%
  map(
    function(x) 
      vcf_data_test %>% 
      dplyr::select(x) %>% 
      separate(x, 
               into = paste0(x, c("_genotype", "aa")), 
               sep = ":")  %>% 
      dplyr::select(paste0(x, "_genotype"))
  ) %>%
  bind_cols()


colnames(vcf_data_final_test) <- c("target",TCGA_combined_surv_genotype$Row.names)
vcf_data_final_test[vcf_data_final_test == '0/0'] <- '0'
vcf_data_final_test[vcf_data_final_test == '0/1' | vcf_data_final_test == '0|1'] <- '1'
vcf_data_final_test[vcf_data_final_test == '1/1' | vcf_data_final_test == '1|1' ] <- '2'
vcf_data_final_test[vcf_data_final_test == './.'] <- '0'
vcf_data_final_test[vcf_data_final_test == '0|0'] <- '0'

vcf_data_final_test <- t(vcf_data_final_test)
colnames(vcf_data_final_test) <- vcf_data_final_test[1,]
vcf_data_final_test <- vcf_data_final_test[-1,]
vcf_data_final_test <- as.data.frame(vcf_data_final_test)


##  8.3 start with survival analysis data preparation ------------------
vcf_data_final_test[vcf_data_final_test == './.'] <- '0'
vcf_data_final_test[vcf_data_final_test == '0|0'] <- '0'
vcf_data_final_test[vcf_data_final_test == '0/0'] <- '0'
vcf_data_final_test[vcf_data_final_test == '0/1' | vcf_data_final_test == '0|1'] <- '1'
vcf_data_final_test[vcf_data_final_test == '1/1' | vcf_data_final_test == '1|1' ] <- '2'
vcf_data_final_test[vcf_data_final_test == './.'] <- '0'
vcf_data_final_test[vcf_data_final_test == '0|0'] <- '0'
View(vcf_data_final_test)

table(colnames(vcf_data_final_test) %in% unique(filtered_target_pairs$genotype))
#### 209 variants 

#### compare the effect of joint calling result with the single-sample calling and merge result: similar, a bit different
vcf_data_final_test$chr10_80032006_T_C == TCGA_combined_surv_genotype$chr10_80032006_T_C
vcf_data_final_test$chr11_206089_G_A == TCGA_combined_surv_genotype$chr11_206089_G_A

#removed KM analysis for both gene and variatns version
### Filter columns of vcf+inal_test based on matching genotypes -------------------

vcf_data_final_test <- vcf_data_final_test[, colnames(vcf_data_final_test) %in% unique(filtered_target_pairs$genotype)]

target_snp <- unique(filtered_target_pairs$genotype)
#intersect(target_snp, colnames(TCGA_combined_surv_genotype))
intersect(target_snp, colnames(vcf_data_final_test))
#209
dim(vcf_data_final_test)
# [1] 287 362
dim(TCGA_combined_surv_genotype)
# [1] 287 130
TCGA_combined_surv_genotype$Row.names
colnames(vcf_data_final_test)
#View(vcf_data_final_test)
rownames(vcf_data_final_test)
#rownames(vcf_data_final_test) <- sapply(strsplit(rownames(vcf_data_final_test), "-"), function(x) paste(x[1:3], collapse = "-"))
# Print the modified row names
print(rownames(vcf_data_final_test))
# Print the modified row names and check whether they are same
rownames(vcf_data_final_test) == clinical_data_TCGA_clean$`sample ID`
# True
table(rownames(vcf_data_final_test) == new_TCGA_combined_surv_genotype$Row.names)
#rm(TCGA_combined_surv_genotype_cox_crop)
#rm(TCGA_combined_surv_genotype_KM_crop)



