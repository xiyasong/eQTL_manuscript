
#combined_surv_genotype <- read.table("ccRCC_genotype_tpm_target_gene_exp.txt",sep="\t",header = TRUE)

##8.4 retrive the TCGA clinical data =======================================================
clinical_data_TCGA <-read.csv("/Users/xiyas/eQTL-continue/eQTL-continue-TCGA_survival_matched-analysis/TCGA_KIRC_clinical_data.txt",fill = TRUE,sep = "\t")
needed_col <- c("bcr_patient_barcode","gender","race","ajcc_pathologic_t","stage_event_tnm_categories","vital_status","age_at_index","days_to_last_follow_up","days_to_death")

clinical_data_TCGA_clean <- clinical_data_TCGA[,colnames(clinical_data_TCGA) %in% needed_col] 
clinical_data_TCGA_clean$LivingDays = clinical_data_TCGA_clean$days_to_last_follow_up
clinical_data_TCGA_clean$LivingDays[clinical_data_TCGA_clean$vital_status == "dead"] = 
  clinical_data_TCGA_clean$days_to_death[clinical_data_TCGA_clean$vital_status == "dead"]

rownames(clinical_data_TCGA_clean) <- clinical_data_TCGA_clean$`sample ID`
new_TCGA_combined_surv_genotype <- merge(clinical_data_TCGA_clean,vcf_data_final_test,by = "row.names",all.x = FALSE)

#> dim(new_TCGA_combined_surv_genotype)
#[1] 287 218


#new_TCGA_combined_surv_genotype <- TCGA_combined_surv_genotype
# Include the tunmor grades, etc
#new_TCGA_combined_surv_genotype <- merge(clinical_data_TCGA_clean,vcf_data_final_test,by = "row.names",all.x = FALSE)



# prepare step : remove snp genotype =========================

t_variant <- t(vcf_data_final_test)
dim(t_variant)
genotype_counts_list_TCGA <- apply(t_variant, 1, function(x) as.data.frame(table(factor(x, levels = c(0, 1, 2)))))
genotype_counts_df_TCGA <- do.call(rbind, lapply(seq_along(genotype_counts_list_TCGA), function(i) {
  df <- genotype_counts_list_TCGA[[i]]
  df$SNP <- rownames(t_variant)[i]
  df
}))
genotype_counts_wide_TCGA <- genotype_counts_df_TCGA %>%
  pivot_wider(names_from = Var1, values_from = Freq, values_fill = list(Freq = 0))
colnames(genotype_counts_wide_TCGA) <- c("SNP", "Genotype_0", "Genotype_1", "Genotype_2")
# View the result
print(genotype_counts_wide_TCGA)

# genotype control
snps_with_few_genotype_TCGA <- genotype_counts_wide_TCGA %>%
  filter(Genotype_1 < 3 & Genotype_2 < 3 ) %>%
  pull(SNP)

length(snps_with_few_genotype_TCGA)
### 3 variant being removed 
#snps_with_few_genotype_TCGA
#[1] "chr16_22374113_T_C" "chr7_6003434_C_T"   "chr9_42850992_A_G" 

## not need to run anymore, this only used for plotting: 8.5 start of KM TCGA_survival_matched analysis based on genotype =======================================================

setwd("/Users/xiyas/eQTL-continue/eQTL-continue_analysis")
source('Cheng_toolbox_beta.R')
#cancerType<-"TCGA_genotype_131"
#colnames(TCGA_combined_surv_genotype)[3]<- "Status"
write.table(new_TCGA_combined_surv_genotype,file= "TCGA_genotype_2024_tpm_target_gene_exp.txt",row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
getwd()
cancerType<-"TCGA_genotype_2024"
#colnames(TCGA_combined_surv_genotype)[3]<- "Status"
setwd("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis")
path_raw<-"/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/"
setwd(path_raw)
#exp<-as.matrix(read.table(paste0(cancerType,"_tpm_target_gene_exp.txt"),header=T,sep="\t"))
exp <- new_TCGA_combined_surv_genotype
genotype_list<-as.matrix(colnames(exp)[10:dim(exp)[2]])
#> length(genotype_list)
#[1] 362
dataDir=path_raw
path_raw
path_out_pdf<-paste0(path_raw,"KM_genotype_pdf_2024_TCGA/",sep = "")
path_out_pdf
dir.create(path_out_pdf)
setwd(path_out_pdf)
p_list<-NULL
cutoff_list<-NULL
coef_list<-NULL
for (j in 1:length(genotype_list)){
  #for (j in 1:2){
  print(j)
  output<-Cheng_generateSurvInputfromTCGA(genotype_list[j,1],cancerType,dataDir)
  setwd(path_out_pdf)
  result<-Cheng_generateKMplot_genoptype(output,outFile=paste0(cancerType,"_",genotype_list[j,1]))
  log_rank_p<-as.matrix(result$logRankP)
  coef<-as.matrix(result$coef)
  p_list<-rbind(p_list,log_rank_p)
  #cutoff_list<-rbind(cutoff_list,cut_off)
  coef_list<-rbind(coef_list,coef)
  rm(log_rank_p,coef)
}
fdr_list<-p.adjust(p_list,method="BH",n=length(p_list))
final_result_TCGA<-cbind(as.matrix(genotype_list),coef_list,p_list,fdr_list)
final_result_TCGA <- as.data.frame(final_result_TCGA)
colnames(final_result_TCGA)<-c("genotype","coef","p","fdr")
#16 variants significant in KM TCGA analysis
TCGA_KM_p_0.05_2024 <- filter(final_result_TCGA,p <= 0.05)

# after remove, 15 left
TCGA_KM_p_0.05_2024 <- TCGA_KM_p_0.05_2024 %>% filter(! genotype %in% snps_with_few_genotype_TCGA)

##8.6 start of Cox survival analysis based on genotype: Repeat JP cohort analysis, different zygosity models =======================================================

setwd("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis")
#survival table didn't remove the lost patients from ccRCC_1 to ccRCC_106

TCGA_survival<-clinical_data_TCGA_clean[,c("Status","LivingDays")]
TCGA_survival[which(TCGA_survival[,1]=="alive"),1]=0
TCGA_survival[which(TCGA_survival[,1]=="dead"),1]=1
TCGA_survival[,2]=gsub(" ","",TCGA_survival[,2])


TCGA_survival$T <- as.numeric(str_extract(clinical_data_TCGA_clean$Stage, "T[0-9]") %>% str_replace("T", ""))
TCGA_survival$N <- as.numeric(str_extract(clinical_data_TCGA_clean$Stage, "N[0-9a-zA-Z]+(?=M)") %>% str_replace("N", ""))
TCGA_survival$M <- as.numeric(str_extract(clinical_data_TCGA_clean$Stage, "M[0-9a-zA-Z]") %>% str_replace("M", ""))

temp <- merge(TCGA_survival,vcf_data_final_test,by = "row.names",all.x = FALSE)
rownames(temp) <- temp$Row.names
# all TRUE
rownames(temp) == colnames(t_variant)
#clean first column because it is alreadly in rownames
temp <- temp[,-1]
#get survival metadata table with matching order with genoptype data
TCGA_survival_matched <- temp[,c(1:5)]
# all TRUE

#can only adjusted by T stage
#check the survival metadata is the numeric number
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

#all true
rownames(TCGA_survival_matched) == colnames(t_variant)
i <- c(1,2)
TCGA_survival_matched[ , i] <- apply(TCGA_survival_matched[ , i], 2,            # Specify own function within apply
                        function(x) as.numeric(as.character(x)))

### 8.5.1 Single-Cox analysis test on single variants: if adjusted for TNM stage ====================

cox_result_adjusted=summary(coxph(Surv(TCGA_survival_matched[,2], TCGA_survival_matched[,1]) ~
                                    factor(t_variant['chr10_27014676_A_T',])))




# 8.5.2 NOT ADJUSTED Get sig_cox_results without adjust tumor grade: A new cox code which test both categorical and numeric value in SNP ========================

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

for (i in 1: dim(t_variant)[1]) {
  additive_model <- coxph(Surv(TCGA_survival_matched[, 2], TCGA_survival_matched[, 1]) ~ as.numeric(t_variant[i, ]))
  additive_result <- summary(additive_model)
  coef_additive <- as.matrix(additive_result$coefficients)[1, "coef"]
  p_additive <- as.matrix(additive_result$coefficients)[1, "Pr(>|z|)"]
  coef_list_additive <- rbind(coef_list_additive, coef_additive)
  p_list_additive <- rbind(p_list_additive, p_additive)
  # Categorical model (factor)
  categorical_model <- tryCatch({
    coxph(Surv(TCGA_survival_matched[, 2], TCGA_survival_matched[, 1]) ~ factor(t_variant[i, ]))
  }, error = function(e) {
    NULL  # Return NULL if only one level is present
  })
  if (!is.null(categorical_model)) {
    categorical_result <- summary(categorical_model)
    # Extract coefficients and p-values for heterozygous and homozygous minor genotypes
    if ("factor(t_variant[i, ])1" %in% rownames(categorical_result$coefficients)) {
      coef_categorical_het <- as.matrix(categorical_result$coefficients)["factor(t_variant[i, ])1", "coef"]
      p_categorical_het <- as.matrix(categorical_result$coefficients)["factor(t_variant[i, ])1", "Pr(>|z|)"]
    } else {
      coef_categorical_het <- NA
      p_categorical_het <- NA
    }
    if ("factor(t_variant[i, ])2" %in% rownames(categorical_result$coefficients)) {
      coef_categorical_hom <- as.matrix(categorical_result$coefficients)["factor(t_variant[i, ])2", "coef"]
      p_categorical_hom <- as.matrix(categorical_result$coefficients)["factor(t_variant[i, ])2", "Pr(>|z|)"]
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
cox_table_TCGA <- data.frame(
  genotype= rownames(t_variant),
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
colnames(cox_table_TCGA)[c(8:10)] = c('LRT_pvalue','Wald_pvalue','Logrank_pvalue')

#### Remove few genotype less than 3 individulas  ========================
##### 209 left 206 snps
cox_table_TCGA <- cox_table_TCGA %>% filter(!genotype %in% snps_with_few_genotype_TCGA)
#? so which p-value in the whole model test should used? > wald test
#### Filter significant results: if take wald test as whole test  ==================

sig_cox_results_TCGA <- cox_table_TCGA[(!is.na(cox_table_TCGA$P_additive) & cox_table_TCGA$P_additive < 0.05) | 
                                    (!is.na(cox_table_TCGA$P_categorical_het) & cox_table_TCGA$P_categorical_het < 0.05) | 
                                    (!is.na(cox_table_TCGA$P_categorical_hom) & cox_table_TCGA$P_categorical_hom < 0.05) |
                                    (!is.na(cox_table_TCGA$Wald_pvalue) & cox_table_TCGA$Wald_pvalue< 0.05),]

#before removing KM results: TCGA_KM_p_0.05_2024$genotype %in% sig_cox_results_TCGA$genotype
# TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## make a plot by plot_3_p_value

# 15 sig_cox_results left in TCGA
dim(sig_cox_results_TCGA)

#sig_cox_results <- cox_table[(!is.na(cox_table$Wald_pvalue) & cox_table$Wald_pvalue< 0.05),]
#sig_cox_results <- cox_table[(!is.na(cox_table$P_additive) & cox_table$P_additive< 0.05) |
# (!is.na(cox_table$Wald_pvalue) & cox_table$Wald_pvalue< 0.05),]
# Count significant SNPs in each model
num_significant_additive <- sum(!is.na(sig_cox_results$P_additive) & sig_cox_results$P_additive < 0.05)
num_significant_categorical_het <- sum(!is.na(sig_cox_results$P_categorical_het) & sig_cox_results$P_categorical_het < 0.05)
num_significant_categorical_hom <- sum(!is.na(sig_cox_results$P_categorical_hom) & sig_cox_results$P_categorical_hom < 0.05)

# 8.5.4 ADJUSTED Get adjusted sig_cox_results which is included tumor grade as adjusted variables ========================

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
for (i in 1:dim(t_variant)[1]) {
  # Additive model (numeric)
  additive_model <- coxph(Surv(TCGA_survival_matched[, 2], TCGA_survival_matched[, 1]) ~ as.numeric(t_variant[i, ])+ TCGA_survival_matched$T)
  additive_result <- summary(additive_model)
  coef_additive <- as.matrix(additive_result$coefficients)[1, "coef"]
  p_additive <- as.matrix(additive_result$coefficients)[1, "Pr(>|z|)"]
  coef_list_additive <- rbind(coef_list_additive, coef_additive)
  p_list_additive <- rbind(p_list_additive, p_additive)
  
  # Categorical model (factor)
  categorical_model <- tryCatch({
    coxph(Surv(TCGA_survival_matched[, 2], TCGA_survival_matched[, 1]) ~ factor(t_variant[i, ]) + TCGA_survival_matched$T)
  }, error = function(e) {
    NULL  # Return NULL if only one level is present
  })
  
  if (!is.null(categorical_model)) {
    categorical_result <- summary(categorical_model)
    # Extract coefficients and p-values for heterozygous and homozygous minor genotypes
    if ("factor(t_variant[i, ])1" %in% rownames(categorical_result$coefficients)) {
      coef_categorical_het <- as.matrix(categorical_result$coefficients)["factor(t_variant[i, ])1", "coef"]
      p_categorical_het <- as.matrix(categorical_result$coefficients)["factor(t_variant[i, ])1", "Pr(>|z|)"]
    } else {
      coef_categorical_het <- NA
      p_categorical_het <- NA
    }
    
    if ("factor(t_variant[i, ])2" %in% rownames(categorical_result$coefficients)) {
      coef_categorical_hom <- as.matrix(categorical_result$coefficients)["factor(t_variant[i, ])2", "coef"]
      p_categorical_hom <- as.matrix(categorical_result$coefficients)["factor(t_variant[i, ])2", "Pr(>|z|)"]
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
cox_table_adjusted_TCGA <- data.frame(
  genotype = rownames(t_variant),
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
colnames(cox_table_adjusted_TCGA)[c(8:10)] = c('LRT_pvalue','Wald_pvalue','Logrank_pvalue')

### 
cox_table_adjusted_TCGA <- cox_table_adjusted_TCGA %>% filter(!genotype %in% snps_with_few_genotype_TCGA)

#### get "sig_cox_results_adjusted" -----------
### in this case, Wald_pvalue is showing all significant results, cause the tumor grade is extremely significant. 
##33 variants
sig_cox_results_adjusted_TCGA <- cox_table_adjusted_TCGA[(!is.na(cox_table_adjusted_TCGA$P_additive) & cox_table_adjusted_TCGA$P_additive < 0.05) | 
                                                 (!is.na(cox_table_adjusted_TCGA$P_categorical_het) & cox_table_adjusted_TCGA$P_categorical_het < 0.05) | 
                                                 (!is.na(cox_table_adjusted_TCGA$P_categorical_hom) & cox_table_adjusted_TCGA$P_categorical_hom < 0.05),]



#  15 unadjusted and 24 sig_cox_results_adjusted_TCGA
dim(sig_cox_results_TCGA)
dim(sig_cox_results_adjusted_TCGA)
#TCGA_KM_p_0.05_2024$genotype %in% sig_cox_results_adjusted_TCGA$genotype
#. [1]  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE

# 9. TCGA NEW Get variants either significant in unadjusted cox or adjusted Cox all models ======================
#union_variant_result_TCGA <- union(TCGA_KM_p_0.05_2024$genotype,sig_cox_results_TCGA$genotype) 
#union_variant_result_TCGA <- union(union_variant_result_TCGA,sig_cox_results_adjusted_TCGA$genotype) %>% as.data.frame()
union_variant_result_TCGA <- union(sig_cox_results_TCGA$genotype,sig_cox_results_adjusted_TCGA$genotype) %>% as.data.frame()
colnames(union_variant_result_TCGA) <- "genotype"
#colnames(cox_p_0.05)[1] <- "genotype"
union_variant_result_TCGA <- merge(union_variant_result_TCGA,sig_cox_results_TCGA,by = "genotype",all.x=TRUE)
#union_variant_result_TCGA <- merge(union_variant_result_TCGA,sig_cox_results_adjusted_TCGA,by = "genotype",all.x=TRUE)
#union_variant_result_TCGA <- merge(union_variant_result_TCGA,TCGA_KM_p_0.05_2024[,c("coef","genotype","p")],by = "genotype",all.x=TRUE)
print(colnames(union_variant_result_TCGA))
#colnames(union_variant_result_TCGA)[(ncol(union_variant_result_TCGA) - 1):ncol(union_variant_result_TCGA)] <- c('KM_coef', 'KM_P')

#variants_with_different_signs <- union_variant_result %>% filter(different_signs == TRUE)
# Merge with sig_cox_results_adjusted
union_variant_result_TCGA <- merge(union_variant_result_TCGA, sig_cox_results_adjusted_TCGA, by = "genotype", all.x = TRUE)

# Rename the columns from sig_cox_results_adjusted to indicate adjusted values
adjusted_columns <- colnames(sig_cox_results_adjusted_TCGA)[-1] # exclude the genotype column
new_adjusted_columns <- paste0(adjusted_columns, "_adjusted")

# Rename columns in the union_variant_result
colnames(union_variant_result_TCGA)[(ncol(union_variant_result_TCGA) - length(adjusted_columns) + 1):ncol(union_variant_result_TCGA)] <- new_adjusted_columns
colnames(union_variant_result_TCGA) <- sub("\\.x$", "", colnames(union_variant_result_TCGA))
print(colnames(union_variant_result_TCGA))
# check signs again 
# 10. TCGA Check the direction of SNP allelic affect ================================
# Define a function to check if a variant has different signs in its coefficients
# allow to have different signs currently
check_signs <- function(row) {
  # Extract coefficients from the row and convert to numeric
  coefs <- c(
    as.numeric(row["COX_coef_additive"]),
    as.numeric(row["COX_coef_categorical_het"]),
    as.numeric(row["COX_coef_categorical_hom"]),
    # as.numeric(row["KM_coef"]),
    as.numeric(row["COX_coef_additive_adjusted"]),
    as.numeric(row["COX_coef_categorical_het_adjusted"]),
    as.numeric(row["COX_coef_categorical_hom_adjusted"])
  )
  # Remove NAs
  coefs <- na.omit(coefs)
  # Check if there are any differences in signs
  return(length(unique(sign(coefs))) > 1)
}

dim(union_variant_result_TCGA)
# Apply the function to each row and add the result as a new column
union_variant_result_TCGA$different_signs <- apply(union_variant_result_TCGA, 1, check_signs)
table(union_variant_result_TCGA$different_signs)
# 16 FALSE 10 TRUE
union_variant_result_TCGA$different_signs <- apply(union_variant_result_TCGA, 1, check_signs)

table(union_variant_result_TCGA$different_signs)
# FALSE  TRUE 
# 16    10 
# 11. Compare allelic effect JP and TCGA ========================================

# Define a function to compare signs of coefficients for each model pair
# Merge the results from the Japanese and TCGA cohorts based on genotype
merged_cohort_results <- merge(union_variant_result, union_variant_result_TCGA, by = "genotype", suffixes = c("_Japan", "_TCGA"))
# where is 42852958 ?
# Define a function to compare signs of coefficients for each model pair only if both p-values are significant
compare_signs_multiple <- function(row, coef_columns_japan, coef_columns_tcga, p_columns_japan, p_columns_tcga, models, significance_threshold = 0.05) {
  satisfied_models <- list()
  
  for (i in seq_along(coef_columns_japan)) {
    # Get p-values for the current model
    p_japan <- as.numeric(row[p_columns_japan[i]])
    p_tcga <- as.numeric(row[p_columns_tcga[i]])
    
    # Check if p-values are significant in both cohorts
    if (!is.na(p_japan) & !is.na(p_tcga) & p_japan <= significance_threshold & p_tcga <= significance_threshold) {
      # Get the signs of the coefficients
      sign_japan <- sign(as.numeric(row[coef_columns_japan[i]]))
      sign_tcga <- sign(as.numeric(row[coef_columns_tcga[i]]))
      
      # If the signs are the same, add model info to the list
      if (sign_japan == sign_tcga) {
        satisfied_models <- append(satisfied_models, models[i])
      }
    }
  }
  
  # Return TRUE if any model satisfies the criteria, along with the model information
  return(list(any_same = length(satisfied_models) > 0, models = satisfied_models))
}

# Example usage with your data
coef_columns_japan <- c("COX_coef_additive_Japan", "COX_coef_categorical_het_Japan", "COX_coef_categorical_hom_Japan", "COX_coef_additive_adjusted_Japan", "COX_coef_categorical_het_adjusted_Japan", "COX_coef_categorical_hom_adjusted_Japan", "KM_coef_Japan")
coef_columns_tcga <- c("COX_coef_additive_TCGA", "COX_coef_categorical_het_TCGA", "COX_coef_categorical_hom_TCGA", "COX_coef_additive_adjusted_TCGA", "COX_coef_categorical_het_adjusted_TCGA", "COX_coef_categorical_hom_adjusted_TCGA", "KM_coef_TCGA")
p_columns_japan <- c("P_additive_Japan", "P_categorical_het_Japan", "P_categorical_hom_Japan", "P_additive_adjusted_Japan", "P_categorical_het_adjusted_Japan", "P_categorical_hom_adjusted_Japan", "KM_P_Japan")
p_columns_tcga <- c("P_additive_TCGA", "P_categorical_het_TCGA", "P_categorical_hom_TCGA", "P_additive_adjusted_TCGA", "P_categorical_het_adjusted_TCGA", "P_categorical_hom_adjusted_TCGA", "KM_P_TCGA")
model_names <- c("Additive", "Categorical Het", "Categorical Hom", "Additive Adjusted", "Categorical Het Adjusted", "Categorical Hom Adjusted", "KM")

results <- apply(merged_cohort_results, 1, compare_signs_multiple, coef_columns_japan, coef_columns_tcga, p_columns_japan, p_columns_tcga, model_names)

# Extract results into separate columns
merged_cohort_results$same_effect_multiple <- sapply(results, function(x) x$any_same)
merged_cohort_results$satisfied_models <- sapply(results, function(x) paste(x$models, collapse = ", "))

# 11.2 Extract the SNPs with the same allelic effect and significant on both TCGA and JP !!! =========================

true_variants <- merged_cohort_results %>% filter(same_effect_multiple == TRUE)

# 12.We should get [final_paris] metaanalysis significant signal on both variant and gene lever right ? --------------------
#final_pairs <- signif_pairs %>%  filter(variant_id %in% true_variants$genotype)
#final_pairs$gene_prognosis <- final_pairs$gene_id %in% filtered_target_pairs$gene_id
#final_pairs <- final_pairs[final_pairs$gene_prognosis == TRUE,]
#dim(final_pairs)

# 25 13 (before moved KM gene and KM variants: 36)
# test method 2:
final_pairs <- filtered_target_pairs %>%  filter(genotype %in% true_variants$genotype)
dim(final_pairs)
#dim(final_pairs_test)
#[1] 25 43
#they are the same meaning 

# Initialize the Ensembl BioMart
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# List of gene IDs to retrieve names for (replace with your actual gene IDs)
#final_pairs$simple_id <- sapply(strsplit(final_pairs$gene_id,'\\.'),'[',1)
#list <- getBM(attributes=c('gene_biotype','ensembl_gene_id','external_gene_name'),
#             filters = 'ensembl_gene_id',
#             values = final_pairs$simple_id, mart = ensembl)
#for (row in 1:nrow(final_pairs)) {
#  if (final_pairs$simple_id[row] %in% list$ensembl_gene_id){
#    final_pairs$gene_biotype[row] <- list[list$ensembl_gene_id ==final_pairs$simple_id[row],]$gene_biotype
#    final_pairs$gene_name[row] <- list[list$ensembl_gene_id ==final_pairs$simple_id[row],]$external_gene_name
#  }
#  else{
#    final_pairs$gene_biotype[row] <- NA
#    final_pairs$gene_name[row] <- NA
#  }
#}

# Add satisfied_models to final_pairs
final_pairs <- final_pairs %>%
  left_join(true_variants %>% select(genotype, satisfied_models), by = c("genotype" = "genotype"))
final_pairs <-  final_pairs[final_pairs$genotype != "chr9_42921598_T_C",]
# 8 snps
length(unique(final_pairs$genotype))
# 
length(unique(final_pairs$gene_id))

final_pairs$gene_km <- final_pairs$gene_id %in% KM_0.05$gene
final_pairs$gene_cox_unadjusted <- final_pairs$gene_id %in% cox_0.05$gene
final_pairs$gene_cox_adjusted <- final_pairs$gene_id %in% cox_0.05_adjusted$gene
final_pairs$gene_cox_all <- final_pairs$gene_id %in% union(cox_0.05_adjusted$gene, cox_0.05$gene)

# Result seems like have KM or not is not affecting the results?
table(final_pairs$gene_km)
table(final_pairs$gene_cox_unadjusted)
table(final_pairs$gene_cox_adjusted)

table(final_pairs$gene_cox_all)


final_pairs %>%  filter(gene_km ==FALSE) %>% select(gene_name) %>% unique()
#only PWP2 is not significant in gene KM analysis, but it is significant in cox analysis

final_pairs %>% filter(gene_biotype == "protein_coding") %>% filter(gene_km ==TRUE) %>% select(gene_name) %>% unique()

#gene_name
#<char>
#  1:   GXYLT1P5
#2: CR848007.2
#3:     ERV3-1
#4:      DISP2

final_pairs %>% filter(gene_cox_all ==TRUE) %>% select(gene_name) %>% unique()
#gene_name
#<char>
#  1:    SEPHS1P1
#2:      ERV3-1 @@
#3:  ANKRD20A7P
#4:        <NA>
#  5:     SNX18P5
#6:    CYP4F60P
#7:            
#  8: TMEM254-AS1
#9:       DISP2 @@
#10:        PWP2 @@

final_pairs %>% filter(gene_biotype == "protein_coding") %>% filter(gene_cox_all ==TRUE) %>% select(gene_name) %>% unique()
#1:    ERV3-1
#2:    DISP2
#3:    PWP2

# 13. conclusion  summary =====================
final_pairs%>% distinct(gene_id,.keep_all = TRUE) %>%select(genotype,gene_name,gene_id,gene_biotype,gene_km,rs_id,
                                                              gene_cox_all,gene_cox_unadjusted,gene_cox_adjusted,satisfied_models) %>% as_tibble()

write.xlsx(final_pairs,file = "/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/Supplementary_Table_9_final_pairs.xlsx")
## remove this varaint because it only same signigicant in heterozygous and tcga heterozygous count is 1


length(unique(final_pairs$gene_id))
length(unique(final_pairs$genotype))
## Table 2 --------------------------
data <- final_pairs%>% distinct() %>%select(rs_id,genotype,gene_id,gene_name,gene_biotype,marker,slope,satisfied_models) %>% as_tibble()
write.table(data,file = "Table2", sep="\t", row.names=FALSE,quote = FALSE)

final_pairs%>% filter(gene_km ==TRUE) %>% select(gene_name) %>% unique()
#10 genes
final_pairs%>% filter(gene_cox_all ==TRUE) %>% select(gene_name) %>% unique()
#10 genes

final_pairs %>% filter(gene_km ==TRUE & gene_cox_all == FALSE) %>% select(gene_name,satisfied_models) %>% unique()
# 5 genes being excluded, if delete KM analysis

#gene_name                            satisfied_models
#<char>                                      <char>
#  1:    BTN3A2                    Categorical Hom Adjusted
#2:   ZNF603P                                          KM
#3:     XRRA1                                Additive, KM
#4:     GATD3                    Categorical Hom Adjusted
#5: MIR3667HG Additive Adjusted, Categorical Hom Adjusted


# let's see the one for JP only analysis
table(filtered_target_pairs$gene_id %in% KM_0.05$gene)
table(filtered_target_pairs$gene_id %in% cox_0.05$gene)
table(filtered_target_pairs$gene_id %in% cox_0.05_adjusted$gene)
table(filtered_target_pairs$gene_id %in% union(cox_0.05_adjusted$gene, cox_0.05$gene))

# LD analysis ------------------------------

# 14. somatic mutation mapping ========================

somatic_table <- read.table('/Users/xiyas/Degree_Project/data/only_somatic_gene.list.txt')
somatic_mutation_unique_gene <- read.table('/Users/xiyas/Degree_Project/data/only_somatic_gene.unique.list.txt',header = FALSE,row.names = NULL)

final_pairs$gene_name %in% somatic_mutation_unique_gene
our_final_genes <- final_pairs%>% select(gene_name,gene_id) %>% unique()
our_final_genes$somatic_status <- our_final_genes$gene_name %in% somatic_mutation_unique_gene$V1
#our_final_genes$gene_name %in% somatic_mutation_unique_gene$V1

somatic_mut <- read.table('/Users/xiyas/Degree_Project/data/somatic.merge.annotated.read.into.R.ann.vcf',header = TRUE)
somatic_mut <- somatic_mut[,c(1:8)]
annotation_cols <- strsplit(as.character(somatic_mut[, 8]), "\\|")
somatic_mut <- somatic_mut %>%
  mutate(annotation_info = somatic_mut[, 8]) %>%  # Create a temporary column to store the original annotation column
  separate(annotation_info, 
           into = c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", 
                    "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", 
                    "HGVS.c", "HGVS.p", "cDNA.pos_len", "CDS.pos_len", "AA.pos_len", 
                    "Distance", "Errors_Warnings_Info"), 
           sep = "\\|", fill = "right")  # Split on the '|' and fill missing values with NA

# View the modified dataframe
head(somatic_mut)

gene_count <- somatic_mut %>%
  group_by(Gene_Name,Gene_ID) %>%              # Group by Gene_ID
  tally() %>%                        # Count the number of occurrences (variants) per gene
  filter(n >=5)    
# remove intergenic 
gene_count <- gene_count %>% filter(!grepl("\\-",Gene_Name))
all_somatic_gene_cou <- somatic_mut %>%
  group_by(Gene_Name,Gene_ID) %>%              # Group by Gene_ID
  tally() %>%                        # Count the number of occurrences (variants) per gene
  filter(n >=1)    
write.xlsx(gene_count,file = "/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/Supplementary_Table 4_somatic_mutated_gene_loads.xlsx")
final_pairs$somatic_mutated <-final_pairs$simple_gene_id %in% all_somatic_gene_cou$Gene_ID

#unique(final_pairs$gene_name) %in% somatic_table$V1
##[1]  TRUE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE
#> table(unique(final_pairs$gene_name) %in% somatic_table$V1)
#FALSE  TRUE 
#11     5 

# 15. mapping  with machine learning results ====================

ml_41_features <- read.table('/Users/xiyas/eQTL-continue/Machine-learning/Xiya_survival_new_Han_simplify/target_41.vcf')
table(filtered_target_pairs$genotype %in% ml_41_features$V3) 
table(final_pairs$genotype %in% ml_41_features$V3) 
# 16. Make it into LDassoc tools =================

eGenes_db <- eGenes_db %>%
  separate(variant_id, into = c("chr", "pos", "ref", "alt"), sep = "_") %>%
  mutate(pos = as.numeric(pos))  # Convert pos to numeric if needed

# Step 2: Group by chr and pos (essentially by each unique SNP)
# Step 3: Find the row with the lowest pval_beta for each SNP
result_table <- eGenes_db %>%
  group_by(chr, pos) %>%
  slice(which.min(pval_beta)) %>%   # Select the row with the minimum pval_beta
  ungroup() %>%
  select(chr, pos, pval_beta)       # Select the desired columns

# Step 4: Rename the columns for clarity
result_table <- result_table %>%
  rename(pvalue = pval_beta)

# Step 5: Display or save the resulting table
dim(result_table)
print(result_table)
write.table(result_table,row.names = FALSE,file = "/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/LDassoc_result_table.txt",quote = FALSE)
