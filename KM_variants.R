#KM analysis for variants 
#rm(list=ls())
#rm(list=ls())
#rm(list=ls())

##read genotype data from df1
##why the read.table miss the # row? Use fread instead it worked
install.packages('R.utils')
library(data.table)
library(purrr)
library(tidyr)
vcf_data<-fread("/Users/xiyas/Degree_Project/data/df1_for_eqtl.vcf.gz",fill = TRUE, header =TRUE)

##remove unneeded data for vcf_data

vcf_data_final <- vcf_data[,c(3,10:109)]
##Method two
vcf_data_final_a = names(vcf_data_final) %>%
  map(
    function(x) 
      vcf_data_final %>% 
      dplyr::select(x) %>% 
      separate(x, 
               into = paste0(x, c("_genotype", "aa")), 
               sep = ":")  %>% 
      dplyr::select(paste0(x, "_genotype"))
  ) %>%
  bind_cols()

#vcf_data_final_a <- cbind(VCF_ID, vcf_data_final_a)
##read significant pairs data table
signif_pairs <- fread("/Users/xiyas/Degree_Project/data/eQTL-ccRCC-final.signifpairs.txt")

##select the snps

#703个 significant variants for KM_cox_0.05_genes(117 prognostic eGenes),基因型只有579个
#sig_variant_id = signif_pairs[signif_pairs$gene_id %in% KM_cox_0.05_genes,]
#new_name <- paste(sig_variant_id$variant_id, sig_variant_id$gene_id, sep = "_")

#variants have same value, to discover:
#n_occur <- data.frame(table(sig_variant_id$variant_id))
#exclude <- n_occur[n_occur$Freq > 1,]
#sig_variant_id_list = list()
#or(i in 1:nrow(sig_variant_id)){
#  print(paste0(i, " ", sig_variant_id$variant_id[i], "_", sig_variant_id$gene_id[i]))
#  sig_variant_id_list[[sig_variant_id$variant_id[i]]] = paste0(sig_variant_id$variant_id[i], "_", sig_variant_id$gene_id[i])
#}

#sig_variant_id$variant_id[!(sig_variant_id$variant_id %in% names(sig_variant_id_list))]

###T his is for all signif pairs km analsis
#sig_variant_genotype_japan = vcf_data_final_a[vcf_data_final_a$VCF_ID %in% sig_variant_id$variant_id,]

###this is for extraction of KM-COX significant Km
sig_variant_genotype_japan_131 <- vcf_data_final_a[vcf_data_final_a$ID_genotype %in% variant_merge$ID,]
## Transform genotypes to "0","1","2"

sig_variant_genotype_japan_131[sig_variant_genotype_japan_131 == '0/0'] <- '0'
sig_variant_genotype_japan_131[sig_variant_genotype_japan_131 == '0/1' | sig_variant_genotype_japan_131 == '0|1'] <- '1'
sig_variant_genotype_japan_131[sig_variant_genotype_japan_131 == '1/1' | sig_variant_genotype_japan_131 == '1|1' ] <- '2'
sig_variant_genotype_japan_131[sig_variant_genotype_japan_131 == './.'] <- '0'
#需要加上》
sig_variant_genotype_japan_131[sig_variant_genotype_japan_131 == '0|0'] <- '0'

## column, row T
sig_variant_genotype_japan_131 <- t(sig_variant_genotype_japan_131)
colnames(sig_variant_genotype_japan_131) <- sig_variant_genotype_japan_131[1,]
#colnames(sig_variant_genotype_japan) <- new_name
sig_variant_genotype_japan_131 <- sig_variant_genotype_japan_131[-1,]

##change row names 
a <-rownames(sig_variant_genotype_japan_131)
a <-strsplit(a,split = '_')

a = sapply(a,"[[",2)
a = gsub("-tumor","",a)
a = gsub("\\-","_",a)
rownames(sig_variant_genotype_japan_131) <- a

###combine survival data with genotype data
combined_surv_genotype_131 <- merge(surv_data,sig_variant_genotype_japan_131,by = "row.names",all.x = FALSE)
combined_surv_genotype_131 <- combined_surv_genotype_131[,-1]
write.table(combined_surv_genotype_131,file= "/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/ccRCC_surv_genotype_combined_131.txt",row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")

#combined_surv_genotype <- fread("ccRCC_surv_genotype_combined_ KM_cox_0.05.txt")
#write.table(combined_surv_genotype,file= "ccRCC_surv_genotype_combined_KM_cox_0.05.txt",row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")


##start of KM survival analysis based on genotype
library(survival)
setwd("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis")
source('Cheng_toolbox_beta.R')

cancerType<-"ccRCC_genotype_586"


write.table(combined_surv_genotype,file= "ccRCC_genotype_586_tpm_target_gene_exp.txt",row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")


path_raw<-"/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/"
setwd(path_raw)
exp<-as.matrix(read.table(paste0(cancerType,"_tpm_target_gene_exp.txt"),header=T,sep="\t"))

exp <- combined_surv_genotype
exp
gene_list<-as.matrix(colnames(exp)[8:dim(exp)[2]])

gnedataDir=path_raw
path_raw
path_out_pdf<-paste0(path_raw,"KM_genotype_pdf_579/",sep = "")
path_out_pdf
dir.create(path_out_pdf)
setwd(path_out_pdf)

p_list<-NULL
cutoff_list<-NULL
coef_list<-NULL
for (j in 1:length(gene_list)){
  print(j)
  output<-Cheng_generateSurvInputfromTCGA(gene_list[j,1],cancerType,dataDir)
  setwd(path_out_pdf)
  result<-Cheng_generateKMplot(output,0.8,outFile=paste0(cancerType,"_",gene_list[j,1]))
  log_rank_p<-as.matrix(result$logRankP)
  cut_off<-0.8
  coef<-as.matrix(result$coef)
  p_list<-rbind(p_list,log_rank_p)
  cutoff_list<-rbind(cutoff_list,cut_off)
  coef_list<-rbind(coef_list,coef)
  rm(log_rank_p,cut_off,coef)
}
fdr_list<-p.adjust(p_list,method="BH",n=length(p_list))
final_result<-cbind(as.matrix(gene_list),coef_list,cutoff_list,p_list,fdr_list)
colnames(final_result)<-c("genotype","coef","cutoff","p","fdr")
write.table(final_result,file="genotype_coef_cutoff_list.txt",sep="\t",row.names=F,col.names=T,quote=F)

###select the variants with fdr <0.05
final_result <- as.data.frame(final_result)
KM_p_0.05 <- filter(final_result,p <= 0.05)
KM_p_0.05 <- KM_p_0.05[order(KM_p_0.05$p), ]
KM_p_0.05 <- KM_p_0.05 %>% filter(KM_p_0.05$p !=0)
KM_fdr_0.05 <- filter(final_result,fdr <= 0.05)


### _________________________Genotype (Japanese) Coxph analysis_____________________

#cox anlysis
rm(list=ls())
library(survival)
library(dplyr)

setwd("/Users/mstermo/eQTL-continue/eQTL-continue-survival-analysis")

target <- rownames(sig_variant_genotype_japan)

exp_tumor <- t_exp_tpm_keep
sig_variant_genotype_japan <- t(sig_variant_genotype_japan)

mean_value<-as.matrix(rowMeans(exp_tumor))

colnames(mean_value)<-"mean_TPM"

survival<-surv_data[,c("Status","LivingDays")]
survival[which(survival[,1]=="alive"),1]=0
survival[which(survival[,1]=="dead"),1]=1
survival[,2]=gsub(" ","",survival[,2])
#?gsub
#class(survival)
#sapply(survival,class)
i <- c(1,2)
survival[ , i] <- apply(survival[ , i], 2,            # Specify own function within apply
                        function(x) as.numeric(as.character(x)))
sapply(survival,class)
#survival table didn't remove the lost patients from ccRCC_1 to ccRCC_106
survival <- survival[-c(33,47,57,63,74,101),]

#as.numeric(unlist(survival))
#mode(survival)="numeric"
#survival
coef_list<-NULL
p_list<-NULL

for (i in 1:dim(sig_variant_genotype_japan)[1]){
  cox_result=summary(coxph(Surv(survival[,2], survival[,1]) ~ as.matrix(sig_variant_genotype_japan[i,])))
  coef=as.matrix(cox_result$coefficients)[1,"coef"]
  p=as.matrix(cox_result$coefficients)[1,"Pr(>|z|)"]
  coef_list<-rbind(coef_list,coef)
  p_list=rbind(p_list,p)
}

#Cox analysis for the first three variants in chr8
cox_result=summary(coxph(Surv(survival[,2], survival[,1]) ~ as.numeric(sig_variant_genotype_japan[1,])+ as.numeric(sig_variant_genotype_japan[2,])+ as.numeric(sig_variant_genotype_japan[3,])))

coef=as.matrix(cox_result$coefficients)[1,"coef"]
p=as.matrix(cox_result$coefficients)[1,"Pr(>|z|)"]

coef_list<-rbind(coef_list,coef)
p_list=rbind(p_list,p)

fdr_list<-p.adjust(p_list,method="BH",n=length(p_list))
cox_list<-cbind(coef_list,p_list,fdr_list)
colnames(cox_list)<-c("COX_coef","P","FDR")

target<-cbind(target,cox_list)
#count the significant genes
#class(target)
target <- as.data.frame(target)
cox_0.05 <- filter(target,P <= 0.05)
cox_0.001 <- filter(target, P <= 0.001)
#write?
write.table(target,file="target_one_cox_in_TCGA.txt",sep="\t",row.names=F,col.names=T,quote=F)

