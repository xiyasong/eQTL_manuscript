# 1. main KM analysis for gene-----------------------------

library(survival)
setwd("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis")
source('Cheng_toolbox_beta_gene.R')

cancerType<-"ccRCC"

#path_raw<-"E:\\L1000_2021\\kidney_drug_reposition\\coexp_network\\based_all_genes\\TPM_1\\Drug_reposition\\Cor_SHvsCP_HA1E\\KM_anaylysis\\target_gene\\"
#path_raw<-paste0(path_raw,"JAP\\")
path_raw<-"/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/"
setwd(path_raw)
exp<-as.matrix(read.table(paste0(cancerType,"_tpm_target_gene_exp.txt"),header=T,sep="\t"))
exp
gene_list<-as.matrix(colnames(exp)[8:dim(exp)[2]])

dataDir=path_raw
path_raw
path_out_pdf<-paste0(path_raw,"KM_pdf/",sep = "")
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
  result<-Cheng_generateKMplot(output,outFile=paste0(cancerType,"_",gene_list[j,1]))
  log_rank_p<-as.matrix(result$logRankP)
  cut_off<-as.matrix(result$EXPcut)
  coef<-as.matrix(result$coef)
  
  p_list<-rbind(p_list,log_rank_p)
  cutoff_list<-rbind(cutoff_list,cut_off)
  coef_list<-rbind(coef_list,coef)
  rm(log_rank_p,cut_off,coef)
}

fdr_list<-p.adjust(p_list,method="BH",n=length(p_list))
final_result<-cbind(as.matrix(gene_list),coef_list,cutoff_list,p_list,fdr_list)
colnames(final_result)<-c("gene","coef","cutoff","p","fdr")
write.table(final_result,file="coef_cutoff_list.txt",sep="\t",row.names=F,col.names=T,quote=F)

final_result <- as.data.frame(final_result)
KM_0.05 <- filter(final_result,p<= 0.05)
fdr_0.05 <- filter(final_result,fdr <=0.05)
KM_0.05_genes <- KM_0.05$gene


#2. main Cox analysis for gene-----------------------------

library(survival)
library(dplyr)

setwd("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis")

#target<-target[over_gene,]
target <- rownames(exp_tumor)
#exp_tumor<-exp_tumor[over_gene,]
#exp_tpm_keep the TPM filter is >0.1
exp_tumor <- t_exp_tpm_keep
exp_tumor <- t(exp_tumor)

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

##  2.1 Test on single gene =====================================

### 2.1.1 Univariate Unadjusted =====================================
cox_result <- summary(coxph(Surv(survival[,2], survival[,1]) ~ as.matrix(exp_tumor['ENSG00000066926.10', ])))
### 2.1.2 Adjusted =====================================
cox_result <- summary(coxph(Surv(survival[,2], survival[,1]) ~ as.matrix(exp_tumor['ENSG00000066926.10', ]) + survival$T + survival$N + survival$M
))
coef <- as.matrix(cox_result$coefficients)[1, "coef"]
p <- as.matrix(cox_result$coefficients)[1, "Pr(>|z|)"]

## 2.2 Univariate Cox survival batch processsing ================================
coef_list<-NULL
p_list<-NULL
for (i in 1:dim(exp_tumor)[1]){
  cox_result=summary(coxph(Surv(survival[,2], survival[,1]) ~ as.matrix(exp_tumor[i,])))
  coef=as.matrix(cox_result$coefficients)[1,"coef"]
  p=as.matrix(cox_result$coefficients)[1,"Pr(>|z|)"]
  coef_list<-rbind(coef_list,coef)
  p_list=rbind(p_list,p)
}
fdr_list<-p.adjust(p_list,method="BH",n=length(p_list))
cox_list<-cbind(coef_list,p_list,fdr_list)
colnames(cox_list)<-c("COX_coef","P","FDR")
target<-cbind(rownames(exp_tumor),cox_list,mean_value)
#count the significant genes
#class(target)
target <- as.data.frame(target)
cox_0.05 <- filter(target,P <= 0.05)

## 2.3 adjusted multivariate  Cox survival batch processsing ================================

coef_list<-NULL
p_list<-NULL
for (i in 1:dim(exp_tumor)[1]){
  cox_result=summary(coxph(Surv(survival[,2], survival[,1]) ~ as.matrix(exp_tumor[i,]) + survival$T + survival$N + survival$M))
  coef=as.matrix(cox_result$coefficients)[1,"coef"]
  p=as.matrix(cox_result$coefficients)[1,"Pr(>|z|)"]
  
  coef_list<-rbind(coef_list,coef)
  p_list=rbind(p_list,p)
}
fdr_list<-p.adjust(p_list,method="BH",n=length(p_list))
cox_list<-cbind(coef_list,p_list,fdr_list)
colnames(cox_list)<-c("COX_coef","P","FDR")
target_adjusted<-cbind(rownames(exp_tumor),cox_list,mean_value)
#count the significant genes
#class(target)
target_adjusted <- as.data.frame(target_adjusted)
colnames(target_adjusted)[1] <- 'gene'
cox_0.05_adjusted <- filter(target_adjusted, P <= 0.05)


#compare the cox and KM ================================

Cox_0.05_genes <- cox_0.05$target
KM_cox_0.05_genes <- intersect(KM_0.05_genes,Cox_0.05_genes)


KM_0.001_genes <- Km_0.001$gene
cox_0.001_genes <- cox_0.001$target
KM_cox_0.001_genes <- intersect(KM_0.001_genes,cox_0.001_genes)
KM_cox_0.001_genes
?intersect

#######
#compare the survival and hub genes
getwd()
hub_genes<- gene_mutation$gene_id[gene_mutation$is_hub==TRUE]

colnames(hub_genes) <- hub_genes[1,] 

hub_genes <- hub_genes[-1,]
hub_genes_id <- hub_genes$gene_id

hub_genes_id
##
KM_cox_0.05_hub <- intersect(KM_cox_0.05_genes,hub_genes_id)
KM_cox_0.001_hub <- intersect(KM_cox_0.001_genes,hub_genes_id)
KM_0.05_hub <- intersect(KM_0.05_genes,hub_genes_id)


