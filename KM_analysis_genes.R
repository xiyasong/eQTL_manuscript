Cheng_generateKMplot <- function(SurvInput,quantCut,outFile,figureTitle) { 
  
#  content = read.table("ASS1.txt",sep = '\t',header = TRUE)
  
#  DeadInd = content$Status %in% c('dead')
#  LivingDays = content$LivingDays
#  SurvInput= as.data.frame(LivingDays)
#  SurvInput$DeadInd = DeadInd
#  EXP = content$ENSG00000130707
#  SurvInput$EXP = EXP

  if (!hasArg("quantCut")) {
    cutoffs = sort(unique(SurvInput$EXP))
    Percentile20 = quantile(SurvInput$EXP,0.2)
    Percentile80 = quantile(SurvInput$EXP,0.8)
    cutoffs = cutoffs[cutoffs>Percentile20&cutoffs<=Percentile80]
    Pvalue = 1
    cutoff = 0
    mini = 0
    Pvalues = cutoffs
    for (i in 1:length(cutoffs)){
      tempdata = SurvInput
      tempdata$EXP[tempdata$EXP<cutoffs[i]]=0
      tempdata$EXP[tempdata$EXP>=cutoffs[i]]=1
      survOut = survfit(SurvObj ~ EXP,tempdata)
      res.km2 <- survdiff(SurvObj~ EXP, tempdata, rho=0)
      icutoff = cutoffs[i]
      iPvalue = pchisq(res.km2$chisq,length(res.km2$n)-1,lower.tail = FALSE)
      Pvalues[i] = iPvalue
      if (iPvalue < Pvalue) {
        cutoff = icutoff
        Pvalue = iPvalue
        mini = i
      }
      
    }
  } else {
    cutoff = quantile(SurvInput$EXP,quantCut)
  }
  SurvInput$EXP[SurvInput$EXP<cutoff]=0
  SurvInput$EXP[SurvInput$EXP>=cutoff]=1
  survOut = survfit(SurvObj ~ EXP,SurvInput)

  pdf(file = paste(outFile,"pdf",sep = "."),width=5.8,height=6)
  plot(survOut, col=c("red","black"), mark.time=T, cex=1.4,xlab="Time (year)",xscale = 365,lty =1, ylab = "Survival Probability",las=1, cex.lab=1.4)
  group1legend = paste("Low expression",paste("(n =",as.character(sum(SurvInput$EXP==0)),")",sep = ""), sep = " ")
  group2legend = paste("High expression",paste("(n =",as.character(sum(SurvInput$EXP==1)),")",sep = ""), sep = " ")
  if (!hasArg("figureTitle")) {
    #title("ASS1")
  } else {
    title(figureTitle)
  }
  legend(
    "bottomleft",
    legend=c(group1legend,group2legend),
    col=c("red","black"),
    horiz=FALSE,
    lty=1,
    bty='n',
    cex=1.4)
#  title("ASS1") 
  res = survdiff(SurvObj ~ EXP, SurvInput)
  logRankP = 1 - pchisq(res$chisq, length(res$n)-1)
  legend("topright",legend =c(paste0("P=",as.character(format(logRankP,scientific = TRUE,digits = 3)))),text.font=2,bty="n",cex=1.4)
  dev.off()
  sum_result<-summary(coxph(SurvObj ~ EXP, SurvInput))
  coef<-sum_result$coefficients[1]
  result = as.data.frame(logRankP)
  result$EXPcut = cutoff
  result$coef=coef
  return(result)
}

Cheng_generateSurvInputfromTCGA <- function(gene,cancerType,dataDir) { 
  #if (!hasArg("dataDir")) {
  #  dataDir = paste0("C:/Pan_cancer/TCGA_transcript_result_store/",cancerType,"/PKLR/")
  #}
  data = Cheng_readFile(paste0(dataDir,paste0(cancerType,"_tpm_target_gene_exp.txt")))
  LivingDays = as.numeric(data$LivingDays)
  SurvInput= as.data.frame(LivingDays)
  SurvInput$EXP = as.numeric(data[,gene])
  SurvInput$DeadInd = data$Status %in% c('dead')
  SurvInput$SurvObj <- with(SurvInput, Surv(LivingDays, DeadInd))
  return(SurvInput)
}


############################################### main KM analysis for gene-----------------------------

library(survival)
setwd("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis")
source('Cheng_toolbox_beta_gene.R')

cancerType<-"ccRCC"
a <- as.matrix(exp[,'ENSG00000066926.10'])
View(exp %>% select('ccRCC_1', ENSG00000066926.10))
     
path_raw<-"/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/"
setwd(path_raw)
exp<-as.matrix(read.table(paste0(cancerType,"_tpm_target_gene_exp.txt"),header=T,sep="\t"))
exp
gene_list<-as.matrix(colnames(exp)[8:dim(exp)[2]])

dataDir=path_raw
path_raw
path_out_pdf<-paste0(path_raw,"KM_pdf_new/",sep = "")
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
#write.table(final_result,file="coef_cutoff_list.txt",sep="\t",row.names=F,col.names=T,quote=F)

KM_final_result <- as.data.frame(final_result)
KM_0.05 <- filter(KM_final_result,p<= 0.05)
fdr_0.05 <- filter(KM_final_result,fdr <=0.05)
KM_0.05_genes <- KM_0.05$gene
