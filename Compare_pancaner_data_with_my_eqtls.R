##########draw the eGenes level comparision among 3 datasets
library(data.table)
library(ggplot2)
library(ggVennDiagram)
install.packages("ggVennDiagram")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library("ggvenn")

Pancancer_eqtl <- read.delim("/Users/xiyas/Degree_Project/data/KIRC_tumor.cis_eQTL.txt")

#8739
length(unique(Pancancer_eqtl$gene))


#####
JP_variants <- read.table("/Users/xiyas/Degree_Project/data/snp_lookup.txt",header = TRUE)
length(which(JP_variants$rs_id_dbsnp == "."))


#410720
which(duplicated(Pancancer_eqtl$SNP))
Pancancer_eqtl$SNP[2929]
unique(Pancancer_eqtl$SNP)
length(unique(Pancancer_eqtl$SNP))

#
eGenes_pancaner <- unique(Pancancer_eqtl$gene)
eGenes_pancaner <- strsplit(eGenes_pancaner,"|",fixed = TRUE)
eGenes_pancaner <- sapply(eGenes_pancaner,"[[",1)


# The eGenes are actually quite similar, in Japanese and TCGA cohort.
intersect(eGenes_pancaner,Japan_eqtl_0.05$gene_name)
intersect(Japan_eqtl_0.05$gene_name,eGenes_pancaner)


eQTL_pancaner <-unique(Pancancer_eqtl$SNP)
eQTL_pancaner[1]
eQTL_pancaner <- strsplit(eQTL_pancaner,":",fixed = TRUE)
eQTL_pancaner <- sapply(eQTL_pancaner,"[[",1)

#intersect(eQTL_pancaner,Japan_eqtl_0.05$rs_id_dbsnp)
intersect(eQTL_pancaner,signif_pairs$variant_id)


####hypergeometric testing
length(intersect(eQTL_pancaner,JP_variants$rs_id_dbsnp))
length(intersect(signif_pairs_GTEx$variant_id,JP_variants$variant_id))





#-----------plot for eGenes comparison ------------------
GTEx_eqtl <-read.table("/Users/xiyas/Degree_Project/data/Kidney_Cortex.v8.egenes.txt",fill = TRUE,header = TRUE)
GTEx_eqtl_0.05 <- GTEx_eqtl[GTEx_eqtl$qval <= 0.05,]


Japan_eqtl <-read.table("/Users/xiyas/Degree_Project/data/eQTL-ccRCC-final.genes.annotated.txt",fill = TRUE,header = TRUE)
Japan_eqtl_0.05 <- Japan_eqtl[Japan_eqtl$qval <= 0.05,]

egene_japan <-unique(Japan_eqtl_0.05$gene_id)

egene_GTEx <- unique(GTEx_eqtl_0.05$gene_id)

hyper_test = function(set1, set2, all = 25508)
  return(phyper(sum(set1 %in% set2)-1, length(set1), all-length(set1), length(set2), lower.tail=F))


egene_japan <-unique(Japan_eqtl_0.05$gene_name)

egene_GTEx <- unique(GTEx_eqtl_0.05$gene_name)

hyper_test(egene_japan,egene_GTEx)
hyper_test(egene_japan,eGenes_pancaner)

x <- list(
  A = Japan_eqtl_0.05$gene_name,
  B = GTEx_eqtl_0.05$gene_name,
  C = eGenes_pancaner
)

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
  #ggsave(my_venn, file="my_venn.png",dpi = 600, width = 14, height =8,device = "png")
}

display_venn(
  x,
  category.names = c("Japanese ccRCC","GTEx kidney","Pancancer ccRCC"),
  output = TRUE ,
  #filename = 'venn_new.png',
  imagetype="png" ,
  height =6000 , 
  width = 6000, 
  resolution = 600,
  # Circles
  lwd = 1.5,
  #lty = 'black',
  #fill = c("#ffce80","#56B4E9","#53c68c"),
  fill=c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  # Numbers
  cex = 3,
  fontface = "italic",
  # Set names
  cat.cex = 3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.15,0.15,0.1),
  margin=0.1
)
names(x) <- c("Japanese ccRCC","GTEx kidney","Pancancer ccRCC")


##method 2 venn.diagram
venn.diagram(
  x,
  category.names = c("Japanese ccRCC","GTEx kidney","Pancancer ccRCC"),
  output = TRUE ,
  filename = 'venn_new.png',
  imagetype="png" ,
  height =8000 , 
  width = 8000, 
  resolution = 600,
  # Circles
  lwd = 1.5,
  #lty = 'black',
  #fill = c("#ffce80","#56B4E9","#53c68c"),
  fill=c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  # Numbers
  cex = 3.5,
  fontface = "italic",
  # Set names
  cat.cex = 4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.15,0.15,0.1),
  margin=0.2
)
dev.off()
ggsave(filename = "Venn_final.pdf",dpi = 600, width = 12, height =8,device = "pdf")



# plot for eQTLs comparison ---------------

signif_pairs <- fread("/Users/xiyas/Degree_Project/data/eQTL-ccRCC-final.signifpairs.txt")
signif_pairs_GTEx <- fread("/Users/xiyas/eQTL-continue/Kidney_Cortex.v8.signif_variant_gene_pairs.txt.gz")
#signif_pairs$variant_id
signif_pairs_GTEx$variant_id <- gsub("_b38","",signif_pairs_GTEx$variant_id)

unique_eqtls_GTEx<- unique(signif_pairs_GTEx$variant_id)
unique_eqtls_Japan <- unique(signif_pairs$variant_id)
U_GTEx_Jp <- intersect(unique_eqtls_GTEx,unique_eqtls_Japan)

hyper_test = function(set1, set2, all = 284774)
  return(phyper(sum(set1 %in% set2)-1, length(set1), all-length(set1), length(set2), lower.tail=F))

##pancancer hypergeometric test
set1 <- unique(signif_pairs$rs_id)
set2 <- unique(intersect(eQTL_pancaner,JP_variants$rs_id_dbsnp))

hyper_test(set1,set2)
#GTEx hypergeometric test
set1 <- unique(signif_pairs$variant_id)
set2 <- unique(intersect(signif_pairs_GTEx$variant_id,JP_variants$variant_id))
hyper_test(set1,set2)


hyper_test(unique_eqtls_pancancer,unique_eqtls_Japan)


signif_pairs <- signif_pairs %>% filter(!rs_id == ".")
dim(signif_pairs)

unique_eqtls_Japan <- unique(signif_pairs$rs_id)
unique_eqtls_pancancer <- unique(eQTL_pancaner)


test <- intersect(unique(eQTL_pancaner),unique(Supp_2$rs_id))



length(unique_eqtls_Japan)
length(unique_eqtls_pancancer)
intersect(unique_eqtls_Japan,unique_eqtls_pancancer)

######adding rs id to GTEx dataset
snp_lookup <- fread("/Users/xiyas/Degree_Project/data/snp_lookup.txt")
for (row in 1:nrow(signif_pairs_GTEx)) {
  signif_pairs_GTEx$rs_id[row] <- snp_lookup[snp_lookup$variant_id == signif_pairs_GTEx$variant_id[row],]$rs_id_dbsnp
}
for (row in 1:nrow(signif_pairs)) {
  signif_pairs$rs_id[row] <- snp_lookup[snp_lookup$variant_id == signif_pairs$variant_id[row],]$rs_id_dbsnp
}

signif_pairs_GTEx$rs_id <- 
x <- list(
  A = unique(Supp_2$rs_id),
  B = unique_eqtls_GTEx,
  C = 
)


display_venn(
  x,
  category.names = c("Japanese ccRCC","GTEx kidney","Pancancer ccRCC"),
  output = TRUE ,
  #filename = 'venn_new.png',
  imagetype="png" ,
  height =6000 , 
  width = 6000, 
  resolution = 600,
  # Circles
  lwd = 1.5,
  #lty = 'black',
  #fill = c("#ffce80","#56B4E9","#53c68c"),
  fill=c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  # Numbers
  cex = 3,
  fontface = "italic",
  # Set names
  cat.cex = 3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.15,0.15,0.1),
  margin=0.1
)
names(x) <- c("Japanese ccRCC","GTEx kidney","Pancancer ccRCC")


#######survival figure 3A and 3B



####meng's code
### n >2 (take n=5 as example)
#install.packages("venn")
library(venn)
library(ggplot2)
x <- list (
  "KM < 0.05" = KM_0.05_genes,
  "Cox < 0.05" = Cox_0.05_genes
)

#Reduce(intersect, list(Japan_eqtl_0.05$gene_name,GTEx_eqtl_0.05$gene_name,eGenes_pancaner))

pdf(file="Figure3A.pdf")
p_up <- venn(x, zcolor = c("#0073C2FF", "#EFC000FF"), ilcs=2, sncs=2,
               box=F, ggplot=F,par =TRUE)
dev.off()

### Figure 3B
x <- list (
  "KM < 0.05" = KM_p_0.05$genotype,
  "Cox < 0.05" = cox_p_0.05$target
)
pdf(file="Figure3B.pdf")
p_up <- venn(x, zcolor = c("#0073C2FF", "#EFC000FF"), ilcs=2, sncs=2,
             box=F, ggplot=F,par =TRUE)
dev.off()

###Figure 2

x <- list(
  "Japanese" = Japan_eqtl_0.05$gene_name,
  "GTEx" = GTEx_eqtl_0.05$gene_name,
  "Pan-cancer" = eGenes_pancaner
)



pdf(file="Figure2D.pdf")
p_up <- venn(x, zcolor = c("#0073C299","#EFC00099","#86868699"), ilcs=2, sncs=2,plotsize = 15,
             box=F, ggplot=F,par =TRUE)
dev.off()



KM_COX_SIG_eqtls$chr <- sapply(strsplit(KM_COX_SIG_eqtls$variant_id,"_",1),"[[",1)

the_61_eqtl$chr <- sapply(strsplit(the_61_eqtl$the_61_eqtl,"_",1),"[[",1)
