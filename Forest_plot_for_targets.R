# 13. Forest plot ========================
library(forestploter)
library(grid)
library(survival)
# 13.1  Prepare the data and define the function ===================
significant_snps<- unique(final_pairs$genotype)
#significant_snps <- c("chr6_26370772_G_C", "chr6_28259826_C_A", "chr7_64992452_C_T", 
#                      "chr7_65005124_G_A", "chr9_42852947_G_A", "chr9_42852958_G_T", 
#                      "chr9_42921598_T_C", "chr10_80032436_A_G", "chr11_74994651_T_C", 
#                      "chr15_40370468_C_T", "chr21_44074302_G_A", "chr21_44092213_G_A", 
#                      "chr22_49553927_C_T")

# new_significant_snps <- c("chr10_80032436_A_G", "chr15_40370468_C_T", "chr21_44074302_G_A", 
# "chr21_44092213_G_A", "chr7_64992452_C_T", "chr7_65005124_G_A", 
# "chr9_42852947_G_A", "chr9_42852958_G_T")
 
# Function to add signs for p-values ===================
add_asterisks <- function(p_value) {
  if (is.na(p_value)) {
    return("")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# Function to fit Cox models for each SNP ===================
fit_cox_models_JP <- function(snp, genotype_data, survival_data, adjust = FALSE) {
  if (adjust) {
    additive_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ as.numeric(genotype_data[snp, ]) + survival_data$T + survival_data$N + survival_data$M)
    categorical_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ factor(genotype_data[snp, ]) + survival_data$T + survival_data$N + survival_data$M)
  } else {
    additive_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ as.numeric(genotype_data[snp, ]))
    categorical_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ factor(genotype_data[snp, ]))
  }
  
  list(additive = additive_model, categorical = categorical_model)
}

fit_cox_models_TCGA <- function(snp, genotype_data, survival_data, adjust = FALSE) {
  if (adjust) {
    additive_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ as.numeric(genotype_data[snp, ]) + survival_data$T )
    categorical_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ factor(genotype_data[snp, ]) + survival_data$T)
  } else {
    additive_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ as.numeric(genotype_data[snp, ]))
    categorical_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ factor(genotype_data[snp, ]))
  }
  
  list(additive = additive_model, categorical = categorical_model)
}


## 13.2 Fit models for JP and TCGA cohorts =========================
cox_models_jp <- lapply(significant_snps, fit_cox_models_JP, genotype_data = sig_variant_genotype_japan, survival_data = survival)
cox_models_tcga <- lapply(significant_snps, fit_cox_models_TCGA, genotype_data = t_variant , survival_data = TCGA_survival_matched)
names(cox_models_jp) <- significant_snps
names(cox_models_tcga) <- significant_snps

cox_models_jp_adjusted <-lapply(significant_snps, fit_cox_models_JP, genotype_data = sig_variant_genotype_japan, survival_data = survival,adjust =TRUE)
cox_models_tcga_adjusted <- lapply(significant_snps, fit_cox_models_TCGA, genotype_data = t_variant , survival_data = TCGA_survival_matched,adjust = TRUE)
names(cox_models_jp_adjusted) <- significant_snps
names(cox_models_tcga_adjusted) <- significant_snps

#cox_result_adjusted=summary(coxph(Surv(TCGA_survival_matched[,2], TCGA_survival_matched[,1]) ~
#                                    +                                     factor(t_variant['chr6_28259826_C_A',])))
# Function to extract results from Cox models
extract_cox_results <- function(model_list, cohort_name) {
  result_list <- lapply(names(model_list), function(snp) {
    models <- model_list[[snp]]
    additive_summary <- summary(models$additive)
    categorical_summary <- summary(models$categorical)
    data.frame(
      SNP = snp,
      model = c("additive", "categorical Het","categorical Hom"),
      #HR = c(exp(coef(models$additive)[1]), exp(coef(models$categorical)[2])),
      HR = c(additive_summary$conf.int[1, "exp(coef)"],categorical_summary$conf.int[1, "exp(coef)"], categorical_summary$conf.int[2, "exp(coef)"]),
      #CI_lower = c(exp(confint(models$additive)[1, 1]), exp(confint(models$categorical)[2, 1])),
      #CI_upper = c(exp(confint(models$additive)[1, 2]), exp(confint(models$categorical)[2, 2])),
      CI_lower = c(additive_summary$conf.int[1, "lower .95"], categorical_summary$conf.int[1, "lower .95"],categorical_summary$conf.int[2, "lower .95"]),
      CI_upper = c(additive_summary$conf.int[1, "upper .95"], categorical_summary$conf.int[1, "upper .95"],categorical_summary$conf.int[2, "upper .95"]),
      #p_value = c(summary(models$additive)$coefficients[1, 5], summary(models$categorical)$coefficients[2, 5]),
      p_value = c(additive_summary$coefficients[1, "Pr(>|z|)"], categorical_summary$coefficients[1, "Pr(>|z|)"],categorical_summary$coefficients[2, "Pr(>|z|)"]),
      cohort = cohort_name
    )
  })
  do.call(rbind, result_list)
}
results_jp <- extract_cox_results(cox_models_jp, "JP")
results_tcga <- extract_cox_results(cox_models_tcga, "TCGA")
results_jp_adjusted <- extract_cox_results(cox_models_jp_adjusted, "JP")
results_tcga_adjusted <- extract_cox_results(cox_models_tcga_adjusted, "TCGA")

#check 
results_jp$SNP == results_jp_adjusted$SNP
results_tcga$SNP == results_tcga_adjusted$SNP
#True
#z <- 1.96  # For 95% confidence interval
#ci_lower_log <- coef - z * se
#ci_upper_log <- coef + z * se

#HR <- exp(coef)
#CI_lower <- exp(ci_lower_log)
#CI_upper <- exp(ci_upper_log)

# Dataframe preparation ==========================================
## JP ==========================================
### Add a blank column for the forest plot to display CI 

results_jp$`HR (95% CI)` <- ifelse(is.na(results_jp$HR), "", 
                                         sprintf("%.2f (%.2f to %.2f)", 
                                                 results_jp$HR, 
                                                 results_jp$CI_lower, 
                                                 results_jp$CI_upper))

###adding the adjusted value
results_jp$`adjusted HR (95% CI)` <- ifelse(is.na(results_jp_adjusted$HR), "", 
                                   sprintf("%.2f (%.2f to %.2f)", 
                                           results_jp_adjusted$HR, 
                                           results_jp_adjusted$CI_lower, 
                                           results_jp_adjusted$CI_upper))

# addin adjusted p-value 
results_jp$p_value_adjusted <- results_jp_adjusted$p_value
results_jp$HR_adjusted <- results_jp_adjusted$HR
results_jp$CI_lower_adjusted <- results_jp_adjusted$CI_lower
results_jp$CI_upper_adjusted <- results_jp_adjusted$CI_upper
### start to modify for the plot =======================
results_jp_modified <- results_jp

# Add two blank columns for CI
results_jp_modified$`HR Range` <- paste(rep(" ", 30), collapse = " ")
# Add genotype count 
# Merge with genotype counts
results_jp_modified <- merge(results_jp_modified, genotype_counts_wide, by.x = "SNP", by.y = "SNP", all.x = TRUE)
colnames(results_jp_modified)[(dim(results_jp_modified)[2]-2):dim(results_jp_modified)[2]] <- c('Ref','Het','Hom')

for(i in 1:nrow(results_jp_modified)) {
  if(i %% 3 !=1) {
    results_jp_modified$SNP[i] <- ""
    results_jp_modified$Ref[i] <- ""
    results_jp_modified$Het[i] <- ""
    results_jp_modified$Hom[i]<- ""
  }
}

# Apply the function to the p_values
results_jp_modified$Unadjusted_P <- mapply(function(p, sig) {
  paste0(sprintf("%.3f", p), " ",sig)
}, results_jp_modified$p_value, sapply(results_jp_modified$p_value, add_asterisks))

results_jp_modified$adjusted_P <- mapply(function(p, sig) {
  paste0(sprintf("%.3f", p), " ",sig)
}, results_jp_modified$p_value_adjusted, sapply(results_jp_modified$p_value_adjusted, add_asterisks))

# 13.2 TCGA cohort plot  ----------------

significant_snps<- unique(final_pairs$genotype)
#significant_snps <- c("chr6_26370772_G_C", "chr6_28259826_C_A", "chr7_64992452_C_T", 
#                      "chr7_65005124_G_A", "chr9_42852947_G_A", "chr9_42852958_G_T", 
#                      "chr9_42921598_T_C", "chr10_80032436_A_G", "chr11_74994651_T_C", 
#                      "chr15_40370468_C_T", "chr21_44074302_G_A", "chr21_44092213_G_A", 
#                      "chr22_49553927_C_T")

# new_significant_snps <- c("chr10_80032436_A_G", "chr15_40370468_C_T", "chr21_44074302_G_A", 
# "chr21_44092213_G_A", "chr7_64992452_C_T", "chr7_65005124_G_A", 
# "chr9_42852947_G_A", "chr9_42852958_G_T", "chr9_42921598_T_C")

# Function to add signs for p-values ===================
add_asterisks <- function(p_value) {
  if (is.na(p_value)) {
    return("")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# Function to fit Cox models for each SNP ===================
fit_cox_models_JP <- function(snp, genotype_data, survival_data, adjust = FALSE) {
  if (adjust) {
    additive_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ as.numeric(genotype_data[snp, ]) + survival_data$T + survival_data$N + survival_data$M)
    categorical_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ factor(genotype_data[snp, ]) + survival_data$T + survival_data$N + survival_data$M)
  } else {
    additive_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ as.numeric(genotype_data[snp, ]))
    categorical_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ factor(genotype_data[snp, ]))
  }
  
  list(additive = additive_model, categorical = categorical_model)
}

fit_cox_models_TCGA <- function(snp, genotype_data, survival_data, adjust = FALSE) {
  if (adjust) {
    additive_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ as.numeric(genotype_data[snp, ]) + survival_data$T )
    categorical_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ factor(genotype_data[snp, ]) + survival_data$T)
  } else {
    additive_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ as.numeric(genotype_data[snp, ]))
    categorical_model <- coxph(Surv(survival_data[, 2], survival_data[, 1]) ~ factor(genotype_data[snp, ]))
  }
  
  list(additive = additive_model, categorical = categorical_model)
}


## 13.2 Fit models for JP and TCGA cohorts =========================
cox_models_jp <- lapply(significant_snps, fit_cox_models_JP, genotype_data = sig_variant_genotype_japan, survival_data = survival)
cox_models_tcga <- lapply(significant_snps, fit_cox_models_TCGA, genotype_data = t_variant , survival_data = TCGA_survival_matched)
names(cox_models_jp) <- significant_snps
names(cox_models_tcga) <- significant_snps

cox_models_jp_adjusted <-lapply(significant_snps, fit_cox_models_JP, genotype_data = sig_variant_genotype_japan, survival_data = survival,adjust =TRUE)
cox_models_tcga_adjusted <- lapply(significant_snps, fit_cox_models_TCGA, genotype_data = t_variant , survival_data = TCGA_survival_matched,adjust = TRUE)
names(cox_models_jp_adjusted) <- significant_snps
names(cox_models_tcga_adjusted) <- significant_snps

#cox_result_adjusted=summary(coxph(Surv(TCGA_survival_matched[,2], TCGA_survival_matched[,1]) ~
#                                    +                                     factor(t_variant['chr6_28259826_C_A',])))
# Function to extract results from Cox models
extract_cox_results <- function(model_list, cohort_name) {
  result_list <- lapply(names(model_list), function(snp) {
    models <- model_list[[snp]]
    additive_summary <- summary(models$additive)
    categorical_summary <- summary(models$categorical)
    data.frame(
      SNP = snp,
      model = c("additive", "categorical Het","categorical Hom"),
      #HR = c(exp(coef(models$additive)[1]), exp(coef(models$categorical)[2])),
      HR = c(additive_summary$conf.int[1, "exp(coef)"],categorical_summary$conf.int[1, "exp(coef)"], categorical_summary$conf.int[2, "exp(coef)"]),
      #CI_lower = c(exp(confint(models$additive)[1, 1]), exp(confint(models$categorical)[2, 1])),
      #CI_upper = c(exp(confint(models$additive)[1, 2]), exp(confint(models$categorical)[2, 2])),
      CI_lower = c(additive_summary$conf.int[1, "lower .95"], categorical_summary$conf.int[1, "lower .95"],categorical_summary$conf.int[2, "lower .95"]),
      CI_upper = c(additive_summary$conf.int[1, "upper .95"], categorical_summary$conf.int[1, "upper .95"],categorical_summary$conf.int[2, "upper .95"]),
      #p_value = c(summary(models$additive)$coefficients[1, 5], summary(models$categorical)$coefficients[2, 5]),
      p_value = c(additive_summary$coefficients[1, "Pr(>|z|)"], categorical_summary$coefficients[1, "Pr(>|z|)"],categorical_summary$coefficients[2, "Pr(>|z|)"]),
      cohort = cohort_name
    )
  })
  do.call(rbind, result_list)
}
results_jp <- extract_cox_results(cox_models_jp, "JP")
results_tcga <- extract_cox_results(cox_models_tcga, "TCGA")
results_jp_adjusted <- extract_cox_results(cox_models_jp_adjusted, "JP")
results_tcga_adjusted <- extract_cox_results(cox_models_tcga_adjusted, "TCGA")

#check 
results_jp$SNP == results_jp_adjusted$SNP
results_tcga$SNP == results_tcga_adjusted$SNP
#True
#z <- 1.96  # For 95% confidence interval
#ci_lower_log <- coef - z * se
#ci_upper_log <- coef + z * se

#HR <- exp(coef)
#CI_lower <- exp(ci_lower_log)
#CI_upper <- exp(ci_upper_log)


# 14. plot drawing==========================================

## JP ==========================================
### Add a blank column for the forest plot to display CI 

results_jp$`HR (95% CI)` <- ifelse(is.na(results_jp$HR), "", 
                                   sprintf("%.2f (%.2f to %.2f)", 
                                           results_jp$HR, 
                                           results_jp$CI_lower, 
                                           results_jp$CI_upper))

###adding the adjusted value
results_jp$`adjusted HR (95% CI)` <- ifelse(is.na(results_jp_adjusted$HR), "", 
                                            sprintf("%.2f (%.2f to %.2f)", 
                                                    results_jp_adjusted$HR, 
                                                    results_jp_adjusted$CI_lower, 
                                                    results_jp_adjusted$CI_upper))

# addin adjusted p-value 
results_jp$p_value_adjusted <- results_jp_adjusted$p_value
results_jp$HR_adjusted <- results_jp_adjusted$HR
results_jp$CI_lower_adjusted <- results_jp_adjusted$CI_lower
results_jp$CI_upper_adjusted <- results_jp_adjusted$CI_upper

# adding rsID and gene info
# Corrected vectorized solution
for (row in 1:nrow(results_jp)) {
    results_jp$rs_id[row] <- final_pairs$rs_id[match(results_jp$SNP[row], final_pairs$genotype)]
    results_jp$gene_name[row] <- final_pairs$gene_name[match(results_jp$SNP[row], final_pairs$genotype)]
}

### start to modify for the plot =======================
results_jp_modified <- results_jp

# Add two blank columns for CI
results_jp_modified$`HR Range` <- paste(rep(" ", 30), collapse = " ")
# Add genotype count 
# Merge with genotype counts
results_jp_modified <- merge(results_jp_modified, genotype_counts_wide, by.x = "SNP", by.y = "SNP", all.x = TRUE)
colnames(results_jp_modified)[(dim(results_jp_modified)[2]-2):dim(results_jp_modified)[2]] <- c('Ref','Het','Hom')

for(i in 1:nrow(results_jp_modified)) {
  if(i %% 3 !=1) {
    results_jp_modified$SNP[i] <- ""
    results_jp_modified$Ref[i] <- ""
    results_jp_modified$Het[i] <- ""
    results_jp_modified$Hom[i]<- ""
    results_jp_modified$rs_id[i]<- ""
    results_jp_modified$gene_name[i]<- ""
  }
}

# Apply the function to the p_values
results_jp_modified$Unadjusted_P <- mapply(function(p, sig) {
  paste0(sprintf("%.3f", p), " ",sig)
}, results_jp_modified$p_value, sapply(results_jp_modified$p_value, add_asterisks))

results_jp_modified$adjusted_P <- mapply(function(p, sig) {
  paste0(sprintf("%.3f", p), " ",sig)
}, results_jp_modified$p_value_adjusted, sapply(results_jp_modified$p_value_adjusted, add_asterisks))


### plot =================

tm <- forest_theme(base_size = 10,
                   #refline_gp = "solid",
                   ci_pch = c(15, 18),
                   ci_col = c("#377eb8", "#4daf4a"),
                   footnote_gp = gpar(col = "blue"),
                   legend_name = "Group",
                   legend_value = c("Unadjusted", "Adjusted"),
                   vertline_lty = c("dashed", "dotted"),
                   vertline_col = c("#d6604d", "#bababa"),
                   # Table cell padding, width 4 and heights 3
                   core = list(bg_params = list(fill = c("#f2f7ff", "white", "white")),
                               padding = unit(c(4, 3), "mm")))

p <- forest(results_jp_modified[, c("SNP", "rs_id","gene_name","Ref","Het","Hom","model", "HR (95% CI)", "adjusted HR (95% CI)","HR Range", 
                                    "Unadjusted_P","adjusted_P")],
            est = list(results_jp_modified$HR,
                       results_jp_modified$HR_adjusted),
            lower = list(results_jp_modified$CI_lower,
                         results_jp_modified$CI_lower_adjusted), 
            upper = list(results_jp_modified$CI_upper,
                         results_jp_modified$CI_upper_adjusted),
            ci_column = 10,
            sizes = 0.4,
            ref_line = 1,
            vert_line = c(0.5, 2),
            arrow_lab = c("Lower Risk", "Higher Risk"),
            xlim = c(0, 6),
            ticks_at = c(0,0.5, 1, 2, 3,4,5,6),
            nudge_y = 0.4,
            theme = tm)

setwd("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/")
pdf("forestplot_JP.pdf", width = 20, height = 14)
par(mar = c(3,3,3,3), oma = c(1, 1, 1, 1))
plot(p)
dev.off()


## TCGA ==========================================
### Add a blank column for the forest plot to display CI 

# adding rsID and gene info
# Corrected vectorized solution
for (row in 1:nrow(results_tcga)) {
  results_tcga$rs_id[row] <- final_pairs$rs_id[match(results_tcga$SNP[row], final_pairs$genotype)]
  results_tcga$gene_name[row] <- final_pairs$gene_name[match(results_tcga$SNP[row], final_pairs$genotype)]
}
results_tcga$`HR (95% CI)` <- ifelse(is.na(results_tcga$HR), "", 
                                   sprintf("%.2f (%.2f to %.2f)", 
                                           results_tcga$HR, 
                                           results_tcga$CI_lower, 
                                           results_tcga$CI_upper))

###adding the adjusted value
results_tcga$`adjusted HR (95% CI)` <- ifelse(is.na(results_tcga_adjusted$HR), "", 
                                            sprintf("%.2f (%.2f to %.2f)", 
                                                    results_tcga_adjusted$HR, 
                                                    results_tcga_adjusted$CI_lower, 
                                                    results_tcga_adjusted$CI_upper))

# addin adjusted p-value 
results_tcga$p_value_adjusted <- results_tcga_adjusted$p_value
results_tcga$HR_adjusted <- results_tcga_adjusted$HR
results_tcga$CI_lower_adjusted <- results_tcga_adjusted$CI_lower
results_tcga$CI_upper_adjusted <- results_tcga_adjusted$CI_upper
### start to modify for the plot =======================
results_tcga_modified <- results_tcga

# Add two blank columns for CI
results_tcga_modified$`HR Range` <- paste(rep(" ", 30), collapse = " ")
# Add genotype count 
# Merge with genotype counts
results_tcga_modified <- merge(results_tcga_modified, genotype_counts_wide_TCGA, by.x = "SNP", by.y = "SNP", all.x = TRUE)
colnames(results_tcga_modified)[(dim(results_tcga_modified)[2]-2):dim(results_tcga_modified)[2]] <- c('Ref','Het','Hom')

for(i in 1:nrow(results_tcga_modified)) {
  if(i %% 3 !=1) {
    results_tcga_modified$SNP[i] <- ""
    results_tcga_modified$Ref[i] <- ""
    results_tcga_modified$Het[i] <- ""
    results_tcga_modified$Hom[i]<- ""
    results_tcga_modified$rs_id[i]<- ""
    results_tcga_modified$gene_name[i]<- ""
  }
}

# Apply the function to the p_values
results_tcga_modified$Unadjusted_P <- mapply(function(p, sig) {
  paste0(sprintf("%.3f", p), " ",sig)
}, results_tcga_modified$p_value, sapply(results_tcga_modified$p_value, add_asterisks))

results_tcga_modified$adjusted_P <- mapply(function(p, sig) {
  paste0(sprintf("%.3f", p), " ",sig)
}, results_tcga_modified$p_value_adjusted, sapply(results_tcga_modified$p_value_adjusted, add_asterisks))



### plot =================
tm <- forest_theme(base_size = 10,
                   #refline_gp = "solid",
                   ci_pch = c(15, 18),
                   ci_col = c("#377eb8", "#4daf4a"),
                   footnote_gp = gpar(col = "blue"),
                   legend_name = "Group",
                   legend_value = c("Unadjusted", "Adjusted"),
                   vertline_lty = c("dashed", "dotted"),
                   vertline_col = c("#d6604d", "#bababa"),
                   # Table cell padding, width 4 and heights 3
                   core = list(bg_params = list(fill = c("#f2f7ff", "white", "white")),
                               padding = unit(c(4, 3), "mm")))

p <- forest(results_tcga_modified[, c("SNP", "rs_id","gene_name","Ref","Het","Hom","model", "HR (95% CI)", "adjusted HR (95% CI)","HR Range", 
                                      "Unadjusted_P","adjusted_P")],
            est = list(results_tcga_modified$HR,
                       results_tcga_modified$HR_adjusted),
            lower = list(results_tcga_modified$CI_lower,
                         results_tcga_modified$CI_lower_adjusted), 
            upper = list(results_tcga_modified$CI_upper,
                         results_tcga_modified$CI_upper_adjusted),
            ci_column = 10,
            sizes = 0.4,
            ref_line = 1,
            vert_line = c(0.5, 2),
            arrow_lab = c("Lower Risk", "Higher Risk"),
            xlim = c(0, 6),
            ticks_at = c(0,0.5, 1, 2, 3,4,5,6),
            nudge_y = 0.4,
            theme = tm)


pdf("forestplot_TCGA.pdf", width = 20, height = 14)
par(mar = c(3,3,3,3), oma = c(1, 1, 1, 1))
plot(p)
dev.off()


