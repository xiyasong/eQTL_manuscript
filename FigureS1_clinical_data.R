library(dplyr)
library(stringr)
## supplementary table 1
## Reaa: cohort patient's statistics ==================

clinical_data_TCGA <-read.csv("/Users/xiyas/eQTL-continue/eQTL-continue-survival-analysis/TCGA_KIRC_clinical_data.txt",fill = TRUE,sep = "\t")
needed_col <- c("bcr_patient_barcode","gender","race_list","ajcc_pathologic_t","stage_event_tnm_categories","vital_status","age_at_initial_pathologic_diagnosis","days_to_last_follow_up","days_to_death")

clinical_data_TCGA <- clinical_data_TCGA[,colnames(clinical_data_TCGA) %in% needed_col] 

# filter to only 287 TCGA samples
clinical_data_TCGA_clean <- clinical_data_TCGA[clinical_data_TCGA$bcr_patient_barcode %in% TCGA_combined_surv_genotype$bcr_patient_barcode,]

# chean and define living day
clinical_data_TCGA_clean$LivingDays = clinical_data_TCGA_clean$days_to_last_follow_up
clinical_data_TCGA_clean$LivingDays[clinical_data_TCGA_clean$vital_status == "dead"] = 
  clinical_data_TCGA_clean$days_to_death[clinical_data_TCGA_clean$vital_status == "dead"]
clinical_data_TCGA_clean$cohort <- 'TCGA'

colnames(clinical_data_TCGA_clean) <- c("sample ID", "Status", "days_to_last_follow_up", 
                                        "days_to_death", "Gender", "Age","Race","Stage", "LivingDays","cohort")

clinical_data_TCGA_clean <- clinical_data_TCGA_clean %>% dplyr::select(`sample ID`,Gender,Race, Stage,Status,Age,LivingDays,cohort)

# ccRCC clinical data 
clinical_data_ccRCC <- surv_data %>% filter (surv_data$`sample ID` %in% combined_surv_genotype_131$`sample ID`)
clinical_data_ccRCC$cohort <- 'ccRCC'

combined_df <- rbind(clinical_data_TCGA_clean, clinical_data_ccRCC)
## Gender customized ==================
combined_df <- combined_df %>%
  mutate(Gender = case_when(
    Gender == "F" ~ "FEMALE",
    Gender == "M" ~ "MALE",
    TRUE ~ Gender  # Keeps the original value for any other cases
  ))

## Race customized ==================
combined_df <- combined_df %>%
  mutate(Race = case_when(
    Race == "Asian" ~ "ASIAN",
    TRUE ~ Race  # Keeps the original value for any other cases
  ))
## Tumor Stage customized ===============
#combined_df$T <- str_extract(combined_df$Stage, "T[0-9a-zA-Z]+(?=N)")

combined_df$T <- str_extract(combined_df$Stage, "T[0-9]")
combined_df$N <- str_extract(combined_df$Stage, "N[0-9a-zA-Z]+(?=M)")
combined_df$M <- str_extract(combined_df$Stage, "M[0-9a-zA-Z]")
head(combined_df)

## Make summary table ==================
library(dplyr)
combined_df$Age <- as.numeric(as.character(combined_df$Age))
summary_table <- combined_df %>%
  group_by(cohort) %>%
  summarise(
    Mean_Age = mean(Age, na.rm = TRUE),
    Median_Age = median(Age, na.rm = TRUE),
    Mean_LivingDays = mean(LivingDays, na.rm = TRUE),
    Median_LivingDays = median(LivingDays, na.rm = TRUE),
    Gender_Distribution = list(table(Gender)),
    Stage_Distribution = list(table(Stage)),
    Status_Distribution = list(table(Status))
  )
## Ploting ==================
# Age Distribution
library(ggplot2)
library(ggsci)

# Age Distribution
p_age <- ggplot(combined_df, aes(x = Age, fill = cohort)) +
  geom_histogram(binwidth = 5) +
  facet_wrap(~ cohort) +
  ggtitle("Age Distribution by Cohort") +
  scale_fill_nejm() +  # Example color scale
  theme_minimal()

# Living Days Distribution
p_living_days <- ggplot(combined_df, aes(x = LivingDays,fill = cohort)) +
  geom_histogram(binwidth = 500) +
  facet_wrap(~ cohort) +
  scale_fill_nejm()+
  ggtitle("Living Days Distribution by Cohort") +
  theme_minimal()

p_living_days
# Gender Distribution
p_gender <- ggplot(combined_df, aes(x = Gender, fill = Gender)) +
  geom_bar() +
  facet_wrap(~ cohort) +
  scale_fill_nejm()+
  ggtitle("Gender Distribution by Cohort") +
  theme_minimal()

# Race Distribution
p_race <- ggplot(combined_df, aes(x = Race, fill = Race)) +
  geom_bar() +
  facet_wrap(~ cohort) +
  ggtitle("Race Distribution by Cohort") +
  scale_fill_nejm()+
  theme_minimal()

# T Stage Distribution
p_T <- ggplot(combined_df, aes(x = T, fill = T)) +
  geom_bar() +
  facet_wrap(~ cohort) +
  scale_fill_nejm()+
  ggtitle("T Stage Distribution by Cohort") +
  theme_minimal()

# N Stage Distribution
p_N <- ggplot(combined_df, aes(x = N, fill = N)) +
  geom_bar() +
  facet_wrap(~ cohort) +
  scale_fill_nejm()+
  ggtitle("N Stage Distribution by Cohort") +
  theme_minimal()

# M Stage Distribution
p_M <- ggplot(combined_df, aes(x = M, fill = M)) +
  geom_bar() +
  facet_wrap(~ cohort) +
  scale_fill_nejm()+
  ggtitle("M Stage Distribution by Cohort") +
  theme_minimal()

library(gridExtra)
library(ggplot2)
# Figure S1 not plotting race now ===============
Figure_S1 <- grid.arrange(p_age, p_living_days, p_gender, p_T, p_N, p_M, ncol = 2)
Figure_S1 <- ggplotGrob(Figure_S1)
ggsave("FigureS1.pdf", plot = Figure_S1, width = 10, height = 8, units = "in")

