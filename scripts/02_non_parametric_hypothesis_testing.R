#Q1: wilcoxon tests

PD1_P<-wilcox.test(protein_data$PD.1_pre, protein_data$PD.1_post, paired=TRUE)
PDL1_P<-wilcox.test(protein_data$PD.L1_pre, protein_data$PD.L1_post, paired=TRUE)
CTLA4_P<-wilcox.test(protein_data$CTLA.4_pre, protein_data$CTLA.4_post, paired=TRUE)
SIRP_P<-wilcox.test(protein_data$SIRP_pre, protein_data$SIRP_post, paired=TRUE)

p_values<-(c(PD1_P$p.value, 
             PDL1_P$p.value,
             CTLA4_P$p.value,
             SIRP_P$p.value))
adj_p_values<-p.adjust(p_values, method="BH")

pval_table <- data.frame(
  Protein = c("PD-1", "PDL1", "CTLA4", "SIRP"),
  p_value = signif(p_values, 3))


print(pval_table)

#Q2: McNemars test

underweight_post<-protein_data$BMI_post_treatment <20  
underweight_pre<-protein_data$BMI_pre_treatment <20 

bmi_table<-table(underweight_pre, underweight_post)
bmi_df <- as.data.frame.matrix(bmi_table)
rownames(bmi_df) <- c("Pre ≥20 (Not Underweight)", "Pre <20 (Underweight)")
colnames(bmi_df) <- c("Post ≥20 (Not Underweight)", "Post <20 (Underweight)")
bmi_df

mcnemar_result<-mcnemar.test(bmi_table)

mcnemar_table <- data.frame(
  Test = "McNemar's Test",
  Chi_Squared = round(unname(mcnemar_result$statistic), 3),
  df = unname(mcnemar_result$parameter),
  p_value = signif(mcnemar_result$p.value, 3)
)

print(mcnemar_table)

#Q3: Friedman test

library(PMCMRplus)
library(FSA)
library(rcompanion)
library(tidyr)
attach(protein_data)

PD1_change<-pivot_longer(protein_data, 
                         cols= c(PD.1_pre, PD.1_during, PD.1_post),
                         names_to="WithinSubjectCondition", 
                         values_to="Expression")
PD1_friedman_result<-friedman.test(Expression~WithinSubjectCondition | Patient_number, data=PD1_change)
posthoc_PD1_result<-frdAllPairsConoverTest(y=PD1_change$Expression, 
                                           groups=PD1_change$WithinSubjectCondition,
                                           blocks=PD1_change$Patient_number)


PDL1_change<-pivot_longer(protein_data, 
                          cols= c(PD.L1_pre, PD.L1_during, PD.L1_post),
                          names_to="WithinSubjectCondition", 
                          values_to="Expression")
PDL1_friedman_result<-friedman.test(Expression~WithinSubjectCondition | Patient_number, data=PDL1_change)
posthoc_PDL1_result<-frdAllPairsConoverTest(y=PDL1_change$Expression, 
                                            groups=PDL1_change$WithinSubjectCondition,
                                            blocks=PDL1_change$Patient_number)


CTLA4_change<-pivot_longer(protein_data, 
                           cols= c(CTLA.4_pre, CTLA.4_during, CTLA.4_post),
                           names_to="WithinSubjectCondition", 
                           values_to="Expression")
CTLA4_friedman_result<-friedman.test(Expression~WithinSubjectCondition | Patient_number, data=CTLA4_change)
posthoc_CTLA4_result<-frdAllPairsConoverTest(y=CTLA4_change$Expression, 
                                             groups=CTLA4_change$WithinSubjectCondition,
                                             blocks=CTLA4_change$Patient_number)


SIRP_change<-pivot_longer(protein_data, 
                          cols= c(SIRP_pre, SIRP_during, SIRP_post),
                          names_to="WithinSubjectCondition", 
                          values_to="Expression")
SIRP_friedman_result<-friedman.test(Expression~WithinSubjectCondition | Patient_number, data=SIRP_change)
posthoc_SIRP_result<-frdAllPairsConoverTest(y=SIRP_change$Expression, 
                                            groups=SIRP_change$WithinSubjectCondition,
                                            blocks=SIRP_change$Patient_number)


friedman_all <- data.frame(
  Protein = c("PD-1", "PDL1", "CTLA4", "SIRP"),
  Statistic = c(PD1_friedman_result$statistic,
                PDL1_friedman_result$statistic,
                CTLA4_friedman_result$statistic,
                SIRP_friedman_result$statistic),
  df = c(PD1_friedman_result$parameter,
         PDL1_friedman_result$parameter,
         CTLA4_friedman_result$parameter,
         SIRP_friedman_result$parameter),
  p_value = c(PD1_friedman_result$p.value,
              PDL1_friedman_result$p.value,
              CTLA4_friedman_result$p.value,
              SIRP_friedman_result$p.value)
)

View(friedman_all)
print(friedman_all) 