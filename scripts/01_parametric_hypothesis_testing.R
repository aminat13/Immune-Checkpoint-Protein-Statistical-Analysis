#Q1: Students ttest 
attach(protein_data)

PD1_ <- t.test(PD.1_pre ~ Treatment_Arm, var.equal = TRUE)
PDL1_ <- t.test(PD.L1_pre ~ Treatment_Arm, var.equal = TRUE)

ttest_table <- data.frame(
  Protein = c("PD-1", "PDL-1"),
  t_value = c(round(PD1_$statistic, 3), round(PDL1_$statistic, 3)),
  df = c(round(PD1_$parameter, 2), round(PDL1_$parameter, 2)),
  p_value = c(signif(PD1_$p.value, 3), signif(PDL1_$p.value, 3)),
  mean_diff = c(round(diff(PD1_$estimate), 3), round(diff(PDL1_$estimate), 3))
)

print(ttest_table)

#Q1: Wilcoxon Signed Rank tests
CTLA4_p<-wilcox.test(CTLA.4_pre~Treatment_Arm)
SIRP_p<-wilcox.test(SIRP_pre~Treatment_Arm)

p_values<-(c(
  CTLA4_p$p.value, 
  SIRP_p$p.value) 
)
p.adjust(p_values)

wilcox_table <- data.frame(
  Protein = c("CTLA4", "SIRP"),
  W_statistic = c(
    CTLA4_p$statistic,
    SIRP_p$statistic),
  p_value = c(
    CTLA4_p$p.value,
    SIRP_p$p.value)
)

wilcox_table$W_statistic <- round(as.numeric(wilcox_table$W_statistic), 3)
wilcox_table$p_value <- signif(as.numeric(wilcox_table$p_value), 3)

print(wilcox_table)

#Q2: fisher's test
library(RBioinf) 
library(ape)
library(rcompanion)

ER_response<-table(ER_status, Treatment_Response) 

ER_response <- matrix(
  c(11, 19,   # ER-negative (nonresponders, responders)
    39, 19),  # ER-positive (nonresponders, responders)
  nrow = 2,
  byrow = TRUE
)

rownames(ER_response) <- c("negative", "positive")
colnames(ER_response) <- c("nonresponders", "responders")

ER_response

fisher.test(ER_response)


#Q3: Kruksal-Wallis test
library(FSA)
rucaparib_treatment<-protein_data[protein_data$Treatment_Arm=="Rucaparib" , ]
protein_values<-c(rucaparib_treatment$PD.1_pre,
                  rucaparib_treatment$PD.L1_pre,
                  rucaparib_treatment$CTLA.4_pre,
                  rucaparib_treatment$SIRP_pre)
protein_groups<-c(rep("PD-1", length(rucaparib_treatment$PD.1_pre)),
                  rep("PDL-1", length(rucaparib_treatment$PD.L1_pre)),
                  rep("CTAL4", length(rucaparib_treatment$CTLA.4_pre)), 
                  rep("SIRP", length(rucaparib_treatment$SIRP_pre)))
kruskal_result<-kruskal.test(protein_values~protein_groups) 

kruskal_table <- data.frame(
  Test = "Kruskal–Wallis",
  Chi_Squared = unname(kruskal_result$statistic),
  df = unname(kruskal_result$parameter),
  p_value = kruskal_result$p.value
)

kruskal_table$Chi_Squared <- round(kruskal_table$Chi_Squared, 3)
kruskal_table$p_value <- signif(kruskal_table$p_value, 3)

print(kruskal_table)

dunn_results<-dunnTest(protein_values~protein_groups, method="bh")
results_table<-dunn_results$res
results_table
