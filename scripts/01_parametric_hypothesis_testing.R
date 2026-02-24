# tests to confirm normality
R_shapiro_PD1 <- shapiro.test(protein_data$PD.1_pre[protein_data$Treatment_Arm == "Rucaparib"])
R_shapiro_PDL1 <- shapiro.test(protein_data$PD.L1_pre[protein_data$Treatment_Arm == "Rucaparib"])
R_shapiro_CTLA4 <- shapiro.test(protein_data$CTLA.4_pre[protein_data$Treatment_Arm == "Rucaparib"])
R_shapiro_SIRP <- shapiro.test(protein_data$SIRP_pre[protein_data$Treatment_Arm == "Rucaparib"])

N_shapiro_PD1 <- shapiro.test(protein_data$PD.1_pre[protein_data$Treatment_Arm == "niraparib"])
N_shapiro_PDL1 <- shapiro.test(protein_data$PD.L1_pre[protein_data$Treatment_Arm == "niraparib"])
N_shapiro_CTLA4 <- shapiro.test(protein_data$CTLA.4_pre[protein_data$Treatment_Arm == "niraparib"])
N_shapiro_SIRP <- shapiro.test(protein_data$SIRP_pre[protein_data$Treatment_Arm == "niraparib"])

-----------------------------------------------------------------------------------------------------

  # test to confirm variance
  library(car)
R_PD1_lev <- leveneTest(protein_data$PD.1_pre ~ protein_data$Treatment_Arm == "Rucaparib")
R_PDL1_lev <- leveneTest(protein_data$PD.L1_pre ~ protein_data$Treatment_Arm == "Rucaparib")
R_CTLA4_lev <- leveneTest(protein_data$CTLA.4_pre ~ protein_data$Treatment_Arm == "Rucaparib")
R_SIRP_lev <- leveneTest(protein_data$SIRP_pre ~ protein_data$Treatment_Arm == "Rucaparib")

N_PD1_lev <- leveneTest(protein_data$PD.1_pre ~ protein_data$Treatment_Arm == "niraparib")
N_PDL1_lev <- leveneTest(protein_data$PD.L1_pre ~ protein_data$Treatment_Arm == "niraparib")
N_CTLA4_lev <- leveneTest(protein_data$CTLA.4_pre ~ protein_data$Treatment_Arm == "niraparib")
N_SIRP_lev <- leveneTest(protein_data$SIRP_pre ~ protein_data$Treatment_Arm == "niraparib")

----------------------------------------------------------------------------------------

  # Q1: Students ttest
  attach(protein_data)

PD1_ <- t.test(PD.1_pre ~ Treatment_Arm, var.equal = TRUE)
PDL1_ <- t.test(PD.L1_pre ~ Treatment_Arm, var.equal = TRUE)

--------------------------------------------------------------------

  # Q1: Wilcoxon Signed Rank tests
  CTLA4_p <- wilcox.test(CTLA.4_pre ~ Treatment_Arm)
SIRP_p <- wilcox.test(SIRP_pre ~ Treatment_Arm)

p_values <- (c(
  CTLA4_p$p.value,
  SIRP_p$p.value
)
)
p.adjust(p_values)

-----------------------------------------------------------------------

  # Q2: fisher's test
  library(RBioinf)
library(ape)
library(rcompanion)

ER_response <- table(ER_status, Treatment_Response)

ER_response <- matrix(
  c(
    11, 19, # ER-negative (nonresponders, responders)
    39, 19
  ), # ER-positive (nonresponders, responders)
  nrow = 2,
  byrow = TRUE
)

rownames(ER_response) <- c("negative", "positive")
colnames(ER_response) <- c("nonresponders", "responders")

ER_response

fisher.test(ER_response)

---------------------------------------------------------------------------

  # Q3: Kruksal-Wallis test
  library(FSA)
rucaparib_treatment <- protein_data[protein_data$Treatment_Arm == "Rucaparib", ]
protein_values <- c(
  rucaparib_treatment$PD.1_pre,
  rucaparib_treatment$PD.L1_pre,
  rucaparib_treatment$CTLA.4_pre,
  rucaparib_treatment$SIRP_pre
)
protein_groups <- c(
  rep("PD-1", length(rucaparib_treatment$PD.1_pre)),
  rep("PDL-1", length(rucaparib_treatment$PD.L1_pre)),
  rep("CTAL4", length(rucaparib_treatment$CTLA.4_pre)),
  rep("SIRP", length(rucaparib_treatment$SIRP_pre))
)
kruskal_result <- kruskal.test(protein_values ~ protein_groups)

dunn_results <- dunnTest(protein_values ~ protein_groups, method = "bh")
