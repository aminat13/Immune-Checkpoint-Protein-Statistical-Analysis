#packages
install.packages("effsize")
install.packages("pwr")
library(effsize)
library(pwr)

#QUESTION 1

#count samples in each group
sum(df$Treatment_Arm == "Rucaparib")
sum(df$Treatment_Arm == "niraparib")

#effect + power calculations
d1<-cohen.d(df$PD.L1_pre[df$Treatment_Arm=="Rucaparib"], df$PD.L1_pre[df$Treatment_Arm=="niraparib"])
d1 
pwr.t.test(d=d1$estimate, sig.level=0.01, power=0.8)

d2<-cohen.d(df$CTLA.4_pre[df$Treatment_Arm=="Rucaparib"], df$CTLA.4_pre[df$Treatment_Arm=="niraparib"])
d2 
pwr.t.test(d=d2$estimate, sig.level=0.01, power=0.8)

d3<-cohen.d(df$SIRP_pre[df$Treatment_Arm=="Rucaparib"], df$SIRP_pre[df$Treatment_Arm=="niraparib"])
d3
pwr.t.test(d=d3$estimate, sig.level=0.01, power=0.8)


#QUESTION 2

#anova tests
df_RR<-df_longer[df_longer$Treatment_Arm=="Rucaparib", ]
df_NN<-df_longer[df_longer$Treatment_Arm=="niraparib", ]

df_RR$Protein <- as.factor(df_RR$Protein)
df_NN$Protein <- as.factor(df_NN$Protein)

aov_R<-aov(Expression ~ Protein, data = df_RR)
aov_N<-aov(Expression ~ Protein, data = df_NN)
summary(aov_R)
summary(aov_N)

library(multcomp)


post_test<-glht(aov_R, linfct = mcp(Protein="Tukey"))
summary(post_test)

post_test1<-glht(aov_N, linfct = mcp(Protein="Tukey"))
summary(post_test1)

# Subset only the proteins with non-significant pairwise differences
df_RR_ns <- df_RR[df_RR$Protein %in% c("PD.1_post", "PD.L1_post", "CTLA.4_post"), ]
df_NN_ns <- df_NN[df_NN$Protein %in% c("PD.1_post", "PD.L1_post", "CTLA.4_post"), ]
aov_R_ns <- aov(Expression ~ Protein, data = df_RR_ns)
aov_N_ns <- aov(Expression ~ Protein, data = df_NN_ns)
summary(aov_R_ns)
summary(aov_N_ns)

#power tests
library(effsize)
library(pwr)

effectsize::cohens_f(aov_R)
pwr.anova.test(f=0.52, k=4, power=0.8, sig.level=0.01) 
effectsize::cohens_f(aov_N)
pwr.anova.test(f=0.50, k=4, power=0.8, sig.level=0.01) 

nrow(df_longer)
table(df_longer$Treatment_Arm=="Rucaparib")
sum(df_longer$Protein == "PD.1_post")

#QUESTION 3
df$Underweight <- df$BMI_post_treatment < 20
contingency_table <- table(df$Treatment_Arm, df$Underweight)
contingency_table
chisq.test(contingency_table) 

w<-ES.w2(contingency_table/sum(contingency_table))
w #0.131
pwr.chisq.test(w=w,df=1, sig.level=0.05, power=0.8) #455

#QUESTION 4.

pop_P<-0.23
data_P<-0.31

h<-ES.h(data_P, pop_P)
h
pwr.2p2n.test(h=h, n2=4000000, sig.level=0.05, power=0.8) 
