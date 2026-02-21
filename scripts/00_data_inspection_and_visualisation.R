#loading in dataset
df<-read.csv("data/protein_data.csv") 

#inspecting dataframe
head(df) 
tail(df) 
colnames(df) 
class(df) 
dim(df)
str(df)


#histograms to visualise distribution of protein expression pre treatment
par(mfrow=c(2,4), oma=c(0,0,3,0))
hist(df$PD.1_pre[df$Treatment_Arm=="Rucaparib"], xlab="Expression Change", ylab="Frequency", main="PD-1 Expression Pre Rucaparib Treatment", col="violet")
hist(df$PD.L1_pre[df$Treatment_Arm=="Rucaparib"], xlab="Expression Change", ylab="Frequency", main="PDL-1 Expression Pre Rucaparib Treatment", col="violet")
hist(df$CTLA.4_pre[df$Treatment_Arm=="Rucaparib"], xlab="Expression Change", ylab="Frequency", main="CTLA4 Expression Pre Rucaparib Treatment", col="violet")
hist(df$SIRP_pre[df$Treatment_Arm=="Rucaparib"],xlab="Expression Change",  ylab="Frequency", main="SIRP Expression Pre Rucaparib Treatment", col="violet") 

hist(df$PD.1_pre[df$Treatment_Arm=="niraparib"], xlab="Expression Change", ylab="Frequency", main="PD-1 Expression Pre Niraparib Treatment", col="purple4")
hist(df$PD.L1_pre[df$Treatment_Arm=="niraparib"], xlab="Expression Change", ylab="Frequency", main="PDL-1 Expression Pre Niraparib Treatment", col="purple4")
hist(df$CTLA.4_pre[df$Treatment_Arm=="niraparib"], xlab="Expression Change", ylab="Frequency", main="CTLA4 Expression Pre Niraparib Treatment", col="purple4")
hist(df$SIRP_pre[df$Treatment_Arm=="niraparib"], xlab="Expression Change", ylab="Frequency", main="SIRP Expression Pre Niraparib Treatment", col="purple4") 
mtext("Distribution of Protein Expression Pre-Treatment", 
      outer = TRUE, cex = 1.5)

#boxplots to visualise spread of protein expression pre treatment
par(mfrow=c(2,2), oma=c(0,0,3,0))
boxplot(df$PD.1_pre~df$Treatment_Arm, xlab="Type of Treatment", ylab="Protein Expression", main="PD-1 Expression Pre Treatment", col="violet")
boxplot(df$PD.L1_pre~df$Treatment_Arm, xlab="Type of Treatment", ylab="Protein Expression", main="PDL-1 Expression Pre Treatment", col="purple4")
boxplot(df$CTLA.4_pre~df$Treatment_Arm, xlab="Type of Treatment", ylab="Protein Expression", main="CTLA-4 Expression Pre Treatment", col="violet")
boxplot(df$SIRP_pre~df$Treatment_Arm, xlab="Type of Treatment", ylab="Protein Expression", main="SIRP Expression Pre Treatment", col="purple4")
mtext("Protein Expression Variability Post-Treatment", 
      outer = TRUE, cex = 1.5)
#
#histograms to visualise distribution of protein expression post treatment
par(mfrow=c(2,4), oma=c(0,0,3,0))
hist(df$PD.1_post[df$Treatment_Arm=="Rucaparib"], xlab="Expression Change", ylab="Frequency", main="PD-1 Expression Post Rucaparib Treatment", col="violet")
hist(df$PD.L1_post[df$Treatment_Arm=="Rucaparib"], xlab="Expression Change", ylab="Frequency", main="PDL-1 Expression Post Rucaparib Treatment", col="violet")
hist(df$CTLA.4_post[df$Treatment_Arm=="Rucaparib"], xlab="Expression Change", ylab="Frequency", main="CTLA4 Expression Post Rucaparib Treatment", col="violet")
hist(df$SIRP_post[df$Treatment_Arm=="Rucaparib"],xlab="Expression Change",  ylab="Frequency", main="SIRP Expression Post Rucaparib Treatment", col="violet") 

hist(df$PD.1_post[df$Treatment_Arm=="niraparib"], xlab="Expression Change", ylab="Frequency", main="PD-1 Expression Post Niraparib Treatment", col="purple4")
hist(df$PD.L1_post[df$Treatment_Arm=="niraparib"], xlab="Expression Change", ylab="Frequency", main="PDL-1 Expression Post Niraparib Treatment", col="purple4")
hist(df$CTLA.4_post[df$Treatment_Arm=="niraparib"], xlab="Expression Change", ylab="Frequency", main="CTLA4 Expression Post Niraparib Treatment", col="purple4")
hist(df$SIRP_post[df$Treatment_Arm=="niraparib"], xlab="Expression Change", ylab="Frequency", main="SIRP Expression Post Niraparib Treatment", col="purple4") 
mtext("Distribution of Protein Expression Post-Treatment", 
      outer = TRUE, cex = 1.5)

#boxplots to visualise spread of protein expression post treatment
par(mfrow=c(2,2), oma=c(0,0,3,0))
boxplot(df$PD.1_post~df$Treatment_Arm, xlab="Type of Treatment", ylab="Protein Expression", main="PD-1 Expression Post Treatment", col="violet")
boxplot(df$PD.L1_post~df$Treatment_Arm, xlab="Type of Treatment", ylab="Protein Expression", main="PDL-1 Expression Post Treatment", col="purple4")
boxplot(df$CTLA.4_post~df$Treatment_Arm, xlab="Type of Treatment", ylab="Protein Expression", main="CTLA-4 Expression Post Treatment", col="violet")
boxplot(df$SIRP_post~df$Treatment_Arm, xlab="Type of Treatment", ylab="Protein Expression", main="SIRP Expression Post Treatment", col="purple4")
mtext("Protein Expression Variability Post-Treatment", 
      outer = TRUE, cex = 1.5)