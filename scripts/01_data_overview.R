#import dataset
protein_data<-read.csv("data/immune_checkpoint_expression.csv")

#inspect structure 
head(protein_data)
tail(protein_data)
colnames(protein_data)
dim(protein_data)