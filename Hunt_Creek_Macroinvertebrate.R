#Begin analyzing macroinvertebrate community data from Hunt Creek Salmon study

#Upload .csv file
Hunt_Creek_Macroinvertebrates<-read.csv("~/Documents/MSU/Research/Hunt_Creek_Salmon/Macroinvertebrates/Hunt_Creek_Macroinvertebrates.csv", sep = ",", header = T )
#Confirm header names
names(Hunt_Creek_Macroinvertebrates)
#Create sample ID name by combining environmental variables into one name in new column
Hunt_Creek_Macroinvertebrates$SampleID<-factor(paste(Hunt_Creek_Macroinvertebrates$Reach, Hunt_Creek_Macroinvertebrates$Subreach, Hunt_Creek_Macroinvertebrates$Year, Hunt_Creek_Macroinvertebrates$Month, Hunt_Creek_Macroinvertebrates$Day, Hunt_Creek_Macroinvertebrates$Type))
#Create new column that combines taxonomic variables to family level
Hunt_Creek_Macroinvertebrates$Taxonomy<-(paste(Hunt_Creek_Macroinvertebrates$Class, Hunt_Creek_Macroinvertebrates$Order, Hunt_Creek_Macroinvertebrates$Family))
#Delete rows with missing values (not IDed yet) and pupae (named Hexapod NA NA in taxonomy column)
HC_Inverts_Subset<-subset(Hunt_Creek_Macroinvertebrates, Value >= 0 & Taxonomy != "Hexapoda NA NA" & Taxonomy !="Hexapoda Trichoptera NA" & Taxonomy != "Hexapoda Plecoptera NA" & Taxonomy !="Hexapoda Diptera NA")
#Convert Taxonomy charater to factor
HC_Inverts_Subset$Taxonomy<-as.factor(HC_Inverts_Subset$Taxonomy)
#Check levels of SampleID and Taxonomy variables
str(HC_Inverts_Subset)
levels(HC_Inverts_Subset$SampleID)
levels(HC_Inverts_Subset$Taxonomy)

#Convert to community table for community analysis
library(reshape)
HC_M_Matrix <- cast(HC_Inverts_Subset, SampleID + Month + Day + Year + Reach + Subreach + Type ~ Taxonomy, value = "Value")
#create file with only community data
str(HC_M_Matrix)
HC_Inverts_Community <- as.matrix(HC_M_Matrix[,8:ncol(HC_M_Matrix)])
#Replace "na" with 0
HC_Inverts_Community[is.na(HC_Inverts_Community)]<-0
str(HC_Inverts_Community)
#Separate environmental data
HC_Inverts_Env<-HC_M_Matrix[,1:7]
str(HC_Inverts_Env)
levels(HC_M_Matrix$Reach)
stream_site=factor(HC_Inverts_Env$Reach, c("Salmon", "Control"))

#NMDS analysis
library(vegan)
HC_Macroinvertebrate_NMDS<-metaMDS(HC_Inverts_Community, distance="bray")
ordiplot(HC_Macroinvertebrate_NMDS, type="n", main="Invertebrate")
two_col_vec <- c("black", "bisque3")
with(HC_Macroinvertebrate_NMDS, points(HC_Macroinvertebrate_NMDS, display="sites", col=two_col_vec[stream_site], pch=21, cex=1))

#Permanova
adonis(HC_Inverts_Community ~ stream_site, data=HC_Inverts_Env, permutations=999)
