#Analyze macroinvertebrate community data from Hunt Creek Salmon study

#load libraries
library(vegan)
library(reshape)
library(ggplot2)

#Upload .csv file
Hunt_Creek_Macroinvertebrates<-read.csv("~/Documents/MSU/Research/Hunt_Creek_Salmon/Macroinvertebrates/Hunt_Creek_Macroinvertebrates.csv", sep = ",", header = T )
#Confirm header names
names(Hunt_Creek_Macroinvertebrates)
#Create sample ID name by combining environmental variables into one name in new column
Hunt_Creek_Macroinvertebrates$SampleID<-factor(paste(Hunt_Creek_Macroinvertebrates$Reach, Hunt_Creek_Macroinvertebrates$Subreach, Hunt_Creek_Macroinvertebrates$Date, Hunt_Creek_Macroinvertebrates$Type))
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
HC_M_Matrix <- cast(HC_Inverts_Subset, SampleID + Date + Reach + Subreach + Type + Days_Since_Carcass_Introduction ~ Taxonomy, value = "Value")
#create file with only community data
str(HC_M_Matrix)
HC_Inverts_Community <- as.matrix(HC_M_Matrix[,7:ncol(HC_M_Matrix)])
#Replace "na" with 0
HC_Inverts_Community[is.na(HC_Inverts_Community)]<-0
str(HC_Inverts_Community)
#Separate environmental data
HC_Inverts_Env<-HC_M_Matrix[,1:6]
str(HC_Inverts_Env)

#Set up factor for Salmon vs Control called stream_site
levels(HC_M_Matrix$Reach)
stream_site=factor(HC_Inverts_Env$Reach, c("Salmon", "Control"))

#Set up factor for Days_Since_Carcass_Introduction
str(HC_M_Matrix$Days_Since_Carcass_Introduction)
Days_Since_Introduction_n=as.numeric(HC_Inverts_Env$Days_Since_Carcass_Introduction)
Days_Since_Introduction_f=factor(HC_Inverts_Env$Days_Since_Carcass_Introduction)
levels(Days_Since_Introduction_f)

#Months since introduction - not done yet
Months_Since_Introduction<-factor(Days_Since_Introduction)
Months_Since_Introduction<-gsub(16, 1, Months_Since_Introduction)
levels(Months_Since_Introduction)

#Set up factor for Days_Since_Carcass_Introduction_f_stream_site
D_S_C_I_f_s_s <- factor(paste(HC_Inverts_Env$Days_Since_Carcass_Introduction, HC_Inverts_Env$Reach))
levels(D_S_C_I_f_s_s)

#NMDS analysis
HC_Macroinvertebrate_NMDS<-metaMDS(HC_Inverts_Community, distance="bray")
ordiplot(HC_Macroinvertebrate_NMDS, type="n", main="Invertebrate")

#Make color schemes
two_col_vec <- c("black", "bisque3")
eleven_col_vec <- c("gray0", "gray9", "gray18", "gray27", "gray36", "gray45", "gray54", "gray63", "gray72", "gray81", "gray90")
  
#NMDS plot for stream_site Salmon vs Control
with(HC_Macroinvertebrate_NMDS, points(HC_Macroinvertebrate_NMDS, display="sites", col=two_col_vec[stream_site], pch=19, pt.bg=two_col_vec))
with(HC_Macroinvertebrate_NMDS, legend("topleft", legend=levels(stream_site), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, stream_site, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Salmon"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, stream_site, kind="se", conf=0.95, lwd=2, col="bisque3", show.groups = "Control"))

#NMDS plot for Days_Since_Introduction_f
with(HC_Macroinvertebrate_NMDS, points(HC_Macroinvertebrate_NMDS, display="sites", col=eleven_col_vec[Days_Since_Introduction_f], pch=19, pt.bg=eleven_col_vec))

#NMDS plot for Days_Since_Introction_f_stream_site

#Permanova for stream_site Salmon vs. Control
adonis(HC_Inverts_Community ~ stream_site, data=HC_Inverts_Env, permutations=999)


