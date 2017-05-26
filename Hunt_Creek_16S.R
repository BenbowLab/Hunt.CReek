#Analysis of 16S microbial communities at Hunt Creek 2014-2016

#Load packages
library(reshape)
library(ggplot2)
library(vegan)

#color vectors
two_col_vec <- c("black", "bisque3")
five_col_vec<- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")

#Get tables from different runs - R1, R2, R3, R4, R5 and R6
Hunt_Creek_16S_R1<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomy_txt_tables/HC_R1_table_tabseparated.txt", sep="\t", header = T)
Hunt_Creek_16S_R2<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomy_txt_tables/HC_R2_table_tabseparated.txt", sep="\t", header = T)
Hunt_Creek_16S_R3<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomy_txt_tables/HC_R3_table_tabseparated.txt", sep="\t", header = T)
Hunt_Creek_16S_R4<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomy_txt_tables/HC_R4_table_tabseparated.txt", sep="\t", header = T)
Hunt_Creek_16S_R5<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomy_txt_tables/HC_R5_table_tabseparated.txt", sep="\t", header = T)
Hunt_Creek_16S_R6<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomy_txt_tables/HC_R6_table_tabseparated.txt", sep="\t", header = T)

#Classify as matrices
Hunt_Creek_16S_R1<-data.frame(Hunt_Creek_16S_R1)
Hunt_Creek_16S_R2<-data.frame(Hunt_Creek_16S_R2)
Hunt_Creek_16S_R3<-data.frame(Hunt_Creek_16S_R3)
Hunt_Creek_16S_R4<-data.frame(Hunt_Creek_16S_R4)
Hunt_Creek_16S_R5<-data.frame(Hunt_Creek_16S_R5)
Hunt_Creek_16S_R6<-data.frame(Hunt_Creek_16S_R6)

#Merge files based on denovo
HC_16S_OTU_R1R2<-merge(Hunt_Creek_16S_R1, Hunt_Creek_16S_R2, by=c("SampleID"), all = TRUE)
HC_16S_OTU_R3R4<-merge(Hunt_Creek_16S_R3, Hunt_Creek_16S_R4, by=c("SampleID"), all = TRUE)
HC_16S_OTU_R1R2R3R4<-merge(HC_16S_OTU_R1R2, HC_16S_OTU_R3R4, by=c("SampleID"), all = TRUE)
HC_16S_OTU_R5R6<-merge(Hunt_Creek_16S_R5, Hunt_Creek_16S_R6, by=c("SampleID"), all = TRUE)
HC_16S_OTU<-merge(HC_16S_OTU_R1R2R3R4, HC_16S_OTU_R5R6, by=c("SampleID"), all = TRUE)

#Format data frame so the denovo is row name
row.names(HC_16S_OTU)<-HC_16S_OTU[,1]

#Delete taxonomy and denovo columns, now that denovo is row name
HC_16S_OTU$SampleID<-NULL
HC_16S_OTU$Consensus.Lineage.x.x<-NULL
HC_16S_OTU$Consensus.Lineage.y.x<-NULL
HC_16S_OTU$Consensus.Lineage.x.y<-NULL
HC_16S_OTU$Consensus.Lineage.y.y<-NULL
HC_16S_OTU$Consensus.Lineage.y<-NULL
HC_16S_OTU$Consensus.Lineage.x<-NULL

#Replace "na" with 0
HC_16S_OTU<-as.matrix(HC_16S_OTU)
HC_16S_OTU[is.na(HC_16S_OTU)] <- 0

#transpose
HC_16S_OTU_t<-t(HC_16S_OTU)
HC_16S_OTU_t<-data.frame(HC_16S_OTU_t)
str(HC_16S_OTU_t)
names(HC_16S_OTU_t)

#Get metadata through mapping file
Hunt_Creek_16S_map <- read.table("/Users/courtneylarson/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Mapping_files/Hunt_Creek_Map.txt", header=T)
row.names(Hunt_Creek_16S_map)<-Hunt_Creek_16S_map[,1]
Hunt_Creek_16S_map<-Hunt_Creek_16S_map[,-c(1)]

#Merge metadata onto data table
HC_16S_OTU_map <-merge(Hunt_Creek_16S_map, HC_16S_OTU_t, by=0)

#make data table for just biofilm samples and delete columns that add up to 0
HC_16S_OTU_map_BF<-subset(HC_16S_OTU_map, Source=="Biofilm")
HC_16S_OTU_map_BF_v<-HC_16S_OTU_map_BF[, colSums(HC_16S_OTU_map_BF != 0) > 0]

#make data table for just carcass samples and delete columns that add up to 0
HC_16S_OTU_map_C<-subset(HC_16S_OTU_map, Source=="Carcass")
HC_16S_OTU_map_C_v<-HC_16S_OTU_map_C[, colSums(HC_16S_OTU_map_C != 0) > 0]

#Make data table for just insect samples
HC_16S_OTU_map_I<-subset(HC_16S_OTU_map, Source!="Carcass" & Source!="Biofilm")
HC_16S_OTU_map_I$Source
HC_16S_OTU_map_I_v<-HC_16S_OTU_map_I[, colSums(HC_16S_OTU_map_I != 0) > 0]

#Create melted file (this step takes a while, so prepare to wait. Skip if just doing community analysis like NMDS)
HC_16S_OTU_map_m<-melt(HC_16S_OTU_map, id=c("Row.names"))
HC_16S_OTU_map_m<-data.frame(HC_16S_OTU_map_m)
names(HC_16S_OTU_map_m)

#Delete zeros
HC_16S_OTU_map_m<-subset(HC_16S_OTU_map_m, value > 0)

#Find top taxa for samples
totals<-rbind(HC_16S_OTU_t, colSums(Hunt_Creek_16S_t))
totals<-totals[-c(1:*number*),]
sort(totals,decreasing=TRUE)[1:6]
#Use names given as output to reduce data frame for visualization
HC_16S_OTU_map_m_sub<-subset(HC_16S_OTU_map_m, variable=="denovo20259" | variable=="denovo5011" | variable=="denovo20477" | variable=="denovo17186" | variable=="denovo17019" | variable=="denovo19917")
HC_16S_OTU_map_m_sub<-data.frame(HC_16S_OTU_map_m_sub)
str(HC_16S_OTU_map_m_sub)
HC_16S_OTU_map_m_sub$Row.names<-as.factor(HC_16S_OTU_map_m_sub$Row.names)
HC_16S_OTU_map_m_sub$Total_Biofilm_Growth_PostCarcass<-as.factor(HC_16S_OTU_map_m_sub$Total_Biofilm_Growth_PostCarcass)

#make stacked bar graphs
ggplot(HC_16S_OTU_map_m_sub, aes(x=HC_16S_OTU_map_m_sub[,1], y=value, fill=variable)) + geom_bar(position = "fill",stat = "identity")
ggplot(HC_16S_OTU_map_m_sub, aes(x=HC_16S_OTU_map_m_sub[,12], y=value, fill=variable)) + geom_bar(position = "fill",stat = "identity")

#create overal community data matrix for community analysis
HC_16S_com<-as.matrix(HC_16S_OTU_map[,15:ncol(HC_16S_OTU_map)])

#biofilm community data matrix
HC_16S_com_BF<-as.matrix(HC_16S_OTU_map_BF[,15:ncol(HC_16S_OTU_map_BF_v)])

#carcass community data matrix
HC_16S_com_C<-as.matrix(HC_16S_OTU_map_C[,15:ncol(HC_16S_OTU_map_C)])

#insect community data matrix
HC_16S_com_I<-as.matrix(HC_16S_OTU_map_I[,15:ncol(HC_16S_OTU_map_I)])

#Create overall environmental data matrix for community analysis
HC_16S_env<-(HC_16S_OTU_map[,1:14])

#biofilm environmental data matrix
HC_16S_env_BF<-(HC_16S_OTU_map_BF_v[,1:14])

#carcass environmental data matrix
HC_16S_env_C<-(HC_16S_OTU_map_C[,1:14])

#Insect environmental data matrix
HC_16S_env_I<-(HC_16S_OTU_map_I[,1:14])

#Make sample type variable for overall analysis
Sample.type<-factor(HC_16S_env$Source)

#Make stream_site variable for biofilm analysis
stream_site=factor(HC_16S_env_BF$Reach, c("Salmon", "Control"))

#Make Year variable for biofilm analysis
Study_Year<-factor(HC_16S_env_BF$Year, c("1", "2"))

#Overall NMDS for all 16S samples
HC_16S_NMDS<-metaMDS(HC_16S_com, distance="bray")
ordiplot(HC_16S_NMDS, type="n", main="16S")

#stressplot for overall nmds
stressplot(HC_16S_NMDS)

#make figure for overall 16S nmds by sample type (biofilm, carcass, insect)
with(HC_16S_NMDS, points(HC_16S_NMDS, display="sites", col=five_col_vec[Sample.type], pch=19, pt.bg=five_col_vec))
with(HC_16S_NMDS, legend("topleft", legend=levels(Sample.type), bty="n", col=five_col_vec, pch=19, pt.bg=two_col_vec, cex=0.75))

#Overall permanova
adonis(HC_16S_com ~ Reach + Month + Year + Source + Reach*Month + Reach*Year + Reach+Source + Month*Year + Month*Source + Year*Source, data=HC_16S_env, permutations=999)

#make NMDS Biofilms
HC_16S_BF_NMDS<-metaMDS(HC_16S_com_BF, distance="bray")
ordiplot(HC_16S_BF_NMDS, type="n", main="16S Biofilms")

#Stressplot for biofilm nmds
stressplot(HC_16S_BF_NMDS)

#Biofilms With site points by reach
with(HC_16S_BF_NMDS, points(HC_16S_BF_NMDS, display="sites", col=two_col_vec[stream_site], pch=19, pt.bg=two_col_vec))
with(HC_16S_BF_NMDS, legend("topleft", legend=levels(stream_site), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))
with(HC_16S_BF_NMDS, ordiellipse(HC_16S_BF_NMDS, stream_site, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Salmon"))
with(HC_16S_BF_NMDS, ordiellipse(HC_16S_BF_NMDS, stream_site, kind="se", conf=0.95, lwd=2, col="bisque3", show.groups = "Control"))

#Biofilms with site points by year
with(HC_16S_BF_NMDS, points(HC_16S_BF_NMDS, display="sites", col=two_col_vec[Study_Year], pch=19, pt.bg=two_col_vec))
with(HC_16S_BF_NMDS, legend("topleft", legend=levels(Study_Year), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))
with(HC_16S_BF_NMDS, ordiellipse(HC_16S_BF_NMDS, Study_Year, kind="se", conf=0.95, lwd=2, col="black", show.groups = "1"))
with(HC_16S_BF_NMDS, ordiellipse(HC_16S_BF_NMDS, Study_Year, kind="se", conf=0.95, lwd=2, col="bisque3", show.groups = "2"))


#biofilms permanova
adonis(HC_16S_com_BF ~ Reach + Date + Reach*Date, data=HC_16S_env_BF, permutations=999)

#make nmds for insects
HC_16S_I_NMDS<-metaMDS(HC_16S_com_I, distance="bray")
#Insufficient data right now, so run once more samples


####################################################################
####################################################################
####################################################################

#Now, work with phyla level data for 16S


#Upload phyla level files for each run
Hunt_Creek_16S_R1_P<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomic_Tables/HC_Phyla_R1_16S_c.txt", sep="\t", header = T)
Hunt_Creek_16S_R2_P<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomic_Tables/HC_Phyla_R2_16S_c.txt", sep="\t", header = T)
Hunt_Creek_16S_R3_P<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomic_Tables/HC_Phyla_R3_16S_c.txt", sep="\t", header = T)
Hunt_Creek_16S_R4_P<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomic_Tables/HC_Phyla_R4_16S_c.txt", sep="\t", header = T)
Hunt_Creek_16S_R5_P<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomic_Tables/HC_Phyla_R5_16S_c.txt", sep="\t", header = T)
Hunt_Creek_16S_R6_P<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomic_Tables/HC_Phyla_R6_16S_c.txt", sep="\t", header = T)

#Clasify as matrices
Hunt_Creek_16S_R1_P<-data.frame(Hunt_Creek_16S_R1_P)
Hunt_Creek_16S_R2_P<-data.frame(Hunt_Creek_16S_R2_P)
Hunt_Creek_16S_R3_P<-data.frame(Hunt_Creek_16S_R3_P)
Hunt_Creek_16S_R4_P<-data.frame(Hunt_Creek_16S_R4_P)
Hunt_Creek_16S_R5_P<-data.frame(Hunt_Creek_16S_R5_P)
Hunt_Creek_16S_R6_P<-data.frame(Hunt_Creek_16S_R6_P)

#Merge files based on Phyla
HC_16S_P_R1R2<-merge(Hunt_Creek_16S_R1_P, Hunt_Creek_16S_R2_P, by=c("SampleID"), all = TRUE)
HC_16S_P_R3R4<-merge(Hunt_Creek_16S_R3_P, Hunt_Creek_16S_R4_P, by=c("SampleID"), all = TRUE)
HC_16S_P_R1R2R3R4<-merge(HC_16S_P_R1R2, HC_16S_P_R3R4, by=c("SampleID"), all = TRUE)
HC_16S_P_R5R6<-merge(Hunt_Creek_16S_R5_P, Hunt_Creek_16S_R6_P, by=c("SampleID"), all = TRUE)
HC_16S_P<-merge(HC_16S_P_R1R2R3R4, HC_16S_P_R5R6, by=c("SampleID"), all = TRUE)

#Format data frame so the denovo is row name
row.names(HC_16S_P)<-HC_16S_P[,1]

#Delete taxonomy column
HC_16S_P$SampleID<-NULL

#Replace "na" with 0
HC_16S_P<-as.matrix(HC_16S_P)
HC_16S_P[is.na(HC_16S_P)] <- 0

#transpose
HC_16S_P_t<-t(HC_16S_P)
HC_16S_P_t<-data.frame(HC_16S_P_t)
str(HC_16S_P_t)
names(HC_16S_P_t)

#Merge metadata onto data table
HC_16S_P_map <-merge(Hunt_Creek_16S_map, HC_16S_P_t, by=0)

#Find most common Phyla
P_totals<-rbind(HC_16S_P_t, colSums(HC_16S_P_t))
P_totals<-P_totals[-c(1:96),]
sort(P_totals,decreasing=TRUE)[1:6]
rowSums(P_totals)

