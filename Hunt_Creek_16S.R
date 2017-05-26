#Analysis of 16S microbial communities at Hunt Creek 2014-2016

#Load packages
library(reshape)
library(ggplot2)
library(vegan)
library(plyr)
library(dplyr)

#color vectors
two_col_vec <- c("black", "bisque3")
three_col_vec<- c("#a6cee3", "#1f78b4", "#b2df8a")
five_col_vec<- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")
six_col_vec<- c("#e41a1c", "#377eb8", "green", "#984ea3", "#ff7f00", "#ffff33")
six_col_vec_cont<-c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")

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

#table for carcass and biofilms before and after
HC_16S_OTU_map$Date<-as.factor(HC_16S_OTU_map$Date)
HC_16S_OTU_map_CBF<-subset(HC_16S_OTU_map, Total_Biofilm_Growth_PostCarcass<=16 & Reach=="Salmon" & Source!="Rhyacophila" & Source!="Baetis" & Source!="Stegopterna")

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

#carcass/b/a community data matrix
HC_16S_com_CBF<-as.matrix(HC_16S_OTU_map_CBF[,15:ncol(HC_16S_OTU_map_CBF)])

#Create overall environmental data matrix for community analysis
HC_16S_env<-(HC_16S_OTU_map[,1:14])

#biofilm environmental data matrix
HC_16S_env_BF<-(HC_16S_OTU_map_BF_v[,1:14])

#carcass environmental data matrix
HC_16S_env_C<-(HC_16S_OTU_map_C[,1:14])

#Insect environmental data matrix
HC_16S_env_I<-(HC_16S_OTU_map_I[,1:14])

#carcassba env data matrix
HC_16S_env_CBF<-(HC_16S_OTU_map_CBF[,1:14])

#Make sample type variable for overall analysis
Sample.type<-factor(HC_16S_env$Source)

#Make reach variable for analysis
stream_site=factor(HC_16S_env$Reach, c("Salmon", "Control"))

#Make year variable for analysis
Study_Year<-factor(HC_16S_env$Year, c("1", "2"))

#Make days since carcass introduction variable for analysis and combine similar values for simplicity
HC_16S_env$Total_Biofilm_Growth_PostCarcass<-as.factor(HC_16S_env$Total_Biofilm_Growth_PostCarcass)
levels(HC_16S_env$Total_Biofilm_Growth_PostCarcass)
HC_16S_env$Total_Biofilm_Growth_PostCarcass<-as.factor(gsub("16","14",HC_16S_env$Total_Biofilm_Growth_PostCarcass))
HC_16S_env$Total_Biofilm_Growth_PostCarcass<-as.factor(gsub("141","142",HC_16S_env$Total_Biofilm_Growth_PostCarcass))
HC_16S_env$Total_Biofilm_Growth_PostCarcass<-as.factor(gsub("265","264",HC_16S_env$Total_Biofilm_Growth_PostCarcass))
Days_Carcass<-factor(HC_16S_env$Total_Biofilm_Growth_PostCarcass, c("0", "14", "162", "222", "264", "306"))

#Make stream_site variable for biofilm analysis
stream_site_bf=factor(HC_16S_env_BF$Reach, c("Salmon", "Control"))

#Make Year variable for biofilm analysis
Study_Year_bf<-factor(HC_16S_env_BF$Year, c("1", "2"))

#variable carcass before after 
HC_16S_env_CBF$before_carcass_after<-paste(HC_16S_env_CBF$Source, HC_16S_env_CBF$Total_Biofilm_Growth_PostCarcass)
HC_16S_env_CBF$before_carcass_after<-gsub("Biofilm 14", "Biofilm After", HC_16S_env_CBF$before_carcass_after)
HC_16S_env_CBF$before_carcass_after<-gsub("Biofilm 16", "Biofilm After", HC_16S_env_CBF$before_carcass_after)
HC_16S_env_CBF$before_carcass_after<-gsub("Biofilm 0", "Biofilm Before", HC_16S_env_CBF$before_carcass_after)
HC_16S_env_CBF$before_carcass_after<-gsub("Carcass 0", "Carcass", HC_16S_env_CBF$before_carcass_after)
carcass<-factor(HC_16S_env_CBF$before_carcass_after)

#variable for insect type
Taxa_I<-factor(HC_16S_env_I$Source, c("Rhyacophila", "Stegopterna", "Baetis"))

#variable for stream site
stream_site_I<-factor(HC_16S_env_I$Reach, c("Salmon", "Control"))

#Make days since carcass introduction variable for analysis and combine similar values for simplicity biofilm
HC_16S_env_BF$Total_Biofilm_Growth_PostCarcass<-as.factor(HC_16S_env_BF$Total_Biofilm_Growth_PostCarcass)
levels(HC_16S_env_BF$Total_Biofilm_Growth_PostCarcass)
HC_16S_env_BF$Total_Biofilm_Growth_PostCarcass<-as.factor(gsub("16","14",HC_16S_env_BF$Total_Biofilm_Growth_PostCarcass))
HC_16S_env_BF$Total_Biofilm_Growth_PostCarcass<-as.factor(gsub("141","142",HC_16S_env_BF$Total_Biofilm_Growth_PostCarcass))
HC_16S_env_BF$Total_Biofilm_Growth_PostCarcass<-as.factor(gsub("265","264",HC_16S_env_BF$Total_Biofilm_Growth_PostCarcass))
Days_Carcass_BF<-factor(HC_16S_env_BF$Total_Biofilm_Growth_PostCarcass, c("0", "14", "142", "222", "264", "306"))

#Overall NMDS for all 16S samples
HC_16S_NMDS<-metaMDS(HC_16S_com, distance="bray")

#stressplot for overall nmds
stressplot(HC_16S_NMDS)

#make figure for overall 16S nmds by sample type (biofilm, carcass, insect) with ellipses
ordiplot(HC_16S_NMDS, type="n", main="16S")
with(HC_16S_NMDS, points(HC_16S_NMDS, display="sites", col=five_col_vec[Sample.type], pch=19, pt.bg=five_col_vec))
with(HC_16S_NMDS, legend("topleft", legend=levels(Sample.type), bty="n", col=five_col_vec, pch=19, pt.bg=two_col_vec, cex=0.75))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Sample.type, kind="se", conf=0.95, lwd=2, col="#7fc97f", show.groups = "Baetis"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Sample.type, kind="se", conf=0.95, lwd=2, col="#beaed4", show.groups = "Biofilm"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Sample.type, kind="se", conf=0.95, lwd=2, col="#fdc086", show.groups = "Carcass"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Sample.type, kind="se", conf=0.95, lwd=2, col="#ffff99", show.groups = "Rhyacophila"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Sample.type, kind="se", conf=0.95, lwd=2, col="#386cb0", show.groups = "Stegopterna"))

#make figure where different colors are different sample types and ellipses are different reaches
ordiplot(HC_16S_NMDS, type="n", main="16S")
with(HC_16S_NMDS, points(HC_16S_NMDS, display="sites", col=five_col_vec[Sample.type], pch=19, pt.bg=five_col_vec))
with(HC_16S_NMDS, legend("topleft", legend=levels(Sample.type), bty="n", col=five_col_vec, pch=19, pt.bg=two_col_vec, cex=0.75))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, stream_site, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Salmon"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, stream_site, kind="se", conf=0.95, lwd=2, col="bisque3", show.groups = "Control"))

#make figure where different colors are different sample types and ellipses are different years
ordiplot(HC_16S_NMDS, type="n", main="16S")
with(HC_16S_NMDS, points(HC_16S_NMDS, display="sites", col=five_col_vec[Sample.type], pch=19, pt.bg=five_col_vec))
with(HC_16S_NMDS, legend("topleft", legend=levels(Sample.type), bty="n", col=five_col_vec, pch=19, pt.bg=two_col_vec, cex=0.75))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Study_Year, kind="se", conf=0.95, lwd=2, col="black", show.groups = "1"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Study_Year, kind="se", conf=0.95, lwd=2, col="bisque", show.groups = "2"))

#make figure where different colors are different sample types and ellipses are different months
ordiplot(HC_16S_NMDS, type="n", main="16S")
with(HC_16S_NMDS, points(HC_16S_NMDS, display="sites", col=five_col_vec[Sample.type], pch=19, pt.bg=five_col_vec))
with(HC_16S_NMDS, legend("topleft", legend=levels(Sample.type), bty="n", col=five_col_vec, pch=19, pt.bg=two_col_vec, cex=0.75))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="black", show.groups = "0"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="black", show.groups = "14"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="black", show.groups = "162"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="black", show.groups = "222"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="black", show.groups = "264"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="black", show.groups = "306"))

#Overal nmds by year of study
ordiplot(HC_16S_NMDS, type="n", main="16S by year")
with(HC_16S_NMDS, points(HC_16S_NMDS, display="sites", col=two_col_vec[Study_Year], pch=19, pt.bg=two_col_vec))
with(HC_16S_NMDS, legend("topleft", legend=levels(Study_Year), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))

#overall nmds by time since carcass introduction
levels(Days_Carcass)
ordiplot(HC_16S_NMDS, type="n", main="16S by time")
with(HC_16S_NMDS, points(HC_16S_NMDS, display="sites", col=six_col_vec_cont[Days_Carcass], pch=19, pt.bg=six_col_vec_cont))
with(HC_16S_NMDS, legend("topleft", legend=levels(Days_Carcass), bty="n", col=six_col_vec_cont, pch=19, pt.bg=six_col_vec_cont, cex=0.75))

#overall nmds by reach
ordiplot(HC_16S_NMDS, type="n", main="16S by site")
with(HC_16S_NMDS, points(HC_16S_NMDS, display="sites", col=two_col_vec[stream_site], pch=19, pt.bg=two_col_vec))
with(HC_16S_NMDS, legend("topleft", legend=levels(stream_site), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))

#Overall permanova
adonis(HC_16S_com ~ Reach + Month + Year + Source + Reach*Month + Reach*Year + Reach+Source + Month*Year + Month*Source + Year*Source, data=HC_16S_env, permutations=999)

#make NMDS Biofilms
HC_16S_BF_NMDS<-metaMDS(HC_16S_com_BF, distance="bray")

#Stressplot for biofilm nmds
stressplot(HC_16S_BF_NMDS)

#Biofilms With site points by reach
ordiplot(HC_16S_BF_NMDS, type="n", main="16S Biofilms")
with(HC_16S_BF_NMDS, points(HC_16S_BF_NMDS, display="sites", col=two_col_vec[stream_site_bf], pch=19, pt.bg=two_col_vec))
with(HC_16S_BF_NMDS, legend("topleft", legend=levels(stream_site_bf), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))
with(HC_16S_BF_NMDS, ordiellipse(HC_16S_BF_NMDS, stream_site_bf, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Salmon"))
with(HC_16S_BF_NMDS, ordiellipse(HC_16S_BF_NMDS, stream_site_bf, kind="se", conf=0.95, lwd=2, col="bisque3", show.groups = "Control"))

#Biofilms with site points by year
ordiplot(HC_16S_BF_NMDS, type="n", main="16S Biofilms Year")
with(HC_16S_BF_NMDS, points(HC_16S_BF_NMDS, display="sites", col=two_col_vec[Study_Year], pch=19, pt.bg=two_col_vec))
with(HC_16S_BF_NMDS, legend("topleft", legend=levels(Study_Year), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))
with(HC_16S_BF_NMDS, ordiellipse(HC_16S_BF_NMDS, Study_Year, kind="se", conf=0.95, lwd=2, col="black", show.groups = "1"))
with(HC_16S_BF_NMDS, ordiellipse(HC_16S_BF_NMDS, Study_Year, kind="se", conf=0.95, lwd=2, col="bisque3", show.groups = "2"))

#biofilms by month
ordiplot(HC_16S_BF_NMDS, type="n", main="16S Biofilms Month")
with(HC_16S_BF_NMDS, points(HC_16S_BF_NMDS, display="sites", col=six_col_vec[Days_Carcass_BF], pch=19, pt.bg=six_col_vec))
with(HC_16S_BF_NMDS, legend("topleft", legend=levels(Days_Carcass_BF), bty="n", col=six_col_vec, pch=19, pt.bg=six_col_vec))

#biofilms permanova
adonis(HC_16S_com_BF ~ Reach + Date + Reach*Date, data=HC_16S_env_BF, permutations=999)

#make nmds for insects
HC_16S_I_NMDS<-metaMDS(HC_16S_com_I, distance="bray")

#Visualize reach for insects
ordiplot(HC_16S_I_NMDS, type="n", main="16S Insect")
with(HC_16S_I_NMDS, points(HC_16S_I_NMDS, display="sites", col=two_col_vec[stream_site_I], pch=19, pt.bg=two_col_vec))

#permanova insects
adonis(HC_16S_com_I ~ stream_site_I*Taxa_I, data=HC_16S_env_I, permutations=999)

#nmds carcass before and after biofilms
HC_16S_CBF_NMDS<-metaMDS(HC_16S_com_CBF, distance="bray")

#nmds plot carcass vs before vs after
ordiplot(HC_16S_CBF_NMDS, type="n", main="16S Carcass and Biofilm")
with(HC_16S_CBF_NMDS, points(HC_16S_CBF_NMDS, display="sites", col=three_col_vec[carcass], pch=19, pt.bg=three_col_vec))
with(HC_16S_CBF_NMDS, legend("topleft", legend=levels(carcass), bty="n", col=three_col_vec, pch=19, pt.bg=three_col_vec))
with(HC_16S_CBF_NMDS, ordiellipse(HC_16S_CBF_NMDS, carcass, kind="se", conf=0.95, lwd=2, col="#a6cee3", show.groups = "Biofilm After"))
with(HC_16S_CBF_NMDS, ordiellipse(HC_16S_CBF_NMDS, carcass, kind="se", conf=0.95, lwd=2, col="#1f78b4", show.groups = "Biofilm Before"))
with(HC_16S_CBF_NMDS, ordiellipse(HC_16S_CBF_NMDS, carcass, kind="se", conf=0.95, lwd=2, col="#b2df8a", show.groups = "Carcass"))

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

#limit phyla data to biofilms
HC_16S_P_map_BF<-subset(HC_16S_P_map, Source == "Biofilm")

#Make line plot for cyanobacteria for salmon vs control
x<-as.vector((HC_16S_P_map_BF$Total_Biofilm_Growth))
sum(HC_16S_P_map_BF$"k__Bacteria.p__Cyanobacteria")
HC_16S_P_map_BF$Cyanobacteria<-as.numeric((HC_16S_P_map_BF$"k__Bacteria.p__Cyanobacteria")/sum(HC_16S_P_map_BF$"k__Bacteria.p__Cyanobacteria"))
HC_16S_P_map_BF$Reach<-as.factor(HC_16S_P_map_BF$Reach)
Cyo_p<-data.frame(HC_16S_P_map_BF[,c(8,11,89)])
tgc <- summarySE(Cyo_p, measurevar="Cyanobacteria", groupvars=c("Total_Biofilm_Growth","Reach"))
gd <- data.frame(Cyo_p %>% 
  group_by(Reach, Total_Biofilm_Growth) %>% 
  summarise(Cyanobacteria = mean(Cyanobacteria)))
ggplot(tgc, aes(x=Total_Biofilm_Growth, y=Cyanobacteria, colour=Reach)) + 
  geom_errorbar(aes(ymin=Cyanobacteria-se, ymax=Cyanobacteria+se), width=.1) +
  geom_line() +
  geom_point()
  
#melt to long format
names(HC_16S_P_map[,1:14])
str(HC_16S_P_map[,1:14])
HC_16S_P_map$Row.names<-as.factor(HC_16S_P_map$Row.names)
HC_16S_P_map$Month<-as.factor(HC_16S_P_map$Month)
HC_16S_P_map$Day<-as.factor(HC_16S_P_map$Day)
HC_16S_P_map$YYYY<-as.factor(HC_16S_P_map$YYYY)
HC_16S_P_map$Total_Biofilm_Growth<-as.factor(HC_16S_P_map$Total_Biofilm_Growth)
HC_16S_P_map$Total_Biofilm_Growth_PostCarcass<-as.factor(HC_16S_P_map$Total_Biofilm_Growth_PostCarcass)
HC_16S_P_map$Year<-as.factor(HC_16S_P_map$Year)
HC_16S_P_map_m<-melt(HC_16S_P_map, 
    id.vars=c("Row.names", "BarcodeSequence", "LinkerPrimerSequence", "Month", "Day", "YYYY", "Date", "Reach", "Subreach", "Source", "Total_Biofilm_Growth", "Total_Biofilm_Growth_PostCarcass", "Year", "Description"))

#Delete zeros in melted file
HC_16S_P_map_m<-subset(HC_16S_P_map_m, value > 0)



