#Analysis of 16S microbial communities at Hunt Creek 2014-2016

#Load packages
library(reshape)
library(ggplot2)
library(vegan)

#Get tables from different runs - R1, R2, R3, R4, R5 and R6
Hunt_Creek_16S_R1<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Taxonomy_txt_tables/HC_R1_table_tabseparated.txt", sep="\t", header = T)
row.names(Hunt_Creek_16S_R1)<-Hunt_Creek_16S_R1[,1]
Hunt_Creek_16S_R1<-Hunt_Creek_16S_R1[,-c(1,14)]

#transpose
Hunt_Creek_16S_R1_t<-t(Hunt_Creek_16S_R1)
Hunt_Creek_16S_R1_t<-data.frame(Hunt_Creek_16S_R1_t)
str(Hunt_Creek_16S_R1_t)
names(Hunt_Creek_16S_R1_t)

#Get metadata through mapping file
Hunt_Creek_16S_R1_map <- read.table("/Users/courtneylarson/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/Mapping_files/Hunt_Creek_Map_R1.txt", header=T)
row.names(Hunt_Creek_16S_R1_map)<-Hunt_Creek_16S_R1_map[,1]
Hunt_Creek_16S_R1_map<-Hunt_Creek_16S_R1_map[,-c(1)]

#Merge metadata onto data table
HC_16S_R1_OTU_map <-merge(Hunt_Creek_16S_R1_map, Hunt_Creek_16S_R1_t, by=0)
HC_16S_R1_OTU_map_m<-melt(HC_16S_R1_OTU_map, id=c("Row.names", "BarcodeSequence", "LinkerPrimerSequence", "Month", "Day", "YYYY", "Date", "Reach", "Subreach", "Source", "Total_Biofilm_Growth", "Total_Biofilm_Growth_PostCarcass", "Year", "Description"))
HC_16S_R1_OTU_map_m<-data.frame(HC_16S_R1_OTU_map_m)
names(HC_16S_R1_OTU_map_m)

#Delete zeros
HC_16S_R1_OTU_map_m<-subset(HC_16S_R1_OTU_map_m, value > 0)

#Reduce to top taxa
totals<-rbind(Hunt_Creek_16S_R1_t, colSums(Hunt_Creek_16S_R1_t))
totals<-totals[-c(1:12),]
sort(totals,decreasing=TRUE)[1:6]
#Use names given as output to reduce data frame for visualization
HC_16S_R1_OTU_map_m_sub<-subset(HC_16S_R1_OTU_map_m, variable=="denovo20259" | variable=="denovo5011" | variable=="denovo20477" | variable=="denovo17186" | variable=="denovo17019" | variable=="denovo19917")
HC_16S_R1_OTU_map_m_sub<-data.frame(HC_16S_R1_OTU_map_m_sub)
str(HC_16S_R1_OTU_map_m_sub)
HC_16S_R1_OTU_map_m_sub$Row.names<-as.factor(HC_16S_R1_OTU_map_m_sub$Row.names)
HC_16S_R1_OTU_map_m_sub$Total_Biofilm_Growth_PostCarcass<-as.factor(HC_16S_R1_OTU_map_m_sub$Total_Biofilm_Growth_PostCarcass)

#color vectors
two_col_vec <- c("black", "bisque3")

#make stacked bar graphs
ggplot(HC_16S_R1_OTU_map_m_sub, aes(x=HC_16S_R1_OTU_map_m_sub[,1], y=value, fill=variable)) + geom_bar(position = "fill",stat = "identity")
ggplot(HC_16S_R1_OTU_map_m_sub, aes(x=HC_16S_R1_OTU_map_m_sub[,12], y=value, fill=variable)) + geom_bar(position = "fill",stat = "identity")

#create community data matrix for community analysis
HC_16S_com<-as.matrix(HC_16S_R1_OTU_map[,15:ncol(HC_16S_R1_OTU_map)])
                      
#Create environmental data matrix for community analysis
HC_16S_env<-(HC_16S_R1_OTU_map[,1:14])

#Make stream_site variable
stream_site=factor(HC_16S_env$Reach, c("Salmon", "Control"))

#make NMDS
HC_16S_NMDS<-metaMDS(HC_16S_com, distance="bray")
ordiplot(HC_16S_NMDS, type="n", main="16S")

#With site points by reach
with(HC_16S_NMDS, points(HC_16S_NMDS, display="sites", col=two_col_vec[stream_site], pch=19, pt.bg=two_col_vec))
with(HC_16S_NMDS, legend("topleft", legend=levels(stream_site), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, stream_site, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Salmon"))
with(HC_16S_NMDS, ordiellipse(HC_16S_NMDS, stream_site, kind="se", conf=0.95, lwd=2, col="bisque3", show.groups = "Control"))




