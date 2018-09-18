#Analysis of 16S microbial communities at Hunt Creek 2014-2016

#########################################
#Load packages
library(reshape)
library(ggplot2)
library(vegan)
library(plyr)
library(dplyr)
library(indicspecies)
library(doBy)
library(ggpubr)
library(lme4)
#Functions
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence 
## interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#color vectors
two_col_vec <- c("black", "bisque3")
two_col_vec_reach<-c("skyblue3", "tomato3")
two_col_vec_taxa<-c("#33a02c", "#1f78b4")
three_col_vec<- c("#1b9e77", "#d95f02", "#7570b3")
four_col_vec<-c("#7fc97f", "#beaed4", "#fdc086", "#ffff99")
five_col_vec<- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")
six_col_vec<- c("#e41a1c", "#377eb8", "green", "#984ea3", "#ff7f00", "#ffff33")
six_col_vec_cont<-c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")
seven_col_vec<-c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")

#First work with rarefication plots to use with presentation
HC_Shannon<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_shannon.txt", sep="\t", header=T)
HC_Shannon$X<-NULL
HC_Shannon$iteration<-NULL
HC_Shannon$sequences.per.sample<-as.double(HC_Shannon$sequences.per.sample)
HC_Sh_m<-melt(HC_Shannon, id.vars="sequences.per.sample")
HC_Sh_c_m<-cast(HC_Sh_m, variable ~ sequences.per.sample, fun.aggregate=mean)
HC_Sh_c_m$Calculation<-rep("mean",172)
HC_Sh_c_var<-cast(HC_Sh_m, variable ~ sequences.per.sample, fun.aggregate=var)
HC_Sh_c_var$Calculation<-rep("variance",172)
HC_Sh_c_m_var<-rbind(HC_Sh_c_m,HC_Sh_c_var)
names(HC_Sh_c_m_var)[names(HC_Sh_c_m_var)=="variable"] <- "SampleID"
HC_Sh_c_m_var$Calculation<-as.factor(HC_Sh_c_m_var$Calculation)
HC_Sh_c_m_var<-as.data.frame(HC_Sh_c_m_var)
HC_Sh_c_m_var_m<-melt(HC_Sh_c_m_var, id.vars=c("Calculation","SampleID"))
HC_Sh_c_m_var_c<-cast(HC_Sh_c_m_var_m, variable + SampleID ~ Calculation)
#Get metadata through mapping file
Hunt_Creek_16S_map <- read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/Hunt_Creek_Map_Filtered_R.txt", header=T)
#Merge metadata onto rarefication file
HC_Sh_map <-merge(Hunt_Creek_16S_map, HC_Sh_c_m_var_c, by="SampleID")
names(HC_Sh_map)[names(HC_Sh_map)=="variable"] <- "Sequences_per_sample"
HC_Sh_map$Sequences_per_sample<-as.numeric(as.character(HC_Sh_map$Sequences_per_sample))
HC_Sh_map$Source<-revalue(HC_Sh_map$Source, c("Baetis"="Insect", "Heptagenia"="Insect", "Rhyacophila"="Insect", "Simulium"="Insect", "Stegopterna"="Insect"))
HC_Sh_map_sum_m <- summarySE(HC_Sh_map, measurevar=c("mean"), groupvars=c("Sequences_per_sample","Source"))
HC_Sh_map_sum_v <- summarySE(HC_Sh_map, measurevar=c("variance"), groupvars=c("Sequences_per_sample","Source"))
HC_Sh_map_sum_v$StandDev<-sqrt(HC_Sh_map_sum_v$variance)
HC_Sh_map_sum_v$StandEr<-HC_Sh_map_sum_v$StandDev/sqrt(HC_Sh_map_sum_v$N)
HC_Sh_sum_m_sd<-merge(HC_Sh_map_sum_m,HC_Sh_map_sum_v, by=0)
#make rarefication plots
ggplot(HC_Sh_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Sequences per sample") +
  ylab("Shannon H' Diversity") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=three_col_vec)

#####################################
#Upload data tables generated in QIIME and manipulate to make "r friendly"
#####################################

#Get 16S OTU table
Hunt_Creek_16S<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_table_tabseparated_R.txt", sep="\t", header = T)

#Format data frame so the denovo is row name
row.names(Hunt_Creek_16S)<-Hunt_Creek_16S[,1]

#Delete taxonomy and denovo columns, now that denovo is row name
Hunt_Creek_16S$ID<-NULL
Hunt_Creek_16S$Consensus.Lineage<-NULL
#transpose
HC_16S_OTU_t<-t(Hunt_Creek_16S)
HC_16S_OTU_t<-data.frame(HC_16S_OTU_t)

#Merge metadata onto data table
Hunt_Creek_16S_map_r<-Hunt_Creek_16S_map
row.names(Hunt_Creek_16S_map_r)<-Hunt_Creek_16S_map_r[,1]
HC_16S_OTU_map <-merge(Hunt_Creek_16S_map_r, HC_16S_OTU_t, by=0)
str(HC_16S_OTU_map)

#UNI-Do the following if you have unifrac distances ready
#UNI-Upload weighted unifrac distance matrix and remove samples with error
Hunt_Creek_16S_uni <- read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/weighted_unifrac_otu_table_HC.txt", header=T)
#UNI-Add metadata to unifrac table
Hunt_Creek_16S_uni_map <-merge(Hunt_Creek_16S_map_r, Hunt_Creek_16S_uni, by=0)
row.names(Hunt_Creek_16S_uni_map)<-Hunt_Creek_16S_uni_map[,1]
Hunt_Creek_16S_uni_map<-Hunt_Creek_16S_uni_map[,-c(1)]

#Create overall environmental data matrix for community analysis with counts
HC_16S_env<-Hunt_Creek_16S_map
row.names(HC_16S_env)<-HC_16S_env[,1]
HC_16S_env<-HC_16S_env[,-c(1)]

#UNI-Create overall environmental data matrix for community analysis with uni distances
HC_16S_uni_env<-Hunt_Creek_16S_uni_map
row.names(HC_16S_uni_env)<-HC_16S_uni_env[,1]
HC_16S_uni_env<-HC_16S_uni_env[,1:16]
HC_16S_uni_env$Total_Biofilm_Growth_PostCarcass<-as.numeric(HC_16S_uni_env$Total_Biofilm_Growth_PostCarcass)
HC_16S_uni_env$Year<-as.factor(HC_16S_uni_env$Year)

#merge insects into one source type
HC_16S_env_i1<-HC_16S_env
HC_16S_env_i1$Source<-gsub("Rhyacophila", "Insect", HC_16S_env_i1$Source)
HC_16S_env_i1$Source<-gsub("Stegopterna", "Insect", HC_16S_env_i1$Source)
HC_16S_env_i1$Source<-gsub("Baetis", "Insect", HC_16S_env_i1$Source)
HC_16S_env_i1$Source<-gsub("Simulium", "Insect", HC_16S_env_i1$Source)
HC_16S_env_i1$Source<-as.factor(gsub("Heptagenia", "Insect", HC_16S_env_i1$Source))
levels(HC_16S_env_i1$Source)
HC_16S_env_i1$Total_Biofilm_Growth_PostCarcass<-as.factor(HC_16S_env_i1$Total_Biofilm_Growth_PostCarcass)
HC_16S_env_i1$Year<-as.factor(HC_16S_env_i1$Year)
HC_16S_env_i1$Reach<-as.factor(HC_16S_env_i1$Reach)

#UNI-merge insects into one source type in unifrac table
HC_16S_uni_env_i1<-HC_16S_uni_env
HC_16S_uni_env_i1$Source<-gsub("Rhyacophila", "Insect", HC_16S_uni_env_i1$Source)
HC_16S_uni_env_i1$Source<-gsub("Stegopterna", "Insect", HC_16S_uni_env_i1$Source)
HC_16S_uni_env_i1$Source<-gsub("Heptagenia", "Insect", HC_16S_uni_env_i1$Source)
HC_16S_uni_env_i1$Source<-gsub("Simulium", "Insect", HC_16S_uni_env_i1$Source)
HC_16S_uni_env_i1$Source<-as.factor(gsub("Baetis", "Insect", HC_16S_uni_env_i1$Source))
levels(HC_16S_uni_env_i1$Source)
HC_16S_uni_env_i1$Total_Biofilm_Growth_PostCarcass<-as.factor(HC_16S_uni_env_i1$Total_Biofilm_Growth_PostCarcass)
HC_16S_uni_env_i1$Year<-as.factor(HC_16S_uni_env_i1$Year)
HC_16S_uni_env_i1$Reach<-as.factor(HC_16S_uni_env_i1$Reach)
#############################################
#Analysis of all 16S community data
#############################################

#Overall permanova
adonis(HC_16S_OTU_t ~ Reach*Year*Source*Total_Biofilm_Growth_PostCarcass, data=HC_16S_env_i1, method="jaccard", permutations=999)
#Significant factors are Year, Source, Totay_Biofilm_Growth_PostCarcass, Year:Source and Reach:Source:Total_Biofilm_Growth_PostCarcass

#UNI-Overall permanova with unifrac distances
adonis(as.dist(Hunt_Creek_16S_uni) ~ Reach*Source*Total_Biofilm_Growth_PostCarcass*Year, data=HC_16S_uni_env, permutations=999)
#Source has the largest impact on community structuring, so do analysis for each source individually

#UNI-Create carcass unifrac distance table with metadata
H_C_16S_uni_map_car<-subset(Hunt_Creek_16S_uni_map, Source=="Carcass")
H_C_16S_uni_car_env<-H_C_16S_uni_map_car[,1:16]
HC_16S_uni_car_Year<-as.factor(H_C_16S_uni_car_env$Year)
HC_16S_uni_car_Year<-gsub("1", "Year 1", HC_16S_uni_car_Year)
HC_16S_uni_car_Year<-gsub("2", "Year 2", HC_16S_uni_car_Year)
HC_16S_uni_car_Year<-as.factor(HC_16S_uni_car_Year)
H_C_16S_uni_car_com<-H_C_16S_uni_map_car[,17:ncol(H_C_16S_uni_map_car)]
H_C_16S_uni_car_com_samples<-as.vector(rownames(H_C_16S_uni_car_com))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_uni_car_com<-as.matrix(H_C_16S_uni_car_com)
H_C_16S_uni_car_com<-subset(H_C_16S_uni_car_com, select=c(H_C_16S_uni_car_com_samples))
#Create carcass otu table with metadata
HC_16S_OTU_map_car<-subset(HC_16S_OTU_map, Source=="Carcass")
#How many OTUs and gene amplicon sequences for carcasses?
#Delete OTUs with no observations
cols_to_drop_car = c(rep(TRUE, 17), colSums(HC_16S_OTU_map_car[,18:ncol(HC_16S_OTU_map_car)]) > 0)
HC_16S_OTU_map_car<-HC_16S_OTU_map_car[,cols_to_drop_car]
#There are 6158 OTUs (number of variables - 17 metadata variables)
sum(colSums(HC_16S_OTU_map_car[,18:ncol(HC_16S_OTU_map_car)]))
#45,600 total sequence reads
#Subset into environmental and community data frames
HC_16S_OTU_car_env<-HC_16S_OTU_map_car[,1:17]
HC_16S_OTU_car_com<-HC_16S_OTU_map_car[,18:ncol(HC_16S_OTU_map_car)]
#Year 1
HC_16S_OTU_map_car_Y1<-subset(HC_16S_OTU_map_car, Year=="1")
#Delete OTUs with no observations
cols_to_drop_car_Y1 = c(rep(TRUE, 17), colSums(HC_16S_OTU_map_car_Y1[,18:ncol(HC_16S_OTU_map_car_Y1)]) > 0)
HC_16S_OTU_map_car_Y1<-HC_16S_OTU_map_car_Y1[,cols_to_drop_car_Y1]
#There are 4854 OTUs in year 1(number of variables - 17 metadata variables)
sum(colSums(HC_16S_OTU_map_car_Y1[,18:ncol(HC_16S_OTU_map_car_Y1)]))
#21,600 total sequence reads
#Year 2
HC_16S_OTU_map_car_Y2<-subset(HC_16S_OTU_map_car, Year=="2")
#Delete OTUs with no observations
cols_to_drop_car_Y2 = c(rep(TRUE, 17), colSums(HC_16S_OTU_map_car_Y2[,18:ncol(HC_16S_OTU_map_car_Y2)]) > 0)
HC_16S_OTU_map_car_Y2<-HC_16S_OTU_map_car_Y2[,cols_to_drop_car_Y2]
#There are 1894 OTUs in year 2(number of variables - 17 metadata variables)
sum(colSums(HC_16S_OTU_map_car_Y2[,18:ncol(HC_16S_OTU_map_car_Y2)]))
#24,000 total sequence reads

#UNI-carcass permanova with unifrac distances
adonis(as.dist(H_C_16S_uni_car_com) ~ Year, data=H_C_16S_uni_car_env, permutations=999)
#NMDS carcass
HC_NMDS_uni_c<-metaMDS(as.dist(H_C_16S_uni_car_com))
HC_NMDS_uni_c
stressplot(HC_NMDS_uni_c)
ordiplot(HC_NMDS_uni_c, type="n", main="Salmon carcass 16S communities")
with(HC_NMDS_uni_c, points(HC_NMDS_uni_c, display="sites", col=two_col_vec[HC_16S_uni_car_Year], pch=19, pt.bg=two_col_vec))
with(HC_NMDS_uni_c, legend("topleft", legend=levels(HC_16S_uni_car_Year), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))
with(HC_NMDS_uni_c, ordiellipse(HC_NMDS_uni_c, HC_16S_uni_car_Year, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Year 1"))
with(HC_NMDS_uni_c, ordiellipse(HC_NMDS_uni_c, HC_16S_uni_car_Year, kind="se", conf=0.95, lwd=2, col="bisque2", show.groups = "Year 2"))
with(HC_NMDS_uni_c, legend("topright", legend="2D Stress: 0.09", bty="n"))

#upload phyla level info and run indicator analysis for carcass/generate box plots
#Upload phyla level files for each run
Hunt_Creek_16S_P<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_Phyla_16S.txt", sep="\t", header = T)
#Clasify as data.frame
Hunt_Creek_16S_P<-data.frame(Hunt_Creek_16S_P)
#Format data frame so the taxonomy is row name
row.names(Hunt_Creek_16S_P)<-Hunt_Creek_16S_P[,1]
#Delete taxonomy column
Hunt_Creek_16S_P$ID<-NULL
#transpose
HC_16S_P_t<-t(Hunt_Creek_16S_P)
HC_16S_P_t<-data.frame(HC_16S_P_t)
str(HC_16S_P_t)
names(HC_16S_P_t)
#Merge metadata onto data table
HC_16S_P_map <-merge(Hunt_Creek_16S_map_r, HC_16S_P_t, by=0)

#subset phyla info into carcass
HC_16S_P_map_car<-subset(HC_16S_P_map, Source=="Carcass")
#determine how many phyla
cols_to_drop_car_P = c(rep(TRUE, 17), colSums(HC_16S_P_map_car[,18:ncol(HC_16S_P_map_car)]) > 0)
HC_16S_P_map_car<-HC_16S_P_map_car[,cols_to_drop_car_P]
#39 total phyla (number of variables - 17 metadata columns)
HC_16S_P_map_car[,1:17]<-sapply(HC_16S_P_map_car[,1:17], as.factor)
HC_16S_P_map_car_env<-HC_16S_P_map_car[,1:17]
HC_16S_P_map_car_com<-HC_16S_P_map_car[,18:ncol(HC_16S_P_map_car)]
#find most abundant phyla
sort(colSums(HC_16S_P_map_car_com),decreasing=TRUE)
#Proteobacteria most abundant
0.9*sum(colSums(HC_16S_P_map_car_com))
#90% of sequence reads are 41040, so top 7 pyla represent 90% of sample (Proteobactera, Bacteroidetes, Firmicutes, Actinobacteria, Acidobacteria, Spirochaetes and Cyanobacteria)
#Indicator analysis to see what phyla are driving this change
HC_p_car_indic<-multipatt(HC_16S_P_map_car_com, HC_16S_P_map_car_env$Year, control = how(nperm=999))
summary(HC_p_car_indic)
#10 phyla with year 1 significant: Spirochaetes, euryarchaota, nitrospirae, tenericutes, chlorobi, armatimonadetes, crenarchaota, elusimicrobia, TM7,OP11

#make stacked bar graphs of phyla for each carcass year
#Create "Other" variable
HC_16S_P_map_car$Other<-rowSums(HC_16S_P_map_car[,-c(1:17,23:24,27,32,38,48:49)])
#simplify and melt dataset
HC_16S_P_map_car_m<-melt(HC_16S_P_map_car)
HC_16S_P_map_car_m<-subset(HC_16S_P_map_car_m, variable=="k__Bacteria.p__Acidobacteria" | variable=="k__Bacteria.p__Actinobacteria" | variable=="k__Bacteria.p__Bacteroidetes" | variable=="k__Bacteria.p__Cyanobacteria" | variable=="k__Bacteria.p__Firmicutes" | variable=="k__Bacteria.p__Proteobacteria"| variable=="k__Bacteria.p__Spirochaetes"| variable=="Other")
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Acidobacteria","Acidobacteria",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Actinobacteria","Actinobacteria",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Bacteroidetes","Bacteroidetes",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Cyanobacteria","Cyanobacteria",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Firmicutes","Firmicutes",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Proteobacteria","Proteobacteria",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Spirochaetes","Spirochaetes",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("Other","Rare Taxa <3%",HC_16S_P_map_car_m$variable)

ggplot(HC_16S_P_map_car_m,aes(x = Year, y = value,fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") +
  xlab("Carcass Year") +
  ylab("Relative abundance") +
  labs(fill = "Phyla") +
  scale_fill_brewer(palette="Dark2") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16))

#############################################
#Separate data based on year of study (1 or 2)
###############################################

#UNI-Create Year 1 unifrac distance table with metadata
H_C_16S_uni_map_Y1<-subset(Hunt_Creek_16S_uni_map, Year=="1")

#UNI-Create Year 2 unifrac distance table with metadata
H_C_16S_uni_map_Y2<-subset(Hunt_Creek_16S_uni_map, Year=="2")

#Create year 1 otu table with metadata
HC_16S_OTU_map_Y1<-subset(HC_16S_OTU_map, Year=="1")

#Create year2 otu table with metadata
HC_16S_OTU_map_Y2<-subset(HC_16S_OTU_map, Year=="2")

####################################################
#Biofilm sample analysis by year
############################################

#Create biofilm data tables by year 

#Biofilm Year 1

#Create biofilm matrix for year 1 with metadata
H_C_16S_map_Y1_BF<-subset(HC_16S_OTU_map_Y1, Source=="Biofilm")
#how many OTUs BF Y1?
#Delete OTUs with no observations
cols_to_drop_BF_Y1 = c(rep(TRUE, 17), colSums(H_C_16S_map_Y1_BF[,18:ncol(H_C_16S_map_Y1_BF)]) > 0)
H_C_16S_map_Y1_BF<-H_C_16S_map_Y1_BF[,cols_to_drop_BF_Y1]
#There are 11051 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_map_Y1_BF[,18:ncol(H_C_16S_map_Y1_BF)]))
#86,400 total sequence reads
#without metadata for biofilms in year 1
H_C_16S_Y1_BF<-H_C_16S_map_Y1_BF[,18:ncol(H_C_16S_map_Y1_BF)]
H_C_16S_Y1_BF_samples<-as.vector(rownames(H_C_16S_Y1_BF))
#Biofilm year 1 environmental variable table
H_C_16S_env_Y1_BF<-H_C_16S_map_Y1_BF[,1:17]
Total_Biofilm_Growth_PostCarcass_Y1<-as.factor(H_C_16S_env_Y1_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y1<-as.factor(H_C_16S_env_Y1_BF$Reach)

#UNI-Create biofilm unifrac distance table for year 1 with metadata
H_C_16S_uni_map_Y1_BF<-subset(H_C_16S_uni_map_Y1, Source=="Biofilm")
#UNI-unifrac distance table without metadata for biofilms in year 1
H_C_16S_uni_Y1_BF<-H_C_16S_uni_map_Y1_BF[,17:ncol(H_C_16S_uni_map_Y1_BF)]
H_C_16S_uni_Y1_BF_samples<-as.vector(rownames(H_C_16S_uni_Y1_BF))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_uni_Y1_BF<-as.matrix(H_C_16S_uni_Y1_BF)
H_C_16S_uni_Y1_BF<-subset(H_C_16S_uni_Y1_BF, select=c(H_C_16S_uni_Y1_BF_samples))
#UNI-Biofilm year 1 environmental variable table
H_C_16S_env_Y1_BF<-H_C_16S_uni_map_Y1_BF[,1:16]
Total_Biofilm_Growth_PostCarcass_Y1<-as.factor(H_C_16S_env_Y1_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y1<-as.factor(H_C_16S_env_Y1_BF$Reach)
BA_Y1<-as.factor(H_C_16S_env_Y1_BF$BA)

#Biofilm year 2
#Create biofilm matrix for year 2 with metadata
H_C_16S_map_Y2_BF<-subset(HC_16S_OTU_map_Y2, Source=="Biofilm")
#how many OTUs BF Y2?
#Delete OTUs with no observations
cols_to_drop_BF_Y2 = c(rep(TRUE, 17), colSums(H_C_16S_map_Y2_BF[,18:ncol(H_C_16S_map_Y2_BF)]) > 0)
H_C_16S_map_Y2_BF<-H_C_16S_map_Y2_BF[,cols_to_drop_BF_Y2]
#There are 9434 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_map_Y2_BF[,18:ncol(H_C_16S_map_Y2_BF)]))
#86,400 total sequence reads
#table without metadata for biofilms in year 2
H_C_16S_Y2_BF<-H_C_16S_map_Y2_BF[,18:ncol(H_C_16S_map_Y2_BF)]
H_C_16S_Y2_BF_samples<-as.vector(rownames(H_C_16S_Y2_BF))
#Biofilm year 2 environmental variable table
H_C_16S_env_Y2_BF<-H_C_16S_map_Y2_BF[,1:17]
Total_Biofilm_Growth_PostCarcass_Y2<-as.factor(H_C_16S_env_Y2_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y2<-as.factor(H_C_16S_env_Y2_BF$Reach)

#UNI-Create biofilm unifrac distance table for year 2 with metadata
H_C_16S_uni_map_Y2_BF<-subset(H_C_16S_uni_map_Y2, Source=="Biofilm")
#UNI_unifrac distance table without metadata for biofilms in year 2
H_C_16S_uni_Y2_BF<-H_C_16S_uni_map_Y2_BF[,17:ncol(H_C_16S_uni_map_Y2_BF)]
H_C_16S_uni_Y2_BF_samples<-as.vector(rownames(H_C_16S_uni_Y2_BF))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_uni_Y2_BF<-as.matrix(H_C_16S_uni_Y2_BF)
H_C_16S_uni_Y2_BF<-subset(H_C_16S_uni_Y2_BF, select=c(H_C_16S_uni_Y2_BF_samples))
#UNI-Biofilm year 2 environmental variable table
H_C_16S_env_Y2_BF<-H_C_16S_uni_map_Y2_BF[,1:16]
Total_Biofilm_Growth_PostCarcass_Y2<-as.factor(H_C_16S_env_Y2_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y2<-as.factor(H_C_16S_env_Y2_BF$Reach)
BA_Y2<-as.factor(H_C_16S_env_Y2_BF$BA)
#Year 1 biofilm permanova using counts
adonis(H_C_16S_Y1_BF ~ Reach_Y1*Total_Biofilm_Growth_PostCarcass_Y1, data=H_C_16S_env_Y1_BF, method="jaccard", permutations=999)

#Year 2 biofilm permanova using counts
adonis(H_C_16S_Y2_BF ~ Reach_Y2*Total_Biofilm_Growth_PostCarcass_Y2, data=H_C_16S_env_Y2_BF, method="jaccard", permutations=999)

#UNI-Year 1 biofilm permanova using unfrac distance
adonis(as.dist(H_C_16S_uni_Y1_BF) ~ Reach_Y1*Total_Biofilm_Growth_PostCarcass_Y1, data=H_C_16S_env_Y1_BF, permutations=999)

#UNI-Year 2 biofilm permanova using unifrac distance
adonis(as.dist(H_C_16S_uni_Y2_BF) ~ Reach_Y2*Total_Biofilm_Growth_PostCarcass_Y2, data=H_C_16S_env_Y2_BF, permutations=999)

#indicator taxa analysis for biofilms
#subset phyla info into biofilms for each year
HC_16S_P_map_bf_y1<-subset(HC_16S_P_map, Source=="Biofilm" & Year=="1")
HC_16S_P_map_bf_y1[,1:17]<-sapply(HC_16S_P_map_bf_y1[,1:17], as.factor)
HC_16S_P_map_bf_y1_env<-HC_16S_P_map_bf_y1[,1:17]
HC_16S_P_map_bf_y1_com<-HC_16S_P_map_bf_y1[,18:ncol(HC_16S_P_map_bf_y1)]
#find most abundant phyla
sort(colSums(HC_16S_P_map_bf_y1_com),decreasing=TRUE)
#Proteobacteria most abundant
0.9*sum(colSums(HC_16S_P_map_bf_y1_com))
#90% of sequence reads are 77760, so top 6 pyla represent 90% of sample (Proteobactera, Cyanobacteria, Bacteroidetes, Verrucomicrobia, Actinobacteria and Acidobacteria)

#subset phyla info into biofilms y2
HC_16S_P_map_bf_y2<-subset(HC_16S_P_map, Source=="Biofilm" & Year=="2")
HC_16S_P_map_bf_y2[,1:17]<-sapply(HC_16S_P_map_bf_y2[,1:17], as.factor)
HC_16S_P_map_bf_y2_env<-HC_16S_P_map_bf_y2[,1:17]
HC_16S_P_map_bf_y2_com<-HC_16S_P_map_bf_y2[,18:ncol(HC_16S_P_map_bf_y2)]
#find most abundant phyla
sort(colSums(HC_16S_P_map_bf_y2_com),decreasing=TRUE)
#Proteobacteria most abundant
0.9*sum(colSums(HC_16S_P_map_bf_y2_com))
#90% of sequence reads are 77760, so top 6 pyla represent 90% of sample (Proteobactera, Cyanobacteria, Bacteroidetes, Verrucomicrobia, Planctomycetes and Actinobacteria)

#create clusters based on time and reach
BF_Y1_Reach<-HC_16S_P_map_bf_y1_env$Reach
BF_Y2_Reach<-HC_16S_P_map_bf_y2_env$Reach
BF_Y1_Time<-HC_16S_P_map_bf_y1_env$Total_Biofilm_Growth_PostCarcass
BF_Y2_Time<-HC_16S_P_map_bf_y2_env$Total_Biofilm_Growth_PostCarcass
BF_Y1_RT<-paste(BF_Y1_Reach, BF_Y1_Time)
BF_Y2_RT<-paste(BF_Y2_Reach, BF_Y2_Time)
#Indicator analysis of biofilms year 1 to see what phyla are driving this change
HC_p_bf_y1_indic<-signassoc(HC_16S_P_map_bf_y1_com, cluster=BF_Y1_RT,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_p_bf_y1_indic[,c(ncol(HC_p_bf_y1_indic),(ncol(HC_p_bf_y1_indic)-1))]
#indicator analysis found 6 phyla: Euryarchaeota (8 (Salmon 14)), Actinobacteria (2 (control 14)), Cyanobacteria (8), Firmicutes (8), Gemmatimonadetes (5 (Control 264)), Spirochaetes (8)
HC_p_bf_y1_indic
#Four indicate salmon reaches 2 weeks after introduction: Euryarchaeota, Cyanobacteria, Firmicutes, and Spirochaetes
#Indicator analysis of biofilms year 2 to see what phyla are driving this change
HC_p_bf_y2_indic<-signassoc(HC_16S_P_map_bf_y2_com, cluster=BF_Y2_RT,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_p_bf_y2_indic[,c(ncol(HC_p_bf_y2_indic),(ncol(HC_p_bf_y2_indic)-1))]
#indicator analysis found 5 phyla: Acidobacteria (12 (Salmon 311)), Armatimonadetes (7 (Salmon 0)), Cyanobacteria (2 (Control 16)), Firmicutes (3 (Control 162)), Proteobacteria (9 (Salmon 162))
HC_p_bf_y2_indic
#Two indicate salmon reaches after introduction: Proteobacteria and Acidobacteria indicate salmon reaches

#Of the six, only plot most abundant: Proteobacteria, Cyanobacteria, and Acidobacteria
HC_16S_P_map_bf<-subset(HC_16S_P_map, Source=="Biofilm")
HC_16S_P_map_bf$k__Bacteria.p__Proteobacteria<-as.numeric(HC_16S_P_map_bf$k__Bacteria.p__Proteobacteria)
ProBF <- summarySE(HC_16S_P_map_bf, measurevar="k__Bacteria.p__Proteobacteria", groupvars=c("Days_Since_Study_Start","Reach"))
names(ProBF)[names(ProBF) == 'k__Bacteria.p__Proteobacteria'] <- 'Mean'
ProBF$Phylum<-rep("Proteobacteria",24)
HC_16S_P_map_bf$k__Bacteria.p__Cyanobacteria<-as.numeric(HC_16S_P_map_bf$k__Bacteria.p__Cyanobacteria)
CyBF <- summarySE(HC_16S_P_map_bf, measurevar="k__Bacteria.p__Cyanobacteria", groupvars=c("Days_Since_Study_Start","Reach"))
names(CyBF)[names(CyBF) == 'k__Bacteria.p__Cyanobacteria'] <- 'Mean'
CyBF$Phylum<-rep("Cyanobacteria",24)
IndBF<-rbind(ProBF,CyBF)
ggplot(IndBF, aes(x=Days_Since_Study_Start, y=Mean, colour=Phylum)) + 
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach)) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean Abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=two_col_vec_taxa) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")
#########internal insect by year by reach########
#Insects

#Baetis bruniecolor
#Create baetis matrix with metadata
H_C_16S_map_Baetis<-subset(Hunt_Creek_16S_uni_map, Source=="Baetis")
#without metadata for baetis
H_C_16S_Baetis<-H_C_16S_map_Baetis[,18:ncol(H_C_16S_map_Baetis)]
#Baetis environmental variable table
H_C_16S_env_Baetis<-H_C_16S_map_Baetis[,1:17]
Total_Biofilm_Growth_PostCarcass_Bae<-as.factor(H_C_16S_env_Baetis$Total_Biofilm_Growth_PostCarcass)
Reach_Baetis<-as.factor(H_C_16S_env_Baetis$Reach)
Source_Baetis<-as.factor(H_C_16S_env_Baetis$Source)

#Create baetis matrix with metadata
H_C_16S_OTU_map_Baetis<-subset(HC_16S_OTU_map, Source=="Baetis")
#without metadata for baetis
H_C_16S_OTU_Baetis<-H_C_16S_OTU_map_Baetis[,18:ncol(H_C_16S_OTU_map_Baetis)]
#Baetis environmental variable table
H_C_16S_OTU_env_Baetis<-H_C_16S_OTU_map_Baetis[,1:17]
Tot_PostCarcass_OTU_Bae<-as.factor(H_C_16S_OTU_env_Baetis$Total_Biofilm_Growth_PostCarcass)
Reach_OTU_Baetis<-as.factor(H_C_16S_OTU_env_Baetis$Reach)
Source_OTU_Baetis<-as.factor(H_C_16S_OTU_env_Baetis$Source)

#split up by year

#Create baetis matrix for year 1 with metadata
H_C_16S_map_Y1_Bae<-subset(H_C_16S_OTU_map_Baetis, Year=="1")
#how many OTUs Baetis Y1?
#Delete OTUs with no observations
cols_to_drop_Bae_Y1 = c(rep(TRUE, 17), colSums(H_C_16S_map_Y1_Bae[,18:ncol(H_C_16S_map_Y1_Bae)]) > 0)
H_C_16S_map_Y1_Bae<-H_C_16S_map_Y1_Bae[,cols_to_drop_Bae_Y1]
#There are 1898 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_map_Y1_Bae[,18:ncol(H_C_16S_map_Y1_Bae)]))
#40,800 total sequence reads

#Create baetis matrix for year 2 with metadata
H_C_16S_map_Y2_Bae<-subset(H_C_16S_OTU_map_Baetis, Year=="2")
#how many OTUs Baetis Y2?
#Delete OTUs with no observations
cols_to_drop_Bae_Y2 = c(rep(TRUE, 17), colSums(H_C_16S_map_Y2_Bae[,18:ncol(H_C_16S_map_Y2_Bae)]) > 0)
H_C_16S_map_Y2_Bae<-H_C_16S_map_Y2_Bae[,cols_to_drop_Bae_Y2]
#There are 2269 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_map_Y2_Bae[,18:ncol(H_C_16S_map_Y2_Bae)]))
#72,000 total sequence reads

#Year 1
H_C_16S_map_Bae_y1<-subset(H_C_16S_map_Baetis, Year=="1")
#without metadata for baetis
H_C_16S_Bae_y1<-H_C_16S_map_Bae_y1[,18:ncol(H_C_16S_map_Bae_y1)]
H_C_16S_Bae_y1_samples<-as.vector(rownames(H_C_16S_Bae_y1))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_Bae_y1<-as.matrix(H_C_16S_Bae_y1)
H_C_16S_Bae_y1<-subset(H_C_16S_Bae_y1, select=c(H_C_16S_Bae_y1_samples))
#Baetis year 1 environmental variable table
H_C_16S_env_Bae_y1<-H_C_16S_map_Bae_y1[,1:17]
Total_Biofilm_Growth_PostCarcass_bae_y1<-as.factor(H_C_16S_env_Bae_y1$Total_Biofilm_Growth_PostCarcass)
Reach_Bae_y1<-as.factor(H_C_16S_env_Bae_y1$Reach)
Source_Bae_y1<-as.factor(H_C_16S_env_Bae_y1$Source)
#Year 2
H_C_16S_map_Bae_y2<-subset(H_C_16S_map_Baetis, Year=="2")
#without metadata for baetis
H_C_16S_Bae_y2<-H_C_16S_map_Bae_y2[,18:ncol(H_C_16S_map_Bae_y2)]
H_C_16S_Bae_y2_samples<-as.vector(rownames(H_C_16S_Bae_y2))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_Bae_y2<-as.matrix(H_C_16S_Bae_y2)
H_C_16S_Bae_y2<-subset(H_C_16S_Bae_y2, select=c(H_C_16S_Bae_y2_samples))
#Baetis year 2 environmental variable table
H_C_16S_env_Bae_y2<-H_C_16S_map_Bae_y2[,1:17]
Total_Biofilm_Growth_PostCarcass_bae_y2<-as.factor(H_C_16S_env_Bae_y2$Total_Biofilm_Growth_PostCarcass)
Reach_Bae_y2<-as.factor(H_C_16S_env_Bae_y2$Reach)
Source_Bae_y2<-as.factor(H_C_16S_env_Bae_y2$Source)

#Baetis permanova using unifrac
adonis(as.dist(H_C_16S_Bae_y1) ~ Reach*Total_Biofilm_Growth_PostCarcass_bae_y1, data=H_C_16S_env_Bae_y1,permutations=999)
adonis(as.dist(H_C_16S_Bae_y2) ~ Reach*Total_Biofilm_Growth_PostCarcass_bae_y2, data=H_C_16S_env_Bae_y2,permutations=999)

HC_Baetis_y1_NMDS<-metaMDS(as.dist(H_C_16S_Bae_y1))
ordiplot(HC_Baetis_y1_NMDS, type="n", main="Baetis brunicolor 16S year 1")
with(HC_Baetis_y1_NMDS, points(HC_Baetis_y1_NMDS, display="sites", col=two_col_vec_reach[Reach_Bae_y1], pch=19, pt.bg=two_col_vec_reach))
with(HC_Baetis_y1_NMDS, legend("topleft", legend=levels(Reach_Bae_y1), bty="n", col=two_col_vec_reach, pch=19, pt.bg=two_col_vec_reach))
with(HC_Baetis_y1_NMDS, ordiellipse(HC_Baetis_y1_NMDS, Reach_Bae_y1, kind="se", conf=0.95, lwd=2, col="skyblue3", show.groups = "Control"))
with(HC_Baetis_y1_NMDS, ordiellipse(HC_Baetis_y1_NMDS, Reach_Bae_y1, kind="se", conf=0.95, lwd=2, col="tomato3", show.groups = "Salmon"))

HC_Baetis_y2_NMDS<-metaMDS(as.dist(H_C_16S_Bae_y2))
ordiplot(HC_Baetis_y2_NMDS, type="n", main="Baetis brunicolor 16S year 2")
with(HC_Baetis_y2_NMDS, points(HC_Baetis_y2_NMDS, display="sites", col=two_col_vec_reach[Reach_Bae_y2], pch=19, pt.bg=two_col_vec_reach))
with(HC_Baetis_y2_NMDS, legend("topleft", legend=levels(Reach_Bae_y2), bty="n", col=two_col_vec_reach, pch=19, pt.bg=two_col_vec_reach))
with(HC_Baetis_y2_NMDS, ordiellipse(HC_Baetis_y2_NMDS, Reach_Bae_y2, kind="se", conf=0.95, lwd=2, col="skyblue3", show.groups = "Control"))
with(HC_Baetis_y2_NMDS, ordiellipse(HC_Baetis_y2_NMDS, Reach_Bae_y2, kind="se", conf=0.95, lwd=2, col="tomato3", show.groups = "Salmon"))

#indicator taxa analysis for baetis
#subset phyla info
HC_16S_P_map_bae<-subset(HC_16S_P_map, Source=="Baetis")
HC_16S_P_map_bae[,1:17]<-sapply(HC_16S_P_map_bae[,1:17], as.factor)

#Year 1
HC_16S_P_map_bae_y1<-subset(HC_16S_P_map_bae, Year=="1")
HC_16S_P_map_bae_y1_com<-HC_16S_P_map_bae_y1[,18:ncol(HC_16S_P_map_bae_y1)]
#find most abundant phyla
sort(colSums(HC_16S_P_map_bae_y1_com),decreasing=TRUE)
#Proteobacteria most abundant
0.9*sum(colSums(HC_16S_P_map_bae_y1_com))
#90% of sequence reads are 36721, so top 7 pyla represent 90% of sample (Proteobactera, Tenericutes, Verrucomicrobia, Firmicutes, Cyanobacteria, Actinobacteria and Bacteroidetes)
HC_16S_P_map_bae_y1_env<-HC_16S_P_map_bae_y1[1:17]

#Year 2
HC_16S_P_map_bae_y2<-subset(HC_16S_P_map_bae, Year=="2")
HC_16S_P_map_bae_y2_com<-HC_16S_P_map_bae_y2[,18:ncol(HC_16S_P_map_bae_y2)]
#find most abundant phyla
sort(colSums(HC_16S_P_map_bae_y2_com),decreasing=TRUE)
#Proteobacteria most abundant
0.9*sum(colSums(HC_16S_P_map_bae_y2_com))
#90% of sequence reads are 64800, so top 6 pyla represent 90% of sample (Proteobactera, Cyanobacteria, Tenericutes, Bacteroidetes, Actinobacteria and Verrucomicrobia)
HC_16S_P_map_bae_y2_env<-HC_16S_P_map_bae_y2[1:17]

#create clusters based on time and reach
Bae_Y1_Reach<-HC_16S_P_map_bae_y1_env$Reach
Bae_Y2_Reach<-HC_16S_P_map_bae_y2_env$Reach
Bae_Y1_Time<-HC_16S_P_map_bae_y1_env$Total_Biofilm_Growth_PostCarcass
Bae_Y2_Time<-HC_16S_P_map_bae_y2_env$Total_Biofilm_Growth_PostCarcass
Bae_Y1_RT<-paste(Bae_Y1_Reach, Bae_Y1_Time)
Bae_Y2_RT<-paste(Bae_Y2_Reach, Bae_Y2_Time)
#Indicator analysis of Baetids year 1 to see what phyla are driving this change
HC_p_bae_y1_indic<-signassoc(HC_16S_P_map_bae_y1_com, cluster=Bae_Y1_RT,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_p_bae_y1_indic[,c(ncol(HC_p_bae_y1_indic),(ncol(HC_p_bae_y1_indic)-1))]
#indicator analysis found ***Acidobacteria*** indicating Salmon 250 days
#Indicator analysis of Baetids year 2 to see what phyla are driving this change
HC_p_bae_y2_indic<-signassoc(HC_16S_P_map_bae_y2_com, cluster=Bae_Y2_RT,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_p_bae_y2_indic[,c(ncol(HC_p_bae_y2_indic),(ncol(HC_p_bae_y2_indic)-1))]
#indicator analysis found three taxa: Chloroflexi (1), Cyanobacteria(3), Tenericutes (1) indicating 
#***Cyanobacteria*** represents Salmon 162 days

#Heptagenia
#Create heptagenia matrix with metadata
H_C_16S_map_Hepta<-subset(Hunt_Creek_16S_uni_map, Source=="Heptagenia")
#without metadata for baetis
H_C_16S_Hepta<-H_C_16S_map_Hepta[,18:ncol(H_C_16S_map_Hepta)]
H_C_16S_Hep_samples<-as.vector(rownames(H_C_16S_Hepta))
H_C_16S_Hepta<-as.matrix(H_C_16S_Hepta)
H_C_16S_Hepta<-subset(H_C_16S_Hepta, select=c(H_C_16S_Hep_samples))
#Heptagenia environmental variable table
H_C_16S_env_Hepta<-H_C_16S_map_Hepta[,1:17]
Total_Biofilm_Growth_PostCarcass_Hep<-as.factor(H_C_16S_env_Hepta$Total_Biofilm_Growth_PostCarcass)
Reach_Hep<-as.factor(H_C_16S_env_Hepta$Reach)

#Create heptagenia matrix with metadata
H_C_16S_OTU_map_Hepta<-subset(HC_16S_OTU_map, Source=="Heptagenia")
#without metadata for heptagenia
H_C_16S_OTU_Hepta<-H_C_16S_OTU_map_Hepta[,18:ncol(H_C_16S_OTU_map_Hepta)]
#Heptagenia environmental variable table
H_C_16S_OTU_env_Hepta<-H_C_16S_OTU_map_Hepta[,1:17]
Tot_PostCarcass_OTU_Hep<-as.factor(H_C_16S_OTU_env_Hepta$Total_Biofilm_Growth_PostCarcass)
Reach_OTU_Hepta<-as.factor(H_C_16S_OTU_env_Hepta$Reach)
Source_OTU_Hepta<-as.factor(H_C_16S_OTU_env_Hepta$Source)

#Don't need to split up by year, because only year 2

#how many OTUs Heptagenia?
#Delete OTUs with no observations
cols_to_drop_Hep = c(rep(TRUE, 17), colSums(H_C_16S_OTU_map_Hepta[,18:ncol(H_C_16S_OTU_map_Hepta)]) > 0)
H_C_16S_OTU_map_Hepta<-H_C_16S_OTU_map_Hepta[,cols_to_drop_Hep]
#There are 979 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_OTU_map_Hepta[,18:ncol(H_C_16S_OTU_map_Hepta)]))
#16,800 total sequence reads

#Heptagenia permanova using unifrac distance
adonis(as.dist(H_C_16S_Hepta) ~ Tot_PostCarcass_OTU_Hep, permutations=999)
#p=0.05

#indicator taxa analysis for heptagenia
#subset phyla info
HC_16S_P_map_hep<-subset(HC_16S_P_map, Source=="Heptagenia")
HC_16S_P_map_hep[,1:17]<-sapply(HC_16S_P_map_hep[,1:17], as.factor)
HC_16S_P_map_hep_com<-HC_16S_P_map_hep[,18:ncol(HC_16S_P_map_hep)]
#find most abundant phyla
sort(colSums(HC_16S_P_map_hep_com),decreasing=TRUE)
#Proteobacteria most abundant
HC_16S_P_map_hep_env<-HC_16S_P_map_hep[1:17]

HC_p_hep_indic<-signassoc(HC_16S_P_map_hep_com, cluster=HC_16S_P_map_hep_env$Total_Biofilm_Growth_PostCarcass,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_p_hep_indic[,c(ncol(HC_p_hep_indic),(ncol(HC_p_hep_indic)-1))]
#No significant phyla

#Simuliidae
#Create simulidae unifrac matrix with metadata
H_C_16S_map_Simuliidae<-subset(Hunt_Creek_16S_uni_map, Source=="Simulium" | Source=="Stegopterna")
#without metadata for simuliidae
H_C_16S_Sim<-H_C_16S_map_Simuliidae[,18:ncol(H_C_16S_map_Simuliidae)]
H_C_16S_Sim_samples<-as.vector(rownames(H_C_16S_Sim))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_Sim<-as.matrix(H_C_16S_Sim)
H_C_16S_Sim<-subset(H_C_16S_Sim, select=c(H_C_16S_Sim_samples))
#Simuliidae environmental variable table
H_C_16S_env_Sim<-H_C_16S_map_Simuliidae[,1:17]
Total_Biofilm_Growth_PostCarcass_Sim<-as.factor(H_C_16S_env_Sim$Total_Biofilm_Growth_PostCarcass)
Reach_Sim<-as.factor(H_C_16S_env_Sim$Reach)
Source_Sim<-as.factor(H_C_16S_env_Sim$Source)

#Create simulidae count matrix with metadata
H_C_16S_OTU_map_Simuliidae<-subset(HC_16S_OTU_map, Source=="Simulium" | Source=="Stegopterna")
#without metadata for simuliidae
H_C_16S_OTU_Sim<-H_C_16S_OTU_map_Simuliidae[,18:ncol(H_C_16S_OTU_map_Simuliidae)]
#Simuliidae environmental variable table
H_C_16S_OTU_env_Sim<-H_C_16S_OTU_map_Simuliidae[,1:17]
Total_Biofilm_Growth_PostCarcass_OTU_Sim<-as.factor(H_C_16S_OTU_env_Sim$Total_Biofilm_Growth_PostCarcass)
Reach_OTU_Sim<-as.factor(H_C_16S_OTU_env_Sim$Reach)
Source_OTU_Sim<-as.factor(H_C_16S_OTU_env_Sim$Source)

#split up by year

#Create simuliidae count matrix for year 1 with metadata
H_C_16S_OTU_map_Y1_Sim<-subset(H_C_16S_OTU_map_Simuliidae, Year=="1")
#how many OTUs Simuliidae Y1?
#Delete OTUs with no observations
cols_to_drop_OTU_Sim_Y1 = c(rep(TRUE, 17), colSums(H_C_16S_OTU_map_Y1_Sim[,18:ncol(H_C_16S_OTU_map_Y1_Sim)]) > 0)
H_C_16S_OTU_map_Y1_Sim<-H_C_16S_OTU_map_Y1_Sim[,cols_to_drop_OTU_Sim_Y1]
#There are 824 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_OTU_map_Y1_Sim[,18:ncol(H_C_16S_OTU_map_Y1_Sim)]))
#12,000 total sequence reads

#Create simuliidae count matrix for year 1 with metadata
H_C_16S_OTU_map_Y2_Sim<-subset(H_C_16S_OTU_map_Simuliidae, Year=="2")
#how many OTUs Simuliidae Y2?
#Delete OTUs with no observations
cols_to_drop_OTU_Sim_Y2 = c(rep(TRUE, 17), colSums(H_C_16S_OTU_map_Y2_Sim[,18:ncol(H_C_16S_OTU_map_Y2_Sim)]) > 0)
H_C_16S_OTU_map_Y2_Sim<-H_C_16S_OTU_map_Y2_Sim[,cols_to_drop_OTU_Sim_Y2]
#There are 2031 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_OTU_map_Y2_Sim[,18:ncol(H_C_16S_OTU_map_Y2_Sim)]))
#43,200 total sequence reads

#Year 1
H_C_16S_map_Sim_y1<-subset(H_C_16S_map_Simuliidae, Year=="1")
#without metadata for baetis
H_C_16S_Sim_y1<-H_C_16S_map_Sim_y1[,18:ncol(H_C_16S_map_Sim_y1)]
H_C_16S_Sim_y1_samples<-as.vector(rownames(H_C_16S_Sim_y1))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_Sim_y1<-as.matrix(H_C_16S_Sim_y1)
H_C_16S_Sim_y1<-subset(H_C_16S_Sim_y1, select=c(H_C_16S_Sim_y1_samples))
#Simuliidae year 1 environmental variable table
H_C_16S_env_Sim_y1<-H_C_16S_map_Sim_y1[,1:17]
Total_Biofilm_Growth_PostCarcass_sim_y1<-as.factor(H_C_16S_env_Sim_y1$Total_Biofilm_Growth_PostCarcass)
Reach_Sim_y1<-as.factor(H_C_16S_env_Sim_y1$Reach)
Source_Sim_y1<-as.factor(H_C_16S_env_Sim_y1$Source)
#Year 2
H_C_16S_map_Sim_y2<-subset(H_C_16S_map_Simuliidae, Year=="2")
#without metadata for baetis
H_C_16S_Sim_y2<-H_C_16S_map_Sim_y2[,18:ncol(H_C_16S_map_Sim_y2)]
H_C_16S_Sim_y2_samples<-as.vector(rownames(H_C_16S_Sim_y2))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_Sim_y2<-as.matrix(H_C_16S_Sim_y2)
H_C_16S_Sim_y2<-subset(H_C_16S_Sim_y2, select=c(H_C_16S_Sim_y2_samples))
#Simuliidae year 1 environmental variable table
H_C_16S_env_Sim_y2<-H_C_16S_map_Sim_y2[,1:17]
Total_Biofilm_Growth_PostCarcass_sim_y2<-as.factor(H_C_16S_env_Sim_y2$Total_Biofilm_Growth_PostCarcass)
Reach_Sim_y2<-as.factor(H_C_16S_env_Sim_y2$Reach)
Source_Sim_y2<-as.factor(H_C_16S_env_Sim_y2$Source)

#Simuliidae permanova using unifrac
adonis(as.dist(H_C_16S_Sim) ~ Reach*Total_Biofilm_Growth_PostCarcass, data=H_C_16S_env_Sim,permutations=999)
adonis(as.dist(H_C_16S_Sim_y1) ~ Reach*Total_Biofilm_Growth_PostCarcass, data=H_C_16S_env_Sim_y1,permutations=999)
adonis(as.dist(H_C_16S_Sim_y2) ~ Reach*Total_Biofilm_Growth_PostCarcass, data=H_C_16S_env_Sim_y2,permutations=999)

#indicator taxa analysis for Simuliidae
#subset phyla info
HC_16S_P_map_sim<-subset(HC_16S_P_map, Source=="Simulium" | Source=="Stegopterna")
HC_16S_P_map_sim[,1:17]<-sapply(HC_16S_P_map_sim[,1:17], as.factor)

#Year 1
HC_16S_P_map_sim_y1<-subset(HC_16S_P_map_sim, Year=="1")
HC_16S_P_map_sim_y1_com<-HC_16S_P_map_sim_y1[,18:ncol(HC_16S_P_map_sim_y1)]
#find most abundant phyla
sort(colSums(HC_16S_P_map_sim_y1_com),decreasing=TRUE)
#Firmicutes most abundant
HC_16S_P_map_sim_y1_env<-HC_16S_P_map_sim_y1[1:17]

#Year 2
HC_16S_P_map_sim_y2<-subset(HC_16S_P_map_sim, Year=="2")
HC_16S_P_map_sim_y2_com<-HC_16S_P_map_sim_y2[,18:ncol(HC_16S_P_map_sim_y2)]
#find most abundant phyla
sort(colSums(HC_16S_P_map_sim_y2_com),decreasing=TRUE)
#Proteobacteria most abundant
HC_16S_P_map_sim_y2_env<-HC_16S_P_map_sim_y2[1:17]

#create clusters based on time and reach
Sim_Y1_Reach<-HC_16S_P_map_sim_y1_env$Reach
Sim_Y2_Reach<-HC_16S_P_map_sim_y2_env$Reach
Sim_Y1_Time<-HC_16S_P_map_sim_y1_env$Total_Biofilm_Growth_PostCarcass
Sim_Y2_Time<-HC_16S_P_map_sim_y2_env$Total_Biofilm_Growth_PostCarcass
Sim_Y1_RT<-paste(Sim_Y1_Reach, Sim_Y1_Time)
Sim_Y2_RT<-paste(Sim_Y2_Reach, Sim_Y2_Time)

#Indicator analysis of Simuliids year 1 to see what phyla are driving this change
HC_p_sim_y1_indic<-signassoc(HC_16S_P_map_sim_y1_com, cluster=Sim_Y1_RT,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_p_sim_y1_indic[,c(ncol(HC_p_sim_y1_indic),(ncol(HC_p_sim_y1_indic)-1))]
#indicator analysis found no phyla
#Indicator analysis of Simuliids year 2 to see what phyla are driving this change
HC_p_sim_y2_indic<-signassoc(HC_16S_P_map_sim_y2_com, cluster=Sim_Y2_RT,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_p_sim_y2_indic[,c(ncol(HC_p_sim_y2_indic),(ncol(HC_p_sim_y2_indic)-1))]
#indicator analysis found Bacteroidetes (4) and Proteobacteria (2)
#Only ***Bacteroidetes*** represents Salmon reach (day 221)

#Create figure with the three indicator taxa faceted with population graphs
HC_16S_P_map_bae$k__Bacteria.p__Acidobacteria<-as.numeric(HC_16S_P_map_bae$"k__Bacteria.p__Acidobacteria")
AciBae <- summarySE(HC_16S_P_map_bae, measurevar="k__Bacteria.p__Acidobacteria", groupvars=c("Days_Since_Study_Start","Reach"))
AciBae[is.na(AciBae)] <- 0
AciBae$Days_Since_Study_Start<-as.numeric(AciBae$Days_Since_Study_Start)
AciBae$Taxa<-rep("Acidobacteria", 15)
names(AciBae)[names(AciBae) == 'k__Bacteria.p__Acidobacteria'] <- 'Mean'
HC_16S_P_map_bae$k__Bacteria.p__Cyanobacteria<-as.numeric(HC_16S_P_map_bae$"k__Bacteria.p__Cyanobacteria")
CyBae <- summarySE(HC_16S_P_map_bae, measurevar="k__Bacteria.p__Cyanobacteria", groupvars=c("Days_Since_Study_Start","Reach"))
CyBae[is.na(CyBae)] <- 0
CyBae$Days_Since_Study_Start<-as.numeric(CyBae$Days_Since_Study_Start)
CyBae$Taxa<-rep("Cyanobacteria", 15)
names(CyBae)[names(CyBae) == 'k__Bacteria.p__Cyanobacteria'] <- 'Mean'
IndBae<-rbind(AciBae,CyBae)
IndBae$Source<-rep("Baetis brunneicolor", 30)
HC_16S_P_map_sim$k__Bacteria.p__Bacteroidetes<-as.numeric(HC_16S_P_map_sim$"k__Bacteria.p__Bacteroidetes")
BacSim <- summarySE(HC_16S_P_map_sim, measurevar="k__Bacteria.p__Bacteroidetes", groupvars=c("Days_Since_Study_Start","Reach"))
BacSim[is.na(BacSim)] <- 0
BacSim$Days_Since_Study_Start<-as.numeric(BacSim$Days_Since_Study_Start)
BacSim$Taxa<-rep("Bacteroidetes", 8)
names(BacSim)[names(BacSim) == 'k__Bacteria.p__Bacteroidetes'] <- 'Mean'
BacSim$Source<-rep("Stegopterna mutata", 8)
IndInv<-rbind(IndBae,BacSim)
IndInv$Type<-rep("Internal Indicator Microbe", 38)
#Add macroinvertebrate population data
#Already run through macroinvertebrate code, so dataframes are uploaded
names(HC_M_Matrix)[names(HC_M_Matrix) == "Hexapoda Ephemeroptera Baetidae"] <- 'Baetidae_count'
btdfac<-summarySE(HC_M_Matrix, measurevar="Baetidae_count", groupvars=c("Days_Since_Study_Start", "Reach"))
btdfac$Taxa<-rep("Baetis brunneicolor", 26)
btdfac$Source<-rep("Baetis brunneicolor", 26)
btdfac$Type<-rep("Macroinvertebrate", 26)
names(btdfac)[names(btdfac) == "Baetidae_count"] <- 'Mean'
names(HC_M_Matrix)[names(HC_M_Matrix) == "Hexapoda Diptera Simuliidae"] <- 'Simuliidae_count'
simfac<-summarySE(HC_M_Matrix, measurevar="Simuliidae_count", groupvars=c("Days_Since_Study_Start", "Reach"))
simfac$Taxa<-rep("Stegopterna mutata", 26)
simfac$Source<-rep("Stegopterna mutata", 26)
simfac$Type<-rep("Macroinvertebrate", 26)
names(simfac)[names(simfac) == "Simuliidae_count"] <- 'Mean'
IntPlot<-rbind(IndInv,btdfac,simfac)
IntPlot$Type <- factor(IntPlot$Type, c("Macroinvertebrate", "Internal Indicator Microbe"))
IntPlot$Taxa <- factor(IntPlot$Taxa, c("Baetis brunneicolor", "Stegopterna mutata", "Cyanobacteria", "Acidobacteria", "Bacteroidetes"))
ggplot(IntPlot, aes(x=Days_Since_Study_Start, y=Mean, colour=Taxa)) + 
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach)) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean Abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_colour_brewer(palette="Dark2") +
  scale_x_continuous(breaks=seq(0,700,200)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment") +
  facet_grid(Type ~ Source, scales="free")

#Figure out what taxa are introduced by the carrion and if they persist in biofilms/internal invertebrates
HC_16S_OTU_map_ag<-HC_16S_OTU_map
HC_16S_OTU_map_ag$DateEnv<-as.factor(paste(HC_16S_OTU_map_ag$Date,HC_16S_OTU_map_ag$Reach,HC_16S_OTU_map_ag$Source))
HC_16S_OTU_map_ag<-aggregate(HC_16S_OTU_map_ag[18:(ncol(HC_16S_OTU_map_ag)-1)], by=list(Date=HC_16S_OTU_map_ag$Date, Reach=HC_16S_OTU_map_ag$Reach, Source=HC_16S_OTU_map_ag$Source, Year=HC_16S_OTU_map_ag$Year, DateEnv=HC_16S_OTU_map_ag$DateEnv), FUN=sum)
str(HC_16S_OTU_map_ag)
HC_16S_OTU_map_ag_m<-melt(HC_16S_OTU_map_ag, id=(c("Date","Reach","Source","DateEnv", "Year")))
levels(HC_16S_OTU_map_ag_m$DateEnv)
#find names of OTUs that are found in control biofilms and delete them
HC_16S_OTU_names_cont<-filter(HC_16S_OTU_map_ag_m, (Reach=="Control" & Source=="Biofilm" & value>=1) | (Date=="10/3/14" & value>=1))$variable
#delete OTUs from this list
HC_16S_OTU_map_ag_m_ucar<-HC_16S_OTU_map_ag_m[!(HC_16S_OTU_map_ag_m$variable %in% HC_16S_OTU_names_cont),]
HC_16S_OTU_map_ag_m_ucar<-subset(HC_16S_OTU_map_ag_m_ucar, (DateEnv=="10/4/14 Salmon Carcass" & value>=1))
str(HC_16S_OTU_map_ag_m_ucar)
length(unique(HC_16S_OTU_map_ag_m_ucar$variable))
HC_16S_OTU_ucar_names<-unique(HC_16S_OTU_map_ag_m_ucar$variable)
#4103 OTUs were introduced via carrion the first year
#keep OTUs from this list that are in downstream biofilms
HC_16S_OTU_dstbf<-subset(HC_16S_OTU_map_ag_m, Reach=="Salmon" & Source=="Biofilm" & Date!="10/3/14" & value>=1)
HC_16S_OTU_map_ag_m_ucabf<-HC_16S_OTU_dstbf[(HC_16S_OTU_dstbf$variable %in% HC_16S_OTU_ucar_names),]
length(unique(HC_16S_OTU_map_ag_m_ucabf$variable))
HC_16S_OTU_ucardsbf_names<-unique(HC_16S_OTU_map_ag_m_ucabf$variable)
#404 OTUs were introduced via carrion the first year and persist in ds biofilms
#keep OTUs from this list that are in downstream invertebrates
HC_16S_OTU_dsti<-subset(HC_16S_OTU_map_ag_m, Reach=="Salmon" & Source!="Biofilm" & Source!="Carcass" & Date!="9/5/14" & value>=1)
HC_16S_OTU_map_ag_m_ucabfi<-HC_16S_OTU_dsti[(HC_16S_OTU_dsti$variable %in% HC_16S_OTU_ucardsbf_names),]
length(unique(HC_16S_OTU_map_ag_m_ucabfi$variable))
HC_16S_OTU_ucardsbfi_names<-unique(HC_16S_OTU_map_ag_m_ucabfi$variable)
#54 OTUs were introduced via carrion the first year and persist in ds biofilms and are found in ds invertebrates
#now make sure these are never found in control reaches
HC_16S_control_OTU_names<-unique(subset(HC_16S_OTU_map_ag_m, Reach=="Control")$variable)
HC_16S_OTU_map_ag_m_ucabfinc<-HC_16S_OTU_map_ag_m_ucabfi[!HC_16S_OTU_map_ag_m_ucabfi$variable %in% HC_16S_control_OTU_names,]
length(unique(HC_16S_OTU_map_ag_m_ucabfinc$variable))
HC_16S_OTU_ucardsbfinc_names<-unique(HC_16S_OTU_map_ag_m_ucabfinc$variable)
#There are no OTUs that are introduced via carrion in the first year, persist in ds biofilms, are found in ds invertebrates and never found in control reaches
#now see if there are any OTUs that are unique to carrion, ds invertebrates and not found in us,control invertebrates 
HC_16S_OTU_dsti<-subset(HC_16S_OTU_map_ag_m, Reach=="Salmon" & Source!="Biofilm" & Source!="Carcass" & Date!="9/5/14" & value>=1)
HC_16S_OTU_map_ag_m_ucai<-HC_16S_OTU_dsti[(HC_16S_OTU_dsti$variable %in% HC_16S_OTU_ucar_names),]
length(unique(HC_16S_OTU_map_ag_m_ucai$variable))
HC_16S_OTU_ucardsi_names<-unique(HC_16S_OTU_map_ag_m_ucai$variable)
HC_16S_OTU_map_ag_m_ucainc<-HC_16S_OTU_map_ag_m_ucai[!HC_16S_OTU_map_ag_m_ucai$variable %in% HC_16S_control_OTU_names,]
length(unique(HC_16S_OTU_map_ag_m_ucainc$variable))
HC_16S_OTU_ucardsinc_names<-unique(HC_16S_OTU_map_ag_m_ucainc$variable)
#No OTUs are introduced via carcasses and persist in ds invertebrates and are not found in US, control reaches
HC_16S_OTU_ucardsbf_Y1_names<-HC_16S_OTU_ucardsbf_names
#now move on to year 2
#find names of OTUs that are found in year 1 and delete them
HC_16S_OTU_names_y1c<-filter(HC_16S_OTU_map_ag_m, (Year=="1" | (Reach=="Control" & Source=="Biofilm") | Date=="10/4/15") & value>=1)$variable
#delete OTUs from this list
HC_16S_OTU_map_ag_m_ucar2<-HC_16S_OTU_map_ag_m[!(HC_16S_OTU_map_ag_m$variable %in% HC_16S_OTU_names_y1c),]
HC_16S_OTU_map_ag_m_ucar2<-subset(HC_16S_OTU_map_ag_m_ucar2, (DateEnv=="10/9/15 Salmon Carcass" & value>=1))
str(HC_16S_OTU_map_ag_m_ucar2)
length(unique(HC_16S_OTU_map_ag_m_ucar2$variable))
HC_16S_OTU_ucar2_names<-unique(HC_16S_OTU_map_ag_m_ucar2$variable)
#857 unique OTUs were introduced via carrion in the second year
#keep OTUs from this list that are in downstream biofilms year 2
HC_16S_OTU_dstbfy2<-subset(HC_16S_OTU_map_ag_m, Reach=="Salmon" & Source=="Biofilm" & Year=="2" & Date!="10/4/15" & value>=1)
HC_16S_OTU_map_ag_m_ucabfy2<-HC_16S_OTU_dstbfy2[(HC_16S_OTU_dstbfy2$variable %in% HC_16S_OTU_ucar2_names),]
length(unique(HC_16S_OTU_map_ag_m_ucabfy2$variable))
HC_16S_OTU_ucardsbfy2_names<-unique(HC_16S_OTU_map_ag_m_ucabfy2$variable)
#42 OTUs were introduced via carrion the second year and persist in ds biofilms
#keep OTUs from this list that are in downstream invertebrates year 2
HC_16S_OTU_dstiy2<-subset(HC_16S_OTU_map_ag_m, Reach=="Salmon" & Source!="Biofilm" & Source!="Carcass" & Date!="10/4/15" & Year=="2" & value>=1)
HC_16S_OTU_map_ag_m_ucabfiy2<-HC_16S_OTU_dstiy2[(HC_16S_OTU_dstiy2$variable %in% HC_16S_OTU_ucardsbfy2_names),]
length(unique(HC_16S_OTU_map_ag_m_ucabfiy2$variable))
HC_16S_OTU_ucardsbfiy2_names<-unique(HC_16S_OTU_map_ag_m_ucabfiy2$variable)
#4 OTUs were introduced via carrion the second year and persist in ds biofilms and are found in ds invertebrates
#now make sure these are never found in control reaches
HC_16S_OTU_map_ag_m_ucabfincy2<-HC_16S_OTU_map_ag_m_ucabfiy2[!HC_16S_OTU_map_ag_m_ucabfiy2$variable %in% HC_16S_control_OTU_names,]
length(unique(HC_16S_OTU_map_ag_m_ucabfincy2$variable))
HC_16S_OTU_ucardsbfincy2_names<-unique(HC_16S_OTU_map_ag_m_ucabfincy2$variable)
#There are no OTUs that are introduced via carrion in the second year, persist in ds biofilms, are found in ds invertebrates and never found in control reaches

#find out the relative abundance of the unique OTUs
str(HC_16S_OTU_map)
HC_16S_OTU_map_uc1<-HC_16S_OTU_map[,c(1:17, which(names(HC_16S_OTU_map) %in% HC_16S_OTU_ucardsbf_Y1_names))]
HC_16S_OTU_map_uc1$RowSum<-rowSums(HC_16S_OTU_map_uc1[18:ncol(HC_16S_OTU_map_uc1)])
HC_16S_uc1_tot<-subset(HC_16S_OTU_map_uc1, Reach=="Salmon" & (Source=="Carcass" | Source=="Biofilm"), select=c(1:17,422))
UY1<-summarySE(HC_16S_uc1_tot, measurevar="RowSum", groupvars=c("Days_Since_Study_Start", "Source"))
UY1$Year<-rep("One", 14)
HC_16S_OTU_map_uc2<-HC_16S_OTU_map[,c(1:17, which(names(HC_16S_OTU_map) %in% HC_16S_OTU_ucardsbfy2_names))]
HC_16S_OTU_map_uc2$RowSum<-rowSums(HC_16S_OTU_map_uc2[18:ncol(HC_16S_OTU_map_uc2)])
HC_16S_uc2_tot<-subset(HC_16S_OTU_map_uc2, Year=="2" & Reach=="Salmon" & (Source=="Carcass" | Source=="Biofilm"), select=c(1:17,60))
UY2<-summarySE(HC_16S_uc2_tot, measurevar="RowSum", groupvars=c("Days_Since_Study_Start", "Source"))
UY2$Year<-rep("Two", 7)
UY<-rbind(UY1,UY2)
#now visualize these unique OTUs
ggplot(UY, aes(x=Days_Since_Study_Start, y=RowSum, shape=Year, color=Source)) + 
  geom_errorbar(aes(ymin=RowSum-se, ymax=RowSum+se), width=.1) +
  geom_point(size=3) +
  geom_line(size=1.5, data=UY[UY$Source!="Carcass", ])+
  xlab("Days since study start") +
  ylab("Mean Introduced OTU Abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=two_col_vec_taxa) +
  scale_x_continuous(breaks=seq(0,700,100))