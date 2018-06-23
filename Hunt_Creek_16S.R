#Analysis of 16S microbial communities at Hunt Creek 2014-2016

#########################################
#Load packages
library(reshape)
library(ggplot2)
library(vegan)
library(plyr)
library(dplyr)
library(indicspecies)
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
three_col_vec<- c("#a6cee3", "#1f78b4", "#b2df8a")
four_col_vec<-c("#7fc97f", "#beaed4", "#fdc086", "#ffff99")
five_col_vec<- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")
six_col_vec<- c("#e41a1c", "#377eb8", "green", "#984ea3", "#ff7f00", "#ffff33")
six_col_vec_cont<-c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")
seven_col_vec<-c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")

#First work with rarefication plots to use with presentation
HC_Shannon<-read.table("~/Desktop/HC_Shannon.txt", sep="\t", header=T)
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
Hunt_Creek_16S_map <- read.table("~/Desktop/Hunt_Creek_Map_Filtered.txt", header=T)
#Merge metadata onto rarefication file
HC_Sh_map <-merge(Hunt_Creek_16S_map, HC_Sh_c_m_var_c, by="SampleID")
names(HC_Sh_map)[names(HC_Sh_map)=="variable"] <- "Sequences_per_sample"
HC_Sh_map$Sequences_per_sample<-as.numeric(as.character(HC_Sh_map$Sequences_per_sample))
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
  ylab("Rarefaction measure: Shannon +/- SE") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=seven_col_vec)

#####################################
#Upload data tables generated in QIIME and manipulate to make "r friendly"
#####################################

#Get 16S OTU table
Hunt_Creek_16S<-read.table("~/Desktop/HC_table_tabseparated.txt", sep="\t", header = T)

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
Hunt_Creek_16S_uni <- read.table("~/Desktop/weighted_unifrac_otu_table_HC.txt", header=T)
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
HC_16S_OTU_car_env<-HC_16S_OTU_map_car[,1:17]
HC_16S_OTU_car_com<-HC_16S_OTU_map_car[,18:ncol(HC_16S_OTU_map_car)]
#UNI-carcass permanova with unifrac distances
adonis(as.dist(H_C_16S_uni_car_com) ~ Year, data=H_C_16S_uni_car_env, permutations=999)
#NMDS carcass
HC_NMDS_uni_c<-metaMDS(as.dist(H_C_16S_uni_car_com))
ordiplot(HC_NMDS_uni_c, type="n", main="Salmon carcass 16S communities")
with(HC_NMDS_uni_c, points(HC_NMDS_uni_c, display="sites", col=two_col_vec[HC_16S_uni_car_Year], pch=19, pt.bg=two_col_vec))
with(HC_NMDS_uni_c, legend("topleft", legend=levels(HC_16S_uni_car_Year), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))
with(HC_NMDS_uni_c, ordiellipse(HC_NMDS_uni_c, HC_16S_uni_car_Year, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Year 1"))
with(HC_NMDS_uni_c, ordiellipse(HC_NMDS_uni_c, HC_16S_uni_car_Year, kind="se", conf=0.95, lwd=2, col="bisque2", show.groups = "Year 2"))
#Indicator analysis to see what groups are driving this change
HC_car_indic<-multipatt(HC_16S_OTU_car_com, HC_16S_OTU_car_env$Year, control = how(nperm=999))
summary(HC_car_indic)
#186 indicator OTUs detected for carcasses in groups year 1 and year 2

#upload phyla level info and run indicator analysis for carcass/generate box plots
#Upload phyla level files for each run
Hunt_Creek_16S_P<-read.table("~/Desktop/HC_Phyla_16S.txt", sep="\t", header = T)
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
HC_16S_P_map_car[,1:17]<-sapply(HC_16S_P_map_car[,1:17], as.factor)
HC_16S_P_map_car_env<-HC_16S_P_map_car[,1:17]
HC_16S_P_map_car_com<-HC_16S_P_map_car[,18:ncol(HC_16S_P_map_car)]
#Indicator analysis to see what phyla are driving this change
HC_p_car_indic<-multipatt(HC_16S_P_map_car_com, HC_16S_P_map_car_env$Year, control = how(nperm=999))
summary(HC_p_car_indic)
#10 phyla with year 1 significant: Spirochaetes, euryarchaota, nitrospirae, tenericutes, chlorobi, armatimonadetes, crenarchaota, elusimicrobia, TM7,OP11

#make stacked bar graphs of phyla for each carcass year
#simplify and melt dataset
HC_16S_P_map_car_m<-melt(HC_16S_P_map_car)
HC_16S_P_map_car_m<-subset(HC_16S_P_map_car_m, variable=="k__Bacteria.p__Spirochaetes" | variable=="k__Archaea.p__Euryarchaeota" | variable=="k__Bacteria.p__Nitrospirae" | variable=="k__Bacteria.p__Chlorobi" | variable=="k__Bacteria.p__Armatimonadetes" | variable=="k__Bacteria.p__Tenericutes")
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Spirochaetes","Spirochaetes",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Archaea.p__Euryarchaeota","Euryarchaeota",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Nitrospirae","Nitrospirae",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Tenericutes","Tenericutes",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Chlorobi","Chlorobi",HC_16S_P_map_car_m$variable)
HC_16S_P_map_car_m$variable<-gsub("k__Bacteria.p__Armatimonadetes","Armatimonadetes",HC_16S_P_map_car_m$variable)
ggplot(HC_16S_P_map_car_m,aes(x = Year, y = value,fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") +
  xlab("Carcass Year") +
  ylab("Relative abundance") +
  labs(fill = "Indicator phyla") +
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
H_C_16S_uni_Y1_BF<-H_C_16S_uni_map_Y1_BF[,14:ncol(H_C_16S_uni_map_Y1_BF)]
H_C_16S_uni_Y1_BF_samples<-as.vector(rownames(H_C_16S_uni_Y1_BF))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_uni_Y1_BF<-as.matrix(H_C_16S_uni_Y1_BF)
H_C_16S_uni_Y1_BF<-subset(H_C_16S_uni_Y1_BF, select=c(H_C_16S_uni_Y1_BF_samples))
#UNI-Biofilm year 1 environmental variable table
H_C_16S_env_Y1_BF<-H_C_16S_uni_map_Y1_BF[,1:13]
Total_Biofilm_Growth_PostCarcass_Y1<-as.factor(H_C_16S_env_Y1_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y1<-as.factor(H_C_16S_env_Y1_BF$Reach)

#Biofilm year 2
#Create biofilm matrix for year 2 with metadata
H_C_16S_map_Y2_BF<-subset(HC_16S_OTU_map_Y2, Source=="Biofilm")
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
H_C_16S_uni_Y2_BF<-H_C_16S_uni_map_Y2_BF[,14:ncol(H_C_16S_uni_map_Y2_BF)]
H_C_16S_uni_Y2_BF_samples<-as.vector(rownames(H_C_16S_uni_Y2_BF))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_uni_Y2_BF<-as.matrix(H_C_16S_uni_Y2_BF)
H_C_16S_uni_Y2_BF<-subset(H_C_16S_uni_Y2_BF, select=c(H_C_16S_uni_Y2_BF_samples))
#UNI-Biofilm year 2 environmental variable table
H_C_16S_env_Y2_BF<-H_C_16S_uni_map_Y2_BF[,1:13]
Total_Biofilm_Growth_PostCarcass_Y2<-as.factor(H_C_16S_env_Y2_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y2<-as.factor(H_C_16S_env_Y2_BF$Reach)

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
#subset phyla info into biofilms y2
HC_16S_P_map_bf_y2<-subset(HC_16S_P_map, Source=="Biofilm" & Year=="2")
HC_16S_P_map_bf_y2[,1:17]<-sapply(HC_16S_P_map_bf_y2[,1:17], as.factor)
HC_16S_P_map_bf_y2_env<-HC_16S_P_map_bf_y2[,1:17]
HC_16S_P_map_bf_y2_com<-HC_16S_P_map_bf_y2[,18:ncol(HC_16S_P_map_bf_y2)]
#create clusters based on time and reach
BF_Y1_Reach<-HC_16S_P_map_bf_y1_env$Reach
BF_Y2_Reach<-HC_16S_P_map_bf_y2_env$Reach
#Indicator analysis of biofilms year 1 to see what phyla are driving this change
HC_p_bf_y1_indic<-multipatt(HC_16S_P_map_bf_y1_com, BF_Y1_Reach, control = how(nperm=999))
summary(HC_p_bf_y1_indic)
#indicator analysis found no phyla
#Indicator analysis of biofilms year 1 to see what phyla are driving this change
HC_p_bf_y2_indic<-multipatt(HC_16S_P_map_bf_y2_com, BF_Y2_Reach, control = how(nperm=999))
summary(HC_p_bf_y2_indic)
#indicator analysis found no phyla

#Find most common Phyla
P_totals<-rbind(HC_16S_P_t, colSums(HC_16S_P_t))
P_totals<-P_totals[-c(1:172),]
sort(P_totals,decreasing=TRUE)[1:6]
rowSums(P_totals)

#limit phyla data to biofilms
HC_16S_P_map_BF<-subset(HC_16S_P_map, Source == "Biofilm")
#Make line plot for cyanobacteria for salmon vs control
HC_16S_P_map_BF$k__Bacteria.p__Cyanobacteria<-as.numeric(HC_16S_P_map_BF$"k__Bacteria.p__Cyanobacteria")
HC_16S_P_map_BF$Treatment<-as.factor(HC_16S_P_map_BF$Reach)
cpy <- summarySE(HC_16S_P_map_BF, measurevar="k__Bacteria.p__Cyanobacteria", groupvars=c("Days_Since_Study_Start","Treatment"))
ggplot(cpy, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Cyanobacteria, colour=Treatment)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Cyanobacteria-se, ymax=k__Bacteria.p__Cyanobacteria+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Days of biofilm growth") +
  ylab("Cyanobacteria relative abundance +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18,margin=margin(40,0,0,0)),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))
kruskal.test(k__Bacteria.p__Cyanobacteria ~ Reach, data = HC_16S_P_map_BF) 
#Make line plot for firmicutes for salmon vs control
HC_16S_P_map_BF$k__Bacteria.p__Firmicutes<-as.numeric(HC_16S_P_map_BF$"k__Bacteria.p__Firmicutes")
HC_16S_P_map_BF$Treatment<-as.factor(HC_16S_P_map_BF$Reach)
fir <- summarySE(HC_16S_P_map_BF, measurevar="k__Bacteria.p__Firmicutes", groupvars=c("Days_Since_Study_Start","Treatment"))
ggplot(fir, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Firmicutes, colour=Treatment)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Firmicutes-se, ymax=k__Bacteria.p__Firmicutes+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Days of biofilm growth") +
  ylab("Firmicutes relative abundance +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))
kruskal.test(k__Bacteria.p__Firmicutes ~ Reach, data = HC_16S_P_map_BF) 
#not significant, split into year
HC_16S_P_map_bf_y1$Reach<-sapply(HC_16S_P_map_bf_y1$Reach,as.factor)
kruskal.test(k__Bacteria.p__Firmicutes ~ Reach, data = HC_16S_P_map_bf_y1) 
HC_16S_P_map_bf_y2$Reach<-sapply(HC_16S_P_map_bf_y2$Reach,as.factor)
kruskal.test(k__Bacteria.p__Firmicutes ~ Reach, data = HC_16S_P_map_bf_y2) 
#Tenericutes
#Make line plot for Tenericutes for salmon vs control
HC_16S_P_map_BF$k__Bacteria.p__Tenericutes<-as.numeric(HC_16S_P_map_BF$"k__Bacteria.p__Tenericutes")
HC_16S_P_map_BF$Treatment<-as.factor(HC_16S_P_map_BF$Reach)
ten <- summarySE(HC_16S_P_map_BF, measurevar="k__Bacteria.p__Tenericutes", groupvars=c("Days_Since_Study_Start","Treatment"))
ggplot(ten, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Tenericutes, colour=Treatment)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Tenericutes-se, ymax=k__Bacteria.p__Tenericutes+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Days of biofilm growth") +
  ylab("Tenericutes relative abundance +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))
kruskal.test(k__Bacteria.p__Tenericutes ~ Reach, data = HC_16S_P_map_BF) 
#not sig
#Proteobacteria
#Make line plot for Proteobacteria for salmon vs control
HC_16S_P_map_BF$k__Bacteria.p__Proteobacteria<-as.numeric(HC_16S_P_map_BF$"k__Bacteria.p__Proteobacteria")
HC_16S_P_map_BF$Treatment<-as.factor(HC_16S_P_map_BF$Reach)
pro <- summarySE(HC_16S_P_map_BF, measurevar="k__Bacteria.p__Proteobacteria", groupvars=c("Days_Since_Study_Start","Treatment"))
ggplot(pro, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria, colour=Treatment)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria-se, ymax=k__Bacteria.p__Proteobacteria+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Days of biofilm growth") +
  ylab("Proteobacteria relative abundance +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))
kruskal.test(k__Bacteria.p__Proteobacteria ~ Reach, data = HC_16S_P_map_BF)
#significant
#Make line plot for Spirochaetes for salmon vs control
HC_16S_P_map_BF$k__Bacteria.p__Spirochaetes<-as.numeric(HC_16S_P_map_BF$"k__Bacteria.p__Spirochaetes")
spi <- summarySE(HC_16S_P_map_BF, measurevar="k__Bacteria.p__Spirochaetes", groupvars=c("Days_Since_Study_Start","Treatment"))
ggplot(spi, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Spirochaetes, colour=Treatment)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Spirochaetes-se, ymax=k__Bacteria.p__Spirochaetes+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Days of biofilm growth") +
  ylab("Spirochaetes relative abundance +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))
kruskal.test(k__Bacteria.p__Spirochaetes ~ Reach, data = HC_16S_P_map_BF) 
HC_16S_P_map_bf_y1$Reach<-sapply(HC_16S_P_map_bf_y1$Reach,as.factor)
kruskal.test(k__Bacteria.p__Spirochaetes ~ Reach, data = HC_16S_P_map_bf_y1) 
###########################################################
#Carcass and biofilm comparisons
##########################################################
#Make nmds plot for carcass and biofilms for year 1
HC_16S_map_Y1_CBF<-subset(HC_16S_OTU_map_Y1, Reach=="Salmon")
HC_16S_map_Y1_CBF<-subset(HC_16S_map_Y1_CBF, Date=="10/3/14" | Date=="10/4/14" | Date=="10/18/14")
HC_16S_map_Y1_CBF_com<-HC_16S_map_Y1_CBF[,18:ncol(HC_16S_map_Y1_CBF)]
HC_16S_map_Y1_CBF_env<-HC_16S_map_Y1_CBF[,1:17]
Y1_CBF<-as.factor(HC_16S_map_Y1_CBF_env$Total_Biofilm_Growth)
Y1_CBF<-gsub("0","Carcass",Y1_CBF)
Y1_CBF<-gsub("14","Biofilms Before",Y1_CBF)
Y1_CBF<-gsub("28","Biofilms After",Y1_CBF)
Y1_CBF<-as.factor(Y1_CBF)
levels(Y1_CBF)
HC_CBF_NMDS_Y1<-metaMDS(HC_16S_map_Y1_CBF_com, distance="jaccard")
ordiplot(HC_CBF_NMDS_Y1, type="n", main="Year 1 Biofilms Before and After Carcass")
with(HC_CBF_NMDS_Y1, points(HC_CBF_NMDS_Y1, display="sites", col=three_col_vec[Y1_CBF], pch=19, pt.bg=three_col_vec))
with(HC_CBF_NMDS_Y1, legend("topleft", legend=levels(Y1_CBF), bty="n", col=three_col_vec, pch=19, pt.bg=three_col_vec))
with(HC_CBF_NMDS_Y1, ordiellipse(HC_CBF_NMDS_Y1, Y1_CBF, kind="se", conf=0.95, lwd=2, col="#a6cee3", show.groups = "Biofilms After"))
with(HC_CBF_NMDS_Y1, ordiellipse(HC_CBF_NMDS_Y1, Y1_CBF, kind="se", conf=0.95, lwd=2, col="#1f78b4", show.groups = "Biofilms Before"))
with(HC_CBF_NMDS_Y1, ordiellipse(HC_CBF_NMDS_Y1, Y1_CBF, kind="se", conf=0.95, lwd=2, col="#b2df8a", show.groups = "Carcass"))

#UNI-Make nmds plot for carcass and biofilms for year 1
HC_16S_uni_map_Y1_CBF<-subset(H_C_16S_uni_map_Y1, Reach=="Salmon")
HC_16S_uni_map_Y1_CBF<-subset(HC_16S_uni_map_Y1_CBF, Date=="10/3/14" | Date=="10/4/14" | Date=="10/18/14")
HC_16S_uni_Y1_CBF_samples<-as.vector(rownames(HC_16S_uni_map_Y1_CBF))
HC_16S_uni_Y1_CBF_com<-as.matrix(HC_16S_uni_map_Y1_CBF)
HC_16S_uni_Y1_CBF_com<-subset(HC_16S_uni_Y1_CBF_com, select=c(HC_16S_uni_Y1_CBF_samples))
HC_16S_uni_map_Y1_CBF_env<-HC_16S_uni_map_Y1_CBF[,1:16]
Y1_CBF<-as.factor(HC_16S_uni_map_Y1_CBF_env$Total_Biofilm_Growth)
Y1_CBF<-gsub("0","Carcass",Y1_CBF)
Y1_CBF<-gsub("14","Biofilms Before",Y1_CBF)
Y1_CBF<-gsub("28","Biofilms After",Y1_CBF)
Y1_CBF<-as.factor(Y1_CBF)
levels(Y1_CBF)
HC_CBF_NMDS_Y1<-metaMDS(as.dist(HC_16S_uni_map_Y1_CBF_com))
ordiplot(HC_CBF_NMDS_Y1, type="n", main="Biofilms Before and After and Carcass")
with(HC_CBF_NMDS_Y1, points(HC_CBF_NMDS_Y1, display="sites", col=three_col_vec[Y1_CBF], pch=19, pt.bg=three_col_vec))
with(HC_CBF_NMDS_Y1, legend("topleft", legend=levels(Y1_CBF), bty="n", col=three_col_vec, pch=19, pt.bg=three_col_vec))
HC_16S_uni_Y1_CBF_BF<-subset(HC_16S_uni_map_Y1_CBF, Source=="Biofilm")
HC_16S_uni_Y1_CBF_BF_samples<-as.vector(rownames(HC_16S_uni_Y1_CBF_BF))
HC_16S_uni_Y1_CBF_BF_com<-as.matrix(HC_16S_uni_Y1_CBF_BF)
HC_16S_uni_Y1_CBF_BF_com<-subset(HC_16S_uni_Y1_CBF_BF_com, select=c(HC_16S_uni_Y1_CBF_BF_samples))
HC_16S_uni_Y1_CBF_BF_env<-HC_16S_uni_Y1_CBF_BF[,1:16]
adonis(as.dist(HC_16S_uni_Y1_CBF_BF_com) ~ SourceSink, data=HC_16S_uni_Y1_CBF_BF_env, permutations=999)
#Make nmds plot for carcass and biofilms for year 2
HC_16S_map_Y2_CBF<-subset(HC_16S_OTU_map_Y2, Reach=="Salmon")
HC_16S_map_Y2_CBF<-subset(HC_16S_map_Y2_CBF, Date=="10/4/15" | Date=="10/9/15" | Date=="10/25/15")
HC_16S_map_Y2_CBF_com<-HC_16S_map_Y2_CBF[,18:ncol(HC_16S_map_Y1_CBF)]
HC_16S_map_Y2_CBF_env<-HC_16S_map_Y2_CBF[,1:17]
Y2_CBF<-as.factor(HC_16S_map_Y2_CBF_env$Total_Biofilm_Growth)
Y2_CBF<-gsub("0","Carcass",Y2_CBF)
Y2_CBF<-gsub("14","Biofilms Before",Y2_CBF)
Y2_CBF<-gsub("35","Biofilms After",Y2_CBF)
Y2_CBF<-as.factor(Y2_CBF)
levels(Y2_CBF)
HC_CBF_NMDS_Y2<-metaMDS(HC_16S_map_Y2_CBF_com, distance="jaccard")
ordiplot(HC_CBF_NMDS_Y2, type="n", main="Year 2 Biofilms Before and After Carcass")
with(HC_CBF_NMDS_Y2, points(HC_CBF_NMDS_Y2, display="sites", col=three_col_vec[Y2_CBF], pch=19, pt.bg=three_col_vec))
with(HC_CBF_NMDS_Y2, legend("topleft", legend=levels(Y2_CBF), bty="n", col=three_col_vec, pch=19, pt.bg=three_col_vec))
with(HC_CBF_NMDS_Y2, ordiellipse(HC_CBF_NMDS_Y2, Y2_CBF, kind="se", conf=0.95, lwd=2, col="#a6cee3", show.groups = "Biofilms After"))
with(HC_CBF_NMDS_Y2, ordiellipse(HC_CBF_NMDS_Y2, Y2_CBF, kind="se", conf=0.95, lwd=2, col="#1f78b4", show.groups = "Biofilms Before"))
with(HC_CBF_NMDS_Y2, ordiellipse(HC_CBF_NMDS_Y2, Y2_CBF, kind="se", conf=0.95, lwd=2, col="#b2df8a", show.groups = "Carcass"))

#########internal insect by year by reach########
#Insects

#Create baetis matrix with metadata
H_C_16S_map_Baetis<-subset(HC_16S_OTU_map, Source=="Baetis")
#without metadata for baetis
H_C_16S_Baetis<-H_C_16S_map_Baetis[,18:ncol(H_C_16S_map_Baetis)]
#Baetis year 1 environmental variable table
H_C_16S_env_Baetis<-H_C_16S_map_Baetis[,1:17]
Total_Biofilm_Growth_PostCarcass<-as.factor(H_C_16S_env_Baetis$Total_Biofilm_Growth_PostCarcass)
Reach_Baetis<-as.factor(H_C_16S_env_Baetis$Reach)
Source_Baetis<-as.factor(H_C_16S_env_Baetis$Source)

#Baetis permanova using counts
adonis(H_C_16S_Baetis ~ Reach*Total_Biofilm_Growth_PostCarcass*Year, data=H_C_16S_env_Baetis, method="jaccard", permutations=999)
HC_Baetis_NMDS<-metaMDS(H_C_16S_Baetis, distance="jaccard")
ordiplot(HC_Baetis_NMDS, type="n", main="Baetis brunicolor 16S")
with(HC_Baetis_NMDS, points(HC_Baetis_NMDS, display="sites", col=two_col_vec_reach[Reach_Baetis], pch=19, pt.bg=two_col_vec_reach))
with(HC_Baetis_NMDS, legend("topleft", legend=levels(Reach_Baetis), bty="n", col=two_col_vec_reach, pch=19, pt.bg=two_col_vec_reach))
with(HC_Baetis_NMDS, ordiellipse(HC_Baetis_NMDS, Reach_Baetis, kind="se", conf=0.95, lwd=2, col="skyblue3", show.groups = "Control"))
with(HC_Baetis_NMDS, ordiellipse(HC_Baetis_NMDS, Reach_Baetis, kind="se", conf=0.95, lwd=2, col="tomato3", show.groups = "Salmon"))

HC_baetis_indic<-multipatt(H_C_16S_Baetis, Reach_Baetis, control = how(nperm=999))
summary(HC_baetis_indic)

#indicator taxa analysis for baetis
#subset phyla info into biofilms
HC_16S_P_map_bae<-subset(HC_16S_P_map, Source=="Baetis")
HC_16S_P_map_bae[,1:17]<-sapply(HC_16S_P_map_bae[,1:17], as.factor)
HC_16S_P_map_bae_env<-HC_16S_P_map_bae[,1:17]
HC_16S_P_map_bae_com<-HC_16S_P_map_bae[,18:ncol(HC_16S_P_map_bae)]
#create clusters based on time and reach
Bae_Reach<-HC_16S_P_map_bae_env$Reach
Bae_Time<-HC_16S_P_map_bae_env$Total_Biofilm_Growth
#Indicator analysis of biofilms year 1 to see what phyla are driving this change
HC_p_bae_indic<-multipatt(HC_16S_P_map_bae_com, Bae_Reach, control = how(nperm=999))
summary(HC_p_bae_indic)
#indicator analysis found no phyla
#Find most common Phyla
P_bae_totals<-rbind(HC_16S_P_map_bae_com, colSums(HC_16S_P_map_bae_com))
P_bae_totals<-P_bae_totals[-c(1:47),]
sort(P_bae_totals,decreasing=TRUE)[1:6]
rowSums(P_bae_totals)
#Make line plot for Tenericutes for salmon vs control
HC_16S_P_map_bae$k__Bacteria.p__Tenericutes<-as.numeric(as.character(HC_16S_P_map_bae$"k__Bacteria.p__Tenericutes"))
HC_16S_P_map_bae$Treatment<-as.factor(HC_16S_P_map_bae$Reach)
ten <- summarySE(HC_16S_P_map_bae, measurevar="k__Bacteria.p__Tenericutes", groupvars=c("Days_Since_Study_Start","Treatment"))
ten$Days_Since_Study_Start<-as.numeric(as.character(ten$Days_Since_Study_Start))
ggplot(ten, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Tenericutes, colour=Treatment)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Tenericutes-se, ymax=k__Bacteria.p__Tenericutes+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Days of biofilm growth") +
  ylab("Tenericutes relative abundance +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))
kruskal.test(k__Bacteria.p__Tenericutes ~ Treatment, data = HC_16S_P_map_bae) 

#Create heatmap
HC_P_Heat<-HC_16S_P_map[,c(10,18:ncol(HC_16S_P_map))]
HC_P_Heat <- HC_P_Heat[order(HC_P_Heat$Source),] 
HC_P_Heat$Source<-NULL
HC_P_Heat<-as.matrix(HC_P_Heat)
heatmap(HC_P_Heat, scale="column")

########here
#Figure out what taxa are introduced by the carrion and if they persist

#Make data table for just biofilm before, biofilm after, and carcass; for each year
#Year 1 before after and carcass data table
HC_16S_OTU_map_Y1_CBF<-subset(HC_16S_OTU_map_Y1, Reach=="Salmon")
HC_16S_OTU_map_Y1_CBF<-subset(HC_16S_OTU_map_Y1_CBF, Date=="10/3/14" | Date=="10/4/14" | Date=="10/18/14")
HC_16S_OTU_map_Y1_CBF$Date<-gsub("10/3/14","Before",HC_16S_OTU_map_Y1_CBF$Date)
HC_16S_OTU_map_Y1_CBF$Date<-gsub("10/4/14","Carcass",HC_16S_OTU_map_Y1_CBF$Date)
HC_16S_OTU_map_Y1_CBF$Date<-as.factor(gsub("10/18/14","After",HC_16S_OTU_map_Y1_CBF$Date))
rownames(HC_16S_OTU_map_Y1_CBF)<-HC_16S_OTU_map_Y1_CBF[,1]
HC_16S_OTU_map_Y1_CBF<-HC_16S_OTU_map_Y1_CBF[,-c(1)]

#Year 2 before after and carcass data table
HC_16S_OTU_map_Y2_CBF<-subset(HC_16S_OTU_map_Y2, Reach=="Salmon")
HC_16S_OTU_map_Y2_CBF<-subset(HC_16S_OTU_map_Y2_CBF, Date=="10/4/15" | Date=="10/9/15" | Date=="10/25/15")
HC_16S_OTU_map_Y2_CBF$Date<-gsub("10/4/15","Before",HC_16S_OTU_map_Y2_CBF$Date)
HC_16S_OTU_map_Y2_CBF$Date<-gsub("10/9/15","Carcass",HC_16S_OTU_map_Y2_CBF$Date)
HC_16S_OTU_map_Y2_CBF$Date<-as.factor(gsub("10/25/15","After",HC_16S_OTU_map_Y2_CBF$Date))
rownames(HC_16S_OTU_map_Y2_CBF)<-HC_16S_OTU_map_Y2_CBF[,1]
HC_16S_OTU_map_Y2_CBF<-HC_16S_OTU_map_Y2_CBF[,-c(1)]

#Subset into OTU's that have 0 in before, but >0 in carcass and after

#Year 1
HC_16S_OTU_map_Y1_CBF_sub<-data.frame(t(HC_16S_OTU_map_Y1_CBF[,18:ncol(HC_16S_OTU_map_Y1_CBF)]))
names(HC_16S_OTU_map_Y1_CBF_sub)
HC_16S_OTU_map_Y1_CBF_sub$HC.1<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.1))
HC_16S_OTU_map_Y1_CBF_sub$HC.2<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.2))
HC_16S_OTU_map_Y1_CBF_sub$HC.3<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.3))
HC_16S_OTU_map_Y1_CBF_sub$HC.13<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.13))
HC_16S_OTU_map_Y1_CBF_sub$HC.14<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.14))
HC_16S_OTU_map_Y1_CBF_sub$HC.15<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.15))
HC_16S_OTU_map_Y1_CBF_sub$HC.16<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.16))
HC_16S_OTU_map_Y1_CBF_sub$HC.17<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.17))
HC_16S_OTU_map_Y1_CBF_sub$HC.18<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.18))
HC_16S_OTU_map_Y1_CBF_sub$HC.19<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.19))
HC_16S_OTU_map_Y1_CBF_sub$HC.20<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.20))
HC_16S_OTU_map_Y1_CBF_sub$HC.21<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.21))
HC_16S_OTU_map_Y1_CBF_sub$HC.4<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.4))
HC_16S_OTU_map_Y1_CBF_sub$HC.5<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.5))
HC_16S_OTU_map_Y1_CBF_sub$HC.6<-as.numeric(as.character(HC_16S_OTU_map_Y1_CBF_sub$HC.6))

HC_16S_OTU_map_Y1_CBF_sub$Before.sum<-rowSums(HC_16S_OTU_map_Y1_CBF_sub[,c(1,9,12)])
HC_16S_OTU_map_Y1_CBF_sub$Carcass.sum<-rowSums(HC_16S_OTU_map_Y1_CBF_sub[,c(2,3,4,5,6,7,8,9,11)])
HC_16S_OTU_map_Y1_CBF_sub$After.sum<-rowSums(HC_16S_OTU_map_Y1_CBF_sub[,13:15])

HC_16S_OTU_map_Y1_CBF_sub<-subset(HC_16S_OTU_map_Y1_CBF_sub, Before.sum==0 & Carcass.sum>0 & After.sum>0)
CBF_Unique_Y1<-as.factor(row.names(HC_16S_OTU_map_Y1_CBF_sub))
#The output of row.names is the names of the OTU's that are unique to carcass and biofilms after, not before

#Year 2
HC_16S_OTU_map_Y2_CBF_sub<-data.frame(t(HC_16S_OTU_map_Y2_CBF[,15:ncol(HC_16S_OTU_map_Y2_CBF)]))
names(HC_16S_OTU_map_Y2_CBF_sub)
HC_16S_OTU_map_Y2_CBF_sub<-as.numeric.factor(HC_16S_OTU_map_Y2_CBF_sub)
HC_16S_OTU_map_Y2_CBF_sub$HC.59<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.59))
HC_16S_OTU_map_Y2_CBF_sub$HC.60<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.60))
HC_16S_OTU_map_Y2_CBF_sub$HC.64<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.64))
HC_16S_OTU_map_Y2_CBF_sub$HC.65<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.65))
HC_16S_OTU_map_Y2_CBF_sub$HC.66<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.66))
HC_16S_OTU_map_Y2_CBF_sub$HC.70<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.70))
HC_16S_OTU_map_Y2_CBF_sub$HC.71<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.71))
HC_16S_OTU_map_Y2_CBF_sub$HC.72<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.72))
HC_16S_OTU_map_Y2_CBF_sub$HC.73<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.73))
HC_16S_OTU_map_Y2_CBF_sub$HC.74<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.74))
HC_16S_OTU_map_Y2_CBF_sub$HC.75<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.75))
HC_16S_OTU_map_Y2_CBF_sub$HC.76<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.76))
HC_16S_OTU_map_Y2_CBF_sub$HC.77<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.77))
HC_16S_OTU_map_Y2_CBF_sub$HC.78<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.78))
HC_16S_OTU_map_Y2_CBF_sub$HC.79<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.79))
HC_16S_OTU_map_Y2_CBF_sub$Before.sum<-rowSums(HC_16S_OTU_map_Y2_CBF_sub[,1:3])
HC_16S_OTU_map_Y2_CBF_sub$After.sum<-rowSums(HC_16S_OTU_map_Y2_CBF_sub[,4:6])
HC_16S_OTU_map_Y2_CBF_sub$Carcass.sum<-rowSums(HC_16S_OTU_map_Y2_CBF_sub[,7:16])
HC_16S_OTU_map_Y2_CBF_sub<-subset(HC_16S_OTU_map_Y2_CBF_sub, Before.sum==0 & Carcass.sum>0 & After.sum>0)
row.names(HC_16S_OTU_map_Y2_CBF_sub)
CBF_Unique_Y2<-as.factor(row.names(HC_16S_OTU_map_Y2_CBF_sub))
#The output of row.names is the names of the OTU's that are unique to carcass and biofilms after, not before

#Exclude unique OTU's that are found in upstream biofilms in year 1

#Create data table with only upstream biofilms for year 1 then find unique
HC_16S_OTU_map_Y1_BF_US<-subset(HC_16S_OTU_map_Y1, Reach=="Control" & Source=="Biofilm")
rownames(HC_16S_OTU_map_Y1_BF_US)<-HC_16S_OTU_map_Y1_BF_US[,1]
HC_16S_OTU_map_Y1_BF_US<-HC_16S_OTU_map_Y1_BF_US[,-c(1)]
HC_16S_OTU_map_Y1_BF_US_unique<-data.frame(t(HC_16S_OTU_map_Y1_BF_US[,15:ncol(HC_16S_OTU_map_Y1_BF_US)]))
HC_16S_OTU_map_Y1_BF_US_unique$OTU<-row.names(HC_16S_OTU_map_Y1_BF_US_unique)
HC_16S_OTU_map_Y1_BF_US_unique<-HC_16S_OTU_map_Y1_BF_US_unique[HC_16S_OTU_map_Y1_BF_US_unique$OTU %in% CBF_Unique_Y1,]
HC_16S_OTU_map_Y1_BF_US_unique$HC.10<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.10))
HC_16S_OTU_map_Y1_BF_US_unique$HC.11<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.11))
HC_16S_OTU_map_Y1_BF_US_unique$HC.25<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.25))
HC_16S_OTU_map_Y1_BF_US_unique$HC.26<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.26))
HC_16S_OTU_map_Y1_BF_US_unique$HC.27<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.27))
HC_16S_OTU_map_Y1_BF_US_unique$HC.12<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.12))
HC_16S_OTU_map_Y1_BF_US_unique$HC.31<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.31))
HC_16S_OTU_map_Y1_BF_US_unique$HC.32<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.32))
HC_16S_OTU_map_Y1_BF_US_unique$HC.33<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.33))
HC_16S_OTU_map_Y1_BF_US_unique$HC.37<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.37))
HC_16S_OTU_map_Y1_BF_US_unique$HC.38<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.38))
HC_16S_OTU_map_Y1_BF_US_unique$HC.39<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.39))
HC_16S_OTU_map_Y1_BF_US_unique$HC.7<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.7))
HC_16S_OTU_map_Y1_BF_US_unique$HC.43<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.43))
HC_16S_OTU_map_Y1_BF_US_unique$HC.44<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.44))
HC_16S_OTU_map_Y1_BF_US_unique$HC.45<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.45))
HC_16S_OTU_map_Y1_BF_US_unique$HC.8<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.8))
HC_16S_OTU_map_Y1_BF_US_unique$HC.9<-as.numeric(as.character(HC_16S_OTU_map_Y1_BF_US_unique$HC.9))
HC_16S_OTU_map_Y1_BF_US_unique$Sums<-rowSums(HC_16S_OTU_map_Y1_BF_US_unique[,1:18])
HC_16S_OTU_map_Y1_BF_US_unique<-subset(HC_16S_OTU_map_Y1_BF_US_unique, Sums==0)
row.names(HC_16S_OTU_map_Y1_BF_US_unique)
CBF_really_unique_Y1<-as.factor(row.names(HC_16S_OTU_map_Y1_BF_US_unique))
#output are OTU's that are unique to carcass and biofilms, not found before or upstream for year 1 only

#Repeat with unique from year 2
#Create data table with only upstream biofilms for year 2 then find unique
HC_16S_OTU_map_Y2_BF_US<-subset(HC_16S_OTU_map_Y2, Reach=="Control" & Source=="Biofilm")
rownames(HC_16S_OTU_map_Y2_BF_US)<-HC_16S_OTU_map_Y2_BF_US[,1]
HC_16S_OTU_map_Y2_BF_US<-HC_16S_OTU_map_Y2_BF_US[,-c(1)]
HC_16S_OTU_map_Y2_BF_US_unique<-data.frame(t(HC_16S_OTU_map_Y2_BF_US[,15:ncol(HC_16S_OTU_map_Y2_BF_US)]))
HC_16S_OTU_map_Y2_BF_US_unique$OTU<-row.names(HC_16S_OTU_map_Y2_BF_US_unique)
HC_16S_OTU_map_Y2_BF_US_unique<-HC_16S_OTU_map_Y2_BF_US_unique[HC_16S_OTU_map_Y2_BF_US_unique$OTU %in% CBF_Unique_Y2,]
HC_16S_OTU_map_Y2_BF_US_unique<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique))
HC_16S_OTU_map_Y2_BF_US_unique$HC.62<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.62))
HC_16S_OTU_map_Y2_BF_US_unique$HC.63<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.63))
HC_16S_OTU_map_Y2_BF_US_unique$HC.67<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.67))
HC_16S_OTU_map_Y2_BF_US_unique$HC.68<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.68))
HC_16S_OTU_map_Y2_BF_US_unique$HC.69<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.69))
HC_16S_OTU_map_Y2_BF_US_unique$HC.80<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.80))
HC_16S_OTU_map_Y2_BF_US_unique$HC.81<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.81))
HC_16S_OTU_map_Y2_BF_US_unique$HC.82<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.82))
HC_16S_OTU_map_Y2_BF_US_unique$HC.94<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.94))
HC_16S_OTU_map_Y2_BF_US_unique$HC.95<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.95))
HC_16S_OTU_map_Y2_BF_US_unique$HC.96<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.96))
HC_16S_OTU_map_Y2_BF_US_unique$Sums<-rowSums(HC_16S_OTU_map_Y2_BF_US_unique[,1:12])
HC_16S_OTU_map_Y2_BF_US_unique<-subset(HC_16S_OTU_map_Y2_BF_US_unique, Sums==0)
row.names(HC_16S_OTU_map_Y2_BF_US_unique)
CBF_really_unique_Y2<-row.names(HC_16S_OTU_map_Y2_BF_US_unique)
#output are OTU's that are unique to carcass and biofilms, not found before or upstream for year 2 only

#Find out if any of the 293 really unique OTUs in year 1 persist in year 2 biofilms, both US and DS
HC_16S_OTU_map_Y2_BF<-subset(HC_16S_OTU_map_Y2, Source=="Biofilm")
rownames(HC_16S_OTU_map_Y2_BF)<-HC_16S_OTU_map_Y2_BF[,1]
HC_16S_OTU_map_Y2_BF<-HC_16S_OTU_map_Y2_BF[,-c(1)]
HC_16S_OTU_map_Y2_BF_r_unique<-data.frame(t(HC_16S_OTU_map_Y2_BF[,15:ncol(HC_16S_OTU_map_Y2_BF)]))
HC_16S_OTU_map_Y2_BF_r_unique$OTU<-row.names(HC_16S_OTU_map_Y2_BF_r_unique)
HC_16S_OTU_map_Y2_BF_r_unique<-HC_16S_OTU_map_Y2_BF_r_unique[HC_16S_OTU_map_Y2_BF_r_unique$OTU %in% CBF_really_unique_Y1,]
HC_16S_OTU_map_Y2_BF_r_unique$HC.58<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.58))
HC_16S_OTU_map_Y2_BF_r_unique$HC.59<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.59))
HC_16S_OTU_map_Y2_BF_r_unique$HC.60<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.60))
HC_16S_OTU_map_Y2_BF_r_unique$HC.61<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.61))
HC_16S_OTU_map_Y2_BF_r_unique$HC.62<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.62))
HC_16S_OTU_map_Y2_BF_r_unique$HC.63<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.63))
HC_16S_OTU_map_Y2_BF_r_unique$HC.64<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.64))
HC_16S_OTU_map_Y2_BF_r_unique$HC.65<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.65))
HC_16S_OTU_map_Y2_BF_r_unique$HC.66<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.66))
HC_16S_OTU_map_Y2_BF_r_unique$HC.67<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.67))
HC_16S_OTU_map_Y2_BF_r_unique$HC.68<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.68))
HC_16S_OTU_map_Y2_BF_r_unique$HC.69<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.69))
HC_16S_OTU_map_Y2_BF_r_unique$HC.80<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.80))
HC_16S_OTU_map_Y2_BF_r_unique$HC.81<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.81))
HC_16S_OTU_map_Y2_BF_r_unique$HC.82<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.82))
HC_16S_OTU_map_Y2_BF_r_unique$HC.83<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.83))
HC_16S_OTU_map_Y2_BF_r_unique$HC.84<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.84))
HC_16S_OTU_map_Y2_BF_r_unique$HC.85<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.85))
HC_16S_OTU_map_Y2_BF_r_unique$HC.91<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.91))
HC_16S_OTU_map_Y2_BF_r_unique$HC.92<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.92))
HC_16S_OTU_map_Y2_BF_r_unique$HC.93<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.93))
HC_16S_OTU_map_Y2_BF_r_unique$HC.94<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.94))
HC_16S_OTU_map_Y2_BF_r_unique$HC.95<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.95))
HC_16S_OTU_map_Y2_BF_r_unique$HC.96<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_r_unique$HC.96))
HC_16S_OTU_map_Y2_BF_r_unique$Sums<-rowSums(HC_16S_OTU_map_Y2_BF_r_unique[,1:24])
HC_16S_OTU_map_Y2_BF_r_unique<-subset(HC_16S_OTU_map_Y2_BF_r_unique, Sums>0)
HC_16S_OTU_map_Y2_BF_r_unique_p<-row.names(HC_16S_OTU_map_Y2_BF_r_unique)

#of the "Y1 really unique persisting" OTU's, find out which ones are not found US biofilms in Y2
HC_16S_OTU_map_Y2_BF_US_unique<-data.frame(t(HC_16S_OTU_map_Y2_BF_US[,15:ncol(HC_16S_OTU_map_Y2_BF_US)]))
HC_16S_OTU_map_Y2_BF_US_unique$OTU<-row.names(HC_16S_OTU_map_Y2_BF_US_unique)
HC_16S_OTU_map_Y2_BF_US_r_unique<-HC_16S_OTU_map_Y2_BF_US_unique[HC_16S_OTU_map_Y2_BF_US_unique$OTU %in% HC_16S_OTU_map_Y2_BF_r_unique_p,]
HC_16S_OTU_map_Y2_BF_US_r_unique$X58<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X58))
HC_16S_OTU_map_Y2_BF_US_r_unique$X59<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X59))
HC_16S_OTU_map_Y2_BF_US_r_unique$X60<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X60))
HC_16S_OTU_map_Y2_BF_US_r_unique$X64<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X64))
HC_16S_OTU_map_Y2_BF_US_r_unique$X65<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X65))
HC_16S_OTU_map_Y2_BF_US_r_unique$X66<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X66))
HC_16S_OTU_map_Y2_BF_US_r_unique$X79<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X79))
HC_16S_OTU_map_Y2_BF_US_r_unique$X80<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X80))
HC_16S_OTU_map_Y2_BF_US_r_unique$X81<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X81))
HC_16S_OTU_map_Y2_BF_US_r_unique$X94<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X94))
HC_16S_OTU_map_Y2_BF_US_r_unique$X95<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X95))
HC_16S_OTU_map_Y2_BF_US_r_unique$X96<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_r_unique$X96))
HC_16S_OTU_map_Y2_BF_US_r_unique$Sums<-rowSums(HC_16S_OTU_map_Y2_BF_US_r_unique[,1:12])
HC_16S_OTU_map_Y2_BF_US_unique<-subset(HC_16S_OTU_map_Y2_BF_US_r_unique, Sums>0)
HC_CBF_U_OTUs<-row.names(HC_16S_OTU_map_Y2_BF_US_unique)

#Find which of the 37 unique OTUs are the highest in abundance for both years
HC_U_sub<-Hunt_Creek_16S
HC_U_sub$OTU<-row.names(HC_U_sub)
HC_U_sub_CBF<-HC_U_sub[HC_U_sub$OTU %in% HC_CBF_U_OTUs,]
HC_U_sub_CBF<-as.numeric(HC_U_sub_CBF)
HC_U_sub_CBF$Sums<-rowSums(HC_U_sub_CBF[,1:96])
#Find top sums to utilize in data table for presentations

#Are any of the 37 unique taxa found in internal insect microbiomes?
#Make data table for just insect samples
Hunt_Creek_16S_OTU_map_I<-subset(HC_16S_OTU_map, Source=="Baetis" | Source=="Rhyacophila" | Source=="Stegopterna")
rownames(Hunt_Creek_16S_OTU_map_I)<-Hunt_Creek_16S_OTU_map_I[,1]
Hunt_Creek_16S_OTU_map_I<-Hunt_Creek_16S_OTU_map_I[,-c(1)]
Hunt_Creek_16S_OTU_map_I_U<-data.frame(t(Hunt_Creek_16S_OTU_map_I[,14:ncol(Hunt_Creek_16S_OTU_map_I)]))
Hunt_Creek_16S_OTU_map_I$OTU<-row.names(Hunt_Creek_16S_OTU_map_I)
Hunt_Creek_16S_OTU_map_I_U<-Hunt_Creek_16S_OTU_map_I[Hunt_Creek_16S_OTU_map_I$OTU %in% HC_CBF_U_OTUs,]
Hunt_Creek_16S_OTU_map_I_U$HC.46<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.46))
Hunt_Creek_16S_OTU_map_I_U$HC.47<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.47))
Hunt_Creek_16S_OTU_map_I_U$HC.48<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.48))
Hunt_Creek_16S_OTU_map_I_U$HC.49<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.49))
Hunt_Creek_16S_OTU_map_I_U$HC.50<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.50))
Hunt_Creek_16S_OTU_map_I_U$HC.51<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.51))
Hunt_Creek_16S_OTU_map_I_U$HC.52<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.52))
Hunt_Creek_16S_OTU_map_I_U$HC.53<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.53))
Hunt_Creek_16S_OTU_map_I_U$HC.54<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.54))
Hunt_Creek_16S_OTU_map_I_U$HC.55<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.55))
Hunt_Creek_16S_OTU_map_I_U$HC.56<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.56))
Hunt_Creek_16S_OTU_map_I_U$HC.57<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.57))
Hunt_Creek_16S_OTU_map_I_U$HC.86<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.86))
Hunt_Creek_16S_OTU_map_I_U$HC.87<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.87))
Hunt_Creek_16S_OTU_map_I_U$HC.88<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.88))
Hunt_Creek_16S_OTU_map_I_U$HC.89<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.89))
Hunt_Creek_16S_OTU_map_I_U$HC.90<-as.numeric(as.character(Hunt_Creek_16S_OTU_map_I_U$HC.90))
Hunt_Creek_16S_OTU_map_I_U$Sums<-rowSums(Hunt_Creek_16S_OTU_map_I[,1:17])
Hunt_Creek_16S_OTU_map_I_U<-subset(Hunt_Creek_16S_OTU_map_I, Sums>0)
HC_unique_CBFI_OTUs<-Hunt_Creek_16S_OTU_map_I_U$OTU

HC_ICBF_I_16S_OTUs<-Hunt_Creek_16S_OTU_map_I[Hunt_Creek_16S_OTU_map_I$OTU %in% HC_unique_CBFI_OTUs,]

#Find abundance of certain OTU's in salmon carcasses in year 1
HC_16S_OTU_map_Y1_S<-subset(HC_16S_OTU_map_Y1, Source=="Carcass")
HC_16S_OTU_map_Y1_S$denovo2027<-as.numeric(as.character(HC_16S_OTU_map_Y1_S$denovo2027))
HC_16S_OTU_map_Y1_S$denovo164486<-as.numeric(as.character(HC_16S_OTU_map_Y1_S$denovo164486))
HC_16S_OTU_map_Y1_S$denovo130121<-as.numeric(as.character(HC_16S_OTU_map_Y1_S$denovo130121))
HC_16S_OTU_map_Y1_S$denovo183000<-as.numeric(as.character(HC_16S_OTU_map_Y1_S$denovo183000))
HC_16S_OTU_map_Y1_S$denovo102338<-as.numeric(as.character(HC_16S_OTU_map_Y1_S$denovo102338))
HC_16S_OTU_map_Y1_S$denovo114115<-as.numeric(as.character(HC_16S_OTU_map_Y1_S$denovo114115))
HC_16S_OTU_map_Y1_S$denovo144456<-as.numeric(as.character(HC_16S_OTU_map_Y1_S$denovo144456))
HC_16S_OTU_map_Y1_S$denovo177870<-as.numeric(as.character(HC_16S_OTU_map_Y1_S$denovo177870))
HC_16S_OTU_map_Y1_S$denovo39627<-as.numeric(as.character(HC_16S_OTU_map_Y1_S$denovo39627))

sum(HC_16S_OTU_map_Y1_S$denovo2027)
sum(HC_16S_OTU_map_Y1_S$denovo164486)
sum(HC_16S_OTU_map_Y1_S$denovo130121)
sum(HC_16S_OTU_map_Y1_S$denovo183000)
sum(HC_16S_OTU_map_Y1_S$denovo39627)
sum(HC_16S_OTU_map_Y1_S$denovo102338)
sum(HC_16S_OTU_map_Y1_S$denovo114115)
sum(HC_16S_OTU_map_Y1_S$denovo144456)
sum(HC_16S_OTU_map_Y1_S$denovo177870)


#Find abundance of certain OTU's in biofilms DS for both years
HC_16S_OTU_map_BF_DS<-subset(HC_16S_OTU_map, Source=="Biofilm" & Reach=="Salmon")
HC_16S_OTU_map_BF_DS$denovo2027<-as.numeric(as.character(HC_16S_OTU_map_BF_DS$denovo2027))
HC_16S_OTU_map_BF_DS$denovo164486<-as.numeric(as.character(HC_16S_OTU_map_BF_DS$denovo164486))
HC_16S_OTU_map_BF_DS$denovo130121<-as.numeric(as.character(HC_16S_OTU_map_BF_DS$denovo130121))
HC_16S_OTU_map_BF_DS$denovo183000<-as.numeric(as.character(HC_16S_OTU_map_BF_DS$denovo183000))
HC_16S_OTU_map_BF_DS$denovo39627<-as.numeric(as.character(HC_16S_OTU_map_BF_DS$denovo39627))
HC_16S_OTU_map_BF_DS$denovo102338<-as.numeric(as.character(HC_16S_OTU_map_BF_DS$denovo102338))
HC_16S_OTU_map_BF_DS$denovo114115<-as.numeric(as.character(HC_16S_OTU_map_BF_DS$denovo114115))
HC_16S_OTU_map_BF_DS$denovo144456<-as.numeric(as.character(HC_16S_OTU_map_BF_DS$denovo144456))
HC_16S_OTU_map_BF_DS$denovo177870<-as.numeric(as.character(HC_16S_OTU_map_BF_DS$denovo177870))

sum(HC_16S_OTU_map_BF_DS$denovo2027)
sum(HC_16S_OTU_map_BF_DS$denovo164486)
sum(HC_16S_OTU_map_BF_DS$denovo130121)
sum(HC_16S_OTU_map_BF_DS$denovo183000)
sum(HC_16S_OTU_map_BF_DS$denovo39627)
sum(HC_16S_OTU_map_BF_DS$denovo102338)
sum(HC_16S_OTU_map_BF_DS$denovo114115)
sum(HC_16S_OTU_map_BF_DS$denovo144456)
sum(HC_16S_OTU_map_BF_DS$denovo177870)


#Create melted file (this step takes a while, so prepare to wait. Skip if just doing community analysis like NMDS)
HC_16S_OTU_map_m<-melt(HC_16S_OTU_map, id=c("Row.names"))
HC_16S_OTU_map_m<-data.frame(HC_16S_OTU_map_m)
names(HC_16S_OTU_map_m)

#Delete zeros
HC_16S_OTU_map_m<-subset(HC_16S_OTU_map_m, value > 0)

####################################################
#Now work with individual populations at the phyla level to visualize time by reach interaction for specific taxa
#######################################################

#Now, work with phyla level data for 16S

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

###################
#Sourcetracker plots
###################


#Upload sourcetracker table obtained using QIIME
HC_16S_SourceTracker<-read.table("~/Desktop/sourcetracker_out/sink_predictions.txt", sep="\t", header = T)
row.names(HC_16S_SourceTracker)<-HC_16S_SourceTracker[,1]
HC_16S_SourceTracker$SampleID<-NULL
#Merge metadata
HC_16S_SourceTracker_map <-merge(Hunt_Creek_16S_map, HC_16S_SourceTracker, by=0)


#Pie graphs for insect control

#subset insect control
HC_16S_SourceTracker_map_ic<-subset(HC_16S_SourceTracker_map, Env=="InsectControl")
HC_16S_SourceTracker_map_ic<-HC_16S_SourceTracker_map_ic[,18:ncol(HC_16S_SourceTracker_map_ic)]
Sums_ic<-colSums(HC_16S_SourceTracker_map_ic)
lables_ic<-c("Control Biofilms", "Treatment Biofilms", "Salmon Carcass", "Unknown")
pie(Sums_ic, labels=lables_ic, main="All Control Insects", col=four_col_vec)
#subset by insect taxa in control
#Rhyacophila
HC_16S_SourceTracker_map_icr<-subset(HC_16S_SourceTracker_map, Env=="InsectControl" & Source=="Rhyacophila")
HC_16S_SourceTracker_map_icr<-HC_16S_SourceTracker_map_icr[,18:ncol(HC_16S_SourceTracker_map_icr)]
Sums_icr<-colSums(HC_16S_SourceTracker_map_icr)
lables_icr<-c("Control Biofilms", "Treatment Biofilms", "Salmon Carcass", "Unknown")
pie(Sums_icr, labels=lables_icr, main="Control Rhyacophila", col=four_col_vec)
#Stegopterna
HC_16S_SourceTracker_map_ics<-subset(HC_16S_SourceTracker_map, Env=="InsectControl" & Source=="Stegopterna")
HC_16S_SourceTracker_map_ics<-HC_16S_SourceTracker_map_ics[,18:ncol(HC_16S_SourceTracker_map_ics)]
Sums_ics<-colSums(HC_16S_SourceTracker_map_ics)
lables_ics<-c("Control Biofilms", "Treatment Biofilms", "Salmon Carcass", "Unknown")
pie(Sums_ics, labels=lables_ics, main="Control Stegopterna", col=four_col_vec)
#Baetis
HC_16S_SourceTracker_map_icb<-subset(HC_16S_SourceTracker_map, Env=="InsectControl" & Source=="Baetis")
HC_16S_SourceTracker_map_icb<-HC_16S_SourceTracker_map_icb[,18:ncol(HC_16S_SourceTracker_map_icb)]
Sums_icb<-colSums(HC_16S_SourceTracker_map_icb)
lables_icb<-c("Control Biofilms", "Treatment Biofilms", "Salmon Carcass", "Unknown")
pie(Sums_icb, labels=lables_icb, main="Control Baetis", col=four_col_vec)

#Pie graphs for insect treatment

#subset insect treatment
HC_16S_SourceTracker_map_it<-subset(HC_16S_SourceTracker_map, Env=="InsectSalmon")
HC_16S_SourceTracker_map_it<-HC_16S_SourceTracker_map_it[,18:ncol(HC_16S_SourceTracker_map_it)]
Sums_it<-colSums(HC_16S_SourceTracker_map_it)
lables_it<-c("Control Biofilms", "Treatment Biofilms", "Salmon Carcass", "Unknown")
pie(Sums_it, labels=lables_it, main="All Treatment Insects", col=four_col_vec)
#subset by insect taxa in control
#Rhyacophila
HC_16S_SourceTracker_map_itr<-subset(HC_16S_SourceTracker_map, Env=="InsectSalmon" & Source=="Rhyacophila")
HC_16S_SourceTracker_map_itr<-HC_16S_SourceTracker_map_itr[,18:ncol(HC_16S_SourceTracker_map_itr)]
Sums_itr<-colSums(HC_16S_SourceTracker_map_itr)
lables_itr<-c("Control Biofilms", "Treatment Biofilms", "Salmon Carcass", "Unknown")
pie(Sums_itr, labels=lables_itr, main="Treatment Rhyacophila", col=four_col_vec)
#Stegopterna
HC_16S_SourceTracker_map_its<-subset(HC_16S_SourceTracker_map, Env=="InsectSalmon" & Source=="Stegopterna")
HC_16S_SourceTracker_map_its<-HC_16S_SourceTracker_map_its[,18:ncol(HC_16S_SourceTracker_map_its)]
Sums_its<-colSums(HC_16S_SourceTracker_map_its)
lables_its<-c("Control Biofilms", "Treatment Biofilms", "Salmon Carcass", "Unknown")
pie(Sums_its, labels=lables_its, main="Treatment Stegopterna", col=four_col_vec)
#Baetis
HC_16S_SourceTracker_map_itb<-subset(HC_16S_SourceTracker_map, Env=="InsectSalmon" & Source=="Baetis")
HC_16S_SourceTracker_map_itb<-HC_16S_SourceTracker_map_itb[,18:ncol(HC_16S_SourceTracker_map_itb)]
Sums_itb<-colSums(HC_16S_SourceTracker_map_itb)
lables_itb<-c("Control Biofilms", "Treatment Biofilms", "Salmon Carcass", "Unknown")
pie(Sums_itb, labels=lables_itb, main="Treatment Baetis", col=four_col_vec)
