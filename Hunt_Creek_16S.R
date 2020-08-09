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
library(VennDiagram)
library(devtools)
library(gplots)
library(venn)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(lattice)
library(gridExtra)
library(pastecs)
library(MASS)
library(extrafont)
font_import()
#Functions
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

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
BFColor<-"#66c2a5"
two_col_vec <- c("black", "bisque3")
two_col_vec_reach<-c("skyblue3", "tomato3")
two_col_vec_taxa<-c("#33a02c", "#1f78b4")
two_col_vec_bfcar<-c("#66c2a5", "#fc8d62")
two_col_vec_baecar<-c("#8da0cb","#fc8d62")
two_col_vec_carsi<-c("#fc8d62","#666666")
three_col_vec<- c("#1b9e77", "#d95f02", "#fc8d62")
three_col_vec_bfcari<-c("#8da0cb","#66c2a5","#fc8d62")
three_col_vec_conttreatcar<-c("#984ea3","#1b9e77","#fc8d62")
three_col_vec_CBaeS<-c("#fc8d62","#8da0cb","#666666")
four_col_vec<-c("#7fc97f", "#beaed4", "#fdc086", "#ffff99")
four_col_vec_BfCBS<-c("#66c2a5","#fc8d62","#8da0cb","#666666")
four_col_vec_BBfCS<-c("#8da0cb","#66c2a5","#fc8d62","#666666")
four_col_vec_CBfBS<-c("#fc8d62","#66c2a5","#8da0cb","#666666")
four_col_vec_y1conttreatcar<-c("#e6ab02","#984ea3","#1b9e77","#fc8d62")
four_col_vec_CBf1SBf2SBf2C<-c("#fc8d62","#034e7b","#66c2a5","#f7fcb9")
four_col_vec_CBa1SBa2SBa2C<-c("#fc8d62","#4a1486","#8da0cb","#005824")
four_col_vec_CS1SS2SS2C<-c("#fc8d62","#386cb0","#666666","#ffff99")
five_col_vec<- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")
five_col_vec_babchs<-c("#66c2a5","#fc8d62","#8da0cb","#f0027f","#666666")
five_col_vec_CBBaSH<-c("#fc8d62","#66c2a5","#8da0cb","#666666","#f0027f")
six_col_vec<- c("#e41a1c", "#377eb8", "green", "#984ea3", "#ff7f00", "#ffff33")
six_col_vec_cont<-c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")
seven_col_vec<-c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")
blue_scale_white_zero<-brewer.pal(9, "Blues")
blue_scale_white_zero[1]<-"white"
#First work with rarefication plots

#Shannon diversity
HC_Shannon_q2<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_shannon_q2.txt", sep="\t", header=T)
rownames(HC_Shannon_q2)<-HC_Shannon_q2[,1]
HC_Shannon_q2[,1]<-NULL
HC_Shannon<-data.frame(t(HC_Shannon_q2))
HC_Shannon$iteration<-NULL
str(HC_Shannon)
HC_Sh_m<-melt(HC_Shannon, id.vars="sequences.per.sample")
any(is.na(HC_Sh_m))
HC_Sh_c_m<-cast(HC_Sh_m, variable ~ sequences.per.sample, fun.aggregate=mean, na.omit=TRUE)
HC_Sh_c_m$Calculation<-rep("mean",175)
HC_Sh_c_var<-cast(HC_Sh_m, variable ~ sequences.per.sample, fun.aggregate=var)
HC_Sh_c_var$Calculation<-rep("variance",175)
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
any(is.na(HC_Sh_map))
HC_Sh_map$Sequences_per_sample<-as.numeric(as.character(HC_Sh_map$Sequences_per_sample))
HC_Sh_map$Source<-factor(revalue(HC_Sh_map$Source, c("Stegopterna"="S. mutata", "Baetis"="B. brunneicolor", "Heptagenia"="H. flavescens")), levels =c("Biofilm", "Carcass", "B. brunneicolor", "H. flavescens", "S. mutata"))
HC_Sh_map_sum_m <- summarySE(HC_Sh_map, measurevar=c("mean"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_Sh_map_sum_v <- summarySE(HC_Sh_map, measurevar=c("variance"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_Sh_map_sum_v$StandDev<-sqrt(HC_Sh_map_sum_v$variance)
HC_Sh_map_sum_v$StandEr<-HC_Sh_map_sum_v$StandDev/sqrt(HC_Sh_map_sum_v$N)
HC_Sh_sum_m_sd<-merge(HC_Sh_map_sum_m,HC_Sh_map_sum_v, by=0)
#make rarefication plots
ggplot(HC_Sh_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Number of reads sampled") +
  ylab("Mean Shannon H' diversity") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=five_col_vec_babchs)

#Faith's diversity
HC_Faith_q2<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_Faith_q2.txt", sep="\t", header=T)
rownames(HC_Faith_q2)<-HC_Faith_q2[,1]
HC_Faith_q2[,1]<-NULL
HC_Faith<-data.frame(t(HC_Faith_q2))
HC_Faith$iteration<-NULL
str(HC_Faith)
HC_F_m<-melt(HC_Faith, id.vars="sequences.per.sample")
HC_F_c_m<-cast(HC_F_m, variable ~ sequences.per.sample, fun.aggregate=mean)
HC_F_c_m$Calculation<-rep("mean",175)
HC_F_c_var<-cast(HC_F_m, variable ~ sequences.per.sample, fun.aggregate=var)
HC_F_c_var$Calculation<-rep("variance",175)
HC_F_c_m_var<-rbind(HC_F_c_m,HC_F_c_var)
names(HC_F_c_m_var)[names(HC_F_c_m_var)=="variable"] <- "SampleID"
HC_F_c_m_var$Calculation<-as.factor(HC_F_c_m_var$Calculation)
HC_F_c_m_var<-as.data.frame(HC_F_c_m_var)
HC_F_c_m_var_m<-melt(HC_F_c_m_var, id.vars=c("Calculation","SampleID"))
HC_F_c_m_var_c<-cast(HC_F_c_m_var_m, variable + SampleID ~ Calculation)
#Merge metadata onto rarefication file
HC_F_map <-merge(Hunt_Creek_16S_map, HC_F_c_m_var_c, by="SampleID")
names(HC_F_map)[names(HC_F_map)=="variable"] <- "Sequences_per_sample"
HC_F_map$Sequences_per_sample<-as.numeric(as.character(HC_F_map$Sequences_per_sample))
HC_F_map$Source<-factor(revalue(HC_F_map$Source, c("Stegopterna"="S. mutata", "Baetis"="B. brunneicolor", "Heptagenia"="H. flavescens")), levels =c("Biofilm", "Carcass", "B. brunneicolor", "H. flavescens", "S. mutata"))
HC_F_map_sum_m <- summarySE(HC_F_map, measurevar=c("mean"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_F_map_sum_v <- summarySE(HC_F_map, measurevar=c("variance"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_F_map_sum_v$StandDev<-sqrt(HC_F_map_sum_v$variance)
HC_F_map_sum_v$StandEr<-HC_F_map_sum_v$StandDev/sqrt(HC_F_map_sum_v$N)
HC_F_sum_m_sd<-merge(HC_F_map_sum_m,HC_F_map_sum_v, by=0)
#make rarefication plots
ggplot(HC_F_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Number of reads sampled") +
  ylab("Mean Faith's phylogenetic diversity") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=five_col_vec_babchs)

#Observed OTUs
HC_Ob_q2<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_Obs_q2.txt", sep="\t", header=T)
rownames(HC_Ob_q2)<-HC_Ob_q2[,1]
HC_Ob_q2[,1]<-NULL
HC_Ob<-data.frame(t(HC_Ob_q2))
HC_Ob$iteration<-NULL
str(HC_Ob)
HC_o_m<-melt(HC_Ob, id.vars="sequences.per.sample")
HC_o_c_m<-cast(HC_o_m, variable ~ sequences.per.sample, fun.aggregate=mean)
HC_o_c_m$Calculation<-rep("mean",175)
HC_o_c_var<-cast(HC_o_m, variable ~ sequences.per.sample, fun.aggregate=var)
HC_o_c_var$Calculation<-rep("variance",175)
HC_o_c_m_var<-rbind(HC_o_c_m,HC_o_c_var)
names(HC_o_c_m_var)[names(HC_o_c_m_var)=="variable"] <- "SampleID"
HC_o_c_m_var$Calculation<-as.factor(HC_o_c_m_var$Calculation)
HC_o_c_m_var<-as.data.frame(HC_o_c_m_var)
HC_o_c_m_var_m<-melt(HC_o_c_m_var, id.vars=c("Calculation","SampleID"))
HC_o_c_m_var_c<-cast(HC_o_c_m_var_m, variable + SampleID ~ Calculation)
#Merge metadata onto rarefication file
HC_o_map <-merge(Hunt_Creek_16S_map, HC_o_c_m_var_c, by="SampleID")
names(HC_o_map)[names(HC_o_map)=="variable"] <- "Sequences_per_sample"
HC_o_map$Sequences_per_sample<-as.numeric(as.character(HC_o_map$Sequences_per_sample))
HC_o_map$Source<-factor(revalue(HC_o_map$Source, c("Stegopterna"="S. mutata", "Baetis"="B. brunneicolor", "Heptagenia"="H. flavescens")), levels =c("Biofilm", "Carcass", "B. brunneicolor", "H. flavescens", "S. mutata"))
HC_o_map_sum_m <- summarySE(HC_o_map, measurevar=c("mean"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_o_map_sum_v <- summarySE(HC_o_map, measurevar=c("variance"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_o_map_sum_v$StandDev<-sqrt(HC_o_map_sum_v$variance)
HC_o_map_sum_v$StandEr<-HC_o_map_sum_v$StandDev/sqrt(HC_o_map_sum_v$N)
HC_o_sum_m_sd<-merge(HC_o_map_sum_m,HC_o_map_sum_v, by=0)
#make rarefication plots
ggplot(HC_o_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Number of reads sampled") +
  ylab("Mean observed OTUs") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=five_col_vec_babchs)

#Chao1
HC_Ch_q2<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_Chao1_q2.txt", sep="\t", header=T)
rownames(HC_Ch_q2)<-HC_Ch_q2[,1]
HC_Ch_q2[,1]<-NULL
HC_Ch<-data.frame(t(HC_Ch_q2))
HC_Ch$iteration<-NULL
str(HC_Ch)
HC_c_m<-melt(HC_Ch, id.vars="sequences.per.sample")
HC_c_c_m<-cast(HC_c_m, variable ~ sequences.per.sample, fun.aggregate=mean)
HC_c_c_m$Calculation<-rep("mean",175)
HC_c_c_var<-cast(HC_c_m, variable ~ sequences.per.sample, fun.aggregate=var)
HC_c_c_var$Calculation<-rep("variance",175)
HC_c_c_m_var<-rbind(HC_c_c_m,HC_c_c_var)
names(HC_c_c_m_var)[names(HC_c_c_m_var)=="variable"] <- "SampleID"
HC_c_c_m_var$Calculation<-as.factor(HC_c_c_m_var$Calculation)
HC_c_c_m_var<-as.data.frame(HC_c_c_m_var)
HC_c_c_m_var_m<-melt(HC_c_c_m_var, id.vars=c("Calculation","SampleID"))
HC_c_c_m_var_c<-cast(HC_c_c_m_var_m, variable + SampleID ~ Calculation)
#Merge metadata onto rarefication file
HC_c_map <-merge(Hunt_Creek_16S_map, HC_c_c_m_var_c, by="SampleID")
names(HC_c_map)[names(HC_c_map)=="variable"] <- "Sequences_per_sample"
HC_c_map$Sequences_per_sample<-as.numeric(as.character(HC_c_map$Sequences_per_sample))
HC_c_map$Source<-factor(revalue(HC_c_map$Source, c("Stegopterna"="S. mutata", "Baetis"="B. brunneicolor", "Heptagenia"="H. flavescens")), levels =c("Biofilm", "Carcass", "B. brunneicolor", "H. flavescens", "S. mutata"))
HC_c_map_sum_m <- summarySE(HC_c_map, measurevar=c("mean"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_c_map_sum_v <- summarySE(HC_c_map, measurevar=c("variance"), groupvars=c("Sequences_per_sample","Source"), na.rm=TRUE)
HC_c_map_sum_v$StandDev<-sqrt(HC_c_map_sum_v$variance)
HC_c_map_sum_v$StandEr<-HC_c_map_sum_v$StandDev/sqrt(HC_c_map_sum_v$N)
HC_c_sum_m_sd<-merge(HC_c_map_sum_m,HC_c_map_sum_v, by=0)
#make rarefication plots
ggplot(HC_c_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Number of reads sampled") +
  ylab("Mean Chao1 estimator") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=five_col_vec_babchs)

#####Created facetted graph for manuscript
HC_Sh_sum_m_sd$alpha<-rep("A. Shannon Diveristy",50)
HC_F_sum_m_sd$alpha<-rep("B. Faith's Phylogenetic Diveristy",50)
HC_o_sum_m_sd$alpha<-rep("C. Observed OTUs",50)
HC_c_sum_m_sd$alpha<-rep("D. Chao1 Estimator",50)
HC_sum_m_sd<-rbind(HC_Sh_sum_m_sd,HC_F_sum_m_sd,HC_o_sum_m_sd,HC_c_sum_m_sd)
ggplot(HC_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Source.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Number of Reads Sampled") +
  ylab("Mean Alpha Diversity Value (+/- SEM)") +
  labs(colour = "Source") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(family="Times New Roman", size=36, margin=margin(t=10,r=0,b=0,l=0)),
        axis.title.y=element_text(family="Times New Roman", size=36, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x=element_text(family="Times New Roman", size=18),
        axis.text.y = element_text(family="Times New Roman", size=18),
        legend.position = 0, 
        strip.text.x = element_text(family="Times New Roman", face="bold", size = 24)) +
  scale_color_manual(values=five_col_vec_CBBaSH,                      
                     limits=c("Carcass","Biofilm","B. brunneicolor","S. mutata","H. flavescens"),
                     labels=c("Carcass","Biofilm",expression(paste(italic('B. brunneicolor'))),expression(paste(italic('S. mutata'))),expression(paste(italic('H. flavescens'))))) +
  facet_wrap(~alpha,scales="free_y")
#####################################
#Upload data tables generated in QIIME and manipulate to make "r friendly"
#####################################

#Get 16S OTU table
Hunt_Creek_16S<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_otu_table_q2.txt", sep="\t", header = T)

#Format data frame so the denovo is row name
row.names(Hunt_Creek_16S)<-Hunt_Creek_16S[,1]

#Delete otu id column, now that otu id is row name
Hunt_Creek_16S$ID<-NULL

#transpose
HC_16S_OTU_t<-t(Hunt_Creek_16S)
HC_16S_OTU_t<-data.frame(HC_16S_OTU_t, check.names = FALSE)

#Merge metadata onto data table
Hunt_Creek_16S_map_r<-Hunt_Creek_16S_map
row.names(Hunt_Creek_16S_map_r)<-Hunt_Creek_16S_map_r[,1]
HC_16S_OTU_map <-merge(Hunt_Creek_16S_map_r, HC_16S_OTU_t, by=0)
str(HC_16S_OTU_map)
HC_16S_OTU_t<-HC_16S_OTU_map[,18:ncol(HC_16S_OTU_map)]

#UNI-Do the following if you have unifrac distances ready
#UNI-Upload weighted unifrac distance matrix and remove samples with error
Hunt_Creek_16S_uni <- read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/weighted_unifrac_HC_q2.txt", header=T)
#UNI-Add metadata to unifrac table
Hunt_Creek_16S_uni_map <-merge(Hunt_Creek_16S_map_r, Hunt_Creek_16S_uni, by=0)
row.names(Hunt_Creek_16S_uni_map)<-Hunt_Creek_16S_uni_map[,1]
Hunt_Creek_16S_uni_map<-Hunt_Creek_16S_uni_map[,-c(1)]
Hunt_Creek_16S_uni<-Hunt_Creek_16S_uni_map[17:ncol(Hunt_Creek_16S_uni_map)]
H_C_16S_uni_com_samples<-as.vector(rownames(Hunt_Creek_16S_uni))
#UNI-use output of names to subset columns into rows to make square matrix
Hunt_Creek_16S_uni<-as.matrix(Hunt_Creek_16S_uni)
Hunt_Creek_16S_uni<-subset(Hunt_Creek_16S_uni, select=c(H_C_16S_uni_com_samples))
str(Hunt_Creek_16S_uni)
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

#Analysis of all 16S community data
#############################################

#Overall permanova
adonis(HC_16S_OTU_t ~ Reach*Year*Source*Total_Biofilm_Growth_PostCarcass, data=HC_16S_env, method="jaccard", permutations=999)
#Significant factors are Reach, Year, Source, time, Reach:Source, Year:Source, Year:time, source:time, year:source:time and reach:year:source:time

#UNI-Overall permanova with unifrac distances
adonis(as.dist(Hunt_Creek_16S_uni) ~ Reach*Source*Total_Biofilm_Growth_PostCarcass*Year, data=HC_16S_uni_env, permutations=999)
#Significant factors are reach, source, time, year, reach:source, reach:time, source:time, reach:year, Source:year, Time:Year, reach:source:time, Source:time:Year
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
str(H_C_16S_uni_car_com)

#Create carcass otu table with metadata
HC_16S_OTU_map_car<-subset(HC_16S_OTU_map, Source=="Carcass")

#How many OTUs and gene amplicon sequences for carcasses?
#Delete OTUs with no observations
cols_to_drop_car = c(rep(TRUE, 17), colSums(HC_16S_OTU_map_car[,18:ncol(HC_16S_OTU_map_car)]) > 0)
HC_16S_OTU_map_car<-HC_16S_OTU_map_car[,cols_to_drop_car]
#There are 2819 OTUs (number of variables - 17 metadata variables)
sum(colSums(HC_16S_OTU_map_car[,18:ncol(HC_16S_OTU_map_car)]))
#37,500 total sequence reads

#Subset into environmental and community data frames
HC_16S_OTU_car_env<-HC_16S_OTU_map_car[,1:17]
HC_16S_OTU_car_com<-HC_16S_OTU_map_car[,18:ncol(HC_16S_OTU_map_car)]

#Year 1
HC_16S_OTU_map_car_Y1<-subset(HC_16S_OTU_map_car, Year=="1")

#Delete OTUs with no observations
cols_to_drop_car_Y1 = c(rep(TRUE, 17), colSums(HC_16S_OTU_map_car_Y1[,18:ncol(HC_16S_OTU_map_car_Y1)]) > 0)
HC_16S_OTU_map_car_Y1<-HC_16S_OTU_map_car_Y1[,cols_to_drop_car_Y1]
#There are 686 OTUs in year 1(number of variables - 17 metadata variables)
sum(colSums(HC_16S_OTU_map_car_Y1[,18:ncol(HC_16S_OTU_map_car_Y1)]))
#12,500 total sequence reads
#Year 2
HC_16S_OTU_map_car_Y2<-subset(HC_16S_OTU_map_car, Year=="2")
#Delete OTUs with no observations
cols_to_drop_car_Y2 = c(rep(TRUE, 17), colSums(HC_16S_OTU_map_car_Y2[,18:ncol(HC_16S_OTU_map_car_Y2)]) > 0)
HC_16S_OTU_map_car_Y2<-HC_16S_OTU_map_car_Y2[,cols_to_drop_car_Y2]
#There are 2196 OTUs in year 2(number of variables - 17 metadata variables)
sum(colSums(HC_16S_OTU_map_car_Y2[,18:ncol(HC_16S_OTU_map_car_Y2)]))
#25,000 total sequence reads

#UNI-carcass permanova with unifrac distances
adonis(as.dist(H_C_16S_uni_car_com) ~ Year, data=H_C_16S_uni_car_env, permutations=99999)
#Year is significant p=0.00027, R2 = .32
#NMDS carcass
HC_NMDS_uni_c<-metaMDS(as.dist(H_C_16S_uni_car_com), permutations=99999)
HC_NMDS_uni_c
#stress = 0.07
stressplot(HC_NMDS_uni_c)
ordiplot(HC_NMDS_uni_c, type="n")
with(HC_NMDS_uni_c, points(HC_NMDS_uni_c, display="sites", col = "black", pch=c(15,17)[HC_16S_uni_car_Year]))
with(HC_NMDS_uni_c, legend("topleft", legend=levels(HC_16S_uni_car_Year), bty="n", col="black",pch=c(15,17)))
with(HC_NMDS_uni_c, ordiellipse(HC_NMDS_uni_c, HC_16S_uni_car_Year, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Year 1"))
with(HC_NMDS_uni_c, ordiellipse(HC_NMDS_uni_c, HC_16S_uni_car_Year, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Year 2"))
with(HC_NMDS_uni_c, legend("topright", legend="2D Stress: 0.07", bty="n"))

#upload phyla level info and run indicator analysis for carcass/generate box plots
#Upload phyla level files for each run
Hunt_Creek_16S_P<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_Phyla_q2.txt", sep="\t", header = T)
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
str(HC_16S_P_map)
head(sort(colSums(HC_16S_P_map[18:ncol(HC_16S_P_map)]),decreasing=TRUE))
#proteobacteria most abundant with 179270
sum((colSums(HC_16S_P_map[18:ncol(HC_16S_P_map)])))
#Total number of reachs = 415000
#subset phyla info into carcass
HC_16S_P_map_car<-subset(HC_16S_P_map, Source=="Carcass")
#determine how many phyla
cols_to_drop_car_P = c(rep(TRUE, 17), colSums(HC_16S_P_map_car[,18:ncol(HC_16S_P_map_car)]) > 0)
HC_16S_P_map_car<-HC_16S_P_map_car[,cols_to_drop_car_P]
HC_16S_P_map_car$k__Bacteria.__<-NULL
HC_16S_P_map_car$k__Bacteria.p__<-NULL
#29 total phyla (number of variables - 17 metadata columns)
HC_16S_P_map_car[,1:17]<-sapply(HC_16S_P_map_car[,1:17], as.factor)
HC_16S_P_map_car_env<-HC_16S_P_map_car[,1:17]
HC_16S_P_map_car_com<-HC_16S_P_map_car[,18:ncol(HC_16S_P_map_car)]
#find most abundant phyla
head(sort(colSums(HC_16S_P_map_car_com),decreasing=TRUE))
#Proteobacteria most abundant 20806
#Indicator analysis to see what phyla are driving this change
HC_p_car_indic<-signassoc(HC_16S_P_map_car_com, cluster=HC_16S_P_map_car_env$Year,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_p_car_indic_sig<-subset(HC_p_car_indic, psidak<=0.05)
#3 phyla significant: k__Bacteria.p__Acidobacteria, k__Bacteria.p__Actinobacteria, k__Bacteria.p__Spirochaetes

#upload family level info and run indicator analysis for carcass
#Upload family level files for each run
Hunt_Creek_16S_F<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_Family_q2.txt", sep="\t", header = T)
#Clasify as data.frame
Hunt_Creek_16S_F<-data.frame(Hunt_Creek_16S_F)
#Format data frame so the taxonomy is row name
row.names(Hunt_Creek_16S_F)<-Hunt_Creek_16S_F[,1]
#Delete taxonomy column
Hunt_Creek_16S_F$ID<-NULL
#transpose
HC_16S_F_t<-t(Hunt_Creek_16S_F)
HC_16S_F_t<-data.frame(HC_16S_F_t)
str(HC_16S_F_t)
names(HC_16S_F_t)
#Merge metadata onto data table
HC_16S_F_map <-merge(Hunt_Creek_16S_map_r, HC_16S_F_t, by=0)

#subset family info into carcass
HC_16S_F_map_car<-subset(HC_16S_F_map, Source=="Carcass")
HC_16S_F_map_car$k__Bacteria.__.__.__.__<-NULL
HC_16S_F_map_car$k__Bacteria.p__.c__.o__.f__<-NULL
#determine how many families
cols_to_drop_car_F = c(rep(TRUE, 17), colSums(HC_16S_F_map_car[,18:ncol(HC_16S_F_map_car)]) > 0)
HC_16S_F_map_car<-HC_16S_F_map_car[,cols_to_drop_car_F]
#349 (number of variables - 17 metadata columns)
HC_16S_F_map_car[,1:17]<-sapply(HC_16S_F_map_car[,1:17], as.factor)
HC_16S_F_map_car_env<-HC_16S_F_map_car[,1:17]
HC_16S_F_map_car_com<-HC_16S_F_map_car[,18:ncol(HC_16S_F_map_car)]
#find most abundant family
head(sort(colSums(HC_16S_F_map_car_com),decreasing=TRUE))
#k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Moraxellaceae most abundant 5837
sum(HC_16S_F_map_car_com)
stat.desc(HC_16S_F_map_car_com$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Moraxellaceae)
#Indicator analysis to see what families are driving this change
HC_f_car_indic<-signassoc(HC_16S_F_map_car_com, cluster=HC_16S_F_map_car_env$Year,  mode=0, alternative = "two.sided",control = how(nperm=9999))
HC_f_car_indic_sig<-subset(HC_f_car_indic, psidak<=0.05)
#20 significant families
write.table(HC_f_car_indic_sig, "~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/16SCarIndSigQ2.txt", sep="\t") 
#figure out which of the most significant are most abundant and or what abundance differences there are between years.
HC_16S_F_map_car_agg<-aggregate(HC_16S_F_map_car[18:ncol(HC_16S_F_map_car)], by=list(Year=HC_16S_F_map_car$Year), FUN=sum)
HC_16S_F_map_car_agg$k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__
#357 found first year, 0 in the second
HC_16S_F_map_car_agg$k__Bacteria.p__Cyanobacteria.c__4C0d.2.o__YS2.f__
#303 first year, 0 second
HC_16S_F_map_car_agg$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae
#448 first year, 0 second
HC_16S_F_map_car_agg$k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfuromonadales.f__Geobacteraceae
#181 first year, 0 second
HC_16S_F_map_car_agg$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Aeromonadales.f__Succinivibrionaceae
#64 first year, 0 second
HC_16S_F_map_car_agg$k__Bacteria.p__Spirochaetes.c__Spirochaetes.o__Spirochaetales.f__Spirochaetaceae
#175 first year, 0 second
#only include indicator taxa in abundances
HC_f_car_indic_abund<-subset(HC_16S_F_map_car_agg, select=c(rownames(HC_f_car_indic_sig)))
tail(sort(colSums(HC_f_car_indic_abund)))
#most abundant is k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae 
HC_16S_F_map_car_agg$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae
#2136 first year, 1574 second year
HC_16S_F_map_car_y1<-subset(HC_16S_F_map_car, Year=="1")
HC_16S_F_map_car_y2<-subset(HC_16S_F_map_car, Year=="2")
stat.desc(HC_16S_F_map_car_y1$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae)
stat.desc(HC_16S_F_map_car_y2$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae)

#############################################
#Separate data based on year of study (1 or 2)
###############################################

#UNI-Create Year 1 unifrac distance table with metadata
H_C_16S_uni_map_Y1<-subset(Hunt_Creek_16S_uni_map, Year=="1")
#Subset into environmental and community data frames
H_C_16S_uni_map_Y1_env<-H_C_16S_uni_map_Y1[,1:16]
H_C_16S_uni_map_Y1_com<-H_C_16S_uni_map_Y1[,17:ncol(H_C_16S_uni_map_Y1)]
H_C_16S_uni_Y1_com_samples<-as.vector(rownames(H_C_16S_uni_map_Y1_com))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_uni_map_Y1_com<-as.matrix(H_C_16S_uni_map_Y1_com)
str(H_C_16S_uni_map_Y1_com)
H_C_16S_uni_map_Y1_com<-subset(H_C_16S_uni_map_Y1_com, select=c(H_C_16S_uni_Y1_com_samples))
H_C_16S_uni_map_Y1_Env<-as.factor(H_C_16S_uni_map_Y1_env$Env)

#UNI-Create Year 2 unifrac distance table with metadata
H_C_16S_uni_map_Y2<-subset(Hunt_Creek_16S_uni_map, Year=="2")
#Subset into environmental and community data frames
H_C_16S_uni_map_Y2_env<-H_C_16S_uni_map_Y2[,1:16]
H_C_16S_uni_map_Y2_com<-H_C_16S_uni_map_Y2[,17:ncol(H_C_16S_uni_map_Y2)]
H_C_16S_uni_Y2_com_samples<-as.vector(rownames(H_C_16S_uni_map_Y2_com))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_uni_map_Y2_com<-as.matrix(H_C_16S_uni_map_Y2_com)
str(H_C_16S_uni_map_Y2_com)
H_C_16S_uni_map_Y2_com<-subset(H_C_16S_uni_map_Y2_com, select=c(H_C_16S_uni_Y2_com_samples))
H_C_16S_uni_map_Y2_Source<-factor(H_C_16S_uni_map_Y2_env$Source)
levels(H_C_16S_uni_map_Y2_Source)<-c("B. brunneicolor","Biofilm","Carcass","H. flavescens","S. mutata")
H_C_16S_uni_map_Y2_Source<- ordered(H_C_16S_uni_map_Y2_Source, levels = c("Biofilm", "Carcass", "B. brunneicolor","H. flavescens","S. mutata"))
H_C_16S_uni_map_Y2_Env<-factor(H_C_16S_uni_map_Y2_env$Env)
levels(H_C_16S_uni_map_Y2_Env)<-c("Control","Salmon","Control","Salmon","Salmon", "Salmon", "Control", "Salmon")
#Create year 1 otu table with metadata
HC_16S_OTU_map_Y1<-subset(HC_16S_OTU_map, Year=="1")
H_C_16S_map_Y1_env<-HC_16S_OTU_map_Y1[,1:17]
H_C_16S_map_Y1_com<-HC_16S_OTU_map_Y1[,18:ncol(HC_16S_OTU_map_Y1)]
#Create year2 otu table with metadata
HC_16S_OTU_map_Y2<-subset(HC_16S_OTU_map, Year=="2")
H_C_16S_map_Y2_env<-HC_16S_OTU_map_Y2[,1:17]
H_C_16S_map_Y2_com<-HC_16S_OTU_map_Y2[,18:ncol(HC_16S_OTU_map_Y2)]

#Visualize overall differences for year 1 and year 2 separately
#UNI-year 1 permanova with unifrac distances
adonis(as.dist(H_C_16S_uni_map_Y1_com) ~ Reach*Source*Total_Biofilm_Growth_PostCarcass, data=H_C_16S_uni_map_Y1_env, permutations=999)
#Reach, source, time, reach:source, reach:time are significant, source is most significant
#NMDS year 1
HC_NMDS_uni_y1<-metaMDS(as.dist(H_C_16S_uni_map_Y1_com))
#stress is nearly 0 - it appears that stegopterna samples are throwing this off because they are so different, so try again without those samples
#UNI-Create Year 1 unifrac distance table with metadata without any simuliidae
H_C_16S_uni_map_Y1ns<-subset(Hunt_Creek_16S_uni_map, Year=="1" & Source!="Stegopterna")
#Subset into environmental and community data frames
H_C_16S_uni_map_Y1ns_env<-H_C_16S_uni_map_Y1ns[,1:16]
H_C_16S_uni_map_Y1ns_com<-H_C_16S_uni_map_Y1ns[,17:ncol(H_C_16S_uni_map_Y1ns)]
H_C_16S_uni_Y1ns_com_samples<-as.vector(rownames(H_C_16S_uni_map_Y1ns_com))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_uni_map_Y1ns_com<-as.matrix(H_C_16S_uni_map_Y1ns_com)
str(H_C_16S_uni_map_Y1ns_com)
H_C_16S_uni_map_Y1ns_com<-subset(H_C_16S_uni_map_Y1ns_com, select=c(H_C_16S_uni_Y1ns_com_samples))
H_C_16S_uni_map_Y1ns_Env<-factor(H_C_16S_uni_map_Y1ns_env$Env)
H_C_16S_uni_map_Y1ns_Source<-factor(H_C_16S_uni_map_Y1ns_env$Source)
levels(H_C_16S_uni_map_Y1ns_Source)<-c("B. brunneicolor", "Biofilm", "Carcass")
levels(H_C_16S_uni_map_Y1ns_Env)<-c("Control","Salmon","Control","Salmon","Salmon")
HC_NMDS_uni_y1ns<-metaMDS(as.dist(H_C_16S_uni_map_Y1ns_com))
HC_NMDS_uni_y1ns
#stress=0.13
stressplot(HC_NMDS_uni_y1ns)
ordiplot(HC_NMDS_uni_y1ns, type="n")
with(HC_NMDS_uni_y1ns, points(HC_NMDS_uni_y1ns, display="sites", col=three_col_vec_bfcari[H_C_16S_uni_map_Y1ns_Source], pch=c(15,17)[H_C_16S_uni_map_Y1ns_Env]))
with(HC_NMDS_uni_y1ns, legend("topleft", legend=c(levels(H_C_16S_uni_map_Y1ns_Source),levels(H_C_16S_uni_map_Y1ns_Env)), bty="n", col=c(three_col_vec_bfcari,"black","black"),pch=c(19,19,19,0,2),pt.bg="black"))
with(HC_NMDS_uni_y1ns, legend("topright", legend="2D Stress: 0.13", bty="n"))

#UNI-year 2 permanova with unifrac distances
adonis(as.dist(H_C_16S_uni_map_Y2_com) ~ Reach*Source*Total_Biofilm_Growth_PostCarcass, data=H_C_16S_uni_map_Y2_env, permutations=999)
#Reach, source, time, source:time reach:source:time significant
#NMDS year 2
HC_NMDS_uni_y2<-metaMDS(as.dist(H_C_16S_uni_map_Y2_com))
HC_NMDS_uni_y2
#stress=0.11
stressplot(HC_NMDS_uni_y2)
ordiplot(HC_NMDS_uni_y2, type="n")
with(HC_NMDS_uni_y2, points(HC_NMDS_uni_y2, display="sites", col=five_col_vec_babchs[H_C_16S_uni_map_Y2_Source], pch=c(15,17)[H_C_16S_uni_map_Y2_Env]))
with(HC_NMDS_uni_y2, legend("topleft", legend=c(levels(H_C_16S_uni_map_Y2_Source),levels(H_C_16S_uni_map_Y2_Env)), bty="n", col=c(five_col_vec_babchs,"black","black"),pch=c(19,19,19,19,19,0,2),pt.bg="black"))
with(HC_NMDS_uni_y2, legend("topright", legend="2D Stress: 0.11", bty="n"))

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
#There are 3118 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_map_Y1_BF[,18:ncol(H_C_16S_map_Y1_BF)]))
#90,000 total sequence reads
#without metadata for biofilms in year 1
H_C_16S_Y1_BF<-H_C_16S_map_Y1_BF[,18:ncol(H_C_16S_map_Y1_BF)]
H_C_16S_Y1_BF_samples<-as.vector(rownames(H_C_16S_Y1_BF))
#Biofilm year 1 environmental variable table
H_C_16S_env_Y1_BF<-H_C_16S_map_Y1_BF[,1:17]
Total_Biofilm_Growth_PostCarcass_Y1<-H_C_16S_env_Y1_BF$Total_Biofilm_Growth_PostCarcass
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
Total_Biofilm_Growth_PostCarcass_Y1<-(H_C_16S_env_Y1_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y1<-as.factor(H_C_16S_env_Y1_BF$Reach)
BA_Y1<-as.factor(H_C_16S_env_Y1_BF$BA)
Env_Y1BF<-factor(revalue(H_C_16S_env_Y1_BF$Env, c("BiofilmSalmon"="Salmon", "BiofilmControl"="Control")), levels =c("Salmon", "Control"))


#Biofilm year 2
#Create biofilm matrix for year 2 with metadata
H_C_16S_map_Y2_BF<-subset(HC_16S_OTU_map_Y2, Source=="Biofilm")
#how many OTUs BF Y2?
#Delete OTUs with no observations
cols_to_drop_BF_Y2 = c(rep(TRUE, 17), colSums(H_C_16S_map_Y2_BF[,18:ncol(H_C_16S_map_Y2_BF)]) > 0)
H_C_16S_map_Y2_BF<-H_C_16S_map_Y2_BF[,cols_to_drop_BF_Y2]
#There are 4132 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_map_Y2_BF[,18:ncol(H_C_16S_map_Y2_BF)]))
#90,000 total sequence reads
#table without metadata for biofilms in year 2
H_C_16S_Y2_BF<-H_C_16S_map_Y2_BF[,18:ncol(H_C_16S_map_Y2_BF)]
H_C_16S_Y2_BF_samples<-as.vector(rownames(H_C_16S_Y2_BF))
#Biofilm year 2 environmental variable table
H_C_16S_env_Y2_BF<-H_C_16S_map_Y2_BF[,1:17]
Total_Biofilm_Growth_PostCarcass_Y2<-H_C_16S_env_Y2_BF$Total_Biofilm_Growth_PostCarcass
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
Total_Biofilm_Growth_PostCarcass_Y2<-H_C_16S_env_Y2_BF$Total_Biofilm_Growth_PostCarcass
Reach_Y2<-as.factor(H_C_16S_env_Y2_BF$Reach)
BA_Y2<-as.factor(H_C_16S_env_Y2_BF$BA)
Env_Y2BF<-factor(revalue(H_C_16S_env_Y2_BF$Env, c("BiofilmSalmon"="Salmon", "BiofilmControl"="Control")), levels =c("Salmon", "Control"))

#Year 1 biofilm permanova using counts
adonis(H_C_16S_Y1_BF ~ Reach_Y1*Total_Biofilm_Growth_PostCarcass_Y1, data=H_C_16S_env_Y1_BF, method="jaccard", permutations=999)
#Reach,time and reach:time interaction all significant

#Year 2 biofilm permanova using counts
adonis(H_C_16S_Y2_BF ~ Reach_Y2*Total_Biofilm_Growth_PostCarcass_Y2, data=H_C_16S_env_Y2_BF, method="jaccard", permutations=999)
#time significant

#UNI-Year 1 biofilm permanova using unfrac distance
adonis(as.dist(H_C_16S_uni_Y1_BF) ~ Reach_Y1*Total_Biofilm_Growth_PostCarcass_Y1, data=H_C_16S_env_Y1_BF, permutations=9999)
#Reach,time and reach:time interaction all significant

#UNI-Year 2 biofilm permanova using unifrac distance
adonis(as.dist(H_C_16S_uni_Y2_BF) ~ Reach_Y2*Total_Biofilm_Growth_PostCarcass_Y2, data=H_C_16S_env_Y2_BF, permutations=999)
#time significant 

#Plot NMDS results for biofilms for each year
HC_NMDS_uni_y1bf<-metaMDS(as.dist(H_C_16S_uni_Y1_BF))
HC_NMDS_uni_y1bf
#stress=0.11
stressplot(HC_NMDS_uni_y1bf)
ordiplot(HC_NMDS_uni_y1bf, type="n")
with(HC_NMDS_uni_y1bf, points(HC_NMDS_uni_y1bf, display="sites", col="black", pch=c(16,5)[Env_Y1BF]))
with(HC_NMDS_uni_y1bf, ordiellipse(HC_NMDS_uni_y1bf, Env_Y1BF, kind="se", conf=0.95, lwd=2, col="black", draw="polygon", show.groups = "Salmon"))
with(HC_NMDS_uni_y1bf, ordiellipse(HC_NMDS_uni_y1bf, Env_Y1BF, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Control"))
with(HC_NMDS_uni_y1bf, legend("topleft", legend=c("Salmon","Control"), bty="n", col="black",pch=c(16,5),pt.bg="black"))
with(HC_NMDS_uni_y1bf, legend("topright", legend="2D Stress: 0.11", bty="n"))

HC_NMDS_uni_y2bf<-metaMDS(as.dist(H_C_16S_uni_Y2_BF))
HC_NMDS_uni_y2bf
#stress=0.12
stressplot(HC_NMDS_uni_y2bf)
ordiplot(HC_NMDS_uni_y2bf, type="n")
with(HC_NMDS_uni_y2bf, points(HC_NMDS_uni_y2bf, display="sites", col="black", pch=c(16,5)[Env_Y2BF]))
with(HC_NMDS_uni_y2bf, ordiellipse(HC_NMDS_uni_y2bf, Env_Y2BF, kind="se", conf=0.95, lwd=2, col="black", draw="polygon", show.groups = "Salmon"))
with(HC_NMDS_uni_y2bf, ordiellipse(HC_NMDS_uni_y2bf, Env_Y2BF, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Control"))
with(HC_NMDS_uni_y2bf, legend("topleft", legend=c("Salmon","Control"), bty="n", col="black",pch=c(16,5),pt.bg="black"))
with(HC_NMDS_uni_y2bf, legend("topright", legend="2D Stress: 0.12", bty="n"))

#indicator taxa analysis for biofilms
#subset phyla info into biofilms for each year
HC_16S_P_map_bf_y1<-subset(HC_16S_P_map, Source=="Biofilm" & Year=="1")
HC_16S_P_map_bf_y1$Unassigned.__<-NULL
HC_16S_P_map_bf_y1$k__Bacteria.__<-NULL
HC_16S_P_map_bf_y1$k__Bacteria.p__<-NULL
HC_16S_P_map_bf_y1[,1:17]<-sapply(HC_16S_P_map_bf_y1[,1:17], as.factor)
HC_16S_P_map_bf_y1_env<-HC_16S_P_map_bf_y1[,1:17]
HC_16S_P_map_bf_y1_com<-HC_16S_P_map_bf_y1[,18:ncol(HC_16S_P_map_bf_y1)]
#find most abundant phyla
head(sort(colSums(HC_16S_P_map_bf_y1_com),decreasing=TRUE))
#Proteobacteria most abundant with 39572

#subset phyla info into biofilms y2
HC_16S_P_map_bf_y2<-subset(HC_16S_P_map, Source=="Biofilm" & Year=="2")
HC_16S_P_map_bf_y2$Unassigned.__<-NULL
HC_16S_P_map_bf_y2$k__Bacteria.__<-NULL
HC_16S_P_map_bf_y2$k__Bacteria.p__<-NULL
HC_16S_P_map_bf_y2[,1:17]<-sapply(HC_16S_P_map_bf_y2[,1:17], as.factor)
HC_16S_P_map_bf_y2_env<-HC_16S_P_map_bf_y2[,1:17]
HC_16S_P_map_bf_y2_com<-HC_16S_P_map_bf_y2[,18:ncol(HC_16S_P_map_bf_y2)]
#find most abundant phyla
head(sort(colSums(HC_16S_P_map_bf_y2_com),decreasing=TRUE))
#Proteobacteria most abundant with 40010

#create clusters based on time and reach
BF_Y1_Reach<-HC_16S_P_map_bf_y1_env$Reach
BF_Y2_Reach<-HC_16S_P_map_bf_y2_env$Reach
BF_Y1_Time<-HC_16S_P_map_bf_y1_env$Total_Biofilm_Growth_PostCarcass
BF_Y2_Time<-HC_16S_P_map_bf_y2_env$Total_Biofilm_Growth_PostCarcass
BF_Y1_RT<-paste(BF_Y1_Reach, BF_Y1_Time)
BF_Y2_RT<-paste(BF_Y2_Reach, BF_Y2_Time)
#Indicator analysis of biofilms year 1 to see what phyla are driving this change
HC_p_bf_y1_indic<-signassoc(HC_16S_P_map_bf_y1_com, cluster=BF_Y1_Reach,  mode=0, alternative = "two.sided",control = how(nperm=9999))
HC_p_bf_y1_indic_sig<-subset(HC_p_bf_y1_indic, psidak<=0.05)
#indicator analysis found 2 phyla: k__Bacteria.p__Proteobacteria, 	k__Bacteria.p__Spirochaetes, 

#subset family info into biofilms for each year
HC_16S_F_map_bf_y1<-subset(HC_16S_F_map, Source=="Biofilm" & Year=="1")
HC_16S_F_map_bf_y1$Unassigned.__.__.__.__<-NULL
HC_16S_F_map_bf_y1$k__Bacteria.__.__.__.__<-NULL
HC_16S_F_map_bf_y1$k__Bacteria.p__.c__.o__.f__<-NULL
HC_16S_F_map_bf_y1[,1:17]<-sapply(HC_16S_F_map_bf_y1[,1:17], as.factor)
HC_16S_F_map_bf_y1_env<-HC_16S_F_map_bf_y1[,1:17]
HC_16S_F_map_bf_y1_com<-HC_16S_F_map_bf_y1[,18:ncol(HC_16S_F_map_bf_y1)]
#find most abundant family
head((sort(colSums(HC_16S_F_map_bf_y1_com),decreasing=TRUE)),n=5)
#k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Stramenopiles.f__   most abundant with 14734 reads
sum(HC_16S_F_map_bf_y1_com)
#87284 sequence reads

#subset family info into biofilms y2
HC_16S_F_map_bf_y2<-subset(HC_16S_F_map, Source=="Biofilm" & Year=="2")
HC_16S_F_map_bf_y2$Unassigned.__.__.__.__<-NULL
HC_16S_F_map_bf_y2$k__Bacteria.__.__.__.__<-NULL
HC_16S_F_map_bf_y2$k__Bacteria.p__.c__.o__.f__<-NULL
HC_16S_F_map_bf_y2[,1:17]<-sapply(HC_16S_F_map_bf_y2[,1:17], as.factor)
HC_16S_F_map_bf_y2_env<-HC_16S_F_map_bf_y2[,1:17]
HC_16S_F_map_bf_y2_com<-HC_16S_F_map_bf_y2[,18:ncol(HC_16S_F_map_bf_y2)]
#find most abundant family
head((sort(colSums(HC_16S_F_map_bf_y2_com),decreasing=TRUE)),n=5)
#k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Stramenopiles.f__  most abundant with 16015 reads
sum(HC_16S_F_map_bf_y2_com)
#88582 sequence reads
HC_16S_F_map_bf<-subset(HC_16S_F_map, Source=="Biofilm")
stat.desc(HC_16S_F_map_bf$k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Stramenopiles.f__)

#create clusters based on time and reach
BF_Y1_Reach<-HC_16S_F_map_bf_y1_env$Reach
BF_Y2_Reach<-HC_16S_F_map_bf_y2_env$Reach
BF_Y1_Time<-HC_16S_F_map_bf_y1_env$Total_Biofilm_Growth_PostCarcass
BF_Y2_Time<-HC_16S_F_map_bf_y2_env$Total_Biofilm_Growth_PostCarcass
BF_Y1_RT<-paste(BF_Y1_Reach, BF_Y1_Time)
BF_Y2_RT<-paste(BF_Y2_Reach, BF_Y2_Time)
BF_Y1_Env<-HC_16S_F_map_bf_y1_env$Env
BF_Y2_Env<-HC_16S_F_map_bf_y2_env$Env
#Indicator analysis of biofilms year 1 to see what families are driving this change
HC_f_bf_y1_indic<-signassoc(HC_16S_F_map_bf_y1_com, cluster=BF_Y1_Env,  mode=0, alternative = "two.sided",control = how(nperm=9999))
HC_f_bf_y1_indic_sig<-subset(HC_f_bf_y1_indic, psidak<=0.05)
write.table(HC_f_bf_y1_indic_sig, "~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/16SBFY1IndSig.txt", sep="\t") 
#indicator analysis found 15 families:

#Plot these families to visualize
HC_16S_F_map_bfcar<-subset(HC_16S_F_map, Source=="Biofilm" | Source== "Carcass")
HC_16S_F_map_bf<-subset(HC_16S_F_map, Source=="Biofilm")
HC_16S_F_map_bf$k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__PK29.f__<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__PK29.f__)
BACP_BF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__PK29.f__", groupvars=c("Days_Since_Study_Start","Reach"))
BACP_BF$Reach<-factor(BACP_BF$Reach, c("Salmon", "Control"))
ggplot(BACP_BF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__PK29.f__)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__PK29.f__-se, ymax=k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__PK29.f__+se), width=.1) +
  geom_point(size=3) +
  geom_line(size=1.5, aes(linetype=Reach))+
  xlab("Days since study start") +
  ylab("Mean unnamed PK29 sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100))

HC_16S_F_map_bf$k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__RB41.f__Ellin6075<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__RB41.f__Ellin6075)
BACREBF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__RB41.f__Ellin6075", groupvars=c("Days_Since_Study_Start","Reach"))
BACREBF$Reach<-factor(BACREBF$Reach, c("Salmon", "Control"))
ggplot(BACREBF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__RB41.f__Ellin6075)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__RB41.f__Ellin6075-se, ymax=k__Bacteria.p__Acidobacteria.c__.Chloracidobacteria..o__RB41.f__Ellin6075+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach)) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean Ellin6075 sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

HC_16S_F_map_bf$k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Chitinophagaceae<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Chitinophagaceae)
BBSSCBF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Chitinophagaceae", groupvars=c("Days_Since_Study_Start","Reach"))
BBSSCBF$Reach<-factor(BBSSCBF$Reach, c("Salmon", "Control"))
ggplot(BBSSCBF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Chitinophagaceae)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Chitinophagaceae-se, ymax=k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Chitinophagaceae+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach)) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean Chitinophagaceae sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

HC_16S_F_map_bf$k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Saprospiraceae<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Saprospiraceae)
BBSSSBF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Saprospiraceae", groupvars=c("Days_Since_Study_Start","Reach"))
BBSSSBF$Reach<-factor(BBSSSBF$Reach, c("Salmon", "Control"))
ggplot(BBSSSBF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Saprospiraceae)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Saprospiraceae-se, ymax=k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Saprospiraceae+se), width=.1, color=BFColor) +
  geom_line(size=1.5, aes(linetype=Reach), color=BFColor) +
  geom_point(size=3, color=BFColor) +
  xlab("Days since study start") +
  ylab("Mean Saprospiraceae sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

HC_16S_F_map_bf$Stramenopiles.f__RA<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Stramenopiles.f__)/2500
BCCS_BF <- summarySE(HC_16S_F_map_bf, measurevar="Stramenopiles.f__RA", groupvars=c("Days_Since_Study_Start","Reach"))
BCCS_BF$Reach<-factor(BCCS_BF$Reach, c("Salmon", "Control"))
ggplot(BCCS_BF, aes(x=Days_Since_Study_Start, y=Stramenopiles.f__RA)) + 
  geom_errorbar(aes(ymin=Stramenopiles.f__RA-se, ymax=Stramenopiles.f__RA+se), width=.1, color=BFColor) +
  geom_line(size=1.5, aes(linetype=Reach), color=BFColor) +
  geom_point(size=3, color=BFColor) +
  xlab("Days since study start") +
  ylab("Mean Stramenopiles Rel. Abund. (±SE)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Hyphomicrobiaceae<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Hyphomicrobiaceae)
BPARHBF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Hyphomicrobiaceae", groupvars=c("Days_Since_Study_Start","Reach"))
BPARHBF$Reach<-factor(BPARHBF$Reach, c("Salmon", "Control"))
ggplot(BPARHBF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Hyphomicrobiaceae)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Hyphomicrobiaceae-se, ymax=k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Hyphomicrobiaceae+se), width=.1, color=BFColor) +
  geom_line(size=1.5, aes(linetype=Reach), color=BFColor) +
  geom_point(size=3, color=BFColor) +
  xlab("Days since study start") +
  ylab("Mean Hyphomicrobiaceae sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

HC_16S_F_map_bfcar$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae<-as.numeric(HC_16S_F_map_bfcar$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae)
BPASSBFC <- summarySE(HC_16S_F_map_bfcar, measurevar="k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae", groupvars=c("Days_Since_Study_Start","Reach","Source"))
BPASSBFC$Reach<-factor(BPASSBFC$Reach, c("Salmon", "Control"))
ggplot(BPASSBFC, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae, color=Source)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae-se, ymax=k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach), data=BPASSBFC[BPASSBFC$Source=="Biofilm",]) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean Sphingomonadaceae sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=two_col_vec_bfcar) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")
HC_16S_F_map_bf_sal<-subset(HC_16S_F_map_bf, Reach=="Salmon")
HC_16S_F_map_bf_cont<-subset(HC_16S_F_map_bf, Reach=="Control")
stat.desc(HC_16S_F_map_bf_sal$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae)
stat.desc(HC_16S_F_map_bf_cont$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae)


HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__.f__<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__.f__)
BPB__BF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__.f__", groupvars=c("Days_Since_Study_Start","Reach"))
BPB__BF$Reach<-factor(BPB__BF$Reach, c("Salmon", "Control"))
ggplot(BPB__BF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__.f__)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__.f__-se, ymax=k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__.f__+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach)) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean unnamed Betaproteobacteria sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

HC_16S_F_map_bfcar$k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfuromonadales.f__Geobacteraceae<-as.numeric(HC_16S_F_map_bfcar$k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfuromonadales.f__Geobacteraceae)
BPDDGBFC <- summarySE(HC_16S_F_map_bfcar, measurevar="k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfuromonadales.f__Geobacteraceae", groupvars=c("Days_Since_Study_Start","Reach","Source"))
BPDDGBFC$Reach<-factor(BPDDGBFC$Reach, c("Salmon", "Control"))
ggplot(BPDDGBFC, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfuromonadales.f__Geobacteraceae, color=Source)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfuromonadales.f__Geobacteraceae-se, ymax=k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfuromonadales.f__Geobacteraceae+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach), data=BPDDGBFC[BPDDGBFC$Source=="Biofilm",]) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean Geobacteraceae sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=two_col_vec_bfcar) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")
stat.desc(HC_16S_F_map_bf_sal$k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfuromonadales.f__Geobacteraceae)
stat.desc(HC_16S_F_map_bf_cont$k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfuromonadales.f__Geobacteraceae)

HC_16S_F_map_bfcar$k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales.f__Helicobacteraceae<-as.numeric(HC_16S_F_map_bfcar$k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales.f__Helicobacteraceae)
BPECHBFC <- summarySE(HC_16S_F_map_bfcar, measurevar="k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales.f__Helicobacteraceae", groupvars=c("Days_Since_Study_Start","Reach","Source"))
BPECHBFC$Reach<-factor(BPECHBFC$Reach, c("Salmon", "Control"))
ggplot(BPECHBFC, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales.f__Helicobacteraceae, color=Source)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales.f__Helicobacteraceae-se, ymax=k__Bacteria.p__Proteobacteria.c__Epsilonproteobacteria.o__Campylobacterales.f__Helicobacteraceae+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach), data=BPECHBFC[BPECHBFC$Source=="Biofilm",]) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean Helicobacteraceae sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=two_col_vec_bfcar) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__)
BPG__BF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__", groupvars=c("Days_Since_Study_Start","Reach"))
BPG__BF$Reach<-factor(BPG__BF$Reach, c("Salmon", "Control"))
ggplot(BPG__BF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__-se, ymax=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.__.__+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach)) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean unnamed Gammaproteobacteria sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae)
BPGAABF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae", groupvars=c("Days_Since_Study_Start","Reach"))
BPGAABF$Reach<-factor(BPGAABF$Reach, c("Salmon", "Control"))
ggplot(BPGAABF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae-se, ymax=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae+se), width=.1, color=BFColor) +
  geom_line(size=1.5, aes(linetype=Reach), color=BFColor) +
  geom_point(size=3, color=BFColor) +
  xlab("Days since study start") +
  ylab("Mean Alteromonadaceae sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")
stat.desc(HC_16S_F_map_bf_sal$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae)
stat.desc(HC_16S_F_map_bf_cont$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae)

HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae)
BPGEEBF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae", groupvars=c("Days_Since_Study_Start","Reach"))
BPGEEBF$Reach<-factor(BPGEEBF$Reach, c("Salmon", "Control"))
ggplot(BPGEEBF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae-se, ymax=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae+se), width=.1, color=BFColor) +
  geom_line(size=1.5, aes(linetype=Reach), color=BFColor) +
  geom_point(size=3, color=BFColor) +
  xlab("Days since study start") +
  ylab("Mean Enterobacteriaceae sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")
stat.desc(HC_16S_F_map_bf_sal$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae)
stat.desc(HC_16S_F_map_bf_cont$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacteriales.f__Enterobacteriaceae)

HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Crenotrichaceae<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Crenotrichaceae)
BPGMCBF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Crenotrichaceae", groupvars=c("Days_Since_Study_Start","Reach"))
BPGMCBF$Reach<-factor(BPGMCBF$Reach, c("Salmon", "Control"))
ggplot(BPGMCBF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Crenotrichaceae)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Crenotrichaceae-se, ymax=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Crenotrichaceae+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach)) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean Crenotrichaceae sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

HC_16S_F_map_bfcar$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales.f__Xanthomonadaceae<-as.numeric(HC_16S_F_map_bfcar$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales.f__Xanthomonadaceae)
BPGXXBFC <- summarySE(HC_16S_F_map_bfcar, measurevar="k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales.f__Xanthomonadaceae", groupvars=c("Days_Since_Study_Start","Reach","Source"))
BPGXXBFC$Reach<-factor(BPGXXBFC$Reach, c("Salmon", "Control"))
ggplot(BPGXXBFC, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales.f__Xanthomonadaceae, color=Source)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales.f__Xanthomonadaceae-se, ymax=k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales.f__Xanthomonadaceae+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach), data=BPGXXBFC[BPGXXBFC$Source=="Biofilm",]) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean Xanthomonadaceae sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=two_col_vec_bfcar) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")
stat.desc(HC_16S_F_map_bf_sal$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales.f__Xanthomonadaceae)
stat.desc(HC_16S_F_map_bf_cont$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales.f__Xanthomonadaceae)


HC_16S_F_map_bf$k__Bacteria.p__Verrucomicrobia.c__.Pedosphaerae..o__.Pedosphaerales..f__<-as.numeric(HC_16S_F_map_bf$k__Bacteria.p__Verrucomicrobia.c__.Pedosphaerae..o__.Pedosphaerales..f__)
BVPP_BF <- summarySE(HC_16S_F_map_bf, measurevar="k__Bacteria.p__Verrucomicrobia.c__.Pedosphaerae..o__.Pedosphaerales..f__", groupvars=c("Days_Since_Study_Start","Reach"))
BVPP_BF$Reach<-factor(BVPP_BF$Reach, c("Salmon", "Control"))
ggplot(BVPP_BF, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Verrucomicrobia.c__.Pedosphaerae..o__.Pedosphaerales..f__)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Verrucomicrobia.c__.Pedosphaerae..o__.Pedosphaerales..f__-se, ymax=k__Bacteria.p__Verrucomicrobia.c__.Pedosphaerae..o__.Pedosphaerales..f__+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach)) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean unnamed Pedoshaerales sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

#Create facetted graph
BBSSSBF$Taxa<-rep("B Saprospiraceae",24)
BBSSSBF$Source<-rep("Biofilm",24)
colnames(BBSSSBF)[colnames(BBSSSBF)=="k__Bacteria.p__Bacteroidetes.c__.Saprospirae..o__.Saprospirales..f__Saprospiraceae"] <- "Mean"
BCCS_BF$Taxa<-rep("C Unnamed Stramenopiles",24)
BCCS_BF$Source<-rep("Biofilm",24)
colnames(BCCS_BF)[colnames(BCCS_BF)=="k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Stramenopiles.f__"] <- "Mean"
BPASSBFC$Taxa<-rep("D Sphingomonadaceae",26)
colnames(BPASSBFC)[colnames(BPASSBFC)=="k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__Sphingomonadaceae"] <- "Mean"
BPGXXBFC$Taxa<-rep("E Xanthomonadaceae",26)
colnames(BPGXXBFC)[colnames(BPGXXBFC)=="k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Xanthomonadales.f__Xanthomonadaceae"] <- "Mean"
BPGAABF$Taxa<-rep("A Alteromonadaceae",24)
BPGAABF$Source<-rep("Biofilm",24)
colnames(BPGAABF)[colnames(BPGAABF)=="k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Alteromonadaceae"] <- "Mean"
BPDDGBFC$Taxa<-rep("F Geobacteraceae",26)
colnames(BPDDGBFC)[colnames(BPDDGBFC)=="k__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__Desulfuromonadales.f__Geobacteraceae"] <- "Mean"
manuscripttaxa<-rbind(BBSSSBF,BCCS_BF,BPASSBFC,BPGXXBFC,BPGAABF,BPDDGBFC)
ggplot(manuscripttaxa, aes(x=Days_Since_Study_Start, y=Mean, color=Source)) + 
  geom_point(size=1.5)+
  geom_errorbar(aes(ymin=Mean-se, ymax=Mean+se), width=30, size=0.5) +
  geom_line(size=0.7, aes(linetype=Reach), data=manuscripttaxa[manuscripttaxa$Source=="Biofilm",]) +
  xlab("Days Since Study Start") +
  ylab("Mean Relative Sequence Abundance (+/- SEM)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(family="Helvetica", size=12),
        axis.title.y=element_text(family="Helvetica", size=12),
        axis.text.x=element_text(family="Helvetica", size=8),
        axis.text.y = element_text(family="Helvetica", size=8),
        legend.position = "none",
        strip.text.x = element_text(family="Helvetica", face="bold", size = 8)) +
  scale_x_continuous(breaks=seq(0,700,200)) + 
  scale_color_manual(values=two_col_vec_bfcar) +
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment") +
  facet_wrap( ~ Taxa, ncol=3, scales="free_y")
ggsave(file="Ind_Tax_BF.tiff",device ="tiff",height=4, width=6, units='in', dpi=300)
#########internal insect by year by reach########
#Insects

#Baetis bruniecolor
#Create baetis matrix with metadata
H_C_16S_map_Baetis<-subset(Hunt_Creek_16S_uni_map, Source=="Baetis")
#without metadata for baetis
H_C_16S_Baetis<-H_C_16S_map_Baetis[,18:ncol(H_C_16S_map_Baetis)]
#Baetis environmental variable table
H_C_16S_env_Baetis<-H_C_16S_map_Baetis[,1:17]
Total_Biofilm_Growth_PostCarcass_Bae<-H_C_16S_env_Baetis$Total_Biofilm_Growth_PostCarcass
Reach_Baetis<-as.factor(H_C_16S_env_Baetis$Reach)
Source_Baetis<-as.factor(H_C_16S_env_Baetis$Source)

#Create baetis matrix with metadata
H_C_16S_OTU_map_Baetis<-subset(HC_16S_OTU_map, Source=="Baetis")
#without metadata for baetis
H_C_16S_OTU_Baetis<-H_C_16S_OTU_map_Baetis[,18:ncol(H_C_16S_OTU_map_Baetis)]
#Baetis environmental variable table
H_C_16S_OTU_env_Baetis<-H_C_16S_OTU_map_Baetis[,1:17]
Tot_PostCarcass_OTU_Bae<-H_C_16S_OTU_env_Baetis$Total_Biofilm_Growth_PostCarcass
Reach_OTU_Baetis<-as.factor(H_C_16S_OTU_env_Baetis$Reach)
Source_OTU_Baetis<-as.factor(H_C_16S_OTU_env_Baetis$Source)

#split up by year

#Create baetis matrix for year 1 with metadata
H_C_16S_map_Y1_Bae<-subset(H_C_16S_OTU_map_Baetis, Year=="1")
#how many OTUs Baetis Y1?
#Delete OTUs with no observations
cols_to_drop_Bae_Y1 = c(rep(TRUE, 17), colSums(H_C_16S_map_Y1_Bae[,18:ncol(H_C_16S_map_Y1_Bae)]) > 0)
H_C_16S_map_Y1_Bae<-H_C_16S_map_Y1_Bae[,cols_to_drop_Bae_Y1]
#There are 1431 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_map_Y1_Bae[,18:ncol(H_C_16S_map_Y1_Bae)]))
#42,500 total sequence reads

#Create baetis matrix for year 2 with metadata
H_C_16S_map_Y2_Bae<-subset(H_C_16S_OTU_map_Baetis, Year=="2")
#how many OTUs Baetis Y2?
#Delete OTUs with no observations
cols_to_drop_Bae_Y2 = c(rep(TRUE, 17), colSums(H_C_16S_map_Y2_Bae[,18:ncol(H_C_16S_map_Y2_Bae)]) > 0)
H_C_16S_map_Y2_Bae<-H_C_16S_map_Y2_Bae[,cols_to_drop_Bae_Y2]
#There are 1335 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_map_Y2_Bae[,18:ncol(H_C_16S_map_Y2_Bae)]))
#75,000 total sequence reads

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
Total_Biofilm_Growth_PostCarcass_bae_y1<-H_C_16S_env_Bae_y1$Total_Biofilm_Growth_PostCarcass
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
Total_Biofilm_Growth_PostCarcass_bae_y2<-H_C_16S_env_Bae_y2$Total_Biofilm_Growth_PostCarcass
Reach_Bae_y2<-as.factor(H_C_16S_env_Bae_y2$Reach)
Source_Bae_y2<-as.factor(H_C_16S_env_Bae_y2$Source)

#Baetis permanova using unifrac
adonis(as.dist(H_C_16S_Bae_y1) ~ Reach*Total_Biofilm_Growth_PostCarcass_bae_y1, data=H_C_16S_env_Bae_y1,permutations=9999)
#nothing significant
adonis(as.dist(H_C_16S_Bae_y2) ~ Reach*Total_Biofilm_Growth_PostCarcass_bae_y2, data=H_C_16S_env_Bae_y2,permutations=999)
#time significant

HC_Baetis_y1_NMDS<-metaMDS(as.dist(H_C_16S_Bae_y1))
#insufficient data 

HC_Baetis_y2_NMDS<-metaMDS(as.dist(H_C_16S_Bae_y2))
ordiplot(HC_Baetis_y2_NMDS, type="n", main="Baetis brunicolor 16S year 2")
with(HC_Baetis_y2_NMDS, points(HC_Baetis_y2_NMDS, display="sites", col=two_col_vec_reach[Reach_Bae_y2], pch=19, pt.bg=two_col_vec_reach))
with(HC_Baetis_y2_NMDS, legend("topleft", legend=levels(Reach_Bae_y2), bty="n", col=two_col_vec_reach, pch=19, pt.bg=two_col_vec_reach))
with(HC_Baetis_y2_NMDS, ordiellipse(HC_Baetis_y2_NMDS, Reach_Bae_y2, kind="se", conf=0.95, lwd=2, col="skyblue3", show.groups = "Control"))
with(HC_Baetis_y2_NMDS, ordiellipse(HC_Baetis_y2_NMDS, Reach_Bae_y2, kind="se", conf=0.95, lwd=2, col="tomato3", show.groups = "Salmon"))

#subset phyla info
HC_16S_P_map_bae<-subset(HC_16S_P_map, Source=="Baetis")
HC_16S_P_map_bae[,1:17]<-sapply(HC_16S_P_map_bae[,1:17], as.factor)

#Year 1
HC_16S_P_map_bae_y1<-subset(HC_16S_P_map_bae, Year=="1")
HC_16S_P_map_bae_y1_com<-HC_16S_P_map_bae_y1[,18:ncol(HC_16S_P_map_bae_y1)]
#find most abundant phyla
head(sort(colSums(HC_16S_P_map_bae_y1_com),decreasing=TRUE))
#Proteobacteria most abundant
HC_16S_P_map_bae_y1_env<-HC_16S_P_map_bae_y1[1:17]

#Year 2
HC_16S_P_map_bae_y2<-subset(HC_16S_P_map_bae, Year=="2")
HC_16S_P_map_bae_y2_com<-HC_16S_P_map_bae_y2[,18:ncol(HC_16S_P_map_bae_y2)]
#find most abundant phyla
head(sort(colSums(HC_16S_P_map_bae_y2_com),decreasing=TRUE))
#Proteobacteria most abundant
HC_16S_P_map_bae_y2_env<-HC_16S_P_map_bae_y2[1:17]

#subset family info
HC_16S_F_map_bae<-subset(HC_16S_F_map, Source=="Baetis")
HC_16S_F_map_bae[,1:17]<-sapply(HC_16S_F_map_bae[,1:17], as.factor)

#Year 1
HC_16S_F_map_bae_y1<-subset(HC_16S_F_map_bae, Year=="1")
HC_16S_F_map_bae_y1_com<-HC_16S_F_map_bae_y1[,18:ncol(HC_16S_F_map_bae_y1)]
#find most abundant family
head(sort(colSums(HC_16S_F_map_bae_y1_com),decreasing=TRUE))
#Unnamed Mollicutes most abundant 5508
sum(HC_16S_F_map_bae_y1_com)
#42500
HC_16S_F_map_bae_y1_env<-HC_16S_F_map_bae_y1[1:17]
stat.desc(HC_16S_F_map_bae_y1_com$k__Bacteria.p__Tenericutes.c__Mollicutes.o__.f__)
#Year 2
HC_16S_F_map_bae_y2<-subset(HC_16S_F_map_bae, Year=="2")
HC_16S_F_map_bae_y2_com<-HC_16S_F_map_bae_y2[,18:ncol(HC_16S_F_map_bae_y2)]
#find most abundant family
head(sort(colSums(HC_16S_F_map_bae_y2_com),decreasing=TRUE))
#Pseudomonadaceae most abundant 17356
sum(HC_16S_F_map_bae_y2_com)
#75000
HC_16S_F_map_bae_y2_env<-HC_16S_F_map_bae_y2[1:17]
stat.desc(HC_16S_F_map_bae_y2_com$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Pseudomonadaceae)

#indicator species analysis for year two
#create clusters based on time and reach
Bae_F_Y2_Env<-HC_16S_F_map_bae_y2_env$Env

#Indicator analysis of baetids year 2 to see what families are driving this change
HC_f_bae_y2_indic<-signassoc(HC_16S_F_map_bae_y2_com, cluster=Bae_F_Y2_Env,  mode=0, alternative = "two.sided",control = how(nperm=99999))
HC_f_bae_y2_indic_sig<-subset(HC_f_bae_y2_indic, psidak<=0.05)
#indicator analysis found six  families
HC_16S_F_map_bae$Env<-as.factor(HC_16S_F_map_bae$Env)
HC_16S_F_map_bae_y2_sal<-subset(HC_16S_F_map_bae, Env=="BaetisSalmon" & Year=="2")
HC_16S_F_map_bae_y2_cont<-subset(HC_16S_F_map_bae, Env=="BaetisControl" & Year=="2")
stat.desc(HC_16S_F_map_bae_y2_sal$k__Bacteria.p__Tenericutes.c__Mollicutes.o__.f__)
stat.desc(HC_16S_F_map_bae_y2_cont$k__Bacteria.p__Tenericutes.c__Mollicutes.o__.f__)
stat.desc(HC_16S_F_map_bae_y2_sal$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__)
stat.desc(HC_16S_F_map_bae_y2_cont$k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__)

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
Total_Biofilm_Growth_PostCarcass_Hep<-H_C_16S_env_Hepta$Total_Biofilm_Growth_PostCarcass
Reach_Hep<-as.factor(H_C_16S_env_Hepta$Reach)

#Create heptagenia matrix with metadata
H_C_16S_OTU_map_Hepta<-subset(HC_16S_OTU_map, Source=="Heptagenia")
#without metadata for heptagenia
H_C_16S_OTU_Hepta<-H_C_16S_OTU_map_Hepta[,18:ncol(H_C_16S_OTU_map_Hepta)]
#Heptagenia environmental variable table
H_C_16S_OTU_env_Hepta<-H_C_16S_OTU_map_Hepta[,1:17]
Tot_PostCarcass_OTU_Hep<-H_C_16S_OTU_env_Hepta$Total_Biofilm_Growth_PostCarcass
Reach_OTU_Hepta<-as.factor(H_C_16S_OTU_env_Hepta$Reach)
Source_OTU_Hepta<-as.factor(H_C_16S_OTU_env_Hepta$Source)

#Don't need to split up by year, because only year 2

#how many OTUs Heptagenia?
#Delete OTUs with no observations
cols_to_drop_Hep = c(rep(TRUE, 17), colSums(H_C_16S_OTU_map_Hepta[,18:ncol(H_C_16S_OTU_map_Hepta)]) > 0)
H_C_16S_OTU_map_Hepta<-H_C_16S_OTU_map_Hepta[,cols_to_drop_Hep]
#There are 299 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_OTU_map_Hepta[,18:ncol(H_C_16S_OTU_map_Hepta)]))
#20,000 total sequence reads

#Heptagenia permanova using unifrac distance for time
adonis(as.dist(H_C_16S_Hepta) ~ Tot_PostCarcass_OTU_Hep, permutations=999)
#p=0.016

#indicator taxa analysis for heptagenia
#subset phyla info
HC_16S_P_map_hep<-subset(HC_16S_P_map, Source=="Heptagenia")
HC_16S_P_map_hep[,1:17]<-sapply(HC_16S_P_map_hep[,1:17], as.factor)
HC_16S_P_map_hep_com<-HC_16S_P_map_hep[,18:ncol(HC_16S_P_map_hep)]
#find most abundant phyla
head(sort(colSums(HC_16S_P_map_hep_com),decreasing=TRUE))
#Proteobacteria most abundant
HC_16S_P_map_hep_env<-HC_16S_P_map_hep[1:17]
Hep_P_Time<-HC_16S_P_map_hep_env$Total_Biofilm_Growth_PostCarcass
HC_p_hep_indic<-signassoc(HC_16S_P_map_hep_com, cluster=Hep_P_Time,  mode=0, alternative = "two.sided",control = how(nperm=9999))
HC_p_hep_indic_sig<-subset(HC_p_hep_indic, psidak<=0.05)
#No significant phyla

#subset family info
HC_16S_F_map_hep<-subset(HC_16S_F_map, Source=="Heptagenia")
HC_16S_F_map_hep$Unassigned.__.__.__.__<-NULL
HC_16S_F_map_hep$k__Bacteria.__.__.__.__<-NULL
HC_16S_F_map_hep$k__Bacteria.p__.c__.o__.f__<-NULL
HC_16S_F_map_hep[,1:17]<-sapply(HC_16S_F_map_hep[,1:17], as.factor)
HC_16S_F_map_hep_com<-HC_16S_F_map_hep[,18:ncol(HC_16S_F_map_hep)]
#find most abundant family
head(sort(colSums(HC_16S_F_map_hep_com),decreasing=TRUE))
#Enterobacteriaceae most abundant
HC_16S_F_map_hep_env<-HC_16S_F_map_hep[1:17]
Hep_F_Time<-HC_16S_F_map_hep_env$Total_Biofilm_Growth_PostCarcass
HC_f_hep_indic<-signassoc(HC_16S_F_map_hep_com, cluster=Hep_F_Time,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_f_hep_indic_sig<-subset(HC_f_hep_indic, psidak<=0.05)
#No significant families

#Simuliidae
#Create simulidae unifrac matrix with metadata
H_C_16S_map_Simuliidae<-subset(Hunt_Creek_16S_uni_map, Source=="Stegopterna")
#without metadata for simuliidae
H_C_16S_Sim<-H_C_16S_map_Simuliidae[,18:ncol(H_C_16S_map_Simuliidae)]
H_C_16S_Sim_samples<-as.vector(rownames(H_C_16S_Sim))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_Sim<-as.matrix(H_C_16S_Sim)
H_C_16S_Sim<-subset(H_C_16S_Sim, select=c(H_C_16S_Sim_samples))
#Simuliidae environmental variable table
H_C_16S_env_Sim<-H_C_16S_map_Simuliidae[,1:17]
Total_Biofilm_Growth_PostCarcass_Sim<-H_C_16S_env_Sim$Total_Biofilm_Growth_PostCarcass
Reach_Sim<-as.factor(H_C_16S_env_Sim$Reach)
Source_Sim<-as.factor(H_C_16S_env_Sim$Source)

#Create simulidae count matrix with metadata
H_C_16S_OTU_map_Simuliidae<-subset(HC_16S_OTU_map, Source=="Stegopterna")
#without metadata for simuliidae
H_C_16S_OTU_Sim<-H_C_16S_OTU_map_Simuliidae[,18:ncol(H_C_16S_OTU_map_Simuliidae)]
#Simuliidae environmental variable table
H_C_16S_OTU_env_Sim<-H_C_16S_OTU_map_Simuliidae[,1:17]
Total_Biofilm_Growth_PostCarcass_OTU_Sim<-H_C_16S_OTU_env_Sim$Total_Biofilm_Growth_PostCarcass
Reach_OTU_Sim<-as.factor(H_C_16S_OTU_env_Sim$Reach)
Source_OTU_Sim<-as.factor(H_C_16S_OTU_env_Sim$Source)

#split up by year

#Create simuliidae count matrix for year 1 with metadata
H_C_16S_OTU_map_Y1_Sim<-subset(H_C_16S_OTU_map_Simuliidae, Year=="1")
#how many OTUs Simuliidae Y1?
#Delete OTUs with no observations
cols_to_drop_OTU_Sim_Y1 = c(rep(TRUE, 17), colSums(H_C_16S_OTU_map_Y1_Sim[,18:ncol(H_C_16S_OTU_map_Y1_Sim)]) > 0)
H_C_16S_OTU_map_Y1_Sim<-H_C_16S_OTU_map_Y1_Sim[,cols_to_drop_OTU_Sim_Y1]
#There are 449 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_OTU_map_Y1_Sim[,18:ncol(H_C_16S_OTU_map_Y1_Sim)]))
#17,500 total sequence reads

#Create simuliidae count matrix for year 1 with metadata
H_C_16S_OTU_map_Y2_Sim<-subset(H_C_16S_OTU_map_Simuliidae, Year=="2")
#how many OTUs Simuliidae Y2?
#Delete OTUs with no observations
cols_to_drop_OTU_Sim_Y2 = c(rep(TRUE, 17), colSums(H_C_16S_OTU_map_Y2_Sim[,18:ncol(H_C_16S_OTU_map_Y2_Sim)]) > 0)
H_C_16S_OTU_map_Y2_Sim<-H_C_16S_OTU_map_Y2_Sim[,cols_to_drop_OTU_Sim_Y2]
#There are 1189 OTUs (number of variables - 17 metadata variables)
sum(colSums(H_C_16S_OTU_map_Y2_Sim[,18:ncol(H_C_16S_OTU_map_Y2_Sim)]))
#42,500 total sequence reads

#Year 1
H_C_16S_map_Sim_y1<-subset(H_C_16S_map_Simuliidae, Year=="1")
#without metadata for sim
H_C_16S_Sim_y1<-H_C_16S_map_Sim_y1[,18:ncol(H_C_16S_map_Sim_y1)]
H_C_16S_Sim_y1_samples<-as.vector(rownames(H_C_16S_Sim_y1))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_Sim_y1<-as.matrix(H_C_16S_Sim_y1)
H_C_16S_Sim_y1<-subset(H_C_16S_Sim_y1, select=c(H_C_16S_Sim_y1_samples))
#Simuliidae year 1 environmental variable table
H_C_16S_env_Sim_y1<-H_C_16S_map_Sim_y1[,1:17]
Total_Biofilm_Growth_PostCarcass_sim_y1<-H_C_16S_env_Sim_y1$Total_Biofilm_Growth_PostCarcass
Reach_Sim_y1<-as.factor(H_C_16S_env_Sim_y1$Reach)
Source_Sim_y1<-as.factor(H_C_16S_env_Sim_y1$Source)
#Year 2
H_C_16S_map_Sim_y2<-subset(H_C_16S_map_Simuliidae, Year=="2")
#without metadata for simulid
H_C_16S_Sim_y2<-H_C_16S_map_Sim_y2[,18:ncol(H_C_16S_map_Sim_y2)]
H_C_16S_Sim_y2_samples<-as.vector(rownames(H_C_16S_Sim_y2))
#UNI-use output of names to subset columns into rows to make square matrix
H_C_16S_Sim_y2<-as.matrix(H_C_16S_Sim_y2)
H_C_16S_Sim_y2<-subset(H_C_16S_Sim_y2, select=c(H_C_16S_Sim_y2_samples))
#Simuliidae year 1 environmental variable table
H_C_16S_env_Sim_y2<-H_C_16S_map_Sim_y2[,1:17]
Total_Biofilm_Growth_PostCarcass_sim_y2<-H_C_16S_env_Sim_y2$Total_Biofilm_Growth_PostCarcass
Reach_Sim_y2<-as.factor(H_C_16S_env_Sim_y2$Reach)
Source_Sim_y2<-as.factor(H_C_16S_env_Sim_y2$Source)

#Simuliidae permanova using unifrac
adonis(as.dist(H_C_16S_Sim) ~ Reach*Total_Biofilm_Growth_PostCarcass*Year, data=H_C_16S_env_Sim,permutations=999)
#reach, time, year, reach:time all significant
adonis(as.dist(H_C_16S_Sim_y1) ~ Reach*Total_Biofilm_Growth_PostCarcass, data=H_C_16S_env_Sim_y1,permutations=999)
#No significant factors
adonis(as.dist(H_C_16S_Sim_y2) ~ Reach*Total_Biofilm_Growth_PostCarcass, data=H_C_16S_env_Sim_y2,permutations=9999)
#Reach:time significant

#indicator taxa analysis for Simuliidae
#subset phyla info
HC_16S_P_map_sim<-subset(HC_16S_P_map, Source=="Stegopterna")
HC_16S_P_map_sim[,1:17]<-sapply(HC_16S_P_map_sim[,1:17], as.factor)

#Year 1
HC_16S_P_map_sim_y1<-subset(HC_16S_P_map_sim, Year=="1")
HC_16S_P_map_sim_y1$k__Bacteria.__<-NULL
HC_16S_P_map_sim_y1$k__Bacteria.p__<-NULL
HC_16S_P_map_sim_y1_com<-HC_16S_P_map_sim_y1[,18:ncol(HC_16S_P_map_sim_y1)]
#find most abundant phyla
head(sort(colSums(HC_16S_P_map_sim_y1_com),decreasing=TRUE))
#Firmicutes most abundant 6353
sum(HC_16S_P_map_sim_y1_com)
#12688
HC_16S_P_map_sim_y1_env<-HC_16S_P_map_sim_y1[1:17]

#Year 2
HC_16S_P_map_sim_y2<-subset(HC_16S_P_map_sim, Year=="2")
HC_16S_P_map_sim_y1$k__Bacteria.__<-NULL
HC_16S_P_map_sim_y1$k__Bacteria.p__<-NULL
HC_16S_P_map_sim_y2_com<-HC_16S_P_map_sim_y2[,18:ncol(HC_16S_P_map_sim_y2)]
#find most abundant phyla
head(sort(colSums(HC_16S_P_map_sim_y2_com),decreasing=TRUE))
#Firmicutes most abundant 13109
sum(HC_16S_P_map_sim_y2_com)
#42500
HC_16S_P_map_sim_y2_env<-HC_16S_P_map_sim_y2[1:17]
stat.desc(HC_16S_P_map_sim$k__Bacteria.p__Firmicutes)

#create clusters based on time and reach
Sim_Y1_Reach<-HC_16S_P_map_sim_y1_env$Reach
Sim_Y2_Reach<-HC_16S_P_map_sim_y2_env$Reach
Sim_Y1_Time<-HC_16S_P_map_sim_y1_env$Total_Biofilm_Growth_PostCarcass
Sim_Y2_Time<-HC_16S_P_map_sim_y2_env$Total_Biofilm_Growth_PostCarcass
Sim_Y1_RT<-paste(Sim_Y1_Reach, Sim_Y1_Time)
Sim_Y2_RT<-paste(Sim_Y2_Reach, Sim_Y2_Time)
Sim_Y1_Env<-HC_16S_P_map_sim_y1_env$Env
Sim_Y2_Env<-HC_16S_P_map_sim_y2_env$Env

#Indicator analysis of Simuliids year 2 to see what phyla are driving this change
HC_p_sim_y2_indic<-signassoc(HC_16S_P_map_sim_y2_com, cluster=Sim_Y2_Env,  mode=0, alternative = "two.sided",control = how(nperm=9999))
HC_p_sim_y2_indic_sig<-subset(HC_p_sim_y2_indic, psidak<=0.05)
#indicator analysis found no phyla

#subset family info
HC_16S_F_map_sim<-subset(HC_16S_F_map, Source=="Stegopterna")
HC_16S_F_map_sim[,1:17]<-sapply(HC_16S_F_map_sim[,1:17], as.factor)

#Year 1
HC_16S_F_map_sim_y1<-subset(HC_16S_F_map_sim, Year=="1")
HC_16S_F_map_sim_y1$k__Bacteria.__.__.__.__<-NULL
HC_16S_F_map_sim_y1$k__Bacteria.p__.c__.o__.f__<-NULL
HC_16S_F_map_sim_y1_com<-HC_16S_F_map_sim_y1[,18:ncol(HC_16S_F_map_sim_y1)]
#find most abundant family
head(sort(colSums(HC_16S_F_map_sim_y1_com),decreasing=TRUE))
#BAcillaceae most abundant
HC_16S_F_map_sim_y1_env<-HC_16S_F_map_sim_y1[1:17]

#Year 2
HC_16S_F_map_sim_y2<-subset(HC_16S_F_map_sim, Year=="2")
HC_16S_F_map_sim_y2$k__Bacteria.__.__.__.__<-NULL
HC_16S_F_map_sim_y2$k__Bacteria.p__.c__.o__.f__<-NULL
HC_16S_F_map_sim_y2_com<-HC_16S_F_map_sim_y2[,18:ncol(HC_16S_F_map_sim_y2)]
#find most abundant family
head(sort(colSums(HC_16S_F_map_sim_y2_com),decreasing=TRUE))
#Bacillaceae most abundant
HC_16S_F_map_sim_y2_env<-HC_16S_F_map_sim_y2[1:17]

#create clusters based on time and reach
Sim_F_Y1_Reach<-HC_16S_F_map_sim_y1_env$Reach
Sim_F_Y2_Reach<-HC_16S_F_map_sim_y2_env$Reach
Sim_F_Y1_Time<-HC_16S_F_map_sim_y1_env$Total_Biofilm_Growth_PostCarcass
Sim_F_Y2_Time<-HC_16S_F_map_sim_y2_env$Total_Biofilm_Growth_PostCarcass
Sim_F_Y1_RT<-paste(Sim_F_Y1_Reach, Sim_F_Y1_Time)
Sim_F_Y2_RT<-paste(Sim_F_Y2_Reach, Sim_F_Y2_Time)
Sim_F_Y1_Env<-HC_16S_F_map_sim_y1_env$Env
Sim_F_Y2_Env<-HC_16S_F_map_sim_y2_env$Env

#Indicator analysis of Simuliids year 2 to see what families are driving this change
HC_f_sim_y2_indic<-signassoc(HC_16S_F_map_sim_y2_com, cluster=Sim_F_Y2_Env,  mode=0, alternative = "two.sided",control = how(nperm=99999))
HC_f_sim_y2_indic_sig<-subset(HC_f_sim_y2_indic, psidak<=0.05)
#indicator analysis found one family: k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Streptophyta.f__
HC_16S_F_map_sim_y2_agg<-aggregate(HC_16S_F_map_sim_y2[18:ncol(HC_16S_F_map_sim_y2)], by=list(Env=HC_16S_F_map_sim_y2$Env), FUN=sum)
HC_16S_F_map_sim_y2_agg$k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Streptophyta.f__
#Streptophyta over 5 greater
HC_16S_F_map_sim_y2_sal<-subset(HC_16S_F_map_sim_y2, Env=="StegopternaSalmon")
HC_16S_F_map_sim_y2_cont<-subset(HC_16S_F_map_sim_y2, Env=="StegopternaControl")
stat.desc(HC_16S_F_map_sim_y2_sal$k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Streptophyta.f__)
stat.desc(HC_16S_F_map_sim_y2_cont$k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Streptophyta.f__)

#Create figure
#Plot these families to visualize
HC_16S_F_map_sim$k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Streptophyta.f__<-as.numeric(HC_16S_F_map_sim$k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Streptophyta.f__)
HC_16S_F_map_sim$Days_Since_Study_Start<-as.numeric(HC_16S_F_map_sim$Days_Since_Study_Start)
BCCS_Sim <- summarySE(HC_16S_F_map_sim, measurevar="k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Streptophyta.f__", groupvars=c("Days_Since_Study_Start","Reach"))
BCCS_Sim$Reach<-factor(BCCS_Sim$Reach, c("Salmon", "Control"))
ggplot(BCCS_Sim, aes(x=Days_Since_Study_Start, y=k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Streptophyta.f__)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Streptophyta.f__-se, ymax=k__Bacteria.p__Cyanobacteria.c__Chloroplast.o__Streptophyta.f__+se), width=.1) +
  geom_line(size=1.5, aes(linetype=Reach)) +
  geom_point(size=3) +
  xlab("Days since study start") +
  ylab("Mean unamed Streptophyta sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  scale_linetype_discrete(name="Treatment")

######Figure out what taxa are introduced by the carrion and if they persist in biofilms/internal invertebrates########

#Create venn diagram for shared OTUs
HC_16S_OTU_map_ag_venn<-HC_16S_OTU_map
str(HC_16S_OTU_map_ag_venn)
HC_16S_OTU_map_ag_venn<-aggregate(HC_16S_OTU_map_ag_venn[18:ncol(HC_16S_OTU_map_ag_venn)], by=list(Env=HC_16S_OTU_map_ag_venn$Env, Year=HC_16S_OTU_map_ag_venn$Year, Date=HC_16S_OTU_map_ag_venn$Date), FUN=sum)
str(HC_16S_OTU_map_ag_venn)
head(colnames(HC_16S_OTU_map_ag_venn))
HC_16S_OTU_map_ag_venn_m<- melt(HC_16S_OTU_map_ag_venn, id=c("Env", "Year", "Date")) 
str(HC_16S_OTU_map_ag_venn_m)
head(colnames(HC_16S_OTU_map_ag_venn_m))
HC_16S_OTU_map_ag_venn_m$EnvYD<-as.factor(paste(HC_16S_OTU_map_ag_venn_m$Env, HC_16S_OTU_map_ag_venn_m$Year, HC_16S_OTU_map_ag_venn_m$Date))
levels(HC_16S_OTU_map_ag_venn_m$EnvYD)
HC_16S_OTU_map_ag_venn_m$EnvYD<-revalue(HC_16S_OTU_map_ag_venn_m$EnvYD, c("BaetisControl 1 3/12/15"     ="BaetisControl1",
                                                                          "BaetisControl 1 5/12/15"     ="BaetisControl1",
                                                                          "BaetisControl 1 8/4/15"      ="BaetisControl1",
                                                                          "BaetisControl 1 9/20/15"     ="BaetisControl1",
                                                                          "BaetisControl 1 9/5/14"      ="BaetisControl1",
                                                                          "BaetisControl 2 5/17/16"     ="BaetisControl2",
                                                                          "BaetisControl 2 6/30/16"     ="BaetisControl2",
                                                                          "BaetisSalmon 1 3/12/15"      ="BaetisSalmon1",
                                                                          "BaetisSalmon 1 5/12/15"      ="BaetisSalmon1",
                                                                          "BaetisSalmon 1 6/23/15"      ="BaetisSalmon1",
                                                                          "BaetisSalmon 1 8/4/15"       ="BaetisSalmon1",
                                                                          "BaetisSalmon 1 9/20/15"      ="BaetisSalmon1",
                                                                          "BaetisSalmon 2 3/19/16"      ="BaetisSalmon2",
                                                                          "BaetisSalmon 2 5/17/16"      ="BaetisSalmon2",
                                                                          "BaetisSalmon 2 6/30/16"      ="BaetisSalmon2",
                                                                          "BiofilmControl 1 10/18/14"   ="BiofilmControl1",
                                                                          "BiofilmControl 1 10/3/14"    ="BiofilmControl1",
                                                                          "BiofilmControl 1 3/12/15"    ="BiofilmControl1",
                                                                          "BiofilmControl 1 5/12/15"    ="BiofilmControl1",
                                                                          "BiofilmControl 1 6/23/15"    ="BiofilmControl1",
                                                                          "BiofilmControl 1 8/4/15"     ="BiofilmControl1",
                                                                          "BiofilmControl 2 10/25/15"   ="BiofilmControl2",
                                                                          "BiofilmControl 2 10/4/15"    ="BiofilmControl1",
                                                                          "BiofilmControl 2 3/19/16"    ="BiofilmControl2",
                                                                          "BiofilmControl 2 5/17/16"    ="BiofilmControl2",
                                                                          "BiofilmControl 2 6/30/16"    ="BiofilmControl2",
                                                                          "BiofilmControl 2 8/15/16"    ="BiofilmControl2",
                                                                          "BiofilmSalmon 1 10/18/14"    ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 1 3/12/15"     ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 1 5/12/15"     ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 1 6/23/15"     ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 1 8/4/15"      ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 2 10/25/15"    ="BiofilmSalmon2",
                                                                          "BiofilmSalmon 2 10/4/15"     ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 2 3/19/16"     ="BiofilmSalmon2",
                                                                          "BiofilmSalmon 2 5/17/16"     ="BiofilmSalmon2",
                                                                          "BiofilmSalmon 2 6/30/16"     ="BiofilmSalmon2",
                                                                          "BiofilmSalmon 2 8/15/16"     ="BiofilmSalmon2",
                                                                          "HeptageniaSalmon 2 3/19/16"  ="HeptageniaSalmon2",
                                                                          "HeptageniaSalmon 2 5/17/16"  ="HeptageniaSalmon2",
                                                                          "SalmonCarcass 1 10/4/14"     ="Carcass1",
                                                                          "SalmonCarcass 2 10/9/15"     ="Carcass2",
                                                                          "StegopternaControl 1 3/12/15"="StegopternaControl1",
                                                                          "StegopternaControl 1 9/5/14" ="StegopternaControl1",
                                                                          "StegopternaControl 2 5/17/16"="StegopternaControl2",
                                                                          "StegopternaControl 2 6/30/16"="StegopternaControl2",
                                                                          "StegopternaSalmon 1 3/12/15" ="StegopternaSalmon1",
                                                                          "StegopternaSalmon 1 6/23/15" ="StegopternaSalmon1",
                                                                          "StegopternaSalmon 2 3/19/16" ="StegopternaSalmon2",
                                                                          "StegopternaSalmon 2 5/17/16" ="StegopternaSalmon2"))
levels(HC_16S_OTU_map_ag_venn_m$EnvYD)
HC_16S_OTU_map_ag_venn_m_c<- cast(HC_16S_OTU_map_ag_venn_m, variable~EnvYD, sum)
str(HC_16S_OTU_map_ag_venn_m_c)
rownames(HC_16S_OTU_map_ag_venn_m_c)<-HC_16S_OTU_map_ag_venn_m_c[,1]
BaetisControl1<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"BaetisControl1"] > 0]
BaetisControl2<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"BaetisControl2"] > 0]
BaetisSalmon1<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"BaetisSalmon1"] > 0]
BaetisSalmon2<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"BaetisSalmon2"] > 0]
BiofilmControl1<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"BiofilmControl1"] > 0]
BiofilmControl2<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"BiofilmControl2"] > 0]
BiofilmSalmon1<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"BiofilmSalmon1"] > 0]
BiofilmSalmon2<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"BiofilmSalmon2"] > 0]
HeptageniaSalmon2<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"HeptageniaSalmon2"] > 0]
Carcass1<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"Carcass1"] > 0]
Carcass2<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"Carcass2"] > 0]
StegopternaControl1<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"StegopternaControl1"] > 0]
StegopternaControl2<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"StegopternaControl2"] > 0]
StegopternaSalmon1<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"StegopternaSalmon1"] > 0]
StegopternaSalmon2<- rownames(HC_16S_OTU_map_ag_venn_m_c)[HC_16S_OTU_map_ag_venn_m_c[,"StegopternaSalmon2"] > 0]
Control1<-unique(c(BaetisControl1,BiofilmControl1,StegopternaControl1))
Salmon1<-unique(c(BaetisSalmon1,BiofilmSalmon1,StegopternaSalmon1))
Y1<-unique(c(Control1,Salmon1,Carcass1))
Control2<-unique(c(BaetisControl2,BiofilmControl2,StegopternaControl2))
Salmon2<-unique(c(BaetisSalmon2,BiofilmSalmon2,StegopternaSalmon2))
BaetisSalmon<-unique(c(BaetisSalmon1,BaetisSalmon2))
StegopternaSalmon<-unique(c(StegopternaSalmon1,StegopternaSalmon2))
BiofilmSalmon<-unique(c(BiofilmSalmon1,BiofilmSalmon2))

source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
quartz()
Y1SalIntVen<-venn_diagram3(Control1, Salmon1, Carcass1,
              "Y1 Control", "Y1 Treatment", "Y1 Carcass",
              colors=three_col_vec_conttreatcar)
#686 carcass OTUs total
#645 novel to stream
Y1SalIntOTUs<-c(Y1SalIntVen$`Y1 Treatment_Y1 Carcass`,Y1SalIntVen$`Y1 Carcass_only`)
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")
quartz()
#year 1 introduced OTUs with OTUs of both years from different sources
Y1SalIntVenDS<-venn_diagram4(Y1SalIntOTUs,BiofilmSalmon,BaetisSalmon,StegopternaSalmon,
                             "Y1 Car Int","Biofilm","B. brunneicolor","S. mutata",
                             color=four_col_vec_CBfBS)
Y1SalIntRet<-c(Y1SalIntVenDS$`Y1 Car Int_Biofilm`,Y1SalIntVenDS$`Y1 Car Int_B. brunneicolor`,Y1SalIntVenDS$`Y1 Car Int_S. mutata`,Y1SalIntVenDS$`Y1 Car Int_Biofilm_B. brunneicolor`)
Y1SalIntRetBF<-c(Y1SalIntVenDS$`Y1 Car Int_Biofilm`,Y1SalIntVenDS$`Y1 Car Int_Biofilm_B. brunneicolor`)
Y1SalIntRetBae<-c(Y1SalIntVenDS$`Y1 Car Int_B. brunneicolor`,Y1SalIntVenDS$`Y1 Car Int_Biofilm_B. brunneicolor`)
Y1SalIntRetSi<-c(Y1SalIntVenDS$`Y1 Car Int_S. mutata`)
#Biofilm with year 1 introduced OTUs
quartz()
Y1SalIntVenBF<-venn_diagram4(Y1SalIntOTUs,BiofilmSalmon1,BiofilmSalmon2,BiofilmControl2,
                             "Y1 Car Int","Y1 Treat Biof","Y2 Treat Biof","Y2 Cont Biof",
                             color=four_col_vec_CBf1SBf2SBf2C)
#Baetis with year 1 introduced OTUs
quartz()
Y1SalIntVenBae<-venn_diagram4(Y1SalIntOTUs,BaetisSalmon1,BaetisSalmon2,BaetisControl2,
                              "Y1 Car Int","Y1 Treat Baet","Y2 Treat Baet","Y2 Cont Baet",
                              color=four_col_vec_CBa1SBa2SBa2C)
#Stegopterna with year 1 introduced OTUs
quartz()
Y1SalIntVenSi<-venn_diagram4(Y1SalIntOTUs,StegopternaSalmon1,StegopternaSalmon2,StegopternaControl2,
                             "Y1 Car Int","Y1 Treat Steg","Y2 Treat Steg","Y2 Cont Steg",
                             color=four_col_vec_CS1SS2SS2C)
#intvertebrates with year 1 introduced OTUs
quartz()
Y1SalIntVenI<-venn_diagram4(Y1SalIntOTUs,StegopternaSalmon,BaetisSalmon,HeptageniaSalmon2,
                          "Y1Car","S. mutata","B. brunneicolor","H. flavescens")
#Year2 Introdcued OTUs
quartz()
Y2SalIntVen<-venn_diagram4(Y1, Control2, Salmon2, Carcass2,
                           "Y1","Y2 Control", "Y2 Treatment", "Y2 Carcass",
                           colors=four_col_vec_y1conttreatcar)
Y2SalIntOTUs<-c(Y2SalIntVen$`Y2 Treatment_Y2 Carcass`,Y2SalIntVen$`Y2 Carcass_only`)
#Y2 introduced OTUs with different sources
quartz()
Y2SalIntVenDS<-venn_diagram4(Y2SalIntOTUs,BiofilmSalmon2,BaetisSalmon2,StegopternaSalmon2,
                             "Y2 Car Int","Biofilm","B. brunneicolor","S. mutata",
                             color=four_col_vec_CBfBS)
Y2SalIntRet<-c(Y2SalIntVenDS$`Y2 Car Int_Biofilm`,Y2SalIntVenDS$`Y2 Car Int_B. brunneicolor`,Y2SalIntVenDS$`Y2 Car Int_S. mutata`,Y2SalIntVenDS$`Y2 Car Int_Biofilm_S. mutata`)
Y2SalIntRetBF<-c(Y2SalIntVenDS$`Y2 Car Int_Biofilm`,Y2SalIntVenDS$`Y2 Car Int_Biofilm_S. mutata`)
Y2SalIntRetBae<-c(Y2SalIntVenDS$`Y2 Car Int_B. brunneicolor`)
Y2SalIntRetSi<-c(Y2SalIntVenDS$`Y2 Car Int_S. mutata`,Y2SalIntVenDS$`Y2 Car Int_Biofilm_S. mutata`)
quartz()
Y2SalIntVenDSi<-venn_diagram4(Y2SalIntOTUs,HeptageniaSalmon2,BaetisSalmon2,BiofilmSalmon2,
                             "Y2Car","Heptagenia","Baetis","Biofilm")
Y2SalIntRetH<-c(Y2SalIntVenDSi$Y2Car_Heptagenia)
SalIntOTUs<-c(Y1SalIntOTUs,Y2SalIntOTUs)
SalIntRet<-c(Y1SalIntRet,Y2SalIntRet)
#venn diagrams for manuscript

#find out the relative abundance of the unique OTUs in salmon and biofilms
str(HC_16S_OTU_map)
HC_16S_OTU_map_uc1<-HC_16S_OTU_map[,c(1:17, which(names(HC_16S_OTU_map) %in% Y1SalIntOTUs))]
HC_16S_OTU_map_uc1$RowSum<-rowSums(HC_16S_OTU_map_uc1[18:ncol(HC_16S_OTU_map_uc1)])/2500
HC_16S_ucbf1_tot<-subset(HC_16S_OTU_map_uc1, Reach=="Salmon" & (Source=="Carcass" | Source=="Biofilm"), select=c(1:17,ncol(HC_16S_OTU_map_uc1)))
UCBFY1<-summarySE(HC_16S_ucbf1_tot, measurevar="RowSum", groupvars=c("Days_Since_Study_Start", "Source"))
UCBFY1$Year<-rep("One", 14)
HC_16S_OTU_map_uc2<-HC_16S_OTU_map[,c(1:17, which(names(HC_16S_OTU_map) %in% Y2SalIntOTUs))]
HC_16S_OTU_map_uc2$RowSum<-rowSums(HC_16S_OTU_map_uc2[18:ncol(HC_16S_OTU_map_uc2)])/2500
HC_16S_ucbf2_tot<-subset(HC_16S_OTU_map_uc2, Year=="2" & Reach=="Salmon" & (Source=="Carcass" | Source=="Biofilm"), select=c(1:17,ncol(HC_16S_OTU_map_uc2)))
UCBFY2<-summarySE(HC_16S_ucbf2_tot, measurevar="RowSum", groupvars=c("Days_Since_Study_Start", "Source"))
UCBFY2$Year<-rep("Two", 7)
UCBFY<-rbind(UCBFY1,UCBFY2)
#now visualize these unique OTUs
ggplot(UCBFY, aes(x=Days_Since_Study_Start, y=RowSum, shape=Year, color=Source)) + 
  geom_errorbar(aes(ymin=RowSum-se, ymax=RowSum+se), width=.1) +
  geom_point(size=3) +
  geom_line(size=1.5, data=UCBFY[UCBFY$Source!="Carcass", ])+
  xlab("Days since study start") +
  ylab("Mean introduced OTU sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none") +
  scale_color_manual(values=two_col_vec_bfcar) +
  scale_x_continuous(breaks=seq(0,700,100)) +
  scale_y_continuous(trans = "log10") +
  scale_shape_manual(values=c(15, 17))

HC_16S_ucbae1_tot<-subset(HC_16S_OTU_map_uc1, Reach=="Salmon" & (Source=="Carcass" | Source=="Baetis"), select=c(1:17,ncol(HC_16S_OTU_map_uc1)))
UCBaeY1<-summarySE(HC_16S_ucbae1_tot, measurevar="RowSum", groupvars=c("Days_Since_Study_Start", "Source"))
UCBaeY1$Year<-rep("One", 11)
HC_16S_ucbae2_tot<-subset(HC_16S_OTU_map_uc2, Year=="2" & Reach=="Salmon" & (Source=="Carcass" | Source=="Baetis"), select=c(1:17,ncol(HC_16S_OTU_map_uc2)))
UCBaeY2<-summarySE(HC_16S_ucbae2_tot, measurevar="RowSum", groupvars=c("Days_Since_Study_Start", "Source"))
UCBaeY2$Year<-rep("Two", 4)
UCBaeY<-rbind(UCBaeY1,UCBaeY2)
#now visualize these unique OTUs
ggplot(UCBaeY, aes(x=Days_Since_Study_Start, y=RowSum, shape=Year, color=Source)) + 
  geom_errorbar(aes(ymin=RowSum-se, ymax=RowSum+se), width=.1) +
  geom_point(size=3) +
  geom_line(size=1.5, data=UCBaeY[UCBaeY$Source!="Carcass", ])+
  xlab("Days since study start") +
  ylab("Mean introduced OTU sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none") +
  scale_color_manual(values=two_col_vec_baecar) +
  scale_x_continuous(breaks=seq(0,700,100)) +
  scale_y_continuous(trans = "log10") +
  scale_shape_manual(values=c(15, 17))

HC_16S_ucsi1_tot<-subset(HC_16S_OTU_map_uc1, Reach=="Salmon" & (Source=="Carcass" | Source=="Stegopterna"), select=c(1:17,ncol(HC_16S_OTU_map_uc1)))
UCSiY1<-summarySE(HC_16S_ucsi1_tot, measurevar="RowSum", groupvars=c("Days_Since_Study_Start", "Source"))
UCSiY1$Year<-rep("One", 8)
HC_16S_ucsi2_tot<-subset(HC_16S_OTU_map_uc2, Year=="2" & Reach=="Salmon" & (Source=="Carcass" | Source=="Stegopterna"), select=c(1:17,ncol(HC_16S_OTU_map_uc2)))
UCSiY2<-summarySE(HC_16S_ucsi2_tot, measurevar="RowSum", groupvars=c("Days_Since_Study_Start", "Source"))
UCSiY2$Year<-rep("Two", 4)
UCSiY<-rbind(UCSiY1,UCSiY2)
#now visualize these unique OTUs
ggplot(UCSiY, aes(x=Days_Since_Study_Start, y=RowSum, shape=Year, color=Source)) + 
  geom_errorbar(aes(ymin=RowSum-se, ymax=RowSum+se), width=.1) +
  geom_point(size=3) +
  geom_line(size=1.5, data=UCSiY[UCSiY$Source!="Carcass", ])+
  xlab("Days since study start") +
  ylab("Mean introduced OTU sequence abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=18),axis.title.y=element_text(size=18),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none") +
  scale_color_manual(values=two_col_vec_carsi) +
  scale_x_continuous(breaks=seq(0,700,100)) +
  scale_y_continuous(trans = "log10") +
  scale_shape_manual(values=c(15, 17))
UCSiY$facet<-rep("paste('F ',italic('S. mutata'))", 12)
UCBaeY$facet<-rep("paste('E ',italic('B. brunneicolor'))", 15)
UCBFY$facet<-rep("paste('D Biofilm')", 21)

UCY<-rbind(UCSiY,UCBFY,UCBaeY)
UCY$facet<-factor(UCY$facet, levels=c("paste('D Biofilm')","paste('E ',italic('B. brunneicolor'))","paste('F ',italic('S. mutata'))"))
ggplot(UCY, aes(x=Days_Since_Study_Start, y=RowSum, shape=Year, color=Source)) + 
  geom_errorbar(aes(ymin=RowSum-se, ymax=RowSum+se), width=.1) +
  geom_point(size=3) +
  geom_line(size=1.5, data=UCY[UCY$Source!="Carcass", ])+
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")+
  xlab("Days Since Study Start") +
  ylab("Mean Introduced OTU Relative Sequence Abundance (+/- SEM)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(family="Helvetica",size=14),
        axis.title.y=element_text(family="Helvetica",size=14),
        axis.text.x=element_text(family="Helvetica",size=10),
        axis.text.y = element_text(family="Helvetica",size=10, face="bold"),
        strip.text.x = element_text(family="Helvetica",size = 12),
        legend.text.align = 0,
        legend.text=element_text(family="Helvetica",size=10),
        legend.title=element_text(family="Helvetica",size=12))+
  scale_color_manual(values=four_col_vec_CBfBS, 
                     limits=c("Carcass","Biofilm","Baetis","Stegopterna"),
                     labels=c("Carcass","Biofilm",expression(paste(italic('B. brunneicolor'))),expression(paste(italic('S. mutata'))))) +
  scale_x_continuous(breaks=seq(0,700,100)) +
  scale_y_continuous(trans = "log10", limits =c(0.1,1700), labels=c(1,10,100,1000), 
                     breaks=c(1,10,100,1000)) +
  scale_shape_manual(values=c(15, 17))+
  facet_wrap(~facet,ncol=1, labeller=label_parsed)
UCY$facet2<-UCY$facet
levels(UCY$facet2)[levels(UCY$facet2)=="paste('D Biofilm')"] <- "Biofilm"
levels(UCY$facet2)[levels(UCY$facet2)=="paste('E ',italic('B. brunneicolor'))"] <- "italic('B. brunneicolor')"
levels(UCY$facet2)[levels(UCY$facet2)=="paste('F ',italic('S. mutata'))"] <- "italic('S. mutata')"
ggplot(UCY, aes(x=Days_Since_Study_Start, y=RowSum, shape=Year, color=Source)) + 
  geom_errorbar(aes(ymin=RowSum-se, ymax=RowSum+se), width=.1) +
  geom_point(size=3) +
  geom_line(size=1.5, data=UCY[UCY$Source!="Carcass", ])+
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")+
  xlab("Days Since Study Start") +
  ylab("Mean Introduced ASV Rel. Sequence Abund. (± SE)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(family="Helvetica",size=20),
        axis.title.y=element_text(family="Helvetica",size=18),
        axis.text.x=element_text(family="Helvetica",size=16),
        axis.text.y = element_text(family="Helvetica",size=16, face="bold"),
        strip.text.x = element_text(family="Helvetica",size = 18),
        legend.text.align = 0,
        legend.text=element_text(family="Helvetica",size=16),
        legend.title=element_text(family="Helvetica",size=18))+
  scale_color_manual(values=four_col_vec_CBfBS, 
                     limits=c("Carcass","Biofilm","Baetis","Stegopterna"),
                     labels=c("Carcass","Biofilm",expression(paste(italic('B. brunneicolor'))),expression(paste(italic('S. mutata'))))) +
  scale_x_continuous(breaks=seq(0,700,200)) +
  scale_y_continuous(trans = "log10", limits =c(0.001,1), labels=c(0.001,0.01,0.1,1), 
                     breaks=c(0.001,0.01,0.1,1)) +
  scale_shape_manual(values=c(15, 17))+
  facet_wrap(~facet2,ncol=3, labeller=label_parsed)

#Create faceted heatmap for all three sources
HC_16S_OTU_map_uc<-HC_16S_OTU_map[,c(1:17, which(names(HC_16S_OTU_map) %in% SalIntRet))]
HC_16S_OTU_map_uc_hm<-subset(HC_16S_OTU_map_uc, Reach=="Salmon" & Source!="Heptagenia" & Source!="Carcass", select=c(10,13,18:ncol(HC_16S_OTU_map_uc)))
HC_16S_OTU_map_uc_hm<-aggregate(HC_16S_OTU_map_uc_hm[,3:ncol(HC_16S_OTU_map_uc_hm)], by=list(Source=HC_16S_OTU_map_uc_hm$Source,Days_Since_Study_Start=HC_16S_OTU_map_uc_hm$Days_Since_Study_Start),FUN=sum)
HC_16S_OTU_map_uc_hm_m<-melt(HC_16S_OTU_map_uc_hm, id=c("Source","Days_Since_Study_Start")) 
HC_16S_OTU_map_uc_hm_m$Source<-factor(HC_16S_OTU_map_uc_hm_m$Source, levels=c("Biofilm","Baetis","Stegopterna"))
HC_16S_OTU_map_uc_hm_m$Source<-revalue(HC_16S_OTU_map_uc_hm_m$Source, c("Baetis"=expression(paste("B. ",italic("B. brunneicolor"))), "Stegopterna"=expression(paste("C. ",italic("S. mutata"))), "Biofilm"=expression(paste("A. Biofilm"))))
HC_16S_OTU_map_uc_hm_m$Days_Since_Study_Start<-factor(HC_16S_OTU_map_uc_hm_m$Days_Since_Study_Start)
HC_16S_OTU_map_uc_hm_m_s<-subset(HC_16S_OTU_map_uc_hm_m, (Days_Since_Study_Start=="29" & Source=='paste(\"A. Biofilm\")' & variable %in% BiofilmSalmon) | value>0)
HC_16S_OTU_map_uc_hm_m_s$value[HC_16S_OTU_map_uc_hm_m_s$value == 0] <- NA
HC_Taxonomy<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_taxonomy.txt", sep="\t", header=T)
UTax<-HC_Taxonomy[HC_Taxonomy$FeatureID %in% HC_16S_OTU_map_uc_hm_m_s$variable,]
UTax<-UTax[match(HC_16S_OTU_map_uc_hm_m_s$variable, UTax$FeatureID),]
UTax$Taxon<-revalue(UTax$Taxon, 
                    c("k__Bacteria"                                                                                                                        ="Unknown Bacteria",
                      "k__Bacteria; p__Acidobacteria; c__[Chloracidobacteria]; o__DS-100; f__; g__; s__"                                                   ="Unknown DS-100",
                      "k__Bacteria; p__Acidobacteria; c__Acidobacteria-6; o__iii1-15; f__; g__; s__"                                                       ="Unknown iii1-15",
                      "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Corynebacteriaceae; g__Corynebacterium; s__kroppenstedtii"="Corynbacterium kroppenstedtii",
                      "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Kineosporiaceae; g__; s__"                                ="Unknown Kineosporiaceae",
                      "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Micrococcaceae; g__Arthrobacter"                          ="Arthrobacter sp.",
                      "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Mycobacteriaceae; g__Mycobacterium; s__vaccae"            ="Mycobacterium vaccae",
                      "k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Pseudonocardiaceae; g__Pseudonocardia; s__"               ="Pseudonocardia spp.",
                      "k__Bacteria; p__Actinobacteria; c__Thermoleophilia; o__Gaiellales; f__Gaiellaceae; g__; s__"                                        ="Unknown Gaiellaceae",
                      "k__Bacteria; p__Actinobacteria; c__Thermoleophilia; o__Solirubrobacterales; f__; g__; s__"                                          ="Unknown Solirurobacterales",
                      "k__Bacteria; p__Bacteroidetes; c__[Saprospirae]; o__[Saprospirales]; f__Chitinophagaceae; g__; s__"                                 ="Unknown Chitinophagaceae",
                      "k__Bacteria; p__Bacteroidetes; c__[Saprospirae]; o__[Saprospirales]; f__Chitinophagaceae; g__Niastella; s__"                        ="Niastella spp.",
                      "k__Bacteria; p__Bacteroidetes; c__[Saprospirae]; o__[Saprospirales]; f__Chitinophagaceae; g__Sediminibacterium; s__"                ="Sediminibacterium spp.",
                      "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__; g__; s__"                                                     ="Unknown Bacteroidales",
                      "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__ML635J-40; g__; s__"                                            ="Unknown ML635J-40",
                      "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Prevotellaceae; g__Prevotella; s__copri"                        ="Prevotella copri",
                      "k__Bacteria; p__Bacteroidetes; c__Cytophagia; o__Cytophagales; f__Cytophagaceae; g__Dyadobacter; s__"                               ="Dyadobacter spp.",
                      "k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__[Weeksellaceae]; g__Chryseobacterium; s__"                ="Chryseobacterium spp.",
                      "k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__[Weeksellaceae]; g__Wautersiella; s__"                    ="Wautersiella spp.",
                      "k__Bacteria; p__Bacteroidetes; c__Sphingobacteriia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Pedobacter; s__"              ="Pedobacter spp.",
                      "k__Bacteria; p__Chloroflexi; c__Ellin6529; o__; f__; g__; s__"                                                                      ="Unknown Ellin6529",
                      "k__Bacteria; p__Chloroflexi; c__TK17; o__mle1-48; f__; g__; s__"                                                                    ="Unknown mle1-48",
                      "k__Bacteria; p__Cyanobacteria; c__Chloroplast; o__Chlorophyta; f__; g__; s__"                                                       ="Unknown Chlorophyta",
                      "k__Bacteria; p__Cyanobacteria; c__Chloroplast; o__Streptophyta; f__; g__; s__"                                                      ="Unknown Streptophyta",
                      "k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus; s__"                                            ="Bacillus spp.",
                      "k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Paenibacillaceae; g__Aneurinibacillus; s__"                               ="Aneurinibacillus spp.",
                      "k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Paenibacillaceae; g__Cohnella; s__"                                       ="Cohnella spp.",
                      "k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Paenibacillaceae; g__Paenibacillus; s__"                                  ="Paenibacillus spp.",
                      "k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Planococcaceae"                                                           ="Unknown Planococcaceae",
                      "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__"                             ="Lactobacillus spp.",
                      "k__Bacteria; p__Firmicutes; c__Bacilli; o__Turicibacterales; f__Turicibacteraceae; g__Turicibacter; s__"                            ="Turicibacter spp.",
                      "k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Ruminococcaceae; g__Oscillospira; s__"                              ="Oscillospira spp.",
                      "k__Bacteria; p__Firmicutes; c__Clostridia; o__MBA08; f__; g__; s__"                                                                 ="Unknown MBA08",
                      "k__Bacteria; p__Nitrospirae; c__Nitrospira; o__Nitrospirales; f__Nitrospiraceae; g__Nitrospira; s__"                                ="Nitrospira spp.",
                      "k__Bacteria; p__Planctomycetes; c__Phycisphaerae; o__Phycisphaerales; f__Phycisphaeraceae; g__; s__"                                ="Unknown Phycisphaeraceae",
                      "k__Bacteria; p__Planctomycetes; c__Planctomycetia; o__Gemmatales; f__Gemmataceae; g__Gemmata; s__"                                  ="Gemmata spp.",
                      "k__Bacteria; p__Planctomycetes; c__Planctomycetia; o__Pirellulales; f__Pirellulaceae; g__; s__"                                     ="Unknown Pirellulaceae",
                      "k__Bacteria; p__Planctomycetes; c__Planctomycetia; o__Planctomycetales; f__Planctomycetaceae; g__Planctomyces; s__"                 ="Planctomyces spp.",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Hyphomicrobiaceae; g__; s__"                             ="Unknown Hyphomicrobiaceae",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Hyphomicrobiaceae; g__Hyphomicrobium; s__"               ="Hyphomicrobium spp.",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Hyphomicrobiaceae; g__Rhodoplanes; s__"                  ="Rhodoplanes spp.",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Phyllobacteriaceae; g__Mesorhizobium; s__"               ="Mesorhizobium spp.",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiaceae; g__; s__"                                  ="Unknown Rhizobiaceae",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae; g__Amaricoccus; s__"               ="Amaricoccus spp.",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Acetobacteraceae; g__; s__"                         ="Unknown Acetobacteraceae",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Rhodospirillaceae; g__Reyranella; s__massiliensis"  ="Reyranella massilienses",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__; g__; s__"                                         ="Unknown Sphingomonadales",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae; g__Kaistobacter; s__"            ="Kaistobacter spp.",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae; g__Novosphingobium"              ="Novoshingobium sp.",
                      "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae; g__Sphingomonas; s__"            ="Sphingomonas spp.",
                      "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__; f__; g__; s__"                                                          ="Unknown Betaproteobacteria",
                      "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Alcaligenaceae; g__; s__"                             ="Unknown Alcaligenacaea",
                      "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Variovorax; s__paradoxus"          ="Variovorax paradoxus",
                      "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Oxalobacteraceae; g__Janthinobacterium; s__"          ="Janthinobacterium spp.",
                      "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Neisseriales; f__Neisseriaceae; g__Iodobacter; s__fluviatilis"            ="Iodobacter fluviatilis",
                      "k__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Desulfuromonadales; f__Geobacteraceae; g__Geobacter; s__"                ="Geobacter spp.",
                      "k__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Myxococcales; f__; g__; s__"                                             ="Unknown Myxococcales",
                      "k__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Syntrophobacterales; f__Syntrophobacteraceae; g__; s__"                  ="Unknown Syntrophobacteraceae",
                      "k__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Syntrophobacterales; f__Syntrophorhabdaceae; g__; s__"                   ="Unknown Syntrophorhabdaceae",
                      "k__Bacteria; p__Proteobacteria; c__Epsilonproteobacteria; o__Campylobacterales; f__Helicobacteraceae; g__Flexispira; s__"           ="Flexispira spp.",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Alteromonadales; f__Shewanellaceae; g__Shewanella; s__"                  ="Shewanella spp.",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae"                                ="Unknown Enterobacteriaceae",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Serratia; s__"              ="Serratia spp.",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter"                     ="Acinetobacter sp.",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter; s__"                ="Acinetobacter spp.", 
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter; s__guillouiae"      ="Acinetobacter guillouiae",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Psychrobacter; s__marincola"       ="Psychrobacter marincola",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Psychrobacter; s__pulmonis"        ="Psychrobacter pulmonis",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__Pseudomonas; s__"               ="Pseudomonas spp.",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Thiotrichales; f__Piscirickettsiaceae; g__; s__"                         ="Unknown Piscirickettsiaceae",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Sinobacteraceae; g__; s__"                           ="Unknown Sinobacteraceae",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Sinobacteraceae; g__Steroidobacter; s__"             ="Steroidobacter spp.",
                      "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae; g__Luteimonas; s__"                ="Luteimonas spp.",
                      "k__Bacteria; p__Spirochaetes; c__Spirochaetes; o__Spirochaetales; f__Spirochaetaceae; g__Treponema; s__"                            ="Treponema spp.",
                      "k__Bacteria; p__Verrucomicrobia; c__[Pedosphaerae]; o__[Pedosphaerales]; f__Ellin517; g__; s__"                                     ="Unknown Ellin517"))
HC_16S_OTU_map_uc_hm_m_s$variable<-UTax$Taxon
str(HC_16S_OTU_map_uc_hm_m_s)
ggplot(HC_16S_OTU_map_uc_hm_m_s, aes(x=Days_Since_Study_Start, y=variable, fill=value)) +
  geom_tile(color="white", size=0.1)+
  geom_vline(xintercept = c(1.5,7.5), linetype = "dotted", colour = "black")+
  facet_wrap(.~Source, ncol=1, scale="free_y",labeller=label_parsed)+
  labs(x="Days Since Study Start", y="Introduced Novel OTUs")+
  scale_fill_gradient(name = "Relative \nSequence \nAbundance", trans = "log", 
                      breaks=c(1,10,100), low="cornflowerblue",high="black", na.value = "white")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=10),axis.text.y = element_blank(),
        strip.text.x = element_text(size = 12),axis.ticks.y = element_blank())

#make table for manuscript
HC_16S_OTU_map_salint<-HC_16S_OTU_map[,c(1:17, which(names(HC_16S_OTU_map) %in% SalIntRet))]
HC_16S_OTU_map_salint_man<-subset(HC_16S_OTU_map_salint, Reach=="Salmon" & Env!="BaetisControl" & Env!="BiofilmControl" & Env!="StegopternaControl" & Source!="Heptagenia", select=c(7,10,14,16,18:ncol(HC_16S_OTU_map_salint)))
HC_16S_OTU_map_salint_man$EnvYD<-as.factor(paste(HC_16S_OTU_map_salint_man$Env, HC_16S_OTU_map_salint_man$Year, HC_16S_OTU_map_salint_man$Date))
levels(HC_16S_OTU_map_salint_man$EnvYD)
HC_16S_OTU_map_salint_man$EnvYD<-revalue(HC_16S_OTU_map_salint_man$EnvYD, c("BaetisSalmon 1 3/12/15"      ="BaetisSalmon1",
                                                                          "BaetisSalmon 1 5/12/15"      ="BaetisSalmon1",
                                                                          "BaetisSalmon 1 6/23/15"      ="BaetisSalmon1",
                                                                          "BaetisSalmon 1 8/4/15"       ="BaetisSalmon1",
                                                                          "BaetisSalmon 1 9/20/15"      ="BaetisSalmon1",
                                                                          "BaetisSalmon 2 3/19/16"      ="BaetisSalmon2",
                                                                          "BaetisSalmon 2 5/17/16"      ="BaetisSalmon2",
                                                                          "BaetisSalmon 2 6/30/16"      ="BaetisSalmon2",
                                                                          "BiofilmSalmon 1 10/18/14"    ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 1 3/12/15"     ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 1 5/12/15"     ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 1 6/23/15"     ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 1 8/4/15"      ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 2 10/25/15"    ="BiofilmSalmon2",
                                                                          "BiofilmSalmon 2 10/4/15"     ="BiofilmSalmon1",
                                                                          "BiofilmSalmon 2 3/19/16"     ="BiofilmSalmon2",
                                                                          "BiofilmSalmon 2 5/17/16"     ="BiofilmSalmon2",
                                                                          "BiofilmSalmon 2 6/30/16"     ="BiofilmSalmon2",
                                                                          "BiofilmSalmon 2 8/15/16"     ="BiofilmSalmon2",
                                                                          "SalmonCarcass 1 10/4/14"     ="Carcass1",
                                                                          "SalmonCarcass 2 10/9/15"     ="Carcass2",
                                                                          "StegopternaSalmon 1 3/12/15" ="StegopternaSalmon1",
                                                                          "StegopternaSalmon 1 6/23/15" ="StegopternaSalmon1",
                                                                          "StegopternaSalmon 2 3/19/16" ="StegopternaSalmon2",
                                                                          "StegopternaSalmon 2 5/17/16" ="StegopternaSalmon2"))
HC_16S_OTU_map_salint_man_ag<-aggregate(HC_16S_OTU_map_salint_man[,5:(ncol(HC_16S_OTU_map_salint_man)-1)], by=list(EnvYD=HC_16S_OTU_map_salint_man$EnvYD),FUN=sum)
HC_16S_OTU_map_salint_man_ag_m<-melt(HC_16S_OTU_map_salint_man_ag, id=c("EnvYD")) 
UTaxman<-HC_Taxonomy[HC_Taxonomy$FeatureID %in% HC_16S_OTU_map_salint_man_ag_m$variable,]
UTaxman<-UTaxman[match(HC_16S_OTU_map_salint_man_ag_m$variable, UTaxman$FeatureID),]
HC_16S_OTU_map_salint_man_ag_m$variable<-UTaxman$Taxon
HC_16S_OTU_map_salint_man_ag_m_c<-cast(HC_16S_OTU_map_salint_man_ag_m, variable~EnvYD, sum)
write.table(HC_16S_OTU_map_salint_man_ag_m_c, "~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/IntOTUs.txt", sep="\t") 

###############################
#Analyze picrust output KEG pathways
##############################
Hunt_Creek_KEG<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/HC_cat_func.txt", sep="\t", header = T)
str(Hunt_Creek_KEG)
#Format data frame so the KEG Code is row name
row.names(Hunt_Creek_KEG)<-Hunt_Creek_KEG[,1]

#Delete otu id column, now that otu id is row name
Hunt_Creek_KEG$ID<-NULL
str(Hunt_Creek_KEG)
#transpose
Hunt_Creek_KEG_t<-t(Hunt_Creek_KEG)
str(Hunt_Creek_KEG_t)
Hunt_Creek_KEG_t<-data.frame(Hunt_Creek_KEG_t, check.names = FALSE)
str(Hunt_Creek_KEG_t)
#Merge metadata onto data table
Hunt_Creek_KEG_map_r<-Hunt_Creek_16S_map
row.names(Hunt_Creek_KEG_map_r)<-Hunt_Creek_KEG_map_r[,1]
HC_16S_KEG_map <-merge(Hunt_Creek_KEG_map_r, Hunt_Creek_KEG_t, by=0)
str(HC_16S_KEG_map)
HC_16S_KEG_map<-subset(HC_16S_KEG_map, Source!="Rhyacophila")
HC_16S_KEG_map$Row.names<-NULL

#Create community matrix
HC_16S_KEG_com<-HC_16S_KEG_map[,17:ncol(HC_16S_KEG_map)]

#Create overall environmental data matrix for community analysis with counts
HC_16S_KEG_env<-HC_16S_KEG_map
row.names(HC_16S_KEG_env)<-HC_16S_KEG_env[,1]
HC_16S_KEG_env<-HC_16S_KEG_env[,-c(1)]

#Overall permanova
adonis(HC_16S_KEG_com ~ Reach*Year*Source*Total_Biofilm_Growth_PostCarcass, data=HC_16S_KEG_env_i1, method="jaccard", permutations=999)
#Significant factors are Source, Total_Biofilm_Growth_PostCarcass, Year:Source, and Source:Total_Biofilm_Growth_PostCarcass

#Create carcass keg table with metadata
HC_16S_KEG_map_car<-subset(HC_16S_KEG_map, Source=="Carcass")
#Subset into environmental and community data frames
HC_16S_KEG_car_env<-HC_16S_KEG_map_car[,1:17]
HC_16S_KEG_car_com<-HC_16S_KEG_map_car[,18:ncol(HC_16S_KEG_map_car)]
HC_16S_KEG_car_Year<-factor(HC_16S_KEG_car_env$Year)
levels(HC_16S_KEG_car_Year)<-c("Year 1","Year 2")
#carcass permanova with KEG codes
adonis(HC_16S_KEG_car_com ~ Year, data=HC_16S_KEG_car_env, method="jaccard", permutations=999)
#Year is significant
#NMDS carcass
HC_NMDS_KEG_c<-metaMDS(HC_16S_KEG_car_com, distance="jaccard")
HC_NMDS_KEG_c
#Stress=0.06
stressplot(HC_NMDS_KEG_c)
ordiplot(HC_NMDS_KEG_c, type="n", xlim=c(-0.3,0.3), ylim=c(-0.2,0.2))
with(HC_NMDS_KEG_c, points(HC_NMDS_KEG_c, display="sites", col = "black", pch=c(15,17)[HC_16S_KEG_car_Year]))
with(HC_NMDS_KEG_c, legend("topleft", legend=levels(HC_16S_KEG_car_Year), bty="n", col="black",pch=c(15,17)))
with(HC_NMDS_KEG_c, ordiellipse(HC_NMDS_KEG_c, HC_16S_KEG_car_Year, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Year 1"))
with(HC_NMDS_KEG_c, ordiellipse(HC_NMDS_KEG_c, HC_16S_KEG_car_Year, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Year 2"))
with(HC_NMDS_KEG_c, legend("topright", legend="2D Stress: 0.06", bty="n"))
#Indicator species analysis
HC_KEG_c_indic<-signassoc(HC_16S_KEG_car_com, cluster=HC_16S_KEG_car_Year,  mode=0, alternative = "two.sided",control = how(nperm=9999))
HC_KEG_c_indic_sig<-subset(HC_KEG_c_indic, psidak<=0.05)
#135 indicator KEG orthologs
write.table(HC_KEG_c_indic_sig, "~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/KEGCarIndSig.txt", sep="\t") 
head(sort(HC_KEG_c_indic_sig$psidak))
HC_16S_KEG_car_agg<-aggregate(HC_16S_KEG_map_car[18:ncol(HC_16S_KEG_map_car)], by=list(Year=HC_16S_KEG_map_car$Year), FUN=sum)
HC_16S_KEG_car_agg_sig<-subset(HC_16S_KEG_car_agg, select=c(rownames(HC_KEG_c_indic_sig)))
tail(sort(colSums(HC_16S_KEG_car_agg_sig)))
HC_16S_KEG_car_Y1<-subset(HC_16S_KEG_map_car, Year==1)
HC_16S_KEG_car_Y2<-subset(HC_16S_KEG_map_car, Year==2)
stat.desc(HC_16S_KEG_car_Y1$DNArepairandrecombinationproteins)
stat.desc(HC_16S_KEG_car_Y2$DNArepairandrecombinationproteins)
stat.desc(HC_16S_KEG_car_Y1$Melanogenesis)
stat.desc(HC_16S_KEG_car_Y2$Melanogenesis)

#Create year 1 keg table with metadata
HC_16S_KEG_map_Y1<-subset(HC_16S_KEG_map, Year=="1")

#Create year2 keg table with metadata
HC_16S_KEG_map_Y2<-subset(HC_16S_KEG_map, Year=="2")

#without metadata year 1
H_C_16S_KEG_Y1<-HC_16S_KEG_map_Y1[,18:ncol(HC_16S_KEG_map_Y1)]
H_C_16S_KEG_Y1_samples<-as.vector(rownames(H_C_16S_KEG_Y1))
#Year 1 environmental variable table
H_C_16S_KEG_env_Y1<-HC_16S_KEG_map_Y1[,1:17]
Total_Biofilm_Growth_PostCarcass_Y1_KEG<-as.factor(H_C_16S_KEG_env_Y1$Total_Biofilm_Growth_PostCarcass)
Reach_Y1_KEG<-factor(H_C_16S_KEG_env_Y1$Reach)
Env_Y1_KEG<-factor(H_C_16S_KEG_env_Y1$Env)
levels(Env_Y1_KEG)<-c("Control","Salmon","Control","Salmon","Salmon")
Source_Y1_KEG<-factor(H_C_16S_KEG_env_Y1$Source)
levels(Source_Y1_KEG)<-c("B. brunneicolor","Biofilm","Carcass","S. mutata")
Source_Y1_KEG<- ordered(Source_Y1_KEG, levels = c("Biofilm", "Carcass", "B. brunneicolor","S. mutata"))
#without metadata year 2
H_C_16S_KEG_Y2<-HC_16S_KEG_map_Y2[,18:ncol(HC_16S_KEG_map_Y2)]
H_C_16S_KEG_Y2_samples<-as.vector(rownames(H_C_16S_KEG_Y2))
#Year 2 environmental variable table
H_C_16S_KEG_env_Y2<-HC_16S_KEG_map_Y2[,1:17]
Total_Biofilm_Growth_PostCarcass_Y2_KEG<-factor(H_C_16S_KEG_env_Y2$Total_Biofilm_Growth_PostCarcass)
Reach_Y2_KEG<-factor(H_C_16S_KEG_env_Y2$Reach)
Env_Y2_KEG<-factor(H_C_16S_KEG_env_Y2$Env)
levels(Env_Y2_KEG)<-c("Control","Salmon","Control","Salmon","Salmon")
Source_Y2_KEG<-factor(H_C_16S_KEG_env_Y2$Source)
levels(Source_Y2_KEG)<-c("B. brunneicolor","Biofilm","Carcass","H. flavescens","S. mutata","S. mutata")
Source_Y2_KEG<- ordered(Source_Y2_KEG, levels = c("Biofilm", "Carcass", "B. brunneicolor","H. flavescens","S. mutata"))
#Year 1 permanova using counts
adonis(H_C_16S_KEG_Y1 ~ Source*Reach_Y1_KEG*Total_Biofilm_Growth_PostCarcass_Y1_KEG, data=H_C_16S_KEG_env_Y1, method="jaccard", permutations=999)
#Source,time,Source:time and ~reach:time interaction all significant

#Year 2 permanova using counts
adonis(H_C_16S_KEG_Y2 ~ Source*Reach_Y2_KEG*Total_Biofilm_Growth_PostCarcass_Y2_KEG, data=H_C_16S_KEG_env_Y2, method="jaccard", permutations=999)
#Source,time and Source:time interaction significant

#NMDS year 1
HC_NMDS_KEG_y1<-metaMDS(H_C_16S_KEG_Y1, distance="jaccard")
HC_NMDS_KEG_y1
stressplot(HC_NMDS_KEG_y1)
#stress 0.17
ordiplot(HC_NMDS_KEG_y1, type="n")
with(HC_NMDS_KEG_y1, points(HC_NMDS_KEG_y1, display="sites", col =c(four_col_vec_BfCBS)[Source_Y1_KEG], pch=c(15,17)[Env_Y1_KEG]))
with(HC_NMDS_KEG_y1, legend("topleft", legend=c(levels(Source_Y1_KEG),levels(Env_Y1_KEG)), bty="n", col=c(four_col_vec_BfCBS,"black","black"),pch=c(16,16,16,16,0,2)))
with(HC_NMDS_KEG_y1, legend("topright", legend="2D Stress: 0.17", bty="n"))

#NMDS year 2
HC_NMDS_KEG_y2<-metaMDS(H_C_16S_KEG_Y2, distance="jaccard")
HC_NMDS_KEG_y2
stressplot(HC_NMDS_KEG_y2)
#stress 0.17
ordiplot(HC_NMDS_KEG_y2, type="n")
with(HC_NMDS_KEG_y2, points(HC_NMDS_KEG_y2, display="sites", col =c(five_col_vec_babchs)[Source_Y2_KEG], pch=c(15,17)[Env_Y2_KEG]))
with(HC_NMDS_KEG_y2, legend("topleft", legend=c(levels(Source_Y2_KEG),levels(Env_Y2_KEG)), bty="n", col=c(five_col_vec_babchs,"black","black"),pch=c(16,16,16,16,16,0,2)))
with(HC_NMDS_KEG_y2, legend("topright", legend="2D Stress: 0.18", bty="n"))

#Create biofilm matrix for year 1 with metadata
H_C_16S_KEG_map_Y1_BF<-subset(HC_16S_KEG_map_Y1, Source=="Biofilm")
#without metadata for biofilms in year 1
H_C_16S_KEG_Y1_BF<-H_C_16S_KEG_map_Y1_BF[,18:ncol(H_C_16S_KEG_map_Y1_BF)]
H_C_16S_KEG_Y1_BF_samples<-as.vector(rownames(H_C_16S_KEG_Y1_BF))
#Biofilm year 1 environmental variable table
H_C_16S_KEG_env_Y1_BF<-H_C_16S_KEG_map_Y1_BF[,1:17]
Total_Biofilm_Growth_PostCarcass_Y1_KEG<-(H_C_16S_KEG_env_Y1_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y1_KEG<-factor(H_C_16S_KEG_env_Y1_BF$Reach)
Env_Y1_KEG<-factor(revalue(H_C_16S_KEG_env_Y1_BF$Env, c("BiofilmSalmon"="Salmon", "BiofilmControl"="Control")), levels =c("Salmon", "Control"))

#Create biofilm matrix for year 2 with metadata
H_C_16S_KEG_map_Y2_BF<-subset(HC_16S_KEG_map_Y2, Source=="Biofilm")
#without metadata for biofilms in year 2
H_C_16S_KEG_Y2_BF<-H_C_16S_KEG_map_Y2_BF[,18:ncol(H_C_16S_KEG_map_Y2_BF)]
H_C_16S_KEG_Y2_BF_samples<-as.vector(rownames(H_C_16S_KEG_Y2_BF))
#Biofilm year 2 environmental variable table
H_C_16S_KEG_env_Y2_BF<-H_C_16S_KEG_map_Y2_BF[,1:17]
Total_Biofilm_Growth_PostCarcass_Y2_KEG<-(H_C_16S_KEG_env_Y2_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y2_KEG<-factor(H_C_16S_KEG_env_Y2_BF$Reach)
Env_Y2_KEG<-factor(revalue(H_C_16S_KEG_env_Y2_BF$Env, c("BiofilmSalmon"="Salmon", "BiofilmControl"="Control")), levels =c("Salmon", "Control"))

#Year 1 biofilm permanova using counts
adonis(H_C_16S_KEG_Y1_BF ~ Reach_Y1_KEG*Total_Biofilm_Growth_PostCarcass_Y1_KEG, data=H_C_16S_KEG_env_Y1_BF, method="jaccard", permutations=9999)
#Reach and time significant

#Year 2 biofilm permanova using counts
adonis(H_C_16S_KEG_Y2_BF ~ Reach_Y2_KEG*Total_Biofilm_Growth_PostCarcass_Y2_KEG, data=H_C_16S_KEG_env_Y2_BF, method="jaccard", permutations=9999)
#nothing significant

#NMDS biofilms year 1
HC_NMDS_KEG_bfy1<-metaMDS(H_C_16S_KEG_Y1_BF, distance="jaccard")
HC_NMDS_KEG_bfy1
stressplot(HC_NMDS_KEG_bfy1)
#stress 0.06
ordiplot(HC_NMDS_KEG_bfy1, type="n", xlim=c(-0.5,0.4), ylim=c(-0.2,0.4))
with(HC_NMDS_KEG_bfy1, points(HC_NMDS_KEG_bfy1, display="sites", col = "black", pch=c(16,5)[Env_Y1_KEG]))
with(HC_NMDS_KEG_bfy1, legend("topleft", legend=levels(Env_Y1_KEG), bty="n", col="black",pch=c(16,5)))
with(HC_NMDS_KEG_bfy1, ordiellipse(HC_NMDS_KEG_bfy1, Env_Y1_KEG, kind="se", conf=0.95, lwd=2, col="black", draw="polygon", show.groups = "Salmon"))
with(HC_NMDS_KEG_bfy1, ordiellipse(HC_NMDS_KEG_bfy1, Env_Y1_KEG, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Control"))
with(HC_NMDS_KEG_bfy1, legend("topright", legend="2D Stress: 0.06", bty="n"))

#NMDS biofilms year 2
HC_NMDS_KEG_bfy2<-metaMDS(H_C_16S_KEG_Y2_BF, distance="jaccard")
HC_NMDS_KEG_bfy2
stressplot(HC_NMDS_KEG_bfy2)
#stress 0.14
ordiplot(HC_NMDS_KEG_bfy2, type="n", xlim=c(-0.4,0.3), ylim=c(-0.2,0.4))
with(HC_NMDS_KEG_bfy2, points(HC_NMDS_KEG_bfy2, display="sites", col = "black", pch=c(16,5)[Env_Y2_KEG]))
with(HC_NMDS_KEG_bfy2, legend("topleft", legend=levels(Env_Y2_KEG), bty="n", col="black",pch=c(16,5)))
with(HC_NMDS_KEG_bfy2, ordiellipse(HC_NMDS_KEG_bfy2, Env_Y2_KEG, kind="se", conf=0.95, lwd=2, col="black", draw="polygon", show.groups = "Salmon"))
with(HC_NMDS_KEG_bfy2, ordiellipse(HC_NMDS_KEG_bfy2, Env_Y2_KEG, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Control"))
with(HC_NMDS_KEG_bfy2, legend("topright", legend="2D Stress: 0.14", bty="n"))

#Indicator species analysis
#year 1
HC_KEG_bfy1_indic<-signassoc(H_C_16S_KEG_Y1_BF, cluster=Env_Y1_KEG,  mode=0, alternative = "two.sided",control = how(nperm=9999))
HC_KEG_bfy1_indic_sig<-subset(HC_KEG_bfy1_indic, psidak<=0.05)
#113 indicator KEG orthologs
write.table(HC_KEG_bfy1_indic_sig, "~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/KEGBFY1IndSig.txt", sep="\t") 
head(sort(HC_KEG_bfy1_indic_sig$psidak))
#Fluorobenzoatedegradation most significant
H_C_16S_KEG_Y1_BF_agg<-aggregate(H_C_16S_KEG_map_Y1_BF[18:ncol(H_C_16S_KEG_map_Y1_BF)], by=list(Env=H_C_16S_KEG_map_Y1_BF$Env), FUN=sum)
H_C_16S_KEG_Y1_BF_agg_sig<-subset(H_C_16S_KEG_Y1_BF_agg, select=c(rownames(HC_KEG_bfy1_indic_sig)))
sort(colSums(H_C_16S_KEG_Y1_BF_agg_sig))
#most abundant Twocomponentsystem
HC_16S_KEG_map_BF_sal<-subset(HC_16S_KEG_map, Env=="BiofilmSalmon")
HC_16S_KEG_map_BF_cont<-subset(HC_16S_KEG_map, Env=="BiofilmControl")
stat.desc(HC_16S_KEG_map_BF_sal$Fluorobenzoatedegradation)
stat.desc(HC_16S_KEG_map_BF_cont$Fluorobenzoatedegradation)
stat.desc(HC_16S_KEG_map_BF_sal$Twocomponentsystem)
stat.desc(HC_16S_KEG_map_BF_cont$Twocomponentsystem)
stat.desc(HC_16S_KEG_map_BF_sal$PhosphotransferasesystemPTS)
stat.desc(HC_16S_KEG_map_BF_cont$PhosphotransferasesystemPTS)
HC_16S_KEG_map_BF_Y1_sal<-subset(H_C_16S_KEG_map_Y1_BF, Env=="BiofilmSalmon")
HC_16S_KEG_map_BF_Y1_cont<-subset(H_C_16S_KEG_map_Y1_BF, Env=="BiofilmControl")
stat.desc(HC_16S_KEG_map_BF_Y1_sal$Melanogenesis)
stat.desc(HC_16S_KEG_map_BF_Y1_cont$Melanogenesis)

#Create baetis matrix for year 1 with metadata
H_C_16S_KEG_map_Y1_Bae<-subset(HC_16S_KEG_map_Y1, Source=="Baetis")
#without metadata for baetids in year 1
H_C_16S_KEG_Y1_Bae<-H_C_16S_KEG_map_Y1_Bae[,18:ncol(H_C_16S_KEG_map_Y1_Bae)]
H_C_16S_KEG_Y1_Bae_samples<-as.vector(rownames(H_C_16S_KEG_Y1_Bae))
#Baetis year 1 environmental variable table
H_C_16S_KEG_env_Y1_Bae<-H_C_16S_KEG_map_Y1_Bae[,1:17]
Total_Biofilm_Growth_PostCarcass_Y1_KEG_Bae<-(H_C_16S_KEG_env_Y1_Bae$Total_Biofilm_Growth_PostCarcass)
Reach_Y1_KEG_Bae<-factor(H_C_16S_KEG_env_Y1_Bae$Reach)
Env_Y1_KEG_Bae<-factor(revalue(H_C_16S_KEG_env_Y1_Bae$Env, c("InsectSalmon"="Salmon", "InsectControl"="Control")), levels =c("Salmon", "Control"))

#Create baetis matrix for year 2 with metadata
H_C_16S_KEG_map_Y2_Bae<-subset(HC_16S_KEG_map_Y2, Source=="Baetis")
#without metadata for baetids in year 2
H_C_16S_KEG_Y2_Bae<-H_C_16S_KEG_map_Y2_Bae[,18:ncol(H_C_16S_KEG_map_Y2_Bae)]
H_C_16S_KEG_Y2_Bae_samples<-as.vector(rownames(H_C_16S_KEG_Y2_Bae))
#Baetis year 2 environmental variable table
H_C_16S_KEG_env_Y2_Bae<-H_C_16S_KEG_map_Y2_Bae[,1:17]
Total_Biofilm_Growth_PostCarcass_Y2_KEG_Bae<-(H_C_16S_KEG_env_Y2_Bae$Total_Biofilm_Growth_PostCarcass)
Reach_Y2_KEG_Bae<-factor(H_C_16S_KEG_env_Y2_Bae$Reach)
Env_Y2_KEG_Bae<-factor(revalue(H_C_16S_KEG_env_Y2_Bae$Env, c("InsectSalmon"="Salmon", "InsectControl"="Control")), levels =c("Salmon", "Control"))

#Year 1 baetis permanova using counts
adonis(H_C_16S_KEG_Y1_Bae ~ Reach_Y1_KEG_Bae*Total_Biofilm_Growth_PostCarcass_Y1_KEG_Bae, data=H_C_16S_KEG_env_Y1_Bae, method="jaccard", permutations=9999)
#Nothing significant

#Year 2 baetis permanova using counts
adonis(H_C_16S_KEG_Y2_Bae ~ Reach_Y2_KEG_Bae*Total_Biofilm_Growth_PostCarcass_Y2_KEG_Bae, data=H_C_16S_KEG_env_Y2_Bae, method="jaccard", permutations=9999)
#Nothing significant

#Create simulids matrix for year 1 with metadata
H_C_16S_KEG_map_Y1_S<-subset(HC_16S_KEG_map_Y1, Source=="Simulium" | Source=="Stegopterna")
#without metadata for simulids in year 1
H_C_16S_KEG_Y1_S<-H_C_16S_KEG_map_Y1_S[,18:ncol(H_C_16S_KEG_map_Y1_S)]
H_C_16S_KEG_Y1_S_samples<-as.vector(rownames(H_C_16S_KEG_Y1_S))
#Simulids year 1 environmental variable table
H_C_16S_KEG_env_Y1_S<-H_C_16S_KEG_map_Y1_S[,1:17]
Reach_Y1_KEG_S<-factor(H_C_16S_KEG_env_Y1_S$Reach)
Env_Y1_KEG_S<-factor(revalue(H_C_16S_KEG_env_Y1_S$Env, c("InsectSalmon"="Salmon", "InsectControl"="Control")), levels =c("Salmon", "Control"))

#Create simulids matrix for year 2 with metadata
H_C_16S_KEG_map_Y2_S<-subset(HC_16S_KEG_map_Y2, Source=="Simulium" | Source=="Stegopterna")
#without metadata for simulids in year 2
H_C_16S_KEG_Y2_S<-H_C_16S_KEG_map_Y2_S[,18:ncol(H_C_16S_KEG_map_Y2_S)]
H_C_16S_KEG_Y2_S_samples<-as.vector(rownames(H_C_16S_KEG_Y2_S))
#Simulids year 2 environmental variable table
H_C_16S_KEG_env_Y2_S<-H_C_16S_KEG_map_Y2_S[,1:17]
Reach_Y2_KEG_S<-factor(H_C_16S_KEG_env_Y2_S$Reach)
Env_Y2_KEG_S<-factor(revalue(H_C_16S_KEG_env_Y2_S$Env, c("InsectSalmon"="Salmon", "InsectControl"="Control")), levels =c("Salmon", "Control"))

#Year 1 simulids permanova using counts
adonis(H_C_16S_KEG_Y1_S ~ Reach*Total_Biofilm_Growth_PostCarcass, data=H_C_16S_KEG_env_Y1_S, method="jaccard", permutations=9999)
#Nothing significant

#Year 2 simulids permanova using counts
adonis(H_C_16S_KEG_Y2_S ~ Reach*Total_Biofilm_Growth_PostCarcass, data=H_C_16S_KEG_env_Y2_S, method="jaccard", permutations=9999)
#Only reach:time interaction significant

#Indicator analysis of Simuliids year 2 to see what KEGG orthologs are driving this change
H_C_16S_KEG_Y2_S_indic<-signassoc(H_C_16S_KEG_Y2_S, cluster=Env_Y2_KEG_S,  mode=0, alternative = "two.sided",control = how(nperm=9999))
H_C_16S_KEG_Y2_indic_sig<-subset(H_C_16S_KEG_Y2_S_indic, psidak<=0.05)
H_C_16S_KEG_Y2_S_agg<-aggregate(H_C_16S_KEG_map_Y2_S[18:ncol(H_C_16S_KEG_map_Y2_S)], by=list(Env=H_C_16S_KEG_map_Y2_S$Env), FUN=sum)
H_C_16S_KEG_Y2_S_agg$Melanogenesis
H_C_16S_KEG_Y2_S_sal<-subset(H_C_16S_KEG_map_Y2_S, Env=="StegopternaSalmon")
H_C_16S_KEG_Y2_S_cont<-subset(H_C_16S_KEG_map_Y2_S, Env=="StegopternaControl")
stat.desc(H_C_16S_KEG_Y2_S_sal$Melanogenesis)
stat.desc(H_C_16S_KEG_Y2_S_cont$Melanogenesis)

#make melanogenesis plot
mel<-summarySE(HC_16S_KEG_map,measurevar="Melanogenesis", groupvars=c("Reach", "Source","Days_Since_Study_Start"))
mel<-subset(mel, Source!="Heptagenia")
mel$Facet<-mel$Source
mel$Facet<-revalue(mel$Facet, c("Carcass"=expression(paste("A. Carcass")),"Biofilm"=expression(paste("B. Biofilm")),"Baetis"=expression(paste("C. ",italic("B. brunneicolor"))),"Stegopterna"=expression(paste("D. ",italic("S. mutata")))))
mel$Facet<-factor(mel$Facet, levels=c(expression(paste("A. Carcass")),expression(paste("B. Biofilm")), expression(paste("C. ",italic("B. brunneicolor"))),expression(paste("D. ",italic("S. mutata")))))
mel$Source<-factor(mel$Source, levels=c("Carcass","Biofilm","Baetis","Stegopterna"))
mellines<- data.frame(z=c(levels(mel$Facet),levels(mel$Facet)), vl=c(NA,30,30,30,NA,400,400,400))
ggplot(mel, aes(y=Melanogenesis, x=Days_Since_Study_Start, color=Source, linetype=Reach)) + 
  geom_pointrange(aes(ymin=Melanogenesis-se, ymax=Melanogenesis+se),size=1.25, linetype=1) +
  geom_line(size=3, data=mel[mel$Source!="Carcass", ])+
  scale_color_manual(values=four_col_vec_CBfBS,
                     limits=c("Carcass","Biofilm","Baetis","Stegopterna"),
                     labels=c("Carcass","Biofilm",expression(paste(italic('B. brunneicolor'))),
                              expression(paste(italic('S. mutata'))))) +
  facet_wrap(.~Facet, scales="free_y", labeller=label_parsed)+
  xlab("Days Since Study Start") +
  ylab("Mean Relative Melanogenesis Abundance (+/- SEM)") +
  scale_x_continuous(breaks=seq(0,700,200)) + 
  scale_linetype_manual(values=c("dashed", "solid"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(family="Times New Roman", size=32),
        axis.title.y=element_text(family="Times New Roman", size=32),
        axis.text.x=element_text(family="Times New Roman", size=24),
        axis.text.y = element_text(family="Times New Roman", size=24),
        strip.text.x = element_text(family="Times New Roman", face="bold", size = 24), legend.position = "none") +
  geom_vline(data=filter(mel, Facet=="paste(\"B. Biofilm\")"), aes(xintercept=c(30)), linetype="dotted", 
             colour="black")+
  geom_vline(data=filter(mel, Facet=="paste(\"B. Biofilm\")"), aes(xintercept=c(400)), linetype="dotted", 
             colour="black")+
  geom_vline(data=filter(mel, Facet=="paste(\"C. \", italic(\"B. brunneicolor\"))"), aes(xintercept=c(30)), 
             linetype="dotted", colour="black")+
  geom_vline(data=filter(mel, Facet=="paste(\"C. \", italic(\"B. brunneicolor\"))"), aes(xintercept=c(400)), 
             linetype="dotted", colour="black")+
  geom_vline(data=filter(mel, Facet=="paste(\"D. \", italic(\"S. mutata\"))"), aes(xintercept=c(30)), 
             linetype="dotted", colour="black")+
  geom_vline(data=filter(mel, Facet=="paste(\"D. \", italic(\"S. mutata\"))"), aes(xintercept=c(400)), 
             linetype="dotted", colour="black")

##############
#Diversity as response 
################
#upload diversity dataset
Hunt_Creek_FPD<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/faith-alpha-diversity.txt", sep="\t", header = T)
row.names(Hunt_Creek_FPD)<-Hunt_Creek_FPD[,1]
Hunt_Creek_FPD_map_r<-Hunt_Creek_16S_map
row.names(Hunt_Creek_FPD_map_r)<-Hunt_Creek_FPD_map_r[,1]
HC_16S_FPD_map <-merge(Hunt_Creek_FPD_map_r, Hunt_Creek_FPD, by=0)
str(HC_16S_FPD_map)
HC_16S_FPD_map$Row.names<-NULL

#subset for carcass
HC_16S_FPD_map_car<-subset(HC_16S_FPD_map, Source=="Carcass")
#determine if diversity differs based on year
t.test(faith_pd~Year, data=HC_16S_FPD_map_car)
#t = -0.49307, df = 12.719, p-value = 0.6304

#subset for biofilms
HC_16S_FPD_map_bf_y1<-subset(HC_16S_FPD_map, Source=="Biofilm" & Year=="1")
HC_16S_FPD_map_bf_y2<-subset(HC_16S_FPD_map, Source=="Biofilm" & Year=="2")
#determine if diversity differs based on treatment
shapiro.test(HC_16S_FPD_map_bf_y1$faith_pd)
#not normally distributed W = 0.9339, p-value = 0.03291
hist(HC_16S_FPD_map_bf_y1$faith_pd)
#reflect
max(HC_16S_FPD_map_bf_y1$faith_pd)
#43.6
HC_16S_FPD_map_bf_y1$faith_pdr<-44-HC_16S_FPD_map_bf_y1$faith_pd
#transform
HC_16S_FPD_map_bf_y1$srfpdr<-sqrt(HC_16S_FPD_map_bf_y1$faith_pdr)
shapiro.test(HC_16S_FPD_map_bf_y1$srfpdr)
#normally distributed W = 0.95359, p-value = 0.1357
#create model year 1
BFDiv_y1.lm = lm(srfpdr ~ Reach*Days_Since_Study_Start, data=HC_16S_FPD_map_bf_y1) 
#Check model assumptions
plot(BFDiv_y1.lm)
#looks okay
summary(BFDiv_y1.lm)
#only time significant decrease Days_Since_Study_Start Estimate: -0.0099848  Std. Error: 0.0016893  t value: -5.911 Pr(>t):1.41e-06 ***

#create model year 2
shapiro.test(HC_16S_FPD_map_bf_y2$faith_pd)
#not normal W = 0.90556, p-value = 0.00486
hist(HC_16S_FPD_map_bf_y2$faith_pd)
HC_16S_FPD_map_bf_y2$srfpd<-sqrt(HC_16S_FPD_map_bf_y2$faith_pd)
shapiro.test(HC_16S_FPD_map_bf_y2$srfpd)
#normal W = 0.95129, p-value = 0.1148
BFDiv_y2.lm = lm(srfpd ~ Reach*Days_Since_Study_Start, data=HC_16S_FPD_map_bf_y2) 
#Check model assumptions
plot(BFDiv_y2.lm)
#looks okay
summary(BFDiv_y2.lm)
#only time significant decrease Days_Since_Study_Start Estimate:-0.005975   Std. Error:0.001930  t value:-3.096  Pr(>t):0.00406 ** 

#subset for baetids
HC_16S_FPD_map_bae_y1<-subset(HC_16S_FPD_map, Source=="Baetis" & Year=="1")
HC_16S_FPD_map_bae_y2<-subset(HC_16S_FPD_map, Source=="Baetis" & Year=="2")
#determine if diversity differs based on treatment
#check normality
shapiro.test(HC_16S_FPD_map_bae_y1$faith_pd)
#not normal W = 0.80957, p-value = 0.002753
#square root transform
HC_16S_FPD_map_bae_y1$srfpd<-sqrt(HC_16S_FPD_map_bae_y1$faith_pd)
shapiro.test(HC_16S_FPD_map_bae_y1$srfpd)
#not normal W = 0.85983, p-value = 0.0152
#log transform
HC_16S_FPD_map_bae_y1$log10fpd<-log10(HC_16S_FPD_map_bae_y1$faith_pd)
shapiro.test(HC_16S_FPD_map_bae_y1$log10fpd)
#normal W = 0.90478, p-value = 0.0817
BaeDiv_y1.lm.log10 = lm(log10fpd ~ Reach*Days_Since_Study_Start, data=HC_16S_FPD_map_bae_y1) 
plot(BaeDiv_y1.lm.log10)
hist(HC_16S_FPD_map_bae_y1$log10fpd)
#looks okay
summary(BaeDiv_y1.lm.log10)
#only time significant decrease Days_Since_Study_Start -0.005763   0.001799  -3.203  0.00693 ** 

#create model year 2
#check normality
shapiro.test(HC_16S_FPD_map_bae_y2$faith_pd)
#normally distributed W = 0.97351, p-value = 0.639
BaeDiv_y2.lm = lm(faith_pd ~ Reach*Days_Since_Study_Start, data=HC_16S_FPD_map_bae_y2) 
#Check model assumptions
plot(BaeDiv_y2.lm)
#looks okay
summary(BaeDiv_y2.lm)
#only time significant decrease
#log transform for consistency between years
HC_16S_FPD_map_bae_y2$log10fpd<-log10(HC_16S_FPD_map_bae_y2$faith_pd)
shapiro.test(HC_16S_FPD_map_bae_y2$log10fpd)
#still normally distributed W = 0.95546, p-value = 0.2362
BaeDiv_y2.lm.log10 = lm(log10fpd ~ Reach*Days_Since_Study_Start, data=HC_16S_FPD_map_bae_y2) 
#Check model assumptions
plot(BaeDiv_y2.lm.log10)
#looks okay
summary(BaeDiv_y2.lm.log10)
#reach significantly decreased in salmon, time significantly decreased, significant interaction

#subset for simulids
HC_16S_FPD_map_si_y1<-subset(HC_16S_FPD_map, (Source=="Stegopterna") & Year=="1")
HC_16S_FPD_map_si_y2<-subset(HC_16S_FPD_map, (Source=="Stegopterna") & Year=="2")
#determine if diversity differs based on treatment
shapiro.test(HC_16S_FPD_map_si_y1$faith_pd)
#normally distributed W = 0.961, p-value = 0.8273
#create model year 1
SiDiv_y1.lm = lm(faith_pd ~ Reach*Days_Since_Study_Start, data=HC_16S_FPD_map_si_y1) 
#Check model assumptions
hist(HC_16S_FPD_map_si_y1$faith_pd)
plot(SiDiv_y1.lm)
#112 outlier, but don't want to delete since such small dataset
summary(SiDiv_y1.lm)
#nothing significant

#create model year 2
shapiro.test(HC_16S_FPD_map_si_y2$faith_pd)
#normally distributed W = 0.93751, p-value = 0.2899
hist(HC_16S_FPD_map_si_y2$faith_pd)
#some possible bimodal distribution
SiDiv_y2.lm = lm(faith_pd ~ Reach*Days_Since_Study_Start, data=HC_16S_FPD_map_si_y2) 
#Check model assumptions
plot(SiDiv_y2.lm)
#58 outlier
HC_16S_FPD_map_si_y2.no<-HC_16S_FPD_map_si_y2[-6,]
SiDiv_y2.lm.no = lm(faith_pd ~ Reach*Days_Since_Study_Start, data=HC_16S_FPD_map_si_y2.no) 
#Check model assumptions
plot(SiDiv_y2.lm.no)
summary(SiDiv_y2.lm.no)
#reach significant

#upload chao1 richness dataset
Hunt_Creek_ch<-read.table("~/Documents/MSU/Research/Hunt_Creek_Salmon/Microbes/16S/chao1-alpha-diversity.txt", sep="\t", header = T)
row.names(Hunt_Creek_ch)<-Hunt_Creek_ch[,1]
Hunt_Creek_ch_map_r<-Hunt_Creek_16S_map
row.names(Hunt_Creek_ch_map_r)<-Hunt_Creek_ch_map_r[,1]
HC_16S_ch_map <-merge(Hunt_Creek_ch_map_r, Hunt_Creek_ch, by=0)
str(HC_16S_ch_map)
HC_16S_ch_map$Row.names<-NULL

#subset for carcass
HC_16S_ch_map_car<-subset(HC_16S_ch_map, Source=="Carcass")
#determine if richness differs based on year
t.test(chao1~Year, data=HC_16S_ch_map_car)
#t = -1.7171, df = 10.462, p-value = 0.1154

#subset for biofilms
HC_16S_ch_map_bf_y1<-subset(HC_16S_ch_map, Source=="Biofilm" & Year=="1")
HC_16S_ch_map_bf_y2<-subset(HC_16S_ch_map, Source=="Biofilm" & Year=="2")
#determine if richness differs based on treatment
shapiro.test(HC_16S_ch_map_bf_y1$chao1)
#not normal W = 0.90584, p-value = 0.004946
hist(HC_16S_ch_map_bf_y1$chao1)
boxcox(chao1 ~ Days_Since_Study_Start, data=HC_16S_ch_map_bf_y1)
HC_16S_ch_map_bf_y1$log10chao1<-log10(HC_16S_ch_map_bf_y1$chao1)
shapiro.test(HC_16S_ch_map_bf_y1$log10chao1)
#not normal, W = 0.84329, p-value = 0.0001339 use original scale
#create model year 1
BFch_y1.lm = lm(log10chao1 ~ Days_Since_Study_Start, data=HC_16S_ch_map_bf_y1) 
#Check model assumptions
plot(BFch_y1.lm)
#looks okay
summary(BFch_y1.lm)
#only time significant increase Days_Since_Study_Start estimate 0.0029878  std error 0.0002377   12.57 2.46e-14 ***


#create model year 2
BFch_y2.lm = lm(chao1 ~ Reach*Days_Since_Study_Start, data=HC_16S_ch_map_bf_y2) 
boxcox(BFch_y2.lm)
hist(HC_16S_ch_map_bf_y2$chao1)
HC_16S_ch_map_bf_y2$log10chao1<-log10(HC_16S_ch_map_bf_y2$chao1)
shapiro.test(HC_16S_ch_map_bf_y2$log10chao1)
hist(HC_16S_ch_map_bf_y2$log10chao1)
BFch_y2.lm.log10 = lm(log10chao1 ~ Days_Since_Study_Start, data=HC_16S_ch_map_bf_y2) 
#Check model assumptions
plot(BFch_y2.lm.log10)
#looks okay
summary(BFch_y2.lm.log10)
#only time significant decrease estimate: Days_Since_Study_Start estimate -0.0021682  standard error 0.0003562  -6.087 6.65e-07 ***

#subset for baetids
HC_16S_ch_map_bae_y1<-subset(HC_16S_ch_map, Source=="Baetis" & Year=="1")
HC_16S_ch_map_bae_y2<-subset(HC_16S_ch_map, Source=="Baetis" & Year=="2")
#determine if richness differs based on treatment
#create model year 1
Baech_y1.lm = lm(chao1 ~ Reach*Days_Since_Study_Start, data=HC_16S_ch_map_bae_y1) 
#Check model assumptions
plot(Baech_y1.lm)
#121, 114 and 47 outliers
HC_16S_ch_map_bae_y1_no<-HC_16S_ch_map_bae_y1[-c(3,9,11),]
#create model without outliers
Baech_y1.lm.no = lm(chao1 ~ Reach*Days_Since_Study_Start, data=HC_16S_ch_map_bae_y1_no) 
#Check model assumptions
plot(Baech_y1.lm.no)
#looks good
summary(Baech_y1.lm.no)
#nothing significant

#create model year 2
Baech_y2.lm = lm(chao1 ~ Reach*Days_Since_Study_Start, data=HC_16S_ch_map_bae_y2) 
#Check model assumptions
plot(Baech_y2.lm)
#15 outlier
HC_16S_ch_map_bae_y2_no<-HC_16S_ch_map_bae_y2[-3,]
#create model year 2 no outlier
Baech_y2.lmno = lm(chao1 ~ Reach*Days_Since_Study_Start, data=HC_16S_ch_map_bae_y2_no) 
#Check model assumptions
plot(Baech_y2.lmno)
#looks good
summary(Baech_y2.lmno)
#nothing significant

#subset for simulids
HC_16S_ch_map_si_y1<-subset(HC_16S_ch_map, (Source=="Stegopterna" | Source=="Simulium") & Year=="1")
HC_16S_ch_map_si_y2<-subset(HC_16S_ch_map, (Source=="Stegopterna" | Source=="Simulium") & Year=="2")
#determine if richness differs based on treatment
#create model year 1
Sich_y1.lm = lm(chao1 ~ Reach*Days_Since_Study_Start, data=HC_16S_ch_map_si_y1) 
#Check model assumptions
plot(Sich_y1.lm)
#normality assumptions not met, possible outlier in 113
HC_16S_ch_map_si_y1_no<-HC_16S_ch_map_si_y1[-4,]
#create model year 1 without outlier
Sich_y1.lm.no = lm(chao1 ~ Reach*Days_Since_Study_Start, data=HC_16S_ch_map_si_y1_no) 
#Check model assumptions
plot(Sich_y1.lm.no)
#looks okay
summary(Sich_y1.lm.no)
#nothing significant

#create model year 2
Sich_y2.lm = lm(chao1 ~ Reach*Days_Since_Study_Start, data=HC_16S_ch_map_si_y2) 
#Check model assumptions
plot(Sich_y2.lm)
#looks good
summary(Sich_y2.lm)
#nothing significant

