#Analysis of 16S microbial communities at Hunt Creek 2014-2016

#########################################
#Load packages
library(reshape)
library(ggplot2)
library(vegan)
library(plyr)
library(dplyr)

#Functions
#Functions
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
three_col_vec<- c("#a6cee3", "#1f78b4", "#b2df8a")
five_col_vec<- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")
six_col_vec<- c("#e41a1c", "#377eb8", "green", "#984ea3", "#ff7f00", "#ffff33")
six_col_vec_cont<-c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15")

#####################################
#Upload data tables generated in QIIME and manipulate to make "r friendly"
#####################################

#Get 16S OTU table
Hunt_Creek_16S<-read.table("~/Desktop/HC_table_tabseparated.txt", sep="\t", header = T)

#Classify as matrices
Hunt_Creek_16S<-data.frame(Hunt_Creek_16S)

#Format data frame so the denovo is row name
row.names(Hunt_Creek_16S)<-Hunt_Creek_16S[,1]

#Delete taxonomy and denovo columns, now that denovo is row name
Hunt_Creek_16S$SampleID<-NULL

#transpose
HC_16S_OTU_t<-t(Hunt_Creek_16S)
HC_16S_OTU_t<-data.frame(HC_16S_OTU_t)

#Get metadata through mapping file
Hunt_Creek_16S_map <- read.table("~/Desktop/Hunt_Creek_Map.txt", header=T)
row.names(Hunt_Creek_16S_map)<-Hunt_Creek_16S_map[,1]
Hunt_Creek_16S_map<-Hunt_Creek_16S_map[,-c(1)]

#Merge metadata onto data table
HC_16S_OTU_map <-merge(Hunt_Creek_16S_map, HC_16S_OTU_t, by=0)

#Upload weighted unifrac distance matrix and remove samples with error
Hunt_Creek_16S_uni <- read.table("~/Desktop/weighted_unifrac_otu_table_HC.txt", header=T)
Hunt_Creek_16S_uni<-Hunt_Creek_16S_uni[1:94,1:94]

#Add metadata to unifrac table
Hunt_Creek_16S_uni_map <-merge(Hunt_Creek_16S_map, Hunt_Creek_16S_uni, by=0)

#create overal community data matrix for community analysis
HC_16S_com<-data.frame(HC_16S_OTU_map[,15:ncol(HC_16S_OTU_map)])

#Create overall environmental data matrix for community analysis
HC_16S_uni_env<-(Hunt_Creek_16S_uni_map[,1:14])
row.names(HC_16S_uni_env)<-HC_16S_uni_env[,1]
HC_16S_uni_env<-HC_16S_uni_env[,-c(1)]
#merge insects into one source type
HC_16S_uni_env_i1<-HC_16S_uni_env
HC_16S_uni_env_i1$Source<-gsub("Rhyacophila", "Insect", HC_16S_uni_env_i1$Source)
HC_16S_uni_env_i1$Source<-gsub("Stegopterna", "Insect", HC_16S_uni_env_i1$Source)
HC_16S_uni_env_i1$Source<-as.factor(gsub("Baetis", "Insect", HC_16S_uni_env_i1$Source))
levels(HC_16S_uni_env_i1$Source)
HC_16S_uni_env_i1$Total_Biofilm_Growth_PostCarcass<-as.factor(HC_16S_uni_env_i1$Total_Biofilm_Growth_PostCarcass)
HC_16S_uni_env_i1$Year<-as.factor(HC_16S_uni_env_i1$Year)
HC_16S_uni_env_i1$Reach<-as.factor(HC_16S_uni_env_i1$Reach)
#############################################
#Analysis of all 16S community data
#############################################

#Overall permanova
adonis(as.dist(Hunt_Creek_16S_uni) ~ Reach*Year*Source*Total_Biofilm_Growth_PostCarcass, data=HC_16S_uni_env_i1, permutations=999)

#############################################
#Separate data based on year of study (1 or 2)
###############################################

#Create Year 1 unifrac distance table with metadata
H_C_16S_uni_map_Y1<-subset(Hunt_Creek_16S_uni_map, Year=="1")
rownames(H_C_16S_uni_map_Y1)<-H_C_16S_uni_map_Y1[,1]
H_C_16S_uni_map_Y1<-H_C_16S_uni_map_Y1[,-c(1)]

#Create Year 2 unifrac distance table with metadata
H_C_16S_uni_map_Y2<-subset(Hunt_Creek_16S_uni_map, Year=="2")
rownames(H_C_16S_uni_map_Y2)<-H_C_16S_uni_map_Y2[,1]
H_C_16S_uni_map_Y2<-H_C_16S_uni_map_Y2[,-c(1)]

#Create year 1 otu table with metadata
HC_16S_OTU_map_Y1<-subset(HC_16S_OTU_map, Year=="1")

#Create year2 otu table with metadata
HC_16S_OTU_map_Y2<-subset(HC_16S_OTU_map, Year=="2")

####################################################
#Biofilm sample analysis by year
############################################

#Create biofilm data tables by year 

#Biofilm Year 1
#Create biofilm unifrac distance table for year 1 with metadata
H_C_16S_uni_map_Y1_BF<-subset(H_C_16S_uni_map_Y1, Source=="Biofilm")
#unifrac distance table without metadata for biofilms in year 1
H_C_16S_uni_Y1_BF<-H_C_16S_uni_map_Y1_BF[,14:ncol(H_C_16S_uni_map_Y1_BF)]
H_C_16S_uni_Y1_BF_samples<-as.vector(rownames(H_C_16S_uni_Y1_BF))
#use output of names to subset columns into rows to make square matrix
H_C_16S_uni_Y1_BF<-as.matrix(H_C_16S_uni_Y1_BF)
H_C_16S_uni_Y1_BF<-subset(H_C_16S_uni_Y1_BF, select=c(H_C_16S_uni_Y1_BF_samples))
#Biofilm year 1 environmental variable table
H_C_16S_env_Y1_BF<-H_C_16S_uni_map_Y1_BF[,1:13]
Total_Biofilm_Growth_PostCarcass_Y1<-as.factor(H_C_16S_env_Y1_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y1<-as.factor(H_C_16S_env_Y1_BF$Reach)

#Biofilm year 2
#Create biofilm unifrac distance table for year 2 with metadata
H_C_16S_uni_map_Y2_BF<-subset(H_C_16S_uni_map_Y2, Source=="Biofilm")
#unifrac distance table without metadata for biofilms in year 2
H_C_16S_uni_Y2_BF<-H_C_16S_uni_map_Y2_BF[,14:ncol(H_C_16S_uni_map_Y2_BF)]
H_C_16S_uni_Y2_BF_samples<-as.vector(rownames(H_C_16S_uni_Y2_BF))
#use output of names to subset columns into rows to make square matrix
H_C_16S_uni_Y2_BF<-as.matrix(H_C_16S_uni_Y2_BF)
H_C_16S_uni_Y2_BF<-subset(H_C_16S_uni_Y2_BF, select=c(H_C_16S_uni_Y2_BF_samples))
#Biofilm year 2 environmental variable table
H_C_16S_env_Y2_BF<-H_C_16S_uni_map_Y2_BF[,1:13]
Total_Biofilm_Growth_PostCarcass_Y2<-as.factor(H_C_16S_env_Y2_BF$Total_Biofilm_Growth_PostCarcass)
Reach_Y2<-as.factor(H_C_16S_env_Y2_BF$Reach)

#Year 1 biofilm permanova
adonis(as.dist(H_C_16S_uni_Y1_BF) ~ Reach_Y1*Total_Biofilm_Growth_PostCarcass_Y1, data=H_C_16S_env_Y1_BF, permutations=999)

#Year 2 biofilm permanova
adonis(as.dist(H_C_16S_uni_Y2_BF) ~ Reach_Y2*Total_Biofilm_Growth_PostCarcass_Y2, data=H_C_16S_env_Y2_BF, permutations=999)

###########################################################
#Carcass and biofilm comparisons
##########################################################

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
HC_16S_OTU_map_Y1_CBF_sub<-data.frame(t(HC_16S_OTU_map_Y1_CBF[,15:ncol(HC_16S_OTU_map_Y1_CBF)]))
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
HC_16S_OTU_map_Y2_CBF_sub$HC.58<-as.numeric(as.character(HC_16S_OTU_map_Y2_CBF_sub$HC.58))
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
#Create data table with only upstream biofilms for year 1 then find unique
HC_16S_OTU_map_Y2_BF_US<-subset(HC_16S_OTU_map_Y2, Reach=="Control" & Source=="Biofilm")
rownames(HC_16S_OTU_map_Y2_BF_US)<-HC_16S_OTU_map_Y2_BF_US[,1]
HC_16S_OTU_map_Y2_BF_US<-HC_16S_OTU_map_Y2_BF_US[,-c(1)]
HC_16S_OTU_map_Y2_BF_US_unique<-data.frame(t(HC_16S_OTU_map_Y2_BF_US[,15:ncol(HC_16S_OTU_map_Y2_BF_US)]))
HC_16S_OTU_map_Y2_BF_US_unique$OTU<-row.names(HC_16S_OTU_map_Y2_BF_US_unique)
HC_16S_OTU_map_Y2_BF_US_unique<-HC_16S_OTU_map_Y2_BF_US_unique[HC_16S_OTU_map_Y2_BF_US_unique$OTU %in% CBF_Unique_Y2,]
HC_16S_OTU_map_Y2_BF_US_unique$HC.61<-as.numeric(as.character(HC_16S_OTU_map_Y2_BF_US_unique$HC.61))
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

#Find out if any of the 20 really unique OTUs in year 1 persist in year 2 biofilms, both US and DS
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


#Make nmds plot for carcass and biofilms for year 1
HC_16S_uni_map_Y1_CBF<-subset(H_C_16S_uni_map_Y1, Reach=="Salmon")
HC_16S_uni_map_Y1_CBF<-subset(HC_16S_uni_map_Y1_CBF, Date=="10/3/14" | Date=="10/4/14" | Date=="10/18/14")
HC_16S_uni_map_Y1_CBF_com<-HC_16S_uni_map_Y1_CBF[,c(14,15,16,42,51,56,57,58,59,60,63,64,84,85,88)]
HC_16S_uni_map_Y1_CBF_env<-HC_16S_uni_map_Y1_CBF[,1:13]
Y1_CBF<-as.factor(HC_16S_uni_map_Y1_CBF_env$Total_Biofilm_Growth)
levels(Y1_CBF)
HC_CBF_NMDS_Y1<-metaMDS(as.dist(HC_16S_uni_map_Y1_CBF_com))
ordiplot(HC_CBF_NMDS_Y1, type="n", main="Biofilms Before and After and Carcass")
with(HC_CBF_NMDS_Y1, points(HC_CBF_NMDS_Y1, display="sites", col=three_col_vec[Y1_CBF], pch=19, pt.bg=three_col_vec))
with(HC_CBF_NMDS_Y1, legend("topleft", legend=levels(Y1_CBF), bty="n", col=three_col_vec, pch=19, pt.bg=three_col_vec))

#nmds plot for carcass and biofilms for year 2
HC_16S_uni_map_Y2_CBF<-subset(H_C_16S_uni_map_Y2, Reach=="Salmon")
HC_16S_uni_map_Y2_CBF<-subset(HC_16S_uni_map_Y2_CBF, Date=="10/4/15" | Date=="10/9/15" | Date=="10/25/15")
HC_16S_uni_map_Y2_CBF_com<-HC_16S_uni_map_Y2_CBF[,c(19,30,31,36,37,40,61,62,65,87,89,91,95,97,100,101)]
HC_16S_uni_map_Y2_CBF_env<-HC_16S_uni_map_Y2_CBF[,1:13]
Y2_CBF<-as.factor(HC_16S_uni_map_Y2_CBF_env$Total_Biofilm_Growth)
levels(Y1_CBF)
HC_CBF_NMDS_Y2<-metaMDS(as.dist(HC_16S_uni_map_Y2_CBF_com))
ordiplot(HC_CBF_NMDS_Y2, type="n",main="Biofilms Before and After and Carcass")
with(HC_CBF_NMDS_Y2, points(HC_CBF_NMDS_Y2, display="sites", col=three_col_vec[Y2_CBF], pch=19, pt.bg=three_col_vec))
with(HC_CBF_NMDS_Y2, legend("topleft", legend=levels(Y2_CBF), bty="n", col=three_col_vec, pch=19, pt.bg=three_col_vec))

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

#Upload phyla level files for each run
Hunt_Creek_16S_P<-read.table("~/Desktop/HC_Phyla_16S.txt", sep="\t", header = T)

#Clasify as data.frame
Hunt_Creek_16S_P<-data.frame(Hunt_Creek_16S_P)

#Format data frame so the taxonomy is row name
row.names(Hunt_Creek_16S_P)<-Hunt_Creek_16S_P[,1]
#Delete taxonomy column
Hunt_Creek_16S_P$SampleID<-NULL

#transpose
HC_16S_P_t<-t(Hunt_Creek_16S_P)
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

#Split into year 1 and year 2 data
HC_16S_P_map_Y1<-subset(HC_16S_P_map, Year== "1")
HC_16S_P_map_Y2<-subset(HC_16S_P_map, Year== "2")

#limit phyla data to biofilms for each year
HC_16S_P_map_Y1_BF<-subset(HC_16S_P_map_Y1, Source == "Biofilm")
HC_16S_P_map_Y2_BF<-subset(HC_16S_P_map_Y2, Source == "Biofilm")

#Make line plot for cyanobacteria for salmon vs control in year 1
HC_16S_P_map_Y1_BF$k__Bacteria.p__Cyanobacteria<-as.numeric(HC_16S_P_map_Y1_BF$"k__Bacteria.p__Cyanobacteria")
HC_16S_P_map_Y1_BF$Reach<-as.factor(HC_16S_P_map_Y1_BF$Reach)
Cyo_p_Y1<-data.frame(HC_16S_P_map_Y1_BF[,c(8,11,35)])
cpy <- summarySE(Cyo_p_Y1, measurevar="k__Bacteria.p__Cyanobacteria", groupvars=c("Total_Biofilm_Growth","Reach"))
ggplot(cpy, aes(x=Total_Biofilm_Growth, y=k__Bacteria.p__Cyanobacteria, colour=Reach)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Cyanobacteria-se, ymax=k__Bacteria.p__Cyanobacteria+se), width=.1) +
  geom_line() +
  geom_point()

#Make line plot for cyanobacteria for salmon vs control in year 2
HC_16S_P_map_Y2_BF$k__Bacteria.p__Cyanobacteria<-as.numeric(HC_16S_P_map_Y2_BF$"k__Bacteria.p__Cyanobacteria")
HC_16S_P_map_Y2_BF$Reach<-as.factor(HC_16S_P_map_Y2_BF$Reach)
Cyo_p_Y2<-data.frame(HC_16S_P_map_Y2_BF[,c(8,11,35)])
cpy2 <- summarySE(Cyo_p_Y2, measurevar="k__Bacteria.p__Cyanobacteria", groupvars=c("Total_Biofilm_Growth","Reach"))
ggplot(cpy2, aes(x=Total_Biofilm_Growth, y=k__Bacteria.p__Cyanobacteria, colour=Reach)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Cyanobacteria-se, ymax=k__Bacteria.p__Cyanobacteria+se), width=.1) +
  geom_line() +
  geom_point()

#Make line plot for firmicutes for salmon vs control in year 1
HC_16S_P_map_Y1_BF$k__Bacteria.p__Firmicutes<-as.numeric(HC_16S_P_map_Y1_BF$"k__Bacteria.p__Firmicutes")
HC_16S_P_map_Y1_BF$Reach<-as.factor(HC_16S_P_map_Y1_BF$Reach)
Fir_p_Y1<-data.frame(HC_16S_P_map_Y1_BF[,c(8,11,41)])
fpy1 <- summarySE(Fir_p_Y1, measurevar="k__Bacteria.p__Firmicutes", groupvars=c("Total_Biofilm_Growth","Reach"))
ggplot(fpy1, aes(x=Total_Biofilm_Growth, y=k__Bacteria.p__Firmicutes, colour=Reach)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Firmicutes-se, ymax=k__Bacteria.p__Firmicutes+se), width=.1) +
  geom_line() +
  geom_point()

#Make line plot for firmicutes for salmon vs control in year 2
HC_16S_P_map_Y2_BF$k__Bacteria.p__Firmicutes<-as.numeric(HC_16S_P_map_Y2_BF$"k__Bacteria.p__Firmicutes")
HC_16S_P_map_Y2_BF$Reach<-as.factor(HC_16S_P_map_Y2_BF$Reach)
Fir_p_Y2<-data.frame(HC_16S_P_map_Y2_BF[,c(8,11,41)])
fpy2 <- summarySE(Fir_p_Y2, measurevar="k__Bacteria.p__Firmicutes", groupvars=c("Total_Biofilm_Growth","Reach"))
ggplot(fpy2, aes(x=Total_Biofilm_Growth, y=k__Bacteria.p__Firmicutes, colour=Reach)) + 
  geom_errorbar(aes(ymin=k__Bacteria.p__Firmicutes-se, ymax=k__Bacteria.p__Firmicutes+se), width=.1) +
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



