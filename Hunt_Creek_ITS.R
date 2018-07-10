#Analysis of ITS microbial communities at Hunt Creek 2014-2016

#########################################
#Load packages
########################################
library(reshape)
library(ggplot2)
library(vegan)
library(indicspecies)

#########################
#Functions
#########################

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

###########################
#color vectors
#############################
two_col_vec_reach<-c("skyblue3", "tomato3")

#First work with rarefication plots
HC_Shannon_ITS<-read.table("~/Desktop/shannon_ITS.txt", sep="\t", header=T)
HC_Shannon_ITS$X<-NULL
HC_Shannon_ITS$iteration<-NULL
HC_Shannon_ITS$sequences.per.sample<-as.double(HC_Shannon_ITS$sequences.per.sample)
HC_Sh_ITS_m<-melt(HC_Shannon_ITS, id.vars="sequences.per.sample")
HC_Sh_ITS_m$value<-as.double(HC_Sh_ITS_m$value)
HC_Sh_ITS_c_m<-cast(HC_Sh_ITS_m, variable ~ sequences.per.sample, fun.aggregate=mean)
HC_Sh_ITS_c_m$Calculation<-rep("mean",36)
HC_Sh_ITS_c_var<-cast(HC_Sh_ITS_m, variable ~ sequences.per.sample, fun.aggregate=var)
HC_Sh_ITS_c_var$Calculation<-rep("variance",36)
HC_Sh_ITS_c_m_var<-rbind(HC_Sh_ITS_c_m,HC_Sh_ITS_c_var)
names(HC_Sh_ITS_c_m_var)[names(HC_Sh_ITS_c_m_var)=="variable"] <- "SampleID"
HC_Sh_ITS_c_m_var$Calculation<-as.factor(HC_Sh_ITS_c_m_var$Calculation)
HC_Sh_ITS_c_m_var<-as.data.frame(HC_Sh_ITS_c_m_var)
HC_Sh_ITS_c_m_var_m<-melt(HC_Sh_ITS_c_m_var, id.vars=c("Calculation","SampleID"))
HC_Sh_ITS_c_m_var_c<-cast(HC_Sh_ITS_c_m_var_m, variable + SampleID ~ Calculation)
#Get metadata through mapping file
Hunt_Creek_ITS_map <- read.table("~/Desktop/Hunt_Creek_ITS_Map.txt", header=T)
#Merge metadata onto rarefication file
HC_Sh_ITS_map <-merge(Hunt_Creek_ITS_map, HC_Sh_ITS_c_m_var_c, by="SampleID")
names(HC_Sh_ITS_map)[names(HC_Sh_ITS_map)=="variable"] <- "Sequences_per_sample"
HC_Sh_ITS_map$Sequences_per_sample<-as.numeric(as.character(HC_Sh_ITS_map$Sequences_per_sample))
HC_Sh_ITS_map_sum_m <- summarySE(HC_Sh_ITS_map, measurevar=c("mean"), groupvars=c("Sequences_per_sample","Reach"))
HC_Sh_ITS_map_sum_v <- summarySE(HC_Sh_ITS_map, measurevar=c("variance"), groupvars=c("Sequences_per_sample","Reach"))
HC_Sh_ITS_map_sum_v$StandDev<-sqrt(HC_Sh_ITS_map_sum_v$variance)
HC_Sh_ITS_map_sum_v$StandEr<-HC_Sh_ITS_map_sum_v$StandDev/sqrt(HC_Sh_ITS_map_sum_v$N)
HC_Sh_ITS_sum_m_sd<-merge(HC_Sh_ITS_map_sum_m,HC_Sh_ITS_map_sum_v, by=0)
#make rarefication plots
ggplot(HC_Sh_ITS_sum_m_sd, aes(x=Sequences_per_sample.x, y=mean, colour=Reach.x)) + 
  geom_errorbar(aes(ymin=mean-StandEr, ymax=mean+StandEr), width=1) +
  geom_line(size=1.5) +
  geom_point(size=1.5) +
  xlab("Sequences per sample") +
  ylab("Rarefaction measure: Shannon +/- SE") +
  labs(colour = "Reach") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=two_col_vec_reach)
#####################################
#Upload data tables generated in QIIME and manipulate to make "r friendly"
#####################################

#Get ITS OTU table
Hunt_Creek_ITS<-read.table("~/Desktop/HC_ITS_OTU_Table.txt", sep="\t", header = T)

#Format data frame so the OTU.ID is row name
row.names(Hunt_Creek_ITS)<-Hunt_Creek_ITS[,1]

#Delete taxonomy and OTU.ID columns, now that OTU.ID is row name
Hunt_Creek_ITS$OTU.ID<-NULL
Hunt_Creek_ITS$taxonomy<-NULL

#transpose
HC_ITS_OTU_t<-t(Hunt_Creek_ITS)
HC_ITS_OTU_t<-data.frame(HC_ITS_OTU_t)

#Merge metadata onto data table
Hunt_Creek_ITS_map_r<-Hunt_Creek_ITS_map
row.names(Hunt_Creek_ITS_map_r)<-Hunt_Creek_ITS_map_r[,1]
HC_ITS_OTU_map <-merge(Hunt_Creek_ITS_map_r, HC_ITS_OTU_t, by=0)
str(HC_ITS_OTU_map)

#Create overall environmental data matrix for community analysis with counts
HC_ITS_env<-Hunt_Creek_ITS_map
row.names(HC_ITS_env)<-HC_ITS_env[,1]
HC_ITS_env<-HC_ITS_env[,-c(1)]

#############################################
#Analysis of all ITS community data
#############################################

#Overall permanova
adonis(HC_ITS_OTU_t ~ Reach*Year*Total_Biofilm_Growth_PostCarcass, data=HC_ITS_env, method="jaccard", permutations=999)
#No significant factors

#Indicator analysis
HC_ITS_indic<-multipatt(HC_ITS_OTU_t, HC_ITS_env$Reach, control = how(nperm=999))
summary(HC_ITS_indic)
#6 indicator OTUs detected for biofilms the salmon reaches
#They are New.CleanUp.ReferenceOTU2973, EU547495, AH008235, New.ReferenceOTU6, New.CleanUp.ReferenceOTU4525, and New.CleanUp.ReferenceOTU4648

#################
#Analysis using phyla level taxonomy table
#################

#upload phyla level info
Hunt_Creek_ITS_P<-read.table("~/Desktop/otu_table_HC_ITS_Phyla.txt", sep="\t", header = T)
#Clasify as data.frame
Hunt_Creek_ITS_P<-data.frame(Hunt_Creek_ITS_P)
#Format data frame so the taxonomy is row name
row.names(Hunt_Creek_ITS_P)<-Hunt_Creek_ITS_P[,1]
#Delete taxonomy column
Hunt_Creek_ITS_P$OTU.ID<-NULL
#transpose
HC_ITS_P_t<-t(Hunt_Creek_ITS_P)
HC_ITS_P_t<-data.frame(HC_ITS_P_t)
str(HC_ITS_P_t)
names(HC_ITS_P_t)
#Merge metadata onto data table
HC_ITS_P_map <-merge(Hunt_Creek_ITS_map_r, HC_ITS_P_t, by=0)

#Indicator analysis for phyla
HC_ITS_p_indic<-multipatt(HC_ITS_P_t, HC_ITS_env$Reach, control = how(nperm=999))
summary(HC_ITS_p_indic)
#none significant

