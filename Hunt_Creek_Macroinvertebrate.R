#Analyze macroinvertebrate community data from Hunt Creek Salmon study

###################################################
#Load packages and functions and get color vectors
##################################################

#load packages
library(vegan)
library(reshape)
library(ggplot2)
library(plyr)
library(dplyr)
library(doBy)
library(scales)
library(RColorBrewer)
library(indicspecies)
library(lme4)
library(ggpubr)
library(nlme)
library(pastecs)

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

#Make color vectors
two_col_vec <- c("black", "bisque3")
two_col_vec_ob<-c("#ef6548", "#3690c0")
eight_col_vec<-c("#f7fcfd", "#e5f5f9", "#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#005824")
eleven_col_vec <- c("gray0", "gray9", "gray18", "gray27", "gray36", "gray45", "gray54", "gray63", "gray72", "gray81", "gray90")
fifteen_col_vec<-c("#fee8c8", "#fff7fb", "#fdd49e", "#ece7f2", "#d0d1e6", "#fdbb84", "#a6bddb", "#fc8d59", "#74a9cf", "#ef6548" ,"#3690c0", "#d7301f", "#0570b0", "#990000", "#034e7b")

##############################################
#Upload data and make "r friendly"
###############################################

#Upload .csv file
Hunt_Creek_Macroinvertebrates<-read.csv("~/Documents/MSU/Research/Hunt_Creek_Salmon/Macroinvertebrates/Hunt_Creek_Macroinvertebrates.csv", sep = ",", header = T )
#Confirm header names
names(Hunt_Creek_Macroinvertebrates)
#Create sample ID name by combining environmental variables into one name in new column
Hunt_Creek_Macroinvertebrates$SampleID<-factor(paste(Hunt_Creek_Macroinvertebrates$Reach, Hunt_Creek_Macroinvertebrates$Subreach, Hunt_Creek_Macroinvertebrates$Date, Hunt_Creek_Macroinvertebrates$Type))
#Create new column that combines taxonomic variables to family level
Hunt_Creek_Macroinvertebrates$Taxonomy<-(paste(Hunt_Creek_Macroinvertebrates$Class, Hunt_Creek_Macroinvertebrates$Order, Hunt_Creek_Macroinvertebrates$Family))
#Delete pupae (named Hexapod NA NA in taxonomy column)
HC_Inverts_Subset<-subset(Hunt_Creek_Macroinvertebrates, Value >= 0 & Taxonomy != "Hexapoda NA NA" & Taxonomy !="Hexapoda Trichoptera NA" & Taxonomy != "Hexapoda Plecoptera NA" & Taxonomy !="Hexapoda Diptera NA" & Taxonomy !="Hexapoda Ephemeroptera NA")
#Convert Taxonomy charater to factor
HC_Inverts_Subset$Taxonomy<-as.factor(HC_Inverts_Subset$Taxonomy)
#create env factor
HC_Inverts_Subset$Env<-as.factor(paste(HC_Inverts_Subset$Date,HC_Inverts_Subset$Reach))
levels(HC_Inverts_Subset$Env)
HC_Inverts_Subset$Env<-revalue(HC_Inverts_Subset$Env, c("9/5/14 Control"="Control",
                                                        "9/5/14 Salmon"="Control",
                                                        "10/4/14 Control"="Control",
                                                        "10/4/14 Salmon"="Salmon",
                                                        "3/12/15 Control"="Control",
                                                        "3/12/15 Salmon"="Salmon",
                                                        "5/12/15 Control"="Control",
                                                        "5/12/15 Salmon"="Salmon",
                                                        "6/23/15 Control"="Control",
                                                        "6/23/15 Salmon"="Salmon",
                                                        "8/4/15 Control"="Control",
                                                        "8/4/15 Salmon"="Salmon",
                                                        "9/20/15 Control"="Control",
                                                        "9/20/15 Salmon"="Salmon",
                                                        "10/4/15 Control"="Control",
                                                        "10/4/15 Salmon"="Salmon",
                                                        "10/25/15 Control"="Control", 
                                                        "10/25/15 Salmon"="Salmon", 
                                                        "3/19/16 Control"="Control", 
                                                        "3/19/16 Salmon"="Salmon",
                                                        "5/17/16 Control"="Control",
                                                        "5/17/16 Salmon"="Salmon",
                                                        "6/30/16 Control"="Control",
                                                        "6/30/16 Salmon"="Salmon",
                                                        "8/15/16 Control"="Control",
                                                        "8/15/16 Salmon"="Salmon"))
#Check levels of SampleID and Taxonomy variables
str(HC_Inverts_Subset)
levels(HC_Inverts_Subset$SampleID)
levels(HC_Inverts_Subset$Taxonomy)
levels(HC_Inverts_Subset$Env)
sum(HC_Inverts_Subset$Value)

#Convert to community table for community analysis
HC_M_Matrix <- cast(HC_Inverts_Subset, Date + Reach + Subreach + Days_Since_Carcass_Introduction + Year + Days_Since_Study_Start + Env ~ Taxonomy, value = "Value", fun.aggregate =sum)
#create file with only community data
str(HC_M_Matrix)
HC_Inverts_Community <- HC_M_Matrix[,8:ncol(HC_M_Matrix)]
colSums(HC_Inverts_Community)
#Replace "na" with 0
HC_Inverts_Community[is.na(HC_Inverts_Community)]<-0
str(HC_Inverts_Community)
#Separate environmental data
HC_Inverts_Env<-HC_M_Matrix[,1:7]
str(HC_Inverts_Env)
HC_Inverts_Env$Sample_number<-as.factor(HC_Inverts_Env$Days_Since_Carcass_Introduction)
HC_Inverts_Env$Year<-as.factor(HC_Inverts_Env$Year)
levels(HC_Inverts_Env$Sample_number)
HC_Inverts_Env$Sample_number<-gsub("161", "162", HC_Inverts_Env$Sample_number)
HC_Inverts_Env$Sample_number<-gsub("265", "264", HC_Inverts_Env$Sample_number)
HC_Inverts_Env$Sample_number<-gsub("312", "306",HC_Inverts_Env$Sample_number)
HC_Inverts_Env$Sample_number<-as.factor(HC_Inverts_Env$Sample_number)
HC_Inverts_Env$Location<-paste(HC_Inverts_Env$Reach,HC_Inverts_Env$Subreach)
################################
#Create summary table for manuscript
###############################
HC_Inverts_SummaryTable_means<-aggregate(HC_M_Matrix[10:ncol(HC_M_Matrix)], by=list(Env=HC_M_Matrix$Env, Year=HC_M_Matrix$Year), FUN=mean)
HC_Inverts_SummaryTable_sums<-aggregate(HC_M_Matrix[10:ncol(HC_M_Matrix)], by=list(Env=HC_M_Matrix$Env, Year=HC_M_Matrix$Year), FUN=sum)
HC_Inverts_SummaryTable_SDs<-aggregate(HC_M_Matrix[10:ncol(HC_M_Matrix)], by=list(Env=HC_M_Matrix$Env, Year=HC_M_Matrix$Year), FUN=sd)
HC_Inverts_Summary_Table_SEs<-HC_Inverts_SummaryTable_SDs[3:ncol(HC_Inverts_SummaryTable_SDs)]/sqrt(HC_Inverts_SummaryTable_sums[3:ncol(HC_Inverts_SummaryTable_SDs)])
##############
#Most abundant invertebrate
##################
sum(HC_Inverts_SummaryTable_sums$`Hexapoda Diptera Chironomidae`)/sum(HC_Inverts_SummaryTable_sums[,2:ncol(HC_Inverts_SummaryTable_sums)])

##########################################
#PERMANOVA for macroinvertebrate data
###########################################

#Permanova for Year, treatment, and days since carcass introduction
adonis(HC_Inverts_Community ~ Reach*Sample_number*Year+Location, data=HC_Inverts_Env, method="bray", permutations=999)

###########################################
#Separate and analyze for year 1 and 2 separately
################################################

#Create year 1 macro data table
HC_M_Matrix_Y1<-subset(HC_M_Matrix, Year==1)
HC_Inverts_Community_Y1 <- HC_M_Matrix_Y1[,8:ncol(HC_M_Matrix_Y1)]
any(is.na(HC_Inverts_Community_Y1))
str(HC_Inverts_Community_Y1)
#Separate environmental data
HC_Inverts_Env_Y1<-HC_M_Matrix_Y1[,1:7]
str(HC_Inverts_Env_Y1)
HC_Inverts_Env_Y1$Location<-paste(HC_Inverts_Env_Y1$Reach,HC_Inverts_Env_Y1$Subreach)
#Create Year 2 macro data table
HC_M_Matrix_Y2<-subset(HC_M_Matrix, Year==2)
HC_Inverts_Community_Y2 <- HC_M_Matrix_Y2[,8:ncol(HC_M_Matrix_Y2)]
any(is.na(HC_Inverts_Community_Y2))
str(HC_Inverts_Community_Y2)
#Separate environmental data
HC_Inverts_Env_Y2<-HC_M_Matrix_Y2[,1:7]
str(HC_Inverts_Env_Y2)
HC_Inverts_Env_Y2_Env<-HC_Inverts_Env_Y2$Env
HC_Inverts_Env_Y2$Location<-paste(HC_Inverts_Env_Y2$Reach,HC_Inverts_Env_Y2$Subreach)
HC_Inverts_ReachDate_Y2<-paste(HC_Inverts_Env_Y2$Env,HC_Inverts_Env_Y2$Days_Since_Study_Start)
#PERMANOVAs for each year
adonis(HC_Inverts_Community_Y1 ~ Reach*Days_Since_Study_Start+Location, data=HC_Inverts_Env_Y1, method="bray", permutations=999)
#Time significant
adonis(HC_Inverts_Community_Y2 ~ Reach*Days_Since_Study_Start+Location, data=HC_Inverts_Env_Y2, method="bray", permutations=999)
#Time significant

#indicator taxa analysis for time
HC_Inverts_Community_Y2_indic<-signassoc(HC_Inverts_Community_Y2, cluster=HC_Inverts_Env_Y2_Env,  mode=0, alternative = "two.sided",control = how(nperm=9999))
HC_Inverts_Community_Y2_indic_sig<-subset(HC_Inverts_Community_Y2_indic, psidak<=0.05)

########################################
#Stacked bar graphs
######################################

#Use melted data set
ggplot(Hunt_Creek_Macroinvertebrates,aes(x = Reach, y = Value,fill = Taxonomy)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = percent_format())

#Set up factor for Salmon vs Control called stream_site
levels(HC_M_Matrix$Reach)
stream_site=factor(HC_Inverts_Env$Reach, c("Salmon", "Control"))
env=factor(HC_Inverts_Env$Reach, c("Salmon","Control"))
#Set up factor for Days_Since_Carcass_Introduction
str(HC_M_Matrix$Days_Since_Carcass_Introduction)
Days_Since_Introduction_n=as.numeric(HC_Inverts_Env$Days_Since_Carcass_Introduction)
Days_Since_Introduction_f=factor(HC_Inverts_Env$Days_Since_Carcass_Introduction)
levels(Days_Since_Introduction_f)
Days_Since_Introduction_f<-as.factor(gsub("161", "162", Days_Since_Introduction_f))
Days_Since_Introduction_f<-as.factor(gsub("265", "264", Days_Since_Introduction_f))
Days_Since_Introduction_f<-as.factor(gsub("312", "306", Days_Since_Introduction_f))
Days_Carcass<-factor(Days_Since_Introduction_f, c("0", "1", "16", "162", "222", "264", "306", "352"))

#Set up factor for Days_Since_Carcass_Introduction_f_stream_site
D_S_C_I_f_s_s <- factor(paste(HC_Inverts_Env$Days_Since_Carcass_Introduction, HC_Inverts_Env$Reach))
levels(D_S_C_I_f_s_s)
D_S_C_I_f_s_s<-as.factor(gsub("161 Control", "162 Control", D_S_C_I_f_s_s))
D_S_C_I_f_s_s<-as.factor(gsub("161 Salmon", "162 Salmon", D_S_C_I_f_s_s))
D_S_C_I_f_s_s<-as.factor(gsub("265 Salmon", "264 Salmon", D_S_C_I_f_s_s))
D_S_C_I_f_s_s<-as.factor(gsub("265 Control", "264 Control", D_S_C_I_f_s_s))
D_S_C_I_f_s_s<-as.factor(gsub("312 Control", "306 Control", D_S_C_I_f_s_s))
D_S_C_I_f_s_s<-as.factor(gsub("312 Salmon", "306 Salmon", D_S_C_I_f_s_s))
Days_Carcass_Reach<-factor(D_S_C_I_f_s_s, c("0 Control", "0 Salmon", "1 Control", "1 Salmon", "16 Salmon", "162 Control", "162 Salmon", "222 Control", "222 Salmon", "264 Control", "264 Salmon", "306 Control", "306 Salmon", "352 Control", "352 Salmon"))

#make factor for year
Year<-factor(HC_Inverts_Env$Date)
levels(Year)
Year<-as.factor(gsub("10/25/15", "2", Year))
Year<-as.factor(gsub("10/4/15", "2", Year))
Year<-as.factor(gsub("3/19/16", "2", Year))
Year<-as.factor(gsub("5/17/16", "2", Year))
Year<-as.factor(gsub("6/30/16", "2", Year))
Year<-as.factor(gsub("10/4/14", "1", Year))
Year<-as.factor(gsub("3/12/15", "1", Year))
Year<-as.factor(gsub("3/19/15", "2", Year)) #typo in data entry will fix master later
Year<-as.factor(gsub("5/12/15", "1", Year))
Year<-as.factor(gsub("6/23/15", "1", Year))
Year<-as.factor(gsub("8/4/15", "1", Year))
Year<-as.factor(gsub("9/20/15", "1", Year))
Year<-as.factor(gsub("9/5/14", "1", Year))
Year<-as.factor(gsub("8/15/16", "2", Year))

#NMDS analysis
HC_Inverts_Community<-as.matrix(HC_Inverts_Community)
HC_Macroinvertebrate_NMDS<-metaMDS(HC_Inverts_Community, distance="bray")

#Stressplot macroinvertebrate Nmds
stressplot(HC_Macroinvertebrate_NMDS)

#NMDS plot for env Salmon vs Control
ordiplot(HC_Macroinvertebrate_NMDS, type="n", main="Invertebrate")
with(HC_Macroinvertebrate_NMDS, points(HC_Macroinvertebrate_NMDS, display="sites", col=two_col_vec[env], pch=19, pt.bg=two_col_vec))
with(HC_Macroinvertebrate_NMDS, legend("topleft", legend=levels(env), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, env, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Salmon"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, env, kind="se", conf=0.95, lwd=2, col="bisque3", show.groups = "Control"))

#NMDS plot for Days_Since_Introduction_f
ordiplot(HC_Macroinvertebrate_NMDS, type="n", main="Macroinvertebrate Days since introduction")
with(HC_Macroinvertebrate_NMDS, points(HC_Macroinvertebrate_NMDS, display="sites", col=eight_col_vec[Days_Carcass], pch=19, pt.bg=eight_col_vec))
with(HC_Macroinvertebrate_NMDS, legend("topleft", legend=levels(Days_Carcass), bty="n", col=eight_col_vec, pch=19, pt.bg=eight_col_vec))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="#f7fcfd", show.groups = "0"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="#e5f5f9", show.groups = "1"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="#ccece6", show.groups = "16"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="#99d8c9", show.groups = "162"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="#66c2a4", show.groups = "222"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="#41ae76", show.groups = "264"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="#238b45", show.groups = "306"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, Days_Carcass, kind="se", conf=0.95, lwd=2, col="#005824", show.groups = "352"))

#NMDS plot by year
ordiplot(HC_Macroinvertebrate_NMDS, type="n", main="Macroinvertebrate Days since introduction")
with(HC_Macroinvertebrate_NMDS, points(HC_Macroinvertebrate_NMDS, display="sites", col=two_col_vec[Year], pch=19, pt.bg=two_col_vec))
with(HC_Macroinvertebrate_NMDS, legend("topleft", legend=levels(Year), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))


#NMDS plot for days carcass reach interaction
ordiplot(HC_Macroinvertebrate_NMDS, type="n", main="Macroinvertebrate Days since introduction Reach")
with(HC_Macroinvertebrate_NMDS, points(HC_Macroinvertebrate_NMDS, display="sites", col=fifteen_col_vec[Days_Carcass_Reach], pch=19, pt.bg=fifteen_col_vec))
with(HC_Macroinvertebrate_NMDS, legend("topleft", legend=levels(stream_site), bty="n", col=two_col_vec_ob, pch=19, pt.bg=two_col_vec_ob))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, Year, kind="se", conf=0.95, lwd=2, col="black", show.groups = "1"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, Year, kind="se", conf=0.95, lwd=2, col="bisque3", show.groups = "2"))

##############
####Diversity metrics
###############
#Calculate richess for year 1
HC_Inverts_Env_Y1$Richness<-rowSums(HC_Inverts_Community_Y1 > 0)
#create model
shapiro.test(HC_Inverts_Env_Y1$Richness)
#normal W = 0.95203, p-value = 0.07634
Rich_y1.lm = glm(Richness ~ Reach, data=HC_Inverts_Env_Y1, family=poisson()) 
#Check model assumptions
#looks okay
summary(Rich_y1.lm)
#significant interaction 

#Calculate richess for year 1
HC_Inverts_Env_Y2$Richness<-rowSums(HC_Inverts_Community_Y2 > 0)
#create model
Rich_y2.lm = glm(Richness ~ Reach*Days_Since_Study_Start, data=HC_Inverts_Env_Y2, family=poisson()) 
#Check model assumptions
#looks okay
summary(Rich_y2.lm)
#nothing significant

#plot richness
HC_Inverts_Env$Richness<-rowSums(HC_Inverts_Community > 0)
rich<-summarySE(HC_Inverts_Env, measurevar="Richness", groupvars=c("Reach", "Days_Since_Study_Start"))
ggplot(rich, aes(x=Days_Since_Study_Start, y=Richness, linetype=Reach)) + 
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), width=.1, linetype=1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days since study start") +
  ylab("Species Richness") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(5,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")

#Calculate Simpson's diversity for year 1
HC_Inverts_Env_Y1$Simp<-diversity(HC_Inverts_Community_Y1, index="simpson")
shapiro.test(HC_Inverts_Env_Y1$Simp)
#normally distributed W = 0.95624, p-value = 0.1082
#create model
Div_y1.lm = lme(Simp ~ Reach*Days_Since_Study_Start, data=HC_Inverts_Env_Y1, random=~1|Location) 
#Check model assumptions
#looks good
summary(Div_y1.lm)
#Nothing significant

#Calculate Simpson's diversity for year 2
HC_Inverts_Env_Y2$Simp<-diversity(HC_Inverts_Community_Y2, index="simpson")
shapiro.test(HC_Inverts_Env_Y2$Simp)
#Not normal W = 0.83367, p-value = 0.0001519
hist(HC_Inverts_Env_Y2$Simp)
#left skewed
max(HC_Inverts_Env_Y2$Simp)
#0.86
HC_Inverts_Env_Y2$Simpr<-0.91-HC_Inverts_Env_Y2$Simp
hist(HC_Inverts_Env_Y2$Simpr)
HC_Inverts_Env_Y2$srSimpr<-sqrt(HC_Inverts_Env_Y2$Simpr)
shapiro.test(HC_Inverts_Env_Y2$srSimpr)
#not normal W = 0.91078, p-value = 0.01026
hist(HC_Inverts_Env_Y2$srSimpr)
HC_Inverts_Env_Y2$log10Simpr<-log10(HC_Inverts_Env_Y2$Simpr)
shapiro.test(HC_Inverts_Env_Y2$log10Simpr)
#normal W = 0.95595, p-value = 0.198
#create model
Div_y2.lm = lme(log10Simpr ~ Reach*Days_Since_Study_Start, data=HC_Inverts_Env_Y2, random=~1|Location) 
#Check model assumptions
#looks okay
summary(Div_y2.lm)
#time and time:reach interaction significant

#plot diversity
HC_Inverts_Env$Simp<-diversity(HC_Inverts_Community, index="simpson")
simp<-summarySE(HC_Inverts_Env, measurevar="Simp", groupvars=c("Reach", "Days_Since_Study_Start"))
ggplot(simp, aes(x=Days_Since_Study_Start, y=Simp, linetype=Reach)) + 
  geom_errorbar(aes(ymin=Simp-se, ymax=Simp+se), width=.1, linetype=1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days Since Study Start") +
  ylab("Simpson's Diversity Index") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(5,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position = "none") +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")

##################
#Density metric
###############
#Calculate density for year 1
HC_Inverts_Env_Y1$Density<-rowSums(HC_Inverts_Community_Y1)
shapiro.test(HC_Inverts_Env_Y1$Density)
# not normal W = 0.79814, p-value = 4.044e-06
hist(HC_Inverts_Env_Y1$Density)
#skewed right
HC_Inverts_Env_Y1$log10density<-log10(HC_Inverts_Env_Y1$Density)
shapiro.test(HC_Inverts_Env_Y1$log10density)
#normal W = 0.94725, p-value = 0.0515
#create model
Dens_y1.lm = lm(log10density ~ Reach*Days_Since_Carcass_Introduction, data=HC_Inverts_Env_Y1) 
#looks okay
summary(Dens_y1.lm)
#only time significant Days_Since_Carcass_Introduction estimate 0.0021641  std error 0.0006257   t value 3.459  p 0.00135 ** 

#Calculate density for year 2
HC_Inverts_Env_Y2$Density<-rowSums(HC_Inverts_Community_Y2)
shapiro.test(HC_Inverts_Env_Y2$Density)
#not normal W = 0.87465, p-value = 0.001248
hist(HC_Inverts_Env_Y2$Density)
#right skewed
HC_Inverts_Env_Y2$log10density<-log10(HC_Inverts_Env_Y2$Density)
shapiro.test(HC_Inverts_Env_Y2$log10density)
#not normal but better W = 0.92841, p-value = 0.03154
#create model
Dens_y2.lm = lm(log10density ~ Reach*Days_Since_Carcass_Introduction, data=HC_Inverts_Env_Y2) 
#Check model assumptions
#looks okay
summary(Dens_y2.lm)
#nothing significant

#Calculate heptageniidae density for year 1
#create model
HC_M_Matrix_Y1$Location<-paste(HC_M_Matrix_Y1$Reach,HC_M_Matrix_Y1$Subreach)
colnames(HC_M_Matrix_Y1)[colnames(HC_M_Matrix_Y1)=="Hexapoda Ephemeroptera Heptageniidae"] <- "Hexapoda_Ephemeroptera_Heptageniidae"
hep_y1.lm = glm(Hexapoda_Ephemeroptera_Heptageniidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y1, family=poisson()) 
#Check model assumptions
#looks okay
summary(hep_y1.lm)
#reach and time significant

HC_M_Matrix_Y2$Location<-paste(HC_M_Matrix_Y2$Reach,HC_M_Matrix_Y2$Subreach)
colnames(HC_M_Matrix_Y2)[colnames(HC_M_Matrix_Y2)=="Hexapoda Ephemeroptera Heptageniidae"] <- "Hexapoda_Ephemeroptera_Heptageniidae"
hep_y2.lm = glm(Hexapoda_Ephemeroptera_Heptageniidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y2, family=poisson()) 
#Check model assumptions
#looks okay
summary(hep_y2.lm)
#reach significant

#Calculate baetidae density for year 1
#create model
colnames(HC_M_Matrix_Y1)[colnames(HC_M_Matrix_Y1)=="Hexapoda Ephemeroptera Baetidae"] <- "Hexapoda_Ephemeroptera_Baetidae"
ba_y1.lm = glm(Hexapoda_Ephemeroptera_Baetidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y1, family=poisson()) 
#Check model assumptions
#looks okay
summary(ba_y1.lm)
#time and time:reach interaction significant significant

#Calculate baetidae density for year 2
#create model
colnames(HC_M_Matrix_Y2)[colnames(HC_M_Matrix_Y2)=="Hexapoda Ephemeroptera Baetidae"] <- "Hexapoda_Ephemeroptera_Baetidae"
ba_y2.lm = glm(Hexapoda_Ephemeroptera_Baetidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y2, family=poisson()) 
#Check model assumptions
#looks okay
summary(ba_y2.lm)
#time and time:reach interaction significant significant

#Calculate simuliidae density for year 1
#create model
colnames(HC_M_Matrix_Y1)[colnames(HC_M_Matrix_Y1)=="Hexapoda Diptera Simuliidae"] <- "Hexapoda_Diptera_Simuliidae"
si_y1.lm = glm(Hexapoda_Diptera_Simuliidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y1, family=poisson()) 
#Check model assumptions
#looks okay
summary(si_y1.lm)
#time and time:reach interaction significant significant, yet doesn't increase by 1 individual in over 500 days so not biologically relevant
si_y1.lmnt = glm(Hexapoda_Diptera_Simuliidae ~ Reach+Reach:Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y1, family=poisson()) 
summary(si_y1.lmnt)
#time:reach interaction significant

#Calculate simuliidae density for year 2
#create model
colnames(HC_M_Matrix_Y2)[colnames(HC_M_Matrix_Y2)=="Hexapoda Diptera Simuliidae"] <- "Hexapoda_Diptera_Simuliidae"
si_y2.lm = glm(Hexapoda_Diptera_Simuliidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y2, family=poisson()) 
#Check model assumptions
#looks okay
summary(si_y2.lm)
#time significant

#now do relative abundances
#Calculate heptageniidae relative abundance for year 1
#create model
HC_M_Matrix_Y1$RowSums<-rowSums(HC_M_Matrix_Y1[,8:(ncol(HC_M_Matrix_Y1)-1)])
HC_M_Matrix_Y1$Heptageniidae<-(HC_M_Matrix_Y1$"Hexapoda_Ephemeroptera_Heptageniidae"/HC_M_Matrix_Y1$RowSums)
hepra_y1.lm = lm(Heptageniidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y1) 
#Check model assumptions
#looks okay
summary(hepra_y1.lm)
#reach:time interaction significant

#Calculate heptageniidae relative abundance for year 2
#create model
HC_M_Matrix_Y2$RowSums<-rowSums(HC_M_Matrix_Y2[,8:(ncol(HC_M_Matrix_Y2)-1)])
HC_M_Matrix_Y2$Heptageniidae<-(HC_M_Matrix_Y2$"Hexapoda_Ephemeroptera_Heptageniidae"/HC_M_Matrix_Y2$RowSums)
hepra_y2.lm = lm(Heptageniidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y2) 
#Check model assumptions
#looks okay
summary(hepra_y2.lm)
#nothing significant

#Calculate baetis relative abundance for year 1
#create model
HC_M_Matrix_Y1$Baetidae<-(HC_M_Matrix_Y1$"Hexapoda_Ephemeroptera_Baetidae"/HC_M_Matrix_Y1$RowSums)
bara_y1.lm = lm(Baetidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y1) 
#Check model assumptions
#looks okay
summary(bara_y1.lm)
#nothing significant

#Calculate baetis relative abundance for year 2
#create model
HC_M_Matrix_Y2$Baetidae<-(HC_M_Matrix_Y2$"Hexapoda_Ephemeroptera_Baetidae"/HC_M_Matrix_Y2$RowSums)
bara_y2.lm = lm(Baetidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y2) 
#Check model assumptions
#looks okay
summary(bara_y2.lm)
#nothing significant

#Calculate simulid relative abundance for year 1
#create model
HC_M_Matrix_Y1$Simuliidae<-(HC_M_Matrix_Y1$"Hexapoda_Diptera_Simuliidae"/HC_M_Matrix_Y1$RowSums)
sira_y1.lm = lm(Simuliidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y1) 
#Check model assumptions
#looks okay
summary(sira_y1.lm)
#nothing significant

#Calculate simulid relative abundance for year 2
#create model
HC_M_Matrix_Y2$Simuliidae<-(HC_M_Matrix_Y2$"Hexapoda_Diptera_Simuliidae"/HC_M_Matrix_Y2$RowSums)
sira_y2.lm = lm(Simuliidae ~ Reach*Days_Since_Carcass_Introduction, data=HC_M_Matrix_Y2) 
#Check model assumptions
#looks okay
summary(sira_y2.lm)
#nothing significant

#############
#Now on to populations rather than communities
#################
#find top taxa
totals<-rbind(HC_Inverts_Community, colSums(HC_Inverts_Community))
totals<-totals[-c(1:58),]
rowSums(totals)
sort(totals,decreasing=TRUE)[1:6]

#Make line plot for Baetids
HC_M_Matrix$RowSums<-rowSums(HC_M_Matrix[,8:ncol(HC_M_Matrix)])
HC_M_Matrix$Treatment<-HC_M_Matrix$Reach
HC_M_Matrix$Baetidae<-(HC_M_Matrix$"Hexapoda Ephemeroptera Baetidae"/HC_M_Matrix$RowSums)
Baetidaeplot<-data.frame(HC_M_Matrix[,c(57,6,58)])
Baetidaeplot$Days_Since_Study_Start<-as.numeric(Baetidaeplot$Days_Since_Study_Start)
btd<-summarySE(Baetidaeplot, measurevar="Baetidae", groupvars=c("Days_Since_Study_Start", "Treatment"))
ggplot(btd, aes(x=Days_Since_Study_Start, y=Baetidae, colour=Treatment)) + 
  geom_errorbar(aes(ymin=Baetidae-se, ymax=Baetidae+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days") +
  ylab("Baetidae relative abundance +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))
colnames(HC_M_Matrix)[colnames(HC_M_Matrix)=="Hexapoda Ephemeroptera Baetidae"] <- "Hexapoda_Ephemeroptera_Baetidae"
btddens<-summarySE(HC_M_Matrix, measurevar="Hexapoda_Ephemeroptera_Baetidae", groupvars=c("Days_Since_Study_Start", "Reach"))
ggplot(btddens, aes(x=Days_Since_Study_Start, y=Hexapoda_Ephemeroptera_Baetidae, colour=Reach)) + 
  geom_errorbar(aes(ymin=Hexapoda_Ephemeroptera_Baetidae-se, ymax=Hexapoda_Ephemeroptera_Baetidae+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days") +
  ylab("Baetidae density +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))
HC_M_Matrix_sal<-subset(HC_M_Matrix, Env=="Salmon")
HC_M_Matrix_cont<-subset(HC_M_Matrix, Env=="Control")
stat.desc(HC_M_Matrix_sal$`Hexapoda_Ephemeroptera_Baetidae`)
stat.desc(HC_M_Matrix_cont$`Hexapoda_Ephemeroptera_Baetidae`)
#Make line plot for Simuliids
HC_M_Matrix$Simuliidae<-(HC_M_Matrix$"Hexapoda Diptera Simuliidae"/HC_M_Matrix$RowSums)
Simuliidaeplot<-data.frame(HC_M_Matrix[,c(57,6,59)])
Simuliidaeplot$Days_Since_Study_Start<-as.numeric(Simuliidaeplot$Days_Since_Study_Start)
sim<-summarySE(Simuliidaeplot, measurevar="Simuliidae", groupvars=c("Days_Since_Study_Start", "Treatment"))
ggplot(sim, aes(x=Days_Since_Study_Start, y=Simuliidae, colour=Treatment)) + 
  geom_errorbar(aes(ymin=Simuliidae-se, ymax=Simuliidae+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days") +
  ylab("Simuliidae relative abundance +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))
colnames(HC_M_Matrix)[colnames(HC_M_Matrix)=="Hexapoda Diptera Simuliidae"] <- "Hexapoda_Diptera_Simuliidae"
sidens<-summarySE(HC_M_Matrix, measurevar="Hexapoda_Diptera_Simuliidae", groupvars=c("Days_Since_Study_Start", "Reach"))
ggplot(sidens, aes(x=Days_Since_Study_Start, y=Hexapoda_Diptera_Simuliidae, colour=Reach)) + 
  geom_errorbar(aes(ymin=Hexapoda_Diptera_Simuliidae-se, ymax=Hexapoda_Diptera_Simuliidae+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days") +
  ylab("Simuliidae density +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))
#Make line plot for Polycentropidae
HC_M_Matrix$Cyrnellusfraternus<-(HC_M_Matrix$"Hexapoda Trichoptera Polycentropidae"/HC_M_Matrix$RowSums)
Cfplot<-data.frame(HC_M_Matrix[,c(6,57,60)])
Cfplot$Days_Since_Study_Start<-as.numeric(Cfplot$Days_Since_Study_Start)
Cf<-summarySE(Cfplot, measurevar="Cyrnellusfraternus", groupvars=c("Days_Since_Study_Start", "Treatment"))
ggplot(Cf, aes(x=Days_Since_Study_Start, y=Cyrnellusfraternus, colour=Treatment)) + 
  geom_errorbar(aes(ymin=Cyrnellusfraternus-se, ymax=Cyrnellusfraternus+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days") +
  ylab("Cyrnellus fraternus relative abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_color_manual(values=c("skyblue3", "tomato3")) +
  scale_x_continuous(breaks=seq(0,700,100))

#Make line plot for Brachycentridae
HC_M_Matrix$Brachycentridae<-(HC_M_Matrix$"Hexapoda Trichoptera Brachycentridae"/HC_M_Matrix$RowSums)
Brachyplot<-data.frame(HC_M_Matrix[,c(57,6,61)])
Brachyplot$Days_Since_Study_Start<-as.numeric(Brachyplot$Days_Since_Study_Start)
bra<-summarySE(Brachyplot, measurevar="Brachycentridae", groupvars=c("Days_Since_Study_Start", "Treatment"))
bra$Treatment<-factor(bra$Treatment, c("Salmon","Control"))
Bralab<-expression(paste("Mean ", italic("Brachycentrus"), " relative abundance"))
ggplot(bra, aes(x=Days_Since_Study_Start, y=Brachycentridae, linetype=Treatment)) + 
  geom_errorbar(aes(ymin=Brachycentridae-se, ymax=Brachycentridae+se), width=.1, linetype=1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days since study start") +
  ylab(Bralab) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(5,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.position=0) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")
HC_M_Matrix_Y2_sal<-subset(HC_M_Matrix_Y2, Reach=="Salmon")
HC_M_Matrix_Y2_cont<-subset(HC_M_Matrix_Y2, Reach=="Control")
stat.desc(HC_M_Matrix_Y2_sal$`Hexapoda Trichoptera Brachycentridae`)
stat.desc(HC_M_Matrix_Y2_cont$`Hexapoda Trichoptera Brachycentridae`)


#Make line plot for Heptageniidae
HC_M_Matrix$Heptageniidae<-(HC_M_Matrix$'Hexapoda Ephemeroptera Heptageniidae'/HC_M_Matrix$RowSums)
Heptageniidplot<-data.frame(HC_M_Matrix[,c(57,6,62)])
Heptageniidplot$Days_Since_Study_Start<-as.numeric(Heptageniidplot$Days_Since_Study_Start)
hep<-summarySE(Heptageniidplot, measurevar="Heptageniidae", groupvars=c("Days_Since_Study_Start", "Treatment"))
hep$Treatment<-factor(hep$Treatment, c("Salmon","Control"))
Hepylab<-expression(paste("Mean ", italic("Heptagenia"), " Rel. Abund. (±SE)"))
ggplot(hep, aes(x=Days_Since_Study_Start, y=Heptageniidae, linetype=Treatment)) + 
  geom_errorbar(aes(ymin=Heptageniidae-se, ymax=Heptageniidae+se), width=.1, linetype=1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days since study start") +
  ylab(Hepylab) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")
colnames(HC_M_Matrix)[colnames(HC_M_Matrix)=="Hexapoda Ephemeroptera Heptageniidae"] <- "Hexapoda_Ephemeroptera_Heptageniidae"
hepdens<-summarySE(HC_M_Matrix, measurevar="Hexapoda_Ephemeroptera_Heptageniidae",groupvars=c("Days_Since_Study_Start","Reach"))
hepdens$Treatment<-factor(hepdens$Reach, c("Salmon","Control"))
Hepdensylab<-expression(paste("Mean ", italic("Heptagenia"), " density"))
ggplot(hepdens, aes(x=Days_Since_Study_Start, y=Hexapoda_Ephemeroptera_Heptageniidae, linetype=Treatment)) + 
  geom_errorbar(aes(ymin=Hexapoda_Ephemeroptera_Heptageniidae-se, ymax=Hexapoda_Ephemeroptera_Heptageniidae+se), width=.1, linetype=1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days since study start") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(5,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  scale_y_continuous(Hepdensylab, expand = c(0,0)) +
  annotate("rect", xmin = 30, xmax = 90, ymin = 0, ymax = 70, alpha=0.5) +
  annotate("rect", xmin = 400, xmax = 480, ymin = 0, ymax = 70, alpha=0.5) +
  geom_vline(xintercept = c(30,400), linetype = "blank", colour = "black", size=1.5)
stat.desc(HC_M_Matrix_sal$`Hexapoda_Ephemeroptera_Heptageniidae`)
stat.desc(HC_M_Matrix_cont$`Hexapoda_Ephemeroptera_Heptageniidae`)

#Make line plot for Elmidae
HC_M_Matrix$Elmidae<-(HC_M_Matrix$"Hexapoda Coleoptera Elmidae"/HC_M_Matrix$RowSums)
Elmidplot<-data.frame(HC_M_Matrix[,c(57,6,63)])
Elmidplot$Days_Since_Study_Start<-as.numeric(Elmidplot$Days_Since_Study_Start)
El<-summarySE(Elmidplot, measurevar="Elmidae", groupvars=c("Days_Since_Study_Start", "Treatment"))
ggplot(El, aes(x=Days_Since_Study_Start, y=Elmidae, linetype=Treatment)) + 
  geom_errorbar(aes(ymin=Elmidae-se, ymax=Elmidae+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days since study start") +
  ylab("Elmidae relative abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")

#Make line plot for Oligochaeta
HC_M_Matrix$Oligochaeta<-(HC_M_Matrix$"Oligochaeta NA NA"/HC_M_Matrix$RowSums)
Oligoplot<-data.frame(HC_M_Matrix[,c(57,6,64)])
Oligoplot$Days_Since_Study_Start<-as.numeric(Oligoplot$Days_Since_Study_Start)
Ol<-summarySE(Oligoplot, measurevar="Oligochaeta", groupvars=c("Days_Since_Study_Start", "Treatment"))
ggplot(Ol, aes(x=Days_Since_Study_Start, y=Oligochaeta, linetype=Treatment)) + 
  geom_errorbar(aes(ymin=Oligochaeta-se, ymax=Oligochaeta+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days since study start") +
  ylab("Oligochaeta relative abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")
#Hydraarachnae
HC_M_Matrix$Hydrachnidae<-(HC_M_Matrix$"Arachnida Trombidiformes Hydrachnidae"/HC_M_Matrix$RowSums)
Hydplot<-data.frame(HC_M_Matrix[,c(57,6,65)])
Hydplot$Days_Since_Study_Start<-as.numeric(Hydplot$Days_Since_Study_Start)
Hy<-summarySE(Hydplot, measurevar="Hydrachnidae", groupvars=c("Days_Since_Study_Start", "Treatment"))
ggplot(Hy, aes(x=Days_Since_Study_Start, y=Hydrachnidae, linetype=Treatment)) + 
  geom_errorbar(aes(ymin=Hydrachnidae-se, ymax=Hydrachnidae+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days since study start") +
  ylab("Hydrachnidae relative abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")
#Chironomidae
HC_M_Matrix$Chironomidae<-(HC_M_Matrix$"Hexapoda Diptera Chironomidae"/HC_M_Matrix$RowSums)
Chiplot<-data.frame(HC_M_Matrix[,c(57,6,66)])
Chiplot$Days_Since_Study_Start<-as.numeric(Chiplot$Days_Since_Study_Start)
Chi<-summarySE(Chiplot, measurevar="Chironomidae", groupvars=c("Days_Since_Study_Start", "Treatment"))
ggplot(Chi, aes(x=Days_Since_Study_Start, y=Chironomidae, linetype=Treatment)) + 
  geom_errorbar(aes(ymin=Chironomidae-se, ymax=Chironomidae+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days since study start") +
  ylab("Chironomidae relative abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")
stat.desc(Chiplot$Chironomidae)

#Make faceted line plot for heptagenia, baetis and stegopterna
hep$Taxa<-rep("paste('B. ',italic('Heptagenia'))",26)
colnames(hep)[colnames(hep)=="Heptageniidae"] <- "Mean"
sim$Taxa<-rep("paste('D. ',italic('Stegopterna'))",26)
colnames(sim)[colnames(sim)=="Simuliidae"] <- "Mean"
btd$Taxa<-rep("paste('C. ',italic('Baetis'))",26)
colnames(btd)[colnames(btd)=="Baetidae"] <- "Mean"
bra$Taxa<-rep("paste('A. ',italic('Brachycentrus'))",26)
colnames(bra)[colnames(bra)=="Brachycentridae"] <- "Mean"
indinv<-rbind(hep,sim,btd,bra)
indinv$Taxa<-factor(indinv$Taxa, levels =c("paste('A. ',italic('Brachycentrus'))","paste('B. ',italic('Heptagenia'))","paste('C. ',italic('Baetis'))","paste('D. ',italic('Stegopterna'))"))
ggplot(indinv, aes(x=Days_Since_Study_Start, y=Mean, linetype=Treatment)) + 
  geom_line(size=6) +
  geom_pointrange(aes(ymin=Mean-se, ymax =Mean+se), size=2.5) +
  xlab("Days Since Study Start") +
  ylab("Mean Relative Abundance (+/- SEM)") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(family="Times New Roman", size=80, margin=margin(t=20,r=0,b=20,l=0)),
        axis.title.y=element_text(family="Times New Roman", size=80),
        axis.text.x=element_text(family="Times New Roman", size=56),
        axis.text.y = element_text(family="Times New Roman", size=56),
        legend.position = "none",
        strip.text.x = element_text(family="Times New Roman", face="bold", size = 56)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black") +
  facet_wrap( ~ Taxa, ncol=1, scales="free_y", labeller=label_parsed)

