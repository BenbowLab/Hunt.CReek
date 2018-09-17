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
#Check levels of SampleID and Taxonomy variables
str(HC_Inverts_Subset)
levels(HC_Inverts_Subset$SampleID)
levels(HC_Inverts_Subset$Taxonomy)
sum(HC_Inverts_Subset$Value)

#Convert to community table for community analysis
HC_M_Matrix <- cast(HC_Inverts_Subset, Date + Reach + Subreach + Days_Since_Carcass_Introduction + Year + Days_Since_Study_Start ~ Taxonomy, value = "Value", fun.aggregate =sum)
#create file with only community data
str(HC_M_Matrix)
HC_Inverts_Community <- HC_M_Matrix[,7:ncol(HC_M_Matrix)]
colSums(HC_Inverts_Community)
#Replace "na" with 0
HC_Inverts_Community[is.na(HC_Inverts_Community)]<-0
str(HC_Inverts_Community)
#Separate environmental data
HC_Inverts_Env<-HC_M_Matrix[,1:6]
str(HC_Inverts_Env)
HC_Inverts_Env$Sample_number<-as.factor(HC_Inverts_Env$Days_Since_Carcass_Introduction)
HC_Inverts_Env$Year<-as.factor(HC_Inverts_Env$Year)
levels(HC_Inverts_Env$Sample_number)
HC_Inverts_Env$Sample_number<-gsub("161", "162", HC_Inverts_Env$Sample_number)
HC_Inverts_Env$Sample_number<-gsub("265", "264", HC_Inverts_Env$Sample_number)
HC_Inverts_Env$Sample_number<-gsub("312", "306",HC_Inverts_Env$Sample_number)
HC_Inverts_Env$Sample_number<-as.factor(HC_Inverts_Env$Sample_number)
##########################################
#PERMANOVA for macroinvertebrate data
###########################################

#Permanova for Year, treatment, and days since carcass introduction
adonis(HC_Inverts_Community ~ Reach*Sample_number*Year, data=HC_Inverts_Env, method="bray", permutations=999)

###########################################
#Separate and analyze for year 1 and 2 separately
################################################

#Create year 1 macro data table
HC_M_Matrix_Y1<-subset(HC_M_Matrix, Year==1)
HC_Inverts_Community_Y1 <- HC_M_Matrix_Y1[,7:ncol(HC_M_Matrix_Y1)]
any(is.na(HC_Inverts_Community_Y1))
str(HC_Inverts_Community_Y1)
#Separate environmental data
HC_Inverts_Env_Y1<-HC_M_Matrix_Y1[,1:6]
str(HC_Inverts_Env_Y1)
HC_Inverts_Env_Y1$Days_Since_Carcass_Introduction<-as.factor(HC_Inverts_Env_Y1$Days_Since_Carcass_Introduction)

#Create Year 2 macro data table
HC_M_Matrix_Y2<-subset(HC_M_Matrix, Year==2)
HC_Inverts_Community_Y2 <- HC_M_Matrix_Y2[,7:ncol(HC_M_Matrix_Y2)]
any(is.na(HC_Inverts_Community_Y2))
str(HC_Inverts_Community_Y2)
#Separate environmental data
HC_Inverts_Env_Y2<-HC_M_Matrix_Y2[,1:6]
str(HC_Inverts_Env_Y2)
HC_Inverts_Env_Y2$Days_Since_Carcass_Introduction<-as.factor(HC_Inverts_Env_Y2$Days_Since_Carcass_Introduction)

#PERMANOVAs for each year
adonis(HC_Inverts_Community_Y1 ~ Reach*Days_Since_Carcass_Introduction, data=HC_Inverts_Env_Y1, method="bray", permutations=999)
adonis(HC_Inverts_Community_Y2 ~ Reach*Days_Since_Carcass_Introduction, data=HC_Inverts_Env_Y2, method="bray", permutations=999)

#Indicator taxa analysis
#Year 1
HC_inv_y1_indic<-signassoc(HC_Inverts_Community_Y1, cluster=HC_Inverts_Env_Y1$Reach,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_inv_y1_indic[,c(ncol(HC_inv_y1_indic),(ncol(HC_inv_y1_indic)-1))]
#Indicator analysis found three taxa: Polycentropidae, Brachycentridae, Heptageniidae
#Year 2
HC_inv_y2_indic<-signassoc(HC_Inverts_Community_Y2, cluster=HC_Inverts_Env_Y2$Reach,  mode=0, alternative = "two.sided",control = how(nperm=999))
HC_inv_y2_indic[,c(ncol(HC_inv_y2_indic),(ncol(HC_inv_y2_indic)-1))]
#Indicator analysis found two taxa: Elmidae and Brachycentridae

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

#NMDS plot for stream_site Salmon vs Control
ordiplot(HC_Macroinvertebrate_NMDS, type="n", main="Invertebrate")
with(HC_Macroinvertebrate_NMDS, points(HC_Macroinvertebrate_NMDS, display="sites", col=two_col_vec[stream_site], pch=19, pt.bg=two_col_vec))
with(HC_Macroinvertebrate_NMDS, legend("topleft", legend=levels(stream_site), bty="n", col=two_col_vec, pch=19, pt.bg=two_col_vec))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, stream_site, kind="se", conf=0.95, lwd=2, col="black", show.groups = "Salmon"))
with(HC_Macroinvertebrate_NMDS, ordiellipse(HC_Macroinvertebrate_NMDS, stream_site, kind="se", conf=0.95, lwd=2, col="bisque3", show.groups = "Control"))

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

#############
#Now on to populations rather than communities

#find top taxa
totals<-rbind(HC_Inverts_Community, colSums(HC_Inverts_Community))
totals<-totals[-c(1:58),]
rowSums(totals)
sort(totals,decreasing=TRUE)[1:6]

#Make line plot for Baetids
HC_M_Matrix$RowSums<-rowSums(HC_M_Matrix[,7:ncol(HC_M_Matrix)])
HC_M_Matrix$Treatment<-HC_M_Matrix$Reach
HC_M_Matrix$Baetidae<-(HC_M_Matrix$"Hexapoda Ephemeroptera Baetidae"/HC_M_Matrix$RowSums)
Baetidaeplot<-data.frame(HC_M_Matrix[,c(56,6,57)])
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

#Make line plot for Simuliids
HC_M_Matrix$Simuliidae<-(HC_M_Matrix$"Hexapoda Diptera Simuliidae"/HC_M_Matrix$RowSums)
Simuliidaeplot<-data.frame(HC_M_Matrix[,c(56,6,58)])
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

#Make line plot for Polycentropidae
HC_M_Matrix$Cyrnellusfraternus<-(HC_M_Matrix$"Hexapoda Trichoptera Polycentropidae"/HC_M_Matrix$RowSums)
Cfplot<-data.frame(HC_M_Matrix[,c(6,56,59)])
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
Brachyplot<-data.frame(HC_M_Matrix[,c(56,6,60)])
Brachyplot$Days_Since_Study_Start<-as.numeric(Brachyplot$Days_Since_Study_Start)
bra<-summarySE(Brachyplot, measurevar="Brachycentridae", groupvars=c("Days_Since_Study_Start", "Treatment"))
ggplot(bra, aes(x=Days_Since_Study_Start, y=Brachycentridae, linetype=Treatment)) + 
  geom_errorbar(aes(ymin=Brachycentridae-se, ymax=Brachycentridae+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days") +
  ylab("Brachycentrus relative abundance +/- SE") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")

#Make line plot for Heptageniidae
HC_M_Matrix$Heptageniidae<-(HC_M_Matrix$"Hexapoda Ephemeroptera Heptageniidae"/HC_M_Matrix$RowSums)
Heptageniidplot<-data.frame(HC_M_Matrix[,c(56,6,61)])
Heptageniidplot$Days_Since_Study_Start<-as.numeric(Heptageniidplot$Days_Since_Study_Start)
hep<-summarySE(Heptageniidplot, measurevar="Heptageniidae", groupvars=c("Days_Since_Study_Start", "Treatment"))
ggplot(hep, aes(x=Days_Since_Study_Start, y=Heptageniidae, linetype=Treatment)) + 
  geom_errorbar(aes(ymin=Heptageniidae-se, ymax=Heptageniidae+se), width=.1) +
  geom_line(size=1.5) +
  geom_point(size=2) +
  xlab("Days since study start") +
  ylab("Heptagenia flavescens relative abundance") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20,margin=margin(40,0,0,0)),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16)) +
  scale_x_continuous(breaks=seq(0,700,100)) + 
  geom_vline(xintercept = c(30,400), linetype = "dotted", colour = "black")

#Make line plot for Elmidae
HC_M_Matrix$Elmidae<-(HC_M_Matrix$"Hexapoda Coleoptera Elmidae"/HC_M_Matrix$RowSums)
Elmidplot<-data.frame(HC_M_Matrix[,c(56,6,62)])
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
