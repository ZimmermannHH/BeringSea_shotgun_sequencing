################################################################################
#R-script 
#author: Heike Zimmermann
#aim: data preparation on family level before resampling
################################################################################
#load libraries
library(vegan)
library(analogue)
library(rioja)
library(dplyr)
library(readr)
library(tidyverse)
library(stats)
library(BiodiversityR)
library(ape)

#read in the data
##################
krakenout<-read_delim("APMG-689_nt0.2.txt", delim="\t", col_names=F)
names(krakenout) <- c("samples", "percent", "clade_count", "tax_count", "rank", "taxid", "name")

##################

#counts for merged/unmerged 
krakenout.counts <- krakenout %>%
  group_by(samples) %>%
  summarise(sumcount = sum(tax_count))
write.table(krakenout.counts, "krakenout.counts.original.txt", sep="\t")

##################
#take only information on family level 
kraken.fam<-subset(krakenout, rank=="F")
fam=as.data.frame(kraken.fam)
dim(kraken.fam)

#remove unnecessary columns
fam[4:6]=NULL
fam[2]=NULL

#rename header
names(fam)[1]<-"sample"
names(fam)[2]<-"counts"
names(fam)[3]<-"family_name"

#counts for merged/unmerged on family level
fam.counts <- fam %>%
  group_by(sample) %>%
  summarise(sumcount = sum(counts))
write.table(fam.counts, "fam.counts.original.txt", sep="\t")

#remove unnecessary string in sample names
fam$sample <-gsub("_dedup_merged_conf0.2_nt.report","",fam$sample)
fam$sample <-gsub("_dedup_unmerged_conf0.2_nt.report","",fam$sample)

#exclude outlier sample 11.0 kyr BP APMG-6-5 and 8-8 with insufficient reads
fam1.1=subset(fam, sample!="APMG-6-5")
fam2=subset(fam1.1, sample!="APMG-8-8")

fam3 <- fam2 %>%
  group_by(sample, family_name) %>%
  summarise(sumcount = sum(counts))

fam3.df=as.data.frame(fam3)
names(fam3.df)
fam.wide.new<-reshape(fam3.df, idvar = "sample",
                      timevar = "family_name", direction = "wide")   

names(fam.wide.new) <-gsub("sumcount.","",names(fam.wide.new))
#write.table(fam.wide.new, "fam.wide.new.txt", sep="\t")

tail(fam.wide.new)
row.names(fam.wide.new)=fam.wide.new$sample
fam.wide.new[1]=NULL
fam.wide.new[is.na(fam.wide.new)] <- 0
write.table(fam.wide.new, "fam.wide.new.txt", sep="\t")

#transpose data
fam.wide.new.t=as.data.frame(t(fam.wide.new))
names(fam.wide.new.t)
sum(colSums(fam.wide.new.t))
#13120676

#check the blanks
##########################################################
#"APMG-6-7" "APMG-6-8" "APMG-8-7" "APMG-9-4" "APMG-9-12" "APMG-9-13"

blanks=as.data.frame(cbind(fam.wide.new.t[6:7],fam.wide.new.t[17],fam.wide.new.t[22:23],fam.wide.new.t[26]))

#remove family names with rowSum=0
blanks=  blanks[apply(blanks[,-1], 1, function(x) !all(x==0)),]
dim(blanks)
#208  6

#number of read counts in blanks
blank.sums=colSums(blanks)
# APMG-6-7  APMG-6-8  APMG-8-7 APMG-9-12 APMG-9-13  APMG-9-4 
#     235       221       540       225       111       187 
blanks.totalsum=sum(blank.sums)
#1519

#number of families in blanks
blank.spec=specnumber(t(blanks))
# APMG-6-7  APMG-6-8  APMG-8-7 APMG-9-12 APMG-9-13  APMG-9-4 
#   35        30       185        43        23        25 

fam.readsums.blanks=sort(rowSums(blanks),decreasing =TRUE)
write.table(fam.readsums.blanks, "fam.readsums.blanks_2022_07_15.txt", sep="\t")

#how many families have only one read counts in blanks
length(which(fam.readsums.blanks == 1))
#98
length(which(fam.readsums.blanks == 2))
#38
blanks.t=t(blanks)
1519/13120676

blanks.t=t(blanks)
ra.blanks=rankabundance(blanks.t)

rankabunplot(ra.blanks, scale='logabun', addit=FALSE, specnames=c(1:7))

par(mfrow=c(1,1),mar=c(8,4,1,1))
barplot(fam.readsums.blanks,ylim=c(0,450), las=2)
write.table(blanks, "blanks_2022_07_15.txt", sep="\t")

barplot(blank.spec, col="lightgrey", ylim=c(0,100),
        ylab="Number of families in blanks", las=2)
barplot(blank.sums, col="lightgrey", ylim=c(0,300),
        ylab="Number of read counts in blanks", las=2)

#assess human contamination
human=subset(blanks, row.names(blanks)=="Hominidae")
human.sum=colSums(human)
# APMG-8-7  APMG-6-8  APMG-6-7 APMG-9-12 APMG-9-13  APMG-9-4 
#      94        91        48        74        54        66 
sum(human.sum)/sum(blank.sums)
#0.281106
prop.human=(human.sum/blank.sums)*100
# APMG-6-7  APMG-6-8  APMG-8-7 APMG-9-12 APMG-9-13  APMG-9-4 
# 40.000000 41.176471  8.888889 32.888889 48.648649 35.294118 

##########################################################


#rename columns according to sample age (cal kyr BP)
##########################################################

names(fam.wide.new.t)[1]="19.9" #"APMG.6.1"
names(fam.wide.new.t)[2]="16.5" #"APMG.6.2"
names(fam.wide.new.t)[3]="14.03" #"APMG.6.3"
names(fam.wide.new.t)[4]="12.61" #"APMG.6.4"
names(fam.wide.new.t)[5]="5.6"   #"APMG.6.6"
names(fam.wide.new.t)[8]="19.30" #"APMG.8.1"
names(fam.wide.new.t)[9]="11.56" #"APMG.8.10"
names(fam.wide.new.t)[10]="11.17" #"APMG.8.11"
names(fam.wide.new.t)[11]="10.73" #"APMG.8.12"

names(fam.wide.new.t)[12]="18.39" #"APMG.8.2"
names(fam.wide.new.t)[13]="17.34" #"APMG.8.3"
names(fam.wide.new.t)[14]="15.91" #"APMG.8.4"
names(fam.wide.new.t)[15]="15.6" #"APMG.8.5"
names(fam.wide.new.t)[16]="13.62" #"APMG.8.6"
names(fam.wide.new.t)[18]="12.13" #"APMG.8.9"

names(fam.wide.new.t)[19]="9.84" #"APMG.9.1"
names(fam.wide.new.t)[20]="1.75" #"APMG.9.10"
names(fam.wide.new.t)[21]="1.08" #"APMG.9.11"

names(fam.wide.new.t)[24]="8.29" #"APMG.9.2"
names(fam.wide.new.t)[25]="7.56" #"APMG.9.3"
names(fam.wide.new.t)[27]="6.26" #"APMG.9.5"
names(fam.wide.new.t)[28]="4.31" #"APMG.9.6"
names(fam.wide.new.t)[29]="3.66" #"APMG.9.7"
names(fam.wide.new.t)[30]="2.99" #"APMG.9.8"
names(fam.wide.new.t)[31]="2.36" #"APMG.9.9"


fam.counts=sort(colSums(fam.wide.new.t))
sum(fam.counts)
#13120676

ori.col.sums=as.data.frame(colSums(fam.wide.new.t))
median(colSums(fam.wide.new.t))
x=seq(0,20,2)
plot(ori.col.sums$`colSums(fam.wide.new.t)`~as.numeric(row.names(ori.col.sums)),
     xlab="Age (cal kyr BP)", ylab="assigned read counts", ylim=c(100000,1700000),
     pch=16)
abline(h=324044, col="grey50", lty=2)
rarecurve(fam.wide.new.t, sample=min(rowSums(fam.wide.new.t)))

#PCoA plots 
##########################################################
names(fam.wide.new.t)
fam.wide.new=as.data.frame(t(fam.wide.new.t))
fam.wide.new.t.bray <- vegdist(fam.wide.new, method = "bray")
str(fam.wide.new.t)
par(mfrow=c(1,1), mar=c(4,4,1,1))
pcoa.ori=pcoa(fam.wide.new.t.bray, correction="none", rn=NULL)
biplot(pcoa.ori, Y=NULL, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1, cex=.6)

#filtering
##########################################################
#remove blanks
fam.wide.new.t[26]=NULL
fam.wide.new.t[22:23]=NULL
fam.wide.new.t[17]=NULL
fam.wide.new.t[6:7]=NULL

#order colnames according to age an calculate sum of counts per sample
fam.wide.new.t2=fam.wide.new.t[,order(as.numeric(names(fam.wide.new.t)))]
fam.wide.new.t.cs=colSums(fam.wide.new.t2)

#remove human contamination 
fam.wide.new.t3<-subset(fam.wide.new.t2, row.names(fam.wide.new.t2)!="Hominidae")
fam.wide.new.t2.cs=colSums(fam.wide.new.t3)

#42090 human read counts removed

#remove other contamination
##########################################################
#add a column called fam.names to filter the data
#col=taxa, row=sample
fam.wide.new.t2$fam.names<-row.names(fam.wide.new.t2)
blanks.list=row.names(blanks)
sort(blanks.list)

#contamination removal of taxa relevant for our data analysis
fam.wide.new.t4<-subset(fam.wide.new.t2, row.names(fam.wide.new.t2)!="Cyprinidae")
fam.wide.new.t5<-subset(fam.wide.new.t4, row.names(fam.wide.new.t4)!="Muridae")
fam.wide.new.t6<-subset(fam.wide.new.t5, row.names(fam.wide.new.t5)!="Echeneidae")

write.table(fam.wide.new.t6, "fam.wide.new.t6.txt", sep="\t", dec=".")
