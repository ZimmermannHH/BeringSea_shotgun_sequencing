#compare pelagic-benthic ratio

pel<-read.table("2_outclean_pelagic_planktonic_2022_07_15_new.txt",header=TRUE, sep="\t", stringsAsFactors = FALSE)
ben<-read.table("2_outclean_benthics_2022_07_15_new.txt",header=TRUE, sep="\t", stringsAsFactors = FALSE)
names(pel) <-gsub("X","",names(pel))
names(ben) <-gsub("X","",names(ben))
str(pel.sum)
pel.sum=as.data.frame(colSums(pel))
ben.sum=as.data.frame(colSums(ben))
pel.ben.ratio=pel.sum$`colSums(pel)`/ben.sum$`colSums(ben)`
dim(pel.sum)

###
pel.ben<-read.table("pel_ben_sums2.txt",header=TRUE, sep="\t", stringsAsFactors = FALSE)
int<-read.table("Interpolations_2022.txt",header=TRUE, sep="\t", stringsAsFactors = FALSE)
int[1]=NULL

p_b_int=cbind(pel.sum,ben.sum, pel.ben.ratio, pel.ben$age, int)
names(p_b_int)[1]="pelagic.sum"
names(p_b_int)[2]="benthic.sum"
names(p_b_int)[3]="pelagic.benthic.ratio"
names(p_b_int)[4]="age"
p_b_int$age=as.numeric(row.names(p_b_int))

#pelagic:benthic ratio
################################

ratio2=(p_b_int$pelagic.sum-p_b_int$benthic.sum)/(p_b_int$pelagic.sum+p_b_int$benthic.sum)
par(mfrow=c(2,1), mar=c(4,4,1,1))
plot(pel.ben.ratio~p_b_int$age, type="l",col="red",
      xlab="Age (cal kyr BP)", ylab="pelagic:benthic ratio")
plot(p_b_int$pelagic.sum~p_b_int$age, type="l",
     xlab="Age (cal kyr BP)", ylab="counts", 
     ylim=c(0,35000))
lines(p_b_int$benthic.sum~p_b_int$age, type="l",
      xlab="Age (cal kyr BP)", ylab="", col="cornflowerblue")
legend("topright", c("pelagic", "benthic") ,lty=c(1,1), bty="n", col=c("black", "cornflowerblue"))

#ratio phototrophic bacteria:phototrophic protists
########################################
library(Hmisc)
library(dplyr)
#load data
list.all.phot.prot<-read.table("list.all.phot.prot.txt",header=FALSE, stringsAsFactors = FALSE)
list.phot_bacteria<-read.table("list.phot_bacteria.txt",header=FALSE, stringsAsFactors = FALSE)
res.plankt_pel<-read.table("resampled_plankton_pelagic_famnumber_sampleeffort6593_aggregated_pcainput_new2022_07_15.txt",header=TRUE, sep="\t", stringsAsFactors = FALSE)
res.plankt_pel.t=as.data.frame(t(res.plankt_pel))

fam.filter2=res.plankt_pel.t
fam.filter2$fam.names=row.names(res.plankt_pel.t)

List.phot_bacteria=list.phot_bacteria$V1
phot_bacteria.filter=fam.filter2 %>% filter_at(vars(fam.names), any_vars(. %in% List.phot_bacteria))
phot_bacteria.filter=as.data.frame(phot_bacteria.filter)
row.names(phot_bacteria.filter)=phot_bacteria.filter$fam.names
phot_bacteria.filter[26]=NULL 

phyto=subset(fam.filter2, (row.names(fam.filter2) %in% list.all.phot.prot$V1))
phyto[26]=NULL

phyto.sum=as.data.frame(colSums(phyto))
bacpl.sum=as.data.frame(colSums(phot_bacteria.filter))#phot.bacteria.fwrm or phot_bacteria.filter
small_large_ratio_df=as.data.frame(bacpl.sum$`colSums(phot_bacteria.filter)`/phyto.sum$`colSums(phyto)`)
small_large_ratio_df2=as.data.frame(cbind(small_large_ratio_df, row.names(phyto.sum)))
names(small_large_ratio_df2)[1]="small_large_ratio"
names(small_large_ratio_df2)[2]="age"

par(mfrow=c(1,1), mar=c(4,4,1,1))
x=seq(0,20,2)
plot(small_large_ratio_df2$small_large_ratio~small_large_ratio_df2$age, type="l", xaxt="n",
     xlab="Age (cal kyr BP)", ylab="phot. bacteria:phot. protists", ylim=c(0,3))
axis(side=1, at=x, labels = TRUE)

