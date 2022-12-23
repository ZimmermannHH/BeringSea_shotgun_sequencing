################################################################################
#R-script 
#authors: Heike Zimmermann, Stefan Kruse https://github.com/StefanKruse/R_Rarefaction
#aim: resampling of benthic taxa 
################################################################################
library(vegan)

#read in datasets

#plankton and pelagic
library(analogue)
plank.pel<-read.table("2_outclean_plankton.pelagic_2020_10_17.txt",header=TRUE, sep="\t",stringsAsFactors = FALSE)
plank.pel[13]=NULL
names(plank.pel) <-gsub("X","",names(plank.pel))
plank.pel=plank.pel[,order(as.numeric(names(plank.pel)))]

rarecurve(t(plank.pel), sample=min(colSums(plank.pel)), xlab="counts", ylab="families")


plank.pel.ct=chooseTaxa(t(plank.pel), max.abun = 10, noccur=3)
plank.pel.ct.t=as.data.frame(t(plank.pel.ct))
########################

#resampling

specseq=plank.pel.ct.t2


str(specseq)# check data is of required format, count data as integer/double, taxa names as strings
specseq[is.na(specseq)]=0# make sure additional rows or NA values are converted to zeros

specseq$scientific_name<-row.names(specseq)
names(specseq)
# define columns that contain raw count data
names(specseq)
colstart=1
colend=25
names(specseq)[colstart:colend]# check
COLUMNNAMESAREYEARS=TRUE

# species must have unique sample names and need to be merged with the family (for technical reasons)
SPECIESNAMECOLUMN=26# position of the species assignments column in the data frame
names(specseq)[SPECIESNAMECOLUMN]="scientific_name"# change species name column to "scientific_name"

FAMILYNAMEPRESENT=FALSE
#FAMILYNAMECOLUMN=2# position of the family assignments column in the data frame
if(FAMILYNAMEPRESENT)
{
  names(specseq)[FAMILYNAMECOLUMN]="family_name"# change family name column to "family_name"
} else
{
  specseq$family_name="NoFamilyName"# add empty family name column
}

#specseq_final_name=paste(specseq$family_name, make.unique(specseq$scientific_name))# make sure to have individual names for each species/taxa entry
specseq_final_name=paste(gsub(" ", "_", specseq$family_name), make.unique(gsub(" ", "_", specseq$scientific_name)))# make sure to have individual names for each species/taxa entry; to avoid issues with multiple present space character in the columns those are replaces by lines

# define method to use, either "replace" for simple resampling based on weights on the occurance of counts or "birksline" for the calculation of the expected number of pollen types following as stated in Birks, H. J. B., & Line, J. M. (1992). The use of Rarefaction Analysis for Estimating Palynological Richness from Quaternary Pollen-Analytical Data. The Holocene, 2(1), 1–10. doi:10.1177/095968369200200101
samplingmethod="replace" 

(nsampleff=sort(apply(specseq[,c(colstart:colend)], 2, sum))[1])# determine min. read counts for rarefaction, here automatic procedure to find the minimum within the data set

### rarefication
# resample loop for each sample/year present in the data table
nsampleff=sort(apply(specseq[,c(colstart:colend)], 2, sum))[1]# determine min. read counts for rarefaction, here automatic procedure to find the minimum within the data set
#5.6
#6172
resamplingnumber=500# set here the number of resamplings, standard==100
genrare=list()
yr=NULL# in case of column/sample names are numbers corresponding to depths, an X was given in front of the number on data table import, in the loop this X is deleted and the year is assigned as the name of the entry in the filled loop
missing=0# this counter is needed to reorganize the list in case samples without any reads are present. These will be skipped (should be not the case otherwise the minimal counts would be zero what will produce in an error so please check before running the script if you have such samples in your data set end exclude them prior analyses).

for(yrcoli in colstart:colend)
{
  print(yrcoli)
  
  allspec=specseq_final_name[which(specseq[,yrcoli]>0)]
  allspec_counts=specseq[which(specseq[,yrcoli]>0),yrcoli]
  
  if(samplingmethod=="birksline"){
    allspecvector=NULL
    for(ni in 1:length(allspec_counts))
      allspecvector=c(allspecvector,rep(allspec[ni],allspec_counts[ni]))
  }
  
  if(length(allspec)>0)
  {
    if(COLUMNNAMESAREYEARS)
    {# convert to number
      yr=c(yr, as.numeric(gsub(",", ".", gsub("X","",names(specseq)[yrcoli]))))
    } else
    {# take the name as it is
      yr=c(yr, names(specseq)[yrcoli])
    }
    
    sampleeffort=list()
    for(nsampleeffi in nsampleff)# would also be possible to process many differrent sampling efforts with the same loop e.g. replace here in the loop call nsampleeff by c(nsampleeff, 100, 1000, ...) but then the script needs to be adapted to take them in the analyses into account
    {
      repeatsample=list()
      for(repi in 1:resamplingnumber)
      {
        if(samplingmethod=="replace"){
          repeatsample[[repi]]=sample(allspec,nsampleeffi,replace=TRUE, prob=allspec_counts/sum(allspec_counts))# weighted resampling
        }else if(samplingmethod=="birksline"){
          repeatsample[[repi]]=sample(allspecvector,nsampleeffi,replace=FALSE)# Birks & Line without replacement
        }	
      }
      sampleeffort[[which(nsampleff==nsampleeffi)]]=repeatsample
    }
    genrare[[yrcoli-missing-(colstart-1)]]=list(allspec,sampleeffort)
  } else
  {
    missing=missing+1
  }
}
names(genrare)=yr

# processing of the resampled data set
# count total species/family number and reads of individual families
famorig=specseq$family_name
familylevels=names(rev(sort(table(famorig))))

totspec=NULL# number of species per sample
totfam=NULL# number of species per family per sample
for(li in 1:length(genrare))
{
  for(li2 in 1:length(genrare[[li]][[2]]))
  {
    print(paste0(li," - ",li2))
    
    spectot=NULL
    spectot4fam=NULL
    for(repi in 1:100)
    {
      pei=unique(genrare[[li]][[2]][[li2]][[repi]])
      spectot=c(spectot,length(pei))
      spectot4fam=rbind(spectot4fam,table(factor(unlist(lapply(strsplit(split=" ",pei),function(x)return(x[1]))),levels=familylevels)))
    }
    totspec=rbind(totspec, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),Nspecies=spectot))
    totfam=rbind(totfam, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),spectot4fam))
  }
}
str(totspec)
str(totfam)

# simple plots of the processed data
# modify column names ifnot years that can be coerced to numbers
if(!COLUMNNAMESAREYEARS)
{
  totspec$T=factor(totspec$T, levels=names(genrare))# ordered levels to contain original sorting of columns
  totfam$T=factor(totfam$T, levels=names(genrare))# ordered levels to contain original sorting of columns
}
png(paste0("resampled_plankton_pelagic_famnumber_totalspecies_sampleeffort_new",nsampleff,"_plot.png"), width=480,height=480)
par(mar=c(8,4,3,1),las=2)
with(totspec,plot(Nspecies~T, main="number of species per sample"))
dev.off()

# save processed data
write.table(totspec, paste0("resampled_plankton_pelagic_famnumber_totalspecies_sampleeffort_new",nsampleff,".txt"), row.names=FALSE, sep="\t")	
write.table(totfam, paste0("resampled_plankton_pelagic_famnumber_totalfamilies_sampleeffort_new",nsampleff,".txt"), row.names=FALSE, sep="\t")	



# aggregate data on species level
famorig=specseq$scientific_name
familylevels=names(rev(sort(table(famorig))))

totfam=NULL
for(li in 1:length(genrare))
{
  for(li2 in 1:length(genrare[[li]][[2]]))
  {
    print(paste0(li," - ",li2))
    
    spectot4fam=NULL
    for(repi in 1:100)
    {
      pei=genrare[[li]][[2]][[li2]][[repi]]
      spectot4fam=rbind(spectot4fam,table(factor(unlist(lapply(strsplit(split=" ",pei),function(x)return(paste(x[-1],collapse = " ")))),levels=familylevels)))
    }
    totfam=rbind(totfam, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),spectot4fam))
  }
}
str(totfam)

# modify column names ifnot years that can be coerced to numbers
if(!COLUMNNAMESAREYEARS)
{
  totfam$T=factor(totfam$T, levels=names(genrare))# ordered levels to contain original sorting of columns
}


totfam$T=as.numeric(as.character(totfam$T))

# calculate mean values for each species/taxa
speciesfamiliesdf_totfam=NULL
pdf(paste0("resampled_plankton_pelagic_famnumber",nsampleff,"_aggregated_new.pdf"))
par(mar=c(8,4,3,1),las=2)
for(fami in names(totfam)[3:dim(totfam)[2]])
{
  plot(totfam[,fami]~totfam$T,col=rainbow(length(names(totfam)[3:dim(totfam)[2]]),s=0.6)[which(names(totfam)[3:dim(totfam)[2]]==fami)],type="n",lwd=2,main=fami,ylab="Sequence counts",xlab="Depth (m)")
  
  mn=aggregate(totfam[,fami],list(as.numeric(totfam$T)),mean)
  ti=mn$Group.1
  mn=mn$x
  sd=aggregate(totfam[,fami],list(as.numeric(totfam$T)),sd)$x*1.96
  
  polygon(y=c(mn+sd,rev(mn-sd)), x=c(ti,rev(ti)), col="gray", border=NA)
  lines(mn~ti,col="steelblue2",lwd=3)
  
  if(!COLUMNNAMESAREYEARS)
  {
    ti=levels(totfam$T)
  }
  speciesfamiliesdf_totfam=rbind(speciesfamiliesdf_totfam, data.frame(Species=fami,TBP=ti,Mean=mn,CI95=sd))
}
dev.off()

str(speciesfamiliesdf_totfam)

# save processed data
write.table(speciesfamiliesdf_totfam, paste0("resampled_plankton_pelagic_famnumber_",nsampleff,"_aggregated_new.txt"), sep="\t", row.names=FALSE)


### post processing	
# reformat for pca input ... rows are samples, cols are species/taxa
ordidf=NULL
for(fi in levels(factor(speciesfamiliesdf_totfam$Species)))
{
  numberi=speciesfamiliesdf_totfam[speciesfamiliesdf_totfam$Species==fi, ]$Mean
  ordidf=rbind(ordidf, data.frame(Spec=fi,t(numberi)))
}
names(ordidf)[2:dim(ordidf)[2]]=speciesfamiliesdf_totfam[speciesfamiliesdf_totfam$Species==fi, ]$TBP
str(ordidf)

# remove species/taxa column
row.names(ordidf)=ordidf$Spec
ordidf=ordidf[,-1]
str(ordidf)
head(ordidf)

# exclude samples or species that have no records
rowsumsnotzero=which(apply(ordidf,1,sum)>0)
colsumsnotzero=which(apply(ordidf,2,sum)>0)
ordidf=ordidf[rowsumsnotzero,colsumsnotzero]

# export data
write.table(t(ordidf), paste0("resampled_plankton_pelagic_famnumber_sampleeffort",nsampleff,"_aggregated_pcainput_new.txt"), sep="\t")

# comparison of original and resampled data set
# prepare the data
ordiorigdf=specseq[colstart:colend]
row.names(ordiorigdf)=make.unique(specseq$scientific_name)
rowsumsnotzero=which(apply(ordiorigdf,1,sum)>0)
colsumsnotzero=which(apply(ordiorigdf,2,sum)>0)
ordiorigdf=ordiorigdf[rowsumsnotzero,colsumsnotzero]

# species count (counts >= 1)
png(paste0("resampled_plankton_pelagic_famnumber_sampleeffort",nsampleff,"_aggregated_comparison_new.png"), width=480,height=480)
par(mar=c(8,4,3,1),las=2)
barplot(apply(ordiorigdf,2,function(x)length(which(x>=1))), col="tomato", border=FALSE)
barplot(apply(ordidf,2,function(x)length(which(x>=1))), add=TRUE, col="skyblue", border=FALSE)
legend("topleft", c("original","resampled"), col=c("tomato","skyblue"), pch=15, pt.cex=1.5, title="Families >= 1 count")
dev.off()