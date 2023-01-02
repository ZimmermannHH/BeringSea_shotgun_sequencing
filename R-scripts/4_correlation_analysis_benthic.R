################################################################################
#R-script 
#author: Heike Zimmermann
#aim: stratigraphic diagram and correlation analysis
################################################################################
#load libraries
library(vegan)
library(analogue)
library(rioja)
library(psych)

#load data: resampled dataset
res.ben<-read.table("resampled_benthic_famnumber_sampleeffort1839_aggregated_pcainput.txt",header=TRUE, sep="\t", stringsAsFactors = FALSE)
names(res.ben)
#interpolated SST and IP25 
interpolations<-read.table("Interpolations_2022.txt",header=TRUE, sep="\t", stringsAsFactors = FALSE)
row.names(interpolations)=interpolations$Age_ka
interpolations[1]=NULL

#make subset of sst.int with only the same times-slices as the samples
res.ct.int = subset(res.ben, (row.names(res.ben) %in% row.names(interpolations)))
res.ct.int.sst=cbind(interpolations, res.ct.int)
res.ct.int.ma=as.matrix(res.ct.int.sst)

ben.fam.sums=as.data.frame(sort(colSums(res.ben)))
age=as.numeric(row.names(res.ben))

sort(colSums(res.ben))

#stratigraphic diagrams
########################################################################
ben.dsqrt=sqrt(sqrt(res.ben))
ben.dsqrt.dist<-vegdist(ben.dsqrt, method="bray", binary=FALSE)
ben.dsqrt.clust<-chclust(ben.dsqrt.dist, method="coniss")

par(mfrow=c(1,2))
ben.dsqrt.bst<-bstick(ben.dsqrt.clust) 
plot(ben.dsqrt.clust, hang=-1, cex=0.7, main="CONISS")

#calculate proportions
all.b=res.ben
all.anteile.b=all.b
for (i in 1:dim(all.b)[1])
{
  for (j in 1:dim(all.b)[2])
  {
    all.anteile.b[i,j]=all.b[i,j]/rowSums(all.b)[i]
  }
}

all.b.a=all.anteile.b
all.b.a2<-as.data.frame(all.b.a)
#make proportions to percentage
all.b.a3=all.b.a2*100

#proportions
res.ben.ct=chooseTaxa(all.b.a3, max.abun = 1.5, noccur=5, type="AND")
write.table(res.ben.ct, "input.stratplot.benthic.txt", sep="\t")

windows(width=500,  height=150) 
par(oma=c(0.25,0.25,0.25,0.25))
x<- strat.plot(res.ben.ct,  yvar  =  age,  y.rev=TRUE,
               scale.percent=FALSE,
               las=2,
               ylabel = "Age (cal kyr BP)",
               cex.axis = 0.8,
               #srt.xlabel=45,
               yBottom=0.07, xLeft=0.05, yTop=0.80, xRight=1, #plotting area margins, relative units, xLeft may require adjustment to make room for y-axes; srt.xlabel=45, #angle of species names
               cex.ylabel=1, cex.xlabel=1, cex.title=1,
               xSpace=0.01, #space between graphs in relative units
               wa.order="none",
               scale.minmax=FALSE, #show only min and max values on x-axes to avoid crowding
               y.axis=TRUE, y.tks=c(0,2,4,6,8,10,12,14,16,18,20), ylim=c(0,20),plot.poly=TRUE,col.poly="peachpuff4",
               plot.line=FALSE, plot.bar=TRUE,
               clust=ben.dsqrt.clust, clust.width=0.05,lwd.bar=2,
               col.bar="black",
               add=TRUE)

#Correlation analysis
########################################################################
benthic.ma=as.matrix(res.ben)

#BH: Benjamini-Hochberg p-value correction
ben.corr=corr.test(benthic.ma, y = NULL, use = "pairwise",method="spearman",
                   adjust="BH",
                   alpha=.1,ci=TRUE,minlength=5)
ben.corr.coeff=ben.corr$r
ben.corr.p=ben.corr$p

#above diagonal are adjusted p-values, below are original p-values
#set lower part of matrix to 0
ben.corr.coeff[ lower.tri(ben.corr.coeff, diag=TRUE) ]<- 0
ben.corr.p[ lower.tri(ben.corr.p, diag=TRUE) ]<- 1

write.table(ben.corr.coeff, "ben.corr.coeff_spear_adjp0.1_uppertriangle_2022_07_19.txt", sep ="\t")
write.table(ben.corr.coeff, "ben.corr.p_spear_adjp0.1_uppertriangle_2022_07_19.txt", sep ="\t")

########################################################################
#bring correlation matrix into long format
#from remotes::install_github("heuselm/mocode")
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
########################################################################
#filter for positive correlations above the threshold of the coefficient 
#rho>0.4 and p<0.1

flat.matrix.ben=flattenCorrMatrix(ben.corr.coeff, ben.corr.p)
flat.matrix.ben.s<-subset(flat.matrix.ben, p < 0.1)
flat.matrix.ben.s.cor<-subset(flat.matrix.ben.s, cor > 0.4)
# row        column        cor       p
# Aplysiidae Oxystominidae 0.7379585 0.05981511

#as only one significant correlation is left, no network was computed 

##############################################################
#pairwise correlations including SST and IP25 with 
#Benjamini-Hochberg p-value correction

set.seed(21)
ben.corr.sst=corr.test(res.ct.int.ma, y = NULL, use = "pairwise",method="spearman",
                       adjust="BH",
                       alpha=.1,ci=TRUE,minlength=5)
ben.corr.sst.coeff.sst=ben.corr.sst$r
ben.corr.sst.p=ben.corr.sst$p

ben.corr.sst.coeff.sst[ lower.tri(ben.corr.sst.coeff.sst, diag=TRUE) ]<- 0
ben.corr.sst.p[ lower.tri(ben.corr.sst.p, diag=TRUE) ]<- 1

#filter for correlations above the threshold of the coefficient 
#rho>0.4 and p<0.1
flat.matrix.ben.sst=flattenCorrMatrix(ben.corr.sst.coeff.sst, ben.corr.sst.p)
flat.matrix.ben.sst.s<-subset(flat.matrix.ben.sst, p < 0.1)
flat.matrix.ben.sst.s.cor<-subset(flat.matrix.ben.sst.s, abs(cor) > 0.4)

#filter which taxa are correlated with SSTs or IP25
flat.matrix.ben.sst.s.sst<-subset(flat.matrix.ben.sst.s.cor, row=="SST_interpol" | column=="SST_interpol")
#          row       column        cor            p
#SST_interpol  Zosteraceae  0.7319346 0.0265436304

flat.matrix.ben.sst.s.ip25<-subset(flat.matrix.ben.sst.s.cor, row=="IP25_interpol" | column=="IP25_interpol")
#          row       column        cor            p
#IP25_interpol  Zosteraceae -0.6854394 0.0969532461
