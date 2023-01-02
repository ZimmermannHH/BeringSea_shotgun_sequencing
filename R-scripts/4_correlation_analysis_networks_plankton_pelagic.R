################################################################################
#R-script 
#author: Heike Zimmermann
#aim: stratigraphic diagram, correlation analysis, network analysis (Spearman and ecoCopula)
################################################################################
#load libraries
library(vegan)
library(analogue)
library(rioja)
library(reshape2)
library(psych)
library(Hmisc)
library(igraph)
library(ecoCopula)
library(data.table)
library(car)

#load data: resampled dataset
res.plankt_pel<-read.table("resampled_plankton_pelagic_famnumber_sampleeffort6593_aggregated_pcainput_new2022_07_15.txt",header=TRUE, sep="\t", stringsAsFactors = FALSE)
res.plankt_pel.t=as.data.frame(t(res.plankt_pel))
dim(res.plankt_pel)
#167 families
 
fam.sums=as.data.frame(rowSums(res.plankt_pel.t))
age=as.numeric(row.names(res.plankt_pel))

#counts per family
families.counts.df=as.data.frame(rowSums(res.plankt_pel.t))
sum(families.counts.df)
#164825

#interpolated SST, IP25
interpolations<-read.table("Interpolations_2022.txt",header=TRUE, sep="\t", stringsAsFactors = FALSE)
row.names(interpolations)=interpolations$Age_ka
interpolations[1]=NULL
#make subset of interpolations with only the same times-slices as the samples
res.ct.int = subset(res.plankt_pel, (row.names(res.plankt_pel) %in% row.names(interpolations)))

#functional annotations
fam_function=read.table("Taxa_habitat_normalized_counts_SO201-2-12KL_22_07_23.txt", sep="\t", header=TRUE,stringsAsFactors = TRUE)

###############

dim(res.plankt_pel)
#167
names(res.plankt_pel.ct)

sum(res.plankt_pel) #164825
sum(res.plankt_pel.ct) #122828

124975-122828 #2147
sort(colSums(res.plankt_pel))

#################

pel.dsqrt=sqrt(sqrt(res.plankt_pel))
pel.dsqrt.dist<-vegdist(pel.dsqrt, method="bray", binary=FALSE)
pel.dsqrt.clust<-chclust(pel.dsqrt.dist, method="coniss")

par(mfrow=c(1,2))
pel.dsqrt.bst<-bstick(pel.dsqrt.clust) 
plot(pel.dsqrt.clust, hang=-1, cex=0.7, main="CONISS")


#stratigraphic diagrams
######################################
##calculate proportions
all=res.plankt_pel
all.anteile=all
for (i in 1:dim(all)[1])
{
  for (j in 1:dim(all)[2])
  {
    all.anteile[i,j]=all[i,j]/rowSums(all)[i]
  }
}

allpp.a=all.anteile
allpp.a2<-as.data.frame(allpp.a)
allpp.a3=allpp.a2*100 #make percentages

res.plankt_pel.ct=chooseTaxa(allpp.a3, max.abun = 1.25, noccur=5, type="AND")
dim(res.plankt_pel.ct)
#34 families
#add Calanidae and Metridinidae to stratigraphic diagram
res.plankt_pel.ct$Calanidae=allpp.a3$Calanidae
res.plankt_pel.ct$Metridinidae=allpp.a3$Metridinidae
write.table(res.plankt_pel.ct, "input.pel.stratplot_2.txt", sep="\t")

names(res.plankt_pel.ct)
windows(width=500,  height=150)
par(oma=c(0.25,0.25,0.25,0.25))
x<- strat.plot(res.plankt_pel.ct,  yvar  =  age,  y.rev=TRUE,
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
               y.axis=TRUE, y.tks=c(0,2,4,6,8,10,12,14,16,18,20), ylim=c(0,20),plot.poly=TRUE,col.poly="palegreen3",
               plot.line=FALSE, plot.bar=TRUE,
               clust=pel.dsqrt.clust, clust.width=0.05,lwd.bar=2,
               col.bar="black",
               add=TRUE)

#Correlation analysis
######################################

#correlation for network computation
set.seed(21)
pel.pp.ma=as.matrix(res.plankt_pel)

#Spearman correlations; BH: Benjamini-Hochberg p-value correction
pel.pp.corr=corr.test(pel.pp.ma, y = NULL, use = "pairwise",method="spearman",
                      adjust="BH",
                      alpha=.1,ci=TRUE,minlength=5)
pel.pp.corr.coeff=pel.pp.corr$r
pel.pp.corr.p=pel.pp.corr$p

#above diagonal are adjusted p-values, below are original p-values
#set lower part of matrix to 0
pel.pp.corr.coeff[ lower.tri(pel.pp.corr.coeff, diag=TRUE) ]<- 0
pel.pp.corr.p[ lower.tri(pel.pp.corr.p, diag=TRUE) ]<- 1

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

#filter for positive correlations above the threshold of the coefficient 
#rho>0.4 and p<0.1
flat.matrix.pp=flattenCorrMatrix(pel.pp.corr.coeff, pel.pp.corr.p)
flat.matrix.pp.s<-subset(flat.matrix.pp, p < 0.1)
flat.matrix.pp.s.cor<-subset(flat.matrix.pp.s, cor > 0.4)
write.table(flat.matrix.pp.s.cor, "flat.matrix.pp.s.cor_2022_05_19.txt", sep ="\t")


#correlation analysis to detect correlations of families
#with SSTs or IP25
##############################################################

#combine composition table with env. variables and convert to matrix
pp.env=as.data.frame(cbind(interpolations,res.ct.int))
pp.env.ma=as.matrix(pp.env)

#pairwise correlations including SST and IP25 with BH = Benjamini-Hochberg correction for false discovery rate
set.seed(21)
pp.corr.env=corr.test(pp.env.ma, y = NULL, use = "pairwise",method="spearman",
                      adjust="BH",
                      alpha=.2,ci=TRUE,minlength=5)
pp.corr.env.coeff=pp.corr.env$r
pp.corr.env.p=pp.corr.env$p
write.table(pp.corr.env.coeff, "pp.corr.new_interpolations.coeff_alpha0.2_2022_10_21.txt", sep ="\t")
write.table(pp.corr.env.p, "pp.corr.new_interpolations.pval_alpha0.2_2022_10_21.txt", sep ="\t")

#the lower triangle contains original p-values, the upper one the adjusted
#for further calculations, we only want the adjusted p-values and set the lower triangle to zero
pp.corr.env.coeff[ lower.tri(pp.corr.env.coeff, diag=TRUE) ]<- 0
pp.corr.env.p[ lower.tri(pp.corr.env.p, diag=TRUE) ]<- 1

flat.matrix.pp.env=flattenCorrMatrix(pp.corr.env.coeff, pp.corr.env.p)
flat.matrix.pp.env.s<-subset(flat.matrix.pp.env, p < 0.05) #additional setting for supplement
flat.matrix.pp.env.s<-subset(flat.matrix.pp.env, p < 0.1)  #setting for manuscript
flat.matrix.pp.env.s<-subset(flat.matrix.pp.env, p < 0.2)  #additional setting for supplement
flat.matrix.pp.env.s.cor<-subset(flat.matrix.pp.env.s, cor > 0.4)
#flat.matrix.pp.env.s.cor<-subset(flat.matrix.pp.env.s, abs(cor) > 0.4) #setting to include negative correlations

flat.matrix.pp.env.s.sst<-subset(flat.matrix.pp.env.s.cor, row=="SST_interpol" | column=="SST_interpol")
#       row             column        cor       p
# 173   SST_interpol       Blenniidae 0.6084615 0.044915470
# 380   SST_interpol Chaetocerotaceae 0.6476923 0.025175899
# 498   SST_interpol    Chlorellaceae 0.6328140 0.031764235
# 530   SST_interpol    Chlorobiaceae 0.5576923 0.081395165
# 3488  SST_interpol       Hyellaceae 0.6483160 0.025122135
# 3572  SST_interpol Klebsormidiaceae 0.5954417 0.053412223
# 4562  SST_interpol   Microcystaceae 0.5776923 0.065857362
# 4853  SST_interpol       Mustelidae 0.7346154 0.004726573
# 9182  SST_interpol   Roseiflexaceae 0.6769231 0.015781813
# 9318  SST_interpol       Salmonidae 0.6184615 0.039319030
# 9732  SST_interpol       Sciaenidae 0.5838462 0.060540471
# 10155 SST_interpol       Sebastidae 0.6178111 0.039750271
# 11630 SST_interpol     Syngnathidae 0.6576923 0.021793585

flat.matrix.pp.env.s.ip25<-subset(flat.matrix.pp.env.s.cor, row=="IP25_interpol" | column=="IP25_interpol")
#       row           column         cor        p
# 92    IP25_interpol Bacillariaceae 0.6860444 0.013339000
# 121   IP25_interpol Bathycoccaceae 0.6913012 0.011995270
# 1954  IP25_interpol        Gadidae 0.7197847 0.006685856
# 4372  IP25_interpol   Metridinidae 0.6221109 0.036842136
# 5996  IP25_interpol   Oxytrichidae 0.7566423 0.002788882
# 6904  IP25_interpol Phaeocystaceae 0.5881449 0.057216952
# 11326 IP25_interpol    Suessiaceae 0.5985177 0.051675669
# 13531 IP25_interpol       Ulvaceae 0.6615296 0.020830579

pp.sst.ip25<-as.data.frame(rbind(flat.matrix.pp.env.s.sst,flat.matrix.pp.env.s.ip25))
write.table(pp.sst.ip25, "newinterpolations_PP_pwcorr_spm_bh_abs0.4_adjp0.2_2022_07_24.txt", sep="\t", dec=".")


#Network analysis
##############################################
#to make an adjacency matrix, the filtered correlations need to be transformed
#back from long format to matrix
flat.matrix.pp.s.ma=acast(flat.matrix.pp.s.cor, row~column, value.var="cor")
flat.matrix.pp.s.ma[is.na(flat.matrix.pp.s.ma)] <- 0

#make square matrix
un1 <- unique(sort(c(colnames(flat.matrix.pp.s.ma), rownames(flat.matrix.pp.s.ma))))
m2 <- matrix(NA, length(un1), length(un1), dimnames = list(un1, un1))
m2[row.names(flat.matrix.pp.s.ma), colnames(flat.matrix.pp.s.ma)] <- flat.matrix.pp.s.ma
m2[is.na(m2)] <- 0

write.table(m2, "spearm_bh_adjp0.01_rho0.4_2022_07_15.txt", sep="\t", dec=".")


#global network computation
################################

g <- graph_from_adjacency_matrix(m2, mode = "undirected", weighted = TRUE, diag = TRUE)
isolated.nodes = which(igraph::degree(g)==0)
#remove isolated nodes (=nodes without edges)
g2=igraph::delete.vertices(g, isolated.nodes)
gsize(g2)
#446 edes in graph
gorder(g2)
#148 nodes in  graph

#calculate node degrees for each family
node.degrees=as.data.frame(igraph::degree(g2))

write.table(node.degrees, "node.degrees_ppfam.txt", sep="\t")

##############################################
#scale node size according to abundance (log) in dataset
res.pp.colSums=as.data.frame(colSums(res.plankt_pel))
names(res.pp.colSums)[1]<-"V1"
taxa.list.size = subset(res.pp.colSums, (row.names(res.pp.colSums) %in% row.names(m2)))

V(g2)$size <- log1p(taxa.list.size$V1)
V(g2)$label <-row.names(taxa.list.size)


#colors for groups: read-in script 2_taxa_prep 
############################################################
####colors pelagic

list.chlorophyta=read.table("list.chlorophyta", header=FALSE)
list.charophyta=read.table("list.charophyta", header=FALSE)
V(g2)$color[names(V(g2))%in% list.all.phot.prot$V1] <- "seagreen3"
V(g2)$color[names(V(g2))%in% list.chlorophyta$V1] <- "lightgreen"
V(g2)$color[names(V(g2))%in% list.charophyta$V1] <- "lightgreen"
V(g2)$color[names(V(g2))%in% list.phot_bacteria$V1] <- "turquoise2"
V(g2)$color[names(V(g2))%in% list.all.het.prot$V1] <- "peachpuff3"
V(g2)$color[names(V(g2))%in% list.fish$V1] <- "orange"
V(g2)$color[names(V(g2))%in% List.chondr] <- "orange"
V(g2)$color[names(V(g2))%in% list.metazoa$V1]<- "steelblue3"
V(g2)$color[names(V(g2))%in% list.cnidaria$V1] <- "rosybrown1"
V(g2)$color[names(V(g2))%in% list.crustaceae$V1] <- "yellow"
V(g2)$color[names(V(g2))%in% list.mollusca$V1] <-"salmon3"

par(bg="white", mar=c(0,0,0,0))
set.seed(21)
plot(g2,
     vertex.label=ifelse(V(g2)$size > 4, V(g2)$label, NA),
     vertex.label.cex=0.5,
     vertex.label.color="Black",
     vertex.label.family="sans",
     vertex.frame.color="white",
     edge.color="lightgrey",
     layout=layout.fruchterman.reingold 
)

legend("topleft",
       legend=c("green algae", "phototrophic protists", "phototrophic bacteria",
                "heterotrophic protists", "Crustaceae", "Cnidaria",
                "Fish, Chondrychthyae", "Marine mammals",
                "Mollusca"
       )  ,
       col = c("lightgreen", "seagreen3", "turquoise2",
               "peachpuff3", "yellow", "rosybrown1",
               "orange", "steelblue3",
               "salmon3"
               
       ) , bty = "n",
       pt.cex = 2, pch=20 , cex = .8 , horiz = FALSE, inset = c(0, 0))


#for sub figure with nodes colored according to positive correlations 
#with IP25 and SSTs
V(g2)$color[names(V(g2))%in% flat.matrix.pp.env.s.sst$column] <-"tomato3"
V(g2)$color[names(V(g2))%in% flat.matrix.pp.env.s.ip25$column] <-"royalblue3"
V(g2)$color[names(V(g2))%nin% c(flat.matrix.pp.env.s.sst$column,flat.matrix.pp.env.s.ip25$column)]<-"lightgrey"

par(fig=c(.75, 1, .5, .95), new = TRUE, mar=c(0,0,0,0))#c(bottom, left, top, right)
set.seed(21)
plot(g2,#margin=-0.1,
     vertex.label=NA,
     edge.width=0.3,
     edge.color="lightgrey",
     vertex.label.cex=0.5,
     vertex.label.color="Black",
     vertex.label.family="sans",
     vertex.frame.color="white",
     layout=layout.fruchterman.reingold)

legend("bottomright", legend=c("SSTs rho > 0.4", "IP25 rho > 0.4")  , col = c("tomato3","royalblue3") , bty = "n",
       pt.cex = 3, pch=20 , cex = 0.8 , horiz = FALSE, inset = c(0, 0))


#links among families and functional groups
##############################################

#which are the neighbors (families) of nodes that are pos. correlated to IP25/SSTs

#SSTs
sst.corr=which(names(V(g2))%in% flat.matrix.pp.env.s.sst$column)
#18  25  28  29  72  73  83  85 117 118 121 124 134

#IP25
ip25.corr=which(names(V(g2))%in% flat.matrix.pp.env.s.ip25$column)
#13  15  55  81  96 103 132 145


group.sst=V(g2)$name[sst.corr]
group.ip25=V(g2)$name[ip25.corr]

#group correlated with sst:
# [1] "Blenniidae"       "Chaetocerotaceae" "Chlorellaceae"    "Chlorobiaceae"   
# [5] "Hyellaceae"       "Klebsormidiaceae" "Microcystaceae"   "Mustelidae"      
# [9] "Roseiflexaceae"   "Salmonidae"       "Sciaenidae"       "Sebastidae"      
# [13] "Syngnathidae"  

#group correlated with ip25:
# [1] "Bacillariaceae" "Bathycoccaceae" "Gadidae"        "Metridinidae"  
# [5] "Oxytrichidae"   "Phaeocystaceae" "Suessiaceae"    "Ulvaceae" 


#adjacent nodes to those which are correlated with SSTS/IP25
adj.vert.sst=adjacent_vertices(g2, group.sst,mode = "total")
adj.vert.ip25=adjacent_vertices(g2, group.ip25,mode = "total")

#convert igraph object to vector of IDs --> produces list of lists
adj.vert.ssts.list <- sapply(adj.vert.sst, as_ids)
adj.vert.ip25.list <- sapply(adj.vert.ip25, as_ids)

#convert to long format
adj.vert.ssts.long=do.call(rbind, lapply(adj.vert.ssts.list, data.frame))
adj.vert.ip25.long=do.call(rbind, lapply(adj.vert.ip25.list, data.frame))

write.table(adj.vert.ssts.long, "mat_adjacent_vertices_SSTs_2022_07_24.txt", sep="\t", dec=".")
write.table(adj.vert.ip25.long, "mat_adjacent_vertices_IP25_2022_07_24.txt", sep="\t", dec=".")

#how often is taxon neighbor to SST- or IP25-correlated nodes
ad.sst.tab=as.data.frame(table(adj.vert.ssts.long))
ad.ip25.tab=as.data.frame(table(adj.vert.ip25.long))

ad.sst.tab.sums=sum(ad.sst.tab$Freq)
#114
ad.ip25.tab.sums=sum(ad.ip25.tab$Freq)
#109

write.table(ad.sst.tab, "adjacent_vertices_SSTs_final_2022_07_24.txt", sep="\t", dec=".")
write.table(ad.ip25.tab, "adjacent_vertices_IP25_final_2022_07_24.txt", sep="\t", dec=".")


#Assign functional groups to correlated- and neighbour-families
names(fam_function)
fam_function[6:9]=NULL
fam_function[5]=NULL

#make subset of families with only the same as in ad.ip25.tab$adj.vert.ip25.long
fam.func.SST = subset(fam_function, fam_function$Family.name %in% ad.sst.tab$adj.vert.ssts.long)
fam.func.IP25 = subset(fam_function, fam_function$Family.name %in% ad.ip25.tab$adj.vert.ip25.long)

#sort columns
fam.func.SST.sorted <- with(fam.func.SST,  fam.func.SST[order(Family.name) , ])
fam.func.IP25.sorted <- with(fam.func.IP25,  fam.func.IP25[order(Family.name) , ])
dim(fam.func.SST.sorted)
dim(ad.sst.tab)

#combine
ad.sst.tab.func<-as.data.frame(cbind(ad.sst.tab,fam.func.SST.sorted))
ad.ip25.tab.fun<-as.data.frame(cbind(ad.ip25.tab,fam.func.IP25.sorted))
#check if all taxa are correctly assigned to functional group, then remove first col
ad.sst.tab.func[1]=NULL
ad.ip25.tab.fun[1]=NULL


#add col with variable
ad.sst.tab.func$Variable="SST"
ad.ip25.tab.fun$Variable="IP25"

#combine
ad.fun.all<-as.data.frame(rbind(ad.sst.tab.func,ad.ip25.tab.fun))

#Families
boxplot(ad.fun.all$Freq ~ ad.fun.all$Variable, xlab=" ", ylab="Number of adjacent links")


#number of adjacent vertices per group
sort(table(ad.sst.tab.func$taxonomic.group))
sort(table(ad.ip25.tab.fun$taxonomic.group))

#SST
fish.av.sst=subset(ad.sst.tab.func, taxonomic.group=="Fish")
sum(fish.av.sst$Freq)/sum(ad.sst.tab.func$Freq)
#0.4298246
centricdiatoms.av.sst=subset(ad.sst.tab.func, taxonomic.group=="Diatoms-centric")
sum(centricdiatoms.av.sst$Freq)/sum(ad.sst.tab.func$Freq)
#0.01754386
Chlorophyta.av.sst=subset(ad.sst.tab.func, taxonomic.group=="Chlorophyta")
sum(Chlorophyta.av.sst$Freq)/sum(ad.sst.tab.func$Freq)
#0.1140351
photbact.av.sst=subset(ad.sst.tab.func, habitat.trophic.status=="plankton, phototrophic bacteria")
sum(photbact.av.sst$Freq)/sum(ad.sst.tab.func$Freq)
#0.254386

#IP25
centricdiatoms.av.ip25=subset(ad.ip25.tab.fun, taxonomic.group=="Diatoms-centric")
sum(centricdiatoms.av.ip25$Freq)/sum(ad.ip25.tab.fun$Freq)
#0.1559633
pennatediatoms.av.ip25=subset(ad.ip25.tab.fun, taxonomic.group=="Diatoms-pennate")
sum(pennatediatoms.av.ip25$Freq)/sum(ad.ip25.tab.fun$Freq)
#0.266055
0.1559633+0.266055
#0.4220183
Chlorophyta.av.ip25=subset(ad.ip25.tab.fun, taxonomic.group=="Chlorophyta")
sum(Chlorophyta.av.ip25$Freq)/sum(ad.ip25.tab.fun$Freq)
#0.1192661
Mammals.av.ip25=subset(ad.ip25.tab.fun, taxonomic.group=="Mammals")
sum(Mammals.av.ip25$Freq)/sum(ad.ip25.tab.fun$Freq)
#0.04587156


#Levene Test (homogenous variances?)
str(ad.fun.all)
ad.fun.all$Variable=as.factor(ad.fun.all$Variable)
leveneTest(ad.fun.all$Freq,ad.fun.all$Variable)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value  Pr(>F)  
# group  1  3.2786 0.07409 .
#       77                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#non-significant result suggests homogenous distribution of variances
t.test(ad.fun.all$Freq ~ ad.fun.all$Variable, var.equal = TRUE, alternative = "two.sided") 
# Two Sample t-test
# 
# data:  ad.fun.all$Freq by ad.fun.all$Variable
# t = 4.4415, df = 77, p-value = 2.948e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.9144252 2.4007009
# sample estimates:
#   mean in group IP25  mean in group SST 
# 3.892857           2.235294 


#Functional groups
#aggregate based on taxonomic.group and get links per functional group
fg.sst=aggregate(ad.sst.tab.func$Freq, by=list(ad.sst.tab.func$taxonomic.group), FUN=sum)
fg.IP25=aggregate(ad.ip25.tab.fun$Freq, by=list(ad.ip25.tab.fun$taxonomic.group), FUN=sum)
fg.sst$Variable<-"SST"
fg.IP25$Variable<-"IP25"
fg.all<-as.data.frame(rbind(fg.sst, fg.IP25))

boxplot(fg.all$x ~ fg.all$Variable, xlab=" ", ylab="Number of adjacent links")

#Levene Test (homogenous variances?)

str(fg.all)
fg.all$Variable=as.factor(fg.all$Variable)
leveneTest(fg.all$x,fg.all$Variable)
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  1  0.1174  0.735
#       23     


#non-significance --> set var.eqal=TRUE
t.test(fg.all$x ~ fg.all$Variable, var.equal = TRUE, alternative = "two.sided") 
# Two Sample t-test
# 
# data:  fg.all$x by fg.all$Variable
# t = 0.07092, df = 23, p-value = 0.9441
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -8.847868  9.476073
# sample estimates:
#   mean in group IP25  mean in group SST 
# 9.083333           8.769231 

table(ad.fun.all$Variable)
# IP25  SST 
# 28   51 
table(fg.all$Variable)
# IP25  SST 
# 12   13 
ad.sst.tab.func
ad.ip25.tab.fun

#sum of links
sum(ad.sst.tab.func$Freq)
#114
13/114
#0.1140351
51/114
#0.4473684

sum(ad.ip25.tab.fun$Freq)
#109
12/109
#0.1100917
28/109
#0.2568807

#######
#groups
boxplot(ad.fun.all$Freq ~ ad.fun.all$Variable, xlab=" ", ylab="Number of adjacent links")
boxplot(fg.all$x ~ fg.all$Variable, xlab=" ", ylab="Number of adjacent links")

#Figure boxplots fr supplement
par(mfrow=c(1,2))
#t = 0.58944, df = 30, p-value = 0.56
#new: t = 0.07092, df = 23, p-value = 0.9441
boxplot(fg.all$x ~ fg.all$Variable, xlab=" ", ylab="Number of neighbors",col="#69b3a2",main="Functional groups")
text(1,49,"t(23)=0.07092 ",cex=0.8,font=3)
text(0.875,47,"p=0.9441",cex=0.8,font=3)
text(1,19,"n=12",cex=0.8,font=3)
text(2,15,"n=13",cex=0.8,font=3)

#t = 4.2636, df = 40, p-value = 0.0001187
#new: t = 4.4415, df = 77, p-value = 2.948e-05
boxplot(ad.fun.all$Freq ~ ad.fun.all$Variable, xlab=" ", ylab="Number of neighbors", col="steelblue3", main="Families", ylim=c(0,8))
text(2,8,"t(77)=4.4415",cex=0.8,font=3)
text(2,7.7,"p=2.948e-05",cex=0.8,font=3)
text(1,7.2,"n=28",cex=0.8,font=3)
text(2,4.2,"n=51",cex=0.8,font=3)


#prepare fig. for supplement
#####################################################
list.sst=group.sst
list.ip25=group.ip25
#
list.freshwater=c("Acanthocerataceae", "Bracteacoccaceae", 
                  "Chloroflexaceae", "Closteriaceae", "Cyanophoraceae", 
                  "Golenkiniaceae",
                  "Haematococcaceae", "Mesotaeniaceae", "Oedogoniaceae", 
                  "Roseiflexaceae", "Sarcinofilaceae", "Scotinosphaeraceae", 
                  "Selenastraceae", "Trebouxiaceae", "Zygnemataceae")

sst.taxa = subset(res.plankt_pel.t, (row.names(res.plankt_pel.t) %in% list.sst))
ip25.taxa = subset(res.plankt_pel.t, (row.names(res.plankt_pel.t) %in% list.ip25))
freshw.taxa = subset(res.plankt_pel.t, (row.names(res.plankt_pel.t) %in% list.freshwater))
sst.taxa.sum=as.data.frame(colSums(sst.taxa))
ip25.taxa.sum=as.data.frame(colSums(ip25.taxa))
freshw.taxa.sum=as.data.frame(colSums(freshw.taxa))
new.plot=cbind(sst.taxa.sum, ip25.taxa.sum, freshw.taxa.sum)

par(mfrow=c(2,1), mar=c(4,4,1,1))
z=seq(0,20,1)#
plot(new.plot$`colSums(ip25.taxa)`~as.numeric(row.names(new.plot)), type="l", xlab="Age (cal kyr BP)", ylab="read counts", ylim=c(200,3500), col="steelblue3");axis(side = 1, at = z,labels = T)
lines(new.plot$`colSums(sst.taxa)`~as.numeric(row.names(new.plot)), type="l", xaxt="none", col="red", ylim=c(200,3500))
legend("topright", legend=c("correlated with IP25", "correlated with SST"),
       bty="n",col=c("steelblue3", "red"), lty=c(1,1))

plot(new.plot$`colSums(freshw.taxa)`~as.numeric(row.names(new.plot)), type="l",  xlab="Age (cal kyr BP)", ylab="read counts", col="seagreen3", ylim=c(0,100));axis(side = 1, at = z,labels = T)
legend("topright", legend=c("freshwater/terrestrial taxa"),
       bty="n",col=c("seagreen3"), lty=c(1,1))


#ecocopula network
##############################################
sstip <-read.table("Interpolations_2022.txt",header=TRUE, sep="\t", stringsAsFactors = FALSE)
sstip$Age <- NULL
pel.data.sst <- cbind(sstip, res.plankt_pel)
res.plankt_pel.ma = as.matrix(res.plankt_pel)

# ROUND UP VALUES FOR ECOCOPULA
pel.data.filt.ma.rnd <- ceiling(res.plankt_pel.ma)
pel_mod <- stackedsdm(pel.data.filt.ma.rnd,~., data = sstip, ncores=4)

ecocopula.pel=cgr(pel_mod, lambda = 0.51)
plot(ecocopula.pel,pad=1)


isolated.nodes.ecocopula = which(igraph::degree(ecocopula.pel$best_graph$igraph_out)==0)
#remove isolated nodes (=nodes without edges)
eg2=igraph::delete.vertices(ecocopula.pel$best_graph$igraph_out, isolated.nodes.ecocopula)
gsize(eg2)
#473 edes in graph
gorder(eg2)
#147 nodes in  graph

#calculate node degrees for each family
node.degrees=as.data.frame(igraph::degree(eg2))

res.pp.colSums=as.data.frame(colSums(res.plankt_pel))
names(res.pp.colSums)[1]<-"V1"
ecocopula.pel$best_graph$graph
taxa.list.size.ec = subset(res.pp.colSums, (row.names(res.pp.colSums) %in% row.names(node.degrees)))
dim(taxa.list.size.ec)
V(eg2)$size <- log1p(taxa.list.size.ec$V1)
V(eg2)$label <-row.names(taxa.list.size.ec)
par(bg="white", mar=c(0,0,0,0))
set.seed(21)
plot(eg2,
     #vertex.label=ifelse(V(g2)$size > 4, V(g2)$label, NA),
     vertex.label.cex=0.5,
     vertex.label.color="Black",
     vertex.label.family="sans",
     vertex.frame.color="white",
     layout=layout.fruchterman.reingold 
)

