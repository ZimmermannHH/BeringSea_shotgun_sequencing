################################################################################
#R-script 
#author: Heike Zimmermann
#aim: preparing taxonomic groups
################################################################################
#load libraries
library(dplyr)
library(tidyr)
library(reshape2)


#load data
fam.filter<-read.table("fam.wide.new.t6.txt",header=TRUE, dec=".", sep="\t",stringsAsFactors = FALSE)
names(fam.filter) <-gsub("X","",names(fam.filter))
fam.filter[26]=NULL
dim(fam.filter)

#sum counts per family and apply filter: at least 10 counts to be kept for analysis
fam.filter$count.sum<-rowSums(fam.filter)
fam.filter=subset(fam.filter, count.sum>9)
fam.filter[26]=NULL
#1694 families with at lest 10 counts left

#add a column called fam.names to filter the data
fam.filter$fam.names<-row.names(fam.filter)

#read lists of family names in taxonomic groups
list.all.phot.prot<-read.table("list.all.phot.prot.txt",header=FALSE, stringsAsFactors = FALSE)
list.phot_bacteria<-read.table("list.phot_bacteria.txt",header=FALSE, stringsAsFactors = FALSE)
list.all.het.prot<-read.table("list.all.het.prot.txt",header=FALSE, stringsAsFactors = FALSE)
list.crustaceae<-read.table("list.crustaceae",header=FALSE, stringsAsFactors = FALSE)
list.cnidaria<-read.table("list.cnidaria",header=FALSE, stringsAsFactors = FALSE)
list.knorpelfische<-read.table("list_knorpelfische.txt",header=FALSE, stringsAsFactors = FALSE)
list.fish<-read.table("list.fish",header=FALSE, stringsAsFactors = FALSE)
list.metazoa<-read.table("list.metazoa",header=FALSE, stringsAsFactors = FALSE)
list.chaetognatha<-read.table("list.chaetognatha",header=FALSE, stringsAsFactors = FALSE)
list.rotifera<-read.table("list.rotifera",header=FALSE, stringsAsFactors = FALSE)
list.echinodermata<-read.table("list.echinodermata",header=FALSE, stringsAsFactors = FALSE)
list.tunicata<-read.table("list.tunicata",header=FALSE, stringsAsFactors = FALSE)
list.crustaceae<-read.table("list.crustaceae",header=FALSE, stringsAsFactors = FALSE)
list.mollusca<-read.table("list.mollusca",header=FALSE, stringsAsFactors = FALSE)
list.cnidaria<-read.table("list.cnidaria",header=FALSE, stringsAsFactors = FALSE)
list.porifera<-read.table("list.porifera",header=FALSE, stringsAsFactors = FALSE)
list.polychaeta<-read.table("list.polychaeta",header=FALSE, stringsAsFactors = FALSE)
list.nematoda<-read.table("nematoda.txt",header=FALSE, stringsAsFactors = FALSE)
list.platyhelminthes<-read.table("list.platyhelminthes.txt",header=FALSE, stringsAsFactors = FALSE)
list.nemertea<-read.table("list.nemertea.txt",header=FALSE, stringsAsFactors = FALSE)
list.chelicerata<-read.table("list.chelicerata.txt",header=FALSE, stringsAsFactors = FALSE)
list.placozoa<-read.table("list.placozoa",header=FALSE, stringsAsFactors = FALSE)
list.rhodophyta<-read.table("list.rhodophyta",header=FALSE, stringsAsFactors = FALSE)
list.phaeophyceae<-read.table("list.phaeophyceae",header=FALSE, stringsAsFactors = FALSE)
list.rhodophyta<-read.table("list.rhodophyta",header=FALSE, stringsAsFactors = FALSE)
list.gastrotricha<-read.table("list.gastrotricha",header=FALSE, stringsAsFactors = FALSE)
list.bryozoa<-read.table("list.bryozoa",header=FALSE, stringsAsFactors = FALSE)
list.brachiopoda<-read.table("list.brachiopoda",header=FALSE, stringsAsFactors = FALSE)

#planktonic+pelagic
####################################################

#phototrophic protists
#######################

List.phot_bacteria=list.phot_bacteria$V1
phot_bacteria.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.phot_bacteria))
phot_bacteria.filter=as.data.frame(phot_bacteria.filter)
row.names(phot_bacteria.filter)=phot_bacteria.filter$fam.names
phot_bacteria.filter[26]=NULL 

List.phot=list.all.phot.prot$V1
phot.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.phot))
phot.filter=as.data.frame(phot.filter)
row.names(phot.filter)=phot.filter$fam.names
phot.filter[26]=NULL
sort(row.names(phot.filter))

#new families in new database: add the following taxa:
list.add.plankt.phot=c("Chromeraceae","Mesodiniidae","Pyrocystaceae")
add.plankt.phot.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% list.add.plankt.phot))
add.plankt.phot.filter[26]=NULL
phototrophic.protists.bacteria=as.data.frame(rbind(phot_bacteria.filter, phot.filter, add.plankt.phot.filter))
write.table(phototrophic.protists.bacteria, "2_outclean_phototrophic_protists_bacteria_2022_07_15.txt", sep="\t")


#heterotrophic protists
#######################

List.all.het.prot=list.all.het.prot$V1
het_prot.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.all.het.prot))
het_prot.filter=as.data.frame(het_prot.filter)
row.names(het_prot.filter)=het_prot.filter$fam.names

het_prot.filter[26]=NULL
#remove benthic heterotrophs, phototrophic taxa and parasitic Amoebae / Apicomplexa
#Neglaeria: thermophilic; Vahlkampfiidae --> most strains are free-living and non-pathogenic
list.benthic.heterotrophs=c("Capitellidae", "Cercomonadidae","Didiniidae","Holostichidae","Maldanidae","Mastigamoebidae",
                            "Nereididae","Philasteridae","Philodinidae","Protaspidae","Raperosteliaceae","Scalibregmatidae",
                            "Serpulidae","Siboglinidae","Spionidae", "Gromiidae", "Trichoplacidae")
list.parasitic.amoebae=c("Acanthamoebidae","Acytosteliaceae","Cavenderiaceae","Flabellulidae","Amoebidae",
                         "Paramoebidae","Physaraceae","Eudubosquellidae","Babesiidae", "Entamoebidae", "Dictyosteliaceae")
list.conoidasida<-read.table("list_conoidasida.txt",header=FALSE, stringsAsFactors = FALSE)
List.conoidasida=list.conoidasida$V1
het_prot.pel = subset(het_prot.filter, !(row.names(het_prot.filter) %in% list.benthic.heterotrophs))
het_prot.pel.2=subset(het_prot.pel, !(row.names(het_prot.pel) %in% list.parasitic.amoebae))
het_prot.pel.3=subset(het_prot.pel.2, !(row.names(het_prot.pel.2) %in% List.conoidasida))
het_prot.pel.4=subset(het_prot.pel.3, !(row.names(het_prot.pel.3) %in% list.add.plankt.phot))

sort(row.names(het_prot.pel.4))

#still includes parasites
write.table(het_prot.pel.4, "2_outclean_heterotrophic_pelagic.protists_2022_07_15_new.txt", sep="\t")


#combine planktonic photo- and heterotrophs
#######################################

plankton<-as.data.frame(rbind(phototrophic.protists.bacteria, het_prot.pel.4))
plankton2=  plankton[apply(plankton[,-1], 1, function(x) !all(x==0)),]
sort(row.names(plankton))

#still includes parasites
write.table(plankton2, "2_outclean_plankton_2022_07_15_new.txt", sep="\t")
dim(plankton2)
#######################################

#only Sagittidae - are included in heterotrophic plankton
# List.chaetognatha=list.chaetognatha$V1
# chaetognatha.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.chaetognatha))
# chaetognatha.filter=as.data.frame(chaetognatha.filter)
# row.names(chaetognatha.filter)=chaetognatha.filter$fam.names
# chaetognatha.filter[26]=NULL

#Rotifera
List.rotifera=list.rotifera$V1
rotifera.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.rotifera))
rotifera.filter=as.data.frame(rotifera.filter)
row.names(rotifera.filter)=rotifera.filter$fam.names
rotifera.filter[26]=NULL

#Crustaceae
List.crustaceae=list.crustaceae$V1
crustaceae.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.crustaceae))
crustaceae.filter=as.data.frame(crustaceae.filter)
row.names(crustaceae.filter)=crustaceae.filter$fam.names
crustaceae.filter[26]=NULL
row.names(crustaceae.filter)

rm.terrestrial.paras=c("Artemiidae","Artemiidae","Cambaridae","Caligidae")
list.crustaceae.pelagic=c("Calanidae","Metridinidae","Temoridae","Eucalanidae","Euphausiidae","Oithonidae","Moinidae",
                          "Cyclopettidae","Acartiidae")
list.crustaceae.benthic=c("Harpacticidae","Daphniidae","Parastacidae","Balanidae","Varunidae","Hyalellidae","Penaeidae",
                          "Lysianassidae","Hyalidae","Pollicipedidae","Allocrangonyctidae","Nephropidae","Tetragonicipitidae",
                          "Palaemonidae","Cypridinidae", "Portunidae")

crustaceae.pelagic.filter2= subset(crustaceae.filter, !(row.names(crustaceae.filter) %in% rm.terrestrial.paras))
crustaceae.pelagic.filter3= subset(crustaceae.pelagic.filter2, !(row.names(crustaceae.pelagic.filter2) %in% list.crustaceae.benthic))

#Cnidara
List.cnidaria=list.cnidaria$V1
cnidaria.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.cnidaria))
cnidaria.filter=as.data.frame(cnidaria.filter)
row.names(cnidaria.filter)=cnidaria.filter$fam.names
cnidaria.filter[26]=NULL
list.cnidaria.parasitic="Myxobolidae"
list.cnidaria.pelagic=c("Hydridae","Pelagiidae","Ulmaridae","Cyaneidae","Eirenidae")
list.cnidaria.benthic=c("Acroporidae","Actiniidae","Aiptasiidae","Edwardsiidae","Hydractiniidae","Merulinidae","Nephtheidae",
                        "Pocilloporidae" ,"Sertulariidae","Alcyoniidae")
cnidaria.filter2= subset(cnidaria.filter, !(row.names(cnidaria.filter) %in% list.cnidaria.parasitic))
cnidaria.filter2.pelagic= subset(cnidaria.filter2, !(row.names(cnidaria.filter2) %in% list.cnidaria.benthic))

rowSums(cnidaria.filter)


#Metazoa
List.metazoa=list.metazoa$V1
metazoa.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.metazoa))
metazoa.filter=as.data.frame(metazoa.filter)
row.names(metazoa.filter)=metazoa.filter$fam.names
metazoa.filter[26]=NULL
#list.terrestrial.mammals=c("Bovidae","Camelidae","Canidae","Bathyergidae","Castoridae","Caviidae","","","","","","","","","","")
list.marine.mammals=c("Balaenopteridae","Monodontidae","Delphinidae","Physeteridae","Balaenidae","Scyliorhinidae",
                      "Phocoenidae","Ziphiidae","Eschrichtiidae","Lipotidae","Phocidae","Otariidae","Odobenidae", "Mustelidae")
metazoa.filter2= subset(metazoa.filter, (row.names(metazoa.filter) %in% list.marine.mammals))

#Chondrychthyes
List.chondr=list.knorpelfische$V1
chondr.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.chondr))
chondr.filter=as.data.frame(chondr.filter)
row.names(chondr.filter)=chondr.filter$fam.names
chondr.filter[26]=NULL

#Fish
List.fish=list.fish$V1
fish.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.fish))
fish.filter=as.data.frame(fish.filter)
row.names(fish.filter)=fish.filter$fam.names

list.fish.clean=c("Acipenseridae","Anarhichadidae", "Anguillidae", "Apogonidae", "Batrachoididae", "Blenniidae", "Bovichtidae", "Carangidae", 
                  "Salmonidae", "Sebastidae", "Gobiidae", "Gadidae", "Sparidae", "Clupeidae", "Fundulidae", "Osmeridae"
                  , "Channichthyidae", "Syngnathidae", "Cottidae", "Lateolabracidae", "Cyclopteridae", "Engraulidae", "Gasterosteidae", 
                  "Myctophidae", "Paralichthyidae", "Holocentridae", "Pleuronectidae", "Tetraodontidae",
                  "Serranidae", "Sciaenidae", "Trichiuridae")
fish.filter2=fish.filter %>% filter_at(vars(fam.names), any_vars(. %in% list.fish.clean))
fish.filter2[26]=NULL

#Mollusca
List.mollusca=list.mollusca$V1
mollusca.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.mollusca))
mollusca.filter=as.data.frame(mollusca.filter)
row.names(mollusca.filter)=mollusca.filter$fam.names
mollusca.filter[26]=NULL
octopodidae<-subset(mollusca.filter, row.names(mollusca.filter)=="Octopodidae")

pelagic<-as.data.frame(rbind(phototrophic.protists.bacteria,
                             het_prot.pel.4,
                             #rotifera.filter, #incl in het prot
                             #chaetognatha.filter, #incl in het prot
                             crustaceae.pelagic.filter3, 
                             cnidaria.filter2.pelagic,
                             octopodidae,
                             fish.filter2,
                             chondr.filter,
                             metazoa.filter2))

pelagic2=  pelagic[apply(pelagic[,-1], 1, function(x) !all(x==0)),]

#remove parasites
#Nematoda
List.nematoda=list.nematoda$V1
nematoda.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.nematoda))
nematoda.filter=as.data.frame(nematoda.filter)
row.names(nematoda.filter)=nematoda.filter$fam.names
nematoda.filter[26]=NULL
nematoda.freeliving<-subset(nematoda.filter, row.names(nematoda.filter)=="Oxystominidae")
nematoda.parasitic<-subset(nematoda.filter, row.names(nematoda.filter)!="Oxystominidae")
list.nematoda.parasitic=row.names(nematoda.parasitic)

parasites=c("Dicrocoeliidae","Hymenolepididae", "Opisthorchiidae", "Schistosomatidae", 
            "Taeniidae", "Onchocercidae", "Trypanosomatidae", "Trichomonadidae", 
            "Theileriidae", "Plasmodiophoridae", "Plasmodiidae", "Dictyosteliaceae")
freshwater.snails=c("Planorbidae","Geoplanidae", "Dichotomosiphonaceae")

parasitic.myxozoa<-subset(cnidaria.filter, row.names(cnidaria.filter)=="Saccosporidae" #parasites of fish and freshwater bryozoans
                          |row.names(cnidaria.filter)=="Myxidiidae" #parasites of freshwater fish
                          |row.names(cnidaria.filter)=="Myxobolidae" #
                          |row.names(cnidaria.filter)=="Kudoidae" #
                          |row.names(cnidaria.filter)=="Sphaerosporidae" #
                          |row.names(cnidaria.filter)=="Ceratomyxidae" #
                          |row.names(cnidaria.filter)=="Monobrachidae" #
)


list.myxozoa.parasitic=row.names(parasitic.myxozoa)

list.parasites=c(list.nematoda.parasitic,list.myxozoa.parasitic,parasites,freshwater.snails)

pelagic3= subset(pelagic2, !(row.names(pelagic2) %in% list.parasites))

sort(row.names(pelagic3))
write.table(pelagic3, "2_outclean_pelagic_planktonic_2022_07_15_new.txt", sep="\t")


#benthic
####################################################
#benthic heterotrophics
benthic.heterotrophs.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% list.benthic.heterotrophs))
benthic.heterotrophs.filter[26]=NULL
#macrophytobenthos

#brown algae
List.phaeophyceae=list.phaeophyceae$V1
phaeophyceae.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.phaeophyceae))
phaeophyceae.filter=as.data.frame(phaeophyceae.filter)
row.names(phaeophyceae.filter)=phaeophyceae.filter$fam.names
phaeophyceae.filter[26]=NULL

#red algae
List.rhodophyta=list.rhodophyta$V1
rhodophyta.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.rhodophyta))
rhodophyta.filter=as.data.frame(rhodophyta.filter)
row.names(rhodophyta.filter)=rhodophyta.filter$fam.names
rhodophyta.filter[26]=NULL

#Zosteraceae
zosteraceae<-subset(fam.filter, fam.names=="Zosteraceae")
zosteraceae[26]=NULL

macrophytobenthos=as.data.frame(rbind(phaeophyceae.filter,rhodophyta.filter,zosteraceae))
write.table(macrophytobenthos, "2_outclean_macrophytobenthos_2022_07_15.txt", sep="\t")

###
List.mollusca=list.mollusca$V1
mollusca.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.mollusca))
mollusca.filter=as.data.frame(mollusca.filter)
row.names(mollusca.filter)=mollusca.filter$fam.names
mollusca.filter[26]=NULL
mollusca.filter2<-subset(mollusca.filter, row.names(mollusca.filter)!="Octopodidae")
mollusca.filter3<-subset(mollusca.filter2, row.names(mollusca.filter2)!="Planorbidae")#remove freshwater snails

#Gastrotricha
List.gastrotricha=list.gastrotricha$V1
gastrotricha.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.gastrotricha))
gastrotricha.filter=as.data.frame(gastrotricha.filter)
row.names(gastrotricha.filter)=gastrotricha.filter$fam.names
gastrotricha.filter[26]=NULL


#Bryozoa
List.bryozoa=list.bryozoa$V1
bryozoa.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.bryozoa))
bryozoa.filter=as.data.frame(bryozoa.filter)
row.names(bryozoa.filter)=bryozoa.filter$fam.names
bryozoa.filter[26]=NULL

#Brachyopoda
List.brachiopoda=list.brachiopoda$V1
brachiopoda.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.brachiopoda))
brachiopoda.filter=as.data.frame(brachiopoda.filter)
row.names(brachiopoda.filter)=brachiopoda.filter$fam.names
brachiopoda.filter[26]=NULL


#Chelicerata
List.chelicerata=list.chelicerata$V1
chelicerata.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.chelicerata))
chelicerata.filter=as.data.frame(chelicerata.filter)
row.names(chelicerata.filter)=chelicerata.filter$fam.names
chelicerata.filter[26]=NULL

#Nematoda
List.nematoda=list.nematoda$V1
nematoda.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.nematoda))
nematoda.filter=as.data.frame(nematoda.filter)
row.names(nematoda.filter)=nematoda.filter$fam.names
nematoda.filter[26]=NULL
nematoda.freeliving<-subset(nematoda.filter, row.names(nematoda.filter)=="Oxystominidae")

#Gromiidae benthic Endomyxa
benthic.endomyx<-subset(het_prot.filter, row.names(het_prot.filter)=="Gromiidae")

#Crustaceae
List.crustaceae=list.crustaceae$V1
crustaceae.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.crustaceae))
crustaceae.filter=as.data.frame(crustaceae.filter)
row.names(crustaceae.filter)=crustaceae.filter$fam.names
crustaceae.filter[26]=NULL
rm.terrestrial.paras=c("Artemiidae","Cambaridae","Caligidae")
list.crustaceae.pelagic=c("Calanidae","Metridinidae","Temoridae","Eucalanidae","Euphausiidae","Oithonidae","Moinidae",
                          "Cyclopettidae","Acartiidae")
list.crustaceae.benthic=c("Harpacticidae","Daphniidae","Parastacidae","Balanidae","Varunidae","Hyalellidae","Penaeidae",
                          "Lysianassidae","Hyalidae","Pollicipedidae","Allocrangonyctidae","Nephropidae","Tetragonicipitidae",
                          "Palaemonidae","Cypridinidae", "Portunidae")
crustaceae.pelagic.filter2= subset(crustaceae.filter, !(row.names(crustaceae.filter) %in% rm.terrestrial.paras))
crustaceae.benthic.filter3= subset(crustaceae.pelagic.filter2, !(row.names(crustaceae.pelagic.filter2) %in% list.crustaceae.pelagic))


#Cnidara
List.cnidaria=list.cnidaria$V1
cnidaria.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.cnidaria))
cnidaria.filter=as.data.frame(cnidaria.filter)
row.names(cnidaria.filter)=cnidaria.filter$fam.names
cnidaria.filter[26]=NULL
list.cnidaria.parasitic="Myxobolidae"
list.cnidaria.pelagic=c("Hydridae","Pelagiidae","Ulmaridae","Cyaneidae","Eirenidae")
list.cnidaria.benthic=c("Acroporidae","Actiniidae","Aiptasiidae","Edwardsiidae","Hydractiniidae","Merulinidae","Nephtheidae",
                        "Pocilloporidae" ,"Sertulariidae","Alcyoniidae")
cnidaria.filter2= subset(cnidaria.filter, !(row.names(cnidaria.filter) %in% list.cnidaria.parasitic))
cnidaria.filter2.benthic= subset(cnidaria.filter2, !(row.names(cnidaria.filter2) %in% list.cnidaria.pelagic))


#Echinodermata
List.echinodermata=list.echinodermata$V1
echinodermata.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.echinodermata))
echinodermata.filter=as.data.frame(echinodermata.filter)
row.names(echinodermata.filter)=echinodermata.filter$fam.names
echinodermata.filter[26]=NULL


#Tunicata
List.tunicata=list.tunicata$V1
tunicata.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.tunicata))
tunicata.filter=as.data.frame(tunicata.filter)
row.names(tunicata.filter)=tunicata.filter$fam.names
tunicata.filter[26]=NULL


#Porifera
List.porifera=list.porifera$V1
porifera.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.porifera))
porifera.filter=as.data.frame(porifera.filter)
row.names(porifera.filter)=porifera.filter$fam.names
porifera.filter[26]=NULL


#polychaeta
List.polychaeta=list.polychaeta$V1
polychaeta.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.polychaeta))
polychaeta.filter=as.data.frame(polychaeta.filter)
row.names(polychaeta.filter)=polychaeta.filter$fam.names
polychaeta.filter[26]=NULL


#Priapulidae
Priapulidae<-subset(fam.filter,row.names(fam.filter)=="Priapulidae")
Priapulidae[26]=NULL

#Plathihelminthes
#Trematoda endoparasitic
List.platyhelminthes=list.platyhelminthes$V1
platyhelminthes.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.platyhelminthes))
platyhelminthes.filter=as.data.frame(platyhelminthes.filter)
row.names(platyhelminthes.filter)=platyhelminthes.filter$fam.names
platyhelminthes.filter[26]=NULL

#Nemertea 
List.nemertea=list.nemertea$V1
nemertea.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% List.nemertea))
nemertea.filter=as.data.frame(nemertea.filter)
row.names(nemertea.filter)=nemertea.filter$fam.names
nemertea.filter[26]=NULL

#Placozoa

placozoa<-subset(fam.filter, row.names(fam.filter)=="Trichoplacidae")
placozoa[26]=NULL

#Oligochaeta nur sehr wenige counts (<10)
list.oligochaeta=c("Lumbricidae", "Megascolecidae")
oligochaeta.filter=fam.filter %>% filter_at(vars(fam.names), any_vars(. %in% list.oligochaeta))
oligochaeta.filter[26]=NULL


benthics<-as.data.frame(rbind(benthic.heterotrophs.filter,
                              macrophytobenthos,
                              #placozoa,
                              brachiopoda.filter,
                              #bryozoa.filter,
                              #gastrotricha.filter,
                              #benthic.endomyx,
                              platyhelminthes.filter,
                              Priapulidae,
                              nematoda.freeliving,
                              porifera.filter,
                              tunicata.filter,
                              crustaceae.benthic.filter3,
                              mollusca.filter3,
                              echinodermata.filter,
                              cnidaria.filter2.benthic,
                              chelicerata.filter,
                              oligochaeta.filter))
benthics2=benthics[apply(benthics[,-1], 1, function(x) !all(x==0)),]

#remove parasites
parasites=c("Dicrocoeliidae","Hymenolepididae", "Opisthorchiidae", "Schistosomatidae", "Taeniidae", "Onchocercidae", "Trypanosomatidae")
freshwater.snails=c("Planorbidae","Geoplanidae")
nematoda.parasitic<-subset(nematoda.filter, row.names(nematoda.filter)!="Oxystominidae")

parasitic.myxozoa<-subset(cnidaria.filter, row.names(cnidaria.filter)=="Saccosporidae" #parasites of fish and freshwater bryozoans
                          |row.names(cnidaria.filter)=="Myxidiidae" #parasites of freshwater fish
                          |row.names(cnidaria.filter)=="Myxobolidae" #
                          |row.names(cnidaria.filter)=="Kudoidae" #
                          |row.names(cnidaria.filter)=="Sphaerosporidae" #
                          |row.names(cnidaria.filter)=="Ceratomyxidae" #
                          |row.names(cnidaria.filter)=="Monobrachidae" #
)
Caligidae<-subset(crustaceae.filter, row.names(crustaceae.filter)=="Caligidae")

list.nematoda.parasitic=row.names(nematoda.parasitic)
list.myxozoa.parasitic=row.names(parasitic.myxozoa)

list.parasites=c(list.nematoda.parasitic,list.myxozoa.parasitic,parasites,freshwater.snails)

benthics3= subset(benthics2, !(row.names(benthics2) %in% list.parasites))
sort(row.names(benthics3))
write.table(benthics3, "2_outclean_benthics_2022_07_15_new.txt", sep="\t")                    
