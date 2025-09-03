setwd("~/OneDrive - UBC/DB_fung_2025.01/")

#Data_scraping
library(devtools)
require(testthat)
library(rusda)
library(rentrez)
library(pestr)
library(taxize)
library(finch)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(Taxonstand)
library(magrittr)
library(WorldFlora)
library(fuzzyjoin)
library(stringr)
library(FUNGuildR)
library(polyglotr)

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ GBIF ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# Downloaded May 12 2022.

#get GBIF associations
# downloaded "Fungi" occurrences. 30GB file.
# Citation for GBIF fungal download : GBIF.org (12 May 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.rb2k7k
#from bash
# grep -v -e 'KINGDOM' -e 'PHYLUM' -e 'CLASS' -e 'ORDER' -e 'FAMILY' occurrence.txt > occurence.no.KPCOF.txt
# { head -1 occurence.no.KPCOF.txt ; grep "SPECIES" occurence.no.KPCOF.txt ;} | cat > species.occurrence.txt
# split -l 11000000 species.occurrence.txt
# mv xaa species.occurrence.1.txt
# head -1 species.occurrence.txt | cat - xab > species.occurrence.2.txt
# { head -1 occurence.no.KPCOF.txt ; grep "GENUS" occurence.no.KPCOF.txt ;} | cat > genus.occurrence.txt

#in R
#####SPECIES#####
# fung1 <- data.table::fread("~/data/GBIF.fungi.occurrence/species.occurrence.1.txt")
# fung2 <- data.table::fread("~/data/GBIF.fungi.occurrence/species.occurrence.2.txt")
# in fung1 and fung2, there were some noted errors in uploading the file. 
# saveRDS(fung1, 'OSF/GBIF/fung.1.GBIF.check.for.errors.rds')
# saveRDS(fung2, 'OSF/GBIF/fung.2.GBIF.check.for.errors.rds')

#I
fung1$associatedTaxa<-dplyr::na_if(fung1$associatedTaxa, "")
fung1$habitat<-dplyr::na_if(fung1$habitat, "")
fung1$associatedOrganisms<-dplyr::na_if(fung1$associatedOrganisms, "")

#II
fung2$associatedTaxa<-dplyr::na_if(fung2$associatedTaxa, "")
fung2$habitat<-dplyr::na_if(fung2$habitat, "")
fung2$associatedOrganisms<-dplyr::na_if(fung2$associatedOrganisms, "")

#2. everything with putative host fields
#I
fung1 <- fung1 %>% filter_at(vars(associatedTaxa, habitat, associatedOrganisms),any_vars(!is.na(.)))
#facts <- sapply(fung1, is.factor)                          
#fung1[facts] <- lapply(fung1[facts], as.character)
#II
fung2 <- fung2 %>% filter_at(vars(associatedTaxa, habitat, associatedOrganisms),any_vars(!is.na(.)))
#facts <- sapply(fung2, is.factor)                          
#fung2[facts] <- lapply(fung2[facts], as.character)

fung.G<-rbind(fung1, fung2)
# save full data, adjusted to remove blank entries for associatedtaxa, etc.
# saveRDS(fung.G, 'OSF/GBIF/fung.G.rds')


#### Keep only relatively pertinent info (e.g., species name, occurrenceID)
fung.G.ltd<- fung.G[,c(1,57,58,68,90,93,96,111,120,125,189,203,239,240,241,242)]
# saveRDS(fung.G.ltd, 'OSF/GBIF/fung.G.ltd.rds')
fung.G.ltd<-data.frame(fung.G.ltd,db="GBIF.db")

#### Filter entries already processed by Mycoportal { Myco.po.ab<-readRDS("~/data2/LPs/rds/Myco.po.ab.rds") }

#some occurrence IDs are not unique. Find the minimum number of characters where there are no duplicates.
unique(duplicated(fung.G.ltd$occurrenceID[nchar(fung.G.ltd$occurrenceID) > 36]))
#the answer is > 36 characters. New column where only known non-duplicated entries are.
#however, it looks like there are useful entries that are 36 characters long. These are long enough to generally trust, if duplicates removed.
#keep those over 35 characters
fung35.GB<-fung.G.ltd[nchar(fung.G.ltd$occurrenceID) > 35,]
#keep non-duplicated for consideration for removal if matches MycoPortal.
non.dups.GB<-fung35.GB[which(!duplicated(fung35.GB$occurrenceID)),]
#remove matching to Mycoportal entries from dataframe
in.Mycop<-which(non.dups.GB$occurrenceID %in% Myco.po.ab$occurrenceID)
fung.G.ltd.nomyco<-fung.G.ltd[-in.Mycop,]

#prepare file for taxonomic name extraction.
associates<- data.frame(associatedTaxa=fung.G.ltd.nomyco$associatedTaxa, habitat=fung.G.ltd.nomyco$habitat, associatedOrganisms=fung.G.ltd.nomyco$associatedOrganisms)
associates<-unique(associates)
# NEW REVISED
# write.csv(associates, 'OSF/GBIF/associations.new.GBIF.csv')


# BASH: 
# split -l 550000 associations.new.GBIF.csv 
# rename:

#### run gnfinder to locate scientific names.
# gnfinder -i ./associations.GBIF.new.1.csv -u > ./associations.GBIF.new.1.gnfinder.results.NOuninomial.csv
# gnfinder -i ./associations.GBIF.new.2.csv -u > ./associations.GBIF.new.2.gnfinder.results.NOuninomial.csv

### navigate back to R
gnf1<-read.csv('OSF/GBIF/associations.GBIF.new.1.gnfinder.results.Nouninomial.csv')
gnf2<-read.csv('OSF/GBIF/associations.GBIF.new.2.gnfinder.results.Nouninomial.csv')
#stick them together.
gnf<-rbind(gnf1,gnf2)
# sort by unique name
gnfName<-sort(unique(gnf$Name))

# probably delete:
#remove 'As'
#gnfName.1<-gnfName[-3490]

# this should speed up loop

fung.G.ltd.hab <- fung.G.ltd.nomyco %>% filter_at(vars(habitat),any_vars(!is.na(.))) 

fung.G.ltd.asOr <- fung.G.ltd.nomyco %>% filter_at(vars(associatedOrganisms),any_vars(!is.na(.))) 

fung.G.ltd.asTa <- fung.G.ltd.nomyco %>%  filter_at(vars(associatedTaxa),any_vars(!is.na(.))) 


#make loops that will search for all possible in habitat, assocOrg, and assocTaxa fields and make new entry for them.
GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in gnfName){
  mm<-ifelse(nn<-grep(i, fung.G.ltd.hab$habitat, ignore.case = TRUE),print(i),0)
  GBIF.plant<-c(GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
hits.G.habitat<-data.frame(GBIF.entry,GBIF.plant)

GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in gnfName){
  mm<-ifelse(nn<-grep(i, fung.G.ltd.asOr$associatedOrganisms, ignore.case = TRUE),print(i),0)
  GBIF.plant<-c(GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
hits.G.associated.Org<-data.frame(GBIF.entry,GBIF.plant)

GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in gnfName){
  mm<-ifelse(nn<-grep(i, fung.G.ltd.asTa$associatedTaxa, ignore.case = TRUE),print(i),0)
  GBIF.plant<-c(GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
hits.G.associated.Taxa<-data.frame(GBIF.entry,GBIF.plant)

#make sure to save in case i need to troubleshoot the individual files.
habGBIF.o<-data.frame(GBIF.entry=hits.G.habitat$GBIF.entry, GBIF.plant=hits.G.habitat$GBIF.plant, fung.G.ltd.hab[hits.G.habitat$GBIF.entry,])
# saveRDS(habGBIF.o, 'OSF/GBIF/habGBIF.o.new.rds')
astaGBIF.o<-data.frame(GBIF.entry=hits.G.associated.Taxa$GBIF.entry, GBIF.plant=hits.G.associated.Taxa$GBIF.plant, fung.G.ltd.asTa[hits.G.associated.Taxa$GBIF.entry,])
# saveRDS(astaGBIF.o, 'OSF/GBIF/astaGBIF.o.new.rds')
asorGBIF.o<-data.frame(GBIF.entry=hits.G.associated.Org$GBIF.entry, GBIF.plant=hits.G.associated.Org$GBIF.plant, fung.G.ltd.asOr[hits.G.associated.Org$GBIF.entry,])
# saveRDS(asorGBIF.o, 'OSF/GBIFs/asorGBIF.o.new.rds')


habGBIF<-habGBIF.o[,-1]
astaGBIF<-astaGBIF.o[,-1]
asorGBIF<-asorGBIF.o[,-1]

# All fungi and plants, at least in same geographic region. Not at species level.
GBIF.hits<-unique(rbind(habGBIF,astaGBIF,asorGBIF))

###the order is mixed up, have to sort by gbif id

#Keep only putative "hosts" at species level and below. 

word.count=NULL
for (i in GBIF.hits$GBIF.plant){
  word.count=c(word.count, print(length(strsplit(i, " ")[[1]])))
}
GBIF.hits<-data.frame(GBIF.hits, word.count)

GBIF.hits.sp <- GBIF.hits %>% 
  filter(!GBIF.hits$word.count=="1")

# These represent association data (does it have a latin name species name in any text field). These imply some level of plant-fungal association. 
# save GBIF.hits to host species level.


#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# I. Here we create GBIF.hits.sp, which will be used as raw data for the "Locally Co-occurring" database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# saveRDS(GBIF.hits.sp, 'OSF/GBIF/GBIF.hits.sp.new.rds')

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


###Part 2: Filter to only fungi ON PLANT TISSUE. This is stricter, and mostly filters out habitat fields that = GBIF plant name fields.

#Remove hits where there is only one habitat or associatedTaxa or associatedOrgs field that is just a plant name (no preposition, eg), per gnfinder.
#make sure all other putative host fields are NA.
GB.hab.na<- GBIF.hits.sp %>% filter(is.na(associatedOrganisms) & is.na(associatedTaxa))
#find direct matches between habitat field and gnf Name
GB.hab.gnf<- data.frame(GB.hab.na[GB.hab.na$habitat %in% gnf$Name,])
#many just have "matrix:" in front of the putative host name. find these
GB.hab.f1<-data.frame(GB.hab.na[which(GB.hab.na$habitat %in% (paste("Matrix:", GB.hab.na$GBIF.plant))),])
GB.hab.f2<-data.frame(GB.hab.na[which(GB.hab.na$habitat %in% (paste(GB.hab.na$GBIF.plant, "forest"))),])
GB.hab.f3<-data.frame(GB.hab.na[which(GB.hab.na$habitat %in% (paste(GB.hab.na$GBIF.plant, "L."))),])
GB.hab.f4<-data.frame(GB.hab.na[which(GB.hab.na$habitat %in% (paste("Matr.:", GB.hab.na$GBIF.plant))),])
GB.hab.f5<-data.frame(GB.hab.na[which(GB.hab.na$habitat %in% (paste(GB.hab.na$GBIF.plant, "L."))),])
GB.hab.f6<-data.frame(GB.hab.na[which(GB.hab.na$habitat %in% (paste("Con", GB.hab.na$GBIF.plant))),])
GB.hab.f7<-data.frame(GB.hab.na[which(GB.hab.na$habitat %in% (paste("Near", GB.hab.na$GBIF.plant))),])
GB.hab.f8<-data.frame(GB.hab.na[which(GB.hab.na$habitat %in% (paste(GB.hab.na$GBIF.plant, "stand"))),])
GB.hab.f9<-data.frame(GB.hab.na[which(GB.hab.na$habitat %in% (paste("with", GB.hab.na$GBIF.plant))),])
GB.hab.f10<-data.frame(GB.hab.na[which(GB.hab.na$habitat %in% (paste("{", GB.hab.na$GBIF.plant, "} forest.", sep =""))),])

###
GB.hab.rem<-rbind(GB.hab.gnf, GB.hab.f1,GB.hab.f2,GB.hab.f3,GB.hab.f4,GB.hab.f5,GB.hab.f6,GB.hab.f7,GB.hab.f8, GB.hab.f9,GB.hab.f10)
saveRDS(GB.hab.rem, "intermediate_rds_csv/GBIF/GB.hab.rem.rds")

#filter out gnf Name direct matches and others.
GB.hab.filt<-GB.hab.na[which(!rownames(GB.hab.na) %in% rownames(GB.hab.rem)),]

#now this has been pared down, but it is still > 1mil entries. 
GB.plant<-unique(GB.hab.filt$GBIF.plant)

# Many plants are included as a species list of the place. Filter these out (grep ", TAXON, ")
GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in GB.plant){
  mm<-ifelse(nn<-grep(paste(", ", i, ", ", sep=""), GB.hab.filt$habitat, ignore.case = TRUE),print(i),0)
  GBIF.plant<-c(GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
sp.in.list<-data.frame(GBIF.entry,GBIF.plant)
# write.csv(sp.in.list, "intermediate_rds_csv/GBIF/sp.in.list.csv")
#rename GBIF.plant column of sp.in.list
names(sp.in.list)[names(sp.in.list) == "GBIF.plant"] <- "sp.in.li.GBIF.plant"
aa<-cbind(GB.hab.filt[sp.in.list$GBIF.entry,], sp.in.list)
#only delete relevant "GBIF.plant"-matching row. 
aa1<-aa %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))

##
#later, created this file to keep these lists of data we don't want to keep
sil.GB<-GB.hab.filt[sp.in.list$GBIF.entry,]
saveRDS(sil.GB, 'intermediate_rds_csv/GBIF/sil.GB.rds')
##


# Filter out "TAXON-". This often precedes "forest"
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in GB.plant){
  mm<-ifelse(nn<-grep(paste(i, "-", sep=""), GB.hab.filt$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
Gp.dash<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
write.csv(Gp.dash, "intermediate_rds_csv/GBIF/GBIF.plant.dash.intermediatestep.csv")
bb<-cbind(GB.hab.filt[Gp.dash$GBIF.entry,], Gp.dash)
bb1<-bb %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))

##
#later, created this file to keep these lists of data we don't want to keep
Gpd.GB<-GB.hab.filt[Gp.dash$GBIF.entry,]
saveRDS(Gpd.GB, 'intermediate_rds_csv/GBIF/Gpd.GB.rds')
##

# Filter out "TAXON forest"
GB.hab.filt.forest<-GB.hab.filt[grep("forest", GB.hab.filt$habitat, ignore.case=TRUE), ]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in GB.plant){
  mm<-ifelse(nn<-grep(paste(i, "forest"), GB.hab.filt.forest$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
Gpf<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
write.csv(Gpf, "intermediate_rds_csv/GBIF/GBIF.plant.forest.intermediatestep.csv")
cc<-cbind(GB.hab.filt.forest[Gpf$GBIF.entry,], Gpf)
cc1<-cc %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))

##
#later, created this file to keep these lists of data we don't want to keep
Gpf.GB<-GB.hab.filt.forest[Gpf$GBIF.entry,]
saveRDS(Gpf.GB, 'intermediate_rds_csv/GBIF/Gpf.GB.rds')
##

# Filter out "TAXON woodland"
GB.hab.filt.woodland<-GB.hab.filt[grep("woodland", GB.hab.filt$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in GB.plant){
  mm<-ifelse(nn<-grep(paste(i, "woodland"), GB.hab.filt.woodland$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
Gpw<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
write.csv(Gpw, "intermediate_rds_csv/GBIF/GBIF.plant.woodland.intermediatestep.csv")
dd<-cbind(GB.hab.filt.woodland[Gpw$GBIF.entry,], Gpw)
dd1<-dd %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))

##
#later, created this file to keep these lists of data we don't want to keep
Gpw.GB<-GB.hab.filt.woodland[Gpw$GBIF.entry,]
saveRDS(Gpw.GB, 'intermediate_rds_csv/GBIF/Gpw.GB.rds')
##

# Filter out "TAXON stand"
GB.hab.filt.stand<-GB.hab.filt[grep("stand", GB.hab.filt$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in GB.plant){
  mm<-ifelse(nn<-grep(paste(i, "stand"), GB.hab.filt.stand$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
Gps<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
write.csv(Gps, "intermediate_rds_csv/GBIF/GBIF.plant.stand.intermediatestep.csv")
ee<-cbind(GB.hab.filt.stand[Gps$GBIF.entry,], Gps)
ee1<-ee %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))

##
#later, created this file to keep these lists of data we don't want to keep
Gps.GB<-GB.hab.filt.stand[Gps$GBIF.entry,]
saveRDS(Gps.GB, 'intermediate_rds_csv/GBIF/Gps.GB.rds')
##

# Filter out "forest with TAXON"
GB.hab.filt.fwith<-GB.hab.filt[grep("forest with", GB.hab.filt$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in GB.plant){
  mm<-ifelse(nn<-grep(paste("forest with", i), GB.hab.filt.fwith$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
fwG<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
write.csv(fwG, "intermediate_rds_csv/GBIF/forest.w.plant.GBIF.intermediatestep.csv")
ff<-cbind(GB.hab.filt.fwith[fwG$GBIF.entry,], fwG)
ff1<-ff %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))

##
#later, created this file to keep these lists of data we don't want to keep
fwG.GB<-GB.hab.filt.fwith[fwG$GBIF.entry,]
saveRDS(fwG.GB, 'intermediate_rds_csv/GBIF/fwG.GB.rds')
##

# Filter out "forest of TAXON"
GB.hab.filt.fof<-GB.hab.filt[grep("forest of", GB.hab.filt$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in GB.plant){
  mm<-ifelse(nn<-grep(paste("forest of", i), GB.hab.filt.fof$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
foG<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
write.csv(foG, "intermediate_rds_csv/GBIF/forest.of.plant.GBIF.intermediatestep.csv")
gg<-cbind(GB.hab.filt.fof[foG$GBIF.entry,], foG)
gg1<-gg %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))

##
#later, created this file to keep these lists of data we don't want to keep
foG.GB<-GB.hab.filt.fof[foG$GBIF.entry,]
saveRDS(foG.GB, 'intermediate_rds_csv/GBIF/foG.GB.rds')
##

# Filter out "bosque de TAXON"
GB.hab.filt.bde<-GB.hab.filt[grep("bosque de", GB.hab.filt$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in GB.plant){
  mm<-ifelse(nn<-grep(paste("bosque de", i), GB.hab.filt.bde$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
bdG<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
#write.csv(bdG, "intermediate_rds_csv/GBIF/bosque.d.plant.GBIF.intermediatestep.csv")
hh<-cbind(GB.hab.filt.bde[bdG$GBIF.entry,], bdG)
hh1<-hh %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))

##
#later, created this file to keep these lists of data we don't want to keep
bdG.GB<-GB.hab.filt.bde[bdG$GBIF.entry,]
saveRDS(bdG.GB, 'intermediate_rds_csv/GBIF/bdG.GB.rds')
##

# Filter out "bosque con TAXON"
GB.hab.filt.bcon<-GB.hab.filt[grep("bosque con", GB.hab.filt$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in GB.plant){
  mm<-ifelse(nn<-grep(paste("bosque con", i), GB.hab.filt.bcon$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
bcG<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
write.csv(bcG, "intermediate_rds_csv/GBIF/bosque.con.plant.GBIF.intermediatestep.csv")
ii<-cbind(GB.hab.filt.bcon[bcG$GBIF.entry,], bcG)
ii1<-ii %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))

##
#later, created this file to keep these lists of data we don't want to keep
bcG.GB<-GB.hab.filt.bcon[bcG$GBIF.entry,]
saveRDS(bcG.GB, 'intermediate_rds_csv/GBIF/bcG.GB.rds')
##

#compile all "flagged"
jj1<-unique(rbind(aa1,bb1,cc1,dd1,ee1,ff1,gg1,hh1,ii1))
#remove them from GB.hab.filt
GB.hab.filt2<-setdiff(GB.hab.filt,jj1[,c(1:19)])

# remove additional ones that may have been missed before 
GB.hab.f11<-data.frame(GB.hab.filt2[which(GB.hab.filt2$habitat %in% (paste(GB.hab.filt2$GBIF.plant, ".", sep =""))),])
GB.hab.filt3<-GB.hab.filt2[which(!rownames(GB.hab.filt2) %in% rownames(GB.hab.f11)),]

#there are many that have this exact string as habitat field.

GB.hab.f12<-grep("Raised bog community with Pinus sylvestris, dwarfshrubs and dominating S. fuscum.", GB.hab.filt3$habitat)
GB.hab.f13<-grep("Old coniferous forests with dominance of P. sibirica", GB.hab.filt3$habitat)
GB.hab.filt4<-GB.hab.filt3[-c(GB.hab.f12,GB.hab.f13),]


#repeat gnffinder Name step ONLY for assocOr and assocTax
#asOr
GB.asOr.na<- GBIF.hits.sp %>% filter(is.na(habitat) & is.na(associatedTaxa))
GB.asOr.gnf<- data.frame(GB.asOr.na[GB.asOr.na$associatedOrganisms %in% gnf$Name,])
GB.asOr.filt<-GB.asOr.na[which(!rownames(GB.asOr.na) %in% rownames(GB.asOr.gnf)),]
#a lot of them just say "GENUS spec." in the associated organism field. Remove these.
GB.asOr.f1<-data.frame(GB.asOr.filt[which(GB.asOr.filt$associatedOrganisms %in% paste(GB.asOr.filt$GBIF.plant, ".", sep="")),])
GB.asOr.f2<-data.frame(GB.asOr.filt[which(GB.asOr.filt$associatedOrganisms %in% paste("Matrix:", GB.asOr.filt$GBIF.plant)),])
GB.asOr.rem<-rbind(GB.asOr.f1,GB.asOr.f2)
GB.asOr.filt.fin<-GB.asOr.filt[which(!rownames(GB.asOr.filt) %in% rownames(GB.asOr.rem)),]

#asTax
GB.asTa.na<- GBIF.hits.sp %>% filter(is.na(associatedOrganisms) & is.na(habitat))
GB.asTa.gnf<- data.frame(GB.asTa.na[GB.asTa.na$associatedTaxa %in% gnf$Name,])
GB.asTa.filt<-GB.asTa.na[which(!rownames(GB.asTa.na) %in% rownames(GB.asTa.gnf)),]

#COMBINE ALL 
GBIF.all.fields.filt<-rbind(GB.hab.filt4, GB.asOr.filt.fin, GB.asTa.filt)

############ Filter to everything that definitely grows on a plant host (ignoring "dead, decay, fallen, etc" search for now)
GBIF.on.plant<-GBIF.all.fields.filt %>%
  filter_at(vars(associatedTaxa,habitat,associatedOrganisms),any_vars(grepl(
    'host|infect| on |^on |"on |[on] |0n |parasit|leaves|leaf|foglie|twig|stem|huesped| sur |^sur |"sur |sulle |sulla | auf |^auf |"auf | ad |^ad |"ad | ^in |in |"in |^an |"an | an |paa |lebenden|anther| bark|ritidoma|branch|
      spot |spots |splotch|rust|streak|causing|inflorescence|fruit |fruits |pathogen|petiole|pod |pods |root |roots|cone|cones|utrículo|utricle|covering |culm|ovar|foli|needle|acícula|caulibus|caule|cortic|in den |flower|
      seed|pe frunze de|pe ramuri de|pe tulpini de|rama|smut |sobre |super |supra | sub |^sub |blättern|blattern|blad |stängeln|kapseln|stam |stamme |trunk|tronc|troco|truncus|living|live |vivant|sôbre|
      Wurzeln|Fruchtstielen|Schoten|Halmen|Rinde | vine|scape| head|stalk|frond|kvist|Balttsticlen|wedel|ramis|caudex|caudice|viva|rama|Epífita|hospedero|encontrado sobre|stolpe |
      Forófito de|onAbies|endophyt|^on:|substr|^på | på |^на | на | bladeren | Op |^Op |^aan | aan |^nos |sobre|substrat|sustrato| ^en |en |"en |bark | коре| ствол| лист| ветв| живы| корнях | корен | стеблях|
     плод| хвое | хвоя | завязи|nährpflanze', (.), ignore.case = TRUE)))
saveRDS(GBIF.on.plant, "OSF/GBIF/GBIF.on.plant.rds")

### Optional
#Go through those that may have good info and revise previous command as needed.
GBIF.not.on.plant<-setdiff(GBIF.all.fields.filt,GBIF.on.plant)
saveRDS(GBIF.not.on.plant, "intermediate_rds_csv/GBIF/GBIF.not.on.plant.rds")
#use Google translate

#### Collect fungi on dead plants and miscellaneous (e.g., soil, base of tree).

bajo<-which(c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa) %in% (paste("bajo", GBIF.on.plant$GBIF.plant)))
humus<-grep('humus', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
soil<-grep('soil', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
suelo<-grep('suelo', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
litter<-grep('litter', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
under<-grep('under', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
dead<-grep('dead', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
decompos<-grep('decompos', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
stump<-grep('stump', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
fallen<-grep('fallen', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
lo<-grep('^log', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
logg<-grep(' log', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
dying<-grep('dying', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
muert<-grep('muert', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
rott<-grep('rott', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
decay<-grep('decay', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
duff<-grep('duff', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
ground<-grep('ground', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
emortu<-grep('emortu', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
baseoftree<-grep('base of tree', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
hypogeous<-grep('hypogeous', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
dung<-grep('dung', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
putr<-grep('putr.', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
dod<-grep('død', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
charred<-grep('charred', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
burned<-grep('burned', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
burnt<-grep('burnt', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
sticks<-grep('sticks', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
voet<-grep('voet', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
wood<-grep('wood', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
ved<-grep('ved ', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
tocon<-grep('tocón', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
stubb<-grep('stubb', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
lignicolous<-grep('lignicolous', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
snag<-grep('snag', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)
madera<-grep('madera', c(GBIF.on.plant$habitat, GBIF.on.plant$associatedOrganisms, GBIF.on.plant$associatedTaxa), ignore.case = TRUE)

#For "on plant tissue" (not final dataset), keep several of those we are removing for the living database.
GBIF.on.plant.remove<-GBIF.on.plant[unique(c(soil, litter, under, duff, dung, ground, bajo, suelo, humus)),]
GBIF.on.plant.FIN<-GBIF.on.plant[-unique(c(soil, litter, under, duff, dung, ground, bajo, suelo, humus)),]

#Note: The following were removed for the GBIF.on.plant.FIN.rds file for the "On Plant" database, but were not initially removed here for the next step (preparation of the "Fungi on Living Plant Tissue database"
# ‘med bakkekontakt’, ‘på bakken’, ‘Growing on bajo’, ‘ on rock\\’, ‘on granit’, ‘låg’ , ‘dead’, ‘liggande’, ‘låg., ‘:låg’, ‘on moss’, lying trunk’, terricolous’, ‘growing on rocks’, ‘on downed’, ‘ mort\\.’, ‘emort’, ‘död’, ‘ sur sol ‘, ‘base of old ‘, ‘on limestone outcrops’,  ‘on sandstone’, ‘on stone’,  ‘, base of trunk’, ‘; on rock’, ‘\\. on rock’, ‘on root’.

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# II. Here we create GBIF.on.plant.fin, which will be used as raw data for the "On plant tissue" database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

saveRDS(GBIF.on.plant.FIN, 'intermediate_rds_csv/GBIF/GBIF.on.plant.FIN.rds')

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


deadplant<-unique(c(humus, bajo,suelo,dead,decompos,stump,fallen,lo,logg,dying,muert,rott,decay,soil,litter,under,duff,ground,emortu,baseoftree,hypogeous,dung,putr,dod, charred, burned, burnt, sticks, 
                    voet, wood, ved, tocon, stubb, lignicolous, snag, madera))

### Not on living plants and some on soil.
GBIF.deadplant<-GBIF.on.plant[deadplant,]
saveRDS(GBIF.deadplant, "~/data2/LPs/rds/GBIF.deadplant.rds")

#### Here is the rough dataframe for living plants
GBIF.liveplant<-GBIF.on.plant[-deadplant,]
saveRDS(GBIF.liveplant, "OSF/GBIF/GBIF.liveplant.rds")


### There are too many live plant entries here. Edit some out. Specifically, those already found from Mycoportal.
myco.edit.liveplant<-read.csv("OSF/Mycoportal/Mycop.edit.liveplant.csv")
# generate list
myco.36<-myco.edit.liveplant[which(nchar(as.character(myco.edit.liveplant$occurrenceID))==36),]
GBIF.myc.filter<- GBIF.liveplant[which(!GBIF.liveplant$occurrenceID %in% myco.36$occurrenceID),]
saveRDS(GBIF.myc.filter, "OSF/GBIF.myc.filter.rds")

#### Filter out fungi (Index Fungorum source)
# NOTE: DOUBLE CHECK GBIF.myc.filter = G.myc.filter
hosts.G.myc.filt<-as.character(unique(G.myc.filter$GBIF.plant))
plnamresolv.F<-gnr_resolve(hosts.G.myc.filt[1:999], preferred_data_sources = 5)
plnamresolv.F1<-gnr_resolve(hosts.G.myc.filt[1000:2000], preferred_data_sources = 5)
plnamresolv.F2<-gnr_resolve(hosts.G.myc.filt[2001:3000], preferred_data_sources = 5)
plnamresolv.F3<-gnr_resolve(hosts.G.myc.filt[3001:4000], preferred_data_sources = 5)
plnamresolv.F4<-gnr_resolve(hosts.G.myc.filt[4001:5000], preferred_data_sources = 5)
plnamresolv.F5<-gnr_resolve(hosts.G.myc.filt[5001:6000], preferred_data_sources = 5)
plnamresolv.F6<-gnr_resolve(hosts.G.myc.filt[6001:7000], preferred_data_sources = 5)
plnamresolv.F7<-gnr_resolve(hosts.G.myc.filt[7001:8000], preferred_data_sources = 5)
plnamresolv.F8<-gnr_resolve(hosts.G.myc.filt[8001:9000], preferred_data_sources = 5)
plnamresolv.F9<-gnr_resolve(hosts.G.myc.filt[9001:10000], preferred_data_sources = 5)
plnamresolv.F10<-gnr_resolve(hosts.G.myc.filt[10001:11000], preferred_data_sources = 5)
plnamresolv.F11<-gnr_resolve(hosts.G.myc.filt[11001:12000], preferred_data_sources = 5)
plnamresolv.F12<-gnr_resolve(hosts.G.myc.filt[12001:13000], preferred_data_sources = 5)
plnamresolv.F13<-gnr_resolve(hosts.G.myc.filt[13001:14000], preferred_data_sources = 5)
plnamresolv.F14<-gnr_resolve(hosts.G.myc.filt[14001:15000], preferred_data_sources = 5)
plnamresolv.F15<-gnr_resolve(hosts.G.myc.filt[15001:16000], preferred_data_sources = 5)
plnamresolv.F16<-gnr_resolve(hosts.G.myc.filt[16001:17000], preferred_data_sources = 5)
plnamresolv.F17<-gnr_resolve(hosts.G.myc.filt[17001:18000], preferred_data_sources = 5)
plnamresolv.F18<-gnr_resolve(hosts.G.myc.filt[18001:19000], preferred_data_sources = 5)
plnamresolv.F19<-gnr_resolve(hosts.G.myc.filt[19001:20000], preferred_data_sources = 5)
plnamresolv.F20<-gnr_resolve(hosts.G.myc.filt[20001:21000], preferred_data_sources = 5)
plnamresolv.F21<-gnr_resolve(hosts.G.myc.filt[21001:22000], preferred_data_sources = 5)
plnamresolv.F22<-gnr_resolve(hosts.G.myc.filt[22001:23000], preferred_data_sources = 5)
plnamresolv.F23<-gnr_resolve(hosts.G.myc.filt[23001:24000], preferred_data_sources = 5)
plnamresolv.F24<-gnr_resolve(hosts.G.myc.filt[24001:25000], preferred_data_sources = 5)
plnamresolv.F25<-gnr_resolve(hosts.G.myc.filt[25001:26000], preferred_data_sources = 5)
plnamresolv.F26<-gnr_resolve(hosts.G.myc.filt[26001:27000], preferred_data_sources = 5)
plnamresolv.F27<-gnr_resolve(hosts.G.myc.filt[27001:28000], preferred_data_sources = 5)
plnamresolv.F28<-gnr_resolve(hosts.G.myc.filt[28001:29000], preferred_data_sources = 5)
plnamresolv.F29<-gnr_resolve(hosts.G.myc.filt[29001:30000], preferred_data_sources = 5)
plnamresolv.F30<-gnr_resolve(hosts.G.myc.filt[30001:31000], preferred_data_sources = 5)
plnamresolv.F31<-gnr_resolve(hosts.G.myc.filt[31001:32000], preferred_data_sources = 5)
plnamresolv.F32<-gnr_resolve(hosts.G.myc.filt[32001:33000], preferred_data_sources = 5)
plnamresolv.F33<-gnr_resolve(hosts.G.myc.filt[33001:34000], preferred_data_sources = 5)
plnamresolv.F34<-gnr_resolve(hosts.G.myc.filt[34001:34868], preferred_data_sources = 5)

plnamresolv.F.all<-rbind.data.frame(plnamresolv.F, plnamresolv.F1, plnamresolv.F2, plnamresolv.F3, plnamresolv.F4, plnamresolv.F5, plnamresolv.F6, 
                                                plnamresolv.F7, plnamresolv.F8, plnamresolv.F9, plnamresolv.F10, plnamresolv.F11, plnamresolv.F12, plnamresolv.F13, 
                                                plnamresolv.F14, plnamresolv.F15, plnamresolv.F16, plnamresolv.F17, plnamresolv.F18, plnamresolv.F19, plnamresolv.F20, 
                                                plnamresolv.F21, plnamresolv.F22, plnamresolv.F23, plnamresolv.F24, plnamresolv.F25, plnamresolv.F26, plnamresolv.F27, 
                                                plnamresolv.F28, plnamresolv.F29, plnamresolv.F30, plnamresolv.F31, plnamresolv.F32, plnamresolv.F33, plnamresolv.F34)
saveRDS(plnamresolv.F.all, "intermediate_rds_csv/GBIF/plnamresolv.Fungi.all.rds")                                     

GBIF.myc.filter$match.number<-match(as.character(GBIF.myc.filter$GBIF.plant), plnamresolv.F.all$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
GBIF.myc.filter$match.name<-print(plnamresolv.F.all$matched_name[GBIF.myc.filter$match.number[1:587872]])
#dataframe without fungi
no.fung<-GBIF.myc.filter[is.na(GBIF.myc.filter$match.name),]
#dataframe of fungi
fung<-GBIF.myc.filter[!is.na(GBIF.myc.filter$match.name),]
#remaining hosts
hosts.nf.liveplant<-as.character(unique(no.fung$GBIF.plant))

#### Assign host names (International Plant Names Database)
plnamresolv.TPL<-gnr_resolve(hosts.nf.liveplant[1:999], preferred_data_sources = 167)
plnamresolv.TPL1<-gnr_resolve(hosts.nf.liveplant[1000:2000], preferred_data_sources = 167)
plnamresolv.TPL2<-gnr_resolve(hosts.nf.liveplant[2001:3000], preferred_data_sources = 167)
plnamresolv.TPL3<-gnr_resolve(hosts.nf.liveplant[3001:4000], preferred_data_sources = 167)
plnamresolv.TPL4<-gnr_resolve(hosts.nf.liveplant[4001:5000], preferred_data_sources = 167)
plnamresolv.TPL5<-gnr_resolve(hosts.nf.liveplant[5001:6000], preferred_data_sources = 167)
plnamresolv.TPL6<-gnr_resolve(hosts.nf.liveplant[6001:7000], preferred_data_sources = 167)
plnamresolv.TPL7<-gnr_resolve(hosts.nf.liveplant[7001:8000], preferred_data_sources = 167)
plnamresolv.TPL8<-gnr_resolve(hosts.nf.liveplant[8001:9000], preferred_data_sources = 167)
plnamresolv.TPL9<-gnr_resolve(hosts.nf.liveplant[9001:10000], preferred_data_sources = 167)
plnamresolv.TPL10<-gnr_resolve(hosts.nf.liveplant[10001:11000], preferred_data_sources = 167)
plnamresolv.TPL11<-gnr_resolve(hosts.nf.liveplant[11001:12000], preferred_data_sources = 167)
plnamresolv.TPL12<-gnr_resolve(hosts.nf.liveplant[12001:13000], preferred_data_sources = 167)
plnamresolv.TPL13<-gnr_resolve(hosts.nf.liveplant[13001:14000], preferred_data_sources = 167)
plnamresolv.TPL14<-gnr_resolve(hosts.nf.liveplant[14001:15000], preferred_data_sources = 167)
plnamresolv.TPL15<-gnr_resolve(hosts.nf.liveplant[15001:16000], preferred_data_sources = 167)
plnamresolv.TPL16<-gnr_resolve(hosts.nf.liveplant[16001:17000], preferred_data_sources = 167)
plnamresolv.TPL17<-gnr_resolve(hosts.nf.liveplant[17001:18000], preferred_data_sources = 167)
plnamresolv.TPL18<-gnr_resolve(hosts.nf.liveplant[18001:19000], preferred_data_sources = 167)
plnamresolv.TPL19<-gnr_resolve(hosts.nf.liveplant[19001:20000], preferred_data_sources = 167)
plnamresolv.TPL20<-gnr_resolve(hosts.nf.liveplant[20001:21000], preferred_data_sources = 167)
plnamresolv.TPL21<-gnr_resolve(hosts.nf.liveplant[21001:22000], preferred_data_sources = 167)
plnamresolv.TPL22<-gnr_resolve(hosts.nf.liveplant[22001:23000], preferred_data_sources = 167)
plnamresolv.TPL23<-gnr_resolve(hosts.nf.liveplant[23001:24000], preferred_data_sources = 167)
plnamresolv.TPL24<-gnr_resolve(hosts.nf.liveplant[24001:25000], preferred_data_sources = 167)
plnamresolv.TPL25<-gnr_resolve(hosts.nf.liveplant[25001:26000], preferred_data_sources = 167)
plnamresolv.TPL26<-gnr_resolve(hosts.nf.liveplant[26001:27000], preferred_data_sources = 167)
plnamresolv.TPL27<-gnr_resolve(hosts.nf.liveplant[27001:28000], preferred_data_sources = 167)
plnamresolv.TPL28<-gnr_resolve(hosts.nf.liveplant[28001:29000], preferred_data_sources = 167)
plnamresolv.TPL29<-gnr_resolve(hosts.nf.liveplant[29001:30000], preferred_data_sources = 167)
plnamresolv.TPL30<-gnr_resolve(hosts.nf.liveplant[30001:31125], preferred_data_sources = 167)
plnamresolv.TPL.all<-rbind.data.frame(plnamresolv.TPL, plnamresolv.TPL1, plnamresolv.TPL2, plnamresolv.TPL3, plnamresolv.TPL4, plnamresolv.TPL5, 
                                            plnamresolv.TPL6, plnamresolv.TPL7, plnamresolv.TPL8, plnamresolv.TPL9, plnamresolv.TPL10, plnamresolv.TPL11, 
                                            plnamresolv.TPL12, plnamresolv.TPL13, plnamresolv.TPL14, plnamresolv.TPL15, plnamresolv.TPL16, plnamresolv.TPL17, 
                                            plnamresolv.TPL18, plnamresolv.TPL19, plnamresolv.TPL20, plnamresolv.TPL21, plnamresolv.TPL22, plnamresolv.TPL23, 
                                            plnamresolv.TPL24, plnamresolv.TPL25, plnamresolv.TPL26, plnamresolv.TPL27, plnamresolv.TPL28, plnamresolv.TPL29, plnamresolv.TPL30)
saveRDS(plnamresolv.TPL.all, "intermediate_rds_csv/GBIF/plnamresolv.TPL.all.rds")  

#reset match name and match number columns for TPL step.
no.fung<-no.fung[,1:19]
no.fung$match.number<-match(no.fung$GBIF.plant, plnamresolv.TPL.all$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
no.fung$match.name<-print(plnamresolv.TPL.all$matched_name[no.fung$match.number[1:567296]])

#remaining ones
not.found.TPL<-no.fung[is.na(no.fung$match.name),]
TPL.G<-no.fung[!is.na(no.fung$match.name),]

#use Tropicos 
host.nf2.LP<-as.character(unique(not.found.TPL$GBIF.plant))
plnamresolv.Tr<-gnr_resolve(host.nf2.LP[1:1000], preferred_data_sources = 165)
plnamresolv.Tr1<-gnr_resolve(host.nf2.LP[1001:2000], preferred_data_sources = 165)
plnamresolv.Tr2<-gnr_resolve(host.nf2.LP[2001:3000], preferred_data_sources = 165)
plnamresolv.Tr3<-gnr_resolve(host.nf2.LP[3001:4000], preferred_data_sources = 165)
plnamresolv.Tr4<-gnr_resolve(host.nf2.LP[4001:5000], preferred_data_sources = 165)
plnamresolv.Tr5<-gnr_resolve(host.nf2.LP[5001:6000], preferred_data_sources = 165)
plnamresolv.Tr6<-gnr_resolve(host.nf2.LP[6001:7000], preferred_data_sources = 165)
plnamresolv.Tr7<-gnr_resolve(host.nf2.LP[7001:8201], preferred_data_sources = 165)
plnamresolv.Tr.all<-rbind.data.frame(plnamresolv.Tr,plnamresolv.Tr1,plnamresolv.Tr2,plnamresolv.Tr3,plnamresolv.Tr4,
                                     plnamresolv.Tr5,plnamresolv.Tr6,plnamresolv.Tr7)
saveRDS(plnamresolv.Tr.all, "intermediate_rds_csv/GBIF/plnamresolv.Tr.all.rds") 

not.found.TPL<-not.found.TPL[,1:19]
not.found.TPL$match.number<-match(not.found.TPL$GBIF.plant, plnamresolv.Tr.all$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
not.found.TPL$match.name<-print(plnamresolv.Tr.all$matched_name[not.found.TPL$match.number[1:100664]])

#remaining ones
not.found.Tr<-not.found.TPL[is.na(not.found.TPL$match.name),]



nf.Tr.hosts<-as.character(unique(not.found.Tr$GBIF.plant))
#found ones

Tr<-not.found.Tr[!is.na(not.found.Tr$match.name),]

plnamresolv.US<-gnr_resolve(nf.Tr.hosts[1:1000], preferred_data_sources = 150)
plnamresolv.US1<-gnr_resolve(nf.Tr.hosts[1001:2000], preferred_data_sources = 150)
plnamresolv.US2<-gnr_resolve(nf.Tr.hosts[2001:3000], preferred_data_sources = 150)
plnamresolv.US3<-gnr_resolve(nf.Tr.hosts[3001:4000], preferred_data_sources = 150)
plnamresolv.US4<-gnr_resolve(nf.Tr.hosts[4001:5000], preferred_data_sources = 150)
plnamresolv.US5<-gnr_resolve(nf.Tr.hosts[5001:6000], preferred_data_sources = 150)
plnamresolv.US6<-gnr_resolve(nf.Tr.hosts[6001:7067], preferred_data_sources = 150)
plnamresolv.US.all<-rbind.data.frame(plnamresolv.US,plnamresolv.US1,plnamresolv.US2,plnamresolv.US3,plnamresolv.US4,plnamresolv.US5,plnamresolv.US6)
saveRDS(plnamresolv.US.all, "intermediate_rds_csv/GBIF/plnamresolv.US.all.rds") 

not.found.Tr<-not.found.Tr[,1:19]
not.found.Tr$match.number<-match(not.found.Tr$GBIF.plant, plnamresolv.US.all$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
not.found.Tr$match.name<-print(plnamresolv.US.all$matched_name[not.found.Tr$match.number[1:94263]])
not.found.US<-not.found.Tr[is.na(not.found.Tr$match.name),]
#found ones

US<-not.found.US[!is.na(not.found.US$match.name),]

######Now compile a list with this host info (TPL, USDA, Tropicos). There will be more to add later.
#GBIF
#Prepare THE PLANT LIST dataframe
TPL.G$new.ID<-"NA"
TPL.G$namecheck.db<-"Intl.plant.names.index"
TPL.G$score<-print(plnamresolv.TPL.all$score[TPL.G$match.number[1:466632]])
GBIF.list.1<-TPL.G[!is.na(TPL.G$match.name),]
#Prepare TROPICOS dataframe
Tr$new.ID<-"NA"
Tr$namecheck.db<-"Tropicos"
Tr$score<-print(plnamresolv.Tr.all$score[Tr$match.number[1:236]])
GBIF.list.2<-Tr[!is.na(Tr$match.name),]

#Prepare usda dataframe
US$new.ID<-"NA"
US$namecheck.db<-"USDA.plant.db"
US$score<-print(plnamresolv.US.all$score[US$match.number[1:1284]])
GBIF.list.3<-US[!is.na(US$match.name),]
list.taxize.GBIF<-rbind.data.frame(GBIF.list.1, GBIF.list.2, GBIF.list.3)
#prepare full file. this can be included as supplementary or for reference.
saveRDS(list.taxize.GBIF, "OSF/GBIF/list.taxize.GBIF.rds")


####Resume finding remaining entries
#eliminate animal entries (Index Animalium)
nf.US.hosts<-as.character(unique(not.found.US$GBIF.plant))
plnamresolv.Z<-gnr_resolve(nf.US.hosts[1:1000], preferred_data_sources = 183)
plnamresolv.Z1<-gnr_resolve(nf.US.hosts[1001:2000], preferred_data_sources = 183)
plnamresolv.Z2<-gnr_resolve(nf.US.hosts[2001:3000], preferred_data_sources = 183)
plnamresolv.Z3<-gnr_resolve(nf.US.hosts[3001:4000], preferred_data_sources = 183)
plnamresolv.Z4<-gnr_resolve(nf.US.hosts[4001:5000], preferred_data_sources = 183)
plnamresolv.Z5<-gnr_resolve(nf.US.hosts[5001:6000], preferred_data_sources = 183)
plnamresolv.Z6<-gnr_resolve(nf.US.hosts[6001:7028], preferred_data_sources = 183)
plnamresolv.Z.all<-rbind.data.frame(plnamresolv.Z,plnamresolv.Z1,plnamresolv.US2,plnamresolv.Z3,plnamresolv.Z4,plnamresolv.Z5,plnamresolv.Z6)
saveRDS(plnamresolv.Z.all, "intermediate_rds_csv/GBIF/plnamresolv.Z.all.rds") 

not.found.US<-not.found.US[,1:19]
not.found.US$match.number<-match(not.found.US$GBIF.plant, plnamresolv.Z.all$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
not.found.US$match.name<-print(plnamresolv.Z.all$matched_name[not.found.US$match.number[1:94027]])
not.found.Z<-not.found.US[is.na(not.found.US$match.name),]

#split these based on if the GBIF plant is coded as abbreviated genus+specific epithet. Taxonstand will not be helpful for these entries.
not.found.Z.abbrev<-not.found.Z[grep("^.\\.", not.found.Z$GBIF.plant),]
not.found.Z.not.abbrev<-not.found.Z[which(!not.found.Z$GBIF.plant %in% not.found.Z.abbrev$GBIF.plant),]

#There are a bunch of entries where GBIF plant has a different genus first letter. Abbreviated sp names might be helpful, but only if they honestly represent the habitat/assOr/asTa.
correct.nf.Z.ab.hab<-not.found.Z.abbrev[mapply(grepl, not.found.Z.abbrev$GBIF.plant, not.found.Z.abbrev$habitat, fixed=TRUE), ]
correct.nf.Z.ab.asta<-not.found.Z.abbrev[mapply(grepl, not.found.Z.abbrev$GBIF.plant, not.found.Z.abbrev$associatedTaxa, fixed=TRUE), ]
correct.nf.Z.ab.asor<-not.found.Z.abbrev[mapply(grepl, not.found.Z.abbrev$GBIF.plant, not.found.Z.abbrev$associatedOrganisms, fixed=TRUE), ]
all.correct.nf.Z.ab<-unique(rbind(correct.nf.Z.ab.hab, correct.nf.Z.ab.asta, correct.nf.Z.ab.asor))

#all marked "incorrect"
incorrect.nf.Z.ab<-anti_join(not.found.Z.abbrev, all.correct.nf.Z.ab)
#these can be saved to verify at the end, but it is doubtful they will be helpful.
write.csv(incorrect.nf.Z.ab, "intermediate_rds_csv/GBIF/incorrect.nf.Z.ab.csv")


###III Taxonstand filtering
#make list of host names for taxonstand.
nf.for.Taxonstand<-as.character(unique(not.found.Z.not.abbrev$GBIF.plant))
#Use Taxonstand to try to interpret some more (from The Plant List). This takes some time.
nf.TPL.taxons<-TPL(nf.for.Taxonstand)
saveRDS(nf.TPL.taxons, "intermediate_rds_csv/GBIF/nf.TPL.taxons.rds")
#save successful Taxonstand hits
nf.Taxon.hits.G<-filter(nf.TPL.taxons, New.Taxonomic.status != "")

#the rest
nf.Taxon.fail.G<-filter(nf.TPL.taxons, New.Taxonomic.status == "")
#Filter fungal genera
#list fung genera (from Index Fungorum, before)
fung.gen.G<-unique(word(fung$GBIF.plant, start=1L))
nf.TfG<-filter(nf.Taxon.fail.G, !New.Genus %in% fung.gen.G)
# Correct genus names. If Genus names are corrected, the TPL and other tools perform better. These are primarily modified latin names, but some include easy-to-spot typos. All apparent genera typos with >1 entry were corrected. 
nf.TfG$sing.genus<-gsub("Populi","Populus",gsub("Abietis","Abies",gsub("Aceris","Acer",
                              gsub("Aconiti","Aconitum",(gsub("Agropyrum","Agropyron",gsub("Agrostidis","Agrostis", gsub("Allii","Allium", gsub("Alni","Alnus",
                              gsub("Baccharidis", "Baccharis", gsub("Bromi","Bromus",gsub("Calamogrostidis","Calamogrostis",gsub("Caricis","Carex", gsub("Carpini","Carpinum",
                              gsub("Clerodendron","Clerodendrum", gsub("Corni","Cornus",gsub("Crataegi","Crataegus", gsub("Cytisi", "Cytisum", gsub("Delphini","Delphinium", 
                              gsub("Desmodii","Desmodium",gsub("Dipsaci", "Dipsacus",gsub("Elymi","Elymus",gsub("Euonymi","Euonymus", gsub("Eupatorii","Eupatorium",
                              gsub("Evonymi","Euonymus",gsub("Evonymus","Euonmyus",gsub("Fagi","Fagus",gsub("Fici","Ficus",gsub("Fraxini","Fraxinus",gsub("Galii","Galium",
                              gsub("Gei","Geum",gsub("Gelsemii","Gelsemium",gsub("Geranii","Gernaium",gsub("Heraclei","Heracleum",gsub("Hordei","Hordeum",gsub("Humuli","Humulus",
                              gsub("Ilicis","Ilex",gsub("Impatientis","Impatiens",gsub("Iridis","Iris",gsub("Juglandis","Juglans",gsub("Junci","Juncus",gsub("Juniperi","Juniperus",
                              gsub("Lappae","Lappa",gsub("Lathyri","Lathyrus",gsub("Lauri","Laurus",gsub("Leontodontis","Leontodon",gsub("Ligustici","Ligusticum",gsub("Lupini","Lupinus",
                              nf.TfG$New.Genus))))))))))))))))))))))))))))))))))))))))))))))))
nf.TfG$sing.genus<-gsub("Lycii","Lycium",gsub("Mori","Morus",gsub("Nasturtii","Nasturtium",gsub("Nerii","Nerium",gsub("Olerodendron","Clerodendron",
                              gsub("Ononidis","Ononis",gsub("Orobi","Orobus",gsub("Osmorrhizae","Osmorhiza",gsub("Panici","Panicum",gsub("Pentstamon","Penstemon",gsub("Pentstemonis","Penstemon",
                              gsub("Peucedani","Peucedanum",gsub("Pini","Pinus",gsub("Plantaginis","Plantago",gsub("Polygoni","Polygonum",gsub("Polypodium","Polypodii",gsub("Populi","Populus",
                              gsub("Pruni","Prunus",gsub("Pteridis","Pteridium",gsub("Pyri","Pyrus",gsub("Rhamni","Rhamnus",gsub("Rhois","Rhus",gsub("Rhoidis","Rhus",gsub("Ribis","Ribes",
                              gsub("Roripa","Rorippa",gsub("Rubi","Rubus",gsub("Rumicis","Rumex",gsub("Salicis","Salex",gsub("Sambuci","Sambucus",gsub("Scirpi","Scirpus",gsub("Senecionis","Senecio",
                              gsub("Silphinium","Silphium",gsub("Sisymbrii","Sisymbrium",gsub("Smilacis","Smilax",gsub("Solani","Solanum",gsub("Sonchi","Sonchus",gsub("Symphyti","Symphytum",gsub("Tami","Tamus",
                              gsub("Tanaceti","Tanacetum",gsub("Teucrii","Teucrium",gsub("Thalictri","Thalictrum",gsub("Thesii","Theisum",gsub("Trifolii","Trifolium",gsub("Tritici","Triticum",gsub("Ulmi","Ulmus",
                              gsub("Vaccinii","Vaccinium",gsub("Veronicae","Veronica",gsub("Viburni","Viburnum", gsub("Zizyphus","Ziziphus",nf.TfG$sing.genus)))))))))))))))))))))))))))))))))))))))))))))))))

#compile corrected genera, paste with original specific epithets.
nf.TfG$Corr.name<-paste(nf.TfG$sing.genus, nf.TfG$New.Species)
#rerun TPL on corrected data, hope for less manual work!
re.TPL.G<-TPL(nf.TfG$Corr.name)
saveRDS(re.TPL.G, "intermediate_rds_csv/GBIF/re.TPL.G.rds")
#successful taxon hits
re.TPL.taxon.hits.G<-filter(re.TPL.G, New.Taxonomic.status != "")

###### IV prepare TPL Taxonstand output.
TPL.list.G<-rbind.data.frame(nf.Taxon.hits.G, re.TPL.taxon.hits.G)
TPL.list.G$match.name<-paste(TPL.list.G$New.Genus, TPL.list.G$New.Hybrid.Marker, TPL.list.G$New.Species,TPL.list.G$New.Authority)
TPL.list.G$match.name<-gsub("  "," ",TPL.list.G$match.name)
TPL.list.G$match.name<-gsub("  "," ",TPL.list.G$match.name)

#add to match.name, new record ID, database name 
TPL.list.G$namecheck.db<-"The.plant.list"

#follow same format as list.taxize.GBIF
not.found.Z.not.abbrev<-not.found.Z.not.abbrev[,1:19]
not.found.Z.not.abbrev$match.number<-match(not.found.Z.not.abbrev$GBIF.plant, TPL.list.G$Taxon, nomatch = NA_integer_, incomparables = NULL)
not.found.Z.not.abbrev$match.name<-print(TPL.list.G$match.name[not.found.Z.not.abbrev$match.number[1:58437]])
not.found.Z.not.abbrev$new.ID<-print(TPL.list.G$New.ID[not.found.Z.not.abbrev$match.number[1:58437]])
list.TPL.GBIF<-not.found.Z.not.abbrev[!is.na(not.found.Z.not.abbrev$match.name),]
#remove extra spaces

list.TPL.GBIF$namecheck.db<-"The.plant.list"
list.TPL.GBIF$score<-"NA"
saveRDS(list.TPL.GBIF, "intermediate_rds_csv/GBIF/list.TPL.GBIF.rds")

##### V WorldFlora Filtering
#taxon fails
re.TPL.taxon.fails.G<-filter(re.TPL.G, New.Taxonomic.status == "")
#WorldFlora pkg
#WFO.download(save.dir = "~/OneDrive - UBC/PhD.Year.Three", WFO.remember = TRUE)
WFO.remember()
re.WFO.G<-WFO.match(re.TPL.taxon.fails.G$Taxon, WFO.data=WFO.data)

#saveRDS(re.WFO.G, "intermediate_rds_csv/GBIF/re.WFO.G.rds")
#write.csv for re.WFO.G as "/re.WFO.G.csv", row.names=FALSE)

#manually edit. 
#editing:
#I scanned through everything looking for obvious errors and ambiguous abbreviations (e.g., under Abietum, Ph. arvense). 
#I sorted by taxonrank and looked through the "NA" for lots of convincing species names, many were beetles.
#I skipped to "species" taxonrank and filtered by fuzzy dist and manually checked. If there were duplicates, I chose one of the "Accepted" taxonomic status entries.
#Then, went through all unknowns, and deleted all except those that I recognized (kept those with "-etum" ending.)


#I started typing all genera that did not come up with old or accepted alternatives in world flora. If there was a clear solution, I added this to the "manual" column. 
#I checked all entries to make sure they made sense, and searched world flora if only a genus determination was made.
#All entries that were clearly incorrect were removed. Verged on over-removal rather than the alternative.
#removed entries where single old name forms 2 currently accepted names (e.g., Centaurea glauca is now Amberboa glauca/Centaurea calocephala) UNLESS there is accepted name that is the same but correct spelling.
#deleted entries where genus clearly did not match.
#deleted entries where tie was convincing ("Lachia pulchella" (sic) is Lawia pulchella or Lechea pulchella)
#changed variety level classification to species level.
#If only genus was detected, I looked in World Flora to see if there that genus had an entry with a similar specific epithet, and added this in "manual".

#Saved as filt.re.WFO.G.csv

re.WFO.G.filt<-read.csv("intermediate_rds_csv/GBIF/filt.re.WFO.G.csv")
re.WFO.G.filt1<-re.WFO.G.filt[which(unique(re.WFO.G.filt$spec.name) %in% re.WFO.G.filt$spec.name),]


#part 3: prepare WorldFlora output
#include all data fields
re.WFO.G.filt.full<-re.WFO.G[which(re.WFO.G$spec.name %in% re.WFO.G.filt1$spec.name),]
re.WFO.G.filt.full1<-re.WFO.G.filt.full[which(unique(re.WFO.G.filt.full$spec.name) %in% re.WFO.G.filt.full$spec.name),]
re.WFO.G.filt.full1$match.number<-match(re.WFO.G.filt.full1$spec.name, re.WFO.G.filt1$spec.name, nomatch = NA_integer_, incomparables = NULL)
re.WFO.G.filt.full1$manual.name<-print(re.WFO.G.filt1$manual.name[re.WFO.G.filt.full1$match.number[1:1283]])
re.WFO.G.filt.full1$manual.auth<-print(re.WFO.G.filt1$manual.auth[re.WFO.G.filt.full1$match.number[1:1283]])
re.WFO.G.filt.full1$match.name<-print(re.WFO.G.filt1$scientificName[re.WFO.G.filt.full1$match.number[1:1283]])
re.WFO.G.filt.full1$match.auth<-print(re.WFO.G.filt1$scientificNameAuthorship[re.WFO.G.filt.full1$match.number[1:1283]])

#replace match name with manual name if a manual name is detected.
re.WFO.G.filt.full1$match.name[!re.WFO.G.filt.full1$manual.name==""]<-re.WFO.G.filt.full1$manual.name[!re.WFO.G.filt.full1$manual.name==""]
re.WFO.G.filt.full1$match.auth[!re.WFO.G.filt.full1$manual.auth==""]<-re.WFO.G.filt.full1$manual.auth[!re.WFO.G.filt.full1$manual.auth==""]
re.WFO.G.filt.full1$match.name<-paste(re.WFO.G.filt.full1$match.name,re.WFO.G.filt.full1$match.auth)

saveRDS(re.WFO.G.filt.full1, "intermediate_rds_csv/GBIF/re.WFO.G.filt.full1.rds")

#replace these in GBIF data
re.WFO.G.filt.full1$tplID[!re.WFO.G.filt.full1$manual.name==""]<-NA
filt.GBIF.WFO<-not.found.Z.not.abbrev[which(not.found.Z.not.abbrev$GBIF.plant %in% re.WFO.G.filt.full1$spec.name),]
filt.GBIF.WFO$match.number<-match(filt.GBIF.WFO$GBIF.plant, re.WFO.G.filt.full1$spec.name, nomatch = NA_integer_, incomparables = NULL)
filt.GBIF.WFO$match.name<-print(re.WFO.G.filt.full1$match.name[filt.GBIF.WFO$match.number[1:3571]])
filt.GBIF.WFO$new.ID<-print(re.WFO.G.filt.full1$tplID[filt.GBIF.WFO$match.number[1:3571]])
filt.GBIF.WFO$namecheck.db<-"WorldFlora"
filt.GBIF.WFO$score<-"NA"

### V Combine them all

On.plant.GBIF.rough<-rbind(list.taxize.GBIF, list.TPL.GBIF, filt.GBIF.WFO)
#the new.ID for the manually changed reWFO will not be there.
saveRDS(On.plant.GBIF.rough, "intermediate_rds_csv/GBIF/On.plant.GBIF.rough.rds")

#Filter out match name only to genus.
On.plant.GBIF.sp<-On.plant.GBIF.rough[grep(" ", On.plant.GBIF.rough$match.name),]

#more habitat filtration
match.n<-unique(On.plant.GBIF.sp$GBIF.plant)
con<-On.plant.GBIF.sp[grep("con ", On.plant.GBIF.sp$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in match.n){
  mm<-ifelse(nn<-grep(paste("con", i), con$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
s.con<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
aaa<-cbind(con[s.con$GBIF.entry,], s.con)
aaaa<-aaa[which(!grepl("sobre ", aaa$habitat, ignore.case=TRUE)),]
aaa1<-aaaa %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))
#remove them from On.plant.GBIF.sp
On.plant.GBIF.sp.2<-setdiff(On.plant.GBIF.sp,aaa1[,c(1:24)])

#filter bracket
bracket<-On.plant.GBIF.sp.2[grep("}", On.plant.GBIF.sp.2$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in match.n){
  mm<-ifelse(nn<-grep(paste(i, "}", sep = ""), bracket$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
s.bracket<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
bbb<-cbind(bracket[s.bracket$GBIF.entry,], s.bracket)
bbba<-bbb[which(!grepl("On \\{", bbb$habitat, ignore.case=TRUE)),]
bbb1<-bbba %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))
#remove them from On.plant.GBIF.sp
On.plant.GBIF.sp.3<-setdiff(On.plant.GBIF.sp.2,bbb1[,c(1:24)])

#filter start bracket
bracket2<-On.plant.GBIF.sp.3[grep("\\{", On.plant.GBIF.sp.3$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in match.n){
  mm<-ifelse(nn<-grep(paste("\\{", i, sep = ""), bracket2$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
s.bracket2<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
ccc<-cbind(bracket2[s.bracket2$GBIF.entry,], s.bracket2)
ccca<-ccc[which(!grepl("On \\{", ccc$habitat, ignore.case=TRUE)),]
ccc1<-ccca %>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))
#remove them from On.plant.GBIF.sp
On.plant.GBIF.sp.4<-setdiff(On.plant.GBIF.sp.3,ccc1[,c(1:24)])

#### NOTE: I don't think these (besides ' on -, see below) actually were removed here.
######

boulder<-On.plant.GBIF.sp.4[grep("on boulder", On.plant.GBIF.sp.4$habitat, ignore.case=TRUE),]
onnn<-On.plant.GBIF.sp.4[grep(" on ", On.plant.GBIF.sp.4$habitat, ignore.case=TRUE),]
onzone<-On.plant.GBIF.sp.4[grep(" on zone", On.plant.GBIF.sp.4$habitat, ignore.case=TRUE),]

#remove duplicates. Sort by longest first to get rid of duplicates like 'Abies spec' from 'Abies spectabilis' in GBIF.plant field.
#first, keep orig row.names.
On.plant.GBIF.sp.4$orig.row.names<-row.names(On.plant.GBIF.sp.4)
On.plant.GBIF.sp.5<- On.plant.GBIF.sp.4 %>%
  group_by(match.name, habitat, associatedTaxa, associatedOrganisms) %>%
  arrange(desc(nchar(as.character(GBIF.plant))), .by_group = TRUE) %>%
  distinct(gbifID, institutionID, collectionID, occurrenceID, associatedTaxa, organismID, associatedOrganisms, habitat, higherGeography, countryCode, scientificName, genus, speciesKey, species, acceptedScientificName, verbatimScientificName, match.name, .keep_all=TRUE)
On.plant.GBIF.sp.5<-data.frame(On.plant.GBIF.sp.5)
saveRDS(On.plant.GBIF.sp.5, "intermediate_rds_csv/GBIF/On.plant.GBIF.sp.5.rds")

#write translated column
#split into different file
hab.Opgs5<-unique(On.plant.GBIF.sp.5$habitat)
hab.translate<-google_translate(hab.Opgs5, target_language="en", source_language="auto")
hab.tr<-unlist(hab.translate)
hab.tr.df<-data.frame(hab.Opgs5, hab.tr)
saveRDS(hab.translate, "intermediate_rds_csv/GBIF/hab.translate.rds")
saveRDS(hab.tr.df, "intermediate_rds_csv/GBIF/hab.tr.df.rds")


#Here, we filtered out habitat information from hab.tr.df that did not infer any direct tissue contact at all.
hab.ed<-read.csv("intermediate_rds_csv/GBIF/hab.tr.df.FILT.csv")
#NOTE: The final form of hab.tr.df.FILT.csv is below. 

hab.ed$hab.Opgs5<-as.character(hab.ed$hab.Opgs5)
hab.ed$hab.tr<-as.character(hab.ed$hab.tr)

#a few R shortcuts found during habitat translation manual cleaning:
hab.ed1<-hab.ed[grep("med bakkekontakt", hab.ed$hab.Opgs5),]
hab.1<-hab.ed[which(!hab.ed$X %in% hab.ed1$X),]
hab.ed2<-hab.1[grep("på bakken", hab.1$hab.Opgs5, ignore.case=TRUE),]
hab.2<-hab.1[which(!hab.1$X %in% hab.ed2$X),]
hab.ed3<-hab.2[grep("Growing on bajo", hab.2$hab.Opgs5, ignore.case=TRUE),]
hab.3<-hab.2[which(!hab.2$X %in% hab.ed3$X),]
hab.ed4<-hab.3[grep(" on rock\\.", hab.3$hab.Opgs5, ignore.case=TRUE),]
hab.4<-hab.3[which(!hab.3$X %in% hab.ed4$X),]
hab.ed5<-hab.4[grep("on granit", hab.4$hab.Opgs5, ignore.case=TRUE),]
hab.5<-hab.4[which(!hab.4$X %in% hab.ed5$X),]
hab.ed6<-hab.5[grep("låg ", hab.5$hab.Opgs5, ignore.case=TRUE),]
hab.6<-hab.5[which(!hab.5$X %in% hab.ed6$X),]
hab.ed7<-hab.6[grep("dead ", hab.6$hab.tr, ignore.case=TRUE),]
hab.7<-hab.6[which(!hab.6$X %in% hab.ed7$X),]
hab.ed8<-hab.7[grep("liggande", hab.7$hab.Opgs5, ignore.case=TRUE),]
hab.8<-hab.7[which(!hab.7$X %in% hab.ed8$X),]
#Note: THESE WERE identified and removed from "intermediate_rds_csv/GBIF/hab.tr.df.FILT.csv"

#After some filtering, I found some patterns I missed. Here I will read in the Filtered Habitat file midway.
#more habitat filtration
match.5<-unique(On.plant.GBIF.sp.5$GBIF.plant)
wi<-On.plant.GBIF.sp.5[grep("with ", On.plant.GBIF.sp.5$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in match.5){
  mm<-ifelse(nn<-grep(paste("with ", i, ",", sep=""), wi$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
s.wi<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
ddd<-cbind(wi[s.wi$GBIF.entry,], s.wi)
dddd<-ddd[which(!grepl(" on ", ddd$habitat, ignore.case=TRUE)),]
ddd1<-dddd%>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))
#remove them from On.plant.GBIF.sp
On.plant.GBIF.sp.6<-setdiff(On.plant.GBIF.sp.5,ddd1[,c(1:25)])

ba<-On.plant.GBIF.sp.6[grep("base ", On.plant.GBIF.sp.6$habitat, ignore.case=TRUE),]
sp.in.li.GBIF.plant<-NULL
GBIF.entry<-NULL
for (i in match.5){
  mm<-ifelse(nn<-grep(paste("on base of", i), ba$habitat, ignore.case = TRUE),print(i),0)
  sp.in.li.GBIF.plant<-c(sp.in.li.GBIF.plant,mm)
  GBIF.entry<-c(GBIF.entry,nn)
}
s.ba<-data.frame(GBIF.entry,sp.in.li.GBIF.plant)
eee<-cbind(ba[s.ba$GBIF.entry,], s.ba)
eee1<-eee%>%
  filter(as.character(GBIF.plant) == as.character(sp.in.li.GBIF.plant))
#remove them from On.plant.GBIF.sp
On.plant.GBIF.sp.7<-setdiff(On.plant.GBIF.sp.6,eee1[,c(1:25)])

#from "Fjordii", habitat sorted
hab.tr.df.FILT<-read.csv("~/Desktop/hab.tr.df.FILT.csv")
wi.filt<-On.plant.GBIF.sp.6[which(On.plant.GBIF.sp.7$habitat %in% hab.tr.df.FILT$hab.Opgs5),]
hab.tr.df.FILT2<-hab.tr.df.FILT[which(hab.tr.df.FILT$hab.Opgs5 %in% wi.filt$habitat),]
hab.ed9<-hab.tr.df.FILT2[grep(" låg.", hab.tr.df.FILT2$hab.Opgs5, ignore.case=TRUE),]
hab.9<-hab.tr.df.FILT2[which(!hab.tr.df.FILT2$X %in% hab.ed9$X),]
hab.ed10<-hab.9[grep(":låg", hab.9$hab.Opgs5, ignore.case=TRUE),]
hab.10<-hab.9[which(!hab.9$X %in% hab.ed10$X),]
hab.ed11<-hab.10[grep("on moss", hab.10$hab.Opgs5, ignore.case=TRUE),]
hab.11<-hab.10[which(!hab.10$X %in% hab.ed11$X),]
hab.ed12<-hab.11[grep("lying trunk", hab.11$hab.Opgs5, ignore.case=TRUE),]
hab.12<-hab.11[which(!hab.11$X %in% hab.ed12$X),]
hab.ed13<-hab.12[grep("terricolous", hab.12$hab.Opgs5, ignore.case=TRUE),]
hab.13<-hab.12[which(!hab.12$X %in% hab.ed13$X),]
hab.ed14<-hab.13[grep("growing on rocks", hab.13$hab.Opgs5, ignore.case=TRUE),]
hab.14<-hab.13[which(!hab.13$X %in% hab.ed14$X),]
hab.ed15<-hab.14[grep("on downed", hab.14$hab.Opgs5, ignore.case=TRUE),]
hab.15<-hab.14[which(!hab.14$X %in% hab.ed15$X),]
hab.ed16<-hab.15[grep(" mort\\.", hab.15$hab.Opgs5, ignore.case=TRUE),]
hab.16<-hab.15[which(!hab.15$X %in% hab.ed16$X),]
hab.ed17<-hab.16[grep("emort", hab.16$hab.Opgs5, ignore.case=TRUE),]
hab.17<-hab.16[which(!hab.16$X %in% hab.ed17$X),]
hab.ed18<-hab.17[grep("död", hab.17$hab.Opgs5, ignore.case=TRUE),]
hab.18<-hab.17[which(!hab.17$X %in% hab.ed18$X),]
hab.ed19<-hab.18[grep(" sur sol ", hab.18$hab.Opgs5, ignore.case=TRUE),]
hab.19<-hab.18[which(!hab.18$X %in% hab.ed19$X),]
hab.ed20<-hab.19[grep("base of old ", hab.19$hab.Opgs5, ignore.case=TRUE),]
hab.20<-hab.19[which(!hab.19$X %in% hab.ed20$X),]
hab.ed21<-hab.20[grep("on limestone outcrops", hab.20$hab.Opgs5, ignore.case=TRUE),]
hab.21<-hab.20[which(!hab.20$X %in% hab.ed21$X),]
hab.ed22<-hab.21[grep("on sandstone", hab.21$hab.Opgs5, ignore.case=TRUE),]
hab.22<-hab.21[which(!hab.21$X %in% hab.ed22$X),]
hab.ed23<-hab.22[grep("on stone", hab.22$hab.Opgs5, ignore.case=TRUE),]
hab.23<-hab.22[which(!hab.22$X %in% hab.ed23$X),]
hab.ed24<-hab.23[grep(", base of trunk", hab.23$hab.Opgs5, ignore.case=TRUE),]
hab.24<-hab.23[which(!hab.23$X %in% hab.ed24$X),]
hab.ed25<-hab.24[grep("; on rock", hab.24$hab.Opgs5, ignore.case=TRUE),]
hab.25<-hab.24[which(!hab.24$X %in% hab.ed25$X),]
hab.ed26<-hab.25[grep("\\. on rock", hab.25$hab.Opgs5, ignore.case=TRUE),]
hab.26<-hab.25[which(!hab.25$X %in% hab.ed26$X),]
hab.ed27<-hab.26[grep("on root", hab.26$hab.Opgs5, ignore.case=TRUE),]
hab.27<-hab.26[which(!hab.26$X %in% hab.ed27$X),]


#This was saved as hab.tr.df.FILT2.csv, and continued to be edited manually as follows:

# Removed : base of trunk UNLESS corticolous, or otherwise indicates it is growing directly on tree.
# Note: Spanish entries (Habitat field starts with "Growing on" seem to be autopopulated... The Growing on should be ignored.)
# Revisit these entries, perhaps discuss need for careful use of language.
## As in Mycoportal, except I coded in removal of "sticks", "burned", "burnt", "charred". I manually deleted entries in all languages at base of tree or trunk AND those that were on rotten/ing plants.
#note: entries kept if they were explicitly stated to be the reason for rot of host plant.
#note: "effective" duplicates (more than one plant host erroneously detected per entry) retained for now.
#note: rhizosphere entries removed.
#note: if host species has ? (to signify uncertainty in identification), it was deleted
#note: decorticated kept. wood kept. host: "" kept, unless clear it was dead tissue. dried tissue of plant kept. 
#note: if host was other fungi on plant, the plant host record was retained, and fungal host record dropped (as in Aecidi(a)(um))
#note: eliminated all badly abbreviated entries (i.e., U. californica for Fragaria californica)
#note: deleted entries with abbreviated genus names "#. canadensis" if no full genus names provided in separate fields.
#note: records retained if associated.taxa had name of host plant but habitat did not, but was still suitable (i.e., "leaves of plant", but not "dead leaves of plant")
#note: If abbreviated genus provided with other full genus (Cornus alternifolia and C. sericea), I inputted full genus name.

hfilt2<-read.csv("~/data2/LPs/rds/hab.tr.df.FILT2.csv")
hfilt.GBIF<-On.plant.GBIF.sp.7[which(On.plant.GBIF.sp.7$habitat %in% hfilt2$hab.Opgs5),]
hfilt.GBIF$trnum<-match(hfilt.GBIF$habitat, hfilt2$hab.Opgs5, nomatch = NA_integer_, incomparables = NULL)
#append habitat translation.
hfilt.GBIF$hab.trans<-print(hfilt2$hab.Opgs5[hfilt.GBIF$trnum])
#saveRDS(hfilt.GBIF, "intermediate_rds_csv/GBIF/hfilt.GBIF.rds"). This file is unaltered.

######~~~~~~######
# Now we have a rough draft of a dataframe. It is imperfect due to a few plants in each habitat field (e.g, one is host, other is in area) 
# and human error from skimming habitat information quickly.
# Now the plan is to draw unique plant-fungal pairings from this.
# This is a good resource for living and dead plant-associating fungi, and fungi on the trunks of trees. Includes lichens.
#####~~~~~#######

#there are a number of entries from Spain that seem to have "Growing on" as a generic field beginner (e.g, "Growing on hojas de Abies"). Remove this text and create new column for original text for these entries.
nonspanish<-hfilt.GBIF %>%
  filter(!countryCode %in% "ES" & !grepl("Growing on", hfilt.GBIF$habitat))
nonspanish$spanish.original.habitat<-NA
spanish<-hfilt.GBIF %>%
  filter(countryCode %in% "ES" & grepl("Growing on", hfilt.GBIF$habitat))
spanish$spanish.original.habitat<-spanish$habitat
spanish$habitat<-gsub("Growing on ", "",spanish$spanish.original.habitat)
hfilt.GBIF.ns<-rbind.data.frame(spanish, nonspanish)

#remove lichens
#read in lichenlist data from "Consortium of Lichen Herbaria" https://lichenportal.org/portal/checklists/checklist.php?clid=1492&pid=558
lichenlist<-read.csv("ext_data/World Checklist of Genera of Lichenized Fungi_1698006802.csv")
#make genus column
lichenlist$genus<-gsub(" sp.", "", lichenlist$ScientificName)
#filter lichenizing fungi
hfilt.GBIF.nl<-hfilt.GBIF.ns[which(!hfilt.GBIF.ns$genus %in% lichenlist$genus),]
#remove fungal taxon entries that contain word root "lichen"
nlaf<-hfilt.GBIF.nl[grep("lichen",hfilt.GBIF.nl$scientificName, ignore.case=TRUE),]
hfilt.GBIF.nl2<-hfilt.GBIF.nl[which(!hfilt.GBIF.nl$scientificName %in% nlaf$scientificName),]

saveRDS(hfilt.GBIF.nl2, "intermediate_rds_csv/GBIF/hfilt.GBIF.nl2.rds")
#### Start by identifying the entries with the clearest evidence of living host association.
#remove dry tissue or trunk entries
cl.hf.G1<-setdiff(hfilt.GBIF.nl2, hfilt.GBIF.nl2[grep(" on trunk of ", hfilt.GBIF.nl2$habitat, ignore.case=TRUE),])
cl.hf.G2<-setdiff(cl.hf.G1, cl.hf.G1[grep("dry leaves of ", cl.hf.G1$habitat, ignore.case=TRUE),])
cl.hf.G3<-setdiff(cl.hf.G2, cl.hf.G2[grep("dry twigs of ", cl.hf.G2$habitat, ignore.case=TRUE),])
cl.hf.G4<-setdiff(cl.hf.G3, cl.hf.G3[grep("dead", cl.hf.G3$habitat, ignore.case=TRUE),])
cl.hf.G5<-setdiff(cl.hf.G4, cl.hf.G4[grep("stump", cl.hf.G4$habitat, ignore.case=TRUE),])
cl.hf.G6<-setdiff(cl.hf.G5, cl.hf.G5[grep("trunk", cl.hf.G5$habitat, ignore.case=TRUE),])
cl.hf.G7<-setdiff(cl.hf.G6, cl.hf.G6[grep("bark", cl.hf.G6$habitat, ignore.case=TRUE),])
cl.hf.G8<-setdiff(cl.hf.G7, cl.hf.G7[grep("leaf litter", cl.hf.G7$habitat, ignore.case=TRUE),])
cl.hf.G9<-setdiff(cl.hf.G8, cl.hf.G8[grep("unknown", cl.hf.G8$associatedTaxa, ignore.case=TRUE),])
cl.hf.G10<-setdiff(cl.hf.G9, cl.hf.G9[grep(" arid", cl.hf.G9$habitat, ignore.case=TRUE),])
cl.hf.G11<-setdiff(cl.hf.G10, cl.hf.G10[grep(" langui ", cl.hf.G10$habitat, ignore.case=TRUE),])
cl.hf.G12<-setdiff(cl.hf.G11, cl.hf.G11[grep("sicc", cl.hf.G11$habitat, ignore.case=TRUE),])
cl.hf.G13<-setdiff(cl.hf.G12, cl.hf.G12[grep("mort", cl.hf.G12$habitat, ignore.case=TRUE),])
cl.hf.G14<-setdiff(cl.hf.G13, cl.hf.G13[grep("dürren", cl.hf.G13$habitat, ignore.case=TRUE),])
cl.hf.G15<-setdiff(cl.hf.G14, cl.hf.G14[grep("secos", cl.hf.G14$habitat, ignore.case=TRUE),])
cl.hf.G16<-setdiff(cl.hf.G15, cl.hf.G15[grep("last year", cl.hf.G15$habitat, ignore.case=TRUE),])
cl.hf.G17<-setdiff(cl.hf.G16, cl.hf.G16[grep(" dry ", cl.hf.G16$habitat, ignore.case=TRUE),])
cl.hf.G18<-setdiff(cl.hf.G17, cl.hf.G17[grep("withering", cl.hf.G17$habitat, ignore.case=TRUE),])
cl.hf.G19<-setdiff(cl.hf.G18, cl.hf.G18[grep("corteza", cl.hf.G18$habitat, ignore.case=TRUE),])
cl.hf.G20<-setdiff(cl.hf.G19, cl.hf.G19[grep("overwintered", cl.hf.G19$habitat, ignore.case=TRUE),])
HGTdf<-setdiff(cl.hf.G20, cl.hf.G20[grep("ramas de", cl.hf.G20$habitat, ignore.case=TRUE),])




### good data: "on leaves of [species]"
HGTdf.1<-data.frame(HGTdf[grep("leaves of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.2<-data.frame(HGTdf[grep("n the leaves of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.3<-data.frame(HGTdf[grep("n living leaves of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.4<-data.frame(HGTdf[grep("n flowers of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.5<-data.frame(HGTdf[grep("n seeds of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.6<-data.frame(HGTdf[grep("n fruits of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.7<-data.frame(HGTdf[grep("n inflorescence of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.8<-data.frame(HGTdf[grep("n culms of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.9<-data.frame(HGTdf[grep("n living twigs of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.10<-data.frame(HGTdf[grep("n living stems of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.11<-data.frame(HGTdf[grep("n ovaries of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.12<-data.frame(HGTdf[grep("ad folia", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.13<-data.frame(HGTdf[grep("An Blättern von", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.14<-data.frame(HGTdf[grep("hojas de", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.15<-data.frame(HGTdf[grep("n leaves and stems of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.16<-data.frame(HGTdf[grep("n leaves and petioles of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.17<-data.frame(HGTdf[grep("n perigyni", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.18<-data.frame(HGTdf[which(HGTdf$habitat %in% (paste("Matrix: In", HGTdf$GBIF.plant))),])
HGTdf.19<-data.frame(HGTdf[which(HGTdf$habitat %in% (paste("Sobre/", HGTdf$GBIF.plant, ", Hoja", sep =""))),])
HGTdf.20<-data.frame(HGTdf[which(HGTdf$habitat %in% (paste("Sobre/", HGTdf$GBIF.plant, ", Planta viva", sep =""))),])
HGTdf.21<-data.frame(HGTdf[grep("den Blättern von", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.22<-data.frame(HGTdf[grep("utrículos de", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.23<-data.frame(HGTdf[grep("Sur les feuilles", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.24<-data.frame(HGTdf[grep("n anthers", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.25<-data.frame(HGTdf[grep("n the anthers", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.26<-data.frame(HGTdf[grep("rust", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.27<-data.frame(HGTdf[grep("uredospore", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.28<-data.frame(HGTdf[grep(" witches", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.29<-data.frame(HGTdf[grep("På blad av", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.30<-data.frame(HGTdf[grep("ovarios de", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.31<-data.frame(HGTdf[grep("levende stam van", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.32<-data.frame(HGTdf[grep("Op blad", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.33<-data.frame(HGTdf[grep("n the lower side of leaves of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.34<-data.frame(HGTdf[grep("n stems of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.35<-data.frame(HGTdf[grep("n stems and leaves of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.36<-data.frame(HGTdf[grep("n stem of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.37<-data.frame(HGTdf[grep("n stalk of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.38<-data.frame(HGTdf[grep("n stalks of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.39<-data.frame(HGTdf[grep("n spikelet of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.40<-data.frame(HGTdf[grep("n seedlings of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.41<-data.frame(HGTdf[grep("n petioles of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.42<-data.frame(HGTdf[grep("n male catkins", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.43<-data.frame(HGTdf[grep("n male spikes", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.44<-data.frame(HGTdf[grep("n live leaves", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.45<-data.frame(HGTdf[grep("n living and wilt", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.46<-data.frame(HGTdf[grep("n leaves, sheaths", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.47<-data.frame(HGTdf[grep("n leaves, stalks", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.48<-data.frame(HGTdf[grep("n leaves, stems", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.49<-data.frame(HGTdf[grep("n leaves, petioles", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.50<-data.frame(HGTdf[grep("n leaves & stems of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.51<-data.frame(HGTdf[grep("n leaves and culms of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.52<-data.frame(HGTdf[grep("n leaves and bracts of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.53<-data.frame(HGTdf[grep("n fruit of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.54<-data.frame(HGTdf[grep("n fronds of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.55<-data.frame(HGTdf[grep("n female catkins of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.56<-data.frame(HGTdf[grep("n culms and leaves of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.57<-data.frame(HGTdf[grep("n canes of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.58<-data.frame(HGTdf[grep("on attached needles", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.59<-data.frame(HGTdf[grep("on attached leaves", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.60<-data.frame(HGTdf[grep("Substrat:Blad", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.61<-data.frame(HGTdf[grep("folia viva", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.62<-data.frame(HGTdf[grep("in foliis", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.63<-data.frame(HGTdf[grep("in utricles of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.64<-data.frame(HGTdf[grep("in floribus", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.65<-data.frame(HGTdf[grep("in antheris", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.66<-data.frame(HGTdf[grep("^leaves of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.67<-data.frame(HGTdf[grep("leaf spot", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.68<-data.frame(HGTdf[grep("leaf lesion", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.69<-data.frame(HGTdf[grep("leaf and stem gall", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.70<-data.frame(HGTdf[grep("inflorescencia", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.71<-data.frame(HGTdf[grep("in spikelets", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.72<-data.frame(HGTdf[grep("in petiolis", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.73<-data.frame(HGTdf[grep("in ovariis", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.74<-data.frame(HGTdf[grep("in internodes", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.75<-data.frame(HGTdf[grep("in internodiis", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.76<-data.frame(HGTdf[grep("in fol.", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.77<-data.frame(HGTdf[grep("In den Bluten von", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.78<-data.frame(HGTdf[grep("In den Antheren von", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.79<-data.frame(HGTdf[grep("flores de ", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.80<-data.frame(HGTdf[grep("espor hojas d", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.81<-data.frame(HGTdf[grep("espigas de", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.82<-data.frame(HGTdf[grep("n hoja de", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.83<-data.frame(HGTdf[grep("n frutas de", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.84<-data.frame(HGTdf[grep("n leaves and bracts of", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.85<-data.frame(HGTdf[grep("Auf stengeln von", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.86<-data.frame(HGTdf[grep("lebenden Blättern", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.87<-data.frame(HGTdf[grep("Auf Asten von", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.88<-data.frame(HGTdf[grep("anteras de", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.89<-data.frame(HGTdf[grep("Aecidien auf", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.90<-data.frame(HGTdf[grep("ad caules", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.91<-data.frame(HGTdf[grep("Auf Blättern", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.92<-data.frame(HGTdf[grep("På blad av", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.93<-data.frame(HGTdf[grep("en tallo de", HGTdf$habitat, ignore.case=TRUE),])
HGTdf.94<-data.frame(HGTdf[grep("en tallos de", HGTdf$habitat, ignore.case=TRUE),])

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~# #1. Kept dataset
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
HGTdf.keep<-rbind.data.frame(HGTdf.1, HGTdf.2, HGTdf.3, HGTdf.4, HGTdf.5, HGTdf.6, HGTdf.7, HGTdf.8, HGTdf.9, 
                                 HGTdf.10, HGTdf.11, HGTdf.12, HGTdf.13, HGTdf.14, HGTdf.15, HGTdf.16, HGTdf.17, HGTdf.18, HGTdf.19, 
                                 HGTdf.20, HGTdf.21, HGTdf.22, HGTdf.23, HGTdf.24, HGTdf.25, HGTdf.26, HGTdf.27, HGTdf.28, HGTdf.29,
                                 HGTdf.30, HGTdf.31, HGTdf.32, HGTdf.33, HGTdf.34, HGTdf.35, HGTdf.36, HGTdf.37, HGTdf.38, HGTdf.39,
                                 HGTdf.40, HGTdf.41, HGTdf.42, HGTdf.43, HGTdf.44, HGTdf.45, HGTdf.46, HGTdf.47, HGTdf.48, HGTdf.49,
                                 HGTdf.50, HGTdf.51, HGTdf.52, HGTdf.53, HGTdf.54, HGTdf.55, HGTdf.56, HGTdf.57, HGTdf.58, HGTdf.59,
                                 HGTdf.60, HGTdf.61, HGTdf.62, HGTdf.63, HGTdf.64, HGTdf.65, HGTdf.66, HGTdf.67, HGTdf.68, HGTdf.69,
                                 HGTdf.70, HGTdf.71, HGTdf.72, HGTdf.73, HGTdf.74, HGTdf.75, HGTdf.76, HGTdf.77, HGTdf.78, HGTdf.79,
                                 HGTdf.80, HGTdf.81, HGTdf.82, HGTdf.83, HGTdf.84, HGTdf.85, HGTdf.86, HGTdf.87, HGTdf.88, HGTdf.89, 
                             HGTdf.90, HGTdf.91, HGTdf.92, HGTdf.93, HGTdf.94)


# which occurrences were not detected between living plant and fungi (sing.HGTdf.all)?
undet.HGTdf<-setdiff(HGTdf, HGTdf.keep)
# which of these were from combos already detected as occurring on live plant tissue?
#1:these combinations have already been seen on living plant tissue, but not necessarily in this occurrence. Note that I am alowing for leaf litter fungi and the
#rest in the hfilt.GBID.nl2 to be included here.
kept.combos<-undet.HGTdf[paste(undet.HGTdf$GBIF.plant, undet.HGTdf$species) %in% 
                           c(paste(hfilt.GBIF.nl2$GBIF.plant, hfilt.GBIF.nl2$species),
                             paste(hfilt.GBIF.nl2$species, hfilt.GBIF.nl2$GBIF.plant)), ]
#2. these combinations have not yet been seen. must be processed still. I will not include leaf litter (and the rest of the hfilt.GBID.nl2 set here) to be further assessed.
unkept.combos<-undet.HGTdf[! paste(undet.HGTdf$GBIF.plant, undet.HGTdf$species) %in% 
      c(paste(HGTdf.keep$GBIF.plant, HGTdf.keep$species),
        paste(HGTdf.keep$species, HGTdf.keep$GBIF.plant)), ]

###########
###Now I will subset the data to just the data with habitat fields that say "On [plant name]"
# I will delete tree species from this data, as these may be corticolous.
on.unkept.match<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("On", unkept.combos$match.name))),])
on.unkept.match2<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("on", unkept.combos$match.name))),])
on.unkept.match3<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("On", unkept.combos$GBIF.plant))),])
on.unkept.match4<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("on", unkept.combos$GBIF.plant))),])
auf.unkept.match<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("Auf", unkept.combos$match.name))),])
auf.unkept.match2<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("auf", unkept.combos$match.name))),])
auf.unkept.match3<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("Auf", unkept.combos$GBIF.plant))),])
auf.unkept.match4<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("Auf ", unkept.combos$GBIF.plant, ".", sep = ""))),])
auf.unkept.match5<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("auf", unkept.combos$GBIF.plant))),])
en.unkept.match<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("En", unkept.combos$match.name))),])
en.unkept.match2<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("en", unkept.combos$match.name))),])
en.unkept.match3<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("En", unkept.combos$GBIF.plant))),])
en.unkept.match4<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("en", unkept.combos$GBIF.plant))),])
sobre.unkept.match<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("Sobre", unkept.combos$match.name))),])
sobre.unkept.match2<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("sobre", unkept.combos$match.name))),])
sobre.unkept.match3<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("Sobre", unkept.combos$GBIF.plant))),])
sobre.unkept.match4<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("sobre", unkept.combos$GBIF.plant))),])
ad.unkept.match<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("Ad", unkept.combos$match.name))),])
ad.unkept.match2<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("ad", unkept.combos$match.name))),])
ad.unkept.match3<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("Ad", unkept.combos$GBIF.plant))),])
ad.unkept.match4<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("ad", unkept.combos$GBIF.plant))),])
an.unkept.match<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("An", unkept.combos$match.name))),])
an.unkept.match2<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("an", unkept.combos$match.name))),])
an.unkept.match3<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("An", unkept.combos$GBIF.plant))),])
an.unkept.match4<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("an", unkept.combos$GBIF.plant))),])
op.unkept.match<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("Op", unkept.combos$GBIF.plant))),])
op.unkept.match2<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("Op ", unkept.combos$GBIF.plant, ".", sep = ""))),])
pa.unkept.match<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("På", unkept.combos$GBIF.plant))),])
pa.unkept.match2<-data.frame(unkept.combos[which(unkept.combos$habitat %in% (paste("På ", unkept.combos$GBIF.plant, ".", sep = ""))),])

unkept.on<-rbind.data.frame(on.unkept.match, on.unkept.match2, on.unkept.match3, on.unkept.match4, auf.unkept.match, auf.unkept.match2,
                             auf.unkept.match3,auf.unkept.match4,auf.unkept.match5,en.unkept.match,en.unkept.match2,en.unkept.match3,en.unkept.match4,
                             sobre.unkept.match,sobre.unkept.match2,sobre.unkept.match3,sobre.unkept.match4,ad.unkept.match,ad.unkept.match2,ad.unkept.match3,
                             ad.unkept.match4,an.unkept.match,an.unkept.match2,an.unkept.match3,an.unkept.match4, op.unkept.match, op.unkept.match2, 
                            pa.unkept.match, pa.unkept.match2)

#I retrieved the csv for the GlobalTreeSearch data on Oct 31 2023. https://rpubs.com/Roeland-KINDT/812716
#read this data in 
globtree<-read.csv("ext_data/global_tree_search_trees_1_7.csv")
#filter by genus in order to eliminate as many shrubby genera, as well.
gltrgenus<-data.frame(str_split_fixed(globtree$TaxonName, " ", 2))
globtree$genus<-gltrgenus$X1
unk.on.genus<-data.frame(str_split_fixed(unkept.on$match.name, " ", 2))
unkept.on$plant.genus<-unk.on.genus$X1
#saveRDS(unkept.on, "intermediate_rds_csv/GBIF/rds/unkept.on.rds")

#filter out tree data from the "on [species]" list

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
######## #2. Kept dataset ####
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
unkept.on.no.tree<-unkept.on[which(!unkept.on$plant.genus %in% globtree$genus),]


#what is leftover? Non ("on-[plant-species]" entries), and excluding "on-[plant-species]" tree species.
not.on.un.com<-setdiff(unkept.combos, unkept.on[,c(1:28)])
#create non-tree list
nouc.genus<-data.frame(str_split_fixed(not.on.un.com$match.name, " ", 2))
not.on.un.com$plant.genus<-nouc.genus$X1

#filter out tree data from the "on [species]" list
nouc.no.tree<-not.on.un.com[which(!not.on.un.com$plant.genus %in% globtree$genus),]

######
#prepare host data of non-trees in associatedTaxa field
nouc.assocta<-nouc.no.tree[which(!is.na(nouc.no.tree$associatedTaxa)),]
nouc.assocta2<-setdiff(nouc.assocta, nouc.assocta[grep("dead ", nouc.assocta$associatedTaxa, ignore.case=TRUE),])
nouc.assocta3<-setdiff(nouc.assocta2, nouc.assocta2[grep("root", nouc.assocta2$associatedTaxa, ignore.case=TRUE),])
### ^ here are more good entries to add (associated Taxa of non-tree species, not on roots or dead tissue.)
#### make df of distinct fungi/plant combos THAT WE WANT TO KEEP ("singles" list) -- (AssociatedTaxa entries - non tree)

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
######## #3. Kept dataset ####
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
nas1<-nouc.assocta3 %>%
  group_by(species, GBIF.plant) %>%
  arrange(str_count(habitat), .by_group=TRUE) %>%
  distinct(species, .keep_all=TRUE)

nouc.noassoc<-setdiff(nouc.no.tree, nouc.assocta)
#edit manually. Sort by habitat and remove entries in Numbers app. delete "old, dried, last year's, straw, cortex/bark, withering". Deleted all entries from Spain that 
# were just "Growing on [plant species]", but kept e.g. "Growing on en el [plant sp.]". "In mountain meadow; on Verbascum thapsus"
write.csv(nouc.noassoc, "intermediate_rds_csv/GBIF/nouc.noassoc.csv")

fil.nouc.nonassoc<-read.csv("intermediate_rds_csv/GBIF/nouc.noassoc.filt.csv")

#### make df of distinct fungi/plant combos THAT WE WANT TO KEEP ("singles" list) -- (final, hand filtered data)

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
######## #4. Kept dataset ####
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
fn.na1<-fil.nouc.nonassoc %>%
  group_by(species, GBIF.plant) %>%
  arrange(str_count(habitat), .by_group=TRUE) %>%
  distinct(species, .keep_all=TRUE)


HGTdf.keep$plant.genus<-word(HGTdf.keep$match.name)
GBIF.wd.dups.raw<-rbind.data.frame(HGTdf.keep, nas1, unkept.on.no.tree, fn.na1[,2:30])
GBIF.wd.dups.raw$match.name.sp<-word(GBIF.wd.dups.raw$match.name, 1, 2)
saveRDS(GBIF.wd.dups.raw, "intermediate_rds_csv/GBIF/GBIF.wd.dups.raw.rds")
match.name.Gwdr<-unique(GBIF.wd.dups.raw$match.name.sp)
#use U.Taxonstand
corr.plantGwdr<-nameMatch_WCVP(match.name.Gwdr)
saveRDS(corr.plantGwdr, "intermediate_rds_csv/GBIF/corr.plantGwdr.rds")


GB1.Gwdr<-GBIF.wd.dups.raw[which(GBIF.wd.dups.raw$match.name.sp %in% corr.plantGwdr.f$Submitted_Name),]

GB1.Gwdr$genus_in_db<-corr.plantGwdr$Genus_in_database[match(GB1.Gwdr$match.name.sp, corr.plantGwdr$Submitted_Name)]
GB1.Gwdr$name_in_db<-corr.plantGwdr$Name_in_database[match(GB1.Gwdr$match.name.sp, corr.plantGwdr$Submitted_Name)]
GB1.Gwdr$author_in_db<-corr.plantGwdr$Author_in_database[match(GB1.Gwdr$match.name.sp, corr.plantGwdr$Submitted_Name)]
GB1.Gwdr$new_name_wcvp<-corr.plantGwdr$New_name[match(GB1.Gwdr$match.name.sp, corr.plantGwdr$Submitted_Name)]
GB1.Gwdr$new_author_wcvp<-corr.plantGwdr$New_author[match(GB1.Gwdr$match.name.sp, corr.plantGwdr$Submitted_Name)]
corrected.Host.sp<-GB1.Gwdr$new_name_wcvp
corrected.Host.sp[which(is.na(corrected.Host.sp))]<-GB1.Gwdr$name_in_db[which(is.na(corrected.Host.sp))]
GB1.Gwdr$correct.genusplant<-word(corrected.Host.sp, 1)

#isolate trees, as the genera may now be detected after taxonomical re-assessment (but weren't before)
GB1.tree<-GB1.Gwdr[which(GB1.Gwdr$correct.genusplant %in% globtree$genus),]
# write.csv(GB1.tree, "intermediate_rds_csv/GBIF/2024-12-21.GB1.tree.csv")
#remove entries that say host
#no tree - these are done.
GB1.notree<-GB1.Gwdr[-which(GB1.Gwdr$correct.genusplant %in% globtree$genus),]

#Edited tree data, as there were still "on stem", "on stems", "host:", "on <plant name>", "Sobre/...<planta viva", and "Host of:" entries here.
GB1.tree.filt<-read.csv("intermediate_rds_csv/GBIF/2024-12-21.GB1.tree.FILT.csv")
GB2<-rbind.data.frame(GB1.tree.filt[,2:37], GB1.notree)
saveRDS(GB2, "intermediate_rds_csv/GBIF/GB2.rds")

GB3<-semi_join(GB1.Gwdr, GB2)
saveRDS(GB3, "intermediate_rds_csv/GBIF/GB3.rds")

#modify data (from the corr.plant taxon sstandardization) so that it is consistent with the other sources
GB3$Accepted_SPNAME<-corr.plantGwdr$Accepted_SPNAME[match(GB3$match.name.sp, corr.plantGwdr$Submitted_Name)]
GB3$Accepted_SPAUTHOR<-GB3$author_in_db
GB3$Accepted_SPAUTHOR[which(!is.na(GB3$new_author_wcvp))]<-GB3$new_author_wcvp[which(!is.na(GB3$new_author_wcvp))]

#remove entries without accepted sp name
GB3.f<-GB3[which(!is.na(GB3$Accepted_SPNAME)),]


#only keep unique fungi-plant pairs
GB4<-GB3.f %>%
  group_by(Accepted_SPNAME, species) %>%
  distinct(Accepted_SPNAME, .keep_all=TRUE)

#make final dataframe with both corrected and uncorrected host names in separate columns.

#
corrected.Host.sp<-GB4$Accepted_SPNAME
corrected.Host.ssp<-word(GB4$Accepted_SPNAME, 3, -1, sep=" ")
corrected.Host.auth<-GB4$Accepted_SPAUTHOR
detected.Host.plant<-as.character(GB4$GBIF.plant)
accepted.Fungi.sp<-GB4$species
accepted.Fungi.sp.full.name<-GB4$acceptedScientificName
verbatim.Fungi.sp<-GB4$verbatimScientificName
habitat.and.plant.part.1<-GB4$habitat
habitat.and.plant.part.2<-GB4$associatedTaxa
habitat.and.plant.part.3<-GB4$associatedOrganisms
DB<-as.character(GB4$db)
Location<-countrycode(GB4$countryCode, origin='iso2c', destination='country.name')
Ref<-GB4$occurrenceID
accepted.Fungi.auth<-word(GB4$acceptedScientificName, 3, -1)

#make final dataframe with both corrected and uncorrected host names in separate columns.
GBIF.df<-cbind.data.frame(corrected.Host.sp, corrected.Host.ssp, corrected.Host.auth, 
                              detected.Host.plant, accepted.Fungi.sp, accepted.Fungi.sp.full.name, accepted.Fungi.auth, verbatim.Fungi.sp, 
                              habitat.and.plant.part.1, habitat.and.plant.part.2, habitat.and.plant.part.3, Location, Ref, DB)


# III. Here we create GBIF.df, which will be used as raw data for the "On living plant tissue"  database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

saveRDS(GBIF.df, "Cleaned_data/2025-01.GBIF.df.rds")

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#











