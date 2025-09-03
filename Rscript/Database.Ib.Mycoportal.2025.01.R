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

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ Mycoportal ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

### PART I ### DOWNLOAD

# Downloaded on May 10 2022.

# Downloaded all data I could find that was not in GBIF via command line.
# Acad.Sci.occurrences.csv Humboldt.fung.occurrences.csv  u.south.carolina.occurrences.csv
# Alab.chytrid.csv    Miami.occurences.csv  u.south.florida.occurrences.csv
# Ariz.myco.herb.spec.csv Santa.barb.occurrences.csv  utah.occurrences.csv
# Arkansas.fung.csv  u.cincinatti.occurrences.csv   u.tennessee.chattanooga.occurrences.csv
# berkeley.occurrences.csv   u.hawaii.occurrences.csv  u.wyoming.occurrences.csv
# CORD.occurrences.csv    u.michigan.occurrences.csv  u.newmexico.occurrences.csv

# Hand downloaded fungi by phylum/subdivision in "Brazil SpeciesLink Fungi data from iDigBio"
# "Ecuador Fungi data from FungiWebEcuador(QCAM)" "Foray Newfoundland and Labrador Fungarium"
# "Institute of the Industrial Ecology Problems of the North of Kola Science Center of the Russian Academy of Sciences."
# "Jewell and Arline Moss Settle Herbarium at SUNY Oneonta" "National Herbarium of Mexico Fungal Collection"
# "Royal Botanic Garden Edinburgh" "Royal Botanic Garden Edinburgh" "University of Maine, Richard Homola Mycological Herbarium"
# "Addis Ababa University, observation-based" "Cryptogamic Russian Information System (CRIS) Observational Data"
# "Foray Newfoundland and Labrador, observation-based" "Fungal Diversity Survey" "Fungal records database of Khanty-Mansi Autonomous Okrug – Yugra"
# "Index of the C.G. Lloyd Mycological Collection Specimens Housed at BPI" "Indian Ascomycetes/Marine/Mushroom/Myxomycetes/Rust Fungal Database(s)"
# "Long Island Mycological Club" "Malta Mycological Association" "Museo Nacional de Costa Rica, observation-based"
# "NAMP - New York Mycological Society: Macrofungi of New York City, New York" "Purdue University Northwest Mycology"
# "University of Arizona, Gilbertson Mycological Herbarium, observation-based" "General Observation and Personal Collections"

# I omitted the MushroomObserver data in Mycoportal because the association information is embedded in a long "Notes" field.
# Perhaps I could parse this, but I am skipping this for now.
# I had to open the u.michigan data in excel and resave it as Mac formatted csv instead of csv, because R could not interpret some of the comments within certain fields and screwed up the file. 

# in bash:

# cat ~/data2/LPs/Fungi.occurence.data/MyCoPortal.curled/*.csv | sort -u -r > Sorted.curled.Mycoportal.data.25May22.csv
# note: the 15may22 file has the old, badly formatted u.michigan file, but 25may has the revised file.
# cat ~/data2/LPs/Fungi.occurence.data/Mycoportal.handdownloaded/*.csv | sort -u -r > Sorted.handdownloaded.Mycoportal.data.15May22.csv

## Myco.po <- read.csv("~/data2/LPs/Fungi.occurence.data/Sorted.curled.Mycoportal.data.25May22.csv")
# Myco.po2<- read.csv("~/data2/LPs/Fungi.occurence.data/Sorted.handdownloaded.Mycoportal.data.15May22.csv")



##### PART II ##### Data wrangling and filtering

Myco.po1.ab <- data.frame(id=Myco.po$id, institutionCode = Myco.po$institutionCode, collectionCode = Myco.po$collectionCode, occurrenceID = Myco.po$occurrenceID, basisOfRecord = Myco.po$basisOfRecord,
                          scientificName = Myco.po$scientificName, genus = Myco.po$genus, taxonID = Myco.po$taxonID, taxonRank = Myco.po$taxonRank, occurrenceRemarks = Myco.po$occurrenceRemarks, habitat = Myco.po$habitat, 
                          associatedTaxa = Myco.po$associatedTaxa, substrate = rep("NA"), recordId = Myco.po$recordId, country = Myco.po$country, db=rep("Mycoportal.db"))
# the handdownloaded data do not have a substrate column, for whatever reason.
Myco.po2.ab <- data.frame(id=Myco.po2$id, institutionCode = Myco.po2$institutionCode, collectionCode = Myco.po2$collectionCode, occurrenceID = Myco.po2$occurrenceID, basisOfRecord = Myco.po2$basisOfRecord, 
                          scientificName = Myco.po2$scientificName, genus = Myco.po2$genus, taxonID = Myco.po2$taxonID, taxonRank = Myco.po2$taxonRank, occurrenceRemarks = Myco.po2$occurrenceRemarks, habitat = Myco.po2$habitat, 
                          associatedTaxa = Myco.po2$associatedTaxa, substrate = Myco.po2$substrate, recordId = Myco.po2$recordId, country = Myco.po2$country, db=rep("Mycoportal.db"))

# 1. everything, unfiltered(contains animal associates
Myco.po.ab<-rbind(Myco.po1.ab, Myco.po2.ab)
###not sure if this is an issue.... Warning messages:
##1: In `[<-.factor`(`*tmp*`, ri, value = c(9851470L, 9851469L, 9851446L,  : invalid factor level, NA generated. 2: In `[<-.factor`(`*tmp*`, ri, value = c(11736L, 187152L, 136862L,  :invalid factor level, NA generated

### saveRDS(Myco.po.ab, 'OSF/Mycoportal/Myco.po.ab.rds')
Myco.po.ab$associatedTaxa<-dplyr::na_if(Myco.po.ab$associatedTaxa, "")
Myco.po.ab$habitat<-dplyr::na_if(Myco.po.ab$habitat, "")
Myco.po.ab$substrate<-dplyr::na_if(Myco.po.ab$substrate, "")
Myco.po.ab$substrate<-dplyr::na_if(Myco.po.ab$substrate, "NA")

#2. everything with putative host fields
Myco.po.ab <- Myco.po.ab %>% filter_at(vars(habitat, substrate, associatedTaxa),any_vars(!is.na(.)))
facts <- sapply(Myco.po.ab, is.factor)                          
Myco.po.ab[facts] <- lapply(Myco.po.ab[facts], as.character) 

#3. filter to all association data (does it have a latin name in any text field). These imply some level of plant-fungal association. 
#note: there is nothing in associated.occurrences, this is later omitted.
associates<- data.frame(associatedTaxa=Myco.po.ab$associatedTaxa, associatedOccurrences=Myco.po.ab$associatedOccurrences, habitat=Myco.po.ab$habitat, substrate=Myco.po.ab$substrate)
associates<-unique(associates)
# write.csv(associates, 'intermediate_rds_csv/Mycoportal/associations.mycoportal.csv')


# BASH: 
# gnfinder -i intermediate_rds_csv/Mycoportal/associations.mycoportal.csv -u > intermediate_rds_csv/Mycoportal/associations.mycoportal.gnfinder.results.NOuninomial.csv
# gnf<-read.csv('intermediate_rds_csv/Mycoportal/associations.mycoportal.gnfinder.results.Nouninomial.csv')

# this leaves us with a set of names we can subset for. Unfortunately a dplyr solution seesm to exhaust R memory. as did {{ subset(habitat, subset = grepl(paste(name,collapse="|"), habitat)) }}
# write.csv(Myco.po.ab, 'OSF/Mycoportal/Myco.po.ab.filt.csv')

Mycop.plant<-NULL
Mycop.entry<-NULL
for (i in gnf$Name){
  mm<-ifelse(nn<-grep(i, Myco.po.ab$habitat, ignore.case = TRUE),print(i),0)
  Mycop.plant<-c(Mycop.plant,mm)
  Mycop.entry<-c(Mycop.entry,nn)
}
hits.habitat<-data.frame(Mycop.entry,Mycop.plant)

Mycop.plant<-NULL
Mycop.entry<-NULL
for (i in gnf$Name){
  mm<-ifelse(nn<-grep(i, Myco.po.ab$substrate, ignore.case = TRUE),print(i),0)
  Mycop.plant<-c(Mycop.plant,mm)
  Mycop.entry<-c(Mycop.entry,nn)
}
hits.substrate<-data.frame(Mycop.entry,Mycop.plant)

Mycop.plant<-NULL
Mycop.entry<-NULL
for (i in gnf$Name){
  mm<-ifelse(nn<-grep(i, Myco.po.ab$associatedTaxa, ignore.case = TRUE),print(i),0)
  Mycop.plant<-c(Mycop.plant,mm)
  Mycop.entry<-c(Mycop.entry,nn)
}
hits.taxa<-data.frame(Mycop.entry,Mycop.plant)

all.hits<-unique(rbind(hits.habitat,hits.substrate,hits.taxa))
all.hits<-all.hits[order(all.hits$Mycop.entry),]

# keep only hosts at species level and below.
word.count=NULL
for (i in all.hits$Mycop.plant){
  word.count=c(word.count, print(length(strsplit(i, " ")[[1]])))
}
all.hits<-data.frame(all.hits, word.count)

all.hits.sp <- all.hits %>% 
  filter(!all.hits$word.count=="1")

# Now we have fungi and associated plant species in one list.
Myco.po.ab.sp<-data.frame(all.hits.sp$Mycop.entry, all.hits.sp$Mycop.plant, Myco.po.ab[all.hits.sp$Mycop.entry,] ) 


#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# I. Here we create Myco.po.ab.sp, which will be used as raw data for the "Locally co-occurring" database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# saveRDS(Myco.po.ab.sp, 'intermediate_rds_csv/Mycoportal/Myco.po.ab.sp.rds')

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



#save everything that didn't detect a latin name at all (genus or species level). There are common names in here 
# that mostly can work for genus, with some for species (ie. balsam fir). There are also Russian names to reassess.
Mycop.undet<-Myco.po.ab[-(sort(unique(all.hits$Mycop.entry))),]
# saveRDS(Mycop.undet, 'intermediate_rds_csv/Mycoportal/Mycop.undet.rds')
# write.csv(Mycop.undet, 'OSF/Mycoportal/Mycoportal.undet.csv')


#4.
#### filter to everything that definitely grows on living plant host (ignoring "dead, decay, fallen, etc" search for now)

Myco.on.plant<-Myco.po.ab.sp %>%
  filter_at(vars(associatedTaxa,habitat,substrate),any_vars(grepl(
    'host|infect| on |^on |parasit|leaves|leaf|twig|stem|huesped| sur |^sur |"sur | auf |^auf |"auf | ad |^ad |"ad | ^in | in |"in |^an |"an | an |lebenden|anther| bark |branch|
      spot |spots |splotch|rust|streak|causing|inflorescence|fruit |fruits |pathogen|petiole|pod |pods |root |roots |cone |cones |covering |culm|ovar|foli|needle|caulibus|cortic|in den |flower|
      seed|pe frunze de|pe ramuri de|pe tulpini de|rama|smut |sobre | sub |^sub |blättern|blattern|stängeln|kapseln|trunk|living|live |vivant|
      Wurzeln|Fruchtstielen|Schoten|Halmen|Rinde | vine|scape| head|stalk|frond|Balttsticlen|wedel|ramis|caudex|caudice|viva|rama', (.), ignore.case = TRUE)))

Mycop.not.on.plant<-setdiff(Myco.po.ab.sp,Myco.on.plant)
write.csv(Mycop.not.on.plant,"~/Desktop/Mycop.not.on.plant.csv")
# look through these for usable data. If a plant name alone was in the substrate field, it was assumed this was growing on the plant itself, and these records were salvaged (lines 8-972).
# entries from Mycop.not.on.plant with question marks did not get coded from Russian correctly, but these were included in the GBIF dataset.
# The INEP:F collection entries from Russia were collected from Mycoportal on July 18, and these were searched by index in the Mycop.not.on.plant. Google Translate
# was used to make sure that these were growing on live or dead plant. Only one entry was removed, as it had a "soil" substrate.

#read in missing INEP.F data as 'inep'

Mycop.not.on.plant$adj.recordID<-gsub("urn:uuid:","",Mycop.not.on.plant$recordId)
russians<-inep[which(inep$recordID %in% Mycop.not.on.plant$adj.recordID),]
#these were located in Mycop.not.on.plant in Excel (59 entries, 1 filtered.)

# All manually checked, incl Russian translated additions apparently not in GBIF, were kept in the salvaged.from.Mycop.not.on.plant.csv file.
# Only a few entries with habitat data implying growth on plant were kept (most had empty habitat fields) and added to Myco.on.plant.

####
####
#######
#####
#####


soil<-c(grep('soil', (Myco.on.plant[,13]), ignore.case = TRUE),grep('soil', (Myco.on.plant[,15]),ignore.case = TRUE),grep('soil', (Myco.on.plant[,16]),ignore.case = TRUE))
litter<-c(grep('litter', (Myco.on.plant[,13]), ignore.case = TRUE),grep('litter', (Myco.on.plant[,15]),ignore.case = TRUE),grep('litter', (Myco.on.plant[,16]),ignore.case = TRUE))
under<-c(grep('under', (Myco.on.plant[,13]), ignore.case = TRUE),grep('under', (Myco.on.plant[,15]),ignore.case = TRUE),grep('under', (Myco.on.plant[,16]),ignore.case = TRUE))
dead<-c(grep('dead', (Myco.on.plant[,13]), ignore.case = TRUE),grep('dead', (Myco.on.plant[,15]),ignore.case = TRUE),grep('dead', (Myco.on.plant[,16]),ignore.case = TRUE))
decomp<-c(grep('decompos', (Myco.on.plant[,13]), ignore.case = TRUE),grep('decompos', (Myco.on.plant[,15]), ignore.case = TRUE),grep('decompos', (Myco.on.plant[,16]), ignore.case = TRUE))
stump<-c(grep('stump', (Myco.on.plant[,13]), ignore.case = TRUE),grep('stump', (Myco.on.plant[,15]), ignore.case = TRUE),grep('stump', (Myco.on.plant[,16]), ignore.case = TRUE))
fallen<-c(grep('fallen', (Myco.on.plant[,13]), ignore.case = TRUE),grep('fallen', (Myco.on.plant[,15]), ignore.case = TRUE),grep('fallen', (Myco.on.plant[,16]), ignore.case = TRUE))
lo<-c(grep('^log', (Myco.on.plant[,13]), ignore.case = TRUE),grep('^log', (Myco.on.plant[,15]), ignore.case = TRUE),grep('^log', (Myco.on.plant[,16]), ignore.case = TRUE))
logg<-c(grep(' log', (Myco.on.plant[,13]), ignore.case = TRUE),grep(' log', (Myco.on.plant[,15]), ignore.case = TRUE),grep(' log', (Myco.on.plant[,16]), ignore.case = TRUE))
dying<-c(grep('dying', (Myco.on.plant[,13]), ignore.case = TRUE),grep('dying', (Myco.on.plant[,15]), ignore.case = TRUE),grep('dying', (Myco.on.plant[,16]), ignore.case = TRUE))
muert<-c(grep('muert', (Myco.on.plant[,13]), ignore.case = TRUE),grep('muert', (Myco.on.plant[,15]), ignore.case = TRUE),grep('muert', (Myco.on.plant[,16]), ignore.case = TRUE))
rott<-c(grep('rott', (Myco.on.plant[,13]), ignore.case = TRUE),grep('rott', (Myco.on.plant[,15]), ignore.case = TRUE),grep('rott', (Myco.on.plant[,16]), ignore.case = TRUE))
decay<-c(grep('decay', (Myco.on.plant[,13]), ignore.case = TRUE),grep('decay', (Myco.on.plant[,15]), ignore.case = TRUE),grep('decay', (Myco.on.plant[,16]), ignore.case = TRUE))
duff<-c(grep('duff', (Myco.on.plant[,13]), ignore.case = TRUE),grep('duff', (Myco.on.plant[,15]), ignore.case = TRUE),grep('duff', (Myco.on.plant[,16]), ignore.case = TRUE))
ground<-c(grep('ground', (Myco.on.plant[,13]), ignore.case = TRUE),grep('ground', (Myco.on.plant[,15]), ignore.case = TRUE),grep('ground', (Myco.on.plant[,16]), ignore.case = TRUE))
emortu<-c(grep('emortu', (Myco.on.plant[,13]), ignore.case = TRUE),grep('emortu', (Myco.on.plant[,15]), ignore.case = TRUE),grep('emortu', (Myco.on.plant[,16]), ignore.case = TRUE))
downtree<-c(grep('down tree', (Myco.on.plant[,13]), ignore.case = TRUE),grep('down tree', (Myco.on.plant[,15]), ignore.case = TRUE),grep('down tree', (Myco.on.plant[,16]), ignore.case = TRUE))
baseoftree<-c(grep('base of tree', (Myco.on.plant[,13]), ignore.case = TRUE),grep('base of tree', (Myco.on.plant[,15]), ignore.case = TRUE),grep('base of tree', (Myco.on.plant[,16]), ignore.case = TRUE))
baseoftrunk<-c(grep('base of trunk', (Myco.on.plant[,13]), ignore.case = TRUE),grep('base of trunk', (Myco.on.plant[,15]), ignore.case = TRUE),grep('base of trunk', (Myco.on.plant[,16]), ignore.case = TRUE))
hypogeous<-c(grep('hypogeous', (Myco.on.plant[,13]), ignore.case = TRUE),grep('hypogeous', (Myco.on.plant[,15]), ignore.case = TRUE),grep('hypogeous', (Myco.on.plant[,16]), ignore.case = TRUE))
faulen<-c(grep('faulen', (Myco.on.plant[,13]), ignore.case = TRUE),grep('faulen', (Myco.on.plant[,15]), ignore.case = TRUE),grep('faulen', (Myco.on.plant[,16]), ignore.case = TRUE))
dung<-c(grep('dung', (Myco.on.plant[,13]), ignore.case = TRUE),grep('dung', (Myco.on.plant[,15]), ignore.case = TRUE),grep('dung', (Myco.on.plant[,16]), ignore.case = TRUE))
putri<-c(grep('putri', (Myco.on.plant[,13]), ignore.case = TRUE),grep('putri', (Myco.on.plant[,15]), ignore.case = TRUE),grep('putri', (Myco.on.plant[,16]), ignore.case = TRUE))
putre<-c(grep('putre', (Myco.on.plant[,13]), ignore.case = TRUE),grep('putre', (Myco.on.plant[,15]), ignore.case = TRUE),grep('putre', (Myco.on.plant[,16]), ignore.case = TRUE))

#For "on plant tissue" (not final dataset), keep several of those we are removing for the living database.
Myco.on.plant.remove<-Myco.on.plant[unique(c(soil, litter, under, duff, dung, ground)),]
Mop.remove2<-Myco.on.plant.remove[grep("^on|host:", Myco.on.plant.remove$associatedTaxa, ignore.case = TRUE, invert = TRUE),]
Myco.on.plant.FIN<-setdiff(Myco.on.plant, Mop.remove2)

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# II. Here we create Myco.on.plant.FIN, which will be used as raw data for the "On plant tissue" database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# saveRDS(Myco.on.plant.FIN, 'intermediate_rds_csv/Mycoportal/Myco.on.plant.FIN.rds')

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# Collect fungi on dead plants and miscellaneous (e.g., soil, base of tree).
deadplant<-unique(c(dead,decomp,stump,fallen,lo,logg,dying,muert,rott,decay,soil,litter,under,duff,ground,emortu,downtree,baseoftree,baseoftrunk,hypogeous,faulen,dung,putri,putre))
Mycop.deadplant<-Myco.on.plant[deadplant,]
Mycop.liveplant<-Myco.on.plant[-deadplant,]
write.csv(Mycop.liveplant, "OSF/Mycoportal/Mycop.liveplant.csv")
write.csv(Mycop.deadplant, "OSF/Mycoportal/Mycop.deadplant.csv")

## Manually edit liveplant entries.
#note: entries kept if they were explicitly stated to be the reason for rot of host plant.
#note: "effective" duplicates (more than one plant host erroneously detected per entry) retained for now.
#note: rhizosphere entries removed.
#note: if host species has ? (to signify uncertainty in identification), it was deleted
#note: decorticated kept. wood kept. host: "" kept, unless clear it was dead tissue. dried tissue of plant kept. sticks removed.
#note: if host was other fungi on plant, the plant host record was retained, and fungal host record dropped (as in Aecidi(a)(um))
#note: eliminated all badly abbreiviated entries (i.e., U. californica for Fragaria californica)
#note: deleted entries with abbreviated genus names "#. canadensis" if no full genus names provided in separate fields.
#note: records retained if associated.taxa had name of host plant but habitat did not, but was still suitable (i.e., "leaves of plant", but not "dead leaves of plant")
#note: removed charred and burn(ed)(t)
#note: If abbrevated genus provided with other full genus (Cornus alternifolia and C. sericea), I inputted full genus name.

#note: faulen, dung, putri, and putre were created on a semi-edited file in order to edit in R (see "deadplant/liveplant.edit" files), but manual editing resumed normally.
#note: as of july 20 if sorting mycop.edit.2.liveplant by substrate, then habitat, then associated taxa, 



##### PART III ##### Resolve names/apply consistency to all taxa. Final formatting.

# search against the International Plant Names Index
edit.liveplant<-read.csv("OSF/Mycoportal/Mycop.edit.liveplant.csv")
#format species lists
hosts.mycop.edit.liveplant<-as.character(unique(edit.liveplant$all.hits.sp.Mycop.plant))
plnamresolv<-gnr_resolve(hosts.mycop.edit.liveplant[1:999], preferred_data_sources = 167)
plnamresolv1<-gnr_resolve(hosts.mycop.edit.liveplant[1000:1999], preferred_data_sources = 167)
plnamresolv2<-gnr_resolve(hosts.mycop.edit.liveplant[2000:2999], preferred_data_sources = 167)
plnamresolv3<-gnr_resolve(hosts.mycop.edit.liveplant[3000:3999], preferred_data_sources = 167)
plnamresolv4<-gnr_resolve(hosts.mycop.edit.liveplant[4000:4999], preferred_data_sources = 167)
plnamresolv5<-gnr_resolve(hosts.mycop.edit.liveplant[5000:5999], preferred_data_sources = 167)
plnamresolv6<-gnr_resolve(hosts.mycop.edit.liveplant[6000:6999], preferred_data_sources = 167)
plnamresolv7<-gnr_resolve(hosts.mycop.edit.liveplant[7000:7999], preferred_data_sources = 167)
plnamresolv8<-gnr_resolve(hosts.mycop.edit.liveplant[8000:8714], preferred_data_sources = 167)
plnamresolv8.5<-gnr_resolve(hosts.mycop.edit.liveplant[8716:8999], preferred_data_sources = 167)
plnamresolv9<-gnr_resolve(hosts.mycop.edit.liveplant[9000:9999], preferred_data_sources = 167)
plnamresolv10<-gnr_resolve(hosts.mycop.edit.liveplant[10000:10999], preferred_data_sources = 167)
plnamresolv11<-gnr_resolve(hosts.mycop.edit.liveplant[11000:11999], preferred_data_sources = 167)
plnamresolv12<-gnr_resolve(hosts.mycop.edit.liveplant[12000:12999], preferred_data_sources = 167)
plnamresolv13<-gnr_resolve(hosts.mycop.edit.liveplant[13000:13999], preferred_data_sources = 167)
plnamresolv14<-gnr_resolve(hosts.mycop.edit.liveplant[14000:14999], preferred_data_sources = 167)
plnamresolv15<-gnr_resolve(hosts.mycop.edit.liveplant[15000:15999], preferred_data_sources = 167)
plnamresolv16<-gnr_resolve(hosts.mycop.edit.liveplant[16000:16999], preferred_data_sources = 167)
plnamresolv17<-gnr_resolve(hosts.mycop.edit.liveplant[17000:17999], preferred_data_sources = 167)
plnamresolv18<-gnr_resolve(hosts.mycop.edit.liveplant[18000:18704], preferred_data_sources = 167)
full.plant.host.name.list<-rbind.data.frame(plnamresolv, plnamresolv1, plnamresolv2, plnamresolv3, plnamresolv4, 
                                            plnamresolv5, plnamresolv6, plnamresolv7, plnamresolv8, plnamresolv8.5, plnamresolv9, 
                                            plnamresolv10, plnamresolv11, plnamresolv12, plnamresolv13, plnamresolv14, 
                                            plnamresolv15, plnamresolv16, plnamresolv17, plnamresolv18)
saveRDS(full.plant.host.name.list, "intermediate_rds_csv/Mycoportal/new.full.plant.host.name.list.Mycop.rds")

edit.liveplant$match.number<-match(edit.liveplant$all.hits.sp.Mycop.plant, full.plant.host.name.list$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
edit.liveplant$match.name<-print(full.plant.host.name.list$matched_name[edit.liveplant$match.number[1:131583]])
#write.csv(edit.liveplant, "OSF/Mycoportal/edit.liveplant.2.csv")

#remaining ones - use Tropicos (Missouri)
not.found<-edit.liveplant[is.na(edit.liveplant$match.name),]
nf.hosts.mycop.edit.liveplant<-as.character(unique(not.found$all.hits.sp.Mycop.plant))
plnamresolv19<-gnr_resolve(nf.hosts.mycop.edit.liveplant[1:1191], preferred_data_sources = 165)
plnamresolv20<-gnr_resolve(nf.hosts.mycop.edit.liveplant[1193:2470], preferred_data_sources = 165)
nf.full.plant.host.name.list<-rbind.data.frame(plnamresolv19,plnamresolv20)
#write.csv(nf.full.plant.host.name.list, "intermediate_rds_csv/Mycoportal/nf.full.plant.host.name.list.csv")
not.found$match.number<-match(not.found$all.hits.sp.Mycop.plant, nf.full.plant.host.name.list$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
not.found$match.name<-print(nf.full.plant.host.name.list$matched_name[not.found$match.number[1:5537]])
not.found2<-not.found[is.na(not.found$match.name),]
nf2.hosts.mycop.edit.liveplant<-as.character(unique(not.found2$all.hits.sp.Mycop.plant))

####double check slowly

#USDA PLANT db
plnamresolv21<-gnr_resolve(nf2.hosts.mycop.edit.liveplant[1:981], preferred_data_sources = 150)
plnamresolv22<-gnr_resolve(nf2.hosts.mycop.edit.liveplant[983:2032], preferred_data_sources = 150)
nf2.full.plant.host.name.list<-rbind.data.frame(plnamresolv21,plnamresolv22)
#write.csv(nf2.full.plant.host.name.list, "intermediate_rds_csv/Mycoportal/nf2.full.plant.host.name.list.csv")
not.found2$match.number<-match(not.found2$all.hits.sp.Mycop.plant, nf2.full.plant.host.name.list$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
not.found2$match.name<-print(nf2.full.plant.host.name.list$matched_name[not.found2$match.number[1:4320]])
not.found3<-not.found2[is.na(not.found2$match.name),]
nf3.hosts.mycop.edit.liveplant<-as.character(unique(not.found3$all.hits.sp.Mycop.plant))
#saveRDS(not.found3, "intermediate_rds_csv/Mycoportal/not.found3.mycop.rds")

#list of all automatically detected plant hosts 
extended.full.plant.list<-rbind(full.plant.host.name.list,plnamresolv19,plnamresolv20,plnamresolv21,plnamresolv22)


#filter out fungi (Index Fungorum source)
plnamresolv23<-gnr_resolve(nf3.hosts.mycop.edit.liveplant[1:966], preferred_data_sources = 5)
plnamresolv24<-gnr_resolve(nf3.hosts.mycop.edit.liveplant[968:1996], preferred_data_sources = 5)
nf3.full.fungi.host.name.list<-rbind.data.frame(plnamresolv23,plnamresolv24)

not.found3$match.number<-match(not.found3$all.hits.sp.Mycop.plant, nf3.full.fungi.host.name.list$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
not.found3$match.name<-print(nf3.full.fungi.host.name.list$matched_name[not.found3$match.number[1:4132]])
#list without fungi
no.fung<-not.found3[is.na(not.found3$match.name),]
#list of fungi
fung<-not.found3[!is.na(not.found3$match.name),]

nf4.hosts.mycop.edit.liveplant<-as.character(unique(no.fung$all.hits.sp.Mycop.plant))

nf<-no.fung %>%
  distinct(all.hits.sp.Mycop.plant, .keep_all = TRUE)

#Use Taxonstand to try to interpret some more (from The Plant List). This takes some time.
nf.TPL<-TPL(nf$all.hits.sp.Mycop.plant)

#save successful Taxonstand hits
nf.Taxon.hits<-filter(nf.TPL, New.Taxonomic.status != "")
#the rest
nf.Taxon.fail<-filter(nf.TPL, New.Taxonomic.status == "")

#Filter fungal genera
#list fung genera (from Index Fungorum, before)
fung.gen<-unique(word(fung$all.hits.sp.Mycop.plant, start=1L))
nf.gen.Taxon.fail<-filter(nf.Taxon.fail, !New.Genus %in% fung.gen)

#Correct genus names. If Genus names are corrected, the TPL and other tools perform better. These are primarily modified latin names, but some include easy-to-spot typos. All apparent genera typos with >1 entry were corrected. 
nf.gen.Taxon.fail$sing.genus<-gsub("Populi","Populus",gsub("Abietis","Abies",gsub("Aceris","Acer",
                                                                                  gsub("Aconiti","Aconitum",(gsub("Agropyrum","Agropyron",gsub("Agrostidis","Agrostis", gsub("Allii","Allium", gsub("Alni","Alnus",
                                                                                                                                                                                                    gsub("Baccharidis", "Baccharis", gsub("Bromi","Bromus",gsub("Calamogrostidis","Calamogrostis",gsub("Caricis","Carex", gsub("Carpini","Carpinum",
                                                                                                                                                                                                                                                                                                                               gsub("Clerodendron","Clerodendrum", gsub("Corni","Cornus",gsub("Crataegi","Crataegus", gsub("Cytisi", "Cytisum", gsub("Delphini","Delphinium", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                     gsub("Desmodii","Desmodium",gsub("Dipsaci", "Dipsacus",gsub("Elymi","Elymus",gsub("Euonymi","Euonymus", gsub("Eupatorii","Eupatorium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  gsub("Evonymi","Euonymus",gsub("Evonymus","Euonmyus",gsub("Fagi","Fagus",gsub("Fici","Ficus",gsub("Fraxini","Fraxinus",gsub("Galii","Galium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              gsub("Gei","Geum",gsub("Gelsemii","Gelsemium",gsub("Geranii","Gernaium",gsub("Heraclei","Heracleum",gsub("Hordei","Hordeum",gsub("Humuli","Humulus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               gsub("Ilicis","Ilex",gsub("Impatientis","Impatiens",gsub("Iridis","Iris",gsub("Juglandis","Juglans",gsub("Junci","Juncus",gsub("Juniperi","Juniperus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              gsub("Lappae","Lappa",gsub("Lathyri","Lathyrus",gsub("Lauri","Laurus",gsub("Leontodontis","Leontodon",gsub("Ligustici","Ligusticum",gsub("Lupini","Lupinus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       nf.gen.Taxon.fail$New.Genus))))))))))))))))))))))))))))))))))))))))))))))))
nf.gen.Taxon.fail$sing.genus<-gsub("Lycii","Lycium",gsub("Mori","Morus",gsub("Nasturtii","Nasturtium",gsub("Nerii","Nerium",gsub("Olerodendron","Clerodendron",
                                                                                                                                 gsub("Ononidis","Ononis",gsub("Orobi","Orobus",gsub("Osmorrhizae","Osmorhiza",gsub("Panici","Panicum",gsub("Pentstamon","Penstemon",gsub("Pentstemonis","Penstemon",
                                                                                                                                                                                                                                                                          gsub("Peucedani","Peucedanum",gsub("Pini","Pinus",gsub("Plantaginis","Plantago",gsub("Polygoni","Polygonum",gsub("Polypodium","Polypodii",gsub("Populi","Populus",
                                                                                                                                                                                                                                                                                                                                                                                                                         gsub("Pruni","Prunus",gsub("Pteridis","Pteridium",gsub("Pyri","Pyrus",gsub("Rhamni","Rhamnus",gsub("Rhois","Rhus",gsub("Rhoidis","Rhus",gsub("Ribis","Ribes",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      gsub("Roripa","Rorippa",gsub("Rubi","Rubus",gsub("Rumicis","Rumex",gsub("Salicis","Salex",gsub("Sambuci","Sambucus",gsub("Scirpi","Scirpus",gsub("Senecionis","Senecio",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       gsub("Silphinium","Silphium",gsub("Sisymbrii","Sisymbrium",gsub("Smilacis","Smilax",gsub("Solani","Solanum",gsub("Sonchi","Sonchus",gsub("Symphyti","Symphytum",gsub("Tami","Tamus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            gsub("Tanaceti","Tanacetum",gsub("Teucrii","Teucrium",gsub("Thalictri","Thalictrum",gsub("Thesii","Theisum",gsub("Trifolii","Trifolium",gsub("Tritici","Triticum",gsub("Ulmi","Ulmus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   gsub("Vaccinii","Vaccinium",gsub("Veronicae","Veronica",gsub("Viburni","Viburnum", gsub("Zizyphus","Ziziphus",nf.gen.Taxon.fail$sing.genus)))))))))))))))))))))))))))))))))))))))))))))))))

#compile corrected genera, paste with original specific epithets.
nf.gen.Taxon.fail$Corr.name<-paste(nf.gen.Taxon.fail$sing.genus, nf.gen.Taxon.fail$New.Species)
#rerun TPS on corrected data, hope for less manual work!
re.TPL<-TPL(nf.gen.Taxon.fail$Corr.name)
saveRDS(re.TPL, "intermediate_rds_csv/Mycoportal/re.TPL.rds")
#successful taxon hits
re.TPL.taxon.hits<-filter(re.TPL, New.Taxonomic.status != "")
#taxon fails
re.TPL.taxon.fails<-filter(re.TPL, New.Taxonomic.status == "")

#WorldFlora pkg
#WFO.download(save.dir = "~/OneDrive - UBC/PhD.Year.Three", WFO.remember = TRUE)
#WFO.remember()
re.WFO<-WFO.match(re.TPL.taxon.fails$Taxon, WFO.data=WFO.data)
re.WFO.read<-data.frame(cbind(re.WFO$spec.name, re.WFO$Fuzzy.dist, re.WFO$Old.name, re.WFO$scientificName, re.WFO$taxonomicStatus))

#write.csv(re.WFO.read, "intermediate_rds_csv/Mycoportal/re.WFO.read.csv")#manually edit. 
#editing:
#I started typing all genera that did not come up with old or accepted alternatives in world flora. If there was a clear solution, I added this to the "manual" column. 
#I checked all entries to make sure they made sense, and searched world flora if only a genus determination was made.
#All entries that were clearly incorrect were removed. Verged on over-removal rather than the alternative.
#removed entries where single old name forms 2 currently accepted names (e.g., Centaurea glauca is now Amberboa glauca/Centaurea calocephala) UNLESS there is accepted name that is the same but correct spelling.
#deleted entries where genus clearly did not match.
#deleted entries where tie was convincing ("Lachia pulchella" (sic) is Lawia pulchella or Lechea pulchella)
#changed variety level classification to species level.
#If only genus was detected, I looked in World Flora to see if there that genus had an entry with a similar specific epithet, and added this in "manual".

re.WFO.r.edit<-read.csv("intermediate_rds_csv/Mycoportal/re.WFO.read.edited.csv")

#Now compile a list with all the hosts

#Mycoportal
edit.liveplant$new.ID<-"NA"
edit.liveplant$namecheck.db<-"Intl.plant.names.index"
edit.liveplant$score<-print(full.plant.name.list$score[edit.liveplant$match.number[1:131583]])
Mycop.list.1<-edit.liveplant[!is.na(edit.liveplant$match.name),]
not.found$new.ID<-"NA"
not.found$namecheck.db<-"Tropicos"
not.found$score<-print(nf.full.plant.host.name.list$score[not.found$match.number[1:5537]])
Mycop.list.2<-not.found[!is.na(not.found$match.name),]
not.found2$new.ID<-"NA"
not.found2$namecheck.db<-"USDA.plant.db"
not.found2$score<-print(nf2.full.plant.host.name.list$score[not.found2$match.number[1:4320]])
Mycop.list.3<-not.found2[!is.na(not.found2$match.name),]
list.taxize.mycop<-rbind.data.frame(Mycop.list.1, Mycop.list.2, Mycop.list.3)
#prepare full file. this can be included as supplementary or for reference.
#write.csv(list.taxize.mycop, "OSF/Mycoportal/list.taxize.mycop.csv")
#prepare abbreviated file to be standardized and added to TPL and WorldFlora files.
list.taxize.mycop$X.1<-NULL

#part 2: prepare TPL output.
TPL.list<-rbind.data.frame(nf.Taxon.hits, re.TPL.taxon.hits)
#write.csv(TPL.list, "intermediate_rds_csv/Mycoportal/TPL.list.csv")
#add to match.name, new record ID, database name 
list.TPL.mycop$namecheck.db<-"The.plant.list"

#follow same format as list.taxize.mycop
no.fung$match.number<-match(no.fung$all.hits.sp.Mycop.plant, TPL.list$Taxon, nomatch = NA_integer_, incomparables = NULL)
TPL.list$match.name<-paste(TPL.list$New.Genus, TPL.list$New.Hybrid.Marker, TPL.list$New.Species,TPL.list$New.Authority)
no.fung$match.name<-print(TPL.list$match.name[no.fung$match.number[1:3385]])
no.fung$new.ID<-print(TPL.list$New.ID[no.fung$match.number[1:3385]])
list.TPL.mycop<-no.fung[!is.na(no.fung$match.name),]
#remove extra spaces
list.TPL.mycop$match.name<-gsub("  "," ",list.TPL.mycop$match.name)
list.TPL.mycop$match.name<-gsub("  "," ",list.TPL.mycop$match.name)
list.TPL.mycop$namecheck.db<-"The.plant.list"
list.TPL.mycop$score<-"NA"
#write.csv(list.TPL.mycop, "intermediate_rds_csv/Mycoportal/list.TPL.mycop.csv")
#filter this. removed errant entries in match name column.
#(From Methods, though appears earlier: "Manual filtering and formatting to match other data sources was conducted on the Taxonstand output." )
#saved as intermediate_rds_csv/Mycoportal/filt.list.TPL.mycop.csv


#part 3: prepare WorldFlora output
#include all data fields
re.WFO.r.edit.full<-re.WFO[re.WFO.r.edit$X,]
#add in manual input
re.WFO.r.edit.full$manual.name[!is.na(re.WFO.r.edit$manual.name)]<-re.WFO.r.edit$manual.name
re.WFO.r.edit.full$manual.auth[!is.na(re.WFO.r.edit$manual.auth)]<-re.WFO.r.edit$manual.auth
remaining.mycop<-no.fung[is.na(no.fung$match.name),]
remaining.mycop$match.number<-match(remaining.mycop$all.hits.sp.Mycop.plant, re.WFO.r.edit.full$spec.name, nomatch = NA_integer_, incomparables = NULL)
remaining.mycop$match.name<-print(re.WFO.r.edit.full$scientificName[remaining.mycop$match.number[1:2163]])
remaining.mycop$match.auth<-print(re.WFO.r.edit.full$scientificNameAuthorship[remaining.mycop$match.number[1:2163]])
remaining.mycop$new.ID<-print(re.WFO.r.edit.full$tplID[remaining.mycop$match.number[1:2163]])
remaining.mycop$namecheck.db<-"WorldFlora"
remaining.mycop$score<-"NA"
remaining.mycop$manual.name<-print(re.WFO.r.edit.full$manual.name[remaining.mycop$match.number[1:2163]])
remaining.mycop$manual.auth<-print(re.WFO.r.edit.full$manual.auth[remaining.mycop$match.number[1:2163]])
#this file retains which were manually changed and which were changed in R.
#saveRDS(remaining.mycop, "intermediate_rds_csv/Mycoportal/remaining.mycop.rds")

remaining.mycop$manual.name[is.na(remaining.mycop$manual.name)] <-''
remaining.mycop$manual.auth[is.na(remaining.mycop$manual.auth)] <-''
remaining.mycop$match.name[is.na(remaining.mycop$match.name)] <-''
remaining.mycop$match.auth[is.na(remaining.mycop$match.auth)] <-''
#combine names and authorities
remaining.mycop$match.name<-paste(remaining.mycop$match.name, remaining.mycop$match.auth)
#replace match name with manual name if a manual name is detected.
remaining.mycop$match.name[!remaining.mycop$manual.name==""]<-remaining.mycop$manual.name[!remaining.mycop$manual.name==""]
remaining.mycop$match.name[!remaining.mycop$manual.name==""]<-paste(remaining.mycop$match.name[!remaining.mycop$manual.name==""],remaining.mycop$manual.auth[!remaining.mycop$manual.name==""])
remaining.mycop$match.auth<-NULL
remaining.mycop$manual.name<-NULL
remaining.mycop$manual.auth<-NULL

#filter this. 

mycop.cleaned<-rbind.data.frame(list.taxize.mycop, filt.list.TPL.mycop, filt.remaining.mycop)
#saveRDS(mycop.cleaned, "intermediate_rds_csv/Mycoportal/mycop.cleaned.rds")






###########PHASE 2 : PREP FOR LINK PREDICTIONS ON LIVE PLANT TISSUE.

#additional cleaning for LIVING TISSUE ONLY!
##IMPORT WORKING LIST FROM Rashmi et al data

mycosphere.W.df<-readRDS("intermediate_rds_csv/Rashmi_and_Globi/myc.as3_MYCOSPHERE.NEW.4.2024.rds")
#for now, change names of columns
mycosphere.W.df$all.hits.sp.Mycop.plant<-mycosphere.W.df$Host
mycosphere.W.df$scientificName<-mycosphere.W.df$Fungi

#which entries are in mycosphere report data?
mycop.cleaned.SPHEREFILT<-anti_join(mycop.cleaned,mycosphere.W.df, by=c("all.hits.sp.Mycop.plant","scientificName"))

# find entries that seem to have duplicated "plant species", so these may be filtered.
mycop.cleaned.SPHEREFILT$first.letters<-substr(mycop.cleaned.SPHEREFILT$match.name, 1, 6)
mycop.cleaned2<-mycop.cleaned.SPHEREFILT %>% 
  group_by(all.hits.sp.Mycop.entry, first.letters) %>% 
  distinct(first.letters, .keep_all=TRUE)
mycop.cleaned3<-data.frame(mycop.cleaned2)


#remove lichens
#read in lichenlist data from "Consortium of Lichen Herbaria" https://lichenportal.org/portal/checklists/checklist.php?clid=1492&pid=558
lichenlist<-read.csv("ext_data/World Checklist of Genera of Lichenized Fungi_1698006802.csv")
#make genus column
lichenlist$genus<-gsub(" sp.", "", lichenlist$ScientificName)
#filter lichenizing fungi
mycop.cleaned.nl<-mycop.cleaned3[which(!mycop.cleaned3$genus %in% lichenlist$genus),]
#remove fungal taxon entries that contain word root "lichen"
nlmy<-mycop.cleaned.nl[grep("lichen",mycop.cleaned.nl$scientificName, ignore.case=TRUE),]
mycop.cleaned.nl2<-mycop.cleaned.nl[which(!mycop.cleaned.nl$scientificName %in% nlmy$scientificName),]

#remove "on [tree species]" entries
#start w/ making plant genus field
mc.nl2.genus<-data.frame(str_split_fixed(mycop.cleaned.nl2$match.name, " ", 2))
mycop.cleaned.nl2$plant.genus<-mc.nl2.genus$X1
#I retrieved the csv for the GlobalTreeSearch data on Oct 31 2023. https://rpubs.com/Roeland-KINDT/812716
#read this data in 
globtree<-read.csv("ext_data/global_tree_search_trees_1_7.csv")
#filter by genus in order to eliminate as many shrubby genera, as well.
gltrgenus<-data.frame(str_split_fixed(globtree$TaxonName, " ", 2))
globtree$genus<-gltrgenus$X1
m.on.unkept.match<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("On", mycop.cleaned.nl2$match.name))),])
m.on.unkept.match2<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("on", mycop.cleaned.nl2$match.name))),])
m.on.unkept.match3<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("On", mycop.cleaned.nl2$all.hits.sp.Mycop.plant))),])
m.on.unkept.match4<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("on", mycop.cleaned.nl2$all.hits.sp.Mycop.plant))),])
m.on.unkept.match5<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("On ", mycop.cleaned.nl2$all.hits.sp.Mycop.plant, ".", sep = ""))),])
m.on.unkept.match6<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("on ", mycop.cleaned.nl2$all.hits.sp.Mycop.plant, ".", sep = ""))),])
m.on.unkept.match7<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("On ", mycop.cleaned.nl2$match.name, ".", sep = ""))),])
m.on.unkept.match8<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("on ", mycop.cleaned.nl2$match.name, ".", sep = ""))),])
m.auf.unkept.match<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("Auf", mycop.cleaned.nl2$match.name))),])
m.auf.unkept.match2<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("auf", mycop.cleaned.nl2$match.name))),])
m.auf.unkept.match3<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("Auf", mycop.cleaned.nl2$all.hits.sp.Mycop.plant))),])
m.auf.unkept.match4<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("Auf ", mycop.cleaned.nl2$all.hits.sp.Mycop.plant, ".", sep = ""))),])
m.auf.unkept.match5<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("auf", mycop.cleaned.nl2$all.hits.sp.Mycop.plant))),])
m.auf.unkept.match6<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("Auf ", mycop.cleaned.nl2$match.name, ".", sep = ""))),])
m.ad.unkept.match<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("Ad", mycop.cleaned.nl2$match.name))),])
m.ad.unkept.match2<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("ad", mycop.cleaned.nl2$match.name))),])
m.ad.unkept.match3<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("Ad", mycop.cleaned.nl2$all.hits.sp.Mycop.plant))),])
m.ad.unkept.match4<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$habitat %in% (paste("ad", mycop.cleaned.nl2$all.hits.sp.Mycop.plant))),])
mas.on.unkept.match<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$associatedTaxa %in% (paste("On", mycop.cleaned.nl2$match.name))),])
mas.on.unkept.match2<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$associatedTaxa %in% (paste("on", mycop.cleaned.nl2$match.name))),])
mas.on.unkept.match3<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$associatedTaxa %in% (paste("On", mycop.cleaned.nl2$all.hits.sp.Mycop.plant))),])
mas.on.unkept.match4<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$associatedTaxa %in% (paste("on", mycop.cleaned.nl2$all.hits.sp.Mycop.plant))),])
mas.on.unkept.match5<-data.frame(mycop.cleaned.nl2[which(mycop.cleaned.nl2$associatedTaxa %in% (paste("on ", mycop.cleaned.nl2$match.name, ".", sep = ""))),])
m.unkept.on<-rbind.data.frame(m.on.unkept.match, m.on.unkept.match2, m.on.unkept.match3, m.on.unkept.match4, m.on.unkept.match5, m.on.unkept.match6, m.on.unkept.match7, m.on.unkept.match8,
                              m.auf.unkept.match, m.auf.unkept.match2,m.auf.unkept.match3, m.auf.unkept.match4, m.auf.unkept.match5, m.auf.unkept.match6, m.ad.unkept.match, m.ad.unkept.match2, 
                              m.ad.unkept.match3, m.ad.unkept.match4, mas.on.unkept.match, mas.on.unkept.match2, mas.on.unkept.match3, mas.on.unkept.match4, mas.on.unkept.match5)

#filter out tree data from the "on [species]" list
m.unkept.on.tree<-m.unkept.on[which(m.unkept.on$plant.genus %in% globtree$genus),]
m.unkept.on.no.tree<-m.unkept.on[which(!m.unkept.on$plant.genus %in% globtree$genus),]
#no "on [tree sp.]" entries data
mycop.cleaned.nt<-anti_join(mycop.cleaned.nl2, m.unkept.on.tree)
# 
mcn.G<-setdiff(mycop.cleaned.nt, mycop.cleaned.nt[grep(" on trunk of ", mycop.cleaned.nt$habitat, ignore.case=TRUE),])
mcn.Ga<-setdiff(mcn.G, mcn.G[grep("langui", mcn.G$habitat, ignore.case=TRUE),])
mcn.Gb<-setdiff(mcn.Ga, mcn.Ga[grep("root", mcn.Ga$habitat, ignore.case=TRUE),])
mcn.Gc<-setdiff(mcn.Gb, mcn.Gb[grep("old ", mcn.Gb$habitat, ignore.case=TRUE),])
mcn.Gd<-setdiff(mcn.Gc, mcn.Gc[grep("wood ", mcn.Gc$habitat, ignore.case=TRUE),])
mcn.Ge<-setdiff(mcn.Gd, mcn.Gd[grep("cortic", mcn.Gd$habitat, ignore.case=TRUE),])
mcn.G0<-setdiff(mcn.Ge, mcn.Ge[grep("saprophyt", mcn.Gd$associatedTaxa, ignore.case=TRUE),])
mcn.G1<-setdiff(mcn.G0, mcn.G0[grep("dried ", mcn.G0$associatedTaxa, ignore.case=TRUE),])
mcn.G2<-setdiff(mcn.G1, mcn.G1[grep("dry leaves of ", mcn.G1$habitat, ignore.case=TRUE),])
mcn.G3<-setdiff(mcn.G2, mcn.G2[grep("dry twigs of ", mcn.G2$habitat, ignore.case=TRUE),])
mcn.G4<-setdiff(mcn.G3, mcn.G3[grep("trunk", mcn.G3$habitat, ignore.case=TRUE),])
mcn.G5<-setdiff(mcn.G4, mcn.G4[grep("bark", mcn.G4$habitat, ignore.case=TRUE),])
mcn.G6<-setdiff(mcn.G5, mcn.G5[grep(" arid", mcn.G5$habitat, ignore.case=TRUE),])
mcn.G7<-setdiff(mcn.G6, mcn.G6[grep("sicc", mcn.G6$habitat, ignore.case=TRUE),])
mcn.G8<-setdiff(mcn.G7, mcn.G7[grep("mort", mcn.G7$habitat, ignore.case=TRUE),])
mcn.G9<-setdiff(mcn.G8, mcn.G8[grep("dürren", mcn.G8$habitat, ignore.case=TRUE),])
mcn.G10<-setdiff(mcn.G9, mcn.G9[grep("last year", mcn.G9$habitat, ignore.case=TRUE),])
mcn.G11<-setdiff(mcn.G10, mcn.G10[grep("dry stems ", mcn.G10$habitat, ignore.case=TRUE),])
mcn.G12<-setdiff(mcn.G11, mcn.G11[grep("withering", mcn.G11$habitat, ignore.case=TRUE),])
mcn.G13<-setdiff(mcn.G12, mcn.G12[grep("overwintered", mcn.G12$habitat, ignore.case=TRUE),])
mcn.G14<-setdiff(mcn.G13, mcn.G13[grep("previous year", mcn.G13$habitat, ignore.case=TRUE),])
mcn.G15<-setdiff(mcn.G14, mcn.G14[grep("dry stalk", mcn.G14$habitat, ignore.case=TRUE),])
mcn.G16<-setdiff(mcn.G15, mcn.G15[grep("dry branch", mcn.G15$habitat, ignore.case=TRUE),])
mcn.G17<-setdiff(mcn.G16, mcn.G16[grep("durren", mcn.G16$habitat, ignore.case=TRUE),])
mcn.G18<-setdiff(mcn.G17, mcn.G17[grep("bark ", mcn.G17$habitat, ignore.case=TRUE),])

#Keep
mcn.k1<-data.frame(mcn.G18[grep("leaves of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k2<-data.frame(mcn.G18[grep("n the leaves of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k3<-data.frame(mcn.G18[grep("n living leaves of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k4<-data.frame(mcn.G18[grep("n flowers of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k5<-data.frame(mcn.G18[grep("n seeds of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k6<-data.frame(mcn.G18[grep("n fruits of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k7<-data.frame(mcn.G18[grep("n inflor", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k8<-data.frame(mcn.G18[grep("n culms of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k9<-data.frame(mcn.G18[grep("n living twigs of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k10<-data.frame(mcn.G18[grep("n living stems of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k11<-data.frame(mcn.G18[grep("n ovaries of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k12<-data.frame(mcn.G18[grep("ad folia", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k13<-data.frame(mcn.G18[grep("n leaves and stems of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k14<-data.frame(mcn.G18[grep("n leaves and petioles of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k15<-data.frame(mcn.G18[grep("n perigyni", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k16<-data.frame(mcn.G18[grep("den Blättern von", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k17<-data.frame(mcn.G18[grep("Sur les feuilles", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k18<-data.frame(mcn.G18[grep("n anthers", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k19<-data.frame(mcn.G18[grep("n the anthers", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k20<-data.frame(mcn.G18[grep("rust", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k23<-data.frame(mcn.G18[grep("n stems of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k24<-data.frame(mcn.G18[grep("n stems and leaves of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k25<-data.frame(mcn.G18[grep("n stem of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k26<-data.frame(mcn.G18[grep("n stalk of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k27<-data.frame(mcn.G18[grep("n stalks of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k29<-data.frame(mcn.G18[grep("n seedlings of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k30<-data.frame(mcn.G18[grep("n petioles of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k33<-data.frame(mcn.G18[grep("n live leaves", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k34<-data.frame(mcn.G18[grep("n living and wilt", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k35<-data.frame(mcn.G18[grep("n leaves, sheaths", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k37<-data.frame(mcn.G18[grep("n leaves, stems", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k38<-data.frame(mcn.G18[grep("n leaves, petioles", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k39<-data.frame(mcn.G18[grep("n leaves & stems of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k40<-data.frame(mcn.G18[grep("n leaves and culms of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k41<-data.frame(mcn.G18[grep("n leaves and bracts of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k42<-data.frame(mcn.G18[grep("n fruit of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k43<-data.frame(mcn.G18[grep("n fronds of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k44<-data.frame(mcn.G18[grep("n female catkins of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k45<-data.frame(mcn.G18[grep("n culms and leaves of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k46<-data.frame(mcn.G18[grep("n canes of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k49<-data.frame(mcn.G18[grep("folia viva", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k50<-data.frame(mcn.G18[grep("in foliis", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k52<-data.frame(mcn.G18[grep("in floribus", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k53<-data.frame(mcn.G18[grep("in antheris", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k54<-data.frame(mcn.G18[grep("^leaves of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k55<-data.frame(mcn.G18[grep("leaf spot", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k57<-data.frame(mcn.G18[grep("leaf and stem gall", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k58<-data.frame(mcn.G18[grep("in spikelets", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k59<-data.frame(mcn.G18[grep("in petiolis", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k60<-data.frame(mcn.G18[grep("in ovariis", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k63<-data.frame(mcn.G18[grep("in fol.", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k64<-data.frame(mcn.G18[grep("In den Bluten von", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k65<-data.frame(mcn.G18[grep("In den Antheren von", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k66<-data.frame(mcn.G18[grep("n leaves and bracts of", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k67<-data.frame(mcn.G18[grep("Auf stengeln von", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k68<-data.frame(mcn.G18[grep("lebenden Blättern", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k69<-data.frame(mcn.G18[grep("Auf Asten von", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k71<-data.frame(mcn.G18[grep("ad caules", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k72<-data.frame(mcn.G18[grep("Auf Blättern", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k73<-data.frame(mcn.G18[grep("Auf Blattern", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k74<-data.frame(mcn.G18[grep("n stammen von", mcn.G18$habitat, ignore.case=TRUE),])
mcn.k75<-data.frame(mcn.G18[grep("living leaves", mcn.G18$substrate, ignore.case=TRUE),])


mcn.keep<-rbind.data.frame(mcn.k1, mcn.k2, mcn.k3, mcn.k4, mcn.k5, mcn.k6, mcn.k7, mcn.k8, mcn.k9, mcn.k10, mcn.k11, mcn.k12, mcn.k13,
                           mcn.k14, mcn.k15, mcn.k16, mcn.k17, mcn.k18, mcn.k19, mcn.k20, mcn.k23, mcn.k24, mcn.k25, mcn.k26, 
                           mcn.k27,mcn.k29, mcn.k30, mcn.k33, mcn.k34, mcn.k35, mcn.k37, mcn.k38, mcn.k39, mcn.k40,
                           mcn.k41, mcn.k42, mcn.k43, mcn.k44, mcn.k45, mcn.k46, mcn.k49, mcn.k50, mcn.k52, mcn.k53,
                           mcn.k54, mcn.k55, mcn.k57, mcn.k58, mcn.k59, mcn.k60, mcn.k63, mcn.k64, mcn.k65, mcn.k66, 
                           mcn.k67, mcn.k68, mcn.k69, mcn.k71, mcn.k72, mcn.k73, mcn.k74, mcn.k75)


#### make df of distinct fungi/plant combos THAT WE WANT TO KEEP ("singles" list)
sm.1<-mcn.keep %>%
  group_by(scientificName, all.hits.sp.Mycop.plant) %>%
  arrange(str_count(habitat), .by_group=TRUE) %>%
  distinct(scientificName, .keep_all=TRUE)
sing.m.1<-data.frame(sm.1)
#list of plants/fungi
sing.m.1b<-sing.m.1[,c(3,9)]
#these are copies of the extras that were "kept"
m.kept.S.extras<-setdiff(mcn.keep, sing.m.1)

# which occurrences were not detected between living plant and fungi (sing.HGTdf.all)?
undet.mcn1<-setdiff(mcn.G18, mcn.keep)

####### New filtration scheme for remaining entries. Restrict tree entries that are not specific, and just say "host".
# subset to only habitat with "NA"
undet.mcn.na.hab<-undet.mcn1[which(is.na(undet.mcn1$habitat)),]
#of these(w/ no other hab info), which entries just say "Host: [plant sp.]" 
undet.mcn.nh1.nohost<-data.frame(undet.mcn.na.hab[which(undet.mcn.na.hab$associatedTaxa %in% (paste("host:", undet.mcn.na.hab$all.hits.sp.Mycop.plant))),])
undet.mcn.nh2.nohost<-data.frame(undet.mcn.na.hab[which(undet.mcn.na.hab$associatedTaxa %in% (paste("host:", undet.mcn.na.hab$match.name))),])
undet.mcn.nhA<-rbind.data.frame(undet.mcn.nh1.nohost, undet.mcn.nh2.nohost)
undet.mcn.tree<-undet.mcn.nhA[which(undet.mcn.nhA$plant.genus %in% globtree$genus),]
#drop these to restrict saprophytes, etc
undet.mcn.2<-setdiff(undet.mcn1, undet.mcn.tree)

# do the same for when habitat matches match.name
undet.mcn.hab.pl<-undet.mcn.2[which(undet.mcn.2$habitat == undet.mcn.2$match.name),]
undet.mcn.as.pl<-undet.mcn.2[which(undet.mcn.2$habitat == undet.mcn.2$all.hits.sp.Mycop.plant),]
undet.mcn.mat<-rbind.data.frame(undet.mcn.hab.pl, undet.mcn.as.pl)
undet.mcn.mat1.nohost<-data.frame(undet.mcn.mat[which(undet.mcn.mat$associatedTaxa %in% (paste("host:", undet.mcn.mat$all.hits.sp.Mycop.plant))),])
undet.mcn.mat2.nohost<-data.frame(undet.mcn.mat[which(undet.mcn.mat$associatedTaxa %in% (paste("host:", undet.mcn.mat$match.name))),])
undet.mcn.matA<-rbind.data.frame(undet.mcn.mat1.nohost, undet.mcn.mat2.nohost)
undet.mcn.mat.tree<-undet.mcn.matA[which(undet.mcn.matA$plant.genus %in% globtree$genus),]
#drop these again to restrict saprophytes, etc
undet.mcn.2b<-setdiff(undet.mcn.2, undet.mcn.mat.tree)

#new section. keep these. 
host.mG18<-rbind(mcn.G18[which(mcn.G18$associatedTaxa %in% (paste("host:", undet.mcn.mat$match.name))),],
                 mcn.G18[which(mcn.G18$associatedTaxa %in% (paste("host:", undet.mcn.mat$all.hits.sp.Mycop.plant))),])
mcn.k76<-data.frame(host.mG18[grep("leaves", host.mG18$habitat, ignore.case=TRUE),])
mcn.k77<-data.frame(host.mG18[grep("leaf", host.mG18$habitat, ignore.case=TRUE),])
mcn.k78<-data.frame(host.mG18[grep("fruit", host.mG18$habitat, ignore.case=TRUE),])
mcn.k79<-data.frame(host.mG18[grep("stem", host.mG18$habitat, ignore.case=TRUE),])
mcn.k80<-data.frame(host.mG18[grep("flower", host.mG18$habitat, ignore.case=TRUE),])
mcn.k81<-data.frame(host.mG18[grep("feuille", host.mG18$habitat, ignore.case=TRUE),])
mcn.k82<-data.frame(host.mG18[grep("petiol", host.mG18$habitat, ignore.case=TRUE),])
mcn.k83<-data.frame(host.mG18[grep("foliis", host.mG18$habitat, ignore.case=TRUE),])
mcn.k84<-data.frame(host.mG18[grep("fructis", host.mG18$habitat, ignore.case=TRUE),])
mcn.k85<-data.frame(host.mG18[grep("seed", host.mG18$habitat, ignore.case=TRUE),])
new.mcn.keep<-rbind.data.frame(mcn.k76, mcn.k77, mcn.k78, mcn.k79, mcn.k80, mcn.k81, mcn.k82,mcn.k83,mcn.k84,mcn.k85)
mcn.h.si<-new.mcn.keep %>%
  group_by(scientificName, all.hits.sp.Mycop.plant) %>%
  arrange(str_count(habitat), .by_group=TRUE) %>%
  distinct(scientificName, .keep_all=TRUE)
mcnh.si.df<-data.frame(mcn.h.si)
#leftovers
undet.mcn.3<-anti_join(undet.mcn.2b, new.mcn.keep)

undetmcn3.singles<-undet.mcn.3 %>%
  distinct(scientificName, match.name, habitat, associatedTaxa, substrate, .keep_all = TRUE)
#split this into tree and nontree sp.
tree.undetmcn3.singles<-undetmcn3.singles[which(undetmcn3.singles$plant.genus %in% globtree$genus),]
#write.csv(tree.undetmcn3.singles, "intermediate_rds_csv/Mycoportal/tree.undetmcn3.singles.csv")
#Keep the good ones. These are good
Filt.tree.uds<-read.csv("intermediate_rds_csv/Mycoportal/FILT.tree.undetmcn3.singles.csv")
Filt.tree.uds<-Filt.tree.uds[,-1]

#nontree
ntree.undetmcn3.singles<-undetmcn3.singles[which(!undetmcn3.singles$plant.genus %in% globtree$genus),]
#still many combinations, so I will search through unique combinations of host/fungi, and delete ones not occurring that are clearly occurring endosymbiotically on living plant tissue.
ntree.undetmcn3.distinct.Pl.F<-ntree.undetmcn3.singles %>%
  distinct(scientificName, all.hits.sp.Mycop.plant, .keep_all = TRUE)

#write.csv(ntree.undetmcn3.distinct.Pl.F, "intermediate_rds_csv/Mycoportal/ntree.undetmcn3.distinct.Pl.F.csv")
#Go through these. Delete the bad ones, keep the good ones for final data. See if the bad ones are duplicated.
#Filtered out branches, wood, dead, burned, dry, etc.
#These are good.
Filt.ntr.undet.dPF<-read.csv("intermediate_rds_csv/Mycoportal/FILT.ntree.undetmcn3.distinct.Pl.F.csv")
Filt.ntr.undet.dPF<-Filt.ntr.undet.dPF[,-1]

#any other entries representing plant fungi combos that were filtered out?
ntree.udmdPF4<-anti_join(ntree.undetmcn3.singles, Filt.ntr.undet.dPF, by = c("scientificName","all.hits.sp.Mycop.plant"))
#Keep the good ones, delete duplicates
Filt.ntr.undet.dPF2<-read.csv("intermediate_rds_csv/Mycoportal/FILT.ntree.udmdPF4.csv")
Filt.ntr.undet.dPF2<-Filt.ntr.undet.dPF2[,-1]
Mycopo.wd<-rbind.data.frame(sing.m.1, Filt.tree.uds, Filt.ntr.undet.dPF, Filt.ntr.undet.dPF2, mcn.h.si)
Mycopo.wdraw<-Mycopo.wd %>%
  group_by(scientificName, all.hits.sp.Mycop.plant) %>%
  distinct(scientificName, .keep_all=TRUE)
Mycopo.wdraw<-data.frame(Mycopo.wdraw)

#there are still match names that need to be fixed. Find and remove entries that have 0-1 words in the match.name field.
Mycopo.no.match.n<-Mycopo.wdraw[which(lengths(strsplit(Mycopo.wdraw$match.name, ' '))<2),]
#remove them from sample, we will add cleaned entries in later.
Mycopo.wdraw<-anti_join(Mycopo.wdraw, Mycopo.no.match.n)

Mycopo.no.match.n$match.cor<-gsub("Aconiti","Aconitum", gsub("Abietis","Abies",gsub("Aceris","Acer",gsub("Acgopodio","Aegopodium", gsub("Aira","Aera",gsub("Aconiti","Aconitum",gsub("Agroparum","Agropyron",gsub("Agropyrum","Agropyron",gsub("Agrostidis","Agrostis", gsub("Aletridis","Aletris",gsub("Allii","Allium", gsub("Alni","Alnus",gsub("Anthirrhini","Antirrhinum", gsub("Antirrhini","Antirrhinum", gsub("Ari ","Arum ",  gsub("Arundinis","Arundo",gsub("Aspidii","Aspidium",gsub("Asparagi","Asparagus",gsub("Baccharidis", "Baccharis", gsub("Betulae", "Betula",gsub("Betulus", "Betula",gsub("Bromi","Bromus",gsub("Calamogrostidis","Calamogrostis",gsub("Caricis","Carex", gsub("Carpini","Carpinum",gsub("Cauliflower-brassica","Brassica", gsub("Cissi","Cissus", gsub("Chrysanthemi","Chrysanthemum",  gsub("Clerodendron","Clerodendrum",gsub("Clerodendri","Clerodendrum", gsub("Corni","Cornus",gsub("Crataegi","Crataegus", gsub("Cynanchi", "Cynanchum", gsub("Cytisi", "Cytisum", gsub("Delphinii","Delphinium",gsub("Desmodii","Desmodium",gsub("Dipsaci", "Dipsacus",gsub("Elymi","Elymus",gsub("Euonymi","Euonymus", gsub("Eupatorii","Eupatorium",gsub("Echii","Echium",gsub("Eleagni","Eleagnus",gsub("Evonymus","Euonmyus",gsub("Fagi","Fagus",gsub("Fici","Ficus",gsub("Fananculi","Ranunculus",gsub("Fraxini","Fraxinus",Mycopo.no.match.n$all.hits.sp.Mycop.plant)))))))))))))))))))))))))))))))))))))))))))))))
Mycopo.no.match.n$match.cor<-gsub("Galii","Galium",gsub("Gei","Geum",gsub("Gelsemii","Gelsemium",gsub("Geranii","Geranium",gsub("Gerantii","Geranium",gsub("Heraclei","Heracleum",gsub("Helianthi","Helianthus",gsub("Hordei","Hordeum",gsub("Humuli","Humulus",gsub("Hyptidus","Hyptis", gsub("Ilicis","Ilex",gsub("Impatientis","Impatiens",gsub("Iridis","Iris",gsub("Juglandis","Juglans",gsub("Junci","Juncus",gsub("Juniperi","Juniperus", gsub("Lappae","Lappa",gsub("Lathyri","Lathyrus",gsub("Lauri","Laurus",gsub("Leontodontis","Leontodon",gsub("Ligustici","Ligusticum",gsub("Lupini","Lupinus",gsub("Lycii","Lycium",gsub("Mori","Morus",gsub("Medicagine","Medicago",gsub("Nasturtii","Nasturtium",gsub("Nerii","Nerium",gsub("Olacis","Olax",gsub("Olerodendron","Clerodendron",gsub("Ononidis","Ononis",gsub("Orobi","Orobus",gsub("Osmorrhizae","Osmorhiza",gsub("Panici","Panicum",gsub("Pentstamon","Penstemon",gsub("Pentstemonis","Penstemon",gsub("Peucedani","Peucedanum",gsub("Pini","Pinus",gsub("Plantaginis","Plantago",gsub("Polygoni","Polygonum",gsub("Polypodium","Polypodii",gsub("Populi","Populus",gsub("Pruni","Prunus",gsub("Pteridis","Pteridium",gsub("Pyri","Pyrus",gsub("Rhamni","Rhamnus",gsub("Rhois","Rhus",gsub("Rhoidis","Rhus",gsub("Ribis","Ribes",gsub("Roripa","Rorippa",Mycopo.no.match.n$match.cor)))))))))))))))))))))))))))))))))))))))))))))))))
Mycopo.no.match.n$match.cor<-gsub("Rubi","Rubus",gsub("Rumicis","Rumex",gsub("Salicis","Salex",gsub("Sambuci","Sambucus",gsub("Scirpi","Scirpus", gsub("Senecionis","Senecio",gsub("Silenis","Silene",gsub("Silphinium","Silphium",gsub("Strobilis","Strobilus",gsub("Sisymbrii","Sisymbrium",gsub("Smilacis","Smilax",gsub("Solani","Solanum",gsub("Sonchi","Sonchus",gsub("Symphyti","Symphytum",gsub("Tami","Tamus",gsub("Tanaceti","Tanacetum",gsub("Teucrii","Teucrium",gsub("Thalictri","Thalictrum",gsub("Thesii","Thesium",gsub("Trifolii","Trifolium",gsub("Tritici","Triticum",gsub("Ulmi","Ulmus",gsub("Vaccinii","Vaccinium",gsub("Veronicae","Veronica",gsub("Viburni","Viburnum", gsub("Zizyphus","Ziziphus",gsub("Zeae","Zea",Mycopo.no.match.n$match.cor)))))))))))))))))))))))))))

Mycopo.no.match<-Mycopo.no.match.n[which(!Mycopo.no.match.n$match.cor %in% Mycopo.no.match.n$all.hits.sp.Mycop.plant),]

#TPL in Taxonstand stopped working so I had to manually change what this did not catch. 
Mycopo.nM.TPL<-gnr_resolve(Mycopo.no.match$match.cor, preferred_data_sources = 167)
write.csv(Mycopo.nM.TPL, "intermediate_rds_csv/Mycoportal/Mycopo.nM.TPL.csv")
#I kept one copy of each, removed bad ones, and manually checked others.
Filt.Mycopo.nM.TPL<-read.csv("intermediate_rds_csv/Mycoportal/FILT.Mycopo.nM.TPL.csv")
Filt.Mycopo.nM.TPL<-Filt.Mycopo.nM.TPL[,2:6]
#only keep the names, not authorities.
Filt.Mycopo.nM.TPL$new.matched.names<-word(Filt.Mycopo.nM.TPL$matched_name, 1,2, sep=" ")
Mycopo.no.match$match.number<-match(Mycopo.no.match$match.cor, Filt.Mycopo.nM.TPL$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
Mycopo.no.match$match.name<-print(Filt.Mycopo.nM.TPL$new.matched.names[Mycopo.no.match$match.number[1:475]])
Filt.Mycopo.no.match<-Mycopo.no.match[which(!is.na(Mycopo.no.match$match.name)),]
Mycopo.wdraw2<-rbind.data.frame(Mycopo.wdraw, Filt.Mycopo.no.match[,1:27])

####last minute removal of on stem/stengeln/stammen
#Find tree genera entries
Mycopo.wdraw2.t<-Mycopo.wdraw2[which(as.character(Mycopo.wdraw2$plant.genus) %in% globtree$genus),]
M.h.t<-Mycopo.wdraw2.t[grep("^stem| stem|^stammen|stammen |^stengeln|stengeln ", Mycopo.wdraw2.t$habitat, ignore.case=TRUE),]
M.at.t<-Mycopo.wdraw2.t[grep("^stem| stem|^stammen|stammen |^stengeln|stengeln ", Mycopo.wdraw2.t$associatedTaxa, ignore.case=TRUE),]
M.all.t<-rbind.data.frame(M.h.t, M.at.t)
#remove these from dataframe
Mycopo.wdraw3<-anti_join(Mycopo.wdraw2, M.all.t)
keep.M.h.all.t<-M.all.t[grep("^lea| lea|petiol|Blattern|gall|canker|fruit|broom", M.all.t$habitat, ignore.case=TRUE),]
keep.M.at.all.t<-M.all.t[grep("^lea| lea|petiol|Blattern|gall|canker|fruit|broom", M.all.t$associatedTaxa, ignore.case=TRUE),]
Mycopo.workingdata<-rbind.data.frame(Mycopo.wdraw3, keep.M.h.all.t, keep.M.at.all.t)

saveRDS(Mycopo.workingdata, "intermediate_rds_csv/Mycoportal/Mycopo.workingdata.2025.05.04.rds") 
write.csv(Mycopo.workingdata, "OSF/Mycoportal/Mycopo.workingdata.2025.05.04.csv")

#keep only fungi IDd to sp.
Mycopo.workingdata.sp<-Mycopo.workingdata[which(!lengths(strsplit(as.character(Mycopo.workingdata$scientificName), ' ')) ==1),]


saveRDS(Mycopo.workingdata.sp, "intermediate_rds_csv/Mycoportal/Mycopo.workingdata.sp.2025.05.04.rds") 
write.csv(Mycopo.workingdata.sp, "OSF/Mycoportal/Mycopo.workingdata.sp.2025.05.04.csv")

#~# UPDATED FILE 4.2024 - Mycoportal data org

### I ORIGINALLY DID NOT INCLUDE FUNGAL AUTH INFORMATION, SO THIS CAN BE RETRIEVED FROM THE RAW DATA.

Myco.po <-read.csv("OSF/Mycoportal/Sorted.curled.Mycoportal.data.25May22.csv")
Myco.po2<- read.csv("OSF/Mycoportal/Sorted.handdownloaded.Mycoportal.data.15May22.csv")

#Do the exact same subsetting I did as before (in Database.Ib.Mycoportal.2023.11.16.R), except adding in scientificNameAuthorship

Myco.po1.ab <- data.frame(id=Myco.po$id, institutionCode = Myco.po$institutionCode, collectionCode = Myco.po$collectionCode, occurrenceID = Myco.po$occurrenceID, basisOfRecord = Myco.po$basisOfRecord,scientificName = Myco.po$scientificName, genus = Myco.po$genus, taxonID = Myco.po$taxonID, taxonRank = Myco.po$taxonRank, occurrenceRemarks = Myco.po$occurrenceRemarks, habitat = Myco.po$habitat, associatedTaxa = Myco.po$associatedTaxa, substrate = rep("NA"), recordId = Myco.po$recordId, country = Myco.po$country, db=rep("Mycoportal.db"), scientificNameAuthorship = Myco.po$scientificNameAuthorship)
# the handdownloaded data do not have a substrate column, for whatever reason.
Myco.po2.ab <- data.frame(id=Myco.po2$id, institutionCode = Myco.po2$institutionCode, collectionCode = Myco.po2$collectionCode, occurrenceID = Myco.po2$occurrenceID, basisOfRecord = Myco.po2$basisOfRecord, scientificName = Myco.po2$scientificName, genus = Myco.po2$genus, taxonID = Myco.po2$taxonID, taxonRank = Myco.po2$taxonRank, occurrenceRemarks = Myco.po2$occurrenceRemarks, habitat = Myco.po2$habitat, associatedTaxa = Myco.po2$associatedTaxa, substrate = Myco.po2$substrate, recordId = Myco.po2$recordId, country = Myco.po2$country, db=rep("Mycoportal.db"), scientificNameAuthorship = Myco.po2$scientificNameAuthorship)

Myco.po.ab<-rbind(Myco.po1.ab, Myco.po2.ab)
#saveRDS(Myco.po.ab, 'OSF/Mycoportal/Myco.po.ab.with.fungiauth.2024.rds')

#match to author names from raw data file.
Mycopo.workingdata.sp$Fungi.auth <- Myco.po.ab$scientificNameAuthorship[match(Mycopo.workingdata.sp$id, Myco.po.ab$id)]


###### added further steps in 2025-01-05 in order to standardize the plant taxa consistently, remove nonplant taxa

#filter out nonplant hosts - ncbi classification with new match.name2 column (taxizedb package).
Mycopo.workingdata.sp$match.name.sp<-word(Mycopo.workingdata.sp$match.name, 1, 2)
my.h<-unique(Mycopo.workingdata.sp$match.name.sp)
class2<-classification(my.h, db="ncbi")

np<-my.h[setdiff(grep("cellular organisms", class2), grep("Viridiplantae", class2))]
#get non-plant genera

Mwd1<-Mycopo.workingdata.sp[which(!Mycopo.workingdata.sp$match.name.sp %in% np),]
my.h2<-unique(Mwd1$match.name.sp)
corr.plant<-nameMatch_WCVP(my.h2)

#these all have accepted_SPNAMEs, so skipping filtration.

Mwd2<-Mwd1[which(Mwd1$match.name.sp %in% corr.plant$Submitted_Name),]
Mwd2$Genus_in_database<-corr.plant$Genus_in_database[match(Mwd2$match.name.sp, corr.plant$Submitted_Name)]
Mwd2$Name_in_database<-corr.plant$Name_in_database[match(Mwd2$match.name.sp, corr.plant$Submitted_Name)]
Mwd2$Author_in_database<-corr.plant$Author_in_database[match(Mwd2$match.name.sp, corr.plant$Submitted_Name)]
Mwd2$New_Author<-corr.plant$New_author[match(Mwd2$match.name.sp, corr.plant$Submitted_Name)]
Mwd2$Accepted_SPNAME<-corr.plant$Accepted_SPNAME[match(Mwd2$match.name.sp, corr.plant$Submitted_Name)]

Mwd2$Accepted_SPAUTHOR<-Mwd2$Author_in_database
Mwd2$Accepted_SPAUTHOR[which(!is.na(Mwd2$New_author))]<-Mwd2$New_author[which(!is.na(Mwd2$New_author))]

saveRDS(Mwd2, "intermediate_rds_csv/Mycoportal/2025.01.Mwd2.rds")

Mwd3<-Mwd2 %>%
  group_by(Accepted_SPNAME, scientificName) %>%
  distinct(Accepted_SPNAME, .keep_all=TRUE)
#
corrected.Host.sp<-Mwd3$Accepted_SPNAME
corrected.Host.ssp<-word(Mwd3$Accepted_SPNAME, 3, -1, sep=" ")
corrected.Host.auth<-Mwd3$Accepted_SPAUTHOR
detected.Host.plant<-Mwd3$all.hits.sp.Mycop.plant
accepted.Fungi.sp<-Mwd3$scientificName
accepted.Fungi.sp.full.name<-paste(Mwd3$scientificName, Mwd3$Fungi.auth)
verbatim.Fungi.sp<-Mwd3$scientificName
habitat.and.plant.part.1<-Mwd3$habitat
habitat.and.plant.part.2<-Mwd3$associatedTaxa
habitat.and.plant.part.3<-Mwd3$substrate
DB<-"Mycoportal"
Location<-Mwd3$country
Ref<-Mwd3$occurrenceID
accepted.Fungi.auth<-Mwd3$Fungi.auth

#

#make final dataframe with both corrected and uncorrected host names in separate columns.
Myco.po.FIN.df<-cbind.data.frame(corrected.Host.sp, corrected.Host.ssp, corrected.Host.auth, 
                                detected.Host.plant, accepted.Fungi.sp, accepted.Fungi.sp.full.name, accepted.Fungi.auth, verbatim.Fungi.sp, 
                                habitat.and.plant.part.1, habitat.and.plant.part.2, habitat.and.plant.part.3, Location, Ref, DB)

Myco.po.FIN.df<-Myco.po.FIN.df[!Myco.po.FIN.df$accepted.Fungi.sp=="",]
saveRDS(Myco.po.FIN.df, "Cleaned_data/2025-05.UPD.Myco.po.FIN.df.rds")
