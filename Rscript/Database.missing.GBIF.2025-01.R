setwd("~/OneDrive - UBC/DB_fung_2025.01/")

# there were some records that need to be dealt with in the GBIF file that have text in 2 fields (associatedTaxa, habitat, associatedOrganisms)
# These records are branched and handled in this file for the Fungi-on-plant and On Living Tissue databases.

GBIF.hits.sp<-readRDS("OSF/GBIF/GBIF.hits.sp.new.rds")
GB.hab.na<- GBIF.hits.sp %>% filter(is.na(associatedOrganisms) & is.na(associatedTaxa))
GB.asOr.na<- GBIF.hits.sp %>% filter(is.na(habitat) & is.na(associatedTaxa))
GB.asTa.na<- GBIF.hits.sp %>% filter(is.na(associatedOrganisms) & is.na(habitat))
GBIF.not.missing<-rbind.data.frame(GB.hab.na, GB.asOr.na, GB.asTa.na)


GBIF.hits.missing<-anti_join(GBIF.hits.sp, GBIF.not.missing, by = "gbifID")
GBIF.hits.missing <- GBIF.hits.missing %>%
  mutate(row_number = row_number())

# The following naming convention and save directory will be used:
saveRDS(GBIF.hits.missing, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2024-12-21.missing.GBIF.hits.R")


#~#~#~## Use taxon adjustment for hosts that I already determined where possible.
#revise name list using manual corrections I already made.
GBIF.hits.missing$genera<-word(GBIF.hits.missing$GBIF.plant, 1)

# Correct genus names. If Genus names are corrected, the TPL and other tools perform better. These are primarily modified latin names, but some include easy-to-spot typos. All apparent genera typos with >1 entry were corrected. 
GBIF.hits.missing$sing.genus<-gsub("Populi","Populus",gsub("Abietis","Abies",gsub("Aceris","Acer",
                                                                               gsub("Aconiti","Aconitum",(gsub("Agropyrum","Agropyron",gsub("Agrostidis","Agrostis", gsub("Allii","Allium", gsub("Alni","Alnus",
                                                                                                                                                                                                 gsub("Baccharidis", "Baccharis", gsub("Bromi","Bromus",gsub("Calamogrostidis","Calamogrostis",gsub("Caricis","Carex", gsub("Carpini","Carpinum",
                                                                                                                                                                                                                                                                                                                            gsub("Clerodendron","Clerodendrum", gsub("Corni","Cornus",gsub("Crataegi","Crataegus", gsub("Cytisi", "Cytisum", gsub("Delphini","Delphinium", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                  gsub("Desmodii","Desmodium",gsub("Dipsaci", "Dipsacus",gsub("Elymi","Elymus",gsub("Euonymi","Euonymus", gsub("Eupatorii","Eupatorium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               gsub("Evonymi","Euonymus",gsub("Evonymus","Euonmyus",gsub("Fagi","Fagus",gsub("Fici","Ficus",gsub("Fraxini","Fraxinus",gsub("Galii","Galium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           gsub("Gei","Geum",gsub("Gelsemii","Gelsemium",gsub("Geranii","Gernaium",gsub("Heraclei","Heracleum",gsub("Hordei","Hordeum",gsub("Humuli","Humulus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            gsub("Ilicis","Ilex",gsub("Impatientis","Impatiens",gsub("Iridis","Iris",gsub("Juglandis","Juglans",gsub("Junci","Juncus",gsub("Juniperi","Juniperus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           gsub("Lappae","Lappa",gsub("Lathyri","Lathyrus",gsub("Lauri","Laurus",gsub("Leontodontis","Leontodon",gsub("Ligustici","Ligusticum",gsub("Lupini","Lupinus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    GBIF.hits.missing$genera))))))))))))))))))))))))))))))))))))))))))))))))
GBIF.hits.missing$sing.genus<-gsub("Lycii","Lycium",gsub("Mori","Morus",gsub("Nasturtii","Nasturtium",gsub("Nerii","Nerium",gsub("Olerodendron","Clerodendron",
                                                                                                                              gsub("Ononidis","Ononis",gsub("Orobi","Orobus",gsub("Osmorrhizae","Osmorhiza",gsub("Panici","Panicum",gsub("Pentstamon","Penstemon",gsub("Pentstemonis","Penstemon",
                                                                                                                                                                                                                                                                       gsub("Peucedani","Peucedanum",gsub("Pini","Pinus",gsub("Plantaginis","Plantago",gsub("Polygoni","Polygonum",gsub("Polypodium","Polypodii",gsub("Populi","Populus",
                                                                                                                                                                                                                                                                                                                                                                                                                      gsub("Pruni","Prunus",gsub("Pteridis","Pteridium",gsub("Pyri","Pyrus",gsub("Rhamni","Rhamnus",gsub("Rhois","Rhus",gsub("Rhoidis","Rhus",gsub("Ribis","Ribes",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   gsub("Roripa","Rorippa",gsub("Rubi","Rubus",gsub("Rumicis","Rumex",gsub("Salicis","Salex",gsub("Sambuci","Sambucus",gsub("Scirpi","Scirpus",gsub("Senecionis","Senecio",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    gsub("Silphinium","Silphium",gsub("Sisymbrii","Sisymbrium",gsub("Smilacis","Smilax",gsub("Solani","Solanum",gsub("Sonchi","Sonchus",gsub("Symphyti","Symphytum",gsub("Tami","Tamus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         gsub("Tanaceti","Tanacetum",gsub("Teucrii","Teucrium",gsub("Thalictri","Thalictrum",gsub("Thesii","Theisum",gsub("Trifolii","Trifolium",gsub("Tritici","Triticum",gsub("Ulmi","Ulmus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                gsub("Vaccinii","Vaccinium",gsub("Veronicae","Veronica",gsub("Viburni","Viburnum", gsub("Zizyphus","Ziziphus",GBIF.hits.missing$sing.genus)))))))))))))))))))))))))))))))))))))))))))))))))

#make new column we can use for Taxonstand comparison
interm.GBIF.plant<-paste(GBIF.hits.missing$sing.genus, word(GBIF.hits.missing$GBIF.plant, 2, -1))


nf.for.Taxonstand<-as.character(unique(interm.GBIF.plant))
#read in Taxonstand conversion sheet from last time (this won't be comprehensive, but a start)
nf.TPL.taxons<-readRDS("intermediate_rds_csv/GBIF/nf.TPL.taxons.rds")
#only keep accepted changes.
nf.TPL.taxons2<-nf.TPL.taxons[!(nf.TPL.taxons$New.Taxonomic.status ==""),]
#Read in the rest of the key
re.TPL.G<-readRDS("intermediate_rds_csv/GBIF/re.TPL.G.rds")
re.TPL.G2<-re.TPL.G[!(re.TPL.G$New.Taxonomic.status ==""),]#
#combine the two keys.
taxon.key<-unique(rbind.data.frame(nf.TPL.taxons2, re.TPL.G2))

#matches
matches.TPL.taxons<-match(nf.for.Taxonstand, taxon.key$Taxon)

GBIF.newname<-paste(taxon.key$New.Genus, taxon.key$New.Hybrid.Marker, taxon.key$New.Species,taxon.key$New.Authority)

#make conversion chart with original data
Ghm.conv<-data.frame(nf.for.Taxonstand)
Ghm.conv$corrected<-NA
Ghm.conv$corrected<-GBIF.newname[matches.TPL.taxons]
gsub("  ", " ", Ghm.conv$corrected)
Ghm.conv$corrected<-gsub("  ", " ", Ghm.conv$corrected)


#replace the corrected ones
matches.T<-match(GBIF.hits.missing$GBIF.plant, Ghm.conv$nf.for.Taxonstand)

#finally, make new column on working data dataframe, keep the normal ones
GBIF.hits.missing$corr.GBIF.plant<-Ghm.conv$corrected[matches.T]

GBIF.hits.missing2 <- GBIF.hits.missing %>%
  mutate(corr.GBIF.plant = if_else(is.na(corr.GBIF.plant), GBIF.plant, corr.GBIF.plant))



saveRDS(GBIF.hits.missing2, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2024-12-21.missing.GBIF.hits.missing2.rds")


#throw out the plant/fungi pairs already kept GBIF:

#NOTE: branches from 2024 R script here.
GB2<-readRDS("intermediate_rds_csv/GBIF/GB2.rds")
GBIF.hits.missing.minus.working<-anti_join(GBIF.hits.missing2, GB2, by = c("species","GBIF.plant"))

###
#For now, we will mostly follow the exact same procedure as before, starting after point in database.1c.GBIF.2023.11.R when all asOr, asTa, and ashab have been combined:


#filter plants with duplicated first letters. This gets rid of a few entries where there are 2 distinct genera in the host field as well, but those are typically descriptions of the area.
GBIF.hits.missing.minus.working$first.letters<-substr(GBIF.hits.missing.minus.working$GBIF.plant, 1, 6)
GBIF.all.fields.filt<-GBIF.hits.missing.minus.working %>% 
  group_by(gbifID, first.letters) %>%
  arrange(desc(nchar(as.character(GBIF.plant)))) %>%
  distinct(first.letters, .keep_all=TRUE)


saveRDS(GBIF.all.fields.filt, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.all.fields.filt.rds")


GBIF.on.plant<-GBIF.all.fields.filt %>%
  filter_at(vars(associatedTaxa,habitat,associatedOrganisms),any_vars(grepl(
    'host|infect| on |^on |"on |[on] |0n |parasit|leaves|leaf|foglie|twig|stem|huesped| sur |^sur |"sur |sulle |sulla | auf |^auf |"auf | ad |^ad |"ad | ^in |in |"in |^an |"an | an |paa |lebenden|anther| bark|ritidoma|branch|
      spot |spots |splotch|rust|streak|causing|inflorescence|fruit |fruits |pathogen|petiole|pod |pods |root |roots|cone|cones|utrículo|utricle|covering |culm|ovar|foli|needle|acícula|caulibus|caule|cortic|in den |flower|
      seed|pe frunze de|pe ramuri de|pe tulpini de|rama|smut |sobre |super |supra | sub |^sub |blättern|blattern|blad |stängeln|kapseln|stam |stamme |trunk|tronc|troco|truncus|living|live |vivant|sôbre|
      Wurzeln|Fruchtstielen|Schoten|Halmen|Rinde | vine|scape| head|stalk|frond|kvist|Balttsticlen|wedel|ramis|caudex|caudice|viva|rama|Epífita|hospedero|encontrado sobre|stolpe |
      Forófito de|onAbies|endophyt|^on:|substr|^på | på |^на | на | bladeren | Op |^Op |^aan | aan |^nos |sobre|substrat|sustrato| ^en |en |"en |bark | коре| ствол| лист| ветв| живы| корнях | корен | стеблях|
     плод| хвое | хвоя | завязи|nährpflanze|on:', (.), ignore.case = TRUE)))

####~#~#~#~ DIFFERENCE FROM other GBIF 
## NOTE: here, we added "on:"
####~#~#~#~ 

saveRDS(GBIF.on.plant, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.on.plant.rds")

### Optional
#Go through those that may have good info and revise previous command as needed.
GBIF.not.on.plant<-setdiff(GBIF.all.fields.filt,GBIF.on.plant)
saveRDS(GBIF.not.on.plant, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.not.on.plant.rds")
#use Google translate

GBIF.remove.on.plant<-GBIF.on.plant[apply(GBIF.on.plant[, c("habitat", "associatedOrganisms", "associatedTaxa")], 1, 
                                         function(row) any(grepl("soil|humus|suelo|litter|under|duff|ground|dung| on rock| on granit|on moss|terricolous| sur sol |on limestone|on sandstone|on stone|in moss|host: along the trail|
                                                                 among moss| along the dried|in dried-up|drying mud|along the creek|along the forest road| open rock outcrop|open mossy rock outcrop|on forest floor|among moss|on grassy", ignore.case=TRUE, row))), ]

# Note: There are likely usable host-plant pairs in GBIF.remove.on.plant that can be salvaged for the "Fungi-on-plant" db, but the data is very messy and many are likely documented elsewhere.
GBIF.m.on.plant.FIN<-data.frame(setdiff(GBIF.on.plant, GBIF.remove.on.plant))


#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# II. Here we create GBIF.m.on.plant.fin, which will be used as raw data for the "On plant tissue" database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

saveRDS(GBIF.m.on.plant.FIN, 'intermediate_rds_csv/GBIF/GBIF_multiple_fields/GBIF.m.on.plant.FIN.rds')

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


GBIF.deadplant <- GBIF.on.plant[apply(GBIF.on.plant[, c("habitat", "associatedOrganisms", "associatedTaxa")], 1, 
                                     function(row) any(grepl("soil|humus|suelo|litter|under|dead|stump|fallen|^log| log|dying|muert|rott|decay|duff|ground|emortu|
                                                             base of tree|hypogeous|dung|putr|død|charred|burned|burnt|sticks|voet|wood|ved|tocon|stubb|lignicolous|snag|
                                                             madera| on rock| on granit|on moss|terricolous| sur sol |base of old |on limestone|on sandstone|on stone|
                                                             base of trunk|in moss|host: along the trail|among moss| along the dried|in dried-up|drying mud|along the creek|
                                                             along the forest road| open rock outcrop|open mossy rock outcrop|on forest floor|among moss|on grassy", ignore.case=TRUE, row))), ]


saveRDS(GBIF.deadplant, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.deadplant.2024-12-21.rds")


###

#### Here are living plants
GBIF.liveplant<-anti_join(GBIF.on.plant, GBIF.deadplant, by = "gbifID")
saveRDS(GBIF.liveplant, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.liveplant.rds")


### There are too many live plant entries here. Edit some out. Specifically, those already found from Mycoportal.
myco.edit.liveplant<-read.csv("OSF/Mycoportal/Mycop.edit.liveplant.csv")
# generate list
myco.36<-myco.edit.liveplant[which(nchar(as.character(myco.edit.liveplant$occurrenceID))==36),]
GBIF.myc.filter<- GBIF.liveplant[which(!GBIF.liveplant$occurrenceID %in% myco.36$occurrenceID),]
saveRDS(GBIF.myc.filter, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.myc.filter.rds")



#### Filter out fungi (Index Fungorum source)

plnamresolv.F.all<-readRDS("intermediate_rds_csv/GBIF/plnamresolv.Fungi.all.rds")
GBIF.myc.filter$match.number<-match(as.character(GBIF.myc.filter$GBIF.plant), plnamresolv.F.all$user_supplied_name, nomatch = NA_integer_, incomparables = NULL)
GBIF.myc.filter$match.name<-print(plnamresolv.F.all$matched_name[GBIF.myc.filter$match.number[1:43365]])

#dataframe without fungi
no.fung<-GBIF.myc.filter[is.na(GBIF.myc.filter$match.name),]
#dataframe of fungi
fung<-GBIF.myc.filter[!is.na(GBIF.myc.filter$match.name),]

#reset match name and match number columns for TPL step.
no.fung<-no.fung[,1:23]

#remove lichens
#read in lichenlist data from "Consortium of Lichen Herbaria" https://lichenportal.org/portal/checklists/checklist.php?clid=1492&pid=558
lichenlist<-read.csv("ext_data/World Checklist of Genera of Lichenized Fungi_1698006802.csv")
#make genus column
lichenlist$genus<-gsub(" sp.", "", lichenlist$ScientificName)
#filter lichenizing fungi
hfilt.GBIF.nl<-no.fung[which(!no.fung$genus %in% lichenlist$genus),]
#remove fungal taxon entries that contain word root "lichen"
nlaf<-hfilt.GBIF.nl[grep("lichen",hfilt.GBIF.nl$scientificName, ignore.case=TRUE),]
hfilt.GBIF.nl2<-hfilt.GBIF.nl[which(!hfilt.GBIF.nl$scientificName %in% nlaf$scientificName),]


# Read in all taxonomically adjusted GBIF data
#NOTE: Difference from earlier draft.
#Here, we will standardize plant host names where possible using WCVP.
#library(U.Taxonstand)
corr.plantGwdr<-nameMatch_WCVP(hfilt.GBIF.nl2$corr.GBIF.plant)
saveRDS(corr.plantGwdr, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.corr.plantGwdr.rds")

hfilt.GBIF.nl2<-data.frame(hfilt.GBIF.nl2)
hfilt.GBIF.nl2$Name_in_database<-corr.plantGwdr$Name_in_database
hfilt.GBIF.nl2$Author_in_database<-corr.plantGwdr$Author_in_database
hfilt.GBIF.nl2$New_name<-corr.plantGwdr$New_name
hfilt.GBIF.nl2$New_author<-corr.plantGwdr$New_author
hfilt.GBIF.nl2$Accepted_SPNAME<-corr.plantGwdr$Accepted_SPNAME

#taxizedb
#look up the ones for which there was no correction in original data.
nonplant.check<-unique(corr.plantGwdr[is.na(corr.plantGwdr$Accepted_SPNAME),])
classif.nonplant.check<-classification(nonplant.check$Submitted_Name, db="ncbi")

# Which are definitely not plants?
np.results<-nonplant.check[setdiff(grep("cellular organisms", classif.nonplant.check), grep("Viridiplantae", classif.nonplant.check)),]
#get non-plant genera
npr.gen<-unique(np.results$Submitted_Genus)

#remove them
GB.rem<-which(word(hfilt.GBIF.nl2$corr.GBIF.plant, 1) %in% npr.gen)
hfilt.GBIF.nl3<-hfilt.GBIF.nl2[-GB.rem,]
saveRDS(hfilt.GBIF.nl3, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/hfilt.GBIF.nl3.rds")


### END of OLD


###~#~#~#~~#~#~#~#~#~#~#~#~###~#~#~#~#~#~#~#~#~#
### Filtering to living tissue (non-conservative)
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##~#~#~#~#~#~#
#here we can start filtering out dead tissue.
#there really aren't any associatedOrganisms fields left
GBIF.deadtissue <-hfilt.GBIF.nl3[apply(hfilt.GBIF.nl3[, c("habitat", "associatedTaxa")], 1, 
                                      function(row) any(grepl("bark|branch| on old |remains|branchlet| daed| remnants|tronc|on debris|in debris|
                                      in leaf debris|in moist debris|on accumulated debris|on remain of|on trunk of|
                                      on grassy|on last year's|on a piece of|fire-killed|base of tree|auf durren|decorticated|
                                      on moist mud|on drying mud|in drying bed|in last years'|on a compost pile|
                                      on an old|on base of tree|on base of living tree|
                                      on uredo |rhizosphere|fol. mort.", ignore.case=TRUE, row))), ]

GBIF.livetissue.rough<-anti_join(hfilt.GBIF.nl3, GBIF.deadtissue, by = "gbifID")

#remove lichens as host remaining in fields
GBIF.livetissue.rough.NL<-GBIF.livetissue.rough[which(!GBIF.livetissue.rough$genera %in% lichenlist$genus),]
GBIF.lichens.rem <- GBIF.livetissue.rough.NL[unlist(sapply(lichenlist$genus, function(x) grep(paste("Host of:", x), GBIF.livetissue.rough.NL$associatedTaxa))),]
GBIF.lichens.rem2 <- GBIF.livetissue.rough.NL[unlist(sapply(lichenlist$genus, function(x) grep(paste("on", x), GBIF.livetissue.rough.NL$habitat))),]
GBIF.lichens.rem3 <- GBIF.livetissue.rough.NL[unlist(sapply(lichenlist$genus, function(x) grep(paste("On", x), GBIF.livetissue.rough.NL$habitat))),]
GBIF.lichens.rem4 <- GBIF.livetissue.rough.NL[unlist(sapply(lichenlist$genus, function(x) grep(paste("on", x), GBIF.livetissue.rough.NL$associatedTaxa))),]
GBIF.lichens.rem5 <- GBIF.livetissue.rough.NL[unlist(sapply(lichenlist$genus, function(x) grep(paste("On", x), GBIF.livetissue.rough.NL$associatedTaxa))),]
GBIF.rem.lichens.all<-rbind.data.frame(GBIF.lichens.rem, GBIF.lichens.rem2, GBIF.lichens.rem3, GBIF.lichens.rem4, GBIF.lichens.rem5)
GBIF.lt.rough.nl<-anti_join(GBIF.livetissue.rough.NL, GBIF.rem.lichens.all, by = "gbifID")

saveRDS(GBIF.lt.rough.nl, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/GBIF.lt.rough.nl.rds")
#set aside surely living tissue

#### Almost final stage.
### a-e blocks are kept and set aside.

######################
##### a. #############
######################
GBIF.livetissue1 <- GBIF.lt.rough.nl[apply(GBIF.lt.rough.nl[, c("habitat", "associatedTaxa")], 1, 
                                      function(row) any(grepl("living leaves|n leaves of|n fertile|follis vivis|n living|n live|n catkins|n ovariis|
                                      n culms|auf lebenden|n foliis|fol.|n the leaves of|n living leaves of|n flowers of|n seeds of|n fruits of|n inflorescence of|n culms of|n living twigs of|
                                      n living stems of|n ovaries of|ad folia|An Blättern von|hojas de|n leaves and stems of|n leaves and petioles of|n perigyni|
                                      den Blättern von|utrículos de|Sur les feuilles|n anthers|n the anthers|rust|uredospore| witches|På blad av|ovarios de|
                                      levende stam van|Op blad|n the lower side of leaves of|n stems of|n stems and leaves of|n stem of|n stalk of|n stalks of|n spikelet of|
                                      n seedlings of|n petioles of|n male catkins|n male spikes|n live leaves|n living and wilt|n leaves, sheaths|n leaves, stalks|n leaves, stems|
                                      n leaves, petioles|n leaves & stems of|n leaves and culms of|n leaves and bracts of|n fruit of|n fronds of|n female catkins of|culms and leaves of|
                                      n canes of|on attached needles|on attached leaves|Substrat:Blad|folia viva|in foliis|in utricles of|in floribus|in antheris|^leaves of|leaf spot|
                                      leaf lesion|leaf and stem gall|inflorescencia|in spikelets|in petiolis|in ovariis|in internodes|in internodiis|in fol.|In den Bluten von|
                                      In den Antheren von|flores de |espor hojas d|espigas de|n hoja de|n frutas de|n leaves and bracts of|Auf stengeln von|lebenden Blättern|
                                      Auf Asten von|anteras de|ad caules|Auf Blättern|På blad av|en tallo de|en tallos de", ignore.case=TRUE, row))), ]

#new running rough dataset
GBIF.livetissue.rough2<-anti_join(GBIF.livetissue.rough, GBIF.livetissue1, by = "gbifID")

#now let's look for entries in which the Host of: and host: do not match the flagged plant species (GBIF.plant). This means that these are non-plant hosts
#ID rows with this format

#subset dataframe to where it reads Host of: then <corresponding host plant>
nonpl.rem <- GBIF.livetissue.rough2 %>%
  rowwise() %>%
  filter(grepl(paste0(".*Host of:.*", GBIF.plant, ".*"), associatedTaxa, ignore.case=TRUE)) %>%
  ungroup ()

#differences in subset dataframe to dataframe with all "Host of" entries.
GBIF.nonpl.rem<-anti_join(GBIF.livetissue.rough2[grep("Host of:", GBIF.livetissue.rough2$associatedTaxa, , ignore.case=TRUE),], nonpl.rem, by="gbifID" )

#Repeat for "Host:
nonpl.rem2 <- GBIF.livetissue.rough2 %>%
  rowwise() %>%
  filter(grepl(paste0(".*Host:.*", GBIF.plant, ".*"), associatedTaxa, ignore.case=TRUE)) %>%
  ungroup ()
GBIF.nonpl.rem2<-anti_join(GBIF.livetissue.rough2[grep("Host:", GBIF.livetissue.rough2$associatedTaxa, ignore.case=TRUE),], nonpl.rem2, by="gbifID" )


GBIF.nonpl.rem.all<-rbind.data.frame(GBIF.nonpl.rem, GBIF.nonpl.rem2)

GBIF.pl<-anti_join(GBIF.livetissue.rough2, GBIF.nonpl.rem.all, by="gbifID")

saveRDS(GBIF.pl, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.pl.rds")

#I retrieved the csv for the GlobalTreeSearch data on Oct 31 2023. https://rpubs.com/Roeland-KINDT/812716
#read this data in 
globtree<-read.csv("ext_data/global_tree_search_trees_1_7.csv")
#filter by genus in order to eliminate as many shrubby genera, as well USING OLD NAMES
gltrgenus<-data.frame(str_split_fixed(globtree$TaxonName, " ", 2))
globtree$genus<-gltrgenus$X1
unk.genus<-data.frame(str_split_fixed(GBIF.pl$sing.genus, " ", 2))
GBIF.pl$plant.genus<-unk.genus$X1


#filter out tree data from the list
### NOTE: This time I am taking the trees out of the remaining data, not just fields that say "On <taxon name>"

# Manually filter this data, but keep "on <plant>"
GBIF.pl.NOTREE<-GBIF.pl[which(!GBIF.pl$plant.genus %in% globtree$genus),]
onGpl<-GBIF.pl.NOTREE[which(GBIF.pl.NOTREE$habitat %in% paste("on", GBIF.pl.NOTREE$GBIF.plant)),]
onGpl2<-GBIF.pl.NOTREE[which(GBIF.pl.NOTREE$habitat %in% paste("On", GBIF.pl.NOTREE$GBIF.plant)),]
onGpl3<-GBIF.pl.NOTREE[which(GBIF.pl.NOTREE$habitat %in% paste("On:", GBIF.pl.NOTREE$GBIF.plant)),]

######################
##### b. #############
######################
#KEEP these
GBIF.pl.NOTREE.ON<-rbind.data.frame(onGpl, onGpl2, onGpl3)



GBIF.pl.NOTREE.FIN<-anti_join(GBIF.pl.NOTREE, GBIF.pl.NOTREE.ON, by="gbifID")
GBIF.pl.NOTREE.FIN<-GBIF.pl.NOTREE.FIN[grep("Vegetación secundaria", GBIF.pl.NOTREE.FIN$GBIF.plant, invert=TRUE),]
write.csv(GBIF.pl.NOTREE.FIN,"intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.pl.NOTREE.FIN.csv" )

# Manually filter this data, but throw out "on <plant>" if no other parts noted.
GBIF.pl.TREE<-GBIF.pl[which(GBIF.pl$plant.genus %in% globtree$genus),]
onGplT<-GBIF.pl.TREE[which(GBIF.pl.TREE$habitat %in% paste("on", GBIF.pl.TREE$GBIF.plant)),]
onGplT2<-GBIF.pl.TREE[which(GBIF.pl.TREE$habitat %in% paste("On", GBIF.pl.TREE$GBIF.plant)),]
onGplT3<-GBIF.pl.TREE[which(GBIF.pl.TREE$habitat %in% paste("On:", GBIF.pl.TREE$GBIF.plant)),]
#DITCH these
GBIF.pl.TREE2<-anti_join(GBIF.pl.TREE,rbind.data.frame(onGplT, onGplT2, onGplT3), by="gbifID")
#does it mention plant part?
gt2.keep<-unique(c(grep("leaves|leaf|twig|Blattern|fitopatógeno|feuille|needle|bract|seed|fruit|samara|pod|aeci|Fruchten|Knollen|capsul|petiol|carpel|ament| bud|ramul|flower|anther|ovar", GBIF.pl.TREE2$habitat, ignore.case = TRUE),
                    grep("leaves|leaf|twig|Blattern|fitopatógeno|feuille|needle|bract|seed|fruit|samara|pod|aeci|Fruchten|Knollen|capsul|petiol|carpel|ament| bud|ramul|flower|anther|ovar", GBIF.pl.TREE2$associatedTaxa, ignore.case=TRUE)))
GBIF.pl.TREE.FIN<-GBIF.pl.TREE2[gt2.keep,]
write.csv(GBIF.pl.TREE.FIN, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.pl.TREE.FIN.csv")
#Manually edited: removed "host: <>" formatted ones, many were mycorrrhizal. Removed all that did not mention living tissue parts somewhere.

######################
##### c. #############
######################
#removed epiphytes and phorophytes, kept twigs and phytopathogens. Removed on cones (usually saprotrophs). Kept aecidial host.
GBIF.pl.TREE.FIN.filt<-read.csv("intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.pl.TREE.FIN.filt.csv")



##MANUAL edits:
#in the non-tree data, there were many plants that did not get detected correctly. 
# Manually removed all tree species (except if on leaves, or similar) and fungal hosts I spotted (not comprehensive), and cleaned up the badly formatted entries.
# I manually changed these and added these into the corr.GBIF.plant column (these can be re-assessed with taxize/U.Taxonstand).
#There were still several duplicate entries. Duplicates removed where possible, and non-plant entries. There are a few animals in here too, needs final clean.
G.notree.filt<-read.csv("intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.pl.NOTREE.FIN.filt.csv")

#filter the trees again
ug.f<-data.frame(str_split_fixed(G.notree.filt$corr.GBIF.plant, " ", 2))
G.notree.filt$plant.genus<-ug.f$X1
G.notree.filt2<-G.notree.filt[which(!G.notree.filt$plant.genus %in% globtree$genus),]
#ditch trees in data. 
xom<-setdiff(G.notree.filt, G.notree.filt2)

######################
##### d. #############
######################
# keep those mentioning leaves, etc.
G.kept.trees<-xom[c(grep("leaves|endofito|fitopatógeno|leaf|feuilles|fruit|capsul|pods|ram.", xom$associatedTaxa),grep("leaves|endofito|fitopatógeno|leaf|feuilles|fruit", xom$habitat)),]

#retry ncbi classification with new match.name2 column (taxizedb package).
cff2<-unique(G.notree.filt2$corr.GBIF.plant)
class2<-classification(cff2, db="ncbi")

np2<-cff2[setdiff(grep("cellular organisms", class2), grep("Viridiplantae", class2))]
#get non-plant genera
np.gen2<-word(np2, 1)

######################
##### e. #############
######################
G.notree.filt2.nf<-G.notree.filt2[which(!G.notree.filt2$plant.genus %in% np.gen2),]


#####

###~#~#~#~~#~#~#~#~#~#~#~#~###~#~#~#~#~#~#~#~#~#
### Filtering to living tissue (conservative+final filter)
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##~#~#~#~#~#~#

#Combine the set aside data (a-e) ("on live plant part", "not on tree on <plant name>", "not on tree, otherwise on plant part", "trees rescued from nontree data", and "on tree plant part" data).
GBIF.miss<-rbind.data.frame(GBIF.livetissue1, GBIF.pl.NOTREE.ON[,1:28], GBIF.pl.TREE.FIN.filt[,2:29], G.kept.trees[,2:29], G.notree.filt2.nf[,2:29])

saveRDS(GBIF.miss, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.miss.rds")

#manually change entries in corr.GBIF.plant field so that they get detected better by U.Taxonstand.
write.csv(GBIF.miss, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.miss.csv")

#Deleted duplicates with worse host name detection (A. macrophyllum, rather than Acer macrophyllum), 
# then fixed the corr.GBIF.plant column if they produced NAs or incorrect entries with previous U.Taxonstand round. 
# Removed errant entries that did not imply direct living tissue contact.
#Removing everything on dead plants or not on plants manually. I fixed the names in the corr.GBIF.plant column if they were not the plant specified as the host.
#I also changed the genera in the corr.GBIF.plant column to non-plural, where I was able to recognize it.

GBIF.miss.filt<-read.csv("intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.GBIF.miss.filt.csv")

#run another name match with U.Taxonstand to aid filtration, now that changes have been made to corr.GBIF.plant.
uni.corr.plant<-unique(GBIF.miss.filt$corr.GBIF.plant)
corr.plantGmf<-nameMatch_WCVP(uni.corr.plant)



#remove the ones that still have no match. these are virtually all nonplants.
corr.plant.GMf.filt<-corr.plantGmf[which(!is.na(corr.plantGmf$Accepted_SPNAME)),]
saveRDS(corr.plant.GMf.filt, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.corr.plantGmf.rds")

#filter main dataset to only those with matches
gmf2<-GBIF.miss.filt[which(GBIF.miss.filt$corr.GBIF.plant %in% corr.plant.GMf.filt$Submitted_Name),]

gmf2$Genus_in_database<-corr.plant.GMf.filt$Genus_in_database[match(gmf2$corr.GBIF.plant, corr.plant.GMf.filt$Submitted_Name)]

gmf2$Name_in_database<-corr.plant.GMf.filt$Name_in_database[match(gmf2$corr.GBIF.plant, corr.plant.GMf.filt$Submitted_Name)]
gmf2$Author_in_database<-corr.plant.GMf.filt$Author_in_database[match(gmf2$corr.GBIF.plant, corr.plant.GMf.filt$Submitted_Name)]
gmf2$New_name<-corr.plant.GMf.filt$New_name[match(gmf2$corr.GBIF.plant, corr.plant.GMf.filt$Submitted_Name)]
gmf2$New_author<-corr.plant.GMf.filt$New_author[match(gmf2$corr.GBIF.plant, corr.plant.GMf.filt$Submitted_Name)]
gmf2$Accepted_SPNAME<-corr.plant.GMf.filt$Accepted_SPNAME[match(gmf2$corr.GBIF.plant, corr.plant.GMf.filt$Submitted_Name)]


#look for tree genera
gmf2gentree<-gmf2[which(gmf2$Genus_in_database %in% globtree$genus),]
gmf2notree<-gmf2[which(!gmf2$Genus_in_database %in% globtree$genus),]
write.csv(gmf2gentree, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/missing.gmf2gentree.csv")

#filtered manually - only live, non-root tissue!
gmf2gentree.filt<-read.csv("intermediate_rds_csv/GBIF/GBIF_multiple_fields/missing.gmf2gentree.filt.csv")
gmf2gentree.filt<-gmf2gentree.filt[,2:31]

live.tissue.GBIF<-rbind.data.frame(gmf2notree, gmf2gentree.filt)
saveRDS(live.tissue.GBIF, "intermediate_rds_csv/GBIF/GBIF_multiple_fields/2025-01.missing.live.tissue.GBIF.rds")

live.tissue.GBIF.org<-live.tissue.GBIF[,c(2,6,7,9,10,12,13,16,17,18,19,24,25,26,27,28,29)]
live.tissue.GBIF.org$Accepted_SPAUTHOR<-live.tissue.GBIF.org$Author_in_database
live.tissue.GBIF.org$Accepted_SPAUTHOR[which(!is.na(live.tissue.GBIF.org$New_author))]<-live.tissue.GBIF.org$New_author[which(!is.na(live.tissue.GBIF.org$New_author))]

corrected.Host.sp<-live.tissue.GBIF.org$Accepted_SPNAME
corrected.Host.ssp<-word(live.tissue.GBIF$Accepted_SPNAME, 3, -1, sep=" ")
corrected.Host.auth<-live.tissue.GBIF.org$Accepted_SPAUTHOR
detected.Host.plant<-live.tissue.GBIF.org$GBIF.plant
accepted.Fungi.sp<-live.tissue.GBIF.org$species
accepted.Fungi.sp.full.name<-live.tissue.GBIF.org$acceptedScientificName
verbatim.Fungi.sp<-live.tissue.GBIF.org$verbatimScientificName
DB<-live.tissue.GBIF.org$db
habitat.and.plant.part.1<-live.tissue.GBIF.org$habitat
habitat.and.plant.part.2<-live.tissue.GBIF.org$associatedTaxa
habitat.and.plant.part.3<-as.character(live.tissue.GBIF.org$associatedOrganisms)
Location<-countrycode(live.tissue.GBIF.org$countryCode, origin='iso2c', destination='country.name')
Ref<-live.tissue.GBIF.org$occurrenceID
accepted.Fungi.auth<-word(live.tissue.GBIF.org$acceptedScientificName, 3, -1, sep=" ")

#make dataframe "GBIF.more" - in same column order as lp.database.TPL
GBIF.more.df<-cbind.data.frame(corrected.Host.sp, corrected.Host.ssp, corrected.Host.auth, 
                                   detected.Host.plant, accepted.Fungi.sp, accepted.Fungi.sp.full.name, accepted.Fungi.auth, verbatim.Fungi.sp, 
                                   habitat.and.plant.part.1, habitat.and.plant.part.2, habitat.and.plant.part.3, Location, Ref, DB)

#only keep unique fungi-plant pairs
GBIF.more.df.FIN<-GBIF.more.df %>%
  group_by(corrected.Host.sp, accepted.Fungi.sp) %>%
  distinct(corrected.Host.sp, .keep_all=TRUE)

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# III. Here we create GBIF.more.FIN.df, which will be used as raw data for the "On living plant tissue"  database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

saveRDS(GBIF.more.df.FIN, "Cleaned_data/2025-01.missing.GBIF.more.FIN.df.rds")

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
