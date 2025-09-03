#Combining data for Fungi-on-plant database

setwd("~/OneDrive - UBC/DB_fung_2025.01/")
library(U.Taxonstand)
library(dplyr)
library(taxizedb)
library(stringr)
library(countrycode)

#Read in raw data for co-occurrence

globi.rawdata<-readRDS("intermediate_rds_csv/Rashmi_and_Globi/2024-12-21.UPD.globi.assoc.rds")
mycosphere.rawdata<-readRDS("intermediate_rds_csv/Rashmi_and_Globi/2024-12-21.UPD.myc.as4_MYCOSPHERE.NEW.4.2024.rds")
Mycoportal.rawdata<-readRDS('intermediate_rds_csv/Mycoportal/Myco.on.plant.FIN.rds')
GBIF.rawdata<-readRDS('intermediate_rds_csv/GBIF/GBIF.on.plant.FIN.rds')
GBIF2.rawdata<-readRDS('intermediate_rds_csv/GBIF/GBIF_multiple_fields/GBIF.m.on.plant.FIN.rds')
USDA<-read.csv('OSF/USDA/Fungus-Host-Data_20211105.csv')
USDA$db<-"USDA"


### I ORIGINALLY DID NOT INCLUDE FUNGAL AUTH INFORMATION for Mycoportal, SO THIS CAN BE RETRIEVED FROM THE RAW DATA.
#read in authority key for mycoportal data (created in database.1b.mycoportal file)
Myco.po.ab<-readRDS('OSF/Mycoportal/Myco.po.ab.with.fungiauth.2024.rds')
#match to author names from raw data file.
Mycoportal.rawdata$Fungi.auth <- Myco.po.ab$scientificNameAuthorship[match(Mycoportal.rawdata$id, Myco.po.ab$id)]


#Combine into giant db.
detected.Host.plant<-c(USDA$host, mycosphere.rawdata$Host, as.character(Mycoportal.rawdata$all.hits.sp.Mycop.plant), as.character(globi.rawdata$globi.fung.V41), GBIF.rawdata$GBIF.plant, as.character(GBIF2.rawdata$GBIF.plant))
accepted.Fungi.sp<-c(USDA$sciName, mycosphere.rawdata$Fungi, Mycoportal.rawdata$scientificName, as.character(globi.rawdata$globi.fung.V3), GBIF.rawdata$species, GBIF2.rawdata$species)
accepted.Fungi.sp.full.name<-c(USDA$sciName, paste(mycosphere.rawdata$Fungi, mycosphere.rawdata$Fungi.auth), paste(Mycoportal.rawdata$scientificName,Mycoportal.rawdata$Fungi.auth), as.character(globi.rawdata$globi.fung.V3), GBIF.rawdata$acceptedScientificName, GBIF2.rawdata$acceptedScientificName)
verbatim.Fungi.sp<-c(USDA$sciName, mycosphere.rawdata$Fungi, Mycoportal.rawdata$scientificName, as.character(globi.rawdata$globi.fung.V3), GBIF.rawdata$verbatimScientificName, GBIF2.rawdata$verbatimScientificName)
habitat.and.plant.part.1<-c(USDA$occurrenceRemarks, mycosphere.rawdata$Host.part, Mycoportal.rawdata$habitat, as.character(globi.rawdata$globi.fung.V70), GBIF.rawdata$habitat, GBIF2.rawdata$habitat)
#add  habitat.and.plant.part.2 and habitat.and.plant.part.3 later
DB<-c(USDA$db, mycosphere.rawdata$db, Mycoportal.rawdata$db, globi.rawdata$db, GBIF.rawdata$db, as.character(GBIF2.rawdata$db))
Ref<-c(USDA$litnum, mycosphere.rawdata$db, Mycoportal.rawdata$occurrenceID, as.character(globi.rawdata$globi.fung.V2), GBIF.rawdata$occurrenceID, GBIF2.rawdata$occurrenceID)
On.raw.df<-cbind.data.frame(detected.Host.plant, accepted.Fungi.sp, accepted.Fungi.sp.full.name, verbatim.Fungi.sp, habitat.and.plant.part.1, DB, Ref)
#Add in additional habitat info to Raw data
On.raw.df$habitat.and.plant.part.2<-""
On.raw.df$habitat.and.plant.part.2[On.raw.df$DB=="Mycoportal.db"]<-Mycoportal.rawdata$associatedTaxa
On.raw.df$habitat.and.plant.part.2[On.raw.df$DB=="Globi.db"]<-as.character(globi.rawdata$globi.fung.V37)
On.raw.df$habitat.and.plant.part.2[On.raw.df$DB=="GBIF.db"]<-c(GBIF.rawdata$associatedTaxa, GBIF2.rawdata$associatedTaxa)
On.raw.df$habitat.and.plant.part.3<-""
On.raw.df$habitat.and.plant.part.3[On.raw.df$DB=="Mycoportal.db"]<-Mycoportal.rawdata$substrate
On.raw.df$habitat.and.plant.part.3[On.raw.df$DB=="Globi.db"]<-as.character(globi.rawdata$globi.fung.V68)
On.raw.df$habitat.and.plant.part.3[On.raw.df$DB=="GBIF.db"]<-c(GBIF.rawdata$associatedOrganisms, GBIF2.rawdata$associatedOrganisms)
On.raw.df$Location<-c(USDA$country, mycosphere.rawdata$Location, Mycoportal.rawdata$country, as.character(globi.rawdata$db), countrycode(c(GBIF.rawdata$countryCode, GBIF2.rawdata$countryCode), origin='iso2c', destination='country.name'))


saveRDS(On.raw.df, "On.DB/intermediate_rds/01-19.On.raw.df.rds")

#subset to GBIF/Myco data
Ord.GBMyc<-On.raw.df[On.raw.df$DB=="GBIF.db"|On.raw.df$DB=="Mycoportal.db",]

#remove duplicate entries with same 6 letter beginning
Ord.GBMyc$first.letters<-substr(Ord.GBMyc$detected.Host.plant, 1, 6)
Ord.GBMyc.f<-Ord.GBMyc%>% 
  group_by(Ref, first.letters) %>%
  arrange(desc(nchar(as.character(detected.Host.plant)))) %>%
  distinct(first.letters, .keep_all=TRUE)

On.raw.df.GM.filt<-rbind.data.frame(setdiff(On.raw.df, Ord.GBMyc[,-11]), Ord.GBMyc.f[,-11])

On.raw.df.f1<-On.raw.df.GM.filt[-grep(" sp.", On.raw.df.GM.filt$detected.Host.plant),]
On.raw.df.f2<-On.raw.df.f1[-grep(" sp.", On.raw.df.f1$accepted.Fungi.sp),]
#remove unlisted genus entries
On.raw.df.f3<-On.raw.df.f2[-grep("^[A-Z]\\.\\s", On.raw.df.f2$detected.Host.plant),]
#remove one-word fungal entries (genera)
On.raw.df.f4 <- On.raw.df.f3[-grep("^\\w+$", On.raw.df.f3$accepted.Fungi.sp),]
On.raw.df.f5 <- On.raw.df.f4[-grep("^\\w+$", On.raw.df.f4$detected.Host.plant),]
On.raw.df.f6 <- On.raw.df.f5[!On.raw.df.f5$detected.Host.plant=="",]
On.raw.df.f7 <- On.raw.df.f6[!On.raw.df.f6$accepted.Fungi.sp=="",]
On.raw.df.f7$detected.Host.plant<-gsub("^ ", "", On.raw.df.f7$detected.Host.plant)
On.raw.df.f7<-On.raw.df.f7[-grep("A--", On.raw.df.f7$detected.Host.plant),]
#Fix input data so that U.Taxonstand is more accurate
On.raw.df.f7$sing.genus<-gsub("Populi","Populus",gsub("Abietis","Abies",gsub("Aceris","Acer",
                                                                       gsub("Aconiti","Aconitum",(gsub("Agropyrum","Agropyron",gsub("Agrostidis","Agrostis", gsub("Allii","Allium", gsub("Alni","Alnus",
                                                                                                                                                                                         gsub("Baccharidis", "Baccharis", gsub("Bromi","Bromus",gsub("Calamogrostidis","Calamogrostis",gsub("Caricis","Carex", gsub("Carpini","Carpinum",
                                                                                                                                                                                                                                                                                                                    gsub("Clerodendron","Clerodendrum", gsub("Corni","Cornus",gsub("Crataegi","Crataegus", gsub("Cytisi", "Cytisum", gsub("Delphini","Delphinium", 
                                                                                                                                                                                                                                                                                                                                                                                                                                          gsub("Desmodii","Desmodium",gsub("Dipsaci", "Dipsacus",gsub("Elymi","Elymus",gsub("Euonymi","Euonymus", gsub("Eupatorii","Eupatorium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       gsub("Evonymi","Euonymus",gsub("Evonymus","Euonmyus",gsub("Fagi","Fagus",gsub("Fici","Ficus",gsub("Fraxini","Fraxinus",gsub("Galii","Galium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   gsub("Gei","Geum",gsub("Gelsemii","Gelsemium",gsub("Geranii","Gernaium",gsub("Heraclei","Heracleum",gsub("Hordei","Hordeum",gsub("Humuli","Humulus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    gsub("Ilicis","Ilex",gsub("Impatientis","Impatiens",gsub("Iridis","Iris",gsub("Juglandis","Juglans",gsub("Junci","Juncus",gsub("Juniperi","Juniperus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   gsub("Lappae","Lappa",gsub("Lathyri","Lathyrus",gsub("Lauri","Laurus",gsub("Leontodontis","Leontodon",gsub("Ligustici","Ligusticum",gsub("Lupini","Lupinus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            On.raw.df.f7$detected.Host.plant))))))))))))))))))))))))))))))))))))))))))))))))
On.raw.df.f7$sing.genus<-gsub("Lycii","Lycium",gsub("Mori","Morus",gsub("Nasturtii","Nasturtium",gsub("Nerii","Nerium",gsub("Olerodendron","Clerodendron",
                                                                                                                      gsub("Ononidis","Ononis",gsub("Orobi","Orobus",gsub("Osmorrhizae","Osmorhiza",gsub("Panici","Panicum",gsub("Pentstamon","Penstemon",gsub("Pentstemonis","Penstemon",
                                                                                                                                                                                                                                                               gsub("Peucedani","Peucedanum",gsub("Pini","Pinus",gsub("Plantaginis","Plantago",gsub("Polygoni","Polygonum",gsub("Polypodium","Polypodii",gsub("Populi","Populus",
                                                                                                                                                                                                                                                                                                                                                                                                              gsub("Pruni","Prunus",gsub("Pteridis","Pteridium",gsub("Pyri","Pyrus",gsub("Rhamni","Rhamnus",gsub("Rhois","Rhus",gsub("Rhoidis","Rhus",gsub("Ribis","Ribes",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           gsub("Roripa","Rorippa",gsub("Rubi","Rubus",gsub("Rumicis","Rumex",gsub("Salicis","Salex",gsub("Sambuci","Sambucus",gsub("Scirpi","Scirpus",gsub("Senecionis","Senecio",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            gsub("Silphinium","Silphium",gsub("Sisymbrii","Sisymbrium",gsub("Smilacis","Smilax",gsub("Solani","Solanum",gsub("Sonchi","Sonchus",gsub("Symphyti","Symphytum",gsub("Tami","Tamus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 gsub("Tanaceti","Tanacetum",gsub("Teucrii","Teucrium",gsub("Thalictri","Thalictrum",gsub("Thesii","Theisum",gsub("Trifolii","Trifolium",gsub("Tritici","Triticum",gsub("Ulmi","Ulmus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        gsub("Vaccinii","Vaccinium",gsub("Veronicae","Veronica",gsub("Viburni","Viburnum", gsub("Zizyphus","Ziziphus",On.raw.df.f7$detected.Host.plant)))))))))))))))))))))))))))))))))))))))))))))))))



saveRDS(On.raw.df.f7, "On.DB/intermediate_rds/01-19.On.raw.df.f7.rds")
#remove entries already found in "living tissue" db - U.Taxonstand will be quicker
LT<-readRDS("Living.DB/intermediate_rds/cleaned_combined_living_tissue_all.rds")
#remove genera that got left here
LT.f <- LT[-grep("^\\w+$", LT$corrected.Host.sp),]
LT.f2<-LT.f %>%
  group_by(corrected.Host.sp, accepted.Fungi.sp) %>%
  distinct(corrected.Host.sp, .keep_all=TRUE)

#Note: the living tissue entries will be added in at end.

#Now, these are not necessarily covered by the living tissue database.
Not.in.living<-anti_join(On.raw.df.f7, LT.f2, by = c("detected.Host.plant", "accepted.Fungi.sp"))
saveRDS(Not.in.living, "On.DB/intermediate_rds/05-05.Not.in.living.rds")
Nil.host<-unique(Not.in.living$sing.genus)

#Use taxizedb classification to remove nonplant hosts (omitted those that cause it to glitch out)
Nil.host.filt<-Nil.host[-c(24658, 44781, 76049, 79635, 3906, 32522, 45992, 10449, 47931, 16553, 38382, 41273, 49771, 34978, 25975, 35650, 52942, 28183, 63222, 33472, 37985, 47598)]

class.Nh<-classification(Nil.host.filt, db="ncbi")

np<-Nil.host.filt[setdiff(grep("cellular organisms", class.Nh), grep("Viridiplantae", class.Nh))]

#get non-plant genera
np.gen<-unique(word(np, 1))

#remove them
Nil.rem<-which(word(Not.in.living$sing.genus, 1) %in% np.gen)

Nil.on.plant<-Not.in.living[-Nil.rem,]
saveRDS(Nil.on.plant, "On.DB/intermediate_rds/05-05.Nil.on.plant.rds")

Nil.op.host<-unique(Nil.on.plant$sing.genus)

#standardize names
Nil.op.host.c<-nameMatch_WCVP(Nil.op.host)
saveRDS(Nil.op.host.c, "On.DB/intermediate_rds/05-05.Nil.op.host.c.rds")
acc.Nil.op.host.c<-Nil.op.host.c[which(!is.na(Nil.op.host.c$Accepted_SPNAME)),]
spname.acc<-Nil.op.host.c[which(!is.na(Nil.op.host.c$Name_spLev)),]
#rerun standardization on species level names
Nil.op.host.c2<-nameMatch_WCVP(spname.acc$Name_spLev)
saveRDS(Nil.op.host.c2, "On.DB/intermediate_rds/05-05.Nil.op.host.c2.rds")

#add the original names of the unaccepted (only due to not at species level) names
Submitted_Name<-Nil.op.host.c$Submitted_Name[match(Nil.op.host.c2$Submitted_Name, Nil.op.host.c$Name_spLev)]
adj.subm<-cbind.data.frame(Nil.op.host.c2[,c(1,3:22)],Submitted_Name)
#I will include the intermediate name I used as the submitted name in the NOTE tab. (generally just the species name without subspecific epithets)
adj.subm$NOTE<-Nil.op.host.c2[,2]
Kept.sp.Nil.op.host<-rbind.data.frame(acc.Nil.op.host.c, adj.subm)

#filter main dataset to only those with matches
On.plant.filt<-Nil.on.plant[which(Nil.on.plant$detected.Host.plant %in% Kept.sp.Nil.op.host$Submitted_Name),]

#keep plants that were able to be taxonomically standardized
saveRDS(On.plant.filt, "On.DB/intermediate_rds/05-05.On.plant.filt.rds")
#some good data gets removed here, but mostly very messy, unusable data. Could potentially return to earlier files to find different ways of adding some data back in.

#add in corrected taxonomic data
On.plant.filt$Name_in_database<-Kept.sp.Nil.op.host$Name_in_database[match(On.plant.filt$detected.Host.plant, Kept.sp.Nil.op.host$Submitted_Name)]
On.plant.filt$Author_in_database<-Kept.sp.Nil.op.host$Author_in_database[match(On.plant.filt$detected.Host.plant, Kept.sp.Nil.op.host$Submitted_Name)]
On.plant.filt$NOTES<-Kept.sp.Nil.op.host$NOTES[match(On.plant.filt$detected.Host.plant, Kept.sp.Nil.op.host$Submitted_Name)]
On.plant.filt$New_author<-Kept.sp.Nil.op.host$New_author[match(On.plant.filt$detected.Host.plant, Kept.sp.Nil.op.host$Submitted_Name)]
On.plant.filt$corrected.Host.sp<-Kept.sp.Nil.op.host$Accepted_SPNAME[match(On.plant.filt$detected.Host.plant, Kept.sp.Nil.op.host$Submitted_Name)]

On.plant.filt$Accepted_SPAUTHOR<-On.plant.filt$Author_in_database
On.plant.filt$Accepted_SPAUTHOR[which(!is.na(On.plant.filt$New_author))]<-On.plant.filt$New_author[which(!is.na(On.plant.filt$New_author))]

#Are there more combos that we found in the living database now that we have standardized host names?
Not.in.living.revisit<-anti_join(On.plant.filt, LT.f2, by = c("corrected.Host.sp", "accepted.Fungi.sp"))
#We can forget the others, as we already have them in the living database.
saveRDS(Not.in.living.revisit, "On.DB/intermediate_rds/05-05.Not.in.living.revisit.rds")

Not.in.living.revisit$Accepted_SPNAME<-Not.in.living.revisit$corrected.Host.sp
#Take only unique combinations. Sort with most reliable data for on-plant fungal occurrences first - 1. USDA 2. Mycosphere 3. Globi 4. Mycoportal 5. GBIF
On.plant.filt2<-Not.in.living.revisit %>%
  arrange(match(DB, c("USDA", "Kushveer.Sarma.19", "Globi.db", "Mycoportal.db", "GBIF.db"))) %>%
  group_by(Accepted_SPNAME, accepted.Fungi.sp) %>%
  distinct(Accepted_SPNAME, .keep_all=TRUE)
saveRDS(On.plant.filt2, "On.DB/intermediate_rds/05-05.On.plant.filt2.rds")

#remove some entries that slipped through cracks before in GBIF data.
sil.GB<-readRDS("intermediate_rds_csv/GBIF/sil.GB.rds")
Gpd.GB<-readRDS("intermediate_rds_csv/GBIF/Gpd.GB.rds")
Gpf.GB<-readRDS("intermediate_rds_csv/GBIF/Gpf.GB.rds")
Gpw.GB<-readRDS("intermediate_rds_csv/GBIF/Gpw.GB.rds")
Gps.GB<-readRDS("intermediate_rds_csv/GBIF/Gps.GB.rds")
foG.GB<-readRDS("intermediate_rds_csv/GBIF/foG.GB.rds")
fwG.GB<-readRDS("intermediate_rds_csv/GBIF/fwG.GB.rds")
bdG.GB<-readRDS("intermediate_rds_csv/GBIF/bdG.GB.rds")
bcG.GB <-readRDS("intermediate_rds_csv/GBIF/bcG.GB.rds")
GB.hab.rem<-readRDS("intermediate_rds_csv/GBIF/GB.hab.rem.rds")
comb<-rbind.data.frame(sil.GB, Gpd.GB, Gpf.GB, Gpw.GB, Gps.GB, foG.GB, fwG.GB, bdG.GB, bcG.GB, GB.hab.rem)
comb$GBIF.plant<-as.character(comb$GBIF.plant)

On.plant.3<-anti_join(On.plant.filt2,comb, by = c("DB" = "db", "accepted.Fungi.sp" = "acceptedScientificName", "habitat.and.plant.part.1" = "habitat", "habitat.and.plant.part.2" = "associatedTaxa", "habitat.and.plant.part.3" = "associatedOrganisms", "detected.Host.plant" = "GBIF.plant", "Ref" = "occurrenceID"))
On.pl.3.GBIF.Myc<-On.plant.3[On.plant.3$DB=="GBIF.db"|On.plant.3$DB=="Mycoportal.db",]

#remove rows with duplicate entries. these usually have too many taxa in one field and are descriptions of areas.
On.pl.3.GBIF.Myc.B <- On.pl.3.GBIF.Myc %>%
  group_by(Ref, DB, accepted.Fungi.sp, habitat.and.plant.part.1) %>%
  filter(n() == 1)

#there are some good ones here, but we cannot keep the ones that match to a tree in the habitat description
On.pl.3.GBIF.Myc.keep.1<-setdiff(On.pl.3.GBIF.Myc, On.pl.3.GBIF.Myc.B)
#there are a bunch to keep in the associatedTaxa field (now hab...2)
On.pl.3.GBIF.Myc.keep.2<-On.pl.3.GBIF.Myc.keep.1[grep("host:| on |^on ", On.pl.3.GBIF.Myc.keep.1$habitat.and.plant.part.2, ignore.case = TRUE),]

#these were sorted by similarity between detected plant and the associated Taxa field
On.pl.3.GBIF.Myc.keep.2$dist<-mapply(adist, On.pl.3.GBIF.Myc.keep.2$sing.genus, On.pl.3.GBIF.Myc.keep.2$habitat.and.plant.part.2)
On.pl.3.GBIF.Myc.keep.2$grepp<-mapply(grep, On.pl.3.GBIF.Myc.keep.2$sing.genus, On.pl.3.GBIF.Myc.keep.2$habitat.and.plant.part.2)

On.pl.3.GBIF.Myc.keep.3<-On.pl.3.GBIF.Myc.keep.2 %>%
  arrange(dist) %>%
  arrange(grepp) %>%
  group_by(Ref, DB, accepted.Fungi.sp) %>%
  distinct(Ref, .keep_all=TRUE)

#keep those with matching names to host:
On.pl.3.GBIF.Myc.keep.4<-On.pl.3.GBIF.Myc.keep.3[On.pl.3.GBIF.Myc.keep.3$grepp=="1",]
#keep only those with string matching score 50 and under.
On.pl.3.GBIF.Myc.keep.5<-On.pl.3.GBIF.Myc.keep.4[which(On.pl.3.GBIF.Myc.keep.4$dist<51),]

#combine GBIF/Mycoportal remaining with the non-GBIF/Myc data with the living tissue data.
cleaned_combined_On_tissue_all<-rbind.data.frame(setdiff(On.plant.3, On.pl.3.GBIF.Myc), On.pl.3.GBIF.Myc.keep.5[,-c(18,19)])

corrected.Host.sp<-cleaned_combined_On_tissue_all$Accepted_SPNAME
corrected.Host.ssp<-word(cleaned_combined_On_tissue_all$Accepted_SPNAME, 3, -1, sep=" ")
corrected.Host.auth<-cleaned_combined_On_tissue_all$Accepted_SPAUTHOR
accepted.Fungi.auth<-word(cleaned_combined_On_tissue_all$accepted.Fungi.sp.full.name, 3, -1)

#make final dataframe with both corrected and uncorrected host names in separate columns.
cleaned_combined_On_tissue_all2<-cbind.data.frame(cleaned_combined_On_tissue_all, corrected.Host.sp, corrected.Host.ssp, corrected.Host.auth, accepted.Fungi.auth)

cleaned_combined_On_tissue_all3<-rbind.data.frame(LT.f2,cleaned_combined_On_tissue_all2[,-c(11:17)])
cleaned_combined_On_tissue_all3$Location[which(cleaned_combined_On_tissue_all3$Location=="Globi.db")]<-""

saveRDS(cleaned_combined_On_tissue_all3, "On.DB/intermediate_rds/cleaned_combined_On_tissue_all3.rds")

#Note, caught the USDA points without corrected Host sp, however the above was used in other databases, filtering before final database.
cleaned_combined_On_tissue_all4<-cleaned_combined_On_tissue_all3[-which(is.na(cleaned_combined_On_tissue_all3$corrected.Host.sp)),]
#remove probably not properly detected plants
cleaned_combined_On_tissue_all5<-cleaned_combined_On_tissue_all4[-grep("^På|^Sida |^Zona|^Sii|^\\w+\\.", cleaned_combined_On_tissue_all4$detected.Host.plant),]

saveRDS(cleaned_combined_On_tissue_all5, "On.DB/05-05.cleaned_combined_On_tissue_all5.rds")

# The following changes were made after standardization using Fungarium in the following section. skip ahead to standardize section first to determine chromists, then return
#capitalize fungi
cleaned_combined_On_tissue_all5$accepted.Fungi.sp<-str_to_sentence(cleaned_combined_On_tissue_all5$accepted.Fungi.sp)

#~#~#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
### Which are standardizable by Fungarium? ######
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#Change order so that fungarium works (it has trouble when there are fungi auth/no fungi auth)
cleaned_combined_On_tissue_all5 <-cleaned_combined_On_tissue_all5 %>% 
  arrange(desc(accepted.Fungi.auth))
cleaned_combined_On_tissue_all5$accepted.Fungi.auth<-gsub("^  ", "", cleaned_combined_On_tissue_all5$accepted.Fungi.auth)
cleaned_combined_On_tissue_all5$accepted.Fungi.auth<-gsub("  ", " ", cleaned_combined_On_tissue_all5$accepted.Fungi.auth)
cleaned_combined_On_tissue_all5$accepted.Fungi.auth<-gsub("^  ", "", cleaned_combined_On_tissue_all5$accepted.Fungi.auth)
cleaned_combined_On_tissue_all5$accepted.Fungi.auth<-gsub("^ ", "", cleaned_combined_On_tissue_all5$accepted.Fungi.auth)

#These have authorship info
t.on.fun1<-taxon_update(head(cleaned_combined_On_tissue_all5, 44339), taxon_col="accepted.Fungi.sp.full.name", authorship_col="accepted.Fungi.auth", force_accepted=TRUE, species_only = TRUE, cores=3)


#these don't
t.on.fun2<-taxon_update(tail(cleaned_combined_On_tissue_all5, (length(cleaned_combined_On_tissue_all5$corrected.Host.sp)-44339)), taxon_col="accepted.Fungi.sp.full.name", authorship_col=NULL, force_accepted=TRUE, species_only = TRUE, cores = 3)

t.on.fun1$authorship<-""
t.on.fun.c<-rbind.data.frame(t.on.fun1, t.on.fun2)

saveRDS(t.on.fun.c, "On.DB/intermediate_rds/05-05.t.on.fun.c.rds")
t.on.fun.c$Fungal.Genus<-word(t.on.fun.c$accepted.Fungi.sp, 1)

#some fungi didnt get formatted correctly (end with na). redo fungarium on these.
fung.2format.repeat<-t.on.fun.c[grep("error", t.on.fun.c$error),1:14]
#fung.2format.repeat$accepted.Fungi.sp <- iconv(fung.2format.repeat$accepted.Fungi.sp, from = "ISO-8859-1", to = "UTF-8")
#some code error to fix
fung.2format.repeat$accepted.Fungi.sp<-gsub("Ã«|âˆšÃ‰Â¬Â´|Г«", "ë", fung.2format.repeat$accepted.Fungi.sp)
fung.2format.repeat$accepted.Fungi.sp<-gsub(" var\\. .*$| f\\. .*$", "", fung.2format.repeat$accepted.Fungi.sp)
fung.form.fix<-taxon_update(fung.2format.repeat, taxon_col="accepted.Fungi.sp", authorship_col=NULL, force_accepted=TRUE, species_only = TRUE, cores=3)

t.on.fun.c2<-rbind.data.frame(t.on.fun.c[-grep("error", t.on.fun.c$error),], fung.form.fix)
t.on.fun.c2$Fungal.Genus<-word(t.on.fun.c2$accepted.Fungi.sp, 1)

saveRDS(t.on.fun.c2, "On.DB/intermediate_rds/05-05.t.on.fun.c2.ONPLANT.rds")

# Delete all "EXACT matches" below 86%.
# Delete all "FUZZY matches" below 94%.

t.on.fun.c2$taxon_conf<-as.numeric(t.on.fun.c2$taxon_conf)
t.on.fun.c2$taxon_conf[which(is.na(t.on.fun.c2$taxon_conf))]<-""
t.on.fun.c2[t.on.fun.c2$taxon_conf < 85, 15:32]<-""
t.on.fun.c2[t.on.fun.c2$taxon_conf < 94 & t.on.fun.c2$taxon_matchtype=="FUZZY", 15:32]<-""





#there are still nonfungi!!
debatable<-t.on.fun.c2[grep("Plantae|Protozoa|Animalia|Bacteria|Chromista", t.on.fun.c2$new_kingdom),]
#but fungarium did not do a perfect job and some species in fungal genera are now Chromists, etc.
genus.to.remove<-setdiff(unique(debatable$Fungal.Genus), c("Botrytis","Asterina","Torula","Taphrina","Sphaerulina","Meliola","Patinella",
                                                           "Sphaeronema","Trichoderma", "Uredo", "Venturia", "Cryptococcus", "Cyathus", "Sphaeropsis",
                                                           "Xyloma", "Hymenula", "Sphaeridium","Poria","Calloria", "Dermea","Aecidium","Septoria", 
                                                           "Puccinia", "Auricularia", "Melanotus", "Uromyces", "Stictis", "Macrosporium","Hexagonia","Asteroma","Dactylina", 
                                                           "Phloeospora","Clavaria", "Cetraria","Aegerita"))





#delete the species themselves and the genera from the filtered database
on_clean_unique_stand_filt<-setdiff(t.on.fun.c2, debatable)
#capitalize fungi
on_clean_unique_stand_filt$accepted.Fungi.sp<-str_to_sentence(on_clean_unique_stand_filt$accepted.Fungi.sp)
on_clean_unique_stand_filt2<-on_clean_unique_stand_filt[which(!on_clean_unique_stand_filt$Fungal.Genus %in% genus.to.remove), ]
#remove nonsp. hosts.
on_clean_unique_stand_filt3<-as.data.frame(on_clean_unique_stand_filt2[-which(str_count(on_clean_unique_stand_filt2$corrected.Host.sp, "\\S+")==1),])

on_clean_unique_stand_filt3<-on_clean_unique_stand_filt3[-grep("Sphagnum", on_clean_unique_stand_filt3$detected.Host.plant),]
on_clean_unique_stand_filt3<-on_clean_unique_stand_filt3[-grep("Musca ", on_clean_unique_stand_filt3$detected.Host.plant),]
Final.on_cleaned_unique_standardized_filtered<-on_clean_unique_stand_filt3[which(!on_clean_unique_stand_filt3$detected.Host.plant=="Rubus spec"),]

saveRDS(Final.on_cleaned_unique_standardized_filtered, "On.DB/05-05.Final.on_cleaned_unique_standardized_filtered.rds")
#345215 interactions

#Without errors from fungal taxonomic standardization.
Final.on_cleaned_unique_standardized_filtered.noerror<-Final.on_cleaned_unique_standardized_filtered[!Final.on_cleaned_unique_standardized_filtered$new_species=="",]
Final.on_cleaned_unique_standardized_filtered.noerror.uniq<-Final.on_cleaned_unique_standardized_filtered.noerror %>%
  group_by(corrected.Host.sp, new_species) %>%
  distinct(corrected.Host.sp, .keep_all=TRUE)
#251224 unique interactions
saveRDS(Final.on_cleaned_unique_standardized_filtered.noerror.uniq[,1:31], "Final_comb_DBs/ON.PLANT.standardized.cleaned.no.error.uniq.rds")

