#Combining data for Locally Co-occurring database


library(dplyr)
library(taxizedb)
library(U.Taxonstand)
library(stringr)
library(countrycode)

setwd("~/OneDrive - UBC/DB_fung_2025.01/")

#Read in raw data for co-occurrence
Globi.raw<-readRDS("intermediate_rds_csv/Rashmi_and_Globi/globi.local.fung.rds")
mycosphere.raw<-readRDS("intermediate_rds_csv/Rashmi_and_Globi/2024-12-21.UPD.myc.as4_MYCOSPHERE.NEW.4.2024.rds")
Mycoportal.raw<-readRDS("intermediate_rds_csv/Mycoportal/Myco.po.ab.sp.rds")
GBIF.raw<-readRDS('OSF/GBIF/GBIF.hits.sp.new.rds')
USDA.rawdata<-read.csv("OSF/USDA/Fungus-Host-Data_20211105.csv")
USDA.rawdata$db<-"USDA"
USDA.rawdata<-USDA.rawdata[-which(USDA.rawdata$host==""),]

### I ORIGINALLY DID NOT INCLUDE FUNGAL AUTH INFORMATION for Mycoportal, SO THIS CAN BE RETRIEVED FROM THE RAW DATA.
#read in authority key for mycoportal data (created in database.1b.mycoportal file)
Myco.po.ab<-readRDS('OSF/Mycoportal/Myco.po.ab.with.fungiauth.2024.rds')
#match to author names from raw data file.
Mycoportal.raw$Fungi.auth <- Myco.po.ab$scientificNameAuthorship[match(Mycoportal.raw$id, Myco.po.ab$id)]


#Combine into giant db.
detected.Host.plant<-c(USDA.rawdata$host, mycosphere.raw$Host, as.character(Mycoportal.raw$all.hits.sp.Mycop.plant), Globi.raw$globi.local.fung.V41, as.character(GBIF.raw$GBIF.plant))
accepted.Fungi.sp<-c(USDA.rawdata$sciName, mycosphere.raw$Fungi, Mycoportal.raw$scientificName, Globi.raw$globi.local.fung.V3, GBIF.raw$species)
accepted.Fungi.sp.full.name<-c(USDA.rawdata$sciName, paste(mycosphere.raw$Fungi, mycosphere.raw$Fungi.auth), paste(Mycoportal.raw$scientificName,Mycoportal.raw$Fungi.auth), Globi.raw$globi.local.fung.V3, GBIF.raw$acceptedScientificName)
verbatim.Fungi.sp<-c(USDA.rawdata$sciName, mycosphere.raw$Fungi, Mycoportal.raw$scientificName, Globi.raw$globi.local.fung.V3, GBIF.raw$verbatimScientificName)
habitat.and.plant.part.1<-c(USDA.rawdata$occurrenceRemarks, mycosphere.raw$Host.part, Mycoportal.raw$habitat, Globi.raw$globi.local.fung.V70, GBIF.raw$habitat)
#add  habitat.and.plant.part.2 and habitat.and.plant.part.3 later
DB<-c(USDA.rawdata$db, mycosphere.raw$db, Mycoportal.raw$db, Globi.raw$db, as.character(GBIF.raw$db))
Ref<-c(USDA.rawdata$litnum, mycosphere.raw$db, Mycoportal.raw$occurrenceID, Globi.raw$globi.local.fung.V2, GBIF.raw$occurrenceID)
L.raw.df<-cbind.data.frame(detected.Host.plant, accepted.Fungi.sp, accepted.Fungi.sp.full.name, verbatim.Fungi.sp, habitat.and.plant.part.1, DB, Ref)
#Add in additional habitat info to Raw data
L.raw.df$habitat.and.plant.part.2<-""
L.raw.df$habitat.and.plant.part.2[L.raw.df$DB=="Mycoportal.db"]<-Mycoportal.raw$associatedTaxa
L.raw.df$habitat.and.plant.part.2[L.raw.df$DB=="Globi.db"]<-Globi.raw$globi.local.fung.V37
L.raw.df$habitat.and.plant.part.2[L.raw.df$DB=="GBIF.db"]<-GBIF.raw$associatedTaxa
L.raw.df$habitat.and.plant.part.3<-""
L.raw.df$habitat.and.plant.part.3[L.raw.df$DB=="Mycoportal.db"]<-Mycoportal.raw$substrate
L.raw.df$habitat.and.plant.part.3[L.raw.df$DB=="Globi.db"]<-Globi.raw$globi.local.fung.V68
L.raw.df$habitat.and.plant.part.3[L.raw.df$DB=="GBIF.db"]<-GBIF.raw$associatedOrganisms
L.raw.df$Location<-c(USDA.rawdata$country, mycosphere.raw$Location, Mycoportal.raw$country, as.character(Globi.raw$db), countrycode(GBIF.raw$countryCode, origin='iso2c', destination='country.name'))



saveRDS(L.raw.df, "Local.DB/intermediate_rds/01-19.L.raw.df.rds")

#subset to GBIF/Myco data
Lrd.GBMyc<-L.raw.df[L.raw.df$DB=="GBIF.db"|L.raw.df$DB=="Mycoportal.db",]

#remove duplicate entries with same 6 letter beginning
Lrd.GBMyc$first.letters<-substr(Lrd.GBMyc$detected.Host.plant, 1, 6)
Lrd.GBMyc.f<-Lrd.GBMyc%>% 
  group_by(Ref, first.letters) %>%
  arrange(desc(nchar(as.character(detected.Host.plant)))) %>%
  distinct(first.letters, .keep_all=TRUE)

L.raw.df.GM.filt<-rbind.data.frame(setdiff(L.raw.df, Lrd.GBMyc[,-11]), Lrd.GBMyc.f[,-11])

L.raw.df.f1<-L.raw.df.GM.filt[-grep(" sp.", L.raw.df.GM.filt$detected.Host.plant),]
L.raw.df.f2<-L.raw.df.f1[-grep(" sp.", L.raw.df.f1$accepted.Fungi.sp),]
#remove unlisted genus entries
L.raw.df.f3<-L.raw.df.f2[-grep("^[A-Z]\\.\\s", L.raw.df.f2$detected.Host.plant),]
#remove one-word fungal entries (genera)
L.raw.df.f4 <- L.raw.df.f3[-grep("^\\w+$", L.raw.df.f3$accepted.Fungi.sp),]
L.raw.df.f5 <- L.raw.df.f4[-grep("^\\w+$", L.raw.df.f4$detected.Host.plant),]
L.raw.df.f6 <- L.raw.df.f5[!L.raw.df.f5$detected.Host.plant=="",]
L.raw.df.f7 <- L.raw.df.f6[!L.raw.df.f6$accepted.Fungi.sp=="",]
L.raw.df.f7$detected.Host.plant<-gsub("^ ", "", L.raw.df.f7$detected.Host.plant)
L.raw.df.f7<-L.raw.df.f7[-grep("A--", L.raw.df.f7$detected.Host.plant),]

L.raw.df.f7<-L.raw.df.f7[-which(L.raw.df.f7$accepted.Fungi.sp=="_ _"),]
#Fix input data so that U.Taxonstand is more accurate
L.raw.df.f7$sing.genus<-gsub("Populi","Populus",gsub("Abietis","Abies",gsub("Aceris","Acer",
                                                                             gsub("Aconiti","Aconitum",(gsub("Agropyrum","Agropyron",gsub("Agrostidis","Agrostis", gsub("Allii","Allium", gsub("Alni","Alnus",
                                                                                                                                                                                               gsub("Baccharidis", "Baccharis", gsub("Bromi","Bromus",gsub("Calamogrostidis","Calamogrostis",gsub("Caricis","Carex", gsub("Carpini","Carpinum",
                                                                                                                                                                                                                                                                                                                          gsub("Clerodendron","Clerodendrum", gsub("Corni","Cornus",gsub("Crataegi","Crataegus", gsub("Cytisi", "Cytisum", gsub("Delphini","Delphinium", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                gsub("Desmodii","Desmodium",gsub("Dipsaci", "Dipsacus",gsub("Elymi","Elymus",gsub("Euonymi","Euonymus", gsub("Eupatorii","Eupatorium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             gsub("Evonymi","Euonymus",gsub("Evonymus","Euonmyus",gsub("Fagi","Fagus",gsub("Fici","Ficus",gsub("Fraxini","Fraxinus",gsub("Galii","Galium",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         gsub("Gei","Geum",gsub("Gelsemii","Gelsemium",gsub("Geranii","Gernaium",gsub("Heraclei","Heracleum",gsub("Hordei","Hordeum",gsub("Humuli","Humulus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          gsub("Ilicis","Ilex",gsub("Impatientis","Impatiens",gsub("Iridis","Iris",gsub("Juglandis","Juglans",gsub("Junci","Juncus",gsub("Juniperi","Juniperus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         gsub("Lappae","Lappa",gsub("Lathyri","Lathyrus",gsub("Lauri","Laurus",gsub("Leontodontis","Leontodon",gsub("Ligustici","Ligusticum",gsub("Lupini","Lupinus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  L.raw.df.f7$detected.Host.plant))))))))))))))))))))))))))))))))))))))))))))))))
L.raw.df.f7$sing.genus<-gsub("Lycii","Lycium",gsub("Mori","Morus",gsub("Nasturtii","Nasturtium",gsub("Nerii","Nerium",gsub("Olerodendron","Clerodendron",
                                                                                                                            gsub("Ononidis","Ononis",gsub("Orobi","Orobus",gsub("Osmorrhizae","Osmorhiza",gsub("Panici","Panicum",gsub("Pentstamon","Penstemon",gsub("Pentstemonis","Penstemon",
                                                                                                                                                                                                                                                                     gsub("Peucedani","Peucedanum",gsub("Pini","Pinus",gsub("Plantaginis","Plantago",gsub("Polygoni","Polygonum",gsub("Polypodium","Polypodii",gsub("Populi","Populus",
                                                                                                                                                                                                                                                                                                                                                                                                                    gsub("Pruni","Prunus",gsub("Pteridis","Pteridium",gsub("Pyri","Pyrus",gsub("Rhamni","Rhamnus",gsub("Rhois","Rhus",gsub("Rhoidis","Rhus",gsub("Ribis","Ribes",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 gsub("Roripa","Rorippa",gsub("Rubi","Rubus",gsub("Rumicis","Rumex",gsub("Salicis","Salex",gsub("Sambuci","Sambucus",gsub("Scirpi","Scirpus",gsub("Senecionis","Senecio",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  gsub("Silphinium","Silphium",gsub("Sisymbrii","Sisymbrium",gsub("Smilacis","Smilax",gsub("Solani","Solanum",gsub("Sonchi","Sonchus",gsub("Symphyti","Symphytum",gsub("Tami","Tamus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       gsub("Tanaceti","Tanacetum",gsub("Teucrii","Teucrium",gsub("Thalictri","Thalictrum",gsub("Thesii","Theisum",gsub("Trifolii","Trifolium",gsub("Tritici","Triticum",gsub("Ulmi","Ulmus",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              gsub("Vaccinii","Vaccinium",gsub("Veronicae","Veronica",gsub("Viburni","Viburnum", gsub("Zizyphus","Ziziphus",L.raw.df.f7$detected.Host.plant)))))))))))))))))))))))))))))))))))))))))))))))))



saveRDS(L.raw.df.f7, "Local.DB/intermediate_rds/01-19.L.raw.df.f7.rds")

#remove entries already found in "on tissue" db - U.Taxonstand will be quicker
#this data we already dealt with in the living tissue dataset
Odb<-readRDS("On.DB/intermediate_rds/cleaned_combined_On_tissue_all3.rds")


#Now, isolate the ones that are not covered by the "pn plant database"
Not.on.plant<-anti_join(L.raw.df.f7, Odb, by = c("detected.Host.plant", "accepted.Fungi.sp"))
saveRDS(Not.on.plant, "Local.DB/intermediate_rds/05-05.Not.on.plant.rds")
Nop.host<-unique(Not.on.plant$sing.genus)

#Use taxizedb classification to remove nonplant hosts (omitted those that cause it to glitch out)
Nop.host.filt<-Nop.host[-c(17263, 51890, 61050, 978, 27082, 2681, 20780, 11598, 15283, 23309, 9085, 47309, 6508, 12632, 14318)]

class.Nh<-classification(Nop.host.filt, db="ncbi")

np<-Nop.host.filt[setdiff(grep("cellular organisms", class.Nh), grep("Viridiplantae", class.Nh))]

#get non-plant genera
np.gen<-unique(word(np, 1))

#remove them
Nop.rem<-which(word(Not.on.plant$sing.genus, 1) %in% np.gen)

Nop2<-Not.on.plant[-Nop.rem,]
saveRDS(Nop2, "Local.DB/intermediate_rds/05-05.Nop2.rds")

Nop.host.unst<-unique(Nop2$sing.genus)

#before standardizing, use the key from the on plant database.
Nil.op.host.c<-readRDS("On.DB/intermediate_rds/05-05.Nil.op.host.c.rds")
#these I can use later. However there are still more.
Nil.op.host.c.common<-Nil.op.host.c[Nil.op.host.c$Submitted_Name %in% Nop.host.unst,]
#here are the rest.
Nil.op.host.c2<-readRDS("On.DB/intermediate_rds/05-05.Nil.op.host.c2.rds")
#add the original names of the unaccepted (only due to not at species level) names
Nil.op.host.c2$Submitted_Name<-Nil.op.host.c$Submitted_Name[match(Nil.op.host.c2$Submitted_Name, Nil.op.host.c$Name_spLev)]


#These needed to be standardized.
Nop.host.unst.uncommon<-setdiff(Nop.host.unst, Nil.op.host.c.common$Submitted_Name)

#standardize names
Nop.host.uncommon.c<-nameMatch_WCVP(Nop.host.unst.uncommon)
saveRDS(Nop.host.uncommon.c, "Local.DB/intermediate_rds/05-05.Nop.host.uncommon.c.rds")
acc.Nop.host.c<-Nop.host.uncommon.c[which(!is.na(Nop.host.uncommon.c$Accepted_SPNAME)),]
spname.acc<-Nop.host.uncommon.c[which(!is.na(Nop.host.uncommon.c$Name_spLev)),]
#rerun standardization on species level names
Nop.host.uncommon.c2<-nameMatch_WCVP(spname.acc$Name_spLev)
saveRDS(Nop.host.uncommon.c2, "Local.DB/intermediate_rds/05-05.Nop.host.uncommon.c2.rds")

#add the original names of the unaccepted (only due to not at species level) names
Submitted_Name<-Nop.host.uncommon.c$Submitted_Name[match(Nop.host.uncommon.c2$Submitted_Name, Nop.host.uncommon.c$Name_spLev)]
adj.subm<-cbind.data.frame(Nop.host.uncommon.c2[,c(1,3:22)],Submitted_Name)
#I will include the intermediate name I used as the submitted name in the NOTE tab. (generally just the species name without subspecific epithets)
adj.subm$NOTE<-Nop.host.uncommon.c2[,2]

Kept.sp.Nop.host<-rbind.data.frame(acc.Nop.host.c, adj.subm, Nil.op.host.c.common, Nil.op.host.c2)
Kept.sp.Nop.host.filt<-Kept.sp.Nop.host[-which(is.na(Kept.sp.Nop.host$Accepted_SPNAME)),]

#filter main dataset to only those with matches
Nop3<-Nop2[which(Nop2$detected.Host.plant %in% Kept.sp.Nop.host.filt$Submitted_Name),]

#keep plants that were able to be taxonomically standardized
saveRDS(Nop3, "Local.DB/intermediate_rds/05-05.Nop3.rds")

#add in corrected taxonomic data
Nop3$Name_in_database<-Kept.sp.Nop.host$Name_in_database[match(Nop3$detected.Host.plant, Kept.sp.Nop.host$Submitted_Name)]
Nop3$Author_in_database<-Kept.sp.Nop.host$Author_in_database[match(Nop3$detected.Host.plant, Kept.sp.Nop.host$Submitted_Name)]
Nop3$NOTES<-Kept.sp.Nop.host$NOTES[match(Nop3$detected.Host.plant, Kept.sp.Nop.host$Submitted_Name)]
Nop3$New_author<-Kept.sp.Nop.host$New_author[match(Nop3$detected.Host.plant, Kept.sp.Nop.host$Submitted_Name)]
Nop3$corrected.Host.sp<-Kept.sp.Nop.host$Accepted_SPNAME[match(Nop3$detected.Host.plant, Kept.sp.Nop.host$Submitted_Name)]

Nop3$Accepted_SPAUTHOR<-Nop3$Author_in_database
Nop3$Accepted_SPAUTHOR[which(!is.na(Nop3$New_author))]<-Nop3$New_author[which(!is.na(Nop3$New_author))]

#Need to fix corrected.Host.sp in Odb to remove NAs
Odb.corr<-Odb[-which(is.na(Odb$corrected.Host.sp)),]

#Are there more combos that we found in the "on-plant"database now that we have standardized host names?
Not.op.revisit<-anti_join(Nop3, Odb.corr, by = c("corrected.Host.sp", "accepted.Fungi.sp"))
#We can forget the others, as we already have them in the living database.
Not.op.revisit$Accepted_SPNAME<-Not.op.revisit$corrected.Host.sp
saveRDS(Not.op.revisit, "Local.DB/intermediate_rds/05-05.Not.op.revisit.rds")




#remove badly formatted entries
Nop4<-Not.op.revisit[-grep("^× |\\+",Not.op.revisit$Accepted_SPNAME),]
#remove probable not plant detections
Nop5<-Nop4[-grep("^På|^Sida |^Zona|^Sii|^\\w+\\.", Nop4$detected.Host.plant),]

colnames(Nop5)[16]<-"corrected.Host.auth"
Nop5$corrected.Host.ssp<-word(Nop5$corrected.Host.sp, 3, -1)
Nop5$accepted.Fungi.auth<-word(Nop5$accepted.Fungi.sp.full.name, 3, -1)
#add back in ODB data
Nop5format<-Nop5[,c(1:10,15,16,18,19)]
Odb2<-data.frame(Odb, stringsAsFactors = FALSE)
combined_data<-rbind.data.frame(Odb2, Nop5format)
combined_data$DB[which(combined_data$DB=="Mycoportal")]<-"Mycoportal.db"

#Take only unique combinations. Sort with most reliable data for on-plant fungal occurrences first - 1. USDA 2. Mycosphere 3. Globi 4. Mycoportal 5. GBIF
comb_local_db<-comlod %>%
  arrange(match(DB, c("USDA", "Kushveer.Sarma.19", "Globi.db", "Mycoportal.db", "GBIF.db"))) %>%
  group_by(corrected.Host.sp, accepted.Fungi.sp) %>%
  distinct(corrected.Host.sp, .keep_all=TRUE)


#Change order so that fungarium works (it has trouble when there are fungi auth/no fungi auth)
comb_local_db <-comb_local_db %>% 
  arrange(desc(accepted.Fungi.auth))

comb_local_db2<-comb_local_db
#note, the row numbers may change if there are any changes to steps. skip the next 2 steps and trigger error to see the correct row numbers.
comb_local_db2$accepted.Fungi.sp.full.name[c(3735,6084,14117,53573)]<-c("Puccinia levis var. panici-sanguinalis (Rangel)","Puccinia substriata var. indica Ramachar","Uromyces sparganii subsp. asiaticus Parmelee","Pluteus atricapillus (Batsch) Fayod")
comb_local_db2$accepted.Fungi.auth[c(3735,6084,14117,53573)]<-c("(Rangel)","Ramachar","Parmelee","(Batsch) Fayod")
comb_local_db2$accepted.Fungi.auth<-gsub("^  ", "", comb_local_db2$accepted.Fungi.auth)
comb_local_db2$accepted.Fungi.auth<-gsub("  ", " ", comb_local_db2$accepted.Fungi.auth)
comb_local_db2$accepted.Fungi.auth<-gsub("^  ", "", comb_local_db2$accepted.Fungi.auth)
comb_local_db2$accepted.Fungi.auth<-gsub("^ ", "", comb_local_db2$accepted.Fungi.auth)
comb_local_db2$accepted.Fungi.auth<-gsub("  ", " ", comb_local_db2$accepted.Fungi.auth)
#note: slight diff in file and object name
saveRDS(comb_local_db2, "Local.DB/intermediate_rds/05-05.comb_local_db.rds")

#~#~#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
### Which are standardizable by Fungarium? ######
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#These have authorship info
t.loc.fun1<-taxon_update(head(comb_local_db2,265032), taxon_col="accepted.Fungi.sp.full.name", authorship_col="accepted.Fungi.auth", force_accepted=TRUE, species_only = TRUE, cores=3)

#these don't
t.loc.fun2<-taxon_update(tail(comb_local_db2, length(comb_local_db2)-265032), taxon_col="accepted.Fungi.sp.full.name", authorship_col=NULL, force_accepted=TRUE, species_only = TRUE, cores = 24)
t.loc.fun1$authorship<-""
t.loc.fun.c<-rbind.data.frame(t.loc.fun1, t.loc.fun2)
t.loc.fun.c$Fungal.Genus<-word(t.loc.fun.c$accepted.Fungi.sp, 1)


saveRDS(t.loc.fun.c, "Local.DB/intermediate_rds/05-05.t.loc.fun.c.rds")


#some fungi didnt get formatted correctly (end with na). redo fungarium on these.
fung.2format.repeat<-t.loc.fun.c[grep("error", t.loc.fun.c$error),1:14]
#fung.2format.repeat$accepted.Fungi.sp <- iconv(fung.2format.repeat$accepted.Fungi.sp, from = "ISO-8859-1", to = "UTF-8")
#some code error to fix
fung.2format.repeat$accepted.Fungi.sp[4946]<-"Uromyces sparganii subsp. asiaticus"
fung.2format.repeat$accepted.Fungi.sp<-gsub("Ã«|âˆšÃ‰Â¬Â´|Г«", "ë", fung.2format.repeat$accepted.Fungi.sp)
fung.2format.repeat$accepted.Fungi.sp<-gsub(" var\\. .*$| f\\. .*$", "", fung.2format.repeat$accepted.Fungi.sp)
fung.form.fix<-taxon_update(fung.2format.repeat, taxon_col="accepted.Fungi.sp", authorship_col=NULL, force_accepted=TRUE, species_only = TRUE, cores=)

t.loc.fun.c2<-rbind.data.frame(t.loc.fun.c[-grep("error", t.loc.fun.c$error),1:31], fung.form.fix)
t.loc.fun.c2$Fungal.Genus<-word(t.loc.fun.c2$accepted.Fungi.sp, 1)

saveRDS(t.loc.fun.c2, "Local.DB/intermediate_rds/05-05.t.loc.fun.c2.LOC.rds")

#remove corrected.host.sp =="NA
t.loc.fun.c2<-t.loc.fun.c2[-which(is.na(t.loc.fun.c2$corrected.Host.sp)),]

# Delete all "EXACT matches" below 85%.
# Delete all "FUZZY matches" below 94%.
t.loc.fun.c2$taxon_conf<-as.numeric(t.loc.fun.c2$taxon_conf)
t.loc.fun.c2$taxon_conf[which(is.na(t.loc.fun.c2$taxon_conf))]<-""
t.loc.fun.c2[t.loc.fun.c2$taxon_conf < 85, 15:32]<-""
t.loc.fun.c2[t.loc.fun.c2$taxon_conf < 94 & t.loc.fun.c2$taxon_matchtype=="FUZZY", 15:32]<-""


#there are still nonfungi!!
debatable<-t.loc.fun.c2[grep("Plantae|Protozoa|Animalia|Bacteria|Chromista", t.loc.fun.c2$new_kingdom),]
#but fungarium did not do a perfect job and some species in fungal genera are now Chromists, etc.
genus.to.remove<-setdiff(unique(debatable$Fungal.Genus), c("Botrytis","Asterina","Torula","Taphrina","Sphaerulina","Meliola","Patinella","Sphaeronema","Trichoderma", 
                                                           "Uredo", "Venturia", "Cryptococcus", "Cyathus", "Calloria", "Poria", "Lycoperdon", "Puccinia", "Uromyces", 
                                                           "Mucor", "Septoria", "Tremella", "Clavaria", "Auricularia", "Absidia", "Valsella","Caeoma", "Thyridium", "Phoma", 
                                                           "Marssonia", "Sphaeropsis", "Xyloma", "Monilia", "Aecidium", "Hexagonia", "Hymenula", "Collophora","Fenestella", 
                                                           "Cetraria", "Volutella", "Melanotus", "Dermea","Aegerita", "Dactylina", "Asteroma", "Helicoma", "Catenulopsora", 
                                                           "Eutypella", "Macrosporium", "Stictis","Phloeospora","Sphaeridium","Lisea"))

#delete the species themselves and the genera from the filtered database
loc_clean_unique_stand_filt<-setdiff(t.loc.fun.c2, debatable)
#capitalize fungi
loc_clean_unique_stand_filt$accepted.Fungi.sp<-str_to_sentence(loc_clean_unique_stand_filt$accepted.Fungi.sp)
loc_clean_unique_stand_filt2<-loc_clean_unique_stand_filt[which(!loc_clean_unique_stand_filt$Fungal.Genus %in% genus.to.remove), ]
#remove nonsp. hosts.
loc_clean_unique_stand_filt3<-as.data.frame(loc_clean_unique_stand_filt2[which(!str_count(loc_clean_unique_stand_filt2$corrected.Host.sp, "\\S+")==1),])

#### Resolve some issues with the text detector. This will not correct all, by any means.
#remove any with a 2 character detected name.

Filt.loc1<-loc_clean_unique_stand_filt3[which(!nchar(word(loc_clean_unique_stand_filt3$detected.Host.plant, 1))<=2),]
Filt.loc2<-Filt.loc1[which(!nchar(word(Filt.loc1$detected.Host.plant, 2))<=2),]

loc_clean_unique_stand_filt3<-Filt.loc2[-grep(" fores$|^Open |^Aca |^Anus |^Aquil |^Canal |^Euc |^Folia |^Folis |^Fra |^Inter |^Kap |^Oak- |^Paa |^Panda |^Plat |^Planted |^S-vendt |^Ser |^Sur | bos$| scrub$| prairie$| camp$| base$| mac$| leave$| pseud$| rec$| area$| small$| vid$| bei$| sur$| vel$| zone$| sub$| stam$| stammen$| sin$| fores$| adjacent$| amer$",
                          Filt.loc2$detected.Host.plant, ignore.case = TRUE),]


loc_clean_unique_stand_filt3<-loc_clean_unique_stand_filt3[-grep("Sphagnum", loc_clean_unique_stand_filt3$detected.Host.plant),]
loc_clean_unique_stand_filt3<-loc_clean_unique_stand_filt3[-grep("Musca ", loc_clean_unique_stand_filt3$detected.Host.plant),]
Final.loc_cleaned_unique_standardized_filtered<-loc_clean_unique_stand_filt3[which(!loc_clean_unique_stand_filt3$detected.Host.plant=="Rubus spec"),]

saveRDS(Final.loc_cleaned_unique_standardized_filtered, "Local.DB/intermediate_rds/05-05.Final.loc_cleaned_unique_standardized_filtered.rds")
#600924 interactions


#Without errors from fungal taxonomic standardization.
Final.loc_cleaned_unique_standardized_filtered.noerror<-Final.loc_cleaned_unique_standardized_filtered[!Final.loc_cleaned_unique_standardized_filtered$new_species=="",]
Final.loc_cleaned_unique_standardized_filtered.noerror.uniq<-Final.loc_cleaned_unique_standardized_filtered.noerror %>%
  group_by(corrected.Host.sp, new_species) %>%
  distinct(corrected.Host.sp, .keep_all=TRUE)
#263489 unique interactions
saveRDS(Final.loc_cleaned_unique_standardized_filtered.noerror.uniq[,1:31], "Final_comb_DBs/LOCAL.loc_cleaned_unique_standardized_filtered.noerror.uniq.rds")



