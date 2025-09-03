# Database.II.Clean.2025.01.R
#updated version of Database.II.Clean.2023.12.R file and Database.II.Clean.2024.4.R file

setwd("~/OneDrive - UBC/DB_fung_2025.01/")

library(countrycode)
library(dplyr)
library(fungarium)
library(jsonlite)
library(stringr)
#import data

sph.Globi<-readRDS("Cleaned_data/FIN.sph.Globi.rds")
Mycoportal<-readRDS("Cleaned_data/2025-05.UPD.Myco.po.FIN.df.rds")
GBIF<-readRDS("Cleaned_data/2025-01.GBIF.df.rds")
GBIF.m<-readRDS("Cleaned_data/2025-01.missing.GBIF.more.FIN.df.rds")
USDA<-readRDS("Cleaned_data/USDA.FIN.df.rds")

all_cleaned_combined<-rbind.data.frame(USDA, GBIF.m, Mycoportal, sph.Globi, GBIF)

##########################################################################################
# We used this to make the other databases:
saveRDS(all_cleaned_combined, "Living.DB/intermediate_rds/cleaned_combined_living_tissue_all.rds")
###########################################################################################


# But we found a few more things to clean up.
#remove USDA points without corrected host species
all_cleaned_combined2<-all_cleaned_combined[-which(is.na(all_cleaned_combined$corrected.Host.sp)),]

#After looking at final data on funguild, there were a number of entries we did not want in there (like lichen hosts that got through)
#remove lichen genera detected in habitat fields
#read in lichenlist data from "Consortium of Lichen Herbaria" https://lichenportal.org/portal/checklists/checklist.php?clid=1492&pid=558
lichenlist<-read.csv("ext_data/World Checklist of Genera of Lichenized Fungi_1698006802.csv")
#make genus column
lichenlist$genus<-gsub(" sp.", "", lichenlist$ScientificName)

#
lichens1<-unlist(sapply(lichenlist$genus, function(genus) {
  which(grepl(genus, all_cleaned_combined2$habitat.and.plant.part.1, fixed=TRUE))
}))

lichens2<-unlist(sapply(lichenlist$genus, function(genus) {
  which(grepl(genus, all_cleaned_combined2$habitat.and.plant.part.2, fixed=TRUE))
}))

lichens3<-unlist(sapply(lichenlist$genus, function(genus) {
  which(grepl(genus, all_cleaned_combined2$habitat.and.plant.part.3, fixed=TRUE))
}))
#no lichens here

acc.lichens2<-all_cleaned_combined2[unique(c(lichens1, lichens2)),]
#some of these are fine, some are not
potential.lichens<-unique(c(lichens1, lichens2))
nolichens<-potential.lichens[-grep("Psoralea|Psoralidium|leaves", acc.lichens2$habitat.and.plant.part.1)]

#ditched the lichens.
all_cleaned_combined3<-all_cleaned_combined2[-potential.lichens,]
all_cleaned_combined3<-rbind.data.frame(all_cleaned_combined3, all_cleaned_combined2[nolichens,])

saveRDS(all_cleaned_combined3, "Living.DB/intermediate_rds/all_cleaned_combined3.rds")

all_cleaned_combined3.unique<-all_cleaned_combined3 %>%
  group_by(corrected.Host.sp, accepted.Fungi.sp) %>%
  distinct(corrected.Host.sp, .keep_all=TRUE)

saveRDS(all_cleaned_combined3.unique, "Living.DB/intermediate_rds/all_cleaned_combined3.unique.rds")


#~#~#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
### Which are standardizable by Fungarium? ######
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#Change order so that fungarium works (it has trouble when there are fungi auth/no fungi auth)
all_cleaned_combined3.unique <-all_cleaned_combined3.unique %>% 
  arrange(desc(accepted.Fungi.auth))

#These have authorship info
t.all.fun<-taxon_update(head(all_cleaned_combined3.unique, 32590), taxon_col="accepted.Fungi.sp.full.name", authorship_col="accepted.Fungi.auth", force_accepted=TRUE, species_only = TRUE, cores=3)
#these don't
t.all.fun2<-taxon_update(tail(all_cleaned_combined3.unique, 20759), taxon_col="accepted.Fungi.sp.full.name", authorship_col=NULL, force_accepted=TRUE, species_only = TRUE, cores = 3)
t.all.fun$authorship<-""
t.all.fun.c<-rbind.data.frame(t.all.fun, t.all.fun2)
#save later
saveRDS(t.all.fun.c, "Living.DB/intermediate_rds_csv/t.all.fun.c.rds")

#some fungi didnt get formatted correctly (end with na). redo fungarium on these.
fung.2format.repeat<-t.all.fun.c[grep("error", t.all.fun.c$error),1:14]
fung.2format.repeat$accepted.Fungi.sp<-gsub("Ã«|âˆšÃ‰Â¬Â´|Г«", "ë", fung.2format.repeat$accepted.Fungi.sp)
fung.2format.repeat$accepted.Fungi.sp<-gsub(" var\\. .*$| f\\. .*$", "", fung.2format.repeat$accepted.Fungi.sp)
fung.form.fix<-taxon_update(fung.2format.repeat, taxon_col="accepted.Fungi.sp", authorship_col=NULL, force_accepted=TRUE, species_only = TRUE, cores=3)

t.all.fun.c2<-rbind.data.frame(t.all.fun.c[-grep("error", t.all.fun.c$error),], fung.form.fix)
t.all.fun.c2$Fungal.Genus<-word(t.all.fun.c2$accepted.Fungi.sp, 1)



#capitalize fungi
t.all.fun.c2$accepted.Fungi.sp<-str_to_sentence(t.all.fun.c2$accepted.Fungi.sp)

saveRDS(t.all.fun.c2, "Living.DB/intermediate_rds_csv/t.all.fun.c2.rds")

#remove nonsp. hosts.
t.all.fun.c3<-as.data.frame(t.all.fun.c2[-which(str_count(t.all.fun.c2$corrected.Host.sp, "\\S+")==1),])
t.all.fun.c4<-t.all.fun.c3[-grep("sect. ", t.all.fun.c3$corrected.Host.sp),]

#check for ultrarare outliers 
t.all.fun.c5<-t.all.fun.c4[-grep("Sphagnum", t.all.fun.c4$detected.Host.plant),]
t.all.fun.c5<-t.all.fun.c5[-grep("Musca ", t.all.fun.c5$detected.Host.plant),]
t.all.fun.c5<-t.all.fun.c5[which(!t.all.fun.c5$detected.Host.plant=="Rubus spec"),]
t.all.fun.c5<-t.all.fun.c5[which(!t.all.fun.c5$new_class=="Geoglossomycetes"),]
t.all.fun.c5<-t.all.fun.c5[which(!t.all.fun.c5$new_class=="Candelariomycetes"),]


#there are still nonfungi!!
debatable<-t.all.fun.c5[grep("Plantae|Protozoa|Animalia|Bacteria|Chromista", t.all.fun.c5$new_kingdom),]

#but fungarium did not do a perfect job and some species in fungal genera are now Chromists, etc.
genus.to.remove<-setdiff(unique(debatable$Fungal.Genus), c("Botrytis","Asterina", "Asteroma", "Torula","Taphrina","Sphaerulina","Meliola","Patinella","Sphaerulina", "Sphaeronema","Trichoderma", "Uredo"))

#remove the non-fungal entries
t.all.fun.c6<-t.all.fun.c5[which(!t.all.fun.c5$Fungal.Genus %in% genus.to.remove), ]

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# QUICK DETOUR: Make copy of unique interactions DB without fungal standardization
Final.unstandardizedFungi.LTDB.uniq<-t.all.fun.c6 %>%
  group_by(corrected.Host.sp, accepted.Fungi.sp) %>%
  distinct(corrected.Host.sp, .keep_all=TRUE)
saveRDS(Final.unstandardizedFungi.LTDB.uniq, "intermediate_rds_csv/Final.unstandardizedFungi.LTDB.uniq.rds")
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#BACK TO STANDARDIZED DB.
#Without errors from fungal taxonomic standardization.
No.error.LT<-t.all.fun.c6[grep("error", t.all.fun.c6$error, invert = TRUE),]

write.csv(No.error.LT, "Living.DB/intermediate_rds_csv/No.error.LT.csv")

#edit the fungi (in "No.error.LT.csv") that did not get properly transcribed (taxon_conf<100 and put in different species than listed)
# For all "EXACT matches" below 85%, I made sure the name of the genus +specific epithet between the query and accepted.fungi.sp matched, otherwise I erased the fungarium data. 
# For all "FUZZY matches" below 85%, I deleted them. 85%-93% I looked over and deleted entries that did not appear to match the correct taxa.
F.No.error.LT<-read.csv("intermediate_rds_csv/FILT.No.error.LT.csv")
F.No.error.LT<-F.No.error.LT[,-1]
F.No.error.LT<-F.No.error.LT[which(F.No.error.LT$new_kingdom=="Fungi"),]

Final.standard.uniq.LT.db<-F.No.error.LT %>%
  group_by(corrected.Host.sp, new_species) %>%
  distinct(corrected.Host.sp, .keep_all=TRUE)
#42191 unique interactions
saveRDS(Final.standard.uniq.LT.db[,1:31], "Final_comb_DBs/LIVING.standardized.cleaned.no.error.uniq.rds")





