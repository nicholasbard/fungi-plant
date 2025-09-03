setwd("~/OneDrive - UBC/DB_fung_2025.01/")
library(dplyr)
library(taxizedb)
library(stringr)
library(U.Taxonstand)
#downloaded on 2023-11-05 from https://data.nal.usda.gov/dataset/united-states-national-fungus-collections-fungus-host-dataset

#includes info up to Nov 21 2021
#only data with occurrence remarks
USDA<-read.csv("OSF/USDA/Fungus-Host-Data_20211105.csv")

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# I & II The USDA data is already ready for the "Locally Co-occurring" and Fungi-on-plant" databases.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# write.csv("OSF/USDA/Fungus-Host-Data_20211105.csv")

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#remove those where occurrence remarks not present (there are many of these)
USDA.occ<-USDA[which(!USDA$occurrenceRemarks == ""),]
remarks<-as.character(unique(USDA.occ$occurrenceRemarks))

u.remov<-remarks[grep("parasitic on scale insects|ponderosa pine|dried|stor|dead|dry|on roots.|wood|bark|base of|fallen|mortu|under|from root|on rhizome|decay|resin|insect|frass|rhizosphere|decline|log|^rot|trunk|on limb|on branch|postharvest|stump|mummified|cutting|^in root|muricata pine|^yellows.|soil|wilting|rott|sapro|type specimen|herbarium|on pines|longleaf pine|hyperparasit|sugar pine|tuber|from necro|rhinanthacearum|dying|on cortic|timber|killed branches|reported on conifers|limb|submerged leaves|Petr|Tulipa edulis|canary island pine|Cat-tail|epiphy|II, III.|leaves old|loblolly pine|lodgepole pine|maritime pine|mycorrhiza|muricate pine|not stated|meliola dieffenbachiae|leaf litter|oidim|patula pine|tanaka |reported on conifers|treated utility poles|on rust", remarks, ignore.case=TRUE)]


filt.rem<-setdiff(remarks,u.remov)
filt.rem.d<-data.frame(filt.rem)
filt.rem.d$occurrenceRemarks<-as.character(filt.rem.d$filt.rem)
USDA.occ$occurrenceRemarks<-as.character(USDA.occ$occurrenceRemarks)
#filtered dataset
filt.rem.df<-semi_join(USDA.occ, filt.rem.d, by="occurrenceRemarks")

#distinct fungi/plant combos
USDA.fr.dis.Pl.F<-filt.rem.df %>%
  distinct(sciName, host, .keep_all = TRUE)
#remove vague entries
USDA.fdPF.nsp<-USDA.fr.dis.Pl.F[grep(" sp.", USDA.fr.dis.Pl.F$host),]
USDA.fdPF.sp<-setdiff(USDA.fr.dis.Pl.F, USDA.fdPF.nsp)
saveRDS(USDA.fdPF.sp, "intermediate_rds_csv/USDA/2025-01.UPD.USDA.fdPF.sp.rds")

#remove fungi not at species level
USDA.ns<-USDA.fdPF.sp[grep(" sp.", USDA.fdPF.sp$sciName, invert=TRUE),]

#remove nonfungi from fungi column (oomycetes) using taxon_update (fungarium)
USfun.r<-as.character(unique(USDA.fdPF.sp$sciName))
#remove non species fungi (sp.)
USfun.rv<-USfun.r[grep(" sp.", USfun.r, invert=TRUE)]
#use taxizedb to detect nonfungi (oomycetes)
class.USfrv<-classification(USfun.rv, db="ncbi")
#what are the non fungi?
class.nonfun<-USfun.rv[grep("Amoebozoa|Stramenopiles", class.USfrv)]
#what are their genera?
nonfun.gen<-unique(word(class.nonfun, 1))

USfun.filt1<-USfun.rv

#grep all the nonfungal genera out
nonfung.index<-NULL
for (i in nonfun.gen){
  mm<-ifelse(nn<-grep(i, USDA.ns$sciName, ignore.case = TRUE),print(i),0)
  nonfung.index<-c(nonfung.index, nn)
}
#filter out non fungal entries
USDA.ns.filt<-USDA.ns[-nonfung.index,]
saveRDS(USDA.ns.filt, "intermediate_rds_csv/USDA/2025-01.UPD.USDA.ns.filt.rds")

#Check WCVP for host plant names

#removes cultivar name
USDA.ns.filt$hostsp<-gsub(" cv. \\w+", "", USDA.ns.filt$host)
my.h2<-unique(USDA.ns.filt$hostsp)
corr.plant<-nameMatch_WCVP(my.h2)

US1<-USDA.ns.filt[which(USDA.ns.filt$hostsp %in% corr.plant$Submitted_Name),]
#all were accepted sp names

US1$Genus_in_database<-corr.plant$Genus_in_database[match(US1$hostsp, corr.plant$Submitted_Name)]
US1$Name_in_database<-corr.plant$Name_in_database[match(US1$hostsp, corr.plant$Submitted_Name)]
US1$Author_in_database<-corr.plant$Author_in_database[match(US1$hostsp, corr.plant$Submitted_Name)]
US1$New_name<-corr.plant$New_name[match(US1$hostsp, corr.plant$Submitted_Name)]
US1$New_author<-corr.plant$New_author[match(US1$hostsp, corr.plant$Submitted_Name)]
US1$Accepted_SPNAME<-corr.plant$Accepted_SPNAME[match(US1$hostsp, corr.plant$Submitted_Name)]

US1$Accepted_SPAUTHOR<-US1$Author_in_database
US1$Accepted_SPAUTHOR[which(!is.na(US1$New_author))]<-US1$New_author[which(!is.na(US1$New_author))]

saveRDS(US1, "intermediate_rds_csv/USDA/2025-1.USDA.US1.rds")


#remove lichens
#read in lichenlist data from "Consortium of Lichen Herbaria" https://lichenportal.org/portal/checklists/checklist.php?clid=1492&pid=558
lichenlist<-read.csv("ext_data/World Checklist of Genera of Lichenized Fungi_1698006802.csv")
#make genus column
lichenlist$genus<-gsub(" sp.", "", lichenlist$ScientificName)
#filter lichenizing fungi
hfilt.US1<-US1[which(!word(US1$sciName) %in% lichenlist$genus),]
#remove fungal taxon entries that contain word root "lichen"
nlaf<-hfilt.US1[grep("lichen",hfilt.US1$sciName, ignore.case=TRUE),]
hfilt.US2<-hfilt.US1[which(!hfilt.US1$sciName %in% nlaf$sciame),]


#Remove "stems" and "shoots" from tree.
globtree<-read.csv("ext_data/global_tree_search_trees_1_7.csv")
#filter by genus in order to eliminate as many shrubby genera, as well.
gltrgenus<-data.frame(str_split_fixed(globtree$TaxonName, " ", 2))
globtree$genus<-gltrgenus$X1
#inspect the occurrence
US1tree.gen<-hfilt.US2[which(word(as.character(hfilt.US2$Genus_in_database), 1) %in% globtree$genus),]
#we kept more tissue types for the woody plants here (e.g., stem), as these all appear to be pathogens that infect living tissue.
US3b<-hfilt.US2[!hfilt.US2$occurrenceRemarks=="ponderosa pine",]
US3c<-US3b[!US3b$occurrenceRemarks=="monterey pine, insignis pine",]
US3d<-US3c[!US3c$occurrenceRemarks=="slash pine",]


#only keep unique fungi-plant pairs
US4<-US3d %>%
  group_by(Accepted_SPNAME, sciName) %>%
  distinct(Accepted_SPNAME, .keep_all=TRUE)

corrected.Host.sp<-US4$Accepted_SPNAME
corrected.Host.ssp<-word(US4$Accepted_SPNAME, 3, -1, sep=" ")
corrected.Host.auth<-US4$Accepted_SPAUTHOR
detected.Host.plant<-US4$host
accepted.Fungi.sp<-US4$sciName
accepted.Fungi.sp.full.name<-US4$sciName
verbatim.Fungi.sp<-US4$sciName
habitat.and.plant.part.1<-US4$occurrenceRemarks
habitat.and.plant.part.2<-""
habitat.and.plant.part.3<-""
Fungi.auth<-""
DB<-"USDA"
Location<-US4$country
Ref<-US4$litnum
accepted.Fungi.auth<-""

#make final dataframe with both corrected and uncorrected host names in separate columns.
USDA.FIN.df<-cbind.data.frame(corrected.Host.sp, corrected.Host.ssp, corrected.Host.auth, 
                               detected.Host.plant, accepted.Fungi.sp, accepted.Fungi.sp.full.name, accepted.Fungi.auth, verbatim.Fungi.sp, 
                               habitat.and.plant.part.1, habitat.and.plant.part.2, habitat.and.plant.part.3, Location, Ref, DB)

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# III. Here we create USDA.FIN.df, which will be used as raw data for the "On living plant tissue"  database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

saveRDS(USDA.FIN.df, "Cleaned_data/USDA.FIN.df.rds")

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#