# new update for 2025.01 (orig 2023.11), which considers newer CABI file, adjusts mycosphere authority spacing.

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


#~#~#~#~#~#~#~#~#~#~#~#~# Mycosphere #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

# get host associations from Rashmi, Kushveer and Sarma 2019.
# I uploaded to smallpdf.com to convert pdf table to a readable excel file, then I needed to reformat slightly in Excel for everything to sit in the right column.
# Then I handtweaked, because there were a few slight formatting errors.
#deleted sponge, lichen, algae substrates. Deleted host information not at the species level and abbreviated species names I could not ascertain from the data provided. 
# Retained root records, as I could be more sure these were not mycorrhizal.
mycosphere <- read.csv("intermediate_rds_csv/Rashmi_and_Globi/FILT.ed.2023.11.06.Mycosphere.data.formatted.via.SMALLPDF.and.unmerge.csv")
mycosph.assoc <- data.frame(Host=mycosphere$Host, Endophytes=mycosphere$Endophytes, Location=mycosphere$Location, Host.part=mycosphere$Host.Part, db=rep("Rashmi.Kushveer.Sarma.19"))

# separate authority from fungi names
mfun<-data.frame(str_split_fixed(mycosph.assoc$Endophytes, ", ", n=2))
myc.as2<-cbind.data.frame(mfun, mycosph.assoc[,-2])
colnames(myc.as2)[1]<-"Fungi"
colnames(myc.as2)[2]<-"Fungi.auth"
#create new rows for each host.
myc.as3<-myc.as2 %>%
  separate_rows(Host, convert = TRUE, sep=", ")
myc.as3<-data.frame(lapply(myc.as3, as.character), stringsAsFactors = FALSE)
myc.as3$hab2<-""
myc.as3$Location<-gsub("C1", "Asia", 
                       gsub("C2", "Africa",
                            gsub("C3", "North America",
                                 gsub("C4", "South America",
                                      gsub("C5", "Antarctica",
                                           gsub("C6", "Europe",
                                                gsub("C7", "Australia",
                                                     myc.as3$Location)))))))

saveRDS(myc.as3, "intermediate_rds_csv/Rashmi_and_Globi/myc.as3_MYCOSPHERE.NEW.4.2024.rds")
myc.as4<-myc.as3[grep(" sp\\.", myc.as3$Fungi, invert=TRUE),]
saveRDS(myc.as4, "intermediate_rds_csv/Rashmi_and_Globi/2024-12-21.UPD.myc.as4_MYCOSPHERE.NEW.4.2024.rds")

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# I & II Here we create myc.as4, which will be used as raw data for the "Locally Co-occurring" and Fungi-on-plant"  databases.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# saveRDS(myc.as4, "intermediate_rds_csv/Rashmi_and_Globi/2024-12-21.UPD.myc.as4_MYCOSPHERE.NEW.4.2024.rds")

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~#~#~#~#~#~#~#~#~#~#~#~# Globi #~#~#~#~#~#~#~#~#~#~#~#~#~#~#~
# retrieved May 7, 2022.

#obtain species interactions from GLOBI
#on website https://www.globalbioticinteractions.org/data
#get "epiphyteOf, livesOn, hasHost, livesInsideOf, parasiteOf, parasitoidOf, mutualistOf, symbiontOf"

# globi <- read.csv("~/Desktop/Fungi.interactions.Globi.csv", header= FALSE)

#Use for locally co-occurring
globi.local.fung <- globi %>% 
  filter(globi$V20=="Fungi" &
           globi$V58=="Plantae")
globi.local.fung <- data.frame(globi.local.fung$V3, globi.local.fung$V4, globi.local.fung$V37, globi.local.fung$V41, globi.local.fung$V42, globi.local.fung$V2, globi.local.fung$V70, globi.local.fung$V68, globi.local.fung$V38, globi.local.fung$V84, globi.local.fung$V82, db=rep("Globi.db"))

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# I. Here we create globi.assoc, which will be used as raw data for the "Locally-Co-occurring"  database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

# saveRDS(globi.local.fung, "intermediate_rds_csv/Rashmi_and_Globi/globi.local.fung.rds")

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#



globi.fung <- globi %>% 
     filter(globi$V20=="Fungi" &
            globi$V58=="Plantae" &
              !globi$V37=="eats" &
              !globi$V37=="hasHabitat" &
              !globi$V37=="interactsWith" &
              !globi$V37=="adjacentTo" &
              !globi$V37=="livesNear" &
              !globi$V37=="hasVector" &
              !globi$V37=="visitsFlowersOf" &
              !globi$V37=="visits" &
              !globi$V37=="coOccursWith" &
              !globi$V4=="kingdom" &
              !globi$V4=="phylum" &
              !globi$V4=="class" &
              !globi$V4=="order" &
              !globi$V4=="family")

# saveRDS(globi.fung, "intermediate_rds_csv/Rashmi_and_Globi/globi.fung.rds")

#subset to meaningful variables
globi.assoc <- data.frame(globi.fung$V3, globi.fung$V4, globi.fung$V37, globi.fung$V41, globi.fung$V42, globi.fung$V2, globi.fung$V70, globi.fung$V68, globi.fung$V38, globi.fung$V84, globi.fung$V82, db=rep("Globi.db"))

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# II. Here we create globi.assoc, which will be used as raw data for the "Fungi-on-plant"  database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
globi.assoc <- data.frame(globi.fung$V3, globi.fung$V4, globi.fung$V37, globi.fung$V41, globi.fung$V42, globi.fung$V2, globi.fung$V70, globi.fung$V68, globi.fung$V38, globi.fung$V84, globi.fung$V82, db=rep("Globi.db"))
# saveRDS(globi.assoc, "intermediate_rds_csv/Rashmi_and_Globi/2024-12-21.UPD.globi.assoc.rds")

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


######for living tissue:
#keep all parasites, parasitoids, and lives-inside-ofs.
globi.assoc3<-globi.assoc2[grep("Inside", globi.assoc2$globi.fung.V37, ignore.case=TRUE),]

#only keep symbionts clearly inside plant tissue.
globi.assoc4<-globi.assoc2[grep("parasit|Host|livesOn|symbiont|mutualist", globi.assoc2$globi.fung.V37, ignore.case=TRUE),]

#eliminate tree data unless specific about plant part, in order to get rid of most ectomycorrhizal.

#I retrieved the csv for the GlobalTreeSearch data on Oct 31 2023. https://rpubs.com/Roeland-KINDT/812716
#read this data in 
globtree<-read.csv("ext_data/global_tree_search_trees_1_7.csv")
#filter by genus in order to eliminate as many shrubby genera, as well.
gltrgenus<-data.frame(str_split_fixed(globtree$TaxonName, " ", 2))
globtree$genus<-gltrgenus$X1
#remove the entries without explicit living plant parts from the tree, but allow on non-woody host plants.
globi.assoc4tree<-globi.assoc4[which(word(as.character(globi.assoc4$globi.fung.V41), 1) %in% globtree$genus),]
globi.assoc5tree<-globi.assoc4tree[which(!globi.assoc4tree$globi.fung.V70 == ""),]
globi.assoc6tree<-globi.assoc5tree[grep("lea|sor|petiole|sheath|anther|fruit|calyx|internode|first year|midrib|cone|peduncle|capitul|acervulus|bud|berry|flower|culm|pod|ovary|appressorium|green parts|spikelet|pycnium|aecium|telium|sepal|pedicel|catkin|flor|fruto|hoja|semilla|pycnothyrium|frond|caespituli|needle|seed|stolon|twig", globi.assoc5$globi.fung.V70, ignore.case=TRUE),]
globi.assoc4notree<-globi.assoc4[which(!word(as.character(globi.assoc4$globi.fung.V41), 1) %in% globtree$genus),]
globi.living<-rbind.data.frame(globi.assoc3, globi.assoc4notree, globi.assoc6tree)


globi.living<-rbind.data.frame(globi.assoc3, globi.assoc6)

globi.living<-data.frame(lapply(globi.living, as.character), stringsAsFactors = FALSE)
#only on fungal sp
globi.living.sp<-globi.living[which(!globi.living$globi.fung.V4 == "genus"),]
globi.living.sp<-globi.living.sp[grep(" sp\\.", globi.living.sp$globi.fung.V3, inver=TRUE),]

#only on host species or lower.
globi.hsp<-globi.living.sp[globi.living.sp$globi.fung.V42=="species"|globi.living.sp$globi.fung.V42=="subspecies"|globi.living.sp$globi.fung.V42=="variety",]

saveRDS(globi.hsp, "intermediate_rds_csv/Rashmi_and_Globi/2024-12-21.UPD.globi.hsp.rds")

#compile only Mycosphere and Globi for Database.
Host.sp <- c(globi.hsp$globi.fung.V41, myc.as4$Host)
Fungi.sp<-c(globi.hsp$globi.fung.V3, myc.as4$Fungi)
DB<-c(globi.hsp$db, myc.as4$db) 
habitat.and.plant.part.1<-c(globi.hsp$globi.fung.V70, myc.as4$Host.part)
habitat.and.plant.part.2<-c(globi.hsp$globi.fung.V37, myc.as4$hab2)
habitat.and.plant.part.3<-c(globi.hsp$globi.fung.V68, myc.as4$hab2)
Location<-c(globi.hsp$db, myc.as4$Location)
Ref<-c(globi.hsp$globi.fung.V2, myc.as4$db)
globi.hsp$Fungi.auth<-""
Fungi.auth<-c(globi.hsp$Fungi.auth, myc.as4$Fungi.auth)

#combine into usable dataframe
sG<-cbind.data.frame(Host.sp, Fungi.sp, DB, habitat.and.plant.part.1, habitat.and.plant.part.2, habitat.and.plant.part.3, Location, Ref, Fungi.auth)
#tibble
sph.Globi.tibble<-sG %>%
  arrange(desc(DB)) %>%
  group_by(Host.sp,Fungi.sp) %>%
  distinct(Fungi.sp, .keep_all = TRUE)
sph.Globi.df<-data.frame(sph.Globi.tibble)
saveRDS(sph.Globi.df, "intermediate_rds_csv/Rashmi_and_Globi/2024-12-21.UPD.sph.Globi.df.UPD.4.2024.rds")

#filter out nonplant hosts - ncbi classification with new match.name2 column (taxizedb package).
sgd.h<-unique(sph.Globi.df$Host.sp)
#Some species names were causing glitches, I selected them out.
class2<-classification(sgd.h[-c(4885,13001,5585)], db="ncbi")

np<-sgd.h[setdiff(grep("cellular organisms", class2), grep("Viridiplantae", class2))]
#get non-plant genera

sph.Globi.df2<-sph.Globi.df[which(!sph.Globi.df$Host.sp %in% np),]

#read in lichenlist data from "Consortium of Lichen Herbaria" https://lichenportal.org/portal/checklists/checklist.php?clid=1492&pid=558
lichenlist<-read.csv("ext_data/World Checklist of Genera of Lichenized Fungi_1698006802.csv")
#make genus column
lichenlist$genus<-gsub(" sp.", "", lichenlist$ScientificName)
#filter lichenizing fungi
sph.Globi.df2$fung.gen<-word(sph.Globi.df2$Fungi.sp, 1)


sph.Globi.df3<-sph.Globi.df2[which(!sph.Globi.df2$fung.gen %in% lichenlist$genus),]
#check if fungal taxon entries that contain word root "lichen"
nlaf<-sph.Globi.df3[grep("lichen",sph.Globi.df3$Host.sp, ignore.case=TRUE),]
#no
#random typo I caught that caused issues.
sph.Globi.df3$Host.sp<-gsub("Abies alna", "Abies alba", sph.Globi.df3$Host.sp)
sph.Globi.df3$Host.sp.o<-gsub("subsp. |var. ", "", sph.Globi.df3$Host.sp)

#library U.Taxonstand
uni.corr.plant<-unique(sph.Globi.df3$Host.sp.o)
corr.plant<-nameMatch_WCVP(uni.corr.plant)
#remove the ones that still have no match. these are virtually all nonplants.
corr.plant.f<-corr.plant[which(!is.na(corr.plant$Accepted_SPNAME)),]



sgd.t<-sph.Globi.df3[which(sph.Globi.df3$Host.sp.o %in% corr.plant.f$Submitted_Name),]
sgd.t$Genus_in_database<-corr.plant.f$Genus_in_database[match(sgd.t$Host.sp, corr.plant.f$Submitted_Name)]
sgd.t$Name_in_database<-corr.plant.f$Name_in_database[match(sgd.t$Host.sp, corr.plant.f$Submitted_Name)]
sgd.t$Author_in_database<-corr.plant.f$Author_in_database[match(sgd.t$Host.sp, corr.plant.f$Submitted_Name)]
sgd.t$Accepted_SPNAME<-corr.plant.f$Accepted_SPNAME[match(sgd.t$Host.sp, corr.plant.f$Submitted_Name)]
sgd.t$Accepted_SPAUTHOR<-sgd.t$Author_in_database
sgd.t$Accepted_SPAUTHOR[which(!is.na(sgd.t$New_author))]<-sgd.t$New_author[which(!is.na(sgd.t$New_author))]

saveRDS(sgd.t, "intermediate_rds_csv/Rashmi_and_Globi/sgd.t.rds")

#keep only distinct combos with the new fungal species data
sph.Globi.df<-sgd.t %>%
  group_by(Accepted_SPNAME,Fungi.sp) %>%
  distinct(Accepted_SPNAME, .keep_all = TRUE)

#Still some fungal genera to get rid of
sph.Globi.df2<-sph.Globi.df[-which(str_count(sph.Globi.df$Fungi.sp, "\\S+") == "1"),]
sph.Globi.df3<-sph.Globi.df2[-grep(" sp\\.", sph.Globi.df2$Fungi.sp, ignore.case = TRUE),]


sph.Globi.df3$accepted.Fungi.sp.full.name<-paste(sph.Globi.df3$Fungi.sp, sph.Globi.df3$Fungi.auth)
sph.Globi.df3$accepted.Fungi.sp<-word(sph.Globi.df3$Fungi.sp, 1, 2, sep=" ")
#Make temp new column for leftover fungal auth/ssp. info that did not get formatted correctly
sph.Globi.df3$fung.sspandauth<-word(sph.Globi.df3$Fungi.sp, 3, -1, sep=" ")

#Manually reformat errant fungi auth.
write.csv(sph.Globi.df3, "intermediate_rds_csv/Rashmi_and_Globi/sph.Globi.df3.csv")
#fixed fungi auth and ssp., deleted "Fungi incertae sedis"
sph.Globi.df3.filt<-read.csv("intermediate_rds_csv/Rashmi_and_Globi/sph.Globi.df3.filt.csv")
sph.Globi.df4<-sph.Globi.df3.filt[,-1]

#there are more trees now that plant species names have been corrected. Filter Globi hits with no accompanying tissue detail, or on woody tissue. Filter phellophytes in Rashmi data (bark)
sph.Globi.df4tree<-sph.Globi.df4[which(word(as.character(sph.Globi.df4$Accepted_SPNAME), 1) %in% globtree$genus),]
sph.Globi.df4treeglobi<-sph.Globi.df4tree[which(sph.Globi.df4tree$DB=="Globi.db"),]
sph.Globi.df4treeglobi.keep<-sph.Globi.df4treeglobi[grep("lea|sor|petiole|sheath|anther|fruit|calyx|internode|first year|midrib|cone|peduncle|capitul|acervulus|bud|berry|flower|culm|pod|ovary|appressorium|green parts|spikelet|pycnium|aecium|telium|sepal|pedicel|catkin|flor|fruto|hoja|semilla|pycnothyrium|frond|caespituli|needle|seed|stolon|twig", sph.Globi.df4treeglobi$habitat.and.plant.part.1, ignore.case=TRUE),]

sph.Globi.df4notree<-sph.Globi.df4[-which(word(as.character(sph.Globi.df4$Accepted_SPNAME), 1) %in% globtree$genus),]
sph.Globi.df4treeKushveer<-sph.Globi.df4tree[which(sph.Globi.df4tree$DB=="Rashmi.Kushveer.Sarma.19"),]
#Combine the ones we will keep
sph.Globi.df5<-rbind.data.frame(sph.Globi.df4treeKushveer, sph.Globi.df4notree, sph.Globi.df4treeglobi.keep)
#remove bark from all
sph.Globi.df6<-sph.Globi.df5[-which(sph.Globi.df5$habitat.and.plant.part.1=="Bark"),]

# some of the Globi records for "hasHost" are not on living tissue and are the same records that were filtered out (e.g., Mycoportal)
# Many of these are probably on living tissue, but hard to say. Remove them.
sph.Globi.df7<-sph.Globi.df6[-which(sph.Globi.df6$habitat.and.plant.part.2 == "hasHost"),]

corrected.Host.sp<-sph.Globi.df7$Accepted_SPNAME
corrected.Host.ssp<-word(sph.Globi.df7$Accepted_SPNAME, 3, -1, sep=" ")
corrected.Host.auth<-sph.Globi.df7$Accepted_SPAUTHOR
detected.Host.plant<-sph.Globi.df7$Host.sp
accepted.Fungi.sp<-sph.Globi.df7$accepted.Fungi.sp
accepted.Fungi.sp.full.name<-sph.Globi.df7$accepted.Fungi.sp.full.name
verbatim.Fungi.sp<-sph.Globi.df7$Fungi.sp
habitat.and.plant.part.1<-sph.Globi.df7$habitat.and.plant.part.1
habitat.and.plant.part.2<-sph.Globi.df7$habitat.and.plant.part.2
habitat.and.plant.part.3<-sph.Globi.df7$habitat.and.plant.part.3
DB<-sph.Globi.df7$DB
Location<-sph.Globi.df7$Location
Ref<-sph.Globi.df7$Ref
accepted.Fungi.auth<-sph.Globi.df7$Fungi.auth

#make final dataframe with both corrected and uncorrected host names in separate columns.
FIN.sph.Globi<-cbind.data.frame(corrected.Host.sp, corrected.Host.ssp, corrected.Host.auth, 
                              detected.Host.plant, accepted.Fungi.sp, accepted.Fungi.sp.full.name, accepted.Fungi.auth, verbatim.Fungi.sp, 
                              habitat.and.plant.part.1, habitat.and.plant.part.2, habitat.and.plant.part.3, Location, Ref, DB)



#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
# III. Here we create FIN.sph.Globi, which will be used as raw data for the "On living plant tissue"  database.
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

saveRDS(FIN.sph.Globi, "Cleaned_data/FIN.sph.Globi.rds")

#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
