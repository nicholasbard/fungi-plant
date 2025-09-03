# Taxonomical Analysis
setwd("~/OneDrive - UBC/DB_fung_2025.01/")

library(countrycode)
library(dplyr)
library(fungarium)
library(jsonlite)
library(stringr)
library(data.table)
library(V.PhyloMaker2)
library(ape)
library(ggplot2)
library(plantlist)
library(forcats)
library(U.Taxonstand)
library(taxizedb)
library(egg)
library(FactoMineR)
library(pals)
library(rWCVP)


#~#~#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
### Which are detected in by Funguild DB? ######
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

Living.plant<-readRDS("Final_comb_DBs/LIVING.standardized.cleaned.no.error.uniq.rds")
#length(unique(Living.plant$new_species)) #13114 fungal species
#length(unique(Living.plant$new_genus)) #2501 fungal species
#length(unique(Living.plant$corrected.Host.sp)) #10945 plant species

#Read in Funguild data (downloaded from website 2025, Jan 20)
funguild <- fromJSON("ext_data/Funguild.data.txt")
funguild$genus<-word(funguild$taxon, 1)

# Which species are in Funguild data?
Funguild.living.sp<-inner_join(Living.plant, funguild, by = c("new_species"="taxon"))
saveRDS(Funguild.living.sp, "Analyses/Results/05-05.Funguild.living.sp.rds")


# Which species are not in Funguild data?
Sps.not.in.Funguild<-setdiff(Living.plant, Funguild.living.sp[,-c(32:42)])
#length(unique(Sps.not.in.Funguild$new_species)) # 12085 fungal species not present in funguild (of 13121)

#Which genera are in Funguild data? Note: use distinct guilds
funguild.dist.genera<-funguild %>%
  arrange(taxon) %>%
  distinct(genus, .keep_all = TRUE)

Funguild.living.gen<-inner_join(Living.plant, funguild.dist.genera, by = c("new_genus"="genus"))

saveRDS(Funguild.living.gen, "Analyses/Results/05-05.Funguild.living.gen.rds")

#which genera are not in Funguild data?
Gens.not.in.Funguild<-setdiff(Living.plant, Funguild.living.gen[,-c(32:42)])
#length(unique(Gens.not.in.Funguild$new_genus)) # 300 fungal genera not present in funguild (of 2505)

#~#~#~#~#~#~#~##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
### Make a phylogenetic heatmap ######
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


#### 1. use standardized fungi


#B. get phylo info for species.

# make a phylogeny for plants

#THIS MAKES a host phylogeny.
# Read in host phylogenetic info
# read in megaphylogeny of host plants 
megatree <- GBOTB.extended.TPL

# read in fungi-host data
Living.plant$Host.sp.phylo<-gsub(" ", "_", Living.plant$corrected.Host.sp)

megatree <- drop.tip(megatree, megatree$tip.label[!megatree$tip.label%in%unique(Living.plant$Host.sp.phylo)])

FACUS.stand.in.meg<-Living.plant[Living.plant$Host.sp.phylo %in% megatree$tip.label,]
FACUS.not<-setdiff(Living.plant, FACUS.stand.in.meg)

#build species list for tips
gen.fam.info<-nameMatch_WCVP(FACUS.not$corrected.Host.sp)
species.list.FACUS.not<-data.frame()
species.list.FACUS.not<-cbind.data.frame("species"=FACUS.not$Host.sp.phylo, "genus"=gen.fam.info$Genus_in_database, "family"=gen.fam.info$Family)


### 
gen.fam.info2<-nameMatch_WCVP(FACUS.stand.in.meg$corrected.Host.sp)
species.list.FACUS.stand.in.meg<-cbind.data.frame("species"=FACUS.stand.in.meg$Host.sp.phylo, "genus"=gen.fam.info2$Genus_in_database, "family"=gen.fam.info2$Family)
species.list.FACUS.stand.in.meg.dist<-species.list.FACUS.stand.in.meg %>% distinct(species, .keep_all = TRUE)
rrr<-rbind.data.frame(species.list.FACUS.stand.in.meg.dist, species.list.FACUS.not[,1:3])
Host.phylo.tree<-phylo.maker(rrr)
saveRDS(Host.phylo.tree, "Analyses/Results/05-05.Host.phylo.tree.rds")
####

#make phylogeny for fungi
FACUSF.stand.f<-Living.plant %>%
  group_by(new_species) %>%
  distinct(new_kingdom, new_phylum, new_class, new_order, new_family, new_genus, new_species, .keep_all = TRUE) 
FACUSF.stand.f.filt<-FACUSF.stand.f %>%
  distinct(new_kingdom, new_phylum, new_class,new_order,new_family,new_genus, new_species) 

FACUSF.stand.f.filt$new_family[which(is.na(FACUSF.stand.f.filt$new_family))]<-""
FACUSF.stand.f.filt$new_order[which(is.na(FACUSF.stand.f.filt$new_order))]<-""
FACUSF.stand.f.filt$new_species<-as.factor(FACUSF.stand.f.filt$new_species)
FACUSF.stand.f.filt$new_genus<-as.factor(FACUSF.stand.f.filt$new_genus)
FACUSF.stand.f.filt$new_family<-as.factor(FACUSF.stand.f.filt$new_family)
FACUSF.stand.f.filt$new_order<-as.factor(FACUSF.stand.f.filt$new_order)
FACUSF.stand.f.filt$new_class<-as.factor(FACUSF.stand.f.filt$new_class)
FACUSF.stand.f.filt$new_phylum<-as.factor(FACUSF.stand.f.filt$new_phylum)
FACUSF.stand.f.filt$new_kingdom<-as.factor(FACUSF.stand.f.filt$new_kingdom)

FACUSF.tree <- as.phylo(~new_kingdom/new_phylum/new_class/new_order/new_family/new_genus/new_species, data =FACUSF.stand.f.filt)


#A. create matrix with fungal inf per host sp

FACUSF.stand.adj<-Living.plant[Living.plant$new_species %in% FACUSF.stand.f.filt$new_species, ]
FACUSF.stand.adj$Host.sp.phylo<-gsub(" ", "_", FACUSF.stand.adj$corrected.Host.sp)
FACUSF.stand.adj2<-FACUSF.stand.adj[FACUSF.stand.adj$Host.sp.phylo %in% Host.phylo.tree$scenario.3$tip.label, ]

FACUSF.stand.adj2$host.gen<-word(FACUSF.stand.adj2$corrected.Host.sp, 1)
saveRDS(FACUSF.stand.adj2, "Analyses/Results/03-03.FACUSF.stand.adj2.rds")

freq.tab<-table(FACUSF.stand.adj2[,c(26,32)])
#saveRDS(freq.tab, "Analyses/int_files/05-05.freq.tab.rds")

Host.plant.tree2 <- drop.tip(Host.phylo.tree$scenario.3, Host.phylo.tree$scenario.3$tip.label[!Host.phylo.tree$scenario.3$tip.label%in%colnames(freq.tab)])
FACUSF.tree2 <-drop.tip(FACUSF.tree, FACUSF.tree$tip.label[!FACUSF.tree$tip.label %in% rownames(freq.tab)])

write.tree(FACUSF.tree2, "Analyses/int_files/05-05.facus.fungal2.tre")
write.tree(Host.plant.tree2, "Analyses/int_files/05-05.Host.plant.tree2.tre")




F.s2.abbrev<-cbind.data.frame("host"=FACUSF.stand.adj2$corrected.Host.sp, "fungi"=FACUSF.stand.adj2$new_species, "value"=1)
Fs2ab.mat<-xtabs(value ~ ., F.s2.abbrev)
rownames(Fs2ab.mat)<-gsub(" ","_", rownames(Fs2ab.mat))

Fs2ab.m_ordered <- Fs2ab.mat[match(Host.plant.tree2$tip.label, rownames(Fs2ab.mat)),match(FACUSF.tree2$tip.label, colnames(Fs2ab.mat))]
Fs2ab.m_ord.df<-as.data.frame(Fs2ab.m_ordered)
Fs2ab.m_ord.df.filt0<-Fs2ab.m_ord.df[-which(Fs2ab.m_ord.df$Freq=="0"),]
Fs2ab.m_ord.df.filt0$one<-"1"

saveRDS(Fs2ab.m_ord.df.filt0, "Analyses/int_files/05-05.Fs2ab.m_ord.df.filt0.rds")

# now retrieve the higher taxonomic info. make order column for fungi
Fs2ab.m_ord.df.filt0$fung.order<-FACUSF.stand.adj2$new_order[match(gsub("_"," ", Fs2ab.m_ord.df.filt0$fungi), FACUSF.stand.adj2$new_species)]
#Add in class when no order provided
Fs2ab.m_ord.df.filt0$fung.order[is.na(Fs2ab.m_ord.df.filt0$fung.order)]<-FACUSF.stand.adj2$new_class[match(gsub("_"," ", Fs2ab.m_ord.df.filt0$fungi)[is.na(Fs2ab.m_ord.df.filt0$fung.order)], FACUSF.stand.adj2$new_species)]
#Add in phylum to ORDER when no class provided
Fs2ab.m_ord.df.filt0$fung.order[which(Fs2ab.m_ord.df.filt0$fung.order=="")]<-FACUSF.stand.adj2$new_phylum[match(gsub("_"," ", Fs2ab.m_ord.df.filt0$fungi)[which(Fs2ab.m_ord.df.filt0$fung.order=="")], FACUSF.stand.adj2$new_species)]
#not sure why this one didnt update, still NA
Fs2ab.m_ord.df.filt0$fung.order[is.na(Fs2ab.m_ord.df.filt0$fung.order)]<-"Ascomycota incertae sedis"
Fs2ab.m_ord.df.filt0$fung.order[Fs2ab.m_ord.df.filt0$fung.order=="Ascomycota"]<-"Ascomycota incertae sedis"
Fs2ab.m_ord.df.filt0$fung.order[Fs2ab.m_ord.df.filt0$fung.order=="Dothideomycetes"]<-"Dothideomycetes incertae sedis"
Fs2ab.m_ord.df.filt0$fung.order[Fs2ab.m_ord.df.filt0$fung.order=="Lecanoromycetes"]<-"Lecanoromycetes incertae sedis"
Fs2ab.m_ord.df.filt0$fung.order[Fs2ab.m_ord.df.filt0$fung.order=="Sordariomycetes"]<-"Sordariomycetes incertae sedis"
Fs2ab.m_ord.df.filt0$fung.order[Fs2ab.m_ord.df.filt0$fung.order=="Leotiomycetes"]<-"Leotiomycetes incertae sedis"
Fs2ab.m_ord.df.filt0$fung.order[Fs2ab.m_ord.df.filt0$fung.order=="Microbotryomycetes"]<-"Microbotryomycetes incertae sedis"
Fs2ab.m_ord.df.filt0$fung.order[Fs2ab.m_ord.df.filt0$fung.order=="Cystobasidiomycetes"]<-"Cystobasidiomycetes incertae sedis"


#make class column for fungi
Fs2ab.m_ord.df.filt0$fung.class<-FACUSF.stand.adj2$new_class[match(gsub("_"," ", Fs2ab.m_ord.df.filt0$fungi), FACUSF.stand.adj2$new_species)]
#Add in phylum (Ascomycota) when no class provided
Fs2ab.m_ord.df.filt0$fung.class[which(Fs2ab.m_ord.df.filt0$fung.class=="")]<-FACUSF.stand.adj2$new_phylum[match(gsub("_"," ", Fs2ab.m_ord.df.filt0$fungi)[which(Fs2ab.m_ord.df.filt0$fung.class=="")], FACUSF.stand.adj2$new_species)] #Change label for these no-class entries
Fs2ab.m_ord.df.filt0$fung.class[which(Fs2ab.m_ord.df.filt0$fung.class=="Ascomycota")]<-"Ascomycota incertae sedis"
#one instance of Phloeosporina fraxini, which seems undefined.
Fs2ab.m_ord.df.filt0$fung.class[which(is.na(Fs2ab.m_ord.df.filt0$fung.class))]<-"Ascomycota incertae sedis"

#make order column for plants

#get family info
# get the gen.fam.info from earlier
#gen.fam.info2<-nameMatch_WCVP(FACUS.stand.in.meg$corrected.Host.sp)
gen.fam.info.complete<-rbind.data.frame(gen.fam.info, gen.fam.info2)
saveRDS(gen.fam.info.complete, "Analyses/int_files/05-05.gen.fam.info.complete.rds")

Fs2ab.m_ord.df.filt0$host.genera<-word(gsub("_", " ", Fs2ab.m_ord.df.filt0$host), 1)
Fs2ab.m_ord.df.filt0$host.families<-gen.fam.info.complete$Family[match(Fs2ab.m_ord.df.filt0$host.genera, gen.fam.info.complete$Genus_in_database)]
#use "plantlist" orders_dat to find orders
Fs2ab.m_ord.df.filt0$host.orders<-orders_dat$ORDER[match(Fs2ab.m_ord.df.filt0$host.families, orders_dat$FAMILY)]

#create new df
Fs4<-Fs2ab.m_ord.df.filt0
Fs4$host.orders[Fs4$host.families=="Icacinaceae"]<-"Icacinales"
Fs4$host.orders[Fs4$host.families=="Metteniusaceae"]<-"Metteniusales"
Fs4$host.orders[Fs4$host=="×_Elyhordeum_macounii"]<-"Poales"
Fs4$host.orders[Fs4$host=="×_Fatshedera_lizei"]<-"Apiales"
Fs4$host.orders[Fs4$host=="×_Comagaria_rosea"]<-"Rosales"
Fs4$host.orders[Fs4$host=="+_Laburnocytisus_adami"]<-"Fabales"
Fs4$host.orders[Fs4$host=="Mazus_pumilus"]<-"Lamiales"
Fs4$host.orders[Fs4$host=="×_Pyronia_veitchii"]<-"Rosales"
Fs4$host.orders[Fs4$host=="×_Elyleymus_aristatus"]<-"Poales"
Fs4$host.orders[Fs4$host=="Mazus_miquelii"]<-"Lamiales"
Fs4$host.orders[Fs4$host=="x_Agropogon_lutosus"]<-"Poales"
Fs4$host.orders[Fs4$host=="×_Agropogon_lutosus"]<-"Poales"
Fs4$host.orders[Fs4$host.orders==" Selaginellales"]<-"Selaginellales"
#will need to assign this to different order 
Fs4$host.orders[Fs4$host.orders=="Sabiales"]<-"Proteales"

ord.info <- classification(unique(Fs4$host.orders, db = "ncbi"))

ord.info2<-data.frame()
for (i in 1:length(ord.info)) {
  ord.info2 <- rbind(ord.info2, ord.info[[i]][13,])
}
ord.info3<-cbind.data.frame(host.ord=unique(Fs4$host.orders), ord.info2)

#Change order using factor with plant phylogeny
Fs4$host.factor <- factor(Fs4$host, levels = Host.plant.tree2$tip.label)
Fs4$fungi.factor<-factor(Fs4$fungi, levels = FACUSF.tree2$tip.label)

#
Fs5.sp<- Fs4 %>%
  arrange(host.factor)
saveRDS(Fs5.sp, "Analyses/int_files/05-05.UNEDITED.Fs5.sp.rds")


#Since the scales are not great, we can make new columns for a better scale
Fs5.sp$plant.hi.t<-Fs5.sp$host.orders
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Acorales"]<-"other Liliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Alismatales"]<-"other Liliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Aquifoliales"]<-"other asterids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Araucariales"]<-"gymnosperms"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Arecales"]<-"other commelinids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Austrobaileyales"]<-"other Magnoliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Buxales"]<-"other Magnoliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Canellales"]<-"other Magnoliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Celastrales"]<-"other rosids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Commelinales"]<-"other commelinids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Cornales"]<-"other asterids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Crossosomatales"]<-"other rosids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Cucurbitales"]<-"other rosids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Cyatheales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Cycadales"]<-"gymnosperms"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Dilleniales"]<-"other Gunneridae"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Dioscoreales"]<-"other Liliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Equisetales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Escalloniales"]<-"other asterids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Garryales"]<-"other asterids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Geraniales"]<-"other rosids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Ginkgoales"]<-"gymnosperms"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Gleicheniales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Gnetales"]<-"gymnosperms"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Gunnerales"]<-"other Gunneridae"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Icacinales"]<-"other asterids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Laurales"]<-"other Magnoliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Lycopodiales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Magnoliales"]<-"other Magnoliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Marattiales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Metteniusales"]<-"other asterids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Nymphaeales"]<-"other Magnoliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Ophioglossales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Osmundales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Oxalidales"]<-"other rosids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Pandanales"]<-"other Liliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Piperales"]<-"other Magnoliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Proteales"]<-"other Magnoliopsida"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Psilotales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Santalales"]<-"other Gunneridae"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Vitales"]<-"other Gunneridae"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Zingiberales"]<-"other commelinids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Zygophyllales"]<-"other rosids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Cupressales"]<-"other rosids"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Selaginellales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Schizaeales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Salviniales"]<-"other tracheophytes"
Fs5.sp$plant.hi.t[Fs5.sp$host.orders=="Pinales"]<-"gymnosperms"

fung.phylum<-FACUSF.stand.adj2$new_phylum[match(gsub("_"," ", Fs5.sp$fungi), FACUSF.stand.adj2$new_species)]
Fs5.sp$fung.hi.t<-fung.phylum
Fs5.sp$fung.hi.t[Fs5.sp$fung.class=="Sordariomycetes"]<-"other Sordariomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.class=="Leotiomycetes"]<-"Leotiomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.class=="Agaricomycetes"]<-"other Agaricomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.class=="Dothideomycetes"]<-"other Dothideomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.class=="Eurotiomycetes"]<-"Eurotiomycetes"

Fs5.sp$fung.hi.t[Fs5.sp$fung.class=="Ustilaginomycetes"]<-"Ustilaginomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.class=="Pucciniomycetes"]<-"Pucciniomycetes"

Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Agaricales"]<-"Agaricales - Agaricomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Polyporales"]<-"Polyporales - Agaricomycetes"
#Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Pucciniales"]<-"Pucciniales - Pucciniomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Hypocreales"]<-"Hypocreales - Sordariomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Xylariales"]<-"Xylariales - Sordariomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Diaporthales"]<-"Diaporthales - Sordariomycetes"
#Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Eurotiales"]<-"Eurotiales - Sordariomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Mycosphaerellales"]<-"Mycosphaerellales - Dothideomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Pleosporales"]<-"Pleosporales - Dothideomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Botryosphaeriales"]<-"Botryosphaeriales - Dothideomycetes"
Fs5.sp$fung.hi.t[Fs5.sp$fung.order=="Capnodiales"]<-"Capnodiales - Dothideomycetes"



Fs5.sp$fung.hi.t[Fs5.sp$fung.hi.t=="Blastocladiomycota"]<-"other phyla"
Fs5.sp$fung.hi.t[Fs5.sp$fung.hi.t=="Chytridiomycota"]<-"other phyla"
Fs5.sp$fung.hi.t[Fs5.sp$fung.hi.t=="Entomophthoromycota"]<-"other phyla"
Fs5.sp$fung.hi.t[Fs5.sp$fung.hi.t=="Mucoromycota"]<-"other phyla"
Fs5.sp$fung.hi.t[Fs5.sp$fung.hi.t=="Glomeromycota"]<-"other phyla"
Fs5.sp$fung.hi.t[Fs5.sp$fung.hi.t=="Ascomycota"]<-"other Ascomycota"
Fs5.sp$fung.hi.t[Fs5.sp$fung.hi.t=="Basidiomycota"]<-"other Basidiomycota"

Fs5.sp$fung.hi.t <- fct_relevel(Fs5.sp$fung.hi.t, "Agaricales - Agaricomycetes", "Polyporales - Agaricomycetes",
                                "other Agaricomycetes", "Ustilaginomycetes", "Pucciniomycetes",
                                "other Basidiomycota", "Hypocreales - Sordariomycetes",
                                "Xylariales - Sordariomycetes", "Diaporthales - Sordariomycetes",
                                "other Sordariomycetes", "Leotiomycetes", "Eurotiomycetes",
                                "Mycosphaerellales - Dothideomycetes","Pleosporales - Dothideomycetes",
                                "Botryosphaeriales - Dothideomycetes","Capnodiales - Dothideomycetes",
                                "other Dothideomycetes", "other Ascomycota", "other phyla")
saveRDS(Fs5.sp, "Analyses/int_files/13-07.Fs5.sp.rds")
#here is the co-occurrence heat map
ggplot(data.frame(Fs5.sp), aes(host, fungi, fill=one)) + 
  geom_tile()+
  scale_fill_discrete(type="black") +
  scale_y_discrete(drop = TRUE, expand = c(0, 0))+
  facet_grid(fung.hi.t~fct_inorder(plant.hi.t), scales="free", space = "free", switch="y")+
  theme(strip.placement = "outside",panel.spacing = unit(0, "in"),strip.background.y = element_rect(fill = "red", color = "gray35"), strip.text.y = element_text(size=8), axis.text.y = element_blank())+
  theme(strip.background.x = element_rect(fill = "white", color = "gray35"),axis.text.x = element_blank())+
  theme(strip.text.y.left = element_text(angle = 0, size =7), strip.text.x.top = element_text(angle = 90, size=7), axis.ticks=element_blank())+
  scale_x_discrete(position="top")+
  guides(fill="none")

#with color, same plot
library(ggh4x)

ggplot(Fs5.sp, aes(host, fungi, fill = one)) +
  geom_tile() +
  scale_fill_discrete(type = "black") +
  scale_y_discrete(drop = TRUE, expand = c(0, 0)) +
  facet_grid2(
    fung.hi.t ~ fct_inorder(plant.hi.t),
    scales = "free",
    space = "free",
    switch = "y",
    strip = strip_themed(
      background_y = list(
        element_rect(fill = alpha("#5A5156",0.3), color = "gray35"),  # agaricale
        element_rect(fill = alpha("#AA0DFE",0.3), color = "gray35"), # polyporale
        element_rect(fill = alpha("#332288",0.3), color = "gray35"), # agaricomyce
        element_rect(fill = alpha("#AA4499",0.3), color = "gray35"), # ustilagino
        element_rect(fill = alpha("#CC6677",0.3), color = "gray35"), # pucciniomycetes
        element_rect(fill = alpha("#0075DC",0.3), color = "gray35"), # basidiomycota
        element_rect(fill = alpha("#1CFFCE",0.3), color = "gray35"), # hypocreales
        element_rect(fill = alpha("#C4451C",0.3), color = "gray35"), # xylariales
        element_rect(fill = alpha("#16FF32",0.3), color = "gray35"), # diaporthales
        element_rect(fill = alpha("#CC6677",0.3), color = "gray35"), # sordariomycetes
        element_rect(fill = alpha("#999933",0.3), color = "gray35"), # leotiomycetes
        element_rect(fill = alpha("#44AA99",0.3), color = "gray35"), # eurotiomycetes
        element_rect(fill = alpha("#90AD1C",0.3), color = "gray35"), # mycosphaerales
        element_rect(fill = alpha("#DEA0FD",0.3), color = "gray35"), # pleoporales
        element_rect(fill = alpha("#F6222E",0.3), color = "gray35"), # boytrosphaeliales
        element_rect(fill = alpha("#FE00FA",0.3), color = "gray35"), # capnodiales
        element_rect(fill = alpha("#88CCEE",0.3), color = "gray35"), # dothideomyctes
        element_rect(fill = alpha("#F0A0FF",0.3), color = "gray35"), # ascomycota
        element_rect(fill = "white", color = "gray35")
      )
    )
  ) +
  theme(
    panel.spacing = unit(0, "in"),
    strip.text.y = element_text(size = 8),
    axis.text.y = element_blank(),
    strip.background.x = element_rect(fill = "white", color = "gray35"),
    axis.text.x = element_blank(),
    strip.text.y.left = element_text(angle = 0, size = 7),
    strip.text.x.top = element_text(angle = 90, size = 7),
    axis.ticks = element_blank()
  ) +
  scale_x_discrete(position = "top") +
  guides(fill = "none")

#saved as "05-05.host.matrix.pdf"
#use human use data
human<-read.csv("Analyses/data/human.use.dataset/raw_data/utilised_plants_species_list.csv")
human$binomial_acc_name<-gsub(" ", "_", human$binomial_acc_name)
nn<-human[which(human$binomial_acc_name %in% Fs5.sp$host),]
mm<-unique(as.character(Fs5.sp$host[which(!Fs5.sp$host %in% human$binomial_acc_name)]))

#add in a bunch of 0s with the rest of the species not in the database.
pp <- data.frame(matrix(ncol = length(colnames(human)), nrow = length(mm)))
colnames(pp)<-colnames(human)
pp$binomial_acc_name<-mm
pp[is.na(pp)]<-0

#This is the complete species dataframe
human.and.non<-rbind.data.frame(nn,pp)




#or we can just add the plant use data alongside the plant column
Fs5.sp$host<-as.character(Fs5.sp$host)
Fs6<-left_join(Fs5.sp, human.and.non, by=c("host"= "binomial_acc_name"))
saveRDS(Fs6, "Analyses/int_files/05-05.Fs6.rds")

#make new df structured for use, instead of plant
Fs6.use<-rbind.data.frame(Fs6[which(Fs6$AnimalFood==1),],Fs6[which(Fs6$EnvironmentalUses==1),],Fs6[which(Fs6$Fuels==1),],Fs6[which(Fs6$GeneSources==1),],Fs6[which(Fs6$HumanFood==1),],
                          Fs6[which(Fs6$InvertebrateFood==1),],Fs6[which(Fs6$Materials==1),],Fs6[which(Fs6$Medicines==1),],Fs6[which(Fs6$Poisons==1),],Fs6[which(Fs6$SocialUses==1),],Fs6[which(Fs6$Totals==0),])

use<-c(rep("AnimalFood", each=length(which(Fs6$AnimalFood==1))),rep("EnvironmentalUses", each=length(which(Fs6$EnvironmentalUses==1))),rep("Fuels", each=length(which(Fs6$Fuels==1))),
       rep("GeneSources", each=length(which(Fs6$GeneSources==1))),rep("HumanFood", each=length(which(Fs6$HumanFood==1))),rep("InvertebrateFood", each=length(which(Fs6$InvertebrateFood==1))),
       rep("Materials", each=length(which(Fs6$Materials==1))),rep("Medicines", each=length(which(Fs6$Medicines==1))),rep("Poisons", each=length(which(Fs6$Poisons==1))),
       rep("SocialUses", each=length(which(Fs6$SocialUses==1))),rep("Unused or unknown", each=length(which(Fs6$Totals==0))))
            
Fs6.use$use<-use
Fs6.use$fungi<-as.character(Fs6.use$fungi)
saveRDS(Fs6.use, "Analyses/int_files/05-05.Fs6.use.rds")

Fs6.use2<-Fs6.use[!Fs6.use$use=="Unused or unknown",]


#Get totals
length(unique(Fs6.use2$host)) #5845
length(unique(Fs6.use2$fungi))#10878

Fs6.use.dist<-Fs6.use2 %>%
  group_by(fungi, host) %>%
  distinct(fungi, .keep_all=TRUE)

Fs6.use.dist$use <- fct_relevel(Fs6.use.dist$use, "AnimalFood","EnvironmentalUses","Fuels","GeneSources","HumanFood",
                                "InvertebrateFood","Materials","Medicines","Poisons","SocialUses")

sum(Fs6.use.dist$AnimalFood)
#8555
sum(Fs6.use.dist$EnvironmentalUses)
sum(Fs6.use.dist$Fuels)
sum(Fs6.use.dist$GeneSources)
sum(Fs6.use.dist$HumanFood)
sum(Fs6.use.dist$InvertebrateFood)
sum(Fs6.use.dist$Materials)
sum(Fs6.use.dist$Medicines)
sum(Fs6.use.dist$Poisons)
sum(Fs6.use.dist$SocialUses)

## To do Pcoa
Fs6.use2_freq<-data.frame(with(Fs6.use2, table(Fs6.use2$use, Fs6.use2$fungi)))
Fs6.use2_freq.t <- reshape2::dcast(Fs6.use2_freq, Fs6.use2_freq[,1] ~ Fs6.use2_freq[,2])
rownames(Fs6.use2_freq.t)<-Fs6.use2_freq.t[,1]
Fs6.use2_freq.t<-Fs6.use2_freq.t[,-1]

#convert to 0s and 1s (presence/absence of fungi in each type)
Fs6.use2_freq.t[!Fs6.use2_freq.t=="0"]<-1
#calculate dissimilarity
dissim.bray <- vegdist(Fs6.use2_freq.t, method = "bray", binary=TRUE) 
#PCOA Bray/Sorenson
pcoa.bray<-cmdscale(dissim.bray, eig = TRUE)
pcoa.bray$species<-wascores(pcoa.bray$points, Fs6.use2_freq.t, expand = TRUE)


rownames(pcoa.bray$points)[8]<-"Medi-\ncines"
rownames(pcoa.bray$points)[2]<-"Envi.\nuses"
dft.pl <- ordiplot(pcoa.bray, type="none", cex=1.3, cex.lab=1.5, cex.axis=1.5, xlab="Dim. 1 (45.54%)", ylab="Dim. 2 (15.38%)",xlim = c(-0.6,1.4), ylim=c(-0.4,0.5)) 
species_use <- as.data.frame(t(Fs6.use2_freq.t))
use_types <- colnames(species_use)
palette <- c("yellow", "lightblue", "violet", "brown3", "yellow4", "purple4", "tan", "grey", "pink1", "darkorange")

for (i in seq_along(use_types)) {
  use <- use_types[i]
  species_in_use <- rownames(species_use)[species_use[[use]] == 1]
  coords <- pcoa.bray$species[species_in_use, ]
  
hull <- chull(coords)
    polygon(coords[hull, ], col = adjustcolor(palette[i], alpha.f = 0.18), border = palette[i], lwd = 2)
}
points(dft.pl, "species", col="white", pch=16, xlim = c(-0.4,0.35), ylim=c(-0.4,0.5)) 
text(dft.pl, "sites", col="black", cex=1.1, cex.lab=1.5, cex.axis=1.5, xlim = c(-0.4,0.35), ylim=c(-0.4,0.5))
legend("right", legend = use_types, fill = adjustcolor(palette, alpha.f = 0.18), border = palette, cex = 1.1, 
       xpd=TRUE,  y.intersp = 0.8, x.intersp = 0.8,text.width = max(strwidth(use_types)) * 0.6,
       title = "Human use type", title.cex = 1.3, title.font=2 )
# save as approximately 1800 X 1200 pixels (prev. 3500 X 1700 pixels )
#saves as "Supp.Fig.1.PCOA.DB.png"

#### To do beta diversity analysis
library(betapart)
library(reshape2)
library(ggplot2)
library(vegan)
#partition species richness as nestedness & capture species composition for turnover
Fs6.beta<-betapart.core(Fs6.use2_freq.t)


Fs6.beta.m<-beta.multi(Fs6.use2_freq.t)
write.csv(Fs6.beta.m, "Analyses/Results/05-05.Fs6.beta.multi.csv")
Fs6.beta.p<-beta.pair(Fs6.use2_freq.t)
saveRDS(Fs6.beta.p, "Analyses/Results/05-05.Fs6.beta.pairwise.rds")
#create matrix from turnover
Fs6.b.turn.m<-as.matrix(Fs6.beta.p$beta.sim)
Fs6.b.turn.m[upper.tri(Fs6.b.turn.m)] <- NA
#create matrix from nestedness
Fs6.b.nest.m<-as.matrix(Fs6.beta.p$beta.sne)
Fs6.b.nest.m[lower.tri(Fs6.b.nest.m)] <- NA
write.csv(Fs6.b.turn.m, "Analyses/Results/05-05.Fs6.beta.turnover.csv")
write.csv(Fs6.b.nest.m, "Analyses/Results/05-05.Fs6.beta.nestedness.csv")

#combine both
Fs6.b.all<-Fs6.b.nest.m
Fs6.b.all[is.na(Fs6.b.all)] <- Fs6.b.turn.m[is.na(Fs6.b.all)]
Fs6.b.all.df<-melt(Fs6.b.all)
Fsg.b.all.a<-Fs6.b.all.df
colnames(Fs6.b.all.df)[3]<-"Sorenson index value"


#make plot of both sorenson indices
#top left is nestedness, bottom right is turnover
#first make comparisons to itself NA so it shows up as a different color scheme, then make plot
Fs6.b.all.df$`Sorenson index value`[Fs6.b.all.df$`Sorenson index value`==0]<-NA
ggplot(Fs6.b.all.df, aes(x = Var1, y = Var2, fill = `Sorenson index value`, direction=-1)) +
  geom_tile()+
  scale_fill_gradient(low="lightblue1", high="darkblue", na.value = "white")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(size = 10.5),axis.text.y=element_text(size = 10.5))

### end of beta diversity analysis

## chi-square to test significant difference between expected interactions and observed
plants=c(sum(human[,3]),sum(human[,4]),sum(human[,5]),sum(human[,6]),sum(human[,7]),sum(human[,8]),sum(human[,9]),sum(human[,10]),sum(human[,11]),sum(human[,12]))
inters=as.vector((table(Fs6.use2$use)))
ch.sq.t.d<-rbind(plants, inters)
colnames(ch.sq.t.d)<-colnames(Fs6.use2[,c(15:24)])
chi.res.hum<-chisq.test(ch.sq.t.d)
#Pearson's Chi-squared test
#data:  ch.sq.t.d
#X-squared = 4014.2, df = 9, p-value < 2.2e-16

chi.res.hum$residuals
#        AnimalFood EnvironmentalUses     Fuels GeneSources HumanFood InvertebrateFood Materials Medicines   Poisons SocialUses
#plants  -12.59131         -19.75506 -3.027729   -7.799455 -3.277991       -10.504148  18.36189  27.86672 -24.07598  1.0869072
#inters   10.14293          15.91369  2.438988    6.282852  2.640586         8.461619 -14.79143 -22.44804  19.39441 -0.8755584
#add proportion changes row
ch.sq.t.d<-rbind(ch.sq.t.d, (plants/sum(plants))/(inters/sum(inters)))
#pairwise post hoc (not sure this tells us anything, except that social uses, invertebrate food, and poisons proportion change compared to every other pairwsie proportion change is significant, while others are insignificant in one pairwise comparison only..)
library(fifer)
chisq.post.hoc(ch.sq.t.d, control = "bonferroni", popsInRows  = FALSE)
# find that all pairwise proportions are significant p<0.01 except for Animalfood/EnviUse; fuels/genesources; fuels/HumanFood; fuels/social uses (p=~0.0141) humanfood/SocialUses (p=~0.0434); InvertebrateFood/Poisons (p=~0.0195); Materials/Medicines
#There is a difference in number of interactions for each human use than what is expected for every category 


# get started accessing geographic data
#WCVP data downloaded 2025-02-15
wcvp.names<-read.csv("Analyses/data/wcvp/wcvp_names.csv", sep = "|")
wcvp.dist<-read.csv("Analyses/data/wcvp/wcvp_distribution.csv", sep = "|")

#updated 2025-04-09
wcvp.names<-wcvp_names
saveRDS(wcvp.names, "ext_data/2025-04-09.wcvp.names.rds")
wcvp.dist<-wcvp_distributions
saveRDS(wcvp.dist, "2025-04-09.wcvp.dist.rds")

#to get total numbers of lifeforms
#need to re-retrieve the ID numbers for host plants
Fs5.sp$host<-gsub("_", " ", Fs5.sp$host)
tm<-nameMatch_WCVP(Fs5.sp$host)
filt.wn<-wcvp.names[which(wcvp.names$plant_name_id %in% tm$ID_in_database),]

filt.wn.uniq<-filt.wn %>%
  arrange(lifeform_description) %>%
  distinct(taxon_name, .keep_all=TRUE)

#now we integrate our data
filt.wnd<-left_join(filt.wn.uniq, wcvp.dist, by="plant_name_id")
Fs5.sp$host<-as.character(gsub("_", " ", Fs5.sp$host))
Fs7<-left_join(Fs5.sp, filt.wnd, by=c("host"="taxon_name"), relationship = "many-to-many")
saveRDS(Fs7, "Analyses/int_files/Fs7.rds")

Fs7.reg.dist<-Fs7 %>%
  group_by(fungi, region) %>%
  distinct(fungi, .keep_all=TRUE)
Fs7.reg.dist$region[is.na(Fs7.reg.dist$region)]<-""
Fs7.reg.dist$region[Fs7.reg.dist$region==""]<-"Unlisted region"

Fs7.habit.dist<-Fs7 %>%
  group_by(fungi, lifeform_description) %>%
  distinct(fungi, .keep_all=TRUE)

Fs7.habit.dist.tree<-Fs7.habit.dist[grep("tree", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.shrub<-Fs7.habit.dist[grep("^shrub| shrub", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.subshrub<-Fs7.habit.dist[grep("subshrub", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.liana<-Fs7.habit.dist[grep("climb|liana", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.annual<-Fs7.habit.dist[grep("annual", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.biennial<-Fs7.habit.dist[grep("biennial", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.perennial<-Fs7.habit.dist[grep("perennial", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.epiphyte<-Fs7.habit.dist[grep("epiphyte", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.lithophyte<-Fs7.habit.dist[grep("lithophyte", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.geophyte<-Fs7.habit.dist[grep("geophyte", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.helophyte<-Fs7.habit.dist[grep("helophyte", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.parasite<-Fs7.habit.dist[grep("parasit", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.succulent<-Fs7.habit.dist[grep("succulent", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.bamboo<-Fs7.habit.dist[grep("bamboo", Fs7.habit.dist$lifeform_description),]
Fs7.habit.dist.rev<-rbind.data.frame(Fs7.habit.dist.tree, Fs7.habit.dist.shrub, Fs7.habit.dist.subshrub, Fs7.habit.dist.liana,Fs7.habit.dist.annual,Fs7.habit.dist.biennial,
                                     Fs7.habit.dist.perennial,Fs7.habit.dist.epiphyte,Fs7.habit.dist.lithophyte,Fs7.habit.dist.geophyte,Fs7.habit.dist.helophyte, Fs7.habit.dist.parasite, Fs7.habit.dist.succulent, Fs7.habit.dist.bamboo)
Fs7.habit.dist.non<-Fs7.habit.dist[which(!Fs7.habit.dist$fungi %in% Fs7.habit.dist.rev$fungi),]


growth.form<-c(rep("tree", each=nrow(Fs7.habit.dist.tree)),rep("shrub", each=nrow(Fs7.habit.dist.shrub)), rep("subshrub", each=nrow(Fs7.habit.dist.subshrub)),
       rep("climber or liana", each=nrow(Fs7.habit.dist.liana)),rep("annual", each=nrow(Fs7.habit.dist.annual)),rep("biennial", each=nrow(Fs7.habit.dist.biennial)),
       rep("perennial", each=nrow(Fs7.habit.dist.perennial)),rep("epiphyte", each=nrow(Fs7.habit.dist.epiphyte)),
       rep("lithophyte", each=nrow(Fs7.habit.dist.lithophyte)),rep("geophyte", each=nrow(Fs7.habit.dist.geophyte)),rep("helophyte", each=nrow(Fs7.habit.dist.helophyte)),
       rep("hemi- or holoparasite", each=nrow(Fs7.habit.dist.succulent)),rep("succulent", each=nrow(Fs7.habit.dist.parasite)),rep("bamboo", each=nrow(Fs7.habit.dist.bamboo)), rep("Unknown or unlisted", each=nrow(Fs7.habit.dist.non)))

Fs7.habit2.dist<-rbind.data.frame(Fs7.habit.dist.rev, Fs7.habit.dist.non)
Fs7.habit2.dist$growth.form<-growth.form
saveRDS(Fs7.habit2.dist, "Analyses/int_files/05-05.Fs7.habit2.dist.rds")
#includes Country frequency information.



################
#For WCVP and Pironon tables
#################
Fs7$fungi<-as.character(Fs7$fungi)
#for continent
Fs7.continent <- Fs7 %>%
  group_by(continent, fungi, host) %>%
  distinct(continent, .keep_all = TRUE)
table(Fs7.continent$continent)
#no blanks, only NAs.
#total number of hosts with continent info
length(unique(Fs7.continent$host[which(!is.na(Fs7.continent$continent))])) #10755
length(unique(Fs7.continent$fungi[which(!is.na(Fs7.continent$continent))])) #13086



#for climate zone
Fs7.climate<-Fs7 %>%
  group_by(climate_description, fungi, host) %>%
  distinct(climate_description, .keep_all = TRUE)
table(Fs7.climate$climate_description)
#blanks and NAs for climate
clim<-Fs7.climate[which(!is.na(Fs7.climate$climate_description)),]
clim2<-clim[which(!clim$climate_description==""),]
length(unique(clim2$host)) #10217
length(unique(clim2$fungi)) #12869

#for lifeform
lifeform<-Fs7[which(!is.na(Fs7$lifeform_description)),]
lifeform2<-lifeform[which(!lifeform$lifeform_description==""),]
table(lifeform2$lifeform_description)

length(unique(lifeform2$host)) #9843
length(unique(lifeform2$fungi))#12690

length(grep("^shrub| shrub", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 209163
length(grep("subshrub", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 125353
length(grep("climb|liana", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 74400
length(grep("tree", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 313739
length(grep("annual", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 670874
length(grep("biennial", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 192531
length(grep("perennial", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 1005968
length(grep("epiphyte", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 4125
length(grep("lithophyte", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 3947
length(grep("geophyte", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 301927
length(grep("helophyte", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 67847
length(grep("parasit", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 3969
length(grep("bamboo", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 5254
length(grep("succulent", lifeform2$lifeform_description, ignore.case = TRUE))
#[1] 6965
# total = 6965+5254+3969+67847+301927+3947+4125+1005968+192531+670874+313739+74400+125353+209163 =2986062
# Pironon totals
sum(Fs6$Totals)
#[1] 105338
sum(Fs6$AnimalFood)
#[1] 8555
sum(Fs6$EnvironmentalUses)
#[1] 17912
sum(Fs6$Fuels)
#[1] 4003
sum(Fs6$GeneSources)
#[1] 8703
sum(Fs6$HumanFood)
#[1] 10640
sum(Fs6$InvertebrateFood)
#[1] 2480
sum(Fs6$Materials)
#[1] 14348
sum(Fs6$Medicines)
#[1] 26750
sum(Fs6$Poisons)
#[1] 8405
sum(Fs6$SocialUses)
#[1] 3542

#total number of hosts in human use
length(unique(Fs6$host[which(Fs6$Totals >=1)])) #5845
length(unique(Fs6$fungi[which(Fs6$Totals >=1)])) #10878
#note: these are the hosts with at least one human use, there were plants with no human use in the orig database.


Fs7.area<-Fs7 %>%
  group_by(fungi, host, area) %>%
  distinct(area, .keep_all=TRUE)
sort(table(Fs7.area$area))

#########
# Map fig
#########
Fs7.area.filt<-Fs7.area[!is.na(Fs7.area$area),]
saveRDS(Fs7.area.filt, "Analyses/int_files/05-05.Fs7.area.filt.rds")
#total count of all fungi per country
Fs7.area.filt2 <-Fs7.area.filt  %>%
  group_by(area_code_l3) %>%
  mutate(frequency = length(fungi))
oo<-Fs7.area.filt2[,c(1:15,36:54)]
ooo<-oo %>%
  distinct(area, .keep_all=TRUE)

map.data<-wgsrpd3  %>%
  left_join(ooo, by=c("LEVEL3_COD"="area_code_l3"))

ggplot(map.data)+
  geom_sf(aes(fill=frequency), col="transparent")+
  theme_void()+
  scale_fill_viridis_c(name = "Number of potential\ninteractions", direction=-1)
# saved as 05-05.potential.inter.freq.map.pdf

