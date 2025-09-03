#Descriptive taxonomical analysis

#Part 1:
#Make a full fungal cladogram for all fungi detected, highlight those on living tissue

setwd("~/OneDrive - UBC/DB_fung_2025.01/")

library(dplyr)
library(fungarium)
library(jsonlite)
library(stringr)
library(data.table)
library(ape)
library(ggtree)
library(ggnewscale)
library(FUNGuildR)
library(viridis)
library(ggplot2)
library(ggtreeExtra)
library(pals)
#1. Read in local fungal data, create cladogram
#make phylogeny for fungi (from in local area db)

Fung.local.raw<-readRDS("Final_comb_DBs/LOCAL.loc_cleaned_unique_standardized_filtered.noerror.uniq.rds")
#55893 unique fungi, 37,869 unique plants
F.loc.1<-Fung.local.raw %>%
  group_by(new_species) %>%
  distinct(new_kingdom, new_phylum, new_class, new_order, new_family, new_genus, new_species, .keep_all = TRUE) 
F.loc.2<-F.loc.1 %>%
  distinct(new_kingdom, new_phylum, new_class,new_order,new_family,new_genus, new_species) 

F.loc.3<- na.omit(F.loc.2)
F.loc.3<-F.loc.3[which(!F.loc.3$new_species==""),]
F.loc.3$new_species<-as.factor(F.loc.3$new_species)
F.loc.3$new_genus<-as.factor(F.loc.3$new_genus)
F.loc.3$new_family<-as.factor(F.loc.3$new_family)
F.loc.3$new_order<-as.factor(F.loc.3$new_order)
F.loc.3$new_class<-as.factor(F.loc.3$new_class)
F.loc.3$new_phylum<-as.factor(F.loc.3$new_phylum)
F.loc.3$new_kingdom<-as.factor(F.loc.3$new_kingdom)


F.loc.3.tree <- as.phylo(~new_kingdom/new_phylum/new_class/new_order/new_family/new_genus/new_species, data =F.loc.3)
#Plot as square phylogeny
ggtree(F.loc.3.tree, layout = "rectangular")

#read in "on plant"
F.on<-readRDS("Final_comb_DBs/ON.PLANT.standardized.cleaned.no.error.uniq.rds")
F.on.1<-F.on %>%
  group_by(new_species) %>%
  distinct(new_kingdom, new_phylum, new_class, new_order, new_family, new_genus, new_species, .keep_all = TRUE) 
F.on.2<-F.on.1 %>%
  distinct(new_kingdom, new_phylum, new_class,new_order,new_family,new_genus, new_species) 
F.on.3<- na.omit(F.on.2)

#read in "living plant"
F.liv<-readRDS("Final_comb_DBs/LIVING.standardized.cleaned.no.error.uniq.rds")
F.liv.1<-F.liv %>%
  group_by(new_species) %>%
  distinct(new_kingdom, new_phylum, new_class, new_order, new_family, new_genus, new_species, .keep_all = TRUE) 
F.liv.2<-F.liv.1 %>%
  distinct(new_kingdom, new_phylum, new_class,new_order,new_family,new_genus, new_species) 
F.liv.3<- na.omit(F.liv.2)




#start a df 
Local_fungi<-F.loc.3.tree$tip.label
Fungal_classes<-as.character(F.loc.3$new_class[match(Local_fungi, F.loc.3$new_species)])
Fungi_on_plant<-F.loc.3.tree$tip.label %in% F.on.3$new_species
Fungi_on_plant[which(Fungi_on_plant=="FALSE")]<-"Not found on plant"
Fungi_on_plant[which(Fungi_on_plant=="TRUE")]<-"Found on plant"
Fungi_on_living_plant<-F.loc.3.tree$tip.label %in% F.liv.3$new_species
Fungi_on_living_plant[which(Fungi_on_living_plant=="FALSE")]<-"Not found on living plant"
Fungi_on_living_plant[which(Fungi_on_living_plant=="TRUE")]<-"Found on living plant"                                
Fungi_on_plant<-data.frame(Fungi_on_plant)
rownames(Fungi_on_plant)<-Local_fungi
Fungi_on_living_plant<-data.frame(Fungi_on_living_plant)
rownames(Fungi_on_living_plant)<-Local_fungi

circb<-ggtree(F.loc.3.tree, layout="circular")
circb$data$x[circb$data$x == "6"] <- 5.03
saveRDS(circb, "Analyses/int_files/05-05.circb.rds")

Flt3.class<-cbind.data.frame(class=circb$data$label[grep("mycetes$", word(circb$data$label))], id=as.numeric(circb$data$node[grep("mycetes$", word(circb$data$label))]))
Flt3.order<-cbind.data.frame(order=circb$data$label[grep("ales$", word(circb$data$label))], id=as.numeric(circb$data$node[grep("ales$", word(circb$data$label))]))
Flt3.phylum<-cbind.data.frame(phylum=circb$data$label[grep("mycota$", word(circb$data$label))], id=as.numeric(circb$data$node[grep("mycota$", word(circb$data$label))]))

F.p<-Flt3.phylum[which(Flt3.phylum$phylum %in% c("Ascomycota", "Basidiomycota")),]
F.c<-Flt3.class[which(Flt3.class$class %in% c("Sordariomycetes", "Leotiomycetes", "Agaricomycetes", "Dothideomycetes", "Ustilaginomycetes", "Pucciniomycetes", "Eurotiomycetes", "Exobasidiomycetes")),]
F.o<-Flt3.order[which(Flt3.order$order %in% c("Agaricales", "Polyporales", "Hypocreales", "Xylariales", "Diaporthales", "Eurotiales", "Mycosphaerellales", "Pleosporales", 
                                              "Botryosphaeriales", "Capnodiales", "Pucciniales", "Helotiales", "Ustilaginales", "Glomerellales", "Amphisphaeriales", "Phyllachorales")),]

hi.fung.id<-c(Flt3.phylum$id, F.c$id, F.o$id)
hi.fung.name<-c(Flt3.phylum$phylum, F.c$class, F.o$order)
hi.fung<-cbind.data.frame(hi.fung.id, hi.fung.name)

circc<-circb+geom_hilight(data=F.o, mapping=aes(fill=order, node=id), alpha=0.3, extend=-2.5)+scale_fill_manual(values=unname(polychrome(16)))+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.title.position = "top", legend.key.size = unit(10, "lines"), legend.position="bottom")#scale_fill_viridis_d(option="H")#geom_cladelab(node=Flt3.class$id, label=Flt3.class$class, align=TRUE,  offset=-1) #angle="auto"

circc1 <- circc + new_scale_fill() 
circd<-circc1+geom_hilight(data=F.c, mapping=aes(fill=class, node=id), alpha=0.3, extend=-3.5)+scale_fill_manual(values=unname(tol(8)))+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.key.size = unit(1, "lines"))#scale_fill_viridis_d(option="B")+ #geom_cladelab(node=Flt3.class$id, label=Flt3.class$class, align=TRUE,  offset=-1) #angle="auto"

circd1 <- circd + new_scale_fill() 
#use this for legend. (should match other, so maybe not important)
circe<-circd1+geom_hilight(data=F.p, mapping=aes(fill=phylum, node=id), alpha=0.3, extend=-4.5)+scale_fill_manual(values=unname(alphabet(2)))+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.key.size = unit(1, "lines"))

#this is the one we will use for final plot
circe.u<-circb+geom_hilight(data=F.o, mapping=aes(fill=order, node=id), alpha=0.3, extend=-1.5)+scale_fill_manual(values=unname(polychrome(16)))+guides(fill="none")+
  new_scale_fill()+geom_hilight(data=F.c, mapping=aes(fill=class, node=id), alpha=0.3, extend=-2.5)+scale_fill_manual(values=unname(tol(8)))+guides(fill="none")+
  new_scale_fill()+geom_hilight(data=F.p, mapping=aes(fill=phylum, node=id), alpha=0.3, extend=-3.5)+scale_fill_manual(values=unname(alphabet(2)))+guides(fill="none")
circf <- circe.u + new_scale_fill() 

#Explore fungal types using funguild data
fg<-get_funguild_db(db = "http://www.stbates.org/funguild_db_2.php")
saveRDS(fg, "ext_data/2025.04.06.Funguild_all.rds")
#break down to genus level
fg.gen<-fg[which(fg$taxonomicLevel=="13"),]

#SIDE CALC: to calculate number of genera in one of our databases vs. endophytes/pathos sensu FUNGuild.
t.liv<-fg.gen[which(fg.gen$taxon %in% F.liv$new_genus ),]
t.liv2<-t.liv[grep("Plant Parasite|Plant Pathogen|Endophyte", t.liv$guild),]
t.liv.nonendopath<-t.liv[which(!t.liv$taxon %in% t.liv2$taxon),]
t.liv.sapro<-t.liv.nonendopath[grep("Saprotroph", t.liv.nonendopath$guild),]
#2501 total genera in living plant db. 2196 also in funguild. 1141 as pathos/endos. 1141/2501=0.4562175 (all taxa). 1141/2196 = 0.5195811 (only genera found in funguild)
#1055 not labeled as pathos/pars/endos. of these, 961 are saprotrophs. 961/2501 (all taxa)=0.3842463. 961/2196=0.4376138 (only genera found in funguild)

t.on<-fg.gen[which(fg.gen$taxon %in% F.on$new_genus ),]
t.on2<-t.on[grep("Plant Parasite|Plant Pathogen|Endophyte", t.on$guild),]
t.on.nonendopath<-t.on[which(!t.on$taxon %in% t.on2$taxon),]
t.on.sapro<-t.on.nonendopath[grep("Saprotroph", t.on.nonendopath$guild),]
#6395 total genera in on plant db. 4966 also in funguild. 1954 as pathos/endos. 1954/6395=0.3055512 (all taxa). 1954/4966=0.3934756 (only genera found in funguild)
#3012 not labeled as pathos/pars/endos. of these, 2496 are saprotrophs. 2496/6395 (all taxa)=0.3903049. 2496/4966=0.5026178 (only genera found in funguild)

t.loc<-fg.gen[which(fg.gen$taxon %in% Fung.local.raw$new_genus ),]
t.loc2<-t.loc[grep("Plant Parasite|Plant Pathogen|Endophyte", t.loc$guild),]
t.loc.nonendopath<-t.loc[which(!t.loc$taxon %in% t.loc2$taxon),]
t.loc.sapro<-t.loc.nonendopath[grep("Saprotroph", t.loc.nonendopath$guild),]
#6581 total genera in local db. 5134 also in funguild. 1960 as pathos/endos. 1960/6581=0.2978271 (all taxa). 1960/5134=0.3817686
#3174 not labeled as pathos/pars/endos. of these, 2547 are saprotrophs. 2547/6581 (all taxa)=0.3870232. 2547/5134=0.4961044 (only genera found in funguild)
#how many fg.gen are plant pathogen/parasite endophyte?
tt<-grep("Plant Parasite|Plant Pathogen|Endophyte", fg.gen$guild)
#



#match to genera 
#first make sure F.loc.3 is usable
F.loc.3$new_genus<-as.character(F.loc.3$new_genus)
#incorporate funguild gen data
funguild.pres<-word(F.loc.3.tree$tip.label, 1) %in% fg.gen$taxon
funguild.pres<-fg.gen[match(word(F.loc.3.tree$tip.label, 1), fg.gen$taxon),]
funguild.pres$present<-word(F.loc.3.tree$tip.label, 1) %in% word(funguild.pres$taxon, 1)
funguild.pres$tip.label<-F.loc.3.tree$tip.label
funguild.pres$guild[which(funguild.pres$present=="FALSE")]<-"No guild data"
fung.guild<-funguild.pres$guild
fung.guild<-gsub("\\|","",fung.guild)
endos<-which(grepl("Endophyte", fung.guild))
pathos<-which(grepl("Plant Parasite|Plant Pathogen", fung.guild))
endos.pathos<-vector(length=length(funguild.pres$guild))
endos.pathos[pathos]<-"Plant Pathogen/Parasite"
endos.pathos[endos]<-"Endophyte"
endos.pathos[which(grepl("Plant Parasite|Plant Pathogen", fung.guild) & grepl("Endophyte", fung.guild))]<-"Endophyte or Plant Pathogen/Parasite"
endos.pathos[endos.pathos=="FALSE"]<-"Not classified as Endophyte or Plant Pathogen/Parasite"
endos.pathos<-data.frame(endos.pathos)
rownames(endos.pathos)<-Local_fungi


p1 <- circf + new_scale_fill() 
p1b <- gheatmap(p1, endos.pathos, offset=-.4, width=.15,
                colnames=FALSE, color = FALSE) + #colnames_angle=30, colnames_offset_y = 0, custom_column_labels = "Guild"
  scale_fill_viridis_d(option="D", name="Guild")+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.title.position = "top", legend.key.size = unit(10, "lines"), legend.position="right")

#p1b

p2b <- p1b + new_scale_fill() 
p3b<- gheatmap(p2b, Fungi_on_living_plant, offset=0.6, width=.15,
         colnames=FALSE, color=FALSE) + # colnames_angle=30, colnames_offset_y = 0, colnames_position = 'top', colnames_offset_y = .0001, custom_column_labels = "On living Plant"
  scale_fill_viridis_d(option="G", name="On living plant")+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.title.position = "top", legend.key.size = unit(10, "lines"))
#p3b

p4b <- p3b + new_scale_fill() 
gheatmap(p4b, Fungi_on_plant, offset=1.6, width=0.15,
          colnames=FALSE, color=FALSE)+ # colnames_angle=30, colnames_offset_y = .25, custom_column_labels = "On plant"
  scale_fill_viridis_d(option="C", name="On plant")+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.title.position = "top", legend.key.size = unit(10, "lines"))+
  annotate("text", x=0, y=0, label="order", size=30, color="black", vjust=-15, hjust=-.1)+
  annotate("text", x=0, y=0, label="class", size=30, color="black", vjust=-10, hjust=-.1)+
  annotate("text", x=0, y=0, label="phylum", size=30, color="black", vjust=-5, hjust=-.1)

#ALT for rectangular cladogram:
gheatmap(p4b, Fungi_on_plant, offset=1.6, width=0.15,
         colnames=FALSE, color=FALSE)+ # colnames_angle=30, colnames_offset_y = .25, custom_column_labels = "On plant"
  scale_fill_viridis_d(option="C", name="On plant")+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.title.position = "top", legend.key.size = unit(10, "lines"))


