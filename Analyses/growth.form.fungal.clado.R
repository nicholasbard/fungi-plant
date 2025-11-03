# Taxonomical Analysis - Growth Form figure.


# Note: Requires "Analyses/int_files/04-10.Fs7.habit2.dist.rds" from Co-occurrence.analysis.R 

setwd("~/OneDrive - UBC/DB_fung_2025.01/")
library(stringr)
library(ape)
library(dplyr)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(pals)
library(purrr)
#1. Read in local fungal data, create cladogram
#make phylogeny for fungi (from in local area db)

Fung.living.raw<-readRDS("Final_comb_DBs/LIVING.standardized.cleaned.no.error.uniq.rds")
F.liv.1<-Fung.living.raw %>%
  group_by(new_species) %>%
  distinct(new_kingdom, new_phylum, new_class, new_order, new_family, new_genus, new_species, .keep_all = TRUE) 
F.liv.1$new_order[which(is.na(F.liv.1$new_order))]<-""
F.liv.1$new_class[which(is.na(F.liv.1$new_class))]<-""
F.liv.1$new_family[which(is.na(F.liv.1$new_family))]<-""

F.liv.3<-F.liv.1 %>%
  distinct(new_kingdom, new_phylum, new_class,new_order,new_family,new_genus, new_species, .keep_all = TRUE) 

F.liv.3$new_species<-as.factor(F.liv.3$new_species)
F.liv.3$new_genus<-as.factor(F.liv.3$new_genus)
F.liv.3$new_family<-as.factor(F.liv.3$new_family)
F.liv.3$new_order<-as.factor(F.liv.3$new_order)
F.liv.3$new_class<-as.factor(F.liv.3$new_class)
F.liv.3$new_phylum<-as.factor(F.liv.3$new_phylum)
F.liv.3$new_kingdom<-as.factor(F.liv.3$new_kingdom)


F.liv.3.tree <- as.phylo(~new_kingdom/new_phylum/new_class/new_order/new_family/new_genus/new_species, data =F.liv.3)

#make initial ggtree
circ.liv3<-ggtree(F.liv.3.tree, layout="circular")
circ.liv3$data$x[circ.liv3$data$x == "6"] <- 5.03



#assess higher taxa
Flt3.class<-cbind.data.frame(class=circ.liv3$data$label[grep("mycetes$", word(circ.liv3$data$label))], id=as.numeric(circ.liv3$data$node[grep("mycetes$", word(circ.liv3$data$label))]))
Flt3.order<-cbind.data.frame(order=circ.liv3$data$label[grep("ales$", word(circ.liv3$data$label))], id=as.numeric(circ.liv3$data$node[grep("ales$", word(circ.liv3$data$label))]))
Flt3.phylum<-cbind.data.frame(phylum=circ.liv3$data$label[grep("mycota$", word(circ.liv3$data$label))], id=as.numeric(circ.liv3$data$node[grep("mycota$", word(circ.liv3$data$label))]))

F.p<-Flt3.phylum[which(Flt3.phylum$phylum %in% c("Ascomycota", "Basidiomycota")),]
F.c<-Flt3.class[which(Flt3.class$class %in% c("Sordariomycetes", "Leotiomycetes", "Agaricomycetes", "Dothideomycetes", "Ustilaginomycetes", "Pucciniomycetes", "Eurotiomycetes", "Exobasidiomycetes")),]
F.o<-Flt3.order[which(Flt3.order$order %in% c("Agaricales", "Polyporales", "Hypocreales", "Xylariales", "Diaporthales", "Eurotiales", "Mycosphaerellales", "Pleosporales", 
                                              "Botryosphaeriales", "Capnodiales", "Pucciniales", "Helotiales", "Ustilaginales", "Glomerellales", "Amphisphaeriales", "Phyllachorales")),]

hi.fung.id<-c(Flt3.phylum$id, F.c$id, F.o$id)
hi.fung.name<-c(Flt3.phylum$phylum, F.c$class, F.o$order)
hi.fung<-cbind.data.frame(hi.fung.id, hi.fung.name)



circ.liv4<-circ.liv3+geom_hilight(data=F.o, mapping=aes(fill=order, node=id), alpha=0.3, extend=-5.5)+scale_fill_manual(values=unname(polychrome(16)))+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.title.position = "top", legend.key.size = unit(10, "lines"), legend.position="bottom")#scale_fill_viridis_d(option="H")#geom_cladelab(node=Flt3.class$id, label=Flt3.class$class, align=TRUE,  offset=-1) #angle="auto"

circ.liv4b <- circ.liv4 + new_scale_fill() 
circ.liv5<-circ.liv4b+geom_hilight(data=F.c, mapping=aes(fill=class, node=id), alpha=0.3, extend=-3.5)+scale_fill_manual(values=unname(tol(8)))+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.key.size = unit(10, "lines"))#scale_fill_viridis_d(option="B")+ #geom_cladelab(node=Flt3.class$id, label=Flt3.class$class, align=TRUE,  offset=-1) #angle="auto"

circ.liv5b <- circ.liv5 + new_scale_fill() 
circ.liv6<-circ.liv5b+geom_hilight(data=F.p, mapping=aes(fill=phylum, node=id), alpha=0.3, extend=-4.5)+scale_fill_manual(values=unname(alphabet(2)))+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.key.size = unit(10, "lines"))#scale_fill_viridis_d(option="E")#name = "higher taxon", geom_cladelab(node=Flt3.class$id, label=Flt3.class$class, align=TRUE,  offset=-1) #angle="auto"
#### save image for legend - 12000 wide
#this one will have no legend.
circ.liv6b<-circ.liv3+geom_hilight(data=F.o, mapping=aes(fill=order, node=id), alpha=0.3, extend=-1.5)+scale_fill_manual(values=unname(polychrome(16)))+guides(fill="none")+
  new_scale_fill()+geom_hilight(data=F.c, mapping=aes(fill=class, node=id), alpha=0.3, extend=-2.5)+scale_fill_manual(values=unname(tol(8)))+guides(fill="none")+
  new_scale_fill()+geom_hilight(data=F.p, mapping=aes(fill=phylum, node=id), alpha=0.3, extend=-3.5)+scale_fill_manual(values=unname(alphabet(2)))+guides(fill="none")
circ.liv6c <- circ.liv6b + new_scale_fill() 

#match to genera 
#first make sure F.liv.3 is usable
F.liv.3$new_genus<-as.character(F.liv.3$new_genus)


# OK, lets do growth habit now.
Fs7.habit2.dist<-readRDS("Analyses/int_files/05-05.Fs7.habit2.dist.rds")
habits <- model.matrix(~ growth.form - 1, data = Fs7.habit2.dist)
colnames(habits) <- sub("growth.form", "", colnames(habits))

habitsdf<-data.frame(habits)
fungi<-Fs7.habit2.dist$fungi
habitsdf2<-cbind.data.frame(fungi, habitsdf)
#compress everything into unique fungi rows
habitsdf3 <- habitsdf2 %>%
  group_by(fungi) %>%
  summarize(across(everything(), ~ ifelse(sum(.) > 0, 1, 0))) 
habitsdf3$fungi<-as.character(habitsdf3$fungi)



habitsdf4<-habitsdf3[match(F.liv.3.tree$tip.label, habitsdf3$fungi),]



habitsdf4$annual[habitsdf4$annual==1]<-"annual"
habitsdf4$annual[habitsdf4$annual==0]<-""
habitsdf4$biennial[habitsdf4$biennial==1]<-"biennial"
habitsdf4$biennial[habitsdf4$biennial==0]<-""
habitsdf4$perennial[habitsdf4$perennial==1]<-"perennial"
habitsdf4$perennial[habitsdf4$perennial==0]<-""
habitsdf4$shrub[habitsdf4$shrub==1]<-"shrub"
habitsdf4$shrub[habitsdf4$shrub==0]<-""
habitsdf4$subshrub[habitsdf4$subshrub==1]<-"subshrub"
habitsdf4$subshrub[habitsdf4$subshrub==0]<-""
habitsdf4$tree[habitsdf4$tree==1]<-"tree"
habitsdf4$tree[habitsdf4$tree==0]<-""
habitsdf4$climber.or.liana[habitsdf4$climber.or.liana==1]<-"climber or liana"
habitsdf4$climber.or.liana[habitsdf4$climber.or.liana==0]<-""
habitsdf4$bamboo[habitsdf4$bamboo==1]<-"bamboo"
habitsdf4$bamboo[habitsdf4$bamboo==0]<-""
habitsdf4$succulent[habitsdf4$succulent==1]<-"succulent"
habitsdf4$succulent[habitsdf4$succulent==0]<-""
habitsdf4$epiphyte[habitsdf4$epiphyte==1]<-"epiphyte"
habitsdf4$epiphyte[habitsdf4$epiphyte==0]<-""
habitsdf4$geophyte[habitsdf4$geophyte==1]<-"geophyte"
habitsdf4$geophyte[habitsdf4$geophyte==0]<-""
habitsdf4$helophyte[habitsdf4$helophyte==1]<-"helophyte"
habitsdf4$helophyte[habitsdf4$helophyte==0]<-""
habitsdf4$lithophyte[habitsdf4$lithophyte==1]<-"lithophyte"
habitsdf4$lithophyte[habitsdf4$lithophyte==0]<-""
habitsdf4$hemi..or.holoparasite[habitsdf4$hemi..or.holoparasite==1]<-"hemi- or holoparasite"
habitsdf4$hemi..or.holoparasite[habitsdf4$hemi..or.holoparasite==0]<-""


habitsdf4<-data.frame(habitsdf4[,-1])
rownames(habitsdf4)<-F.liv.3.tree$tip.label
habitsdf4 <- habitsdf4[, c("annual", "biennial", "perennial", "shrub", "subshrub", "tree", "climber.or.liana", "bamboo", 
                           "succulent","epiphyte", "geophyte", "helophyte", "lithophyte", "hemi..or.holoparasite")]
gheatmap(circ.liv6c, habitsdf4, offset=0.1, width=0.8, font.size=3, colnames=FALSE, hjust=0) +
  scale_fill_manual(breaks=c("annual", "biennial", "perennial", "shrub", "subshrub", "tree", "climber or liana", "bamboo", 
                             "succulent","epiphyte", "geophyte", "helophyte", "lithophyte", "hemi- or holoparasite"), 
                    values=c("yellow4","grey50", "lightsalmon2", "maroon", "steelblue", "darkgreen", "orange3", "mediumorchid4", "brown", "tan4", "slateblue4","darkseagreen4",
                             "mediumblue",  "grey0", "purple4", "tomato3", "turquoise4"), 
                    name="Growth Form", na.value = "white")+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.key.size = unit(10, "lines"))+
  annotate("text", x=0, y=0, label="order", size=30, color="black", vjust=-11.5)+
  annotate("text", x=0, y=0, label="class", size=30, color="black", vjust=-7.5)+
  annotate("text", x=0, y=0, label="phylum", size=30, color="black", vjust=-3.5)
#saved at 12000 pixel width X 8500 height

#Alt rectangular

gheatmap(circ.liv6c, habitsdf4, offset=0.1, width=0.8, font.size=3, colnames=FALSE, hjust=0) +
  scale_fill_manual(breaks=c("annual", "biennial", "perennial", "shrub", "subshrub", "tree", "climber or liana", "bamboo", 
                             "succulent","epiphyte", "geophyte", "helophyte", "lithophyte", "hemi- or holoparasite"), 
                    values=c("yellow4","grey50", "lightsalmon2", "maroon", "steelblue", "darkgreen", "orange3", "mediumorchid4", "brown", "tan4", "slateblue4","darkseagreen4",
                             "mediumblue",  "grey0", "purple4", "tomato3", "turquoise4"), 
                    name="Growth Form", na.value = "white")+
  theme(legend.title = element_text(size = 85, face="bold"), legend.text = element_text(size = 85), legend.key.size = unit(10, "lines"))
