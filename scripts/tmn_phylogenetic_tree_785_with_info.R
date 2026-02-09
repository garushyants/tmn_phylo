library(ggtree)
library(treeio)
library(ggtreeExtra)
library(ggplot2)
library(ape)
library(phytools)
library(stringr)
library(tidyr)
library(dplyr)
library(ggnewscale)
library(scales)
library(readxl)
library(openxlsx)

mainpath<-"/Volumes/garushyants/tmn_antitmn/20260128_tmn_padloc_rebuild_tree-relaxed/"
setwd(mainpath)
FigDir<-"tree_figures"

if (!dir.exists(FigDir)){
  dir.create(FigDir)
} else {
  print("Directory already exists!")
}
#Read tree
#I use this method from treio to uncover bootstraps
tree<-read.iqtree("./Tmn_withActiveWalker.IQTree.treefile")


###getting bootsrap info
BootstrapValues<-get.data(tree)
#subset UFboot > 80
BootstrapValuesA80<-subset(BootstrapValues,
                           BootstrapValues$UFboot >= 80)
#I can do additional filtering by length if neccesary

#reroot in the middle
#this is also modified in order to deal with tree
TreeMidRoot<-midpoint.root(tree@phylo)

#########
#Read Info about genomes and PADLOC predictions
ct <- rep("guess", 25)
ct[24] <- "text"
LeavesInfo<-read_xlsx("SupplementaryTable1.xlsx", col_types = ct)

########
#Manipulate basic data
CommonGenusInSet<-LeavesInfo %>% group_by(`Genus name`) %>%
  count()
CommonGenusInSet<-CommonGenusInSet[order(-CommonGenusInSet$n),]

sum(pull(CommonGenusInSet[c(1:12),2]))/length(LeavesInfo$Proteinid)
#I am taking 12 most common genus, because they together account for ~85% of all genomes, and tmn there found at least in >50 genomes per genus
GenusToKeep<-c(pull(CommonGenusInSet[c(1:12),1]))
LeavesInfo$TopGenus<-ifelse(LeavesInfo$`Genus name` %in% GenusToKeep,
                            LeavesInfo$`Genus name`, "Other")

HitsCountsPerLeaf<-LeavesInfo %>%
  group_by(TreeRepresentative)%>%
  count()
TaxonomyCountsPerLeaf<-LeavesInfo %>%
  group_by(TreeRepresentative, TopGenus) %>% 
  count(TopGenus) %>%
  group_by(TreeRepresentative) %>%
  mutate(percent = n / sum(n) * 100)
TaxonomyCountsPerLeaf$TopGenus<-factor(TaxonomyCountsPerLeaf$TopGenus,
                                       levels=c(GenusToKeep,"Other"))

DomainCountsPerLeaf<-LeavesInfo %>%
  group_by(TreeRepresentative, `Domain/Realm name`) %>% 
  count(`Domain/Realm name`) %>%
  group_by(TreeRepresentative) %>%
  mutate(percent = n / sum(n) * 100)
#Class level
#Select most common class
CommonClassInSet<-LeavesInfo %>% group_by(`Class name`) %>%
  count()%>%
  arrange(desc(n))
ClassDf<-LeavesInfo[,c(23,19)]
ClassDf$CommonClass<-ifelse(ClassDf$`Class name` %in% CommonClassInSet$`Class name`[1:11],
                            ClassDf$`Class name`,
                            "Other")
#I pick all Classes with > 15 genomes
ClassCountsPerLeaf<-ClassDf %>%
  group_by(TreeRepresentative, CommonClass) %>% 
  count(CommonClass) %>%
  group_by(TreeRepresentative) %>%
  mutate(percent = n / sum(n) * 100)
unique(ClassCountsPerLeaf$CommonClass)
ClassCountsPerLeaf$CommonClass<-factor(ClassCountsPerLeaf$CommonClass,
                                              levels=c(CommonClassInSet$`Class name`[1:11],
                                                       "Other"))
######Vizualize genome counts better (2025/01/02) !!!!!!
#adding a small value so it will be visually seen which genus is there
HitsCountsPerLeaf$logCount<-log10(HitsCountsPerLeaf$n)+.1

PreExpStrains<-subset(LeavesInfo, !is.na(LeavesInfo$Tmn_variant))
ExpStrains<-unique(PreExpStrains[,c("TreeRepresentative","Tmn_variant")])
#######################################################################################
#Draw basic phylogenetic tree

BasicTreePlot<-ggtree(TreeMidRoot,
                      layout = 'fan',
                      open.angle = 10,
                      size=0.2,
                      color="#636363") 
BasicTreePlot

BasicTreePlotWithBoot<-BasicTreePlot %<+% BootstrapValuesA80 +
  geom_nodepoint(aes(size = UFboot),
                 color = '#4292c6',
                 alpha=.3)+
  scale_size_continuous(range = c(0.01,1))+
  guides(size = guide_legend(nrow = 4))
#BasicTreePlotWithBoot

BasicTreePlotWithLabels<-BasicTreePlotWithBoot%<+% ExpStrains+
  geom_tiplab2(aes(label=Tmn_variant),
               face="ArielMT",
               size=3,
               hjust=-0.4)+
  geom_treescale(y=1, x=4.5, fontsize=3, linesize=0.7, offset=1.2)
BasicTreePlotWithLabels

# ###Let's test
# ggtree(TreeMidRoot,
#        size=0.2,
#        color="#636363") %<+% ExpStrains+
#   geom_tiplab(aes(label=Tmn_variant),
#                face="ArielMT",
#                size=3,
#                hjust=-0.4)+
#   geom_text2(aes(label = node), hjust = -0.3)+
#   geom_treescale(y=1, x=4.5, fontsize=3, linesize=0.7, offset=1.2)
#######################################################################################
###Add clusters
ClusterAncestryNodeOfInterest<-data.frame(ClusterID = c("Ia","Ib","II","III","IV","V","VI","VII","VIII","IX"),
                                          ClMRCA = c(1250,1312,1384,895,799,1053,794,1212,1396,1196))
ClusterAncestryNodeOfInterest$ClusterID<-factor(ClusterAncestryNodeOfInterest$ClusterID,
                                                levels = c("Ia","Ib","II","III","IV","V","VI","VII","VIII","IX"))
#plot with clusters

BasicTreePlotWithClusters <-BasicTreePlotWithLabels + geom_hilight(data = ClusterAncestryNodeOfInterest,
                                                                   aes(node=ClMRCA,
                                                                       fill=ClusterID),
                                                                   alpha=.3,
                                                                   align="right",
                                                                   to.bottom=T)+
  scale_fill_manual(values = c("#8dd3c7","#ffffb3","#bebada","#fb8072",
                               "#80b1d3","#fdb462","#b3de69","#fccde5",
                               "#d9d9d9","#bc80bd"))+
  guides(fill = guide_legend(nrow = 4))
BasicTreePlotWithClusters

###Add ClusterIDs to SupplementaryTable1
get_leaves <- function(tree, node) {
  desc <- getDescendants(tree, node)
  tree$tip.label[desc[desc <= length(tree$tip.label)]]
}

TmnClusterAssignments <- ClusterAncestryNodeOfInterest %>%
  rowwise() %>%
  mutate(TreeRepresentative = list(get_leaves(TreeMidRoot, ClMRCA))) %>%
  unnest(TreeRepresentative)

TableToSaveWithClusters<-merge(LeavesInfo[,c(1:25)],TmnClusterAssignments[,c(1,3)], all.x =T,
                               by= "TreeRepresentative")
# #Save the complete table 1 with clusters
# write.xlsx(TableToSaveWithClusters,
#            file="SupplementaryTable1_all_info.xlsx",
#            rowNames = F)
# write.table(TableToSaveWithClusters,
#            file="SupplementaryTable1_all_info.tsv",
#            sep="\t",
#            row.names = F,
#            quote = F)

#######################################################################################
#######################################################################################
#Adding data on MGE from genomad
LeavesForMGE<-subset(TableToSaveWithClusters, TableToSaveWithClusters$RepresentativeGenome == "Y")
##Plasmids
pathtogenomadfolder<-"./genomad_output/"
filelistgenomadplasmid = list.files(pattern="\\plasmid_summary.tsv$",
                                    recursive = T,
                                    path = pathtogenomadfolder)
setwd(paste0(mainpath,"/",pathtogenomadfolder))
GenomadResultsPlasmids<-readr::read_tsv(filelistgenomadplasmid, id="file_name")
setwd(mainpath)
GenomadResultsPlasmids<-separate(data =  GenomadResultsPlasmids,
                                 col = file_name,
                                 into=c("RefSeqGenomeID",NA,NA),
                                 sep="/")
colnames(GenomadResultsPlasmids)[2]<-"Contig"



PlasmidInfoForPlot<-merge(LeavesForMGE[,c("RefSeqGenomeID",
                                        "Contig",
                                        "TreeRepresentative",
                                        "Proteinid",
                                        "ClusterID")],
                             GenomadResultsPlasmids,
                             by= c("RefSeqGenomeID","Contig"),
                             all.x =T)
nrow(subset(PlasmidInfoForPlot, !is.na(PlasmidInfoForPlot$length)))

#########################
#viruses
filelistgenomadvirus = list.files(pattern="\\virus_summary.tsv$",
                                  recursive = T,
                                  path = pathtogenomadfolder)
setwd(paste0(mainpath,"/",pathtogenomadfolder))
GenomadResultsViruses<-readr::read_tsv(filelistgenomadvirus, id="file_name")
setwd(mainpath)

GenomadResultsViruses<-separate(data = GenomadResultsViruses,
                                col = file_name,
                                into=c("RefSeqGenomeID",NA,NA),
                                sep="/")
GenomadResultsViruses<-separate(data = GenomadResultsViruses,
                                col = seq_name,
                                into=c("Contig","Provirus"),
                                sep="\\|", remove = F)
GenomadResultsViruses<-separate(data = GenomadResultsViruses,
                                col = coordinates,
                                into=c("Start","End"),
                                sep="\\-")

#it is essential to convert to numbers here, because otherwise it is interpreted as string
GenomadResultsViruses$Start<-as.integer(GenomadResultsViruses$Start)
GenomadResultsViruses$End<-as.integer(GenomadResultsViruses$End)
#add coordinates where prophage is a whole contig
GenomadResultsViruses$Start<-ifelse(is.na(GenomadResultsViruses$Start),
                                    1,
                                    GenomadResultsViruses$Start)
GenomadResultsViruses$End<-ifelse(is.na(GenomadResultsViruses$End),
                                  as.integer(GenomadResultsViruses$length),
                                    GenomadResultsViruses$End)

GenomadResultsVirusesCons<-subset(GenomadResultsViruses,
                                  GenomadResultsViruses$virus_score > .8)


PreGenomadVirusesOnTmnContigs<-merge(LeavesForMGE[,c(1:8,26)],
                                     GenomadResultsVirusesCons, by=c("RefSeqGenomeID","Contig"))
PreGenomadVirusesOnTmnContigs$InProphage<-ifelse((PreGenomadVirusesOnTmnContigs$Start.x >= PreGenomadVirusesOnTmnContigs$Start.y) & 
                                                   (PreGenomadVirusesOnTmnContigs$End.x <= PreGenomadVirusesOnTmnContigs$End.y) &
                                                   (PreGenomadVirusesOnTmnContigs$Start.x <= PreGenomadVirusesOnTmnContigs$End.y) & 
                                                   (PreGenomadVirusesOnTmnContigs$End.x >= PreGenomadVirusesOnTmnContigs$Start.y),
                                                 1,0)

GenomadVirusesOnTmnContigs<-subset(PreGenomadVirusesOnTmnContigs,
                                   InProphage == 1)


##merging plasmid and Viral data
TmnVirPlas<-merge(PlasmidInfoForPlot[,c(1:3,6,10)],
                  GenomadVirusesOnTmnContigs[,c(1:3,18,23)],
                  all =T,
                  by = c("RefSeqGenomeID","Contig","TreeRepresentative"))
#####
TmnVirPlas$InProphage<-ifelse(TmnVirPlas$virus_score>0, "V", NA)
TmnVirPlas$InPlasmid<-ifelse(TmnVirPlas$length>0, "P", NA)
MGEDataForPlot<-TmnVirPlas[,c("TreeRepresentative","InPlasmid","InProphage")]
# #save to file, so I don't have to reload it again later
# write.table(MGEDataForPlot,
#             file="Tmn_785_MGE_data_genomad.summary.tsv",
#             sep="\t",
#             row.names = F,
#             quote = F)
colnames(MGEDataForPlot)[2:3]<-c("Plasmid","Prophage")
#MGEDataForPlot[is.na(MGEDataForPlot)]<-0
MGEDataForPlotLong<-gather(MGEDataForPlot,
                           key = "Location", value = "Prediction", Plasmid:Prophage)

#######################################################################################
#######################################################################################

######################
###Plot on the phylogenetic tree
#Colors
# colors<-c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628",
#           "#ffffbf","#a6cee3",
#           "#d9d9d9")
genuscolors<-c("#e31a1c","#33a02c","#1f78b4","#6a3d9a",
               "#ff7f00","#b15928","#ffff99","#fb9a99",
               "#fdbf6f","#b2df8a","#a6cee3","#cab2d6", "#d9d9d9")
names(genuscolors)<-c(GenusToKeep,"Other")

classcolors<-c("#8dd3c7","#ffffb3","#bebada","#fb8072",
               "#80b1d3","#fdb462","#b3de69","#fccde5",
               "#bc80bd","#ccebc5","#ffed6f","#d9d9d9")
names(classcolors)<-c(CommonClassInSet$`Class name`[1:11],"Other")
#Add MGE info
TreeWithMGE<-BasicTreePlotWithClusters +
  new_scale_fill()+
  geom_fruit(data=MGEDataForPlotLong,
             geom = geom_tile,
             mapping = aes(y=TreeRepresentative, x= Location,
                           fill = as.character(Prediction)),
             color="#bdbdbd",
             pwidth=.05,
             offset = 0.18,
             axis.params=list(axis="x",
                              text.size=2,
                              text.angle=60,
                              vjust=0,
                              hjust=1))+
  scale_fill_manual(values=c("#238443","#0570b0"), guide="none", na.value = "#ffffff")
TreeWithMGEAndTax<-TreeWithMGE +
  new_scale_fill()+
  geom_fruit(data=DomainCountsPerLeaf,
             geom = geom_col,
             mapping = aes(y=TreeRepresentative,
                           fill = `Domain/Realm name`,
                           x = percent),
             offset =0.05, pwidth =.03)+
  scale_fill_manual(values=c("#b2182b","#9ecae1"), name = "Domain")+
  guides(fill = guide_legend(nrow = 3))+
  new_scale_fill()+
  geom_fruit(data=ClassCountsPerLeaf,
             geom = geom_col,
             mapping = aes(y=TreeRepresentative,
                           fill = CommonClass,
                           x = percent),
             offset =0.01, pwidth =.05)+
  scale_fill_manual(values=classcolors, name = "Class")+
  guides(fill = guide_legend(nrow = 4))+
  new_scale_fill()+
  geom_fruit(data=TaxonomyCountsPerLeaf,
           geom = geom_col,
           mapping = aes(y=TreeRepresentative,
                         fill = TopGenus,
                         x = percent),
           offset =0.01, pwidth =.2)+
  scale_fill_manual(values=genuscolors, name = "Genus")+
  guides(fill = guide_legend(nrow = 4))

MainTreeToSave<-TreeWithMGEAndTax +
  geom_fruit(data = HitsCountsPerLeaf,
             geom = geom_col,
             mapping = aes(y=TreeRepresentative,
                           x=logCount), fill= "#878787",
             pwidth =.4,
             axis.params=list(axis="x",
                              text.size=3,
                              line.size=.3),
             grid.params = list(size=.3,
                                alpha=.3))+
  theme(legend.position = "bottom")
MainTreeToSave  
  
####################
#Save main tree figure
ggsave("tmn_785_tree_with_MGE_and_Tax_v3.png",
       plot=MainTreeToSave,
       path=FigDir,
       width=30,
       height=30,
       dpi=300,
       units="cm")

ggsave("tmn_785_tree_with_MGE_and_Tax_v3.svg",
       plot=MainTreeToSave,
       path=FigDir,
       width=30,
       height=30,
       dpi=300,
       units="cm")





#############################################################
#############################################################
##Get protein lengths per cluster

clustercolors<-c("#8dd3c7","#ffffb3","#bebada","#fb8072",
                 "#80b1d3","#fdb462","#b3de69","#fccde5",
                 "#d9d9d9","#bc80bd")

ProteinLengthsDf<-subset(TableToSaveWithClusters, TableToSaveWithClusters$RepresentativeGenome == "Y" &
                           !is.na(TableToSaveWithClusters$ClusterID))

ProteinLengthsDf$length<-(ProteinLengthsDf$End - ProteinLengthsDf$Start+1)/3

min(ProteinLengthsDf$length)
max(ProteinLengthsDf$length)
median(ProteinLengthsDf$length)

ProtLengthPerCluster<-ggplot(ProteinLengthsDf,
                             aes(x = length, fill = ClusterID, group = ClusterID))+
  geom_histogram(color="#bdbdbd", binwidth =20)
#######
###add Median
hist_data <- ggplot_build(ProtLengthPerCluster)$data[[1]]
hist_data$ClusterID <- rep(levels(ProteinLengthsDf$ClusterID), each = length(unique(hist_data$x)))
#Compute max count per ClusterID
max_counts <- hist_data %>%
  group_by(ClusterID) %>%
  summarize(max_count = max(count))

medians <- ProteinLengthsDf %>% 
  group_by(ClusterID) %>% 
  summarize(median_value = median(length))
medians_max<-merge(medians, max_counts, by="ClusterID")
####
ProtLengthWithMedian<-ProtLengthPerCluster+ 
  geom_vline(data = medians_max, aes(xintercept = median_value), color = "#cb181d", linewidth =1.2)+
  geom_text(data = medians_max, 
          aes(x = median_value, y = max_count*.9,
              label = paste0("", round(median_value,0))), hjust = -0.1)+
  facet_wrap(~ClusterID, ncol =1, scales = "free_y")+
  scale_fill_manual(values = clustercolors,
                    guide="none")+
  xlim(900,1500)+
  xlab("Protein length (aa)")+
  ylab("")+
  theme_minimal()+
  theme(strip.text = element_text(size=14),
        axis.text = element_text(size=12))
ProtLengthWithMedian
ggsave("tmn_cluster_lengths.png",
       plot=ProtLengthWithMedian,
       path=FigDir,
       width=8,
       height=25,
       dpi=300,
       units="cm")

ggsave("tmn_cluster_lengths.svg",
       plot=ProtLengthWithMedian,
       path=FigDir,
       width=8,
       height=25,
       dpi=300,
       units="cm")
