library(rstudioapi)
library(readxl)
library(dplyr)
library(tidyr)
library(treeio)
library(phytools)
library(ape)
library(ggtree)
library(ggplot2)
library(Biostrings)
library(ggmsa)
library(stringr)
#library(cowplot)
library(ggpubr)
library(ggplotify)

##############
mainpath<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mainpath)
##############
#Read info about leaves required for tree construction
ct <- rep("guess", 26)
ct[24] <- "text"
LeavesInfo<-read_xlsx("../data/SupplementaryTable1_all_info.xlsx", col_types = ct)

######
LeavesRepresentativesInfo<-subset(LeavesInfo,
                                  LeavesInfo$RepresentativeGenome == "Y")
#######get tree
tree<-read.iqtree("../data/Tmn_withActiveWalker.IQTree.treefile")
###getting bootsrap info
BootstrapValues<-get.data(tree)
#subset UFboot > 80
BootstrapValuesA80<-subset(BootstrapValues,
                           BootstrapValues$UFboot >= 80)
#I can do additional filtering by length if neccesary

#reroot in the middle
TreeMidRoot<-midpoint.root(tree@phylo)
####################################
####################################
#Plotting tmn that occur in more than one species
GenomesWithMultTmns<-LeavesInfo %>% group_by(RefSeqGenomeID) %>%
  count() %>%
  arrange(desc(n))%>%
  filter( n > 1)
#I have 25 of those
MultipleCopiesDf<-subset(LeavesInfo,
                         LeavesInfo$RefSeqGenomeID %in% GenomesWithMultTmns$RefSeqGenomeID)

MultiCopiesForPlot<-MultipleCopiesDf[,c("RefSeqGenomeID","TreeRepresentative",
                                        "Genus name")] %>% 
  arrange(RefSeqGenomeID) %>%
  mutate(names = rep(c("from","to"),length(MultipleCopiesDf$Proteinid)/2)) %>%
  pivot_wider(
    id_cols = c(RefSeqGenomeID, `Genus name`),
    names_from = names,
    values_from = TreeRepresentative) %>%
  mutate(label = paste(RefSeqGenomeID, `Genus name`))%>%
  arrange(from,to)%>%
  mutate(size = c(seq(1,3.1,by =0.3),rep(3,8),8,4,2,2.2,3,10,7,2,4))

MultiCopiesForPlotF<-merge(MultiCopiesForPlot,
                           LeavesRepresentativesInfo[,c("RepresentativeGenome","RefSeqGenomeID")],
                           by = "RefSeqGenomeID", all.x =T)  
###get all clusters
ClusterNodesAll<-LeavesRepresentativesInfo %>% 
  filter(!is.na(ClusterID)) %>%
  group_by(ClusterID) %>%
  summarise(ClMRCA=getMRCA(TreeMidRoot, TreeRepresentative))
# labelcolors<-qualitative_hcl(25, palette = "Dark 3")
# names(labelcolors)<- MultiCopiesForPlotF[order(MultiCopiesForPlotF$from,
#                                              MultiCopiesForPlotF$to),]$label

ClusterNodesAll$ClusterID<-factor(ClusterNodesAll$ClusterID,
                                  levels = c("Ia","Ib","II","III","IV","V","VI","VII","VIII","IX"))

clustercolors<-c("#8dd3c7","#ffffb3","#bebada","#fb8072",
                 "#80b1d3","#fdb462","#b3de69","#fccde5",
                 "#d9d9d9","#bc80bd")
names(clustercolors)<-c("Ia","Ib","II","III","IV","V","VI","VII","VIII","IX")
#I also want to add experimental strains to this plot
PreExpStrains<-subset(LeavesInfo, !is.na(LeavesInfo$Tmn_variant))
ExpStrains<-unique(PreExpStrains[,c("TreeRepresentative","Tmn_variant")])
#
InvertedTreeWithDupls<-ggtree(TreeMidRoot,
                              layout="inward_circular",
                              color = "#525252",
                              xlim = c(12,0)) + 
  geom_hilight(data = ClusterNodesAll,
               aes(node=ClMRCA,
                   fill=ClusterID),
               alpha=.3)+
  scale_fill_manual(values = clustercolors)+
  geom_taxalink(data = MultiCopiesForPlotF,
                mapping=aes(taxa1=from,taxa2=to,
                            #color = RepresentativeGenome,
                            hratio=size),
                #ncp=10,
                curvature = 0.5,
                color="#3690c0",
                offset =0.02) +
  guides(color="none")

InvertedTreeWithDuplsS<-InvertedTreeWithDupls%<+% ExpStrains+
  geom_tiplab2(aes(label=Tmn_variant),
               face="ArielMT",
               size=3,
               hjust=1)

ggsave("Tmn_multiple_copies_per_genome.pdf",
       plot = InvertedTreeWithDuplsS,
       path = "../figures/tree_figures/",
       width = 30,
       height=30,
       units="cm",
       dpi=300)

write.table(MultipleCopiesDf,
            file=paste0("../data/Tmn_multiple_copies_per_genome.tsv"),
            sep="\t",
            row.names = F,
            quote = F)

unique(MultipleCopiesDf$Proteinid)

############################################################
#####read fastas and look at alignments
TmnRepAlignment <- readAAStringSet("../data/Tmn_multiple/Tmn_padloc_1421_proteins.faa")
names(TmnRepAlignment) <- sub(" .*", "", names(TmnRepAlignment))
#############
#ploting function

getAliPlot<-function(aliobj,path_to_ali,chunk_size =100,zoomcoeff = 14)
{
  writeXStringSet(as(aliobj, "AAStringSet"), file = path_to_ali)
  alilength<-unique(width(as(aliobj, "AAStringSet")))
  #do plotting
  msa_plots <- list()
  final_plot<-list()
  i<-1
  for(start in seq(1,alilength, by=chunk_size)) {
    end <- min(start + chunk_size - 1, alilength)
    p <- ggmsa(path_to_ali, seq_name = TRUE, start = start, end = end) +
      geom_msaBar() 
    if ((start + chunk_size - 1) > alilength)
    {
      final_plot[[1]]<- as.ggplot(p)
    }
    else{
      msa_plots[[i]] <- as.ggplot(p)
    }
    
    i<-i+1
  }
  
  #combine all together
  fullchunk_plot<-ggarrange(plotlist = msa_plots,
                            ncol=1,align = "none")
  emptyplot<-ggplot() + theme_void()
  terminal_plot<-ggarrange(final_plot[[1]],
                           emptyplot,
                           widths = c(2,zoomcoeff), nrow =1)
  combined_plot<-ggarrange(fullchunk_plot,
                           terminal_plot,
                           heights = c(alilength%/%chunk_size,1),
                           ncol=1,
                           align = "v")
  
  return(combined_plot)
}

##get variable positions
GetVariablePositionsAli<-function(aln)
{
  aln_matrix <- as.matrix(aln)
  
  # Find columns where not all residues are identical
  variable_positions <- which(apply(aln_matrix, 2, function(col) {
    length(unique(col)) > 1
  }))
  
  return(variable_positions)
}
#############
##Case 1 Cluster 5
Clade5group<-subset(MultipleCopiesDf, MultipleCopiesDf$RefSeqGenomeID == "GCF_014596695.2")
Clade5ToKeep<-Clade5group$Proteinid
Clade5_aln_sub <- TmnRepAlignment[names(TmnRepAlignment) %in% Clade5ToKeep]
#realign pair
Clade5_realn <- msa(Clade5_aln_sub, method = "Muscle")

#####get basic stat
Clade5_AllChanges<-GetVariablePositionsAli(Clade5_realn)
Clade5ShArm<-seq(417,571)
Clade5LongArm<-seq(612,unique(width(as(Clade5_realn, "AAStringSet"))))
#calculate portion of substitutions per element
Clade5ChngOutsideArms<-length(Clade5_AllChanges[!(Clade5_AllChanges %in% c(Clade5ShArm,Clade5LongArm))])/(alilength - length(c(Clade5ShArm,Clade5LongArm)))
paste("Changes outside arms:",Clade5ChngOutsideArms)
Clade5ChngShArm<-length(Clade5_AllChanges[Clade5_AllChanges%in% Clade5ShArm])/length(Clade5ShArm)
paste("Changes in short arm:",Clade5ChngShArm)
Clade5ChngLongArm<-length(Clade5_AllChanges[Clade5_AllChanges%in% Clade5LongArm])/length(Clade5LongArm)
paste("Changes in long arm:",Clade5ChngLongArm)
#print list of residues for ChimeraX
toString(Clade5_AllChanges)

##

Clade5AliPlot<-getAliPlot(Clade5_realn,"../data/Tmn_multiple/Clade5_tmn.fasta",zoomcoeff =15)
ggsave(plot = Clade5AliPlot,
       filename = "Clade5_copies_alignment.pdf",
       path="../figures/Tmn_multiple/",
       height = 40,
       width = 36,
       dpi = 300,
       units = "cm")
       

#Case 2 Clade 4 neighboring leaves
Clade4group<-subset(MultipleCopiesDf, MultipleCopiesDf$RefSeqGenomeID == "GCF_019931755.1")
Clade4ToKeep<-Clade4group$Proteinid
Clade4_aln_sub <- TmnRepAlignment[names(TmnRepAlignment) %in% Clade4ToKeep]
Clade4_realn <- msa(Clade4_aln_sub, method = "Muscle")
#####get basic stat
Clade4_AllChanges<-GetVariablePositionsAli(Clade4_realn)
Clade4ShArm<-seq(424,599)
Clade4LongArm<-seq(634,unique(width(as(Clade4_realn, "AAStringSet"))))
#calculate portion of substitutions per element
Clade4ChngOutsideArms<-length(Clade4_AllChanges[!(Clade4_AllChanges %in% c(Clade4ShArm,Clade4LongArm))])/(alilength - length(c(Clade4ShArm,Clade4LongArm)))
paste("Changes outside arms:",Clade4ChngOutsideArms)
Clade4ChngShArm<-length(Clade4_AllChanges[Clade4_AllChanges%in% Clade4ShArm])/length(Clade4ShArm)
paste("Changes in short arm:",Clade4ChngShArm)
Clade4ChngLongArm<-length(Clade4_AllChanges[Clade4_AllChanges%in% Clade4LongArm])/length(Clade4LongArm)
paste("Changes in long arm:",Clade4ChngLongArm)
#print list of residues for ChimeraX
Clade4_AllChangesAdj<-ifelse(Clade4_AllChanges<1024, Clade4_AllChanges,
                             Clade4_AllChanges -2)
toString(Clade4_AllChangesAdj)

##plot alignment and save
Clade4AliPlot<-getAliPlot(Clade4_realn,"../data/Tmn_multiple/Clade4_tmn.fasta", zoomcoeff = 3.5)
ggsave(plot = Clade4AliPlot,
       filename = "Clade4_copies_alignment.pdf",
       path="../figures/Tmn_multiple/",
       height = 40,
       width = 36,
       dpi = 300,
       units = "cm")
# #Case 3 Clade 4 more distant leaves
# Clade42group<-subset(MultipleCopiesDf, MultipleCopiesDf$RefSeqGenomeID == "GCF_014645115.1")
# Clade42ToKeep<-Clade42group$Proteinid
# Clade42_aln_sub <- TmnRepAlignment[names(TmnRepAlignment) %in% Clade42ToKeep]
# Clade42_realn <- msa(Clade42_aln_sub, method = "Muscle")
# 
# Clade42AliPlot<-getAliPlot(Clade42_realn,"../data/Tmn_multiple/Clade4_tmn_GCF_014645115.1.fasta", zoomcoeff = 3.6)
# 
# ggsave(plot = Clade42AliPlot,
#        filename = "Clade4_GCF_014645115.1_copies_alignment.pdf",
#        path="../figures/Tmn_multiple/",
#        height = 40,
#        width = 36,
#        dpi = 300,
#        units = "cm")
