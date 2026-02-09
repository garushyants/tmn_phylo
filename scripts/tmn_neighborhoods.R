library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(gggenes)
library(ggtree)
library(treeio)
library(ape)
library(phytools)
library(tidyr)
library(ggnewscale)
library(ggtreeExtra)
library(aplot)
library(scales)
library(readxl)
library(RColorBrewer)
library(colorspace)


#this parameter defines the length of the neighborhood to study
neighLength<-15000
###
mainpath<-"/Volumes/garushyants/tmn_antitmn/20260128_tmn_padloc_rebuild_tree-relaxed/"
setwd(mainpath)
###Set figures folder###
Figures<-"./genomic_neighborhoods"
if (!dir.exists(Figures)) {
  dir.create(Figures)
  cat("Folder created.\n")
} else {
  cat("Folder already exists.\n")
}
###
######################
#########Load and draw phylogenetic tree
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
##############
#Read info about leaves required for tree construction
ct <- rep("guess", 26)
ct[24] <- "text"
LeavesInfo<-read_xlsx("SupplementaryTable1_all_info.xlsx", col_types = ct)
######
LeavesRepresentativesInfo<-subset(LeavesInfo,
                                  LeavesInfo$RepresentativeGenome == "Y")
#Add desired borders around Tmn
LeavesRepresentativesInfo$LeftBorder<-ifelse((LeavesRepresentativesInfo$Start - neighLength) < 0,
                                             1,
                                             LeavesRepresentativesInfo$Start - neighLength)
LeavesRepresentativesInfo$RightBorder<-LeavesRepresentativesInfo$End + neighLength
##Also get info about prophages and plasmids
MGEinfo<-read.csv("Tmn_785_MGE_data_genomad.summary.tsv",
                  header=T, sep="\t")
colnames(MGEinfo)[2:3]<-c("Plasmid","Prophage")
####################################################################################
####################################################################################
###Now let's dig into genomic annotation

#############Part 1.Loading GFFs####################################################
GFFPath<-"./representative_genomes/gff"
GFFFiles<-list.files(pattern="\\_genomic.gff.gz$",
                     path = GFFPath)
setwd(paste(mainpath,GFFPath,sep="/"))
#Reading archived GFFs takes some time
GFFdata<-readr::read_tsv(GFFFiles, id="file_name", skip = 9, col_names = F)

#I select here both CDS and some non-coding genes
GFFdataCDS<-subset(GFFdata, !(GFFdata$X3 %in% c("gene","pseudogene","exon","region")) & !is.na(GFFdata$X2))
GFFdataCDS$GenomeID<-str_replace(GFFdataCDS$file_name,"_genomic.gff.gz","")
setwd(mainpath)

#######Select the required annotation fields
extractField<-function(pattern, column){
  #pattern<-'\\;product=([^;]+)\\;'
  #column<-TmnNeigGFF$V10
  r<-regexpr(pattern,column)
  out <- rep(NA,length(column))
  out[r!=-1] <- regmatches(column, r)
  out<-str_replace_all(out,";","")
  return(out)
}
#####
GFFdataCDS$X9<-paste0(GFFdataCDS$X9,";")#I do that in case the product is last record
GFFdataCDS$ID<-extractField('ID=[^-]+-([^;]+)\\;',GFFdataCDS$X9)
GFFdataCDS$ID<-str_split(GFFdataCDS$ID,"-",simplify = T, n = 5)[,2]
GFFdataCDS$gene<-extractField('\\;gene=([^;]+)\\;',GFFdataCDS$X9)
GFFdataCDS$gene<-str_replace(GFFdataCDS$gene,"gene=","")
GFFdataCDS$note<-extractField('\\;Note=([^;]+)\\;',GFFdataCDS$X9)
GFFdataCDS$note<-str_replace(GFFdataCDS$note,"Note=","")
GFFdataCDS$product<-extractField('\\;product=([^;]+)\\;',GFFdataCDS$X9)
GFFdataCDS$product<-str_replace(GFFdataCDS$product,"product=","")
names(GFFdataCDS)[2]<-"seqid"

######Extract parts of genomes in the vecinity of Tmn
GFFdataCDSwTmn<-merge(GFFdataCDS,
                      LeavesRepresentativesInfo[c(1:7,26:28)],
                      by.x = c("GenomeID","seqid"),
                      by.y = c("RefSeqGenomeID","Contig"))
####
GFFofInterest<-GFFdataCDSwTmn %>%
  filter(X4 >= LeftBorder,
         X5   <= RightBorder)

# ####save genes for PFAM run
# GFFIDsToSave<-unique(GFFofInterest[,c("GenomeID","ID")])
# #I get 25424 records that I save
# write.table(GFFIDsToSave, file = paste0(Figures,"/PFAM_GeneIDs_",neighLength,".tsv"),
#             sep="\t",quote = F, row.names = F, col.names = F)
#############Loading PFAM data######################
##Reading PFAM hmmscan output
PFAMpath<-"./PFAM_search_20260201/"
PFAMselectedFiles<-list.files(pattern="\\.csv$",
                              recursive = T,
                              path = PFAMpath)
setwd(paste0(mainpath,"/",PFAMpath))
PFAMdata<-readr::read_csv(PFAMselectedFiles, id="file_name")

PFAMdata$GenomeID<-str_replace(PFAMdata$file_name,"_pfamscan.csv","")
setwd(mainpath)

PFAMdata$PFAMID<-str_split_i(PFAMdata$hmm_acc,"\\.",1)

PFAMdatashort<-PFAMdata %>% group_by(GenomeID,seq_id) %>%
  summarise(PFAM_IDs = paste(sort(unique(PFAMID)), collapse = ";"),
            PHMM_names = paste(unique(hmm_name), collapse = ";"),
            .groups = "drop")
#############Merge GFF with PFAM
GFFwPFAMofInterest<-merge(GFFofInterest,
                          PFAMdatashort,
                          by.x = c("GenomeID","ID"),
                          by.y = c("GenomeID","seq_id"),
                          all.x=T)
#############Part 2.Loading annotation from Genomad###################################
#I do this annotation on the same selected proteoomes as PFAM search
GenomadAnnotation<-read.csv("./genomic_neighborhoods/annotate_genomaddb/GenomadAnnotateResults_20260201.tsv",
                            sep="\t", header=F)
#Filter results
GenomadAnnotationFiltered<-subset(GenomadAnnotation,
                                  GenomadAnnotation$V3 < 0.001 & #E-value filter as in default genomad
                                  GenomadAnnotation$V6 > 0.75) #V6 is target coverage
##
GenomadAnnotationMeta<-read.csv("./genomic_neighborhoods/annotate_genomaddb/genomad_metadata_v1.9.tsv", 
                                sep="\t")
GenomadAnnotationFull<-merge(GenomadAnnotationFiltered, 
                             GenomadAnnotationMeta,
                             by.x="V2",
                             by.y="MARKER",
                             all.x =T)
#This is the short data that I get out of it
collapse_safe <- function(x) {
  x <- sort(unique(na.omit(x)))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

GenomadAnnoCollapsed<-GenomadAnnotationFull %>% group_by(V1) %>%
  summarise(
    annotation_accessions        = collapse_safe(ANNOTATION_ACCESSIONS),
    annotation_descriptions      = collapse_safe(ANNOTATION_DESCRIPTION),
    genomad_hmm_ids               = collapse_safe(V2),
    genomad_specificity_class     = collapse_safe(SPECIFICITY_CLASS),
    .groups = "drop"
  )
########Merge with GFF and PFAM above
GFFwPGofInterest<-merge(GFFwPFAMofInterest,
                        GenomadAnnoCollapsed,
                          by.x = c("ID"),
                          by.y = c("V1"),
                          all.x=T)
#############Part 3.Loading PADLOC data###############################################
PADLOCdatapath<-"../20231218_padlocdb_2.0/defence_systems_tables/"
PADLOCfiles<-paste0(sort(unique(LeavesRepresentativesInfo$RefSeqGenomeID)),".csv")

setwd(paste0(mainpath,"/",PADLOCdatapath))

PADLOCdata<-readr::read_csv(PADLOCfiles, id="file_name")
PADLOCdata$GenomeID<-str_replace(PADLOCdata$file_name,".csv","")
setwd(mainpath)

##subset only records in the intervals of interest
PADLOCdatawTmn<-merge(PADLOCdata[,-1],
                      LeavesRepresentativesInfo[c(1:7,26:28)],
                      by.x = c("GenomeID","seqid"),
                      by.y = c("RefSeqGenomeID","Contig"))
###
PADLOCofInterest<-PADLOCdatawTmn %>%
  filter(start >= LeftBorder,
         end   <= RightBorder)
######################################################################################
#########################Merging annotations together#################################

####PADLOC and GFF
GFFwPADaGen<-merge(PADLOCofInterest[,c(1:26)], GFFwPGofInterest,
      by.x = c("GenomeID","seqid","start","end","TreeRepresentative",
               "Strand","ClusterID","Start","End","strand"),
      by.y = c("GenomeID","seqid","X4","X5","TreeRepresentative",
               "Strand","ClusterID","Start","End","X7"), all = T)

GFFwPADaGen$HMM_Annot_IDs<-ifelse(!is.na(GFFwPADaGen$annotation_accessions),
                                  GFFwPADaGen$annotation_accessions,
                                  ifelse(!is.na(GFFwPADaGen$PFAM_IDs),
                                         GFFwPADaGen$PFAM_IDs,
                                         NA))
GFFwPADaGen$product_description<-ifelse(!is.na(GFFwPADaGen$protein.name),
                                        GFFwPADaGen$protein.name,
                                        ifelse(!is.na(GFFwPADaGen$annotation_descriptions),
                                               GFFwPADaGen$annotation_descriptions,
                                               ifelse(is.na(GFFwPADaGen$product),
                                                      ifelse(is.na(GFFwPADaGen$note),
                                                             GFFwPADaGen$X3,
                                                             GFFwPADaGen$note),
                                                      GFFwPADaGen$product)))
####created combined annotation
GFFwPADaGen$PlotLabel<-ifelse(!is.na(GFFwPADaGen$protein.name),
                              GFFwPADaGen$protein.name,
                              ifelse(!is.na(GFFwPADaGen$gene),
                                     GFFwPADaGen$gene,
                                     ifelse(!is.na(GFFwPADaGen$PFAM_IDs),
                                            GFFwPADaGen$PFAM_IDs,
                                            ifelse(!is.na(GFFwPADaGen$annotation_descriptions),
                                                   GFFwPADaGen$annotation_accessions,
                                                   ifelse(GFFwPADaGen$X3 != "CDS",
                                                          ifelse(GFFwPADaGen$X3 =="sequence_feature",
                                                                 GFFwPADaGen$note,
                                                                 GFFwPADaGen$X3),
                                                          GFFwPADaGen$product)))))

####revert coordinates
GFFwPADaGen$LocalStart<-ifelse(GFFwPADaGen$Strand == "+",
                               GFFwPADaGen$start - (GFFwPADaGen$Start+(GFFwPADaGen$End-GFFwPADaGen$Start)/2) -1,
                               (GFFwPADaGen$Start+(GFFwPADaGen$End-GFFwPADaGen$Start)/2) - GFFwPADaGen$end +1)
GFFwPADaGen$LocalEnd<-ifelse(GFFwPADaGen$Strand == "+",
                             GFFwPADaGen$end - (GFFwPADaGen$Start+(GFFwPADaGen$End-GFFwPADaGen$Start)/2) +1,
                             (GFFwPADaGen$Start+(GFFwPADaGen$End-GFFwPADaGen$Start)/2) - GFFwPADaGen$start -1)
GFFwPADaGen$LocalStrand<-ifelse(GFFwPADaGen$Strand == "+", #Strand is for representative tmn
                                ifelse(GFFwPADaGen$strand == "+", #strand is for particular gene
                                       T,F),
                                ifelse(GFFwPADaGen$strand == "+",
                                       F,T))
#here I was before looking at the wrong column
GFFwPADaGen$Integrase<-grepl("integrase|recombinase", GFFwPADaGen$product_description, ignore.case=T)

GFFwPADaGen$fill<-ifelse(!is.na(GFFwPADaGen$system.number),
                        ifelse(startsWith(GFFwPADaGen$system,"PDC"),
                               "else",GFFwPADaGen$system),
                         ifelse(GFFwPADaGen$X3 != "CDS",
                                ifelse(GFFwPADaGen$X3 !="sequence_feature",
                                       GFFwPADaGen$X3,
                                       "else"),
                                ifelse(!is.na(GFFwPADaGen$genomad_specificity_class),
                                       ifelse(grepl("VV",GFFwPADaGen$genomad_specificity_class),
                                              "phage",
                                              ifelse(grepl("PP",GFFwPADaGen$genomad_specificity_class),
                                                     "plasmid",
                                                     "else")),
                                       ifelse(GFFwPADaGen$Integrase,
                                              "integrase","else"))))
###Getting rid of duplicates
GFFwPADaGen<-GFFwPADaGen[order(GFFwPADaGen$GenomeID,GFFwPADaGen$seqid,GFFwPADaGen$start),]
GFFwPADaGen$duplic<-ifelse((GFFwPADaGen$LocalEnd == lag(GFFwPADaGen$LocalEnd, default=-100) |
                                    GFFwPADaGen$LocalEnd == lead(GFFwPADaGen$LocalEnd, default=-100) |
                                    GFFwPADaGen$LocalStart == lag(GFFwPADaGen$LocalStart, default=-100) |
                                    GFFwPADaGen$LocalStart == lead(GFFwPADaGen$LocalStart, default=-100)) &
                                   is.na(GFFwPADaGen$system),
                                 F,T)
GFFwPADaGenNoDupl<-subset(GFFwPADaGen, GFFwPADaGen$duplic)
# #######Clean Table To Save
# ##Supplementary Table 2
# SupTable2<-GFFwPADaGenNoDupl[,c(5,7,1,2,27,10,3,4,12,30,34,45:48,53)]#60:62,67,51,52)]
# names(SupTable2)[c(4,9:10)]<-c("Contig","DefenseSystem","Feature_type")
# write.table(SupTable2, paste0(Figures,"/Supplementary_table_2_genomic_contexts_20k.tsv"),
#             sep="\t",
#             row.names = F,
#             quote = F)
#######Clean up for Plot
DataForPlotClean<-GFFwPADaGenNoDupl[,c(5,1,2,7,48:51,53)]
######################################################################################
Clades<-unique(na.omit(LeavesRepresentativesInfo$ClusterID))
cladecolors<-c("#8dd3c7","#ffffb3","#bebada","#fb8072",
                "#80b1d3","#fdb462","#b3de69","#fccde5",
                "#d9d9d9","#bc80bd")
names(cladecolors)<-c("Ia","Ib","II","III","IV","V","VI","VII","VIII","IX")
drawNeigborhoodPlot<-function(Clade)
{
  # ###Testing on clade Ia
  #Clade<-"Ia"#Clades[1]
  CladeLeavesRepresentativesInfo<-subset(LeavesRepresentativesInfo,
                                         LeavesRepresentativesInfo$ClusterID == Clade)
  subtree<-keep.tip(TreeMidRoot,CladeLeavesRepresentativesInfo$TreeRepresentative)
  subanno<-subset(DataForPlotClean, DataForPlotClean$ClusterID == Clade)
  
  ###Draw tree
  BasicTreePlot<-ggtree(subtree,
                        size=0.5,
                        color="#989898") +
    geom_hilight(node=getMRCA(subtree,unique(subanno$TreeRepresentative)),
                 fill=cladecolors[Clade],
                 alpha=.3,
                 align="right",
                 to.bottom=T)
  #BasicTreePlot
  ###get leaves order for subsequent plots
  ###It is essential for arranging all the neighborhoods plots with the tree
  leaf_order<-BasicTreePlot$data %>%
    filter(isTip) %>% arrange (y)
  ####Add MGE data
  subMGEinfo<-subset(MGEinfo,
                     MGEinfo$TreeRepresentative %in% subanno$TreeRepresentative)
  subMGElong<-gather(subMGEinfo,
                     key = "Location", value = "Prediction", Plasmid:Prophage)
  subMGElong$TreeRepresentative<-factor(subMGElong$TreeRepresentative, levels = leaf_order$label)
  HGTPLot<-ggplot(data = subMGElong, 
                  aes(y=TreeRepresentative,
                      x=Location,
                      fill = Prediction),
                  color="#bdbdbd")+
    geom_tile()+
    scale_fill_manual(values=c("#238443","#0570b0"), guide="none", na.value = "#ffffff")+
    theme_tree2()+
    theme(legend.position = 'none',
          axis.text.x = element_text(angle=90))
  #HGTPLot
  
  #########AddGenes
  ##setting color scheme
  allfillgroups<-unique(subanno$fill)
  #invariable colors
  essential_colors<-c("#377eb8","#4daf4a","#4daf4a","#ec7014","#e41a1c","#ffffff","#c6dbef")
  names(essential_colors)<-c("tmn","tmRNA","tRNA","phage","integrase","else","plasmid")
  #colors for various defense
  other_groups <- setdiff(allfillgroups, names(essential_colors))
  base_palette<-brewer.pal(9,"Set3")
  extended_palette <- colorRampPalette(base_palette)(length(other_groups))
  names(extended_palette)<-other_groups
  #combining palletes together
  combinedcolors<-c(essential_colors,extended_palette)
  
  #Ploting genes
  subanno$TreeRepresentative<-factor(subanno$TreeRepresentative, levels = leaf_order$label)
  
  GenesPlot<-ggplot(data = subanno,
                    aes(y =  TreeRepresentative,
                        xmin = LocalStart,
                        xmax = LocalEnd))+ 
    geom_hline(aes(yintercept =TreeRepresentative),
               linewidth =.5,
               color="#bdbdbd")+
    geom_gene_arrow(aes(fill = fill,
                        forward=LocalStrand))+
    geom_gene_label(aes(label=PlotLabel))+
    scale_fill_manual(values = combinedcolors, na.value="white", name ="")+
    theme_tree2()+
    theme(legend.position = "right")+
    guides(fill=guide_legend(ncol=1))
  
  NeighPlotToSave<-HGTPLot %>% insert_left(BasicTreePlot, width = 4) %>% insert_right(GenesPlot, width=30) 
  ggsave(file = paste0("Clade_",Clade,"_genomic_neighborhood_",neighLength,".pdf"),
         path = Figures,
         plot = NeighPlotToSave,
         width=50,
         height = ifelse(length(leaf_order$label)/2 <=10,
                         12,
                         length(leaf_order$label)/2),
         limitsize = F,
         units="cm",
         dpi=300)
  #return(NeighPlotToSave)
}
for(cl in Clades){
  print(cl)
  drawNeigborhoodPlot(Clade = cl)
}


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
  scale_fill_manual(values = c("#8dd3c7","#ffffb3","#bebada","#fb8072",
                               "#80b1d3","#fdb462","#b3de69","#fccde5",
                               "#d9d9d9","#bc80bd"))+
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
       path = Figures,
       width = 30,
       height=30,
       units="cm",
       dpi=300)

write.table(MultipleCopiesDf,
            file=paste0(Figures,"/Tmn_multiple_copies_per_genome.tsv"),
            sep="\t",
            row.names = F,
            quote = F)

unique(MultipleCopiesDf$Proteinid)

# ##########################################################
# ##########################################################
# ###Looking at the genomes with more than one Tmn available
# ###Starting from the ones for which I already have neighborhoods
# CountTmns845<-TmnCoordFP %>% group_by(GenomeID) %>%
#   summarise(count=n())
# TmnNotSingl<-subset(CountTmns845,
#                     CountTmns845$count >1)
# GenomesOfInterest<-unique(TmnNotSingl$GenomeID)
# 
# GenesToPlotNotSingl<-subset(GenesToPlotNC,
#                             GenesToPlotNC$GenomeID %in% GenomesOfInterest)
# 
# DublTmnsPlot<-ggplot(data = GenesToPlotNotSingl,
#        aes(y =  molecule,
#            xmin = pstart,
#            xmax = pend))+ 
#   geom_hline(aes(yintercept = molecule, color = GenomeID))+
#   geom_gene_arrow(aes(fill = fill,
#                       forward=pstrand))+
#   geom_gene_label(aes(label=Agene))+
#   scale_fill_manual(values = GenesColors, na.value="white")+
#   theme_tree2()+
#   theme(legend.position = "right")+
#   guides(fill=guide_legend(ncol=1))
# 
# Tmn3GenomesDupl<-HGTPLot %>% insert_left(BasicTreePlotWithLabels, width = 4) %>% insert_right(DublTmnsPlot, width=15)
# Tmn3GenomesDupl
# 
# ggsave("Tmn_3genomes_dupl_Tree_with_neighborhoods.pdf",
#        plot=Tmn3GenomesDupl,
#        path=Figures,
#        width=30, height=30,
#        limitsize = F,
#        units="in",
#        dpi=300)
# 
# #############################################
# ####After that I can work with all duplicates
# ###For that I need full Tmn table from PADLOC
# AllTmnsPADLOC<-read.csv("./20251017_tmn_neighborhoods/tmn_padloc20_in_refseq_filtered_hmmcov0.8.noPseudogenes.rfriendly.csv",
#          header = F)
# #subset to only the ones on the tree
# #but actually the data in the previous file is already filtered
# TmnsPADLOCOnTree<-subset(AllTmnsPADLOC,
#                          AllTmnsPADLOC$V5 %in% PADLOCTmn$target.name)
# ###Find genomes that have more than one copy
# TmnsPADLOCOnTreeCounts<-TmnsPADLOCOnTree %>% group_by(V1) %>%
#   summarize(count=n())
# TmnAllDupl<-subset(TmnsPADLOCOnTreeCounts,
#                    TmnsPADLOCOnTreeCounts$count > 1)
# TmnAllDuplProteins<-subset(TmnsPADLOCOnTree,TmnsPADLOCOnTree$V1 %in% TmnAllDupl$V1)
# ###
# #Load assembly summary
# AllTmnsAssemblySummary<-readr::read_tsv("./20251017_tmn_neighborhoods/tmn_padloc_genome_with_assembly_summary.tsv",
#                                  col_names = F)
# GenomesOfInterestSummary<-subset(AllTmnsAssemblySummary,
#                                  AllTmnsAssemblySummary$X1 %in% TmnAllDupl$V1)
# GenomesOfInterestSummary$genus<-str_split_fixed(GenomesOfInterestSummary$X8," ",2)[,1]
# #1,8,9,39
# 
# TmnAllDuplProteinsWithTax<-merge(TmnAllDuplProteins[,c(1,3,5)],
#                                  GenomesOfInterestSummary[,c(1,8,9,39)],
#                                  by.x ="V1",
#                                  by.y ="X1")
# #MergeWithCl
# TmnAllDuplProteinsWithTaxWithCl<-merge(TmnAllDuplProteinsWithTax,
#                                        ClustersOfInterestDf,
#                                        by.x="V5", by.y = "SequenceName", all.x =T)
# #####Those only belong to clusters Ia, Ib and III
# ##Subseting Tree
# LeavesToKeep<-subset(ClustersOfInterestDf,
#                      ClustersOfInterestDf$ClusterID %in% c('Ia','Ib',"II",'III'))
# 
# TmnSubTree<-keep.tip(TreeMidRoot,LeavesToKeep$SequenceName)
# #Get MRCA for clade IDs
# SubTreeClusterAncestryNodeOfInterest<-LeavesToKeep %>% group_by(ClusterID) %>%
#   summarise(ClMRCA=getMRCA(TmnSubTree, SequenceName))
# 
# TmnSubTreePlot<-ggtree(TmnSubTree,
#                       size=0.2,
#                       color="#636363") + geom_hilight(data =SubTreeClusterAncestryNodeOfInterest,
#                                                     aes(node=ClMRCA,
#                                                         fill=ClusterID),
#                                                     alpha=.3,
#                                                     align="right",
#                                                     to.bottom=T,
#                                                     extend=.05)+
#   scale_fill_manual(values = c("#8dd3c7","#ffffb3","#bebada","#fb8072"))
# TmnSubTreePlot
# 
# #get order
# SubTreeleaf_order<-TmnSubTreePlot$data %>%
#   filter(isTip) %>% arrange (y)
# names(TmnAllDuplProteinsWithTaxWithCl)[1]<-"molecule"
# ##adding other leaves
# TmnAllDuplProteinsWithTaxWithClFull<-merge(TmnAllDuplProteinsWithTaxWithCl,
#                                            data.frame(molecule = LeavesToKeep[,c(2)]), by = "molecule",
#                                            all = T)
# TmnAllDuplProteinsWithTaxWithClFull$V1<-ifelse(is.na(TmnAllDuplProteinsWithTaxWithClFull$V1),
#                                                "",TmnAllDuplProteinsWithTaxWithClFull$V1)
# 
# TmnAllDuplProteinsWithTaxWithClFull$molecule<-factor(TmnAllDuplProteinsWithTaxWithClFull$molecule, levels = SubTreeleaf_order$label)
# 
# AllDuplPlot<-ggplot(data = TmnAllDuplProteinsWithTaxWithClFull, aes(y = molecule,x = V1,fill=genus))+
#   geom_tile()+ 
#   scale_y_discrete(position="right") +
#   scale_fill_manual(values = c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00"),
#                     na.value = "white",
#                     na.translate = FALSE)+
#   theme_tree2() + 
#   theme(axis.text.x = element_text(angle=90)) +
#   xlab(NULL)+ylab(NULL)
# AllDuplPlot
# 
# AllDuplicatesToSave<-AllDuplPlot %>% insert_left(TmnSubTreePlot)
# AllDuplicatesToSave
# 
# 
# ggsave("Tmn_All_dupl_with_Tree.pdf",
#        plot=AllDuplicatesToSave,
#        path=Figures,
#        width=10, height=10,
#        limitsize = F,
#        units="in",
#        dpi=300)
