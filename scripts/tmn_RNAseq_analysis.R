library(Rsubread)
library(stringr)
#library(limma)
library(edgeR)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ape)
library(ggrepel)

setwd("/net/intdev/cbb01/garushyants/tmn_antitmn/20241210_RNAseq_data")

# buildindex(basename="BL21_AI_tmn_index",reference="reference_genome_yi/tmn_reference_multiline.fa")

conditions<-read.csv("tmn_list_of_samples", sep="/", header=F)
conditions$cond<-str_replace(conditions$V2, "_15","")
condition<-conditions$cond
#Trying to align reads
AlignRNAseqReads<-function(prefix)
{
  #prefix<-"tmn_T2"
  if( !file.exists(paste0(prefix,"_all.bam")))
  {
    align(index="BL21_AI_tmn_index",
          readfile1=paste0("20241210_fastp/",prefix,"_15_1.fq.gz"),
          readfile2=paste0("20241210_fastp/",prefix,"_15_2.fq.gz"),
          type="dna", 
          output_file=paste0(prefix,"_all.bam"),
          minFragLength=50,
          maxFragLength=600,
          nthreads=8)
  }
  
  if( file.exists(paste0(prefix,"_all.bam")))
  {
    ExpCounts<-featureCounts(files=paste0(prefix,"_all.bam"),
                             annot.ext="20241216_T2_genbank/tmn_reference_20241216.gff", 
                             isGTFAnnotationFile=TRUE,
                             GTF.featureType="CDS",GTF.attrType="ID",
                             isPairedEnd = T,
                             requireBothEndsMapped = T,
                             countMultiMappingReads = T,
                             nthreads=8)
    return(ExpCounts)
  }
  else
  {
    message("No BAM file")
  }

}

#Get counts for each sample
YFPT2<-AlignRNAseqReads(condition[1])
YFPControl<-AlignRNAseqReads(condition[2])
tmnT2<-AlignRNAseqReads(condition[3])
tmnControl<-AlignRNAseqReads(condition[4])

YFPCounts<-merge(YFPT2$counts, YFPControl$counts, 
      by = 'row.names', all = TRUE) 
TmnCounts<-merge(tmnT2$counts, tmnControl$counts, 
                 by = 'row.names', all = TRUE) 
AllCounts<-merge(TmnCounts, YFPCounts, 
                 by = c('Row.names'), all = TRUE) 
rownames(AllCounts)<-AllCounts$Row.names

AllAnnotation<-YFPControl$annotation

##get gene annotation
GFFAnnotationIn<-read.gff("20241216_T2_genbank/tmn_reference_20241216.gff",
                          GFF3=T)
get_attribute_by_name<-function(attribute){
  #attribute<-"product"
  tmpdf<-separate(GFFAnnotationIn, col=attributes, sep=paste0(attribute,"="),
                  into=c("p1","p2"))
  tmpdf<-separate(tmpdf,col=p2,sep=";",
                  into=c("attribute"),
                  extra="drop")
  return(tmpdf$attribute)
}
GFFAnnotationIn$product<- get_attribute_by_name("product")
GFFAnnotationIn$gene<- get_attribute_by_name("gene")
GFFAnnotationIn$ID<- get_attribute_by_name("ID")

###look how close all the samples are

x<-DGEList(counts=AllCounts[,c(2:5)],
           genes=AllAnnotation[,c("GeneID","Length")])
isexpr <- rowSums(cpm(x) > 10) >= 2 
x <- x[isexpr,]

y<-voom(x,plot=T)

plotMDS(y)
################
#Look at Tmn expression
AllAnnotationTmn<-subset(AllAnnotation,
                        AllAnnotation$Chr == "pACYC_pBAD_tmnECOR25_nativePromoter")



CountsTmn<-subset(AllCounts[,c(2:5)],
                 row.names(AllCounts) %in% AllAnnotationTmn$GeneID)
xTmn <- DGEList(counts=CountsTmn[,c(1,2)], genes=AllAnnotationTmn[,c("GeneID","Length")])

#FilterLowExpressed
isexpr <- rowSums(cpm(xTmn) > 10) >= 2 
xTmn <- xTmn[isexpr,]
RPKMTmn<-data.frame(rpkm(xTmn))

RPKMTmn$foldchange<-log2(RPKMTmn$tmn_T2_all.bam/RPKMTmn$tmn_control_all.bam)
################
#Look at T2 only
AllAnnotationT2<-subset(AllAnnotation,
                          AllAnnotation$Chr == "T2_in_house")



CountsT2<-subset(AllCounts[,c(2:5)],
                   row.names(AllCounts) %in% AllAnnotationT2$GeneID)

xT2 <- DGEList(counts=CountsT2[,c(1,3)], genes=AllAnnotationT2[,c("GeneID","Length")])

#FilterLowExpressed
isexpr <- rowSums(cpm(xT2) > 10) >= 2 
xT2 <- xT2[isexpr,]

#ConvertToRPKMs
RPKMT2<-data.frame(rpkm(xT2))

RPKMT2$foldchange<-log2(RPKMT2$tmn_T2_all.bam/RPKMT2$YFP_T2_all.bam)

###Plot correlation between YFP and tmn
T2RPKMP<-ggplot(RPKMT2, aes(x=log10(tmn_T2_all.bam),
                            y=log10(YFP_T2_all.bam)))+
  geom_point(aes(colour = foldchange))+
  geom_smooth(method='lm', se =T)+
  ylab("log10(RPKM T2 YFP)")+
  xlab("log10(RPKM T2 tmn)")+
  annotate(x=2, y=4.5, 
           label=paste("R = ", round(cor(log10(RPKMT2$tmn_T2_all.bam), 
                                         log10(RPKMT2$YFP_T2_all.bam)),2)), 
           geom="text", size=5)+
  theme_bw()
T2RPKMP
####Plot genes along genome
RPKMT2WithAnnotInit<-merge(RPKMT2,AllAnnotationT2,
                         by.x='row.names',
                         all.x =T,
                         by.y="GeneID")
RPKMT2WithAnnotFull<-merge(RPKMT2WithAnnotInit,
                           GFFAnnotationIn[,c("ID","gene","product")],
                           by.x="Row.names",
                           all.x=T,
                           by.y="ID")
#########
write.table(RPKMT2WithAnnotFull,
            file="T2_RPKM_along_genome_20241216.tsv",
            sep="\t",
            quote=F,
            row.names = F)

ForPlot<-RPKMT2WithAnnotFull
ForPlot$label<-ifelse(abs(ForPlot$foldchange) > 2 & (ForPlot$product !="hypothetical protein"), 
                      ForPlot$product, NA)
ForPlot$alpha<-ifelse(abs(ForPlot$foldchange) > 2, 1,.6)

T2RPKMAlongGenomePlot<-ggplot(data=ForPlot,
       aes(x=Start, 
           y=foldchange, 
           fill=foldchange))+
  geom_point(aes(alpha=alpha),
             shape=21,
             size=4)+
  geom_text_repel(aes(label = label),
            nudge_y = -0.2,
            segment.size =0.09,
            na.rm = T,
            min.segment.length = 0,
            max.overlaps = 10,
            size=3)+
  ylab("log2 fold change")+
  xlab("Genome position")+
 geom_hline(yintercept = c(2,-2),
            colour="#3690c0",
            linetype='dashed')+
  theme_bw()+
  theme(legend.position = "none",
        axis.text =element_text(size=12),
        axis.title = element_text(size=14))
T2RPKMAlongGenomePlot

ggsave("T2_RPKM_along_genome_20241216.png",
       plot=T2RPKMAlongGenomePlot,
       width =40,
       height =15,
       dpi=300,
       units = "cm")

ggsave("T2_RPKM_along_genome_20241216.pdf",
       plot=T2RPKMAlongGenomePlot,
       width =40,
       height =15,
       dpi=300,
       units = "cm")

#######################################################################
##################New part of RNAseq analysis (2024/12/28)
###########################################################
####Fold change along bacterial chromosome

AnnotationChr<-subset(AllAnnotation,
                        AllAnnotation$Chr == "BL21_AI")

CountsChr<-subset(AllCounts[,c(2:5)],
                 row.names(AllCounts) %in% AnnotationChr$GeneID)

xChr <- DGEList(counts=CountsChr, genes=AnnotationChr[,c("GeneID","Length")])

# #Filter Out Low Expressed
# isexpr <- rowSums(cpm(xChr) > 10) >= 2
# xChr <- xChr[isexpr,]

#I use here modified function because I am considering the sets of genes that I can compare for tmn+T2 vs YFP+T2
isexpr <- rowSums(cpm(xChr)[,c(1,3)] > 10) >= 2
xChr <- xChr[isexpr,]
#ConvertToRPKMs
RPKMChr<-data.frame(rpkm(xChr))
RPKMChr$OddsRatio<-((RPKMChr$tmn_T2_all.bam+.5)/(RPKMChr$tmn_control_all.bam+.5))/((RPKMChr$YFP_T2_all.bam+.5)/(RPKMChr$YFP_control_all.bam+.5))
# RPKMChr$tmnT2YFPT2<-(RPKMChr$tmn_T2_all.bam+.5)/(RPKMChr$YFP_T2_all.bam+.5)
RPKMChr$log2Ratio<-log2(RPKMChr$OddsRatio)

# RPKMChr$fctmnT2tmn<-log2((RPKMChr$tmn_T2_all.bam+.5)/(RPKMChr$tmn_control_all.bam+.5))
# RPKMChr$fcYFPT2YFP<-log2((RPKMChr$YFP_T2_all.bam+.5)/(RPKMChr$YFP_control_all.bam+.5))
# RPKMChr$fctmnT2YFPT2<-log2((RPKMChr$tmn_T2_all.bam+.5)/(RPKMChr$YFP_T2_all.bam+.5))
RPKMChr$fctmnYFP<-log2((RPKMChr$tmn_control_all.bam+.5)/(RPKMChr$YFP_control_all.bam+.5))
RPKMChrWithAnnot<-merge(RPKMChr,AnnotationChr,
                           by.x='row.names',
                           all.x =T,
                           by.y="GeneID")
RPKMChrWithAnnotFull<-merge(RPKMChrWithAnnot,
                           GFFAnnotationIn[,c("ID","gene","product")],
                           by.x="Row.names",
                           all.x=T,
                           by.y="ID")
RPKMChrWithAnnotFull$PredictionConfidence<-ifelse((RPKMChrWithAnnotFull$tmn_T2_all.bam < 1 |
                                                     RPKMChrWithAnnotFull$tmn_control_all.bam < 1 |
                                                     RPKMChrWithAnnotFull$YFP_T2_all.bam < 1 |
                                                     RPKMChrWithAnnotFull$YFP_control_all.bam < 1),
                                                  "Low","High")
########################
write.table(RPKMChrWithAnnotFull,
            file="RPKM_along_genome_BL21_AI_in_presense_of_T2_20241230.tsv",
            sep="\t",
            quote=F,
            row.names = F)
########################
RPKMChrForPlot<- subset(RPKMChrWithAnnotFull,
                       RPKMChrWithAnnotFull$PredictionConfidence == 'High')
RPKMChrForPlot$label<-ifelse(abs(RPKMChrForPlot$log2Ratio) > 2,
                             RPKMChrForPlot$gene, NA)
RPKMChrForPlot$color<-ifelse(abs(RPKMChrForPlot$log2Ratio) > 2, "#2b8cbe","#bdbdbd")

ChrRPKMAlongGenomePlot<-ggplot(data=RPKMChrForPlot,
                              aes(x=Start,
                                  y=log2Ratio,
                                  fill=color))+
  geom_point(shape=21,
             size=3)+
  geom_text(aes(label=label),
            check_overlap = T,
            nudge_x=0,
            nudge_y=.12)+
  ylab("log2 Odds ratio")+
  xlab("Genome position")+
  scale_fill_identity()+
  geom_hline(yintercept = c(2,-2),
             colour="#3690c0",
             linetype='dashed')+
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14))
ChrRPKMAlongGenomePlot

ggsave("RPKM_along_genome_BL21_AI_in_presense_of_T2_20241230.png",
       plot=ChrRPKMAlongGenomePlot,
       width =40,
       height =15,
       dpi=300,
       units = "cm")
ggsave("RPKM_along_genome_BL21_AI_in_presense_of_T2_20241230.pdf",
       plot=ChrRPKMAlongGenomePlot,
       width =40,
       height =15,
       dpi=300,
       units = "cm")

