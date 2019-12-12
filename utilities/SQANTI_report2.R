#####################################
##### SQANTI report generation ######
#####################################

### Author: Lorena de la Fuente
### Date: 14/10/2016

#********************** Taking arguments from python script

args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1]
junc.file <- args[2]

if (length(args)>2) {
  param.file <- args[3]
}

report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1];
report.file <- paste(report.prefix, "sqanti_report.pdf", sep="_");
class.file2 <- paste(report.prefix, "_classification_TPM.txt", sep='');


#********************** Packages (install if not found)

list_of_packages <- c("ggplot2", "scales", "reshape", "gridExtra", "grid", "dplyr")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

library(ggplot2)
library(scales)
library(reshape)
library(gridExtra)
library(grid)
library(dplyr)

# ***********************
# PLOTS
# p0: Splicing complexity (X) Isoforms per Gene (Y) Number of Genes
# p1: Distribution of Isoform Classification
# p2: Distribution of Ref Lengths (FSM ISM only)
# p3: Distribution of Ref Exons   (FSM ISM only)
# p4: Distribution of Isoform Lengths, by Classification
# p5: Distribution of Exon Counts, by Classification
# p6: Distribution of Mono- vs Multi-Exons, Novel vs Annotated Genes
# p7: Splicing complexity, Novel vs Annotated Genes
# p.classByLen.a: Structural categories with increasing transcript length, absolute
# p.classByLen.b: Structural categories with increasing transcript length, normalized
#
# p21.a: Distance to polyA site for FSM, absolute
# p21.b: Distance to polyA site for FSM, percentage
# p21.dist3.ISM.a: Distance to polyA site for ISM, absolute
# p21.dist3.ISM.b: Distance to polyA site for ISM, percentage
# p22.a: Distance to start site for FSM, absolute
# p22.b: Distance to start site for FSM, percentage

# p23.a: Splice Junctions by Classification (known canonical, known non-canonical, novel canonical, novel non-canonical)
# p23.b: Splice Junctions by Classification (canonical vs non-canonical)
#
# p29.a: Splice Junction, % of RT switching, all junctions
# p29.b: Splice Junction, % of RT switching, unique junctions
#
# p30: intra-priming, by Classification
# p31: intra-priming, Mono- vs Multi-Exons
# p32: intra-priming, Coding vs Non-Coding
# ***********************

#********************** Reading Information

########## Classification information

data.class = read.table(class.file, header=T, as.is=T, sep="\t")
rownames(data.class) <- data.class$isoform

# (Liz) not sorting by expression
#if (!all(is.na(data.class$iso_exp))){
#  sorted <- data.class[order(data.class$iso_exp, decreasing = T),]
#  FSMhighestExpIsoPerGene <- sorted[(!duplicated(sorted$associated_gene) & sorted$structural_category=="full-splice_match"),"isoform"]
#  data.class[which(data.class$isoform%in%FSMhighestExpIsoPerGene),"RTS_stage"] <- FALSE
#  write.table(data.class, file=class.file, row.names=FALSE, quote=F, sep="\t")
#}

xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")

legendLabelF1 <- levels(as.factor(data.class$coding));

data.class$structural_category = factor(data.class$structural_category,
                                       labels = xaxislabelsF1, 
                                       levels = xaxislevelsF1,
                                       ordered=TRUE)

data.FSMISM <- subset(data.class, structural_category %in% c('FSM', 'ISM'))
data.FSM <- subset(data.class, (structural_category=="FSM" & exons>1))
data.ISM <- subset(data.class, (structural_category=="ISM" & exons>1))


########### Junction information

data.junction <- read.table(junc.file, header=T, as.is=T, sep="\t")

# (Liz) don't sort junction file by expression
#if (!all(is.na(data.class$iso_exp))){
#  data.junction[which(data.junction$isoform%in%FSMhighestExpIsoPerGene),"RTS_junction"] <- FALSE
#  write.table(data.junction, file=junc.file, row.names=FALSE, quote=F, sep="\t")
#}

# make a unique identifier using chrom_strand_start_end
data.junction$junctionLabel = with(data.junction, paste(chrom, strand, genomic_start_coord, genomic_end_coord, sep="_"))

data.junction$SJ_type <- with(data.junction, paste(junction_category,canonical,"SJ", sep="_"))
data.junction$SJ_type <- factor(data.junction$SJ_type, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                       labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nNon-canonical "), order=T)

data.junction$structural_category = data.class[data.junction$isoform, "structural_category"]

uniqJunc <- unique(data.junction[,c("junctionLabel", "SJ_type", "total_coverage")]);
uniqJuncRTS <- unique(data.junction[,c("junctionLabel","SJ_type", "RTS_junction")]);


########## Generating plots


#*** Global plot parameters

myPalette = c("#6BAED6","#FC8D59","#78C679","coral2","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")

mytheme <- theme_classic(base_family = "Palatino") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=14),
        axis.text.x  = element_text(size=13),
        axis.title.y = element_text(size=14),
        axis.text.y  = element_text(vjust=0.5, size=13) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=11), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=13)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm")) 



# Create a new attribute called "novelGene"

data.class$novelGene <- "Annotated Genes"
data.class[grep("novelGene", data.class$associated_gene), "novelGene"] <- "Novel Genes"
data.class$novelGene = factor(data.class$novelGene,
                              levels = c("Novel Genes","Annotated Genes"),
                              ordered=TRUE)

# Create a new attribute called "exonCat"

data.class[which(data.class$exons>1), "exonCat"] <- "Multi-Exon"
data.class[which(data.class$exons==1), "exonCat"] <- "Mono-Exon"
data.class$exonCat = factor(data.class$exonCat,
                      levels = c("Multi-Exon","Mono-Exon"),
                      ordered=TRUE)

data.class$all_canonical = factor(data.class$all_canonical,
                                  levels = c("canonical","non_canonical"),
                                  ordered=TRUE)

# ----------------------------------------------------------
# Make "isoPerGene" which is aggregated information by gene
#  $associatedGene - either the ref gene name or novelGene_<index>
#  $novelGene      - either "Novel Genes" or "Annotated Genes"
#  $FSM_class      - "A", "B", or "C"
#  $geneExp        - gene expression info
#  $nIso           - number of isoforms associated with this gene
#  $nIsoCat        - splicing complexity based on number of isoforms
# ----------------------------------------------------------

if (!all(is.na(data.class$gene_exp))){
  isoPerGene = aggregate(data.class$isoform,
                         by = list("associatedGene" = data.class$associated_gene,
                                   "novelGene" = data.class$novelGene,
                                   "FSM_class" = data.class$FSM_class,
                                   "geneExp"=data.class$gene_exp),
                         length)
} else {
  isoPerGene = aggregate(data.class$isoform,
                         by = list("associatedGene" = data.class$associated_gene,
                                   "novelGene" = data.class$novelGene,
                                   "FSM_class" = data.class$FSM_class),
                         length)
}
# assign the last column with the colname "nIso" (number of isoforms)
colnames(isoPerGene)[ncol(isoPerGene)] <- "nIso"


isoPerGene$FSM_class2 = factor(isoPerGene$FSM_class, 
                               levels = c("A", "B", "C"), 
                               labels = c("MonoIsoform Gene", "MultiIsoform Genes\nwithout expression\nof a FSM", "MultiIsoform Genes\nexpressing at least\none FSM"), 
                               ordered=TRUE)

isoPerGene$novelGene = factor(isoPerGene$novelGene, 
                              levels = c("Annotated Genes", "Novel Genes"), 
                              ordered=TRUE)

isoPerGene$nIsoCat =cut(isoPerGene$nIso, breaks = c(0,1,3,5,max(isoPerGene$nIso)+1), labels = c("1", "2-3", "4-5", ">=6"))


# single FL count file provided
if (!all(is.na(data.class$FL))){
    total_fl <- sum(data.class$FL, na.rm=T)
    data.class$FL_TPM <- round(data.class$FL*(10**6)/total_fl)
}


# see if there are multple FL samples
FL_multisample_indices <- which(substring(colnames(data.class), 1,3)=="FL.")
if (length(FL_multisample_indices)>0)
{
    FL_multisample_names <- substring(colnames(data.class)[FL_multisample_indices],4)
    FL_TPM_multisample_names <- c();

    for (i in 1:length(FL_multisample_indices)) {
        j <- FL_multisample_indices[i]
        name <- paste("FL_TPM", FL_multisample_names[i], sep='.')
        name2 <- paste(name, "_log10", sep='')
        FL_TPM_multisample_names <- c(FL_TPM_multisample_names, name)
        total_fl <- sum(data.class[j])
        data.class[,name] <- data.class[j]*(10**6)/total_fl
        data.class[,name2] <- log10(data.class[j]*(10**6)/total_fl + 1)
    }
    write.table(data.class, class.file2, quote=F, sep='\t');
}


# PLOT length of isoforms
# p.length.all: length of all isoforms, regardless of category
# p.length.cat: length of isoforms, by category
# p.length.exon: length of isoforms, mono- vs mult-exon
# (optional) p.length.all.sample
# (optional) p.length.exon.sample

p.length.all <- ggplot(data.class, aes(x=length)) +
  geom_histogram(binwidth=100) +
  guides(fill=FALSE) +
  scale_fill_manual(values = myPalette[c(2:5)]) +
  labs(x="Transcript Length", y="Count", title="Transcript Lengths, all transcripts") +
  theme(legend.position="top") +
  mytheme

if (length(FL_multisample_indices)>0) {  # has multiple samples
    df.length_by_sample <- data.frame();
    for (i in 1:length(FL_multisample_indices)) {
        j <- FL_multisample_indices[i];
        df_new <- data.class[data.class[j]>0,c("isoform", "length", "exonCat")];
        df_new$sample <- FL_multisample_names[i];
        df.length_by_sample <- rbind(df.length_by_sample, df_new);
    }

    p.length.all.sample <- ggplot(df.length_by_sample, aes(x=length, color=sample)) +
        geom_freqpoly(binwidth=100) +
        guides(fill=FALSE) +
        scale_fill_manual(values = myPalette[c(2:5)]) +
        labs(x="Transcript Length", y="Count", title="Transcript Lengths, By Sample") +
        theme(legend.position="top") +
        mytheme

    p.length.exon.sample <- ggplot(df.length_by_sample, aes(x=length, color=sample, lty=exonCat)) +
        geom_freqpoly(binwidth=100) +
        guides(fill=FALSE) +
        scale_fill_manual(values = myPalette[c(2:5)]) +
        labs(x="Transcript Length", y="Count", title="Transcript Lengths, Mono- vs Multi-Exons, By Sample") +
        theme(legend.position="top") +
        mytheme
}

p.length.cat <- ggplot(data.class, aes(x=length, color=structural_category)) +
  geom_freqpoly(binwidth=100) +
  guides(fill=FALSE) +
  scale_fill_manual(values = myPalette[c(2:5)]) +
  labs(x="Transcript Length", y="Count", title="Transcript Lengths, by structural category") +
  theme(legend.position="top") +
  mytheme

p.length.exon <- ggplot(data.class, aes(x=length, color=exonCat)) +
  geom_freqpoly(binwidth=100) +
  guides(fill=FALSE) +
  scale_fill_manual(values = myPalette[c(2:5)]) +
  labs(x="Transcript Length", y="Count", title="Transcript Lengths, Mono- vs Multi-Exons") +
  theme(legend.position="top") +
  mytheme


# p0: Distribution of Number of Isoforms

p0 <- ggplot(isoPerGene, aes(x=nIsoCat, fill=nIsoCat)) +
  geom_bar(stat="count", aes(y= (..count..)/sum(..count..)), color="black", size=0.3, width=0.7) +
  guides(fill=FALSE) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(2:5)]) +
  labs(x ="Isoforms Per Gene", title="Number of Isoforms per Gene\n\n\n", y = "% Genes") +
  mytheme

# p7: Distribution of Number of Isoforms, separated by Novel vs Annotated Genes

p7 <- ggplot(data=isoPerGene, aes(x=novelGene)) +
  geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=nIsoCat), color="black", size=0.3, width=0.5) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(name = "Isoforms Per Gene",
                    values = myPalette[c(2:5)]) +
  ylab("% Genes ") +
  xlab("Gene Type") +
  mytheme +
  theme(axis.title.x=element_blank()) +
  theme(legend.position="bottom") +
  guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
  ggtitle("Number of Isoforms per Gene, Novel vs Known Geness\n\n\n\n" )


#**** PLOT 1: Structural Classification

p1 <- ggplot(data=data.class, aes(x=structural_category)) +
  geom_bar(aes(y = (..count..)/sum(..count..), alpha=coding, fill=structural_category), color="black", size=0.3, width=0.7) +
  scale_y_continuous(labels = percent, expand = c(0,0), limits = c(0,1)) +
  geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25)  +
  scale_x_discrete(drop=FALSE) +
  scale_alpha_manual(values=c(1,0.3), 
                     name = "Coding prediction", 
                     labels = legendLabelF1)+
  scale_fill_manual(values = myPalette, guide='none') + 
  xlab("") + 
  ylab("% Transcripts") +
  mytheme + 
  geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Isoform distribution across structural categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1)) 


#**** PLOTS 2-3: refLength and refExons for ISM and FSM transcripts. Plot if any ISM or FSM transcript

if (nrow(data.FSMISM) > 0) {

  p2 <- ggplot(data=data.FSMISM, aes(x=structural_category, y=ref_length/1000, fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
    scale_fill_manual(values = myPalette) +
    scale_x_discrete(drop = TRUE) +
    guides(fill=FALSE) +
    xlab("") +  
    ylab("Matched Reference Length (in kb)") +
    labs(title="Length Distribution of Matched Reference Transcripts\n\n\n",
         subtitle="Applicable only to FSM and ISM categories\n\n")

  p3 <- ggplot(data=data.FSMISM, aes(x=structural_category, y=ref_exons, fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +   
    scale_x_discrete(drop = TRUE) +
    xlab("") +  
    ylab("Matched reference exon count") +
    scale_fill_manual(values = myPalette) +
    guides(fill=FALSE) +
    mytheme +
    labs(title="Exon Count Distribution of Matched Reference Transcripts\n\n\n",
         subtitle="Applicable only to FSM and ISM categories\n\n")

}


#****  PLOT 4: Transcript lengths by category

p4 <- ggplot(data=data.class, aes(x=structural_category, y=length, fill=structural_category)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  scale_x_discrete(drop=FALSE) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  ggtitle("Transcript Lengths by Structural Classification\n\n" ) +
  theme(axis.title.x=element_blank()) 


##**** PLOT 5: Exon counts by category

p5 <- ggplot(data=data.class, aes(x=structural_category, y=exons, fill=structural_category)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("Number of exons") +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(values = myPalette) +
  guides(fill=FALSE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) + 
  theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
  ggtitle("Exon Counts by Structural Classification\n\n" ) +
  theme(axis.title.x=element_blank())


##**** PLOT 6: Mono vs Multi-exon distribution for Known vs Novel Genes

p6 <- ggplot(data=data.class, aes(x=novelGene)) +
  geom_bar(position="fill",aes(y = (..count..)/sum(..count..), fill=exonCat), color="black", size=0.3, width=0.5) +
  scale_x_discrete(drop=FALSE) +
  scale_y_continuous(labels = percent, expand = c(0,0)) +
  scale_fill_manual(name = "Transcript type", 
                    values = myPalette[c(2:5)]) +
  ylab("% Transcripts ") +  
  mytheme +
  theme(axis.title.x=element_blank()) + 
  theme(legend.position="bottom") +   
  ggtitle("Distribution of Mono- vs Multi-Exon Transcripts\n\n" )



##**** PLOT  absolute and normalized % of different categories with increasing transcript length
# requires use of dplyr package
data.class$lenCat <- as.factor(as.integer(data.class$length %/% 1000))
data.class.byLen <- data.class %>% group_by(lenCat, structural_category) %>% summarise(count=n()) %>% mutate(perc=count/sum(count))
data.class.byLen$structural_category <- factor(data.class.byLen$structural_category, levels=rev(xaxislabelsF1), order=TRUE)

p.classByLen.a <- ggplot(data.class.byLen, aes(x=lenCat, y=count, fill=factor(structural_category))) +
    geom_bar(stat='identity') +
    labs(x="Transcript Length (in kb)", y="Percentages", title="Classifications by Transcript Length")


p.classByLen.b <- ggplot(data.class.byLen, aes(x=lenCat, y=perc*100, fill=factor(structural_category))) +
    geom_bar(stat='identity') +
    labs(x="Transcript Length (in kb)", y="Percentages", title="Classifications by Transcript Length, normalized")



##**** PLOT 8: Expression, if isoform expression provided (iso_exp is in TPM)
if (!all(is.na(data.class$iso_exp))){
  p8 <- ggplot(data=data.class, aes(x=structural_category, y=log2(iso_exp+1))) +
    geom_violin(aes(fill=structural_category), draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_x_discrete(drop=FALSE) +
    ylab("log2(TPM+1)") +
    scale_fill_manual(values = myPalette) +
    guides(fill=FALSE) +
    mytheme  + theme(axis.text.x = element_text(angle = 45)) + 
    theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Transcript Expression by Structural Category\n\n" )
}


# PLOT 9: FL number, if FL count provided
# convert FL count to TPM
if (!all(is.na(data.class$FL))){
    p9 <- ggplot(data=data.class, aes(x=structural_category, y=log2(FL_TPM+1))) +
    geom_violin(aes(fill=structural_category), draw_quantiles = c(0.25, 0.5, 0.75)) +
    ylab("log2(FL_TPM+1)") +
    scale_x_discrete(drop=FALSE) +
    scale_fill_manual(values = myPalette) +
    guides(fill=FALSE) +
    mytheme +
    theme(axis.text.x = element_text(angle = 45)) + 
    theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) +
    ggtitle("FL Count (normalized) by Structural Category\n\n" )
}



# PLOT 10: Gene Expression, if expresion provided
if (!all(is.na(data.class$iso_exp))){
  p10 <- ggplot(data=isoPerGene, aes(x=novelGene, y=log2(geneExp+1))) +
    geom_violin(aes(fill=novelGene), draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_x_discrete(drop=FALSE) +
    xlab("Structural Classification") +  
    ylab("log2(Gene_TPM+1)") +
    scale_fill_manual(values = myPalette[c(3:4)]) +
    guides(fill=FALSE) +
    mytheme +
    theme(axis.title.x=element_blank()) +
    ggtitle("Gene Expression, Annotated vs Novel\n\n" )
}


# PLOT 11: Gene FL number, if FL count provided
if (!all(is.na(data.class$FL))){
    FL_gene <- aggregate(as.integer(data.class$FL), by = list("associatedGene" = data.class$associated_gene), sum)
    colnames(FL_gene)[ncol(FL_gene)] <- "FL_gene"
    isoPerGene <- merge(isoPerGene, FL_gene, by="associatedGene")
    isoPerGene$FL_gene_TPM <- isoPerGene$FL_gene*(10**6)/total_fl

    p11 <- ggplot(data=isoPerGene, aes(x=novelGene, y=log2(FL_gene_TPM+1))) +
    geom_violin(aes(fill=novelGene), draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_x_discrete(drop=FALSE) +
    ylab("log2(FL_TPM+1)") +
    scale_fill_manual(values = myPalette[c(3:4)]) +
    guides(fill=FALSE) +
    mytheme +
    theme(axis.title.x=element_blank()) +
    ggtitle("Number of FL reads per Gene by type of gene annotation\n\n" )
}



# PLOT 12: NNC expression genes vs not NNC expression genes
# NNC expression genes vs not NNC expression genes

if (!all(is.na(data.class$gene_exp))){
  if (nrow(data.class[data.class$structural_category=="NNC",])!=0){
    
    NNC_genes <- unique(data.class[data.class$structural_category=="NNC","associated_gene"])
    notNNC_genes <- unique(data.class[!data.class$associated_gene%in%NNC_genes,"associated_gene"])
    isoPerGene[isoPerGene$associatedGene %in% notNNC_genes, "NNC_class"] <- "Genes without\n NNC isoforms"
    isoPerGene[isoPerGene$associatedGene %in% NNC_genes, "NNC_class"] <- "Genes with\n NNC isoforms"
    
    isoPerGene$NNC_class <- factor(isoPerGene$NNC_class, levels=c("Genes with\n NNC isoforms","Genes without\n NNC isoforms"),
                            labels=c("Genes with\n NNC isoforms","Genes without\n NNC isoforms"), order=T)
    
    p12 <- ggplot(data=isoPerGene[!is.na(isoPerGene$NNC_class),], aes(x=NNC_class, y=log2(geneExp+1))) +
      geom_violin(aes(fill=NNC_class), draw_quantiles=c(0.25, 0.5, 0.75)) +
      xlab("") +  
      ylab("log2(Gene_TPM+1)") +
      scale_x_discrete(drop=FALSE) +
      scale_fill_manual(values = c(myPalette[4],"grey38")) +
      guides(fill=FALSE) +
      mytheme +
      theme(axis.title.x=element_blank()) + 
      ggtitle("Gene Expression between NNC and not NNC containing Genes\n\n" )
    }
}

# ToDO: deal with expression data later
# (Liz) not plotting p13, p13.c for now
# PLOT 13: Genes expression to only FSM Genes, only NNC Genes and both containing genes

if (!all(is.na(data.class$gene_exp))){
  if (nrow(data.class[data.class$structural_category=="NNC",])!=0 & nrow(data.class[data.class$structural_category=="FSM",])!=0 ){
    FSM_just_genes = unique(data.class[data.class$FSM_class=="A" & data.class$structural_category=="FSM","associated_gene"])
    NNC_just_genes = unique(data.class[data.class$FSM_class=="A" & data.class$structural_category=="NNC","associated_gene"])
    FSMandNNCgenes = unique(data.class[data.class$FSM_class=="C" & data.class$structural_category=="NNC","associated_gene"])
    isoPerGene[isoPerGene$associatedGene %in% FSMandNNCgenes, "FSM_NNC_class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
    isoPerGene[isoPerGene$associatedGene %in% NNC_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only NNC isoforms"
    isoPerGene[isoPerGene$associatedGene %in% FSM_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only FSM isoforms"
    data.class[data.class$associated_gene %in% FSMandNNCgenes, "class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
    data.class[data.class$associated_gene %in% NNC_just_genes, "class"] <- "Genes expressing\n only NNC isoforms"
    data.class[data.class$associated_gene %in% FSM_just_genes, "class"] <- "Genes expressing\n only FSM isoforms"
    
    isoPerGene$FSM_NNC_class = factor(isoPerGene$FSM_NNC_class, levels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"),
                                labels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"), order=T)
    
    p13 <- ggplot(data=isoPerGene[!is.na(isoPerGene$FSM_NNC_class),], aes(x=FSM_NNC_class, y=log2(geneExp+1))) +
      geom_violin(aes(fill=FSM_NNC_class), draw_quantiles = c(0.25, 0.5, 0.75)) +
      ylab("log2( # Short reads per gene + 1)") +
      theme(axis.title.x=element_blank()) +
      #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
      scale_fill_manual(values = c("grey38",myPalette[[4]],myPalette[[1]])) +
      guides(fill=FALSE) +
      mytheme +
      theme(axis.title.x=element_blank()) +
      ggtitle("Gene expression level in NNC/FSM containing genes\n\n" ) +
      scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
                                "Genes expressing\n only FSM isoforms",
                                "Genes expressing\n only NNC isoforms"),
                       labels=c("NNC/FSM genes",
                                "FSM genes",
                                "NNC genes"), drop=FALSE) 
    
    
    p13.c <- ggplot(data=data.class[!is.na(data.class$class),], aes(x=class, y=log2(iso_exp+1))) +
      geom_violin(aes(fill=structural_category), draw_quantiles = c(0.25, 0.5, 0.75)) +
      ylab("log2( # Short reads per transcript + 1)") +
      theme(axis.title.x=element_blank()) +
      #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
      scale_fill_manual(values = myPalette) +
      guides(fill=FALSE) +
      mytheme +
      theme(axis.title.x=element_blank()) +
      ggtitle("Transcript expression level in NNC/FSM containing genes\n\n" ) +
      scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
                                "Genes expressing\n only FSM isoforms",
                                "Genes expressing\n only NNC isoforms"),
                       labels=c("NNC/FSM genes",
                                "FSM genes",
                                "NNC genes"), drop=F) 
    
  }
}



# PLOT 23: Junction categories
# PLOT 24-26: Junction distance to TSS

if (nrow(data.junction) > 0){

    p23.a <- ggplot(data.junction, aes(x=structural_category)) +
      geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=SJ_type), color="black",  size=0.3, width = 0.7) +
      scale_y_continuous(labels = percent, expand = c(0,0)) +
      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
      ylab("% of Splice Junctions") +
      mytheme +
      guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
      theme(legend.position="bottom", legend.title=element_blank())  +
      theme(axis.text.x = element_text(angle = 45)) +
      theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
      theme(axis.title.x=element_blank()) +
      ggtitle("Distribution of Splice Junctions by Structural Classification\n\n\n")

    t <- subset(data.class, exons > 1)  # select only multi-exon isoforms

    p23.b <- ggplot(data=t, aes(x=structural_category)) +
      geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=all_canonical), color="black", size=0.3, width = 0.7) +
      scale_y_continuous(labels = percent, expand = c(0,0)) +
      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
      xlab("") +
      ylab("% of Transcripts ") +
      mytheme +
      guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
      theme(legend.position="bottom", legend.title=element_blank())  +
      theme(axis.text.x = element_text(angle = 45)) +
      theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
      theme(axis.title.x=element_blank()) +
      ggtitle("Distribution of Transcripts by Splice Junctions\n\n\n")


#   p24 <- ggplot(data=data.junction, aes(x=transcript_coord, fill = canonical_known)) +
#     geom_density(alpha=0.7,  size=0.3) +
#     scale_y_continuous(expand = c(0,0))  +
#     scale_x_continuous(expand = c(0,0)) +
#     scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
#     #geom_vline(xintercept = 200, linetype = "longdash", color="red") +
#     mytheme +
#     guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7))+
#     theme(legend.position="bottom" ) +
#     labs(list(fill = "Junction Type", x= "Distance to TSS (bp)",
#               title = "Distribution of splice junctions distance to TSS\n\n") )
#
#   p25 <- ggplot(data=data.junction, aes(x=TSSrange, fill=canonical_known)) +
#     geom_bar(aes(y = (..count..)/sum(..count..)), color="black", size=0.3, width=0.7, position="fill") +
#     scale_y_continuous(labels = percent, expand = c(0,0)) +
#     scale_x_discrete(drop=FALSE) +
#     scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=F) +
#     mytheme +
#     theme(legend.position="bottom") +
#     guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.7))+
#     labs(list(fill = "Junction Type", x= "TSS distance range (bp)", y="% Junctions",
#               title = "Splice junction distance to TSS across junction type\n\n\n") ) +
#     theme(axis.title.x  = element_text(margin=margin(10,0,0,0), size=12))
#
#   p26 <- ggplot(data=data.junction, aes(x=canonical_known, fill=TSSrange)) +
#     geom_bar(aes(y = (..count..)/sum(..count..)), color="black", size=0.3, width=0.7, position="fill") +
#     scale_y_continuous(labels = percent, expand = c(0,0)) +
#     scale_fill_manual(values = myPalette, drop=F) +
#     scale_x_discrete(drop=FALSE) +
#     ylab("% Junctions") +
#     mytheme +
#     theme(legend.position="bottom") +
#     guides(fill = guide_legend(title = "TSS distance range (bp)")) +
#     ggtitle( "Splice junction distance to TSS across junction type\n\n\n") +
#     theme(axis.title.x=element_blank())
#
}


# PLOT 29: RT-switching


if (sum(data.junction$RTS_junction=='TRUE') > 0) {

    a <- data.frame(table(data.junction$SJ_type));
    b <- data.frame(table(subset(data.junction, RTS_junction=="TRUE")$SJ_type));

    df.RTS <- merge(a, b, by="Var1");
    df.RTS$perc <- df.RTS$Freq.y/df.RTS$Freq.x *100
    df.RTS[is.na(df.RTS$perc),"perc"] <- 0

    maxH <- min(100, (max(df.RTS$perc) %/% 5) * 5 + 5);

    p29.a <- ggplot(data=df.RTS, aes(x=Var1, y=perc, fill=Var1)) +
       geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
       geom_text(label=paste(round(df.RTS$perc),"%",sep=''), nudge_y=0.3) +
       scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=F) +
       labs(x="", y="% RT-switching junctions") +
       ggtitle("RT-switching, all junctions\n\n" ) +
       mytheme +
       guides(fill=FALSE) +
       scale_y_continuous(expand = c(0,0), limits = c(0,maxH)) +
       theme(axis.text.x = element_text(size=11))

    c <- data.frame(table(uniqJuncRTS$SJ_type));
    d <- data.frame(table(subset(uniqJuncRTS, RTS_junction=='TRUE')$SJ_type));


    df.uniqRTS <- merge(c, d, by="Var1");
    df.uniqRTS$perc <- df.uniqRTS$Freq.y/df.uniqRTS$Freq.x *100
    df.uniqRTS[is.na(df.uniqRTS$perc),"perc"] <- 0


    p29.b <- ggplot(data=df.uniqRTS, aes(x=Var1, y=perc, fill=Var1)) +
       geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
       geom_text(label=paste(round(df.uniqRTS$perc),"%",sep=''), nudge_y=0.3) +
       scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=F) +
       labs(x="", y="% RT-switching junctions") +
       ggtitle("RT-switching, unique junctions\n\n" ) +
       mytheme +
       guides(fill=FALSE) +
       scale_y_continuous(expand = c(0,0), limits = c(0,maxH)) +
       theme(axis.text.x = element_text(size=11))

}


# PLOT pn4-5: Splice Junction Coverage (if coverage provided)

if (!all(is.na(data.junction$total_coverage))){

    uniqJuncCov <- unique(data.junction[,c("junctionLabel","SJ_type", "total_coverage")])

    e <- data.frame(table(uniqJuncCov$SJ_type))
    f <- data.frame(table(uniqJuncCov[which(uniqJuncCov$total_coverage>0),"SJ_type"]))


    df.juncSupport <- data.frame(type=e$Var1, count=e$Freq-f$Freq, name='Unsupported')
    df.juncSupport <- rbind(df.juncSupport, data.frame(type=f$Var1, count=f$Freq, name='Supported'))

    pn4.a <- ggplot(df.juncSupport, aes(x=type, y=count, fill=name)) +
             geom_bar(stat='identity') +
             scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
             scale_y_continuous( expand = c(0,0)) +
             labs(x='', y='# of Junctions', title='Unique junctions w/ or w/out short read coverage\n\n\n') +
             mytheme +
             theme(legend.position="bottom", legend.title=element_blank()) +
             guides(fill = guide_legend(title = "") )


    df.SJcov <- merge(e, f, by="Var1")
    # calculate the percentage of junctions that have zero short read junction coverage
    df.SJcov$perc <- 100-df.SJcov$Freq.y / df.SJcov$Freq.x * 100;
    df.SJcov[is.na(df.SJcov$perc), "perc"] <- 0

    pn4.b <- ggplot(df.SJcov, aes(x=Var1,fill=Var1, y=perc)) +
      geom_bar(stat="identity", position = position_dodge(), color="black", size=0.3, width=0.7) +
        geom_text(label=paste(round(df.SJcov$perc,1),"%",sep=''), nudge_y=0.3) +
      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
      scale_y_continuous( expand = c(0,0)) +
        labs(x='', y='# of Junctions', title='Unique junctions w/out short read coverage (percentage)\n\n\n') +
      ylab("% of Junctions") +
      mytheme +
      guides(fill=FALSE)


}


# # PLOT pn1.2: Splice Junction relative coverage (if coverage and expression provided)
#
# if (nrow(data.junction) > 0){
#
#   if (!all(is.na(data.junction$total_coverage)) & !all(is.na(data.class$iso_exp))){
#
#     data.junction$isoExp = data.class[data.junction$isoform, "iso_exp"]
#
#     total = aggregate(cbind(total_coverage,isoExp,transcript_coord) ~ junctionLabel, data = data.junction,
#                       FUN = function(x) c(mn = sum(x), n = min(x) ) )
#
#     total$relCov = total$total_coverage[,"n"] / total$isoExp[,"mn"]
#     total$minTSS = total$transcript_coord[,"n"]
#
#     uniqJunc = unique(data.junction[,c("junctionLabel", "canonical_known", "total_coverage")])
#     uniqJunc$notCov = uniqJunc$total_coverage == 0
#
#     uniqueJunc_nonCov = as.data.frame(table(uniqJunc[uniqJunc$totalCoverage==0,"canonical_known"])/table(uniqJunc$canonical_known)*100)
#
#     uniqJunc2 = merge(total, uniqJunc, by=1)
#     uniqJunc2$TSSrange =cut(uniqJunc2$minTSS, breaks = c(0,40,80,120,160,200,10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))
#
#
#
#     # calculate total expression associated to each unique junction
#     sumExpPerJunc = tapply(data.junction$isoExp, data.junction$junctionLabel, sum)
#
#     data.junction$sumIsoExp = sumExpPerJunc[data.junction$junctionLabel]
#
#     data.junction$relCov = data.junction$total_coverage / data.junction$sumIsoExp
#
#     max_dist = max(data.junction$transcript_coord) +1
#
#     data.junction$TSSrange = cut(data.junction$transcript_coord, breaks = c(0,20,40,60,80,100,120,140,160,180,200,max_dist), labels = c("0-20", "21-40","41-80","61-80", "81-100","101-120", "121-140","141-160", "161-180", "181-200", ">200"))
#
#     pn1.2 <-ggplot(data=data.junction[data.junction$relCov<1,], aes(y=relCov,x=TSSrange,fill=canonical_known)) +
#       geom_boxplot(outlier.size = 0.2, size=0.3) +
#       scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
#       ylab("Relative coverage") +
#       xlab("# TSS distance range") +
#       mytheme_bw +
#       theme(legend.position="bottom", legend.title=element_blank())  +
#       ggtitle( "Relative Coverage of junctions\n\n\n") +
#       theme(axis.text.x = element_text(angle = 45,margin=margin(15,0,0,0), size=12))
#
#
#   }else{    uniqJunc = unique(data.junction[,c("junctionLabel", "canonical_known")])
# }
#
# }
#
#
#
#
# PLOT p21-22: Full-lengthness (if FSM/ISM transcripts)

if (nrow(data.FSM) > 0) {

    diff_max <- max(max(abs(data.FSM$diff_to_TSS)), max(abs(data.FSM$diff_to_TTS)));
    diff_breaks <- c(-(diff_max+1), seq(-200, 200, by = 20), diff_max+1);
    data.FSM$diffTTSCat = cut(-(data.FSM$diff_to_TTS), breaks = diff_breaks);
    data.FSM$diffTSSCat = cut(-(data.FSM$diff_to_TSS), breaks = diff_breaks);

    max_height <- max(max(table(data.FSM$diffTSSCat)), max(table(data.FSM$diffTTSCat)));
    max_height <- (max_height %/% 10+1) * 10;

    # plot histogram of distance to polyA site, Y-axis absolute count
    p21.a <- ggplot(data=data.FSM, aes(x=diffTTSCat)) +
      geom_bar(fill=myPalette[4], color="black", size=0.3) +
      scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
      mytheme +
      scale_x_discrete(drop=F) +
      ylab("Number of FSM Transcripts")+
      xlab("Distance to Annotated Polyadenylation Site (bp)")+
      labs(     title="Distance to Annotated Polyadenylation Site, FSM only\n\n",
               subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # plot histogram of distance to polyA site, Y-axis percentages
    p21.b <- ggplot(data=data.FSM, aes(x=diffTTSCat)) +
      geom_bar(aes(y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3)+
      scale_y_continuous(labels = percent_format(), limits = c(0,1), expand = c(0,0))+
      scale_x_discrete(drop=F) +
      mytheme +
      ylab("Percent of FSM Transcripts")+
      xlab("Distance to Annotated Polyadenylation Site (bp)")+
      labs(     title="Distance to Annotated Polyadenylation Site, FSM only\n\n",
               subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # plot histogram of distance to start site, Y-axis absolute count
    p22.a <- ggplot(data=data.FSM, aes(x=diffTSSCat)) +
      geom_bar(fill=myPalette[6], color="black", size=0.3)+
      scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
      scale_x_discrete(drop=F) +
      mytheme +
      ylab("Number of FSM Transcripts")+
      xlab("Distance to Annotated Transcription Start Site (bp)")+
      labs(     title="Distance to Annotated Transcription Start Site, FSM only\n\n",
               subtitle="Negative values indicate downstream of annotated TSS\n\n") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # plot histogram of distance to start site, Y-axis absolute count
    p22.b <- ggplot(data=data.FSM, aes(x=diffTSSCat)) +
      geom_bar(aes(y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3)+
      scale_y_continuous(labels = percent_format(), limits = c(0,1), expand = c(0,0))+
      scale_x_discrete(drop=F) +
      mytheme +
      ylab("Percent of FSM Transcripts")+
      xlab("Distance to Annotated Transcription Start Site (bp)")+
      labs(title="Distance to Annotated Transcription Start Site, FSM only\n\n",
           subtitle="Negative values indicate downstream of annotated TSS\n\n") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


if (nrow(data.ISM) > 0) {

    diff_max <- max(max(abs(data.ISM$diff_to_TSS)), max(abs(data.ISM$diff_to_TTS)));
    diff_breaks <- c(-(diff_max+1), seq(-10000, 10000, by = 1000), diff_max+1);
    data.ISM$diffTTSCat = cut(-(data.ISM$diff_to_TTS), breaks = diff_breaks);
    data.ISM$diffTSSCat = cut(-(data.ISM$diff_to_TSS), breaks = diff_breaks);

    max_height <- max(max(table(data.ISM$diffTSSCat)), max(table(data.ISM$diffTTSCat)));
    max_height <- (max_height %/% 10+1) * 10;

    p21.dist3.ISM.a <- ggplot(data=data.ISM, aes(x=diffTTSCat)) +
      geom_bar(fill=myPalette[4], color="black", size=0.3) +
      scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
      mytheme +
      scale_x_discrete(drop=F) +
      ylab("Number of ISM Transcripts")+
      xlab("Distance to Annotated Polyadenylation Site (bp)")+
      labs(     title="Distance to Annotated Polyadenylation Site, ISM only\n\n",
               subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # plot histogram of distance to polyA site, Y-axis percentages
    p21.dist3.ISM.b <- ggplot(data=data.ISM, aes(x=diffTTSCat)) +
      geom_bar(aes(y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3)+
      scale_y_continuous(labels = percent_format(), limits = c(0,1), expand = c(0,0))+
      scale_x_discrete(drop=F) +
      mytheme +
      ylab("Percent of ISM Transcripts")+
      xlab("Distance to Annotated Polyadenylation Site (bp)")+
      labs(     title="Distance to Annotated Polyadenylation Site, ISM only\n\n",
               subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))


    # plot histogram of distance to start site, Y-axis absolute count
    p22.dist5.ISM.a <- ggplot(data=data.ISM, aes(x=diffTSSCat)) +
      geom_bar(fill=myPalette[6], color="black", size=0.3)+
      scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
      scale_x_discrete(drop=F) +
      mytheme +
      ylab("Number of FSM Transcripts")+
      xlab("Distance to Annotated Transcription Start Site (bp)")+
      labs(     title="Distance to Annotated Transcription Start Site, ISM only\n\n",
               subtitle="Negative values indicate downstream of annotated TSS\n\n") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # plot histogram of distance to start site, Y-axis absolute count
    p22.dist5.ISM.b <- ggplot(data=data.ISM, aes(x=diffTSSCat)) +
      geom_bar(aes(y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3)+
      scale_y_continuous(labels = percent_format(), limits = c(0,1), expand = c(0,0))+
      scale_x_discrete(drop=F) +
      mytheme +
      ylab("Percent of FSM Transcripts")+
      xlab("Distance to Annotated Transcription Start Site (bp) ")+
      labs(title="Distance to Annotated Transcription Start Site, ISM only\n\n",
           subtitle="Negative values indicate downstream of annotated TSS\n\n") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


if (sum(!is.na(data.class$polyA_dist)) > 10) {
p.polyA_dist <- ggplot(data.class, aes(x=polyA_dist, color=structural_category)) +
    geom_freqpoly(binwidth=1) +
    xlab("Distance of polyA motif from 3' end (bp)") +
    ylab("Count") +
    labs(title="Distance of detected polyA motif from 3' end")
}

# PLOT p28: Attribute summary if junctions

if (nrow(data.junction) > 0){

        # ToDo: USE COVERAGE DATA LATER

        # for FSM, ISM, NIC, and NNC, plot the percentage of RTS and non-canonical junction
        x <- filter(data.class, structural_category %in% c("FSM", "ISM", "NIC", "NNC" ) & exons > 1)

        t1.RTS <- group_by(x, structural_category, RTS_stage) %>% summarise(count=n())
        t2.RTS <- group_by(x, structural_category) %>% summarise(count=n())
        t3.RTS <- merge(t1.RTS, t2.RTS, by="structural_category")
        t3.RTS$perc <- t3.RTS$count.x / t3.RTS$count.y * 100
        t3.RTS <- subset(t3.RTS, RTS_stage=='TRUE');

        t1.SJ <- group_by(x, structural_category, all_canonical) %>% summarise(count=n())
        t3.SJ <- merge(t1.SJ, t2.RTS, by="structural_category")
        t3.SJ$perc <- t3.SJ$count.x / t3.SJ$count.y * 100
        t3.SJ <- subset(t3.SJ, all_canonical=='non_canonical');

        p28.RTS <- ggplot(t3.RTS, aes(x=structural_category, y=perc)) +
            geom_col(position='dodge', width = 0.7,  size=0.3, fill='darkred', color="black") +
            geom_text(label=paste(round(t3.RTS$perc, 1),"%",sep=''), nudge_y=0.5) +
            scale_fill_manual(values = myPalette[9:11]) +
            ylab("% of Isoforms") +
            xlab("") +
            mytheme +
            theme(legend.position="bottom", axis.title.x = element_blank()) +
            ggtitle("Incidence of RT-switching\n\n") +
            guides(fill = guide_legend(title = "QC Attributes") )

        p28.SJ <- ggplot(t3.SJ, aes(x=structural_category, y=perc)) +
            geom_col(position='dodge', width = 0.7,  size=0.3, fill='lightblue', color="black") +
            geom_text(label=paste(round(t3.SJ$perc, 1),"%",sep=''), nudge_y=0.5) +
            scale_fill_manual(values = myPalette[9:11]) +
            ylab("% of Isoforms") +
            xlab("") +
            mytheme +
            theme(legend.position="bottom", axis.title.x = element_blank()) +
            ggtitle("Incidence of Non-Canonical Junctions\n\n") +
            guides(fill = guide_legend(title = "QC Attributes") )

}



# PLOT p30,p31,p32: percA by subcategory
p30 <- ggplot(data=data.class, aes(y=perc_A_downstream_TTS, x=structural_category, fill=subcategory)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
    mytheme +
    xlab("Structural Category") +
    ylab("Percent 'A's (%) ") +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(legend.position="bottom", legend.title=element_blank(), legend.direction = "horizontal", legend.box = "vertical") +
    guides(fill=guide_legend(nrow=5,byrow=TRUE)) +
    theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12)) +
    labs(title="Possible Intra-Priming by Structural Category\n\n",
         subtitle="Percent of genomic 'A's in downstream 20 bp\n\n") +
    theme(axis.title.x=element_blank()) +
    scale_fill_manual(values=myPalette, breaks=c("intron_retention", "3prime_fragment", "internal_fragment", "5prime_fragment",
                             "mono-exon", "multi-exon", "combination_of_known_junctions",
                             "no_combination_of_known_junctions", "mono-exon_by_intron_retention/s",
                            "not any annotated donor/acceptor", "any annotated donor/acceptor"),
                      labels=c("Intron retention", "3' fragment", "Internal fragment", "5' fragment",
                             "Mono-exon", "Multi-exon", "Combination of annotated junctions",
                             "Not combination of annotated junctions", "Mono-exon by intron retention",
                             "Without annotated donors/acceptors", "At least one annotated donor/acceptor"), drop=F)

p31 <- ggplot(data=data.class, aes(y=perc_A_downstream_TTS, x=structural_category, fill=exonCat)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
    mytheme +
    scale_fill_manual(breaks=c("Mono-Exon", "Multi-Exon"),
                      labels=c("Mono-Exon Isoforms", "Multi-Exon Isoforms"), values=myPalette) +
    ylab("Percent 'A's (%) ") +
    theme(legend.position="bottom", legend.title=element_blank()) +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
    labs(title = "Possible Intra-Priming, Mono- vs Multi-Exon\n\n",
         subtitle = "Percent of genomic 'A's in downstream 20 bp\n\n") +
    theme(axis.title.x=element_blank())

p32 <- ggplot(data=data.class, aes(y=perc_A_downstream_TTS, x=structural_category, fill=coding)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
    scale_fill_manual(breaks=c("coding", "non_coding"),
                      labels=c("Coding Isoforms", "Non-Coding Isoforms"), values=myPalette[3:4]) +
    ylab("Percent 'A's (%) ") +
    theme(legend.position="bottom", legend.title=element_blank() ) +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x  = element_text(margin=margin(17,0,0,0), size=12))+
    labs(title = "Possible Intra-Priming, Coding vs Non-Coding\n\n",
         subtitle = "Percent of genomic 'A's in downstream 20 bp\n\n") +
    theme(axis.title.x=element_blank())









###** Output plots

pdf(file=report.file, width = 6.5, height = 6.5)


# cover
grid.newpage()
cover <- textGrob("SQANTI2 report",
    gp=gpar(fontface="italic", fontsize=40, col="orangered"))
grid.draw(cover)


# document the parameters, if available
grid.newpage()
if (length(args)>2) {
  param <- read.table(param.file, sep='\t', header=F)
  table.param <- tableGrob(param, rows = NULL, cols = NULL)
  grid.draw(table.param)
}


# TABLE 1: Number of isoforms in each structural category
freqCat <- group_by(data.class, by=structural_category) %>% summarise(num_iso=n(), num_gene=length(unique(associated_gene)))
table1 <- tableGrob(freqCat, rows = NULL, cols = c("Category","# Isoforms", "# Genes"))
title1 <- textGrob("Characterization of transcripts\n based on splice junctions", gp=gpar(fontface="italic", fontsize=17), vjust = -3.5)
gt1 <- gTree(children=gList(table1, title1))

# TABLE 2: Number of Novel vs Known Genes
freqCat = as.data.frame(table(isoPerGene$novelGene))
table2 <- tableGrob(freqCat, rows = NULL, cols = c("Category","# Genes"))
title2 <- textGrob("Gene classification", gp=gpar(fontface="italic", fontsize=17), vjust = -3.5)
gt2 <- gTree(children=gList(table2, title2))


# TABLE 3: Junction Classification

uniq_sj_count <- nrow(uniqJunc)

freqCat <- as.data.frame(table(uniqJunc$SJ_type))
freqCat$Var1 <- gsub(" ", "", freqCat$Var1);
freqCat$Var1 <- gsub("\n", " ", freqCat$Var1);
freqCat$Frac <- round(freqCat$Freq*100 / uniq_sj_count, 2)
table2 <- tableGrob(freqCat, rows = NULL, cols = c("Category","# SJs","Percent"))
title2 <- textGrob("Splice Junction Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
gt3 <- gTree(children=gList(table2, title2))


# TABLE 4: Summary number of Unique Isoforms and Unique Genes
nGenes = nrow(isoPerGene)
nIso = nrow(data.class)
sn = paste("Unique Genes: ", nGenes, "\n", "Unique Isoforms: ", nIso)
gt4 <- textGrob(sn, gp=gpar(fontface="italic", fontsize=17), vjust = 0)


# Plot Table 1 and Table 2
grid.arrange(gt4,gt2,gt3,gt1, layout_matrix = cbind(c(1,2,3),c(1,4,4)))


s <- textGrob("Gene Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
print(p0)
print(p7)
print(p6)
print(p.classByLen.a)
print(p.classByLen.b)

if (!all(is.na(data.class$iso_exp))){
   print(p10)
}
if (!all(is.na(data.class$FL))){
   print(p11)
}

# PLOT length of isoforms
# p.length.all: length of all isoforms, regardless of category
# p.length.cat: length of isoforms, by category
# p.length.exon: length of isoforms, mono- vs mult-exon
# (optional) p.length.all.sample: length of all isoforms by sample
print(p.length.all)
print(p.length.cat)
print(p.length.exon)
if (length(FL_multisample_indices)>0) {
    print(p.length.all.sample)
    print(p.length.exon.sample)
}

# 2. general parameters by structual categories
s <- textGrob("Structural Isoform Characterization\nby Splice Junctions", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
print(p1)
print(p4)
print(p5)
if (!all(is.na(data.class$iso_exp))){
    print(p8)
}
if (!all(is.na(data.class$FL))){
    print(p9)
}

# (optional) table of FL counts by structural category
# case 1: single FL sample
if (!all(is.na(data.class$FL)))
{
    m1 <- group_by(data.class, by=structural_category) %>% summarise(count=n(), fl=sum(FL))
    cols <- c("category", "isoforms", "FL")
    table.FL <- tableGrob(m1, rows = NULL, cols = cols)
    title.FL <- textGrob("FL counts by category", gp=gpar(fontface="italic", fontsize=17), vjust = -10)
    gt.FL <- gTree(children=gList(table.FL, title.FL))
    grid.arrange(gt.FL, ncol=1)
}

# case 2: multi FL sample
if (length(FL_multisample_indices)>0)
{
    j <- FL_multisample_indices[1]
    name <- colnames(data.class)[j]
    cols <- c("category", "isoforms", name)

    m1 <- group_by(data.class, by=structural_category) %>% summarise(count=n(), fl=sum(!!sym(name)))

    for (i in 2:length(FL_multisample_indices)) {
        j <- FL_multisample_indices[i]
        name <- colnames(data.class)[j]

        m2 <- group_by(data.class, by=structural_category) %>% summarise(fl=sum(!!sym(name)))
        m1 <- merge(m1, m2, by="by")
        cols <- c(cols, name)
    }
    colnames(m1) <- cols

    table.FL <- tableGrob(m1, rows = NULL, cols = cols)
    title.FL <- textGrob("FL counts by category", gp=gpar(fontface="italic", fontsize=17), vjust = -10)
    gt.FL <- gTree(children=gList(table.FL, title.FL))
    grid.arrange(gt.FL, ncol=1)
}


# (optional) p.FL_TPM_sample.by_cat
# (optional) p.FL_TMP_sample.by_length
if (length(FL_multisample_indices)>0) {
    data.class$length_cat <- "<1kb"
    data.class[data.class$length>=1000,"length_cat"] <- "1-3kb"
    data.class[data.class$length>=3000&data.class$length<5000,"length_cat"] <- "3-5kb"
    data.class[data.class$length>=5000&data.class$length<10000,"length_cat"] <- "5-10kb"
    data.class[data.class$length>10000,"length_cat"] <- ">10kb"

    for (i in 1:(length(FL_TPM_multisample_names)-1)) {
        j1 <- FL_TPM_multisample_names[i]
        j1_log10 <- paste(FL_TPM_multisample_names[i], "_log10", sep='')
        n1 <- FL_multisample_names[i]
        for (i2 in (i+1):length(FL_multisample_names)) {
            j2 <- FL_TPM_multisample_names[i2]
            j2_log10 <- paste(FL_TPM_multisample_names[i2], "_log10", sep='')
            n2 <- FL_multisample_names[i2]

            print(paste("Printing FL TPM for sample", j1, "vs", j2, "...."))

            max_j1j2 <- floor(max(data.class[,j1_log10], data.class[,j2_log10])) + 1
            pearson <- round(cor(data.class[,j1], data.class[,j2], method="pearson"), 2)
            p.tmp <- ggplot(data.class, aes_string(j1_log10, j2_log10, color="structural_category")) +
                geom_point(alpha=0.3) +
                annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
                xlim(c(0, max_j1j2)) +
                ylim(c(0, max_j1j2)) +
                labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2)) +
                guides(fill=FALSE) +
                mytheme
            print(p.tmp)
            ggsave(paste("Rplot.",n1,"vs",n2,".by_cat.png",sep=''), p.tmp, width=8, height=6)

            data.class.gene <- group_by(data.class, by=associated_gene) %>% summarise(n=n(), sum1=sum(!!sym(j1)), sum2=sum(!!sym(j2)))
            pearson <- round(cor(data.class.gene$sum1, data.class.gene$sum2, method="pearson"), 2)
            p.tmp.gene <- ggplot(data.class.gene, aes(x=log10(sum1), y=log10(sum2))) +
                geom_point(alpha=0.3, color='orange') +
                annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
                xlim(c(0, max_j1j2)) +
                ylim(c(0, max_j1j2)) +
                xlab(j1) +
                ylab(j2) +
                labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, ", grouped by gene")) +
                guides(fill=FALSE) +
                mytheme
            print(p.tmp.gene)
            ggsave(paste("Rplot.",n1,"vs",n2,".summed_by_gene.png",sep=''), p.tmp.gene, width=8, height=6)


            data.class.tmp <- subset(data.class,length_cat=="<1kb")
            pearson <- round(cor(data.class.tmp[,j1], data.class.tmp[,j2], method="pearson"), 2)
            p.tmp.le1k <- ggplot(data.class.tmp, aes_string(j1_log10, j2_log10)) +
                geom_point(alpha=0.3, color='orange') +
                annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
                xlim(c(0, max_j1j2)) +
                ylim(c(0, max_j1j2)) +
                labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, "< 1kb only"))+
                guides(fill=FALSE) +
                mytheme
            print(p.tmp.le1k)
            ggsave(paste("Rplot.",n1,"vs",n2,".le1k.png",sep=''), p.tmp.le1k, width=8, height=6)

            data.class.tmp <- subset(data.class,length_cat=="1-3kb")
            pearson <- round(cor(data.class.tmp[,j1], data.class.tmp[,j2], method="pearson"), 2)
            p.tmp.1to3k <- ggplot(data.class.tmp, aes_string(j1_log10, j2_log10)) +
                geom_point(alpha=0.3, color='purple') +
                annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
                xlim(c(0, max_j1j2)) +
                ylim(c(0, max_j1j2)) +
                labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, "1-3kb only")) +
                guides(fill=FALSE) +
                mytheme
            print(p.tmp.1to3k)
            ggsave(paste("Rplot.",n1,"vs",n2,".1to3k.png",sep=''), p.tmp.1to3k, width=8, height=6)

            data.class.tmp <- subset(data.class,length_cat=="3-5kb")
            pearson <- round(cor(data.class.tmp[,j1], data.class.tmp[,j2], method="pearson"), 2)
            p.tmp.3to5k <- ggplot(data.class.tmp, aes_string(j1_log10, j2_log10)) +
                geom_point(alpha=0.3, color='royalblue4') +
                annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
                xlim(c(0, max_j1j2)) +
                ylim(c(0, max_j1j2)) +
                labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, "3-5kb only")) +
                guides(fill=FALSE) +
                mytheme
            print(p.tmp.3to5k)
            ggsave(paste("Rplot.",n1,"vs",n2,".3to5k.png",sep=''), p.tmp.3to5k, width=8, height=6)

            data.class.tmp <- subset(data.class,length_cat=="5-10kb")
            pearson <- round(cor(data.class.tmp[,j1], data.class.tmp[,j2], method="pearson"), 2)
            p.tmp.5to10k <- ggplot(data.class.tmp, aes_string(j1_log10, j2_log10)) +
                geom_point(alpha=0.3, color='hotpink4') +
                annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
                xlim(c(0, max_j1j2)) +
                ylim(c(0, max_j1j2)) +
                labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, "5-10 kb only")) +
                guides(fill=FALSE) +
                mytheme
            print(p.tmp.5to10k)
            ggsave(paste("Rplot.",n1,"vs",n2,".5to10k.png",sep=''), p.tmp.5to10k, width=8, height=6)

            data.class.tmp <- subset(data.class,length_cat==">10kb")
            pearson <- round(cor(data.class.tmp[,j1], data.class.tmp[,j2], method="pearson"), 2)
            p.tmp.ge10k <- ggplot(data.class.tmp, aes_string(j1_log10, j2_log10)) +
                geom_point(alpha=0.3, color='darkolivegreen4') +
                annotate("text", x=max_j1j2-0.5, y=max_j1j2-0.5, label=paste("Pearson:", pearson)) +
                xlim(c(0, max_j1j2)) +
                ylim(c(0, max_j1j2)) +
                labs(title=paste("FL TPM (log10 scale)", n1, "vs", n2, ">10 kb only")) +
                guides(fill=FALSE) +
                mytheme
            print(p.tmp.ge10k)
            ggsave(paste("Rplot.",n1,"vs",n2,".ge10k.png",sep=''), p.tmp.ge10k, width=8, height=6)


            #grid.arrange(p.tmp.le1k, p.tmp.1to3k, p.tmp.3to5k, p.tmp.5to10k, p.tmp.ge10k, ncol=2)

        }
    }
}

#
if (nrow(data.FSM) > 0 ) {
    print(p2)
    print(p3)
}
if (!all(is.na(data.class$gene_exp))){
   if (nrow(data.class[data.class$structural_category=="NNC",])!=0){
     print(p12)
   }
}
#if (!all(is.na(data.class$gene_exp))){
#    if (nrow(data.class[data.class$structural_category=="NNC",])!=0 & nrow(data.class[data.class$structural_category=="FSM",])!=0 ){]
#        print(p13)
#        print(p13.c)
#    }
#}


#3. splice junction

s <- textGrob("Splice Junction Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
print(p23.a)
print(p23.b)
#   print(p24)
#   print(p25)
#   print(p26)
#
#   if (!all(is.na(data.junction$total_coverage)) & !all(is.na(data.class$iso_exp))){
#     print(pn1.2)
#   }
#

if (!all(is.na(data.junction$total_coverage))) {
    print(pn4.a)
    print(pn4.b)
}

if (sum(data.junction$RTS_junction=='TRUE') > 0) {
    print(p29.a)
    print(p29.b)
}


s <- textGrob("Comparison with Annotated TSS and PolyA Sites", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
if (nrow(data.FSM) > 0) {
    print(p21.a)
    print(p21.b)
    print(p22.a)
    print(p22.b)
}

if (nrow(data.ISM) > 0) {
    print(p21.dist3.ISM.a)
    print(p21.dist3.ISM.b)
    print(p22.dist5.ISM.a)
    print(p22.dist5.ISM.b)
}

if (sum(!is.na(data.class$polyA_dist)) > 10) {
    print(p.polyA_dist)

    # PLOT polyA motif ranking, distance from 3' end
    df.polyA <- as.data.frame(group_by(data.class, by=structural_category) %>%
                summarise(count=n(),
                          polyA_detected=sum(!is.na(polyA_motif)),
                          polyA_detected_perc=round(polyA_detected*100/count)))

    table.polyA <- tableGrob(df.polyA, rows = NULL, cols = c("Category","Count","polyA\nDetected","%"))
    title.polyA <- textGrob("Number of polyA Motifs Detected", gp=gpar(fontface="italic", fontsize=15), vjust=-10)
    gt.polyA <- gTree(children=gList(table.polyA, title.polyA))


    df.polyA_freq <- as.data.frame(sort(table(data.class$polyA_motif),decreasing=T))
    df.polyA_freq$perc <- round(df.polyA_freq$Freq*100/sum(df.polyA_freq$Freq),1)
    table.polyA_freq <- tableGrob(df.polyA_freq, rows = NULL, cols = c("Motif", "Count", "%"))
    title.polyA_freq <- textGrob("Frequency of polyA motifs", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
    gt.polyA_freq <- gTree(children=gList(title.polyA_freq, table.polyA_freq))

    grid.arrange(gt.polyA, gt.polyA_freq, ncol=2)
}

if (sum(!is.na(data.class$dist_to_cage_peak)) > 10) {

    df.cage <- group_by(data.class, by=structural_category) %>%
        summarise(count=n(), cage=sum(abs(dist_to_cage_peak)<=50,na.rm=T), freq=round(cage*100/count))
    table.cage <- tableGrob(df.cage, rows=NULL, cols=c("Category", "Count", "Has CAGE peak\nwithin 50bp", "%"))
    title.cage <- textGrob("Number of close by CAGE Peaks Detected" ,gp=gpar(fontface="italic", fontsize=15), vjust=-10)
    gt.cage <- gTree(children=gList(table.cage, title.cage))
    grid.arrange(gt.cage, ncol=1)
}

s <- textGrob("Intra-Priming Quality Check", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
print(p30)
print(p31)
print(p32)


s <- textGrob("Quality Controls", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
print(p28.RTS)
print(p28.SJ)

dev.off()


print("SQANTI2 report successfully generated!")

















