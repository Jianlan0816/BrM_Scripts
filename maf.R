library(maftools)
library(stringr)
library(readxl)
library(dplyr)
refgenes = read.table("/Users/Jianlan/Downloads/gene_RefSeqGene.txt",sep = "\t")
censusgenes = read.table("~/Downloads/Census_allSat Dec 11 10_36_25 2021.csv",sep = "," , header = T)
genelist=c("TTN","MUC16","FLG","RYR2","TP53",
           "USH2A","OBSCN","DST","LRP1",
           "PTEN","PIK3R1","STK11","MTOR",
           "PIK3CA","AKT1","TSC1","TSC2","EGFR",
           "ERBB2","MET","ALK","RET",
           "KRAS","NRAS","NF1",
           "BRAF"
           )

genelist2 = c("TTN","TP53","MUC16","RYR2","LRP1B",
              "USH2A","RYR3","MUC17","PCDH15","OBSCN",
              "KMT2D","FAT3","NEB","PCLO","RELN",
              "SYNE1","DST","KRAS","CSMD2","CUBN",
              "DNAH11","FLG","HERC2","RYR1","COL11A1",
              "DNAH9","FBN2","HCN1","RP1L1","SPEG",
              "SYNE2","ANK2","APOB","CACNA1A","FAT1",
              "FSIP2","HUWE1","PTPRD","TCHH","VPS13B",
              "ABCA13","KMT2C","NF1","NRXN1","SPTA1",
              "TG","UNC79","VCAN","HMCN1","PLEC")

genes_brain = read.table("~/Desktop/Lab/Brain_meta/Lung/snvpathway/braingene.txt")


vc_cols = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'Amp',
  'Del'
)

vc_cols2 = RColorBrewer::brewer.pal(n = 7, name = 'Paired')
names(vc_cols2) = c("Nonsense_Mutation",
                    "Missense_Mutation",
                    "Frame_Shift_Del",
                    "Splice_Site",
                    "Frame_Shift_Ins",
                    "Multi_Hit",
                    "In_Frame_Del")

all.lesions <- "/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalbrain_results/all_lesions.conf_90.txt"
amp.genes <- "/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalbrain_results/amp_genes.conf_90.txt"
del.genes <- "/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalbrain_results/del_genes.conf_90.txt"
scores.gis <- "/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalbrain_results/scores.gistic"


clinic = read_excel("/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/Lung_Clinical_Total.xlsx")

colnames(clinic)[31] = "Cohorts"
#merge_brain <- read.maf("/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/Brain_new/maf2/mergebrain2.maf",
#                        removeDuplicatedVariants = TRUE,clinicalData = clinic
#                        )

merge_brain <- read.maf("/Users/Jianlan/Desktop/Lab/Brain_meta/Lung_own/brainmaf/brain_total.maf",
                        removeDuplicatedVariants = TRUE,clinicalData = clinic,
                        #gisticAllLesionsFile = all.lesions,
                        #gisticAmpGenesFile = amp.genes,
                        #gisticDelGenesFile = del.genes,
                        #gisticScoresFile = scores.gis
                        )

merge_brain2 = subsetMaf(merge_brain,genes = refgenes$V3)
#merge_brain3 = subsetMaf(merge_brain2,genes = censusgenes$Gene.Symbol)

oncoplot(merge_brain2,genes = genelist)
oncoplot(merge_primary,top=50)

brain.gistic = readGistic(gisticAllLesionsFile = all.lesions,
                          gisticAmpGenesFile = amp.genes,
                          gisticDelGenesFile = del.genes,
                          gisticScoresFile = scores.gis)

gisticOncoPlot(brain.gistic)

brain.gistic = subsetMaf(merge_brain,genes = refgenes$V3)
plotmafSummary(merge_brain)

fabcolors = RColorBrewer::brewer.pal(n = 3,name = 'Set1')
names(fabcolors) = c("Brastianos 2015", "Fukumura 2021", 
                     "SYSUCC")
sexcolors = c("#FF5733","#8AAFDE")
names(sexcolors) = c("F","M")
smokecolors = c("#D4F9EE","#37BA9C")
names(smokecolors) = c("no","yes")
histologycolor = c("#DAF7A6","#FFC300","#FF5733","#C70039","#C5D1CE")
names(histologycolor) = c("Adenocarcinoma","Lung carcinoma",
                          "Small cell carcinoma","Squamous carcinoma",NA)
pathologycolor = c("#F2D7D5","#E6B0AA","#D98880","#C0392B","#C5D1CE")
names(pathologycolor) = c("I","II","III","IV",NA)
chemocolors = c("#D4F9EE","#37BA9C","#C5D1CE")
names(chemocolors) = c("0","1",NA)
targetcolors = c("#D4F9EE","#37BA9C","#C5D1CE")
names(targetcolors) = c("0","1",NA)
immunocolors = c("#D4F9EE","#37BA9C","#C5D1CE")
names(immunocolors) = c("0","1",NA)
radiocolors = c("#D4F9EE","#37BA9C","#C5D1CE")
names(radiocolors) = c("0","1",NA)
fabcolors = list(Cohorts = fabcolors,Sex=sexcolors,Smoking=smokecolors,
                 Histology = histologycolor,
                 Pathology = pathologycolor,
                 Chemotherapy = chemocolors,
                 Targeted_therapy = targetcolors,
                 Immunotherapy = immunocolors,
                 Radiation_therapy = radiocolors)

tmb_brainnum=merge_brain2@variants.per.sample
tmb_brainnum[,2] = tmb_brainnum$Variants/50
colnames(tmb_brainnum)[2]="TMB"
tmb_brainnum=data.frame(tmb_brainnum)
tmb_brainnum = rbind(tmb_brainnum, c("PP3",0))

tmb_brainnum2 = arrange(tmb_brainnum,Tumor_Sample_Barcode,TMB)
SYSUCC <- tmb_brainnum2[1:23,]
Fuku <- tmb_brainnum2[24:37,]
Bras <- tmb_brainnum2[38:75,]

SYSUCC_order <- arrange(SYSUCC,desc(TMB))
Fuku_order <- arrange(Fuku,desc(TMB))
Bras_order <- arrange(Bras,desc(TMB))
ordersample = c(as.character(Bras_order$Tumor_Sample_Barcode),as.character(SYSUCC_order$Tumor_Sample_Barcode),
                as.character(Fuku_order$Tumor_Sample_Barcode))

par(mar=c(5.1,4.1,4.1,2.1))

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/primary_onco.tiff",width=6000,height=3000,res=300)
oncoplot(maf = merge_primary2,
         colors = vc_cols,
         #top = 20,
         genes = genelist,
         #draw_titv = TRUE,
         #topBarData = tmb_brain,
         leftBarData = NULL,
         removeNonMutated = FALSE,
         clinicalFeatures = c("Cohorts"),
         #clinicalFeatures = names(fabcolors),
         #annotationColor = fabcolors,
         titleText = "Primary",
         sortByAnnotation = TRUE,
         logColBar = TRUE,
         includeColBarCN = FALSE,
         #showTumorSampleBarcodes = TRUE,
         #barcodeSrt = 90,
         sampleOrder = ordersample,
         #rightBarData = NULL,
         drawRowBar = FALSE,
         drawColBar = FALSE,
         fontSize = 1.2,
         #SampleNamefontSize = 1,
         legendFontSize = 1.5,
         sepwd_samples = 3,
         sepwd_genes = 3,
         anno_height = 0.5,
         #legend_height = 8,
         annotationFontSize = 1.5,
         gene_mar = 8,
         barcode_mar = 5,
         bgCol = "white"
         #altered = TRUE,
)

dev.off()

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/lungcnvchrom.tiff",width=2000,height=2000,res=300)
gisticOncoPlot(gistic = brain.gistic,top=60,removeNonAltered = FALSE,
               sortByAnnotation = TRUE,#colors = col,
               sepwd_genes = 4,sepwd_samples = 4)
gisticBubblePlot(gistic = brain.gistic)
gisticChromPlot(gistic=lung.gistic,y_lims = c(3,-1))
dev.off()

laml.titv = titv(maf = merge_brain, useSyn = TRUE)
plotTiTv(laml.titv)

pall.lesions <- "/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalprimary_results/all_lesions.conf_90.txt"
pamp.genes <- "/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalprimary_results/amp_genes.conf_90.txt"
pdel.genes <- "/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalprimary_results/del_genes.conf_90.txt"
pscores.gis <- "/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalprimary_results/scores.gistic"

#merge_primary <- read.maf("/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/Lung_new/maf2/lungmerge2.maf",
#                       removeDuplicatedVariants = TRUE,clinicalData = clinic)
primary.gistic = readGistic(gisticAllLesionsFile = pall.lesions,
                          gisticAmpGenesFile = pamp.genes,
                          gisticDelGenesFile = pdel.genes,
                          gisticScoresFile = pscores.gis)

merge_primary <- read.maf("/Users/Jianlan/Desktop/Lab/Brain_meta/Lung_own/lungmaf/lung_total.maf",
                       removeDuplicatedVariants = TRUE,clinicalData = clinic,
                       #gisticAllLesionsFile = pall.lesions,
                       #gisticAmpGenesFile = pamp.genes,
                       #gisticDelGenesFile = pdel.genes,
                       #gisticScoresFile = pscores.gis
                       )

merge_primary2 = subsetMaf(merge_primary,genes = refgenes$V3)

#merge_primary3 = subsetMaf(merge_primary2,genes = censusgenes$Gene.Symbol)
lung.gistic = readGistic(gisticAllLesionsFile = pall.lesions,
                          gisticAmpGenesFile = pamp.genes,
                          gisticDelGenesFile = pdel.genes,
                          gisticScoresFile = pscores.gis)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/lung_tmb.tiff",width=2500,height=1200,res=300)
laml.mutload = tcgaCompare(maf = c(merge_primary2, merge_brain2), cohortName = c('Primary','BrM'), logscale = TRUE, capture_size = 50,
                           col = c("#33B5FF", "red"),
                           bg_col = c("white", "#2C7FB8"),)
dev.off()


plotmafSummary(merge_primary)
oncoplot(maf = merge_primary,genes = genelist2)
lung_freq = data.frame(genes = genes_brain$Hugo_Symbol,
                  lung_freq = c(0,0,0.02,0.02,0.03,0.02,0.02,
                                0,0,0,0.01,0.01,0.02,0.02,0.02))


tmb_lung=merge_primary2@variants.per.sample
tmb_lung[,2] = tmb_lung$Variants/50
colnames(tmb_lung)[2]="TMB"
tmb_lung=data.frame(tmb_lung)


brain_pvalue = genes_brain[,c(1,4)]

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/primonco.tiff",width=3000,height=1600,res=300)
oncoplot(maf = merge_primary2,
         colors = vc_cols,
         top = 20,
         #genes = genes_primary$Hugo_Symbol,
         #draw_titv = TRUE,
         #topBarData = tmb_brain,
         removeNonMutated = FALSE,
         #clinicalFeatures = names(fabcolors),
         clinicalFeatures = c("Cohorts"),
         annotationColor = fabcolors,
         #titleText = "Primary",
         sortByAnnotation = TRUE,
         logColBar = TRUE,
         includeColBarCN = FALSE,
         #showTumorSampleBarcodes = TRUE,
         #barcodeSrt = 90,
         sampleOrder = ordersample,
         #leftBarData = brain_pvalue,
         fontSize = 0.8,
         #SampleNamefontSize = 1,
         legendFontSize = 1,
         #sepwd_samples = 1,
         #sepwd_genes = 1,
         anno_height = .5,
         legend_height = 5,
         annotationFontSize = 1,
         gene_mar = 8,
         #barcode_mar = 5,
         bgCol = "white"
         #altered = TRUE,
)         
dev.off()

mafSurvival(maf = merge_brain2, genes = 'SMAD4', time = 'OS_Craniotomy_to_censored', Status = 'Survival_Status', isTCGA = FALSE)
prog_geneset = survGroup(maf = merge_brain2, top = 20, geneSetSize = 2, 
                         time = "OS_Craniotomy_to_censored", Status = "Survival_Status", verbose = FALSE)

pt.vs.rt <- mafCompare(m1 = merge_primary2, m2 = merge_brain2,
                       m1Name = 'Primary', m2Name = 'BrM', minMut = 3)
print(pt.vs.rt)
diff=pt.vs.rt$results
diff=diff[order(diff$pval),]
diff = diff[which(diff$pval < 0.1),]
diff$lfc = log(diff$Primary/diff$BrM)
genes_brain = diff[which(diff$lfc<0),]
genes_primary = diff[which(diff$lfc>0),]

write.table(genes_brain$Hugo_Symbol,"/Users/Jianlan/Desktop/braingene.txt",quote = FALSE,
            row.names = FALSE,col.names = FALSE)
write.table(genes_primary$Hugo_Symbol,"/Users/Jianlan/Desktop/lunggene.txt",quote = FALSE,
            row.names = FALSE,col.names = FALSE)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/brainpathway.tiff",width=3000,height=2000,res=300)
oncoplot(maf=merge_brain3,
         pathways = pathway_table,
         colors = vc_cols,
         gene_mar = 30,
         fontSize = 0.6,
         SampleNamefontSize = 0.7,
         legendFontSize = 1,
         sepwd_samples = 2,
         sepwd_genes = 2,
         titleText = "")
dev.off()

OncogenicPathways(maf = merge_brain2,pathways = pathway_table)
OncogenicPathways(maf = merge_primary2)

genes_react = c()
pathways_react = c()

for (i in 1:length(names(pathways.hallmark))){
  genes_react = c(genes_react,pathways.hallmark[[i]])
  pathways_react = c(pathways_react,rep(names(pathways.hallmark)[i],length(pathways.hallmark[[i]])))
}

pathway_table = data.frame(genes_react,pathways_react)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/forestplot.tiff",
     width=2000,height=2000,res=300)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.01,geneFontSize = 0.8)
dev.off()

total_genes = c(genes_primary$Hugo_Symbol[1:5],"EGFR")
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/brain_cooncoplot.tiff",
     width=4000,height=1000,res=300)
coOncoplot(m1 = merge_primary2, 
           m2 = merge_brain2, 
           m1Name = 'Primary', 
           m2Name = 'BrM', 
           #genes = total_genes, 
           removeNonMutated = TRUE,
           #colors = vc_cols,
           sepwd_genes1 = 1,
           sepwd_samples1 = 1,
           sepwd_genes2 = 1,
           sepwd_samples2 = 1,
           legendFontSize = 1,
           #clinicalFeatures1 = c("Cohort"),
           #clinicalFeatures2 = c("Cohort"),
           #annotationColor1 = fabcolors,
           #annotationColor2 = fabcolors,
           #sortByAnnotation1 = TRUE,
           #sortByAnnotation2 = TRUE,
           #anno_height = 1,
           keepGeneOrder = TRUE,
           bgCol = "white"
           )
dev.off()

library("mclust")
brain.hetero = inferHeterogeneity(maf = merge_brain2, tsb="DHP18",vafCol = 'tumor_f',minVaf = 0.32,
                                  #segFile = "/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/brain/brain_beforemerge.seg_copy.txt"
                                  )
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/brain_ithdhp18.tiff",
     width=1000,height=800,res=200)
plotClusters(clusters = brain.hetero, 
             #genes = 'CN_altered',
             tsb="DHP18")
dev.off()

primary.hetero = inferHeterogeneity(maf = merge_primary2, tsb="DHP18",vafCol = 'tumor_f',minVaf = 0.32,
                                    #segFile = "/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/primary/primary_beforemerge.seg_copy.txt"
                                    )
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/lung_ithdhp18.tiff",
     width=1000,height=800,res=200)
plotClusters(clusters = primary.hetero,
             #genes = 'CN_altered',
             tsb = "DHP18")
dev.off()

brain.math = math.score(merge_brain2,vafCol = 'tumor_f')
primary.math = math.score(merge_primary2,vafCol = 'tumor_f')

brain.math = brain.math[order(brain.math$Tumor_Sample_Barcode),]
primary.math = primary.math[order(primary.math$Tumor_Sample_Barcode),]

brain.math$fc = brain.math$MATH/primary.math$MATH

total.math = rbind(primary.math,brain.math[,c(1:3)])
total.math$type = c(rep("Primary",75),rep("BrM",75))
total.math = total.math[,-1]
median = data.frame(Sample = primary.math$Tumor_Sample_Barcode,
                    Primary = primary.math$MedianAbsoluteDeviation,
                    BrM = brain.math$MedianAbsoluteDeviation)
mathtable = data.frame(Sample = primary.math$Tumor_Sample_Barcode,
                       Primary = primary.math$MATH,
                       BrM = brain.math$MATH)

library(ggpubr)
tiff(filename = "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/math_score.tiff",width = 1600,height = 1600,res = 300)
ggpaired(mathtable,cond1 = "Primary",cond2="BrM",
         color = "condition",line.color = "gray",
         line.size = 0.4,palette = c("#1b98e0", "red")) +
  stat_compare_means(paired=TRUE) + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "None") +
  ylab("Mutant-allele Tumor Heterogeneity")
dev.off()

tiff(filename = "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/median_score.tiff",width = 1600,height = 1600,res = 300)
ggpaired(median,cond1 = "Primary",cond2="BrM",
         color = "condition",line.color = "gray",
         line.size = 0.4,palette = c("#1b98e0", "red")) +
  stat_compare_means(paired=TRUE) + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "None") +
  ylab("Median Absolute Deviation") + ylim(c(0,24))
dev.off()

library(dplyr)
library(data.table)
total.math2 = melt(total.math,id="type")
total.math$type <- factor(total.math$type,
                               levels = c('Primary','BrM'))
tiff(filename = "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/math_score.tiff",width = 2000,height = 1600,res = 300)
ggplot(total.math2, aes(x = variable, y = value, color = type)) +
  geom_boxplot(width = 0.7) + 
  stat_compare_means(label = "p.signif",paired=TRUE) + 
  theme_classic2() +
  ylab("Mutant-allele Tumor Heterogeneity") +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position="top",
        axis.text.x = element_text(angle = 0, hjust=0.5,size = 15),
        axis.text.y = element_text(size = 15)) +
  scale_color_manual(values = c("#1b98e0", "red")) + scale_y_log10()
dev.off()

write.csv(brain.math,file = "brainmath.csv")
write.csv(primary.math,file="primarymath.csv")
library(ggpubr)
mad_table <- data.frame(primary = primary.math$MedianAbsoluteDeviation, BrM = brain.math$MedianAbsoluteDeviation)
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/brain_ith1.tiff",
     width=1000,height=1000,res=200)
ggpaired(mad_table,cond1 = "primary",cond2 = "BrM", fill = "condition", palette = "jco") +
  stat_compare_means() + ylab("MAF Median Absolute Deviation")
dev.off()
math_table <- data.frame(primary = primary.math$MATH, BrM = brain.math$MATH)
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/brain_ith2.tiff",
     width=1000,height=1000,res=200)
ggpaired(math_table,cond1 = "primary",cond2 = "BrM", fill = "condition", palette = "jco") +
  stat_compare_means() + ylab("Mutant-allele Tumor Heterogeneity")
dev.off()

#threshhold 0.3, confidence level 99%
all.lesions <- "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/cnv/354479/brain.all_lesions.conf_99.txt"
amp.genes <- "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/cnv/354479/brain.amp_genes.conf_99.txt"
del.genes <- "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/cnv/354479/brain.del_genes.conf_99.txt"
scores.gis <- "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/cnv/354479/brain.scores.gistic"

clinical = read.table("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/Lung_Clinical_Total.txt",sep="\t",header=TRUE)
clinical$Overall_Survival_Status <- 1
clinical$Overall_Survival_Status[which(is.na(clinical$Age.of.Death))] <- 0
clinical$time = 365*(clinical$Age.of.Death-clinical$Age.of.Diagnosis.of.Primary.Tumor)
  
brain.gistic = readGistic(gisticAllLesionsFile = all.lesions, 
                         gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, 
                         gisticScoresFile = scores.gis, isTCGA = FALSE)

gisticChromPlot(gistic = brain.gistic,
                markBands = "all",
                fdrCutOff = 0.1,
                maf = merge_brain,
                mutGenes = c("KRAS","ERBB2","EGFR","PIK3CA","CDKN2A"),
                y_lims = c(-0.5,0.5))

#for i in *maf; do j=$(basename "$i" .funcotated.maf); sed 's/__UNKNOWN__/'$j'/g' $i > $j.new.maf; done
#cat *new.maf >> merge.maf
merge_brain <-NULL
merge_brain <- read.maf("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/Brain/merge.maf",
                  removeDuplicatedVariants = TRUE,
                  clinicalData = clinical,
                  gisticAllLesionsFile = all.lesions,
                  gisticAmpGenesFile = amp.genes,
                  gisticDelGenesFile = del.genes,
                  gisticScoresFile = scores.gis)
#droplevels()
#threshold 2.0 to 0.7
all.lesionsp <- "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/cnv/354398/primary.all_lesions.conf_99.txt"
amp.genesp <- "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/cnv/354398/primary.amp_genes.conf_99.txt"
del.genesp <- "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/cnv/354398/primary.del_genes.conf_99.txt"
scores.gisp <- "/Users/jianlan/Desktop/Lab/Brain_meta/Lung/cnv/354398/primary.scores.gistic"
primary.gistic = readGistic(gisticAllLesionsFile = all.lesionsp, 
                         gisticAmpGenesFile = amp.genesp, gisticDelGenesFile = del.genesp, 
                         gisticScoresFile = scores.gisp, isTCGA = FALSE)

gisticChromPlot(gistic = primary.gistic, markBands = "all")
gisticBubblePlot(gistic = primary.gistic)
gisticOncoPlot(gistic = primary.gistic,
               sortByAnnotation = TRUE, top = 20)


merge_primary <- NULL
merge_primary <- read.maf("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/Primary/merge.maf",
                          removeDuplicatedVariants = TRUE,
                          clinicalData = clinical,
                          gisticAllLesionsFile = all.lesionsp,
                          gisticAmpGenesFile = amp.genesp,
                          gisticDelGenesFile = del.genesp,
                          gisticScoresFile = scores.gisp)

merge_adeno <- read.maf("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/Primary/adenocarcinoma/adenomerge.maf",
                        removeDuplicatedVariants = TRUE)


#plot summary


dev.off()
#cnv analysis
brain.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, 
                         gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = FALSE)
primary.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, 
                            gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = FALSE)
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/brain_CNV.tif",
     width=2000,height=1500,res=300)
gisticChromPlot(gistic = brain.gistic)
dev.off()

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/primary_CNV.tif",
     width=2000,height=1500,res=300)
gisticChromPlot(gistic = primary.gistic, markBands = "all",cytobandTxtSize = 1,mutGenesTxtSize = 0.2,
                txtSize = 0.5)
dev.off()

#signiture analysis
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
primary.tnm = trinucleotideMatrix(maf = merge_primary2, 
                               prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
brain.tnm = trinucleotideMatrix(maf = merge_brain2, 
                                  prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
plotApobecDiff(tnm = primary.tnm, maf = merge_brain2, pVal = 0.2)
library('NMF')
primary.sign = estimateSignatures(mat = primary.tnm, nTry = 6)
primary.sig = extractSignatures(mat = primary.tnm, n = 3)
primary.og30.cosm = compareSignatures(nmfRes = primary.sig, sig_db = "legacy")
primary.v3.cosm = compareSignatures(nmfRes = primary.sig, sig_db = "SBS")

brain.sign = estimateSignatures(mat = brain.tnm, nTry = 6)
brain.sig = extractSignatures(mat = brain.tnm, n = 4)
brain.og30.cosm = compareSignatures(nmfRes = brain.sig, sig_db = "legacy")
brain.v3.cosm = compareSignatures(nmfRes = brain.sig, sig_db = "SBS")

library('pheatmap')
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/primary_sig_sim.tif",
     width=2000,height=1000,res=300)
pheatmap::pheatmap(mat = primary.v3.cosm$cosine_similarities, 
                   cluster_rows = FALSE, main = "Primary cosine similarity against COSMIC signatures")
dev.off()

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/brain_sig.tiff",
     width=2000,height=1500,res=300)
maftools::plotSignatures(nmfRes = brain.sign, title_size = 1.2, sig_db = "SBS")
dev.off()

library("barplot3d")
#Visualize first signature
brainsig1 = brain.sig$signatures[,1]
barplot3d::legoplot3d(contextdata = brainsig1, labels = FALSE, scalexy = 0.01, sixcolors = "sanger", alpha = 0.5)

#mutation signiture analysis
library("gridExtra")
library("MutationalPatterns")
library(RColorBrewer)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome)
library(BSgenome.Celegans.UCSC.ce2)
library("NMF")
brain_mat <- t(brain.tnm[["nmf_matrix"]])
brain_mat = brain_mat + 0.0001
brain_estimate = nmf(brain_mat, rank=2:6, method="brunet", nrun=1000, seed=123456)
plot(brain_estimate)