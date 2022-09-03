setwd
setwd("~/Desktop/Lab/Brain_meta/Lung_own/RNA/count/empty/")
files<- dir(pattern="*.txt")
files

count<-read.table(files[1], sep='\t',header = F)
count1<-read.table(files[2], sep='\t',header = F)
total = NULL
total<-merge(count, count1, by="V1")

for(i in 3:46){
  c<-read.table(files[i], sep='\t',header = F)
  total<-merge(total, c, by="V1")
}

library(stringr)
samplename<-str_remove_all(files, ".count.txt")
colnames(total)[2:47]=samplename[1:46]

colnames(total)[1]="ENSEMBL"

library(edgeR)
library(cpm)
library(dplyr)
library(openxlsx)
library(lumi)



total <- total[-c(1:5),]

total$ENSEMBL = sub("\\..*","",total$ENSEMBL)
total1=total
write.xlsx(total,"Count.xlsx")

library("readxl")
my_data <- read_excel("Count.xlsx")

library(org.Hs.eg.db)
annotation<-AnnotationDbi::select(org.Hs.eg.db, keys = as.character(my_data$ENSEMBL),columns = "SYMBOL", keytype = "ENSEMBL")

annotation2<-AnnotationDbi::select(org.Hs.eg.db, keys = as.character(total$GENE_ID) ,columns = "GENENAME", keytype = "ENSEMBL")
colnames(my_data)[1]="ENSEMBL"
totalnew1<-merge(annotation,my_data,by="ENSEMBL")
totalnew2<-merge(annotation2,totalnew1,by="ENSEMBL")

k=3
while (k<37){
  total[,k]=as.numeric(as.character(totalnew1[,k]))
  k=k+1
}

totalnew3 <- totalnew1[,-c(1)]
totalnew3 <- aggregate(.~SYMBOL,data=totalnew3, FUN=sum)

l=c(1:47)
k=l[grepl("B",colnames(totalnew3))]
m=l[!grepl("B",colnames(totalnew3))]
BrMvsP<-totalnew3[,c(k,m)]

BrM <- BrMvsP[,c(1:9)]
P <- BrMvsP[,c(10:14)]

row.names(BrMvsP)=BrMvsP$SYMBOL
BrMvsP = BrMvsP[,-1]

BrM = as.matrix(BrM)
P = as.matrix(P)
dataNormtest <- TCGAanalyze_Normalization(tabDF = cbind(BrM, P),
                                      geneInfo = TCGAbiolinks::geneInfo,
                                      method = "gcContent") #18323   672
dataFilttest <- TCGAanalyze_Filtering(tabDF = dataNormtest,
                                  method = "quantile",
                                  qnt.cut =  0.25)  # 13742   672
save(dataFilttest, file = paste0("BrM_P_Norm_IlluminaHiSeq.rda"))

dataDEGstest <- TCGAanalyze_DEA(mat1 = BrM2,
                            mat2 = P,
                            Cond1type = "BrM",
                            Cond2type = "P",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")
# Number of differentially expressed genes (DEG)

myCPM <- cpm(BrMvsP)

thresh <- myCPM > 0.5

keep <- rowSums(thresh) >= 2

counts.keep <- BrMvsP[keep,]
group<-c(rep(1,17),rep(2,17))
dgeObj <- DGEList(counts.keep,group=factor(group))
dgeObj <- calcNormFactors(dgeObj)

dgeObj <- estimateCommonDisp(dgeObj)
dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj)
ex <- exactTest(dgeObj, pair = c(2,1), dispersion = 0.1)

result <- ex$table
result$FDR <- p.adjust(result$PValue,method = "fdr")
result.cpm <- merge(result,myCPM[,c(1:34)],by="row.names")
colnames(result.cpm)[c(6:7)] <- paste(colnames(result.cpm)[c(6:7)],".CPM",by="")
result.cpm <- result.cpm[order(result.cpm$PValue),]

tiff(filename = "MD_OralTvsC.tiff", width = 6000, height = 5000, res = 600,compression = "lzw")
plotMD(ex)
abline(h=0,col='red')
dev.off()

colnames(result.cpm)[1]="SYMBOL"
bothdfs <- merge(result.cpm, totalnew3,by="SYMBOL")
bothdfs2 <- merge(bothdfs,totalnew2[,c(1,2,3)],by="SYMBOL")

write.table(myCPM,"totalcpm2.txt",sep="\t",quote = F,row.names = T,col.names = T)

result.cpm$FDR<-result.cpm$PValue
row.names(result.cpm)<-result.cpm$SYMBOL


library(ggplot2)
tiff(filename = "Volcano_OralTvsC.tiff", width = 6000, height = 5000, res = 600,compression = "lzw")
EnhancedVolcanoEdgeR(toptable=result.cpm, NominalCutoff=0.05, AdjustedCutoff=0.01, LabellingCutoff=6, FCCutoff=3, main="Oral_TumorvsControl")
dev.off()


logcounts <- cpm(dgeObj,log=TRUE)

library(ggplot2)
library(DESeq2)
dataset = ExpressionSet(assayData = data.matrix(logcounts))
boxplot(dataset, color="red")

tiff(filename = "Oral_boxplot.tiff", width = 3000, height = 5000, res = 600)
boxplot(dataset)
dev.off()

#SampleRelation
group<-c(1,2,3,4,5,6)
tiff(filename = "Oral_samplerelation.tiff", width = 4000, height = 4000, res = 600)
plot(dataset,what='sampleRelation', method="mds",pch = 15, col=c("black","red","green","yellow","blue","purple"))
legend("bottom",
       c("2033-Oral_A5","2033-Oral_A6","Bots_HN_D29","Bots_HN_F103","Oral_con_100C","Oral_con_102C"),
       pch = 15,
       col = c("black","red","green","yellow","blue","purple"))
dev.off()

#Heatmap
library(pheatmap)
rna.dist <- data.matrix(dist(t(logcounts),diag = T,upper = T))
png(filename="Oral_heatmap.png",width = 3000, height = 2000, res=600)
pheatmap(rna.dist,show_rownames = T,show_colnames = T,fontsize = 6)
dev.off()

#deseq2
#OralTvsC<-totalnew2[c(5,6,11,12,13,14)]
total_Oral <- read.xlsx("/Users/jianlan/Desktop/Lab/Hiseq4000/htseq/Oral/OralTvsC.xlsx")
OralTvsC <- total_Oral[,c(1,10:15)]
OralTvsC2 <- OralTvsC[!duplicated(OralTvsC$SYMBOL), ]
rownames(OralTvsC2)<-OralTvsC2[,1]
OralTvsC2<-OralTvsC2[,-1]

#cervical
total_Cerv <- read.xlsx("/Users/jianlan/Desktop/Lab/Hiseq4000/htseq/Cervical/CervicalTvsC.xlsx")
CervTvsC <- total_Cerv[,c(1,10:15)]
CervTvsC2 <- CervTvsC[!duplicated(CervTvsC$SYMBOL), ]
rownames(CervTvsC2)<-CervTvsC2[,1]
CervTvsC2<-CervTvsC2[,-1]

library(GSVA)
#------------ssgsea---------------
score <- gsva(myCPM, pathways.hallmark, method="ssgsea", min.sz=1)
write.csv(score, "reactome_ssgsea_score.csv", row.names=T)

#only for plotting
score.t = t(scale(t(score), center=T, scale=T))
score.t = data.frame(score.t, check.names = F)
write.csv(score.t, "reactome_ssgsea_scorescale.csv", row.names=T)
#based on difference in score between avg of B and avg of L
score.t2=score.t
score.t2$avg1=rowMeans(score.t2[1:21])
score.t2$avg2=rowMeans(score.t2[22:43])
score.t2$diff=score.t2$avg2-score.t2$avg1
score.t2=score.t2[order(score.t2$diff,decreasing = T),]

annotation <- data.frame(Group=rep(c("BrM","Primary"),c(23,23)))
rownames(annotation) <- colnames(score.t2)[1:46]

annotation_row = data.frame(
  pathway = factor(rep(c("Path-BrM", "Path-Primary"), c(20, 20))))
rownames(annotation_row) = rownames(score.t2[c(1585:1604,1:20),])

tiff(file="ssgsea_heatmap2.tiff", width = 5000, height = 3600, res=400, compression = "lzw")

library(pheatmap)
pheatmap(
  score.t2[c(1:20,1585:1604),c(43:1)],
  #cluster_rows=F,
  cluster_cols=F,
  annotation_row = annotation_row,
  annotation_col = annotation,
  annotation_legend = TRUE,
  #annotation_colors = list("Group"=c("BrM/Primary"="Black", "Primary/BrM"="grey")),
  show_colnames=T,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  border_color=NA,
  fontsize=6,
  cellwidth = 8,
  cellheight = 10,
  main = "BrM to Primary ssGSEA Score Difference in Reactome"
)

dev.off()

library(DESeq2)
#design <- data.frame(Intercept=rep(1,12), Group = c("Cancer","AD","nonAD","Cancer","Cancer","AD","nonAD","AD","nonAD","Cancer","AD","nonAD"))
design<-data.frame(Intercept=rep(1,46),Group = c(rep("B",23),rep("A",23)))
                                         


dds <- DESeqDataSetFromMatrix(countData = BrMvsP,
                              colData = design,
                              design= ~ Group)

dds <- DESeq(dds)


resultsNames(dds) # lists the coefficients


res <- results(dds, name="Group_B_vs_A")

res[,7]=rownames(res)
colnames(res)[7]="SYMBOL"

#https://stephenturner.github.io/deseq-to-fgsea/
#fgsea
res2 <- as.data.frame(res) %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))

library(fgsea)
library(dplyr)
library(tibble)
#deframe() converts two-column data frames to a named vector or list, using the first column as name and the second column as value.
ranks <- deframe(res2)
head(ranks, 20)
res3 = res2[which(res2$stat < -2 | res2$stat>2),]
rankss = deframe(res3)
#-------------------  4.2 EA: enrichment analysis             --------------------
ansEAtest <- TCGAanalyze_EAcomplete(TFname="DEA genes BrM vs P", RegulonList = names(rankss))
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/GO.tiff",width=2000,height=2000,res=200)
TCGAvisualize_EAbarplot(tf = rownames(ansEAtest$ResBP),
                        filename = NULL,
                        GOBPTab = ansEAtest$ResBP, 
                        nRGTab = names(rankss),
                        nBar = 20)
dev.off()
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/CC.tiff",width=2000,height=2000,res=200)
TCGAvisualize_EAbarplot(tf = rownames(ansEAtest$ResBP),
                        filename = NULL,
                        GOCCTab = ansEAtest$ResCC,
                        nRGTab = rownames(rankss),
                        nBar = 20)
dev.off()
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/MF.tiff",width=2000,height=2000,res=200)
TCGAvisualize_EAbarplot(tf = rownames(ansEAtest$ResBP),
                        filename = NULL,
                        GOMFTab = ansEAtest$ResMF, 
                        nRGTab = rownames(rankss),
                        nBar = 20)
dev.off()
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/Pathway.tiff",width=2000,height=2000,res=200)
TCGAvisualize_EAbarplot(tf = rownames(ansEAtest$ResBP),
                        filename = NULL,
                        PathTab = ansEAtest$ResPat,
                        nRGTab = rownames(rankss),
                        nBar = 20)
dev.off()
#The gmtPathways() function will take a GMT file you downloaded from MSigDB and turn it into a list. 
pathways.hallmark <- gmtPathways("/Users/jianlan/Downloads/c7.immunesigdb.v7.4.symbols.gmt")

# Look at them all if you want (uncomment)
# pathways.hallmark

# Show the first few pathways, and within those, show only the first few genes. 
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange((padj))

set_lists_to_chars <- function(x) {
  if(class(x) == 'list') {
    y <- paste(unlist(x[1]), sep='', collapse=', ')
  } else {
    y <- x 
  }
  return(y)
}
fgseaRes2=fgseaRes[,-9]
for (i in 1:nrow(fgseaRes2)){fgseaRes2$leadingEdge2[i] = paste(unlist(fgseaRes2$leadingEdge[i]), sep='', collapse=', ')}

new_frame <- data.frame(lapply(fgseaResTidy, set_lists_to_chars), stringsAsFactors = F)
new_frame=new_frame[,-8]

write.csv(new_frame, file='fgsea.csv')

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

#Plot the normalized enrichment scores. Color the bar indicating whether or not the pathway was significant:
library(ggplot2)
tiff(file="fgsea_barplot_NES.tiff",width = 2000, height = 900, res=150,compression = "lzw")
ggplot(fgseaResTidy[c(1:40),], aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES<0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=NULL) + 
  theme_minimal()
dev.off()

topPathwaysUp <- fgseaRes[NES > 0][head(order(NES,padj,decreasing = T), n=20), pathway]
topPathwaysDown <- fgseaRes[NES < 0][head(order(NES,padj), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))


tiff(file="topgseatable.tiff",width = 2000, height = 800, res=110,compression = "lzw")
plotGseaTable(pathways.hallmark[topPathwaysUp], ranks, fgseaRes, 
              gseaParam = 0.5,colwidths = c(19, 2, 0.6, 1, 1),)
dev.off()

plotEnrichment(pathways.hallmark[["REACTOME_DISEASE"]], ranks)
plotEnrichment(pathways.hallmark[["REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION"]],ranks)
plotEnrichment(pathways.hallmark[["REACTOME_POST_TRANSLATIONAL_PROTEIN_MODIFICATION"]],ranks)

plotEnrichment(pathways.hallmark[["REACTOME_RELEASE_OF_HH_NP_FROM_THE_SECRETING_CELL"]],ranks)

for (i in 1:20) {
  tiff(filename=paste("./fgsea_plots/",topPathways[i],".tiff",sep = ""),width = 3000, height = 2000, res=600,compression = "lzw")
  print(plotEnrichment(pathways.hallmark[[topPathways[i]]], ranks) + labs(title=topPathways[i]))
  dev.off()
}

#plotenrichment sample
data("examplePathways")
data("exampleRanks")
plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],exampleRanks)

# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="Group_groupB_vs_groupA", type="apeglm")

res1 <- lfcShrink(dds2, coef="Group_groupBD_vs_groupAC", type="apeglm")

#MAPlot VolcanoPlot
write.csv(as.data.frame(res), 
          file="GroupAvsB.csv")

write.csv(as.data.frame(res2),
          file="GroupACvsBD.csv")
dev.new()

plotMA(res,ylim=c(-0.000005,0.000005))

plotMA(res2,ylim=c(-5,5))

library(EnhancedVolcano)
# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  res$log2FoldChange < -2.5, 'royalblue',
  ifelse(res$log2FoldChange > 2.5, 'red',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'

tiff(filename="volcano_BrMvsP2.0.tiff",width = 2000, height = 2000, res=300)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = rownames(res)[which(names(keyvals) %in% c('high', 'low'))],
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 1.5,
                FCcutoff = 2.5, 
                pCutoff = 10e-4,
                colCustom = keyvals,
                xlim = c(-8, 8),
                title = NULL,
                subtitle=NULL)
dev.off()

tiff(filename="volcano_ACvsBD.tiff",width = 5000, height = 5000, res=600)
EnhancedVolcano(res2,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                xlim = c(-5, 8))
dev.off()

#p-value
qwe=which(is.na(ex$table$PValue))
new_res_p=result$PValue[-qwe]
head(new_res_p)
tiff(filename="Oral_pvalue.tiff",width = 4000, height = 3000, res=600)
hist(new_res_p,breaks = 50,col="red")
hist(result$PValue,breaks = 50,col="red")
dev.off()

library(clusterProfiler)
egtest = as.data.frame(bitr(res3$SYMBOL,
                        fromType="SYMBOL",
                        toType="ENTREZID",
                        OrgDb="org.Hs.eg.db"))
egtest <- egtest[!duplicated(egtest$SYMBOL),]
egtest <- egtest[order(egtest$SYMBOL,decreasing=FALSE),]
all(egtest$SYMBOL == res3$SYMBOL)
res_df = as.data.frame(res)
res_df = subset(res_df,res_df$SYMBOL %in% egtest$SYMBOL)
egtest2 = as.data.frame(bitr(res_df$SYMBOL,
                            fromType="SYMBOL",
                            toType="ENTREZID",
                            OrgDb="org.Hs.eg.db"))
res_df = cbind(res_df,egtest2$ENTREZID)
colnames(res_df)[8] = "GeneID"
dataDEGsFiltLevel_subtest <- subset(res_df, select = c("GeneID", "log2FoldChange"))
genelistDEGstest <- as.numeric(dataDEGsFiltLevel_subtest$log2FoldChange)
names(genelistDEGstest) <- dataDEGsFiltLevel_subtest$GeneID
# pathway.id: hsa05214 is the glioma pathway
# limit: sets the limit for gene expression legend and color
library(pathview)
hsa05214test <- pathview::pathview(gene.data  = genelistDEGstest,
                               pathway.id = "hsa05214",
                               species    = "hsa",
                               limit = list(gene=as.integer(max(abs(genelistDEGstest)))))
hsa8569 <- pathview::pathview(gene.data  = genelistDEGstest,
                                   pathway.id = "hsa04151",
                                   species    = "hsa",
                                   limit = list(gene=as.integer(max(abs(genelistDEGstest)))))


