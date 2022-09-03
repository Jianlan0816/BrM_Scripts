setwd("~/Desktop/Lab/Brain_meta/Lung/RNA/totalrna/")
data = read.table("combatTotalRNA.txt", header = TRUE)

####DE ANALYSIS####
library(DESeq2)
design<-data.frame(Intercept=rep(1,64),Group = c(rep("B",32),rep("A",32)))
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = design,
                              design= ~ Group)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Group_B_vs_A")
res[,7]=rownames(res)
colnames(res)[7]="SYMBOL"
write.table(res,file = "TotalDEanalysis.txt",sep="\t",quote = F,row.names = F,col.names = T)

#######volcano plot########
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
tiff(filename="volcano_BrMvsP.tiff",width = 2000, height = 2000, res=300)
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

library(fgsea)
library(dplyr)
library(tibble)
############fgsea#############
res2 <- as.data.frame(res) %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)
res3 = res2[which(res2$stat < -2 | res2$stat>2),]
rankss = deframe(res3)

#The gmtPathways() function will take a GMT file you downloaded from MSigDB and turn it into a list. 
pathways.hallmark <- gmtPathways("/Users/jianlan/Downloads/c2.cp.reactome.v7.4.symbols.gmt")
pathways.hallmark2 <- gmtPathways("/Users/jianlan/Downloads/c2.cgp.v7.4.symbols.gmt")

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange((padj))
fgseaRes2=fgseaRes
for (i in 1:nrow(fgseaRes2)){fgseaRes2$leadingEdge2[i] = paste(unlist(fgseaRes2$leadingEdge[i]), 
                                                               sep='', collapse=', ')}
fgseaRes2 = fgseaRes2[,-8]
fgseaRes3 = fgseaRes2[order(fgseaRes2$padj),]
write.csv(fgseaRes3, file='~/Desktop/Lab/Brain_meta/Lung/RNA/totalrna/fgsea_cpg.csv')

topPathwaysUp <- fgseaRes[NES > 0][head(order(NES,padj,decreasing = T), n=20), pathway]
topPathwaysDown <- fgseaRes[NES < 0][head(order(NES,padj), n=20), pathway]

tiff(file="reactome_topgseatable.tiff",width = 2000, height = 800, res=110,compression = "lzw")
plotGseaTable(pathways.hallmark[topPathwaysUp], ranks, fgseaRes, 
              gseaParam = 0.5,
              colwidths = c(16, 2, 0.6, 1, 1)
)
dev.off()

tiff(file="cpg_topgseatable.tiff",width = 2000, height = 800, res=110,compression = "lzw")
plotGseaTable(pathways.hallmark2[topPathwaysUp], ranks, fgseaRes, 
              gseaParam = 0.5,
              colwidths = c(4, 2, 0.6, 1, 1)
)
dev.off()
