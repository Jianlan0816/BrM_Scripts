library(GenVisR)
braindf = data.frame(merge_brain2@data)
colnames(braindf)[3] = "Tumor_Sample_Barcode"
waterfall(braindf, mainRecurCutoff = 0.06, 
          maxGenes = 20)

braincnvdf = data.frame(brain.gistic@data)
braincnvdf = braincnvdf[!grepl("001",braincnvdf$Tumor_Sample_Barcode),]
braincnvdf = braincnvdf[which(braincnvdf$Hugo_Symbol %in% refgenes$V3),]

######https://cancer.sanger.ac.uk/census
#######https://cancer.sanger.ac.uk/cosmic/census/tables?name=dels
censuscosmic= read.table("~/Downloads/Census_allSat Dec 11 10_36_25 2021.csv",
                         sep = ",",header = TRUE)
census_amp = read.table("~/Downloads/Census_ampSat Dec 11 16_58_51 2021.csv",
                        sep=",", header = TRUE)
census_del = read.table("~/Downloads/Census_delsSat Dec 11 16_59_19 2021.csv",
                        sep=",", header = TRUE)
cosmicgenes = union(census_amp$Gene.Symbol, census_del$Gene.Symbol)
braincnvdf = braincnvdf[which(braincnvdf$Hugo_Symbol %in% censuscosmic$Gene.Symbol),]
braincnvdf = braincnvdf[which(braincnvdf$Hugo_Symbol %in% cosmicgenes),]
braincnvdf2 = braincnvdf[,c(1:3)]
braincnvdf2 = dplyr::distinct(braincnvdf2)
braincnvdf2 = braincnvdf2[which(braincnvdf2$Hugo_Symbol %in% refgenes$V3),]
######plot heatmap of cnv######
library(reshape2)
braincnvdf2$Variant_Classification = sub("^(Amp).*",1,braincnvdf2$Variant_Classification)
braincnvdf2$Variant_Classification = sub("^(Del).*",-1,braincnvdf2$Variant_Classification)
braincnvdf2 = braincnvdf2[-c(129,188,291),]
braincnvdf2 = spread(braincnvdf2, Hugo_Symbol, Variant_Classification)
braincnvdf2[is.na(braincnvdf2)] <- 0
rownames(braincnvdf2) = braincnvdf2$Tumor_Sample_Barcode
braincnvdf2 = braincnvdf2[,-1]                                          
i <- c(1, 19) 
braincnvdf2[ , i] <- apply(braincnvdf2[ , i], 2,            # Specify own function within apply
                          function(x) as.numeric(as.character(x)))
braincnvdf3 = as.matrix(sapply(braincnvdf2, as.numeric))
rownames(braincnvdf3) = rownames(braincnvdf2)

clinicsub = clinic[which(clinic$Tumor_Sample_Barcode %in% rownames(braincnvdf3)), ]
mix_ha <- HeatmapAnnotation(Cohort = clinicsub$Cohort, 
                            col = list(Cohort = fabcolors))
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/brmcnvheatmap.tiff",width=3500,height=2000,res=300)
Heatmap(t(braincnvdf3),cluster_rows = T,
        cluster_columns = T,
        rect_gp = gpar(type="none"),
        cell_fun = cell_fun,
        col = col_fun,
        top_annotation = mix_ha
)
dev.off()

lungcnvdf = data.frame(lung.gistic@data)
lungcnvdf = lungcnvdf[!grepl("001",lungcnvdf$Tumor_Sample_Barcode),]
lungcnvdf = lungcnvdf[which(lungcnvdf$Hugo_Symbol %in% refgenes$V3),]
lungcnvdf = lungcnvdf[which(lungcnvdf$Hugo_Symbol %in% censuscosmic$Gene.Symbol),]
lungcnvdf = lungcnvdf[which(lungcnvdf$Hugo_Symbol %in% cosmicgenes),]

######plot heatmap of cnv######
lungcnvdf2 = lungcnvdf[,c(1:3)]
lungcnvdf2 = dplyr::distinct(lungcnvdf2)
library(reshape2)
lungcnvdf2$Variant_Classification = sub("^(Amp).*", 1,lungcnvdf2$Variant_Classification)
lungcnvdf2$Variant_Classification = sub("^(Del).*", -1,lungcnvdf2$Variant_Classification)

#lungcnvdf2 = acast(lungcnvdf2, Hugo_Symbol~Tumor_Sample_Barcode, value.var='Variant_Classification',fill=0)
library(tidyr)
lungcnvdf2 = lungcnvdf2[-151,]
lungcnvdf2 = spread(lungcnvdf2, Hugo_Symbol, Variant_Classification)
lungcnvdf2[is.na(lungcnvdf2)] <- 0
rownames(lungcnvdf2) = lungcnvdf2$Tumor_Sample_Barcode
lungcnvdf2 = lungcnvdf2[,-1]
i <- c(1, 21) 
lungcnvdf2[ , i] <- apply(lungcnvdf2[ , i], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))
lungcnvdf3 = as.matrix(sapply(lungcnvdf2, as.numeric))
rownames(lungcnvdf3) = rownames(lungcnvdf2)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
fabcolors = RColorBrewer::brewer.pal(n = 5,name = 'Set1')
names(fabcolors) = c("Brastianos 2015", "Fukumura 2021", "Liu 2021", 
                     "SYSUCC","Paik 2015")
clinicsub = clinic[which(clinic$Tumor_Sample_Barcode %in% rownames(lungcnvdf2)), ]
mix_ha <- HeatmapAnnotation(Cohort = clinicsub$Cohort, 
                            col = list(Cohort = fabcolors))

cell_fun = function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y, width = width * 0.95, height = height*0.95, 
            gp = gpar(col = "grey", fill = fill, lty = 1, lwd = 0.5))
}
col_fun = colorRamp2(c(-1,0,1),c("blue","white","red"))
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/lungcnvheatmap.tiff",width=3500,height=2000,res=300)
Heatmap(t(lungcnvdf3),cluster_rows = T,
        cluster_columns = T,
        rect_gp = gpar(type="none"),
        cell_fun = cell_fun,
        col = col_fun,
        top_annotation = mix_ha
        )
dev.off()

clinical <- melt(clinic, id.vars = c("Tumor_Sample_Barcode"))

cnSpec(LucCNseg, genome = "hg19")
brain_cn = read.table("~/Desktop/Lab/Brain_meta/Lung/cnv/brain/brain_beforemerge.seg.txt",
                      header = T)
brain_cn = brain_cn[!grepl("001",brain_cn$ID),]
brain_cn$ID <- sub("^.*", "BrM", brain_cn$ID)

brain_cn$ID <- sub("^(Lung).*", "Fukumura", brain_cn$ID)
brain_cn$ID <- sub("^(DHP).*", "SYSUCC", brain_cn$ID)
brain_cn$ID <- sub("^(PB).*", "Brastinos", brain_cn$ID)
write.table(brain_cn,"~/Desktop/Lab/Brain_meta/Lung/cnv/brain/brainall_beforemerge_liuremove.seg.txt",
            col.names = TRUE,quote = FALSE,row.names = FALSE,sep = "\t")

brain_cn = brain_cn[,c(2:6,1)]
colnames(brain_cn) = colnames(LucCNseg)
brain_cn$segmean = brain_cn$segmean + 1
pdf(file="~/Desktop/Lab/Brain_meta/Lung/plots/cnv.pdf")
cnSpec(brain_cn, genome = "hg19",CNscale = "absolute")
dev.off()

primary_cn = read.table("~/Desktop/Lab/Brain_meta/Lung/cnv/primary/primary_beforemerge.seg.txt",
                        header = T,sep = "\t")

primary_cn = primary_cn[!grepl("001",primary_cn$Sample),]
#primary_cn$Sample <- sub("^.*", "Primary", primary_cn$Sample)
write.table(primary_cn,"~/Desktop/Lab/Brain_meta/Lung/cnv/primary/primary_beforemerge_liuremove.seg.txt",
            col.names = TRUE,quote = FALSE,row.names = FALSE,sep = "\t")
