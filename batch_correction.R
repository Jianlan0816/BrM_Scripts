setwd("~/Desktop/Lab/Brain_meta/Lung/RNA/Fukumara/")
library(ggplot2)
library(preprocessCore)
library(sva)
library(DESeq2)

data = read.csv("~/Desktop/Lab/Brain_meta/Lung/RNA/Fukumara/totalRNA.csv")
rownames(data) = data$X
data = data[,-1]

file = file[which(file$GeneID %in% data$X),]
data = data[which(data$X %in% file$GeneID),]

file = file[order(file$GeneID),]
data = data[order(data$X),]

total = cbind(data,file[,-1])
total = total[,-1]

l=c(1:150)
k=l[grepl("B",colnames(total))]
m=l[!grepl("B",colnames(total))]
BrMvsP<-total[,c(k,m)]

i=1
while (i<151){
  BrMvsP[,i]=as.numeric(as.character(BrMvsP[,i]))
  i=i+1
}

Tumor_type = c(rep("Brm",75), rep("Primary", 75))
Batch = c(rep("Batch1", 9), rep("Batch2", 23), rep("Batch3",43),
          rep("Batch1",9), rep("Batch2", 23), rep("Batch3", 43))

design<-data.frame(Tumor_type=Tumor_type,
                   Batch = Batch)

dds <- DESeqDataSetFromMatrix(countData = BrMvsP,
                              colData = design,
                              design= ~ Batch)

######DESEQ2 normalization and transformation######
dds_norm <- vst(dds)
pca_results <-
  plotPCA(
    dds_norm,
    intgroup = c("Batch", "Tumor_type"),
    returnData = TRUE # This argument tells R to return the PCA values
  )
# Plot using `ggplot()` function and save to an object
annotated_pca_plot <- ggplot(
  pca_results,
  aes(
    x = PC1,
    y = PC2,
    # plot points with different colors for each `refinebio_treatment` group
    color = Batch,
    # plot points with different shapes for each `refinebio_disease` group
    shape = Tumor_type
  )
) +
  # Make a scatter plot
  geom_point()

# display annotated plot
tiff("Prebatch.tiff", compression = "zip", 
     res = 200, height = 1000, width = 1000)
annotated_pca_plot
dev.off()

data_combat = ComBat_seq(as.matrix(BrMvsP), batch = Batch, group = Tumor_type)

dds <- DESeqDataSetFromMatrix(countData = data_combat,
                              colData = design,
                              design= ~ Batch)

######DESEQ2 normalization and transformation######
dds_norm <- vst(dds)

pca_results <-
  plotPCA(
    dds_norm,
    intgroup = c("Batch", "Tumor_type"),
    returnData = TRUE # This argument tells R to return the PCA values
  )
# Plot using `ggplot()` function and save to an object
annotated_pca_plot <- ggplot(
  pca_results,
  aes(
    x = PC1,
    y = PC2,
    # plot points with different colors for each `refinebio_treatment` group
    color = Batch,
    # plot points with different shapes for each `refinebio_disease` group
    shape = Tumor_type
  )
) +
  # Make a scatter plot
  geom_point()

# display annotated plot
tiff("Postbatch.tiff", compression = "zip", 
     res = 200, height = 1000, width = 1000)
annotated_pca_plot
dev.off()

write.table(BrMvsP, "~/Desktop/Lab/Brain_meta/Lung/RNA/totalrna/TotalRNA_3.txt", sep = "\t", 
            col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(data_combat, "~/Desktop/Lab/Brain_meta/Lung/RNA/totalrna/combatTotalRNA_3.txt", sep = "\t", 
            col.names = TRUE, row.names = TRUE, quote = FALSE)
dds_norm_norm = assay(dds_norm)
write.table(dds_norm_norm, "~/Desktop/Lab/Brain_meta/Lung/RNA/totalrna/combatTotalRNA_norm_3.txt", sep = "\t", 
            col.names = TRUE, row.names = TRUE, quote = FALSE)

