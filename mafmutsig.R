data("ToothGrowth")
head(ToothGrowth)

ggbarplot(ToothGrowth, x = "dose", y = "len", add = "mean_se",
          color = "supp", palette = "jco", 
          position = position_dodge(0.8))+
  stat_compare_means(aes(group = supp), label = "p.signif", label.y = 29) +
  coord_flip() +scale_y_continuous(position="right")

df <- data.frame(dose=c("D0.5", "D1", "D2"),
                 len=c(4.2, 10, 29.5))
print(df)

# Basic plot with label outsite
# +++++++++++++++++++++++++++
ggbarplot(df, x = "dose", y = "len",
          label = TRUE, label.pos = "out")
ggbarplot(ToothGrowth, x = "dose", y = "len",
          add = c("mean_se", "dotplot"))

df2 = data.frame(color = c(rep("#845EC2",2),rep("#F6C8FF",2)))

########mut signature########
library(MutationalPatterns)
library(BSgenome)
head(available.genomes())
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
brm.tnm = trinucleotideMatrix(maf = merge_brain2, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
brm_snv = brm.tnm[["APOBEC_scores"]][,c(1:13)]
brm_snv$`C>A_new` = brm_snv$`G>A` + brm_snv$`C>A`
brm_snv$`C>G_new` = brm_snv$`C>G` + brm_snv$`G>C`
brm_snv$`C>T_new` = brm_snv$`G>T` + brm_snv$`C>T`
brm_snv$`T>A_new` = brm_snv$`A>T` + brm_snv$`T>A`
brm_snv$`T>C_new` = brm_snv$`T>C` + brm_snv$`A>C`
brm_snv$`T>G_new` = brm_snv$`T>G` + brm_snv$`A>G`
brm_snv2 = brm_snv[,c(1,14:19)]
colnames(brm_snv2) = str_remove(colnames(brm_snv2),"_new")
brm_snv2 = melt(brm_snv2)
brm_snv2 <- within(brm_snv2, Tumor_Sample_Barcode <- factor(Tumor_Sample_Barcode, levels = names(sort(table(Tumor_Sample_Barcode), decreasing = TRUE))))

p1 = ggplot(brm_snv2, aes(fill = variable, y=value, x=Tumor_Sample_Barcode)) +
  geom_bar(position = "fill", stat = "identity",width = 0.6) +
  theme_classic2() + 
  ylab("Proportion of \nbase change (%)") +
  scale_fill_manual(values=c("#43CA89", "#F3E482", "#B9AD4E",
                              "#817919", "#6A9B7F", "#D9EDDF")) +
  scale_y_continuous(labels = function(x) paste0(x*100, "")) +
  ggtitle("BrM") + 
  theme(axis.title.x=element_blank(), 
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    #axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    text=element_text(size=8),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 0,  # Top margin
                         r = 20,  # Right margin
                         b = 20,  # Bottom margin
                         l = 20)
  )  
p1 

prim.tnm = trinucleotideMatrix(maf = merge_primary2, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
prim_snv = prim.tnm[["APOBEC_scores"]][,c(1:13)]
prim_snv$`C>A_new` = prim_snv$`G>A` + prim_snv$`C>A`
prim_snv$`C>G_new` = prim_snv$`C>G` + prim_snv$`G>C`
prim_snv$`C>T_new` = prim_snv$`G>T` + prim_snv$`C>T`
prim_snv$`T>A_new` = prim_snv$`A>T` + prim_snv$`T>A`
prim_snv$`T>C_new` = prim_snv$`T>C` + prim_snv$`A>C`
prim_snv$`T>G_new` = prim_snv$`T>G` + prim_snv$`A>G`
prim_snv2 = prim_snv[,c(1,14:19)]
colnames(prim_snv2) = str_remove(colnames(prim_snv2),"_new")
prim_snv2 = melt(prim_snv2)
prim_snv2 <- within(prim_snv2, Tumor_Sample_Barcode <- factor(Tumor_Sample_Barcode, 
                                                            levels = names(sort(table(Tumor_Sample_Barcode), decreasing = TRUE))))
p2 = ggplot(prim_snv2, aes(fill = variable, y=value, x=Tumor_Sample_Barcode)) +
  geom_bar(position = "fill", stat = "identity",width = 0.6) +
  theme_classic2() + 
  ylab("Proportion of \nbase change (%)") +
  scale_y_continuous(labels = function(x) paste0(x*100, "")) +
  ggtitle("Primary") +
  scale_fill_manual(values=c("#43CA89", "#F3E482", "#B9AD4E",
                             "#817919", "#6A9B7F", "#D9EDDF")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        text=element_text(size=8),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 0,  # Bottom margin
                             l = 20)
  )  
p2

p_leg = ggplot(prim_snv2, aes(fill = variable, y=value, x=Tumor_Sample_Barcode)) +
  geom_bar(position = "fill", stat = "identity",width = 0.6) +
  theme_classic2() + 
  ylab("Primary") +
  scale_fill_manual(values=c("#43CA89", "#F3E482", "#B9AD4E",
                             "#817919", "#6A9B7F", "#D9EDDF")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 20)
  ) 
leg1 <- get_legend(p_leg)
blank_p <- plot_spacer() + theme_void()
# combine legend 3 & blank plot
leg10 <- plot_grid(leg1, blank_p,
                   nrow = 2,
                   rel_heights = c(1,0.1)
)
final_p <- plot_grid(p2,
                     p1, 
                     nrow = 2,
                     ncol = 1,
                     align = "v",
                     #axis = "t",
                     rel_heights = c(1,1)
)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-2Cdraft.tif",
     width=1500,height=800,res=300)
print(final_p)
dev.off()

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-2C_legend.tif",
     width=300,height=500,res=300)
print(leg10)
dev.off()

library('NMF')
#brm.sign = estimateSignatures(mat = brm.tnm, nTry = 6)
brm.sig = extractSignatures(mat = brm.tnm, n = 10)
brm.v3.cosm = compareSignatures(nmfRes = brm.sig, sig_db = "legacy")
brm_cont = brm.sig[["contributions"]]
brm_cont = data.frame(t(brm_cont))
brm_cont$Tumor_Sample_Barcode = rownames(brm_cont)
brm_cont$COSMIC_1 = brm_cont$Signature_1
brm_cont$COSMIC_2 = brm_cont$Signature_2 
brm_cont$COSMIC_4 = brm_cont$Signature_3 + brm_cont$Signature_6 + brm_cont$Signature_7
brm_cont$COSMIT_7 = brm_cont$Signature_4
brm_cont$COSMIT_6 = brm_cont$Signature_5
brm_cont$Unknown = brm_cont$Signature_8 + brm_cont$Signature_9
brm_cont$COSMIT_13 = brm_cont$Signature_10

brm_cont2 = brm_cont[,c(11,12,13,14,16,15,18,17)]
brm_cont2 = melt(brm_cont2)

brm_cont2 <- within(brm_cont2, Tumor_Sample_Barcode <- factor(Tumor_Sample_Barcode, 
                                                              levels = names(sort(table(Tumor_Sample_Barcode), decreasing = TRUE))))
p1 = ggplot(brm_cont2, aes(fill = variable, y=value, x=Tumor_Sample_Barcode)) +
  geom_bar(position = "fill", stat = "identity",width = 0.6) +
  theme_classic2() + 
  ylab("Signature \n contribution (%)") +
  scale_y_continuous(labels = function(x) paste0(x*100, "")) +
  ggtitle("BrM") +
  scale_fill_manual(
    labels = c("COSMIC_1 (spontaneous)",
               "COSMIC_2 (APOBEC C>T)",
               "COSMIC_4 (smoking)",
               "COSMIC_6 (defective DNA mismatch repair)",
               "COSMIC_7 (UV exposure)",
               "COSMIC_13 (APOBEC C>G)",
               "Unknown"),
                    values=c("#845EC2", "#D65DB1", "#FF6F91",
                             "#FF9671", "#FFC75F", "#F9F871",
                             "#B39CD0")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        text=element_text(size=8),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 0,  # Top margin
                             r = 20,  # Right margin
                             b = 20,  # Bottom margin
                             l = 20)
  ) 
p1

prim.sig = extractSignatures(mat = prim.tnm, n = 10)
prim.v3.cosm = compareSignatures(nmfRes = prim.sig, sig_db = "legacy")
prim_cont = prim.sig[["contributions"]]
prim_cont = data.frame(t(prim_cont))
prim_cont$Tumor_Sample_Barcode = rownames(prim_cont)
prim_cont$COSMIC_1 = prim_cont$Signature_9 + prim_cont$Signature_6
prim_cont$COSMIC_4 = prim_cont$Signature_3 + prim_cont$Signature_7
prim_cont$COSMIT_6 = prim_cont$Signature_5
prim_cont$COSMIT_7 = prim_cont$Signature_4
prim_cont$COSMIT_13 = prim_cont$Signature_10
prim_cont$Unknown = prim_cont$Signature_1 + prim_cont$Signature_8 + prim_cont$Signature_2
prim_cont2 = prim_cont[,c(11:17)]
prim_cont2 = melt(prim_cont2)

prim_cont2 <- within(prim_cont2, Tumor_Sample_Barcode <- factor(Tumor_Sample_Barcode, 
                                                              levels = names(sort(table(Tumor_Sample_Barcode), decreasing = TRUE))))
p2 = ggplot(prim_cont2, aes(fill = variable, y=value, x=Tumor_Sample_Barcode)) +
  geom_bar(position = "fill", stat = "identity",width = 0.6) +
  theme_classic2() + 
  ylab("Signature \n contribution (%)") +
  scale_y_continuous(labels = function(x) paste0(x*100, "")) +
  ggtitle("Primary") +
  scale_fill_manual(
    labels = c("COSMIC_1 (spontaneous)",
               "COSMIC_4 (smoking)",
               "COSMIC_6 (defective DNA mismatch repair)",
               "COSMIC_7 (UV exposure)",
               "COSMIC_13 (APOBEC C>G)",
               "Unknown"),
    values=c("#845EC2", "#FF6F91",
             "#FF9671", "#FFC75F", "#F9F871",
             "#B39CD0")) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        text=element_text(size=8),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 0,  # Bottom margin
                             l = 20)
  ) 
p2

leg1 <- get_legend(p1)
blank_p <- plot_spacer() + theme_void()
# combine legend 3 & blank plot
leg10 <- plot_grid(leg1, blank_p,
                   nrow = 2,
                   rel_heights = c(1,0.1)
)
final_p <- plot_grid(p2,
                     p1, 
                     nrow = 2,
                     ncol = 1,
                     align = "v",
                     #axis = "t",
                     rel_heights = c(1,1)
)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-2Ddraft.tif",
     width=1500,height=800,res=300)
print(final_p)
dev.off()

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-2D_legend.tif",
     width=1000,height=500,res=300)
print(leg10)
dev.off()
