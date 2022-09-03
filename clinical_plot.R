#plot Muts per Mb
tmb_brain <- merge_brain2@variant.type.summary
tmb_brain <- tmb_brain[order(tmb_brain$Tumor_Sample_Barcode),]

tmb_primary <- merge_primary2@variant.type.summary
tmb_primary <- tmb_primary[order(tmb_primary$Tumor_Sample_Barcode),]

tmb_total =  data.frame("Tumor_Sample_Barcode"=tmb_primary$Tumor_Sample_Barcode,
                        "BrainSNP"=tmb_brain$SNP,
                        "PrimarySNP"=tmb_primary$SNP,
                        'BrainIndel'=tmb_brain$DEL + tmb_brain$INS,
                        'PrimaryIndel'=tmb_primary$DEL + tmb_primary$INS,
                        'BrainCNV'=tmb_brain$CNV,
                        'PrimaryCNV'=tmb_primary$CNV)

tmb_total[,c(2:7)] = tmb_total[,c(2:7)]/50

library(ggplot2)
library("gridExtra")
library("cowplot")

ylimits <- c(0,max(c(tmb_total$BrainSNP,tmb_total$PrimarySNP)))

Brastianos_tmb <- tmb_total[c(38:75),]
Brastianos_tmb$Tumor_Sample_Barcode <- factor(Brastianos_tmb$Tumor_Sample_Barcode, 
                          levels = Brastianos_tmb$Tumor_Sample_Barcode[order(Brastianos_tmb$BrainSNP,decreasing = TRUE)])

Brastianos_tmb_plot <- ggplot(Brastianos_tmb,aes(x=Tumor_Sample_Barcode)) + 
  geom_point(aes(y=BrainSNP,color = "BrM"),shape=2,size=2) + 
  geom_point(aes(y=PrimarySNP,,color = "Primary"),shape=4,size=2) + 
  labs(y="Muts per Mb") + 
  #geom_hline(yintercept=17,linetype="dashed") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y=element_blank(), 
        #axis.ticks.y=element_blank(),
        text=element_text(size=10),
        #legend.position = "none",
        plot.margin=unit(c(1,0,0,1), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.2, linetype = 'dashed',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  scale_colour_discrete("TMB") + 
  ylim(ylimits) + 
  ggtitle("Brastianos 2015 \n (n=38)") +
  theme(plot.title = element_text(size=10, hjust = 0.5))
Brastianos_tmb_plot

Own_tmb <- tmb_total[c(1:23),]
Own_tmb$Tumor_Sample_Barcode <- factor(Own_tmb$Tumor_Sample_Barcode, 
                                       levels = Own_tmb$Tumor_Sample_Barcode[order(Own_tmb$BrainSNP,decreasing = TRUE)])

Own_tmb_plot <- ggplot(Own_tmb,aes(x=Tumor_Sample_Barcode)) + 
  geom_point(aes(y=BrainSNP,color = "BrM"),shape=2,size=2) + 
  geom_point(aes(y=PrimarySNP,color = "Primary"),shape=4,size=2) + 
  labs(y="Muts per Mb") + 
  #geom_hline(yintercept=17,linetype="dashed") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        text=element_text(size=10),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin=unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  scale_colour_discrete("TMB") + 
  ylim(ylimits) + 
  ggtitle("SYSUCC \n (n=23)") +
  theme(plot.title = element_text(size=10, hjust = 0.5))
Own_tmb_plot

Fuku_tmb <- tmb_total[c(24:37),]
Fuku_tmb$Tumor_Sample_Barcode <- factor(Fuku_tmb$Tumor_Sample_Barcode, 
                                       levels = Fuku_tmb$Tumor_Sample_Barcode[order(Fuku_tmb$BrainSNP,decreasing = TRUE)])
Fuku_tmb_plot <- ggplot(Fuku_tmb,aes(x=Tumor_Sample_Barcode)) + 
  geom_point(aes(y=BrainSNP,color = "BrM"),shape=2,size=2) + 
  geom_point(aes(y=PrimarySNP,color = "Primary"),shape=4,size=2) + 
  labs(y="Muts per Mb") + 
  #geom_hline(yintercept=17,linetype="dashed") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        text=element_text(size=10),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin=unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  scale_colour_discrete("TMB") + 
  ylim(ylimits) + 
  ggtitle("Fukumara 2021 \n (n=14)") +
  theme(plot.title = element_text(size=10, hjust = 0.5))
Fuku_tmb_plot

Bras_tmb <- tmb_total[c(51:88),]
Bras_tmb$Tumor_Sample_Barcode <- factor(Bras_tmb$Tumor_Sample_Barcode, 
                                        levels = Bras_tmb$Tumor_Sample_Barcode[order(Bras_tmb$BTMB,decreasing = TRUE)])

Bras_tmb_plot <- ggplot(Bras_tmb,aes(x=Tumor_Sample_Barcode)) + 
  geom_point(aes(y=BTMB,color = "BrM"),shape=15,size=1) + 
  geom_point(aes(y=PTMB,,color = "Primary"),shape=10,size=1) + 
  labs(y="Muts per Mb") + 
  geom_hline(yintercept=17,linetype="dashed") +
  theme_linedraw() +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        text=element_text(size=6),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin=unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "lightyellow",
                                        colour = "lightyellow",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightyellow"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "lightyellow")
  ) +
  scale_colour_discrete("TMB") + 
  ylim(ylimits) + 
  annotate("text", x = 32, y = 100, label = "Hypermut.cutoff\n17.0 per Mb",size = 2) +
  ggtitle("Brastianos cohort \n (n=38)") +
  theme(plot.title = element_text(size=10, hjust = 0.5))

legend_tmb <- get_legend(Brastianos_tmb_plot)

mut_all = grid.arrange(Brastianos_tmb_plot,Own_tmb_plot,Fuku_tmb_plot,nrow=1,
                       widths = c(39,23,14))
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/Mutspermb.tif",
     width=2000,height=600,res=200)
grid.arrange(Brastianos_tmb_plot,Own_tmb_plot,Fuku_tmb_plot,nrow=1,
             widths = c(39,23,14))

dev.off()

tmb_plot <- plot_grid(Liu_tmb_plot,Fuku_tmb_plot,Bras_tmb_plot,nrow=1,
                         rel_widths = c(13,14,38))
# Indels per Mb
#DHP24 PRIMARYINDEL = 5 (58)
indel_ylimits <- c(0,max(c(tmb_total$BrainIndel,tmb_total$PrimaryIndel)))

Brastianos_indel <- tmb_total[c(38:75),]
Brastianos_indel$Tumor_Sample_Barcode <- factor(Brastianos_indel$Tumor_Sample_Barcode, 
                                              levels = Brastianos_indel$Tumor_Sample_Barcode[order(Brastianos_tmb$BrainSNP,decreasing = TRUE)])

Brastianos_indel_plot <- ggplot(Brastianos_indel,aes(x=Tumor_Sample_Barcode)) + 
  geom_point(aes(y=BrainIndel,color = "BrM"),shape=2,size=2) + 
  geom_point(aes(y=PrimaryIndel,,color = "Primary"),shape=4,size=2) + 
  labs(y="Indels per Mb") + 
  #geom_hline(yintercept=17,linetype="dashed") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(size=10),
        legend.position = "none",
        plot.margin=unit(c(1,0,0,1), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  scale_colour_discrete("TMB") + 
  ylim(indel_ylimits)
  #ggtitle("Brastianos 2015 \n (n=38)") +
  #theme(plot.title = element_text(size=10, hjust = 0.5))
Brastianos_indel_plot

Own_indel <- tmb_total[c(1:23),]
Own_indel$Tumor_Sample_Barcode <- factor(Own_indel$Tumor_Sample_Barcode, 
                                       levels = Own_indel$Tumor_Sample_Barcode[order(Own_tmb$BrainSNP,decreasing = TRUE)])

Own_indel_plot <- ggplot(Own_indel,aes(x=Tumor_Sample_Barcode)) + 
  geom_point(aes(y=BrainIndel,color = "BrM"),shape=2,size=2) + 
  geom_point(aes(y=PrimaryIndel,color = "Primary"),shape=4,size=2) + 
  labs(y="Muts per Mb") + 
  #geom_hline(yintercept=17,linetype="dashed") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        text=element_text(size=10),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin=unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  scale_colour_discrete("Indel") + 
  ylim(indel_ylimits) 
  #ggtitle("SYSUCC \n (n=23)") +
  #theme(plot.title = element_text(size=10, hjust = 0.5))
Own_indel_plot


Fuku_indel <- tmb_total[c(24:37),]
Fuku_indel$Tumor_Sample_Barcode <- factor(Fuku_indel$Tumor_Sample_Barcode, 
                                        levels = Fuku_indel$Tumor_Sample_Barcode[order(Fuku_tmb$BrainSNP,decreasing = TRUE)])
Fuku_indel_plot <- ggplot(Fuku_indel,aes(x=Tumor_Sample_Barcode)) + 
  geom_point(aes(y=BrainIndel,color = "BrM"),shape=2,size=2) + 
  geom_point(aes(y=PrimaryIndel,color = "Primary"),shape=4,size=2) + 
  labs(y="Muts per Mb") + 
  #geom_hline(yintercept=17,linetype="dashed") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        text=element_text(size=10),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin=unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  scale_colour_discrete("TMB") + 
  ylim(ylimits)
  #ggtitle("Fukumara 2021 \n (n=14)") +
  #theme(plot.title = element_text(size=10, hjust = 0.5))
Fuku_indel_plot

Fuku_tmb_indel <- plot_grid(Fuku_tmb_plot, Fuku_indel_plot,ncol=1,nrow=2,align="v",rel_heights = c(1,0.5))
legend_indel = get_legend(Bras_indel_plot)

indel_all = plot_grid(Brastianos_indel_plot,Own_indel_plot,Fuku_indel_plot,nrow=1,
                      rel_widths = c(39,23,14))

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/Indelspermb.tif",
     width=2000,height=400,res=200)
plot_grid(Brastianos_indel_plot,Own_indel_plot,Fuku_indel_plot,nrow=1,
          rel_widths = c(39,23,14))
dev.off()

####CNV
cnv_ylimits <- c(0,max(c(tmb_total$BrainCNV,tmb_total$PrimaryCNV)))

Brastianos_cnv <- tmb_total[c(38:75),]
Brastianos_cnv$Tumor_Sample_Barcode <- factor(Brastianos_cnv$Tumor_Sample_Barcode, 
                                                levels = Brastianos_cnv$Tumor_Sample_Barcode[order(Brastianos_tmb$BrainSNP,decreasing = TRUE)])

Brastianos_cnv_plot <- ggplot(Brastianos_cnv,aes(x=Tumor_Sample_Barcode)) + 
  geom_point(aes(y=BrainCNV,color = "BrM"),shape=2,size=2) + 
  geom_point(aes(y=PrimaryCNV,,color = "Primary"),shape=4,size=2) + 
  labs(y="CNVs per Mb") + 
  #geom_hline(yintercept=17,linetype="dashed") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(size=10),
        legend.position = "none",
        plot.margin=unit(c(1,0,0,1), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  scale_colour_discrete("CNV") + 
  ylim(cnv_ylimits)
#ggtitle("Brastianos 2015 \n (n=38)") +
#theme(plot.title = element_text(size=10, hjust = 0.5))
Brastianos_cnv_plot

Own_cnv <- tmb_total[c(1:23),]
Own_cnv$Tumor_Sample_Barcode <- factor(Own_cnv$Tumor_Sample_Barcode, 
                                         levels = Own_cnv$Tumor_Sample_Barcode[order(Own_tmb$BrainSNP,decreasing = TRUE)])

Own_cnv_plot <- ggplot(Own_cnv,aes(x=Tumor_Sample_Barcode)) + 
  geom_point(aes(y=BrainCNV,color = "BrM"),shape=2,size=2) + 
  geom_point(aes(y=PrimaryCNV,color = "Primary"),shape=4,size=2) + 
  labs(y="CNVs per Mb") + 
  #geom_hline(yintercept=17,linetype="dashed") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        text=element_text(size=10),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin=unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  scale_colour_discrete("CNV") + 
  ylim(cnv_ylimits) 
#ggtitle("SYSUCC \n (n=23)") +
#theme(plot.title = element_text(size=10, hjust = 0.5))
Own_cnv_plot


Fuku_cnv <- tmb_total[c(24:37),]
Fuku_cnv$Tumor_Sample_Barcode <- factor(Fuku_cnv$Tumor_Sample_Barcode, 
                                          levels = Fuku_cnv$Tumor_Sample_Barcode[order(Fuku_tmb$BrainSNP,decreasing = TRUE)])
Fuku_cnv_plot <- ggplot(Fuku_cnv,aes(x=Tumor_Sample_Barcode)) + 
  geom_point(aes(y=BrainCNV,color = "BrM"),shape=2,size=2) + 
  geom_point(aes(y=PrimaryCNV,color = "Primary"),shape=4,size=2) + 
  labs(y="CNVs per Mb") + 
  #geom_hline(yintercept=17,linetype="dashed") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        text=element_text(size=10),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.margin=unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")) +
  scale_colour_discrete("CNV") + 
  ylim(cnv_ylimits)
#ggtitle("Fukumara 2021 \n (n=14)") +
#theme(plot.title = element_text(size=10, hjust = 0.5))
Fuku_cnv_plot

cnv_all = plot_grid(Brastianos_cnv_plot,Own_cnv_plot,Fuku_cnv_plot,nrow=1,
                    rel_widths = c(39,23,14))

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/CNVspermb.tif",
     width=2000,height=400,res=200)
plot_grid(Brastianos_cnv_plot,Own_cnv_plot,Fuku_cnv_plot,nrow=1,
          rel_widths = c(39,23,14))
dev.off()

SYSUCC <- tmb_total[1:23,]
Fuku <- tmb_total[24:37,]
Bras <- tmb_total[38:75,]

SYSUCC_order <- arrange(SYSUCC,desc(BrainSNP))
Fuku_order <- arrange(Fuku,desc(BrainSNP))
Bras_order <- arrange(Bras,desc(BrainSNP))
ordersample = c(as.character(Bras_order$Tumor_Sample_Barcode),as.character(SYSUCC_order$Tumor_Sample_Barcode),
                as.character(Fuku_order$Tumor_Sample_Barcode))

######Coverage######
library(readxl)
library(reshape2)
coverage_brm = read_excel("~/Desktop/Lab/Brain_meta/total_coverage.xlsx",sheet = "brm")
coverage_lung = read_excel("~/Desktop/Lab/Brain_meta/total_coverage.xlsx",sheet = "lung")
coverage_normal = read_excel("~/Desktop/Lab/Brain_meta/total_coverage.xlsx",sheet = "normal")

coverage_total = data.frame("Sample" = coverage_brm$...1,
                            'BrM' = coverage_brm$MEAN_TARGET_COVERAGE,
                            'Primary' = coverage_lung$MEAN_TARGET_COVERAGE,
                            'Normal' = coverage_normal$MEAN_TARGET_COVERAGE)

brat_coverage = coverage_total[c(15:52),]
#brat_coverage$Sample = coverage_total$Sample[15:52]
#coverage_long = melt(coverage_total, id="Sample")
brat_coverage$Sample <- factor(brat_coverage$Sample,
                                levels = as.character(Bras_order$Tumor_Sample_Barcode))

brat_coverage_plot <- ggplot(brat_coverage, aes(Sample,group=1)) + 
  geom_point(aes(y = BrM, colour = "BrM"),shape=20,size=1) + 
  geom_point(aes(y = Primary, colour = "Primary"),shape=20,size=1) +
  geom_point(aes(y = Normal, colour = "Normal"),shape=20,size=1) +
  ylab("Coverage") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        text=element_text(size=10),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        #plot.margin=unit(c(-0.1,0,0,0), "cm"),
        plot.margin=unit(c(1,0,0,1), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")
  ) +
  ylim(c(0,400))
brat_coverage_plot

own_coverage = coverage_total[c(53:75),]
#brat_coverage$Sample = coverage_total$Sample[15:52]
#coverage_long = melt(coverage_total, id="Sample")
own_coverage$Sample <- factor(own_coverage$Sample,
                               levels = as.character(SYSUCC_order$Tumor_Sample_Barcode))

own_coverage_plot <- ggplot(own_coverage, aes(Sample,group=1)) + 
  geom_point(aes(y = BrM, colour = "BrM"),shape=20,size=1) + 
  geom_point(aes(y = Primary, colour = "Primary"),shape=20,size=1) +
  geom_point(aes(y = Normal, colour = "Normal"),shape=20,size=1) +
  ylab("Coverage") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=10),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        #plot.margin=unit(c(-0.1,0,0,0), "cm"),
        plot.margin=unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")
  ) +
  ylim(c(0,400))
own_coverage_plot

fuku_coverage = coverage_total[c(1:14),]
#brat_coverage$Sample = coverage_total$Sample[15:52]
#coverage_long = melt(coverage_total, id="Sample")
fuku_coverage$Sample <- factor(fuku_coverage$Sample,
                              levels = as.character(Fuku_order$Tumor_Sample_Barcode))

fuku_coverage_plot <- ggplot(fuku_coverage, aes(Sample,group=1)) + 
  geom_point(aes(y = BrM, colour = "BrM"),shape=20,size=1) + 
  geom_point(aes(y = Primary, colour = "Primary"),shape=20,size=1) +
  geom_point(aes(y = Normal, colour = "Normal"),shape=20,size=1) +
  ylab("Coverage") +
  theme_linedraw()+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=10),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        #plot.margin=unit(c(-0.1,0,0,0), "cm"),
        plot.margin=unit(c(1,0,0,0), "cm"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white")
  )+ylim(c(0,400))
fuku_coverage_plot

coverage_legend = get_legend(brat_coverage_plot)
coverage_all = plot_grid(brat_coverage_plot,own_coverage_plot,fuku_coverage_plot,nrow=1,
                         rel_widths = c(39,23,14))

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/coverage.tif",
     width=2000,height=400,res=200)
plot_grid(brat_coverage_plot,own_coverage_plot,fuku_coverage_plot,nrow=1,
          rel_widths = c(39,23,14))
dev.off()

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/clinical2.tif",width=2000,height=2000,res=300)
plot_grid(mut_all,indel_all,cnv_all,coverage_all,
          ncol=1,nrow=4,align="v",axis = "lr",
          rel_heights = c(3,2,2,2)
) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
dev.off()

#######boxplot########
library(ggpubr)
library(tidyverse)
library(rstatix)
tmb_total2 = melt(tmb_total)
tmb_total$cohort = c(rep("SYSUCC",23),rep("Fukumura 2021",14),rep("Brastianos 2015",38))
tmb_total$cohort = as.factor(tmb_total$cohort)
stat.test <- tmb_total2 %>%
  group_by(cohort) %>%
  t_test(value~variable) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")


bp <- ggbarplot(
  tmb_total2, x = "variable", y = "value", fill = "cohort",
  palette = "npg", add = "mean_sd",
  position = position_dodge(0.8)
) + stat_pvalue_manual(
  stat.test, x = "variable",  label = "p.adj.signif", 
  tip.length = 0.01, vjust = 2
)
bp

###brainSNP####
par(mar=c(5,6,4,2)+0.1,mgp=c(5,1,0))

df_snp = tmb_total2[c(1:150),]
df_snp$SNP = rep("SNP",150)

ggbarplot(df_snp, x = "SNP", y="value",
          color = "variable",
          add = c("mean_se", "point"),
          lab.size = 0.5,
          width = 0.5,
          palette = c("#1b98e0", "red"),
          #fill = "variable",
          position = position_dodge(),
          order = c("PrimarySNP","BrainSNP")) +
  stat_compare_means(aes(group = variable),
                     method = "wilcox.test",
                     paired=TRUE,
                     label.y = 60,
                     label = "p.signif") +
  coord_flip() +scale_y_continuous(position="right")

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/snp_sig.tif",width=2000,height=800,res=300)
p1 = ggplot(data = df_snp, aes(x=SNP,y=value,fill=variable)) +
  geom_boxplot(outlier.size = 1,position = position_dodge(1.0),width = 0.5) +
  stat_compare_means(aes(group = variable),
                     method = "wilcox.test",
                     paired=TRUE,
                     label.x = 0.1, 
                     label.y = 62,
                     label = "p.format") +
  scale_color_manual(values=c("red", "#1b98e0")) +
  theme_classic2() + 
  ylab("TMB") +
  theme(#axis.title.x=element_blank(), 
        #axis.text.x=element_blank(), 
        #axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        text=element_text(size=15),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 20)
  ) +
  coord_flip() +
  scale_y_continuous(position = "right",limits = c(0,70))

dev.off()

p1
df_indel = tmb_total2[c(151:300),]
df_indel$Indel = rep("Indel",150)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/indel_sig.tif",width=2000,height=700,res=300)
p2 = ggplot(data = df_indel, aes(x=Indel,y=value,fill=variable)) +
  geom_boxplot(outlier.size = 1,position = position_dodge(1.0),width = 0.5) +
  stat_compare_means(aes(group = variable),
                     method = "wilcox.test",
                     paired=TRUE,
                     label.x = 0.1, 
                     label.y = 4.3,
                     label = "p.format"
                     ) +
  scale_color_manual(values=c("red", "#1b98e0")) +
  theme_classic2() + 
  ylab("TMB") +
  theme(axis.title.x=element_blank(), 
    #axis.text.x=element_blank(), 
    #axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    text=element_text(size=15),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 20,  # Top margin
                         r = 20,  # Right margin
                         b = 10,  # Bottom margin
                         l = 20)
  ) +
  coord_flip() +
  scale_y_continuous(position = "right",limits = c(0,5))

dev.off()

df_cnv = tmb_total2[c(301:450),]
df_cnv$CNV = rep("CNV",150)
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/cnv_sig.tif",width=2000,height=800,res=300)
p3 = ggplot(data = df_cnv, aes(x=CNV,y=value,fill=variable)) +
  geom_boxplot(outlier.size = 1,position = position_dodge(1.0),width = 0.5) +
  stat_compare_means(aes(group = variable),
                     method = "wilcox.test",
                     paired=TRUE,
                     label.x = 0.1, 
                     label.y = 70,
                     label = "p.format") +
  scale_color_manual(values=c("red", "#1b98e0")) +
  theme_classic2() + 
  ylab("TMB") +
  theme(axis.title.x=element_blank(), 
        #axis.text.x=element_blank(), 
        #axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        text=element_text(size=15),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 10,  # Bottom margin
                             l = 20)
  ) +
  coord_flip() +
  scale_y_continuous(position = "right",limits = c(0,80))

dev.off()

p_leg = ggplot(data = df_snp, aes(x=SNP,y=value,fill=variable)) +
  geom_boxplot(outlier.size = 1,position = position_dodge(1.0),width = 0.5) +
  stat_compare_means(aes(group = variable),
                     method = "wilcox.test",
                     paired=TRUE,
                     label.x = 0.1, 
                     label.y = 62,
                     label = "p.format") +
  scale_color_manual(values=c("red", "#1b98e0")) +
  scale_fill_discrete(labels = c("BrM","Primary")) +
  theme_classic2() + 
  ylab("TMB") +
  theme(#axis.title.x=element_blank(), 
    #axis.text.x=element_blank(), 
    #axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    text=element_text(size=15),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = margin(t = 20,  # Top margin
                         r = 20,  # Right margin
                         b = 10,  # Bottom margin
                         l = 20)
  ) +
  coord_flip() +
  scale_y_continuous(position = "right",limits = c(0,70))

p_leg
leg1 <- get_legend(p_leg)
# create a blank plot for legend alignment 
blank_p <- plot_spacer() + theme_void()
# combine legend 3 & blank plot
leg10 <- plot_grid(leg1, blank_p,
                   nrow = 2
)

leg12 = plot_grid(legend_tmb,coverage_legend,blank_p,nrow = 3)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-2legend.tif",
     width=500,height=500,res=200)
leg12
dev.off()

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/coverage.tif",
     width=2000,height=400,res=200)
plot_grid(brat_coverage_plot,own_coverage_plot,fuku_coverage_plot,nrow=1,
          rel_widths = c(39,23,14))
dev.off()


final_p <- plot_grid(p1,
                     p2, 
                     p3,
                     leg10,
                     nrow = 4,
                     ncol = 1,
                     align = "v",
                     #axis = "t",
                     rel_heights = c(1,0.9,0.8,0.3)
)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-2draft.tif",
     width=1000,height=1500,res=200)
print(final_p)
dev.off()
