library(stringr)

primary_cyto = primary.gistic@cytoband.summary
brm_cyto = brain.gistic@cytoband.summary

brm.scores.gis <- read.table("/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalbrain_results/scores.gistic",
                         sep = "\t")
prim.scores.gis <- read.table("/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalprimary_results/scores.gistic",
                          sep = "\t")
colnames(prim.scores.gis) = prim.scores.gis[1,]
prim.scores.gis = prim.scores.gis[-1,]

colnames(brm.scores.gis) = brm.scores.gis[1,]
brm.scores.gis = brm.scores.gis[-1,]

brm.scores.gis$cytoband = NA

brm_freq = data.frame(group = "BrM", value = brm.scores.gis$frequency)
prim_freq = data.frame(group = "Primary", value = prim.scores.gis$frequency)
df<-data.frame(X=c(brm.scores.gis$`G-score`,prim.scores.gis$`G-score`),
               Grp=rep(c("BrM","Primary"),times=c(28309,29961)))
df$X = as.numeric(df$X)
df$Grp <- factor(df$Grp,     # Reorder factor levels
                         c("Primary", "BrM"))
library(ggplot2)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-S4C.tif",
     width=600,height=1000,res=200)
ggplot(data=df, aes(x=Grp, y=X, color = Grp)) + geom_boxplot(outlier.shape = NA) +
  #stat_compare_means(method = "wilcox.test") +
  annotate("text", x = 1.5, y = 0.28, label="p = 0.18", size = 5) +
  ylab("GISTIC Score") +
  scale_color_manual(values=c("#1b98e0","red")) + 
  theme_bw() +
  theme(axis.title.x=element_blank(), 
        #axis.text.x=element_blank(), 
        #axis.ticks.x=element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 20,  # Bottom margin
                             l = 20),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,0.3))

dev.off()
for (k in 1:28309){
  chr1 = brm.scores.gis$Chromosome[k]
  start1 = strtoi(brm.scores.gis$Start[k])
  end1 = strtoi(brm.scores.gis$End[k])
  for (i in 1:154){
    stri = brm_cyto$Wide_Peak_Limits[i]
    chr = str_remove(str_split(stri,":")[[1]][1],"chr")
    start = strtoi(str_split(str_split(stri,":")[[1]][2],"-")[[1]][1])
    end = strtoi(str_split(str_split(stri,":")[[1]][2],"-")[[1]][2])
    if ((chr == chr1) & (start < start1) & (end > end1)){
      brm.scores.gis$cytoband[k] = brm_cyto$Cytoband[i]
    }
  }
}

prim.scores.gis$cytoband = NA
for (k in 1:29961){
  chr1 = prim.scores.gis$Chromosome[k]
  start1 = strtoi(prim.scores.gis$Start[k])
  end1 = strtoi(prim.scores.gis$End[k])
  for (i in 1:178){
    stri = primary_cyto$Wide_Peak_Limits[i]
    chr = str_remove(str_split(stri,":")[[1]][1],"chr")
    start = strtoi(str_split(str_split(stri,":")[[1]][2],"-")[[1]][1])
    end = strtoi(str_split(str_split(stri,":")[[1]][2],"-")[[1]][2])
    if ((chr == chr1) & (start < start1) & (end > end1)){
      prim.scores.gis$cytoband[k] = primary_cyto$Cytoband[i]
    }
  }
}

brm.scores.gis <- na.omit(brm.scores.gis) 
prim.scores.gis <- na.omit(prim.scores.gis)

brm.scores.gis$diff = NA

for (i in 1:3253){
  type = brm.scores.gis$Type[i]
  chr = brm.scores.gis$Chromosome[i]
  start = strtoi(brm.scores.gis$Start[i])
  end = strtoi(brm.scores.gis$End[i])
  for (k in 1:7330){
    type1 = prim.scores.gis$Type[i]
    chr1 = prim.scores.gis$Chromosome[i]
    start1 = strtoi(prim.scores.gis$Start[i])
    end1 = strtoi(prim.scores.gis$End[i])
    if ((type == type1) & (chr == chr1) & ((abs(start - start1)) < 1000000)){
      brm.scores.gis$diff[i] = abs(as.numeric(brm.scores.gis$`G-score`[i]) - as.numeric(prim.scores.gis$`G-score`[k]))
    }
  }
}
brm.scores.gis2 <- aggregate(list(brm.scores.gis$`G-score`,brm.scores.gis$`-log10(q-value)`), 
                             by = list(brm.scores.gis$Type, brm.scores.gis$cytoband), max)
prim.scores.gis2 <- aggregate(list(prim.scores.gis$`G-score`,prim.scores.gis$`-log10(q-value)`), 
                             by = list(prim.scores.gis$Type, prim.scores.gis$cytoband), max)
colnames(brm.scores.gis2) = c("Type","Cytoband","G-score","-log10(q-value)")
colnames(prim.scores.gis2) = c("Type","Cytoband","G-score","-log10(q-value)")
brm.scores.gis2$diff = NA

for (i in 1:114){
  type = brm.scores.gis2$Type[i]
  cyto = brm.scores.gis2$Cytoband[i]
  for (k in 1:99){
    type1 = prim.scores.gis2$Type[k]
    cyto1 = prim.scores.gis2$Cytoband[k]
    if ((type == type1) && (cyto ==cyto1)){
      brm.scores.gis2$diff[i] = abs(as.numeric(brm.scores.gis2$`G-score`[i]) - as.numeric(prim.scores.gis2$`G-score`[k]))
    }
  }
}
brm.scores.gis2$`-log10(q-value)` = as.numeric(brm.scores.gis2$`-log10(q-value)`)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-4B.tif",
     width=1200,height=1200,res=200)
ggplot(data = brm.scores.gis2, aes(x = diff, y = `-log10(q-value)`,color = Type)) +
  geom_point() +
  geom_label(
    data=brm.scores.gis2 %>% filter(diff>0.2 & `-log10(q-value)`>2), # Filter data first
    aes(label=Cytoband), show.legend = FALSE, 
    nudge_x = 0.08, nudge_y = 0.05,check_overlap = T) +
  scale_color_manual(labels = c("Amplication", "Deletion"),
                     values = c("Amp" = "red", "Del" = "#1b98e0")) +
  xlab("GISTIC score difference") + 
  geom_hline(yintercept=2,linetype=2) + 
  geom_vline(xintercept=0.2,linetype=2) + 
  theme_linedraw() +
  theme(#axis.title.x=element_blank(), 
        #axis.text.x=element_blank(), 
        #axis.ticks.x=element_blank(),
        #axis.title.y=element_blank(),
        #\axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        legend.position = c(0.15,0.9),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 20,  # Bottom margin
                             l = 20),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )
dev.off()

brm_chr3q25 = chr3.1[which(chr3.1$midpoint > 1 & chr3.1$midpoint < 18022430),]
brm_chr3q25$mid = brm_chr3q25$midpoint / 1e+06
prim_chr3q25 = primchr3.1[which(primchr3.1$midpoint > 1 & primchr3.1$midpoint < 18022430),]
prim_chr3q25$mid = prim_chr3q25$midpoint / 1e+06

p3q25 = ggplot(data = brm_chr3q25, aes(x=mid,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = prim_chr3q25, aes(x=mid,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  annotate("text", x = 15, y = 0.18, label="3q25", color="#1b98e0", size = 5) + 
  annotate("text", x = 15, y = 0.17, label="q = 0.2", size = 5) +
  #annotate("text", x = 10, y = 0.32, label="CD69, CLEC2D, CLECL1", size = 5,fontface = "italic") +
  theme_linedraw() + 
  ylab("GISTIC score") +
  xlab("Genomic position (Mb)") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(#axis.title.x=element_blank(), 
    #axis.text.x=element_blank(), 
    #axis.ticks.x=element_blank(),
    #axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    text=element_text(size=12),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = margin(t = 20,  # Top margin
                         r = 20,  # Right margin
                         b = 20,  # Bottom margin
                         l = 20),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
    #panel.border = element_blank()
  )
#+ ylim(c(0,1.2))
p3q25

chr11.1$Start = as.integer(chr11.1$Start)
primchr11.1$Start = as.integer(primchr11.1$Start)
brm_chr11q13 = chr11.1[which(chr11.1$Start > 7.5e+07 & chr11.1$End < 5e+08),]
brm_chr11q13$mid = brm_chr11q13$midpoint / 1e+06
prim_chr11q13 = primchr11.1[which(primchr11.1$Start > 7.5e+07 & primchr11.1$End < 5e+08),]
prim_chr11q13$mid = prim_chr11q13$midpoint / 1e+06

p11q13 = ggplot(data = brm_chr11q13, aes(x=mid,y=`G-score`,group=1)) + 
  geom_line(color = "red",size = 2, alpha = 0.6) +
  geom_line(data = prim_chr11q13, aes(x=mid,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  annotate("text", x = 128, y = 0.6, label="11q13", color="red", size = 5) + 
  annotate("text", x = 128, y = 0.55, label="q = 0.00015", size = 5) +
  annotate("text", x = 110, y = 0.45, label="SCYL1", size = 5,fontface = "italic") +
  theme_linedraw() + 
  ylab("GISTIC score") +
  xlab("Genomic position (Mb)") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(#axis.title.x=element_blank(), 
        #axis.text.x=element_blank(), 
        #axis.ticks.x=element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 20,  # Bottom margin
                             l = 20),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )
#+ ylim(c(0,1.2))
p11q13

chr14.1$Start = as.integer(chr14.1$Start)
primchr14.1$Start = as.integer(primchr14.1$Start)
brm_chr14q11 = chr14.1[which(chr14.1$midpoint > 12591978 & chr14.1$midpoint < 32938171),]
brm_chr14q11$mid = brm_chr14q11$midpoint / 1e+06
prim_chr14q11 = primchr14.1[which(primchr14.1$midpoint > 12591978 & primchr14.1$midpoint < 32938171),]
prim_chr14q11$mid = prim_chr14q11$midpoint / 1e+06

p14q11 = ggplot(data = brm_chr14q11, aes(x=mid,y=`G-score`,group=1)) + 
  geom_line(color = "red",size = 2, alpha = 0.6) +
  geom_line(data = prim_chr14q11, aes(x=mid,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  annotate("text", x = 30, y = 0.7, label="14q11", color="red", size = 5) + 
  annotate("text", x = 30, y = 0.65, label="q = 0.00025", size = 5) +
  annotate("text", x = 25, y = 0.55, label="DAD1", size = 5,fontface = "italic") +
  theme_linedraw() + 
  ylab("GISTIC score") +
  xlab("Genomic position (Mb)") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(#axis.title.x=element_blank(), 
    #axis.text.x=element_blank(), 
    #axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    text=element_text(size=12),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = margin(t = 20,  # Top margin
                         r = 20,  # Right margin
                         b = 20,  # Bottom margin
                         l = 20),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
    #panel.border = element_blank()
  )
#+ ylim(c(0,1.2))
p14q11

brm_chr12p13 = chr12.1[which(chr12.1$midpoint > 8760437 & chr12.1$midpoint < 10980181),]
brm_chr12p13$mid = brm_chr12p13$midpoint / 1e+06
prim_chr12p13 = primchr12.1[which(primchr12.1$midpoint > 8760437 & primchr12.1$midpoint < 10980181),]
prim_chr12p13$mid = prim_chr12p13$midpoint / 1e+06

p12p13 = ggplot(data = brm_chr12p13, aes(x=mid,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = prim_chr12p13, aes(x=mid,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  annotate("text", x = 10.6, y = 0.37, label="12p13", color="#1b98e0", size = 5) + 
  annotate("text", x = 10.6, y = 0.35, label="q = 0.0012", size = 5) +
  annotate("text", x = 10, y = 0.32, label="CD69, CLEC2D, CLECL1", size = 5,fontface = "italic") +
  theme_linedraw() + 
  ylab("GISTIC score") +
  xlab("Genomic position (Mb)") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(#axis.title.x=element_blank(), 
    #axis.text.x=element_blank(), 
    #axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    text=element_text(size=12),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_blank(),
    #legend.position = "none",
    plot.margin = margin(t = 20,  # Top margin
                         r = 20,  # Right margin
                         b = 20,  # Bottom margin
                         l = 20),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
    #panel.border = element_blank()
  )
#+ ylim(c(0,1.2))
p12p13

final_p <- plot_grid(p3q25,p11q13,p14q11,p12p13,
                     nrow = 1,
                     ncol = 4,
                     align = "h",
                     #axis = "t",
                     rel_widths = c(1,1,1,0.9)
)
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-4C.tif",
     width=5500,height=1500,res=300)
print(final_p)
dev.off()
