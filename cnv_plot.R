####cnv analysis####
scores.gis <- read.table("/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalbrain_results/scores.gistic",
                         sep = "\t")
colnames(scores.gis) = scores.gis[1,]
scores.gis = scores.gis[-1,]
chr1 = scores.gis[which(scores.gis$Type == "Del"),]
chr1.1 = chr1[which(chr1$Chromosome == " 1"),]
chr1.1$midpoint = (as.integer(chr1.1$Start) + as.integer(chr1.1$End)) /2
chr1.1$`G-score` = as.numeric(chr1.1$`G-score`)

chr2 = scores.gis[which(scores.gis$Type == "Del"),]
chr2.1 = chr2[which(chr2$Chromosome == " 2"),]
chr2.1$midpoint = (as.integer(chr2.1$Start) + as.integer(chr2.1$End)) /2
chr2.1$`G-score` = as.numeric(chr2.1$`G-score`)

chr3 = scores.gis[which(scores.gis$Type == "Del"),]
chr3.1 = chr3[which(chr3$Chromosome == " 3"),]
chr3.1$midpoint = (as.integer(chr3.1$Start) + as.integer(chr3.1$End)) /2
chr3.1$`G-score` = as.numeric(chr3.1$`G-score`)

chr4 = scores.gis[which(scores.gis$Type == "Del"),]
chr4.1 = chr4[which(chr4$Chromosome == " 4"),]
chr4.1$midpoint = (as.integer(chr4.1$Start) + as.integer(chr4.1$End)) /2
chr4.1$`G-score` = as.numeric(chr4.1$`G-score`)

chr5 = scores.gis[which(scores.gis$Type == "Del"),]
chr5.1 = chr5[which(chr5$Chromosome == " 5"),]
chr5.1$midpoint = (as.integer(chr5.1$Start) + as.integer(chr5.1$End)) /2
chr5.1$`G-score` = as.numeric(chr5.1$`G-score`)

chr6 = scores.gis[which(scores.gis$Type == "Del"),]
chr6.1 = chr6[which(chr6$Chromosome == " 6"),]
chr6.1$midpoint = (as.integer(chr6.1$Start) + as.integer(chr6.1$End)) /2
chr6.1$`G-score` = as.numeric(chr6.1$`G-score`)

chr7 = scores.gis[which(scores.gis$Type == "Del"),]
chr7.1 = chr7[which(chr7$Chromosome == " 7"),]
chr7.1$midpoint = (as.integer(chr7.1$Start) + as.integer(chr7.1$End)) /2
chr7.1$`G-score` = as.numeric(chr7.1$`G-score`)

chr8 = scores.gis[which(scores.gis$Type == "Del"),]
chr8.1 = chr8[which(chr8$Chromosome == " 8"),]
chr8.1$midpoint = (as.integer(chr8.1$Start) + as.integer(chr8.1$End)) /2
chr8.1$`G-score` = as.numeric(chr8.1$`G-score`)

chr9 = scores.gis[which(scores.gis$Type == "Del"),]
chr9.1 = chr9[which(chr9$Chromosome == " 9"),]
chr9.1$midpoint = (as.integer(chr9.1$Start) + as.integer(chr9.1$End)) /2
chr9.1$`G-score` = as.numeric(chr9.1$`G-score`)

chr9 = scores.gis[which(scores.gis$Type == "Del"),]
chr9.1 = chr9[which(chr9$Chromosome == " 9"),]
chr9.1$midpoint = (as.integer(chr9.1$Start) + as.integer(chr9.1$End)) /2
chr9.1$`G-score` = as.numeric(chr9.1$`G-score`)

chr10 = scores.gis[which(scores.gis$Type == "Del"),]
chr10.1 = chr10[which(chr10$Chromosome == "10"),]
chr10.1$midpoint = (as.integer(chr10.1$Start) + as.integer(chr10.1$End)) /2
chr10.1$`G-score` = as.numeric(chr10.1$`G-score`)

chr11 = scores.gis[which(scores.gis$Type == "Amp"),]
chr11.1 = chr11[which(chr11$Chromosome == "11"),]
chr11.1$midpoint = (as.integer(chr11.1$Start) + as.integer(chr11.1$End)) /2
chr11.1$`G-score` = as.numeric(chr11.1$`G-score`)

chr12 = scores.gis[which(scores.gis$Type == "Del"),]
chr12.1 = chr12[which(chr12$Chromosome == "12"),]
chr12.1$midpoint = (as.integer(chr12.1$Start) + as.integer(chr12.1$End)) /2
chr12.1$`G-score` = as.numeric(chr12.1$`G-score`)

chr13 = scores.gis[which(scores.gis$Type == "Del"),]
chr13.1 = chr13[which(chr13$Chromosome == "13"),]
chr13.1$midpoint = (as.integer(chr13.1$Start) + as.integer(chr13.1$End)) /2
chr13.1$`G-score` = as.numeric(chr13.1$`G-score`)

chr14 = scores.gis[which(scores.gis$Type == "Amp"),]
chr14.1 = chr14[which(chr14$Chromosome == "14"),]
chr14.1$midpoint = (as.integer(chr14.1$Start) + as.integer(chr14.1$End)) /2
chr14.1$`G-score` = as.numeric(chr14.1$`G-score`)

chr15 = scores.gis[which(scores.gis$Type == "Del"),]
chr15.1 = chr15[which(chr15$Chromosome == "15"),]
chr15.1$midpoint = (as.integer(chr15.1$Start) + as.integer(chr15.1$End)) /2
chr15.1$`G-score` = as.numeric(chr15.1$`G-score`)

chr16 = scores.gis[which(scores.gis$Type == "Del"),]
chr16.1 = chr16[which(chr16$Chromosome == "16"),]
chr16.1$midpoint = (as.integer(chr16.1$Start) + as.integer(chr16.1$End)) /2
chr16.1$`G-score` = as.numeric(chr16.1$`G-score`)

chr17 = scores.gis[which(scores.gis$Type == "Del"),]
chr17.1 = chr17[which(chr17$Chromosome == "17"),]
chr17.1$midpoint = (as.integer(chr17.1$Start) + as.integer(chr17.1$End)) /2
chr17.1$`G-score` = as.numeric(chr17.1$`G-score`)

chr18 = scores.gis[which(scores.gis$Type == "Del"),]
chr18.1 = chr18[which(chr18$Chromosome == "18"),]
chr18.1$midpoint = (as.integer(chr18.1$Start) + as.integer(chr18.1$End)) /2
chr18.1$`G-score` = as.numeric(chr18.1$`G-score`)

chr19 = scores.gis[which(scores.gis$Type == "Del"),]
chr19.1 = chr19[which(chr19$Chromosome == "19"),]
chr19.1$midpoint = (as.integer(chr19.1$Start) + as.integer(chr19.1$End)) /2
chr19.1$`G-score` = as.numeric(chr19.1$`G-score`)

chr20 = scores.gis[which(scores.gis$Type == "Del"),]
chr20.1 = chr20[which(chr20$Chromosome == "20"),]
chr20.1$midpoint = (as.integer(chr20.1$Start) + as.integer(chr20.1$End)) /2
chr20.1$`G-score` = as.numeric(chr20.1$`G-score`)

chr21 = scores.gis[which(scores.gis$Type == "Del"),]
chr21.1 = chr21[which(chr21$Chromosome == "21"),]
chr21.1$midpoint = (as.integer(chr21.1$Start) + as.integer(chr21.1$End)) /2
chr21.1$`G-score` = as.numeric(chr21.1$`G-score`)

chr22 = scores.gis[which(scores.gis$Type == "Del"),]
chr22.1 = chr22[which(chr22$Chromosome == "22"),]
chr22.1$midpoint = (as.integer(chr22.1$Start) + as.integer(chr22.1$End)) /2
chr22.1$`G-score` = as.numeric(chr22.1$`G-score`)

scores.gis2 <- read.table("/Users/Jianlan/Desktop/Lab/Brain_meta/Lung/cnv/totalprimary_results/scores.gistic",
                         sep = "\t")
colnames(scores.gis2) = scores.gis2[1,]
scores.gis2 = scores.gis2[-1,]
primchr1 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr1.1 = primchr1[which(primchr1$Chromosome == " 1"),]
primchr1.1$midpoint = (as.integer(primchr1.1$Start) + as.integer(primchr1.1$End)) /2
primchr1.1$`G-score` = as.numeric(primchr1.1$`G-score`)

primchr2 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr2.1 = primchr2[which(primchr2$Chromosome == " 2"),]
primchr2.1$midpoint = (as.integer(primchr2.1$Start) + as.integer(primchr2.1$End)) /2
primchr2.1$`G-score` = as.numeric(primchr2.1$`G-score`)

primchr3 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr3.1 = primchr3[which(primchr3$Chromosome == " 3"),]
primchr3.1$midpoint = (as.integer(primchr3.1$Start) + as.integer(primchr3.1$End)) /2
primchr3.1$`G-score` = as.numeric(primchr3.1$`G-score`)

primchr4 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr4.1 = primchr4[which(primchr4$Chromosome == " 4"),]
primchr4.1$midpoint = (as.integer(primchr4.1$Start) + as.integer(primchr4.1$End)) /2
primchr4.1$`G-score` = as.numeric(primchr4.1$`G-score`)

primchr5 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr5.1 = primchr5[which(primchr5$Chromosome == " 5"),]
primchr5.1$midpoint = (as.integer(primchr5.1$Start) + as.integer(primchr5.1$End)) /2
primchr5.1$`G-score` = as.numeric(primchr5.1$`G-score`)

primchr6 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr6.1 = primchr6[which(primchr6$Chromosome == " 6"),]
primchr6.1$midpoint = (as.integer(primchr6.1$Start) + as.integer(primchr6.1$End)) /2
primchr6.1$`G-score` = as.numeric(primchr6.1$`G-score`)

primchr7 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr7.1 = primchr7[which(primchr7$Chromosome == " 7"),]
primchr7.1$midpoint = (as.integer(primchr7.1$Start) + as.integer(primchr7.1$End)) /2
primchr7.1$`G-score` = as.numeric(primchr7.1$`G-score`)

primchr8 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr8.1 = primchr8[which(primchr8$Chromosome == " 8"),]
primchr8.1$midpoint = (as.integer(primchr8.1$Start) + as.integer(primchr8.1$End)) /2
primchr8.1$`G-score` = as.numeric(primchr8.1$`G-score`)

primchr9 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr9.1 = primchr9[which(primchr9$Chromosome == " 9"),]
primchr9.1$midpoint = (as.integer(primchr9.1$Start) + as.integer(primchr9.1$End)) /2
primchr9.1$`G-score` = as.numeric(primchr9.1$`G-score`)

primchr10 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr10.1 = primchr10[which(primchr10$Chromosome == "10"),]
primchr10.1$midpoint = (as.integer(primchr10.1$Start) + as.integer(primchr10.1$End)) /2
primchr10.1$`G-score` = as.numeric(primchr10.1$`G-score`)

primchr11 = scores.gis2[which(scores.gis2$Type == "Amp"),]
primchr11.1 = primchr11[which(primchr11$Chromosome == "11"),]
primchr11.1$midpoint = (as.integer(primchr11.1$Start) + as.integer(primchr11.1$End)) /2
primchr11.1$`G-score` = as.numeric(primchr11.1$`G-score`)

primchr12 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr12.1 = primchr12[which(primchr12$Chromosome == "12"),]
primchr12.1$midpoint = (as.integer(primchr12.1$Start) + as.integer(primchr12.1$End)) /2
primchr12.1$`G-score` = as.numeric(primchr12.1$`G-score`)

primchr13 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr13.1 = primchr13[which(primchr13$Chromosome == "13"),]
primchr13.1$midpoint = (as.integer(primchr13.1$Start) + as.integer(primchr13.1$End)) /2
primchr13.1$`G-score` = as.numeric(primchr13.1$`G-score`)

primchr14 = scores.gis2[which(scores.gis2$Type == "Amp"),]
primchr14.1 = primchr14[which(primchr14$Chromosome == "14"),]
primchr14.1$midpoint = (as.integer(primchr14.1$Start) + as.integer(primchr14.1$End)) /2
primchr14.1$`G-score` = as.numeric(primchr14.1$`G-score`)

primchr15 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr15.1 = primchr15[which(primchr15$Chromosome == "15"),]
primchr15.1$midpoint = (as.integer(primchr15.1$Start) + as.integer(primchr15.1$End)) /2
primchr15.1$`G-score` = as.numeric(primchr15.1$`G-score`)

primchr16 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr16.1 = primchr16[which(primchr16$Chromosome == "16"),]
primchr16.1$midpoint = (as.integer(primchr16.1$Start) + as.integer(primchr16.1$End)) /2
primchr16.1$`G-score` = as.numeric(primchr16.1$`G-score`)

primchr17 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr17.1 = primchr17[which(primchr17$Chromosome == "17"),]
primchr17.1$midpoint = (as.integer(primchr17.1$Start) + as.integer(primchr17.1$End)) /2
primchr17.1$`G-score` = as.numeric(primchr17.1$`G-score`)

primchr18 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr18.1 = primchr18[which(primchr18$Chromosome == "18"),]
primchr18.1$midpoint = (as.integer(primchr18.1$Start) + as.integer(primchr18.1$End)) /2
primchr18.1$`G-score` = as.numeric(primchr18.1$`G-score`)

primchr19 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr19.1 = primchr19[which(primchr19$Chromosome == "19"),]
primchr19.1$midpoint = (as.integer(primchr19.1$Start) + as.integer(primchr19.1$End)) /2
primchr19.1$`G-score` = as.numeric(primchr19.1$`G-score`)

primchr20 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr20.1 = primchr20[which(primchr20$Chromosome == "20"),]
primchr20.1$midpoint = (as.integer(primchr20.1$Start) + as.integer(primchr20.1$End)) /2
primchr20.1$`G-score` = as.numeric(primchr20.1$`G-score`)

primchr21 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr21.1 = primchr21[which(primchr21$Chromosome == "21"),]
primchr21.1$midpoint = (as.integer(primchr21.1$Start) + as.integer(primchr21.1$End)) /2
primchr21.1$`G-score` = as.numeric(primchr21.1$`G-score`)

primchr22 = scores.gis2[which(scores.gis2$Type == "Del"),]
primchr22.1 = primchr22[which(primchr22$Chromosome == "22"),]
primchr22.1$midpoint = (as.integer(primchr22.1$Start) + as.integer(primchr22.1$End)) /2
primchr22.1$`G-score` = as.numeric(primchr22.1$`G-score`)

p1 = ggplot(data = chr1.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr1.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("1") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
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
                         r = 0,  # Right margin
                         b = 20,  # Bottom margin
                         l = 20),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
    #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p2 = ggplot(data = chr2.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr2.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("2") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p3 = ggplot(data = chr3.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr3.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("3") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p4 = ggplot(data = chr4.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr4.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("4") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p5 = ggplot(data = chr5.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr5.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("5") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p6 = ggplot(data = chr6.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr6.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("6") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p7 = ggplot(data = chr7.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr7.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("7") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p8 = ggplot(data = chr8.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr8.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("8") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p9 = ggplot(data = chr9.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr9.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("9") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p10 = ggplot(data = chr10.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr10.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("10") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p11 = ggplot(data = chr11.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr11.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("11") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p12 = ggplot(data = chr12.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr12.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("12") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p13 = ggplot(data = chr13.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr13.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("13") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p14 = ggplot(data = chr14.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr14.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("14") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p15 = ggplot(data = chr15.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr15.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("15") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p16 = ggplot(data = chr16.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr16.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("16") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p17 = ggplot(data = chr17.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr17.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("17") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p18 = ggplot(data = chr18.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr18.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("18") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p19 = ggplot(data = chr19.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr19.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("19") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p20 = ggplot(data = chr20.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr20.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("20") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p21 = ggplot(data = chr21.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr21.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("21") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  )+ ylim(c(0,1.2))

p22 = ggplot(data = chr22.1, aes(x=midpoint,y=`G-score`,group=1)) + 
  geom_line(color = "#1b98e0",size = 2, alpha = 0.6) +
  geom_line(data = primchr22.1, aes(x=midpoint,y=`G-score`,group=1), 
            color = "black",size = 2,alpha = 0.6) +
  theme_linedraw() + 
  ggtitle("22") + 
  ylab("GISTIC score") +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=12),
        legend.key.size = unit(0.5, "cm"),
        legend.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
        #panel.border = element_blank()
  ) + ylim(c(0,1.2))

final_p <- plot_grid(p1,p2,p3,
                     p4,p5,p6, 
                     p7,p8,p9,
                     p10,p11,p12,
                     p13,p14,p15,
                     p16,p17,p18,
                     p19,p20,p21,
                     p22,
                     nrow = 1,
                     ncol = 22,
                     align = "h",
                     #axis = "t",
                     rel_widths = c(52.67,34.83,47.67,
                                     18,28.33,39,
                                     20.17,24.5,44.5,
                                     29.5,29.33,19.17,
                                     32.33,36,20.83,26.5,
                                     27,18.17,23.17,
                                     19.17,24.17,24.5)
)
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-4A_del_draft.tif",
     width=5000,height=1000,res=300)
print(final_p)
dev.off()

dat_leg = data.frame(chr1 = chr1$Start[1:100], BrM = chr1$`G-score`[1:100],Primary = primchr1$`G-score`[1:100])
dat_leg = melt(dat_leg, id.vars = "chr1")
dat_leg$variable <- factor(dat_leg$variable, levels = c("Primary", "BrM"))


p_leg = ggplot(dat_leg, aes(x=chr1, y = value, group = variable, colour = variable)) + 
  geom_line(alpha = 0.6,size = 2) +
  scale_color_manual(values = c("black","#1b98e0")) +
  theme(legend.position="bottom",
        legend.title = element_blank())
p_leg
leg1 <- get_legend(p_leg)
# create a blank plot for legend alignment 
blank_p <- plot_spacer() + theme_void()
# combine legend 3 & blank plot
leg10 <- plot_grid(leg1, blank_p,
                   nrow = 2
)

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/figure-4A_del_legend.tif",
     width=500,height=500,res=200)
leg10
dev.off()
