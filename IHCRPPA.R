library(reshape2)
library(ggpubr)

data1=NULL
data1=data.frame("condition"=c(rep("Primary",75),rep("BrM",75)),
                 "TMB"=c(tmb_total$PTMB[1:75],tmb_total$BTMB[1:75]))
tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/tmbcompare.tiff",
     width=1000,height=1000,res=200)
ggpaired(data1, x = "condition", y = "TMB",
         color = "condition", line.color = "gray", line.size = 0.2,
         palette = c("#1b98e0", "red"))+
  stat_compare_means(paired = TRUE,label.x = 1.2,label.y = 80) +
  labs(y="TMB \n (Muts per Mb)") + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = "none")
dev.off()

ihc_1 = read_excel("~/Desktop/Lab/Brain_meta/Lung/IHC plots.xlsx",sheet = 1)
ihc_2 = read_excel("~/Desktop/Lab/Brain_meta/Lung/IHC plots.xlsx",sheet = 2)
ihc_3 = read_excel("~/Desktop/Lab/Brain_meta/Lung/IHC plots.xlsx",sheet = 3)

ihc_1.1 = melt(ihc_1)
colnames(ihc_1.1) = c("Sample","Score")
ihcplot1 <- ggpaired(ihc_1.1, x = "Sample", y = "Score",
         color = "Sample", line.color = "gray", line.size = 0.2,
         palette = c("#1b98e0", "red"))+
  #ggtitle("UQCRC2") + 
  stat_compare_means(paired = TRUE,label.x = 1.2,label.y = 15,label = "p.signif") +
  labs(y="UQCRC2") + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = "none")
ihcplot1

ihc_2.1 = melt(ihc_2)
colnames(ihc_2.1) = c("Sample","Score")
ihcplot2 <- ggpaired(ihc_2.1, x = "Sample", y = "Score",
                     color = "Sample", line.color = "gray", line.size = 0.2,
                     palette = c("#1b98e0", "red"))+
  #ggtitle("MTCO1") + 
  stat_compare_means(paired = TRUE,label.x = 1.2,label.y = 15,label = "p.signif") +
  labs(y="MTCO1") + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = "none")
ihcplot2

ihc_3.1 = melt(ihc_3)
colnames(ihc_3.1) = c("Sample","Score")
ihcplot3 <- ggpaired(ihc_3.1, x = "Sample", y = "Score",
                     color = "Sample", line.color = "gray", line.size = 0.2,
                     palette = c("#1b98e0", "red"))+
  #ggtitle("COXIV") + 
  stat_compare_means(paired = TRUE,label.x = 1.2,label.y = 15,label = "p.signif") +
  labs(y="COXIV") + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = "none")
ihcplot3

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/IHCplot3.tiff",
     width=600,height=1500,res=300)
ihcplot3

dev.off()

rppa1 = read_excel("~/Desktop/Lab/Brain_meta/Lung/RPPA plots.xlsx",sheet = 1)
rppa2 = read_excel("~/Desktop/Lab/Brain_meta/Lung/RPPA plots.xlsx",sheet = 2)
rppa3 = read_excel("~/Desktop/Lab/Brain_meta/Lung/RPPA plots.xlsx",sheet = 3)
rppa4 = readxl::read_excel("~/Desktop/Lab/Brain_meta/Lung/RPPA plots.xlsx",sheet = 4)

rppa5 = readxl::read_excel("~/Downloads/Copy of CD3.xlsx")
rppa5 = rppa5[,c(2,5)]

rppa6 = readxl::read_excel("~/Desktop/Lab/Brain_meta/Copy of CD3CD68-LC.xlsx",sheet = 1)
rppa6 = rppa6[-1,9]
rppa6.1 = data.frame("Primary" = rppa6$...9[1:28], 
                     "BrM" = rppa6$...9[29:56])
rppa6.1 <- as.data.frame(sapply(rppa6.1, as.numeric))

colnames(rppa5) = colnames(rppa4)

rppa1.1 = melt(rppa1)
rppa2.1 = melt(rppa2)
rppa3.1 = melt(rppa3)
rppa4.1 = melt(rppa4)
rppa5.1 = melt(rppa5)
rppa6.2 = melt(rppa6.1)
``
colnames(rppa1.1) = c("Sample","Score")
colnames(rppa2.1) = c("Sample","Score")
colnames(rppa3.1) = c("Sample","Score")
colnames(rppa4.1) = c("Sample","Score")
colnames(rppa5.1) = c("Sample","Score")
colnames(rppa6.2) = c("Sample","Score")

ippaplot1 <- ggpaired(rppa1.1, x = "Sample", y = "Score",
                     color = "Sample", line.color = "gray", line.size = 0.2,
                     palette = c("#1b98e0", "red"))+
  stat_compare_means(paired = TRUE,label.x = 1.2,label.y = 1.5,label = "p.signif",method = "t.test") +
  labs(y="UQCRC2") + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = "none")
ippaplot1

ippaplot2 <- ggpaired(rppa2.1, x = "Sample", y = "Score",
                      color = "Sample", line.color = "gray", line.size = 0.2,
                      palette = c("#1b98e0", "red"))+
  stat_compare_means(paired = TRUE,label.x = 1.2,label.y = 4,label = "p.signif") +
  labs(y="MTCO1") + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = "none")
ippaplot2

ippaplot3 <- ggpaired(rppa3.1, x = "Sample", y = "Score",
                      color = "Sample", line.color = "gray", line.size = 0.2,
                      palette = c("#1b98e0", "red"))+
  stat_compare_means(paired = TRUE,label.x = 1.2,label.y = 4,label = "p.signif") +
  labs(y="COXIV") + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = "none")
ippaplot3

ippaplot4 <- ggpaired(rppa4.1, x = "Sample", y = "Score",
                      color = "Sample", line.color = "gray", line.size = 0.2,
                      palette = c("#1b98e0", "red"))+
  stat_compare_means(paired = TRUE,label.x = 1.2,label.y = 2,label = "p.signif") +
  labs(y="SDHA") + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = "none")

compare_means(Score~Sample, data = rppa1.1, method="t.test",paired = TRUE)

ippaplot4

ippaplot5 <- ggpaired(rppa5.1, x = "Sample", y = "Score",
                      color = "Sample", line.color = "gray", line.size = 0.2,
                      palette = c("#1b98e0", "red"))+
  stat_compare_means(paired = TRUE,label.x = 1.2,label.y = 3000,label = "p.signif",method = "t.test") +
  labs(y = expression('CD3'*~ '(Num'~'positive'~'per'~'mm'^2*')')) + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = "none")
ippaplot5

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/RPPAplot5.tiff",
     width=650,height=1500,res=300)
ippaplot5

dev.off()

ippaplot6 <- ggpaired(rppa6.2, x = "Sample", y = "Score",
                      color = "Sample", line.color = "gray", line.size = 0.2,
                      palette = c("#1b98e0", "red"))+
  stat_compare_means(paired = FALSE,label.x = 1.2,label.y = 3000,label = "p.signif",method = "t.test") +
  labs(y= expression('CD68'*~ '('*'N'*'u'*'m'~'p'*'o'*'s'*'i'*'t'*'i'*'v'*'e'~'p'*'e'*'r'~'m'*'m'^2*')')) + 
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        legend.title = element_blank()) +
  theme(legend.position = "none")

ippaplot6

tiff("/Users/jianlan/Desktop/Lab/Brain_meta/Lung/plots/RPPAplot6.tiff",
     width=650,height=1500,res=300)
ippaplot6

dev.off()
