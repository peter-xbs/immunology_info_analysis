#source("https://bioconductor.org/biocLite.R")
#library(BiocInstaller)
# 
# install.packages("ggrepel")
library(Biobase)
# library(GEOquery)
library(limma)
library(ggrepel)
fit2 = read.table('GSE53986.process.diff', header = TRUE, sep = '\t')
#tT <- topTable(fit2, adjust="fdr", sort.by="logFC", number=1000000)

library(ggplot2)
## Volcano plot
require(ggplot2)
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
gene_list = fit2
no_of_genes = dim(gene_list)[1]
#gene_list$threshold = as.factor(abs(gene_list$logFC) >=1 & gene_list$adj.P.Val < 0.05)
gene_list$threshold = as.factor(gene_list$adj.P.Val < 0.05)

gene_list$color_flag <- ifelse(gene_list$threshold==TRUE & gene_list$logFC > 1.8, "upregulation", ifelse(gene_list$threshold==TRUE & gene_list$logFC < -1.8, "downregulation", "no-significance"))

##Construct the plot object
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(adj.P.Val), colour=color_flag)) +
  geom_point(alpha=1, size=1) +
  scale_colour_manual(values=c("blue","grey","red"))+
  theme_bw(base_size = 16) + theme(legend.position = "bottom") +
  xlim(c(-5, 5)) + ylim(c(0, 8)) +
  xlab("log2 fold change") + ylab("-log10 adj.p.value")+
  geom_label_repel(
    # data = subset(gene_list, Gene.symbol %in% c("Zc3h12a","Mmp9","Il1b","Ifitm2","Gapdh","Cybb","Cd14","Calcrl","Ass1","Aoah","Tlr4","Sirt3")),
    #data = subset(gene_list, key3=="Y" & color_flag %in% c("downregulation", "upregulation")),
    data = subset(gene_list, key12=="Y" & threshold==TRUE),
    aes(label = Gene.symbol),
    size = 5,
    segment.colour = "black",
    color="black",
    #direction = 'x',
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.5, "lines")
  )
  #theme_classic(base_rect_size = 16)
g

#anti_inflam = c("Tnfrsf1b","Tnfaip3","Tecpr1","Stat1","Pparg",
#                "Nr1h3","Nos2","Il1rn","Gbp7","Bst2")

#pro_inflam = c("Zc3h12a","Mmp9","Il1b","Ifitm2","Gapdh",
#               "Cybb","Cd14","Calcrl","Ass1","Aoah","Tlr4")

#glocolysis = c("Tpi1","Pgk1","Pfkp","Pfkl","Hif1a","Gpi1","Gapdh","Aldoa",
#               "Ldhb","Glut1","Eno1","Aco1","Aldoc","Slc2a1","Slc2a6","Hk1",
#               "Hk2","Pfkp","Pgm1","Pgm2","Pkm2","Ldha")

#3mitocondrial = c("Mfn1","Mfn2","Drp1","Opa1","Fis1","Gdap1","Dnm1")
