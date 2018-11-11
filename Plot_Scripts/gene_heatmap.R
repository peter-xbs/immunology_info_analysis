#library(BiocInstaller)
#biocLite("GEOquery")
library(GEOquery)
library(gplots)

# gse_id <- "GSE60290"
# gse_matrix <- getGEO(gse_id,GSEMatrix=TRUE)
# gse_eset = exprs(gse_matrix[[1]])
# eset <- cbind(gse_eset[,1:2], gse_eset[,8:10])
data = read.table('GSE53986.heatmap', header = TRUE, sep = '\t')

patientcolors = c("blue", "blue", "blue", "blue", "red", "red", "red", "red")
# patientcolors = c("blue", "blue", "blue", "blue", "yellow", "yellow", "yellow", "yellow")

g1 = subset(data, key3=="Y")
rownames(g1) = g1$Gene
g2 = data.matrix(g1[,2:9])
heatmap.2(g2, col=greenred, scale="row", ColSideColors = patientcolors,
          key=T,lwid=c(0.45,1), keysize=0.6, key.par = list(cex=0.62),
          density.info="none",dendrogram = 'row', Colv = F, 
          trace="none", cexRow=0.0001, cexCol = 0.01,
          margins = c(5,10), #sepwidth=c(0.00000001,0.00000001) ,sepcolor=c("grey"),
          #colsep = 1:ncol(g2), rowsep = 1:nrow(g2),
          lhei = c(1,8))
# legend("topleft",horiz=F,   # location of the legend on the heatmap plot
#        legend = c("LPS", "Ctrl"), 
#        box.lwd = 0,
#        inset = 0.01,
#        col = c("red", "blue"),  # color key
#        lty = 1,            # line style
#        lwd = 6         # line width
# )


#dev.off()

