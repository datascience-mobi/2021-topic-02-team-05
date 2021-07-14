#heatmap of all our genes of interest on all 14 chips
#with the pheatmap function
library(pheatmap)
setwd("/Users/Line/Documents/tra_project/thyroid_cancer_sessions/plots")
pheatmap(thyroid.expr.matrix.sub, 
         main = "gene expression of thyroid TRA's \nin PTC and healthy thyroid tissue (GSE35570) \n(Handkiewicz-Junak et all., 2016)",
         cluster_rows = F, 
         cluster_cols = F, 
         cellheight = 3, 
         cellwidth = 3,
         legend =TRUE,
         fontsize = 3,
         filename = "pheatmap_thyroid_TRA_expression.pdf",
         gaps_col = c(5,9)
        )



#distribution our genes over chromosomes (out of human TRA lists=> List e and f can be ignored!)
#preparation => I need a matrix with the information "chromosome location" for each gene
tiss1=unique(as.character(b[,11])) # "thyroid"
th.genes1= b [b$Max_tissue %in% tiss1[31],]
dim(th.genes1)

tiss2= unique(as.character(c[,11])) 
length(tiss2) 
th.genes2= c[c$max.tissue %in% tiss2[35],  ]
dim(th.genes2)

tiss3=unique(as.character(d[,11])) #"thyroid_gland"
th.genes3= d[d$max.tissue %in% tiss3[20],]
dim(th.genes3)

tiss6=unique(as.character(g[,10])) 
th.genes6= g[g$max.tissue %in% tiss6[44],] # "Thyroid"
dim(th.genes6)

th.genes6[,11]= as.factor("no information") # we have no information about chromosome location of these genes

genes1= data.frame(th.genes1[,6],th.genes1[,2])
genes2= data.frame(th.genes2[,3],th.genes2[,7])
genes3= data.frame(th.genes3[,3],th.genes3[,7])
genes6= data.frame(th.genes6[,3],th.genes6[,11])

colnames(genes1)= c("gene.symbol", "chromosome")
colnames(genes2)= c("gene.symbol", "chromosome")
colnames(genes3)= c("gene.symbol", "chromosome")
colnames(genes6)= c("gene.symbol", "chromosome")

allgenesandchromosomes= rbind(genes1,genes2,genes3,genes6) # all genes of all TRA lists together
uniquegenes= duplicated(allgenesandchromosomes$gene.symbol) # selecting unique genes only
Index= which(uniquegenes==FALSE )
uniqueinfos=allgenesandchromosomes[Index,]

indexx=which(uniqueinfos$gene.symbol %in% symbol) # which genes are examined on our microarrays?
finalinfos= uniqueinfos[indexx,]
Ychrom= data.frame("inventedGene", "Y")
colnames(Ychrom)= c("gene.symbol", "chromosome")
finalinfos= rbind(finalinfos, Ychrom)
View(finalinfos) 

#final plots
chrnumbers= table(finalinfos[,2]) #wie ift ein chromosme vorkommt 
chromosomes= chrnumbers[c(2,13,17:23,3:12,14:16,26)]
chromswith0= c(chromosomes,chrnumbers[25])
chromswithY= c(chromosomes, chrnumbers[27])
chromswithY[24]= 0
View(chromosomes)

library(viridisLite)
par(mai= c(2,2,2,2))
pdf(file="piechart_chromosomes.pdf")
pie(chromswith0, clockwise= TRUE,border= FALSE, col= c(plasma(23),"lightgoldenrodyellow"), init.angle=180, radius=1, main=" Distribution of \n thyroid specific genes \n over chromosomes \n Dinkelacker, 2019, PhD thesis, University of Heidelberg")
dev.off()
pdf(file="barplot_chromosomes.pdf")
par(las=1)
barplot(chromswithY, col=plasma(23), cex.names=0.5,main=" Distribution of thyroid specific genes over chromosomes \n Dinkelacker, 2019, PhD thesis, University of Heidelberg", xlab= "chromosomes", ylab= "number of genes", ylim=c(0,20))
abline(h=c(5,10,15,20), lty=3) # makes the grey lines. Tiping the barplot-code again, but with "add=TRUE" puts the lines in the background
barplot(chromswithY, col=plasma(23), add=TRUE, cex.names=0.5,main=" Distribution of thyroid specific genes over chromosomes \n Dinkelacker, 2019, PhD thesis, University of Heidelberg", xlab= "chromosomes", ylab= "number of genes", ylim=c(0,20))
dev.off()
