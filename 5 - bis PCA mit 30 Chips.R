###kompletter angepasster Code mit 30 Chips und verändertem Working Directory, 
###um vorherige Plots nicht zu überschreiben 
###(clustering und pca weg gelassen wenn nicht sinnvoll)

###quality control

#packages 
library(affy)
library(vsn)
library(AnnotationDbi)
library(hgu133plus2hsenstcdf)
library(hgu133plus2hsenstprobe)

#read in .CEL files
setwd("/Users/Line/Documents/tra_project/thyroid_cancer_rawdata")
data.thyroid = ReadAffy()
data.thyroid@cdfName <- "HGU133Plus2_Hs_ENST"

setwd("/Users/Line/Documents/tra_project/with_30_chips/rda")
save.image(file = "rawdata.thyroid.rda")

#quality control

#single chip control 
setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
pdf(file="single_chip_control.pdf")
image(data.thyroid, col = rainbow(14, start = 0, end = 0.75)[10:1])
dev.off()
#alle chips scheinen okay 

#normalization
thyroid.vsnrma<-vsnrma(data.thyroid)
setwd("/Users/Line/Documents/tra_project/with_30_chips/rda")
save.image(file = "normalized_data.rda")

#meanSdPlot
setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
install.packages("hexbin")
library(hexbin)
pdf(file = 'meanSdPlot_thyroid_vsnrma_normalized.pdf')
meanSdPlot(thyroid.vsnrma)
dev.off()

#Boxplot
#before normalization
pdf(file = 'boxplot_thyroid_rawdata.pdf')
par(mar=c(9,5,3,5))
boxplot(data.thyroid,
        col=rainbow(14), 
        cex.axis=0.5, 
        ylab="intensity", 
        las=2, 
        main="Gene expression in pappillary thyroid cancer (GSE35570) \nbefore normalization (Handkiewicz-Junak et all., 2016)")
abline(h=c(6,8,10,12,14), lty=3)
boxplot(data.thyroid,
        col=rainbow(14), 
        cex.axis=0.5, 
        ylab="intensity", 
        las=2, 
        add=TRUE,
        main="Gene expression in pappillary thyroid cancer (GSE35570) \nbefore normalization (Handkiewicz-Junak et all., 2016)")
dev.off()

#after normalization
pdf(file = 'boxplot_thyroid_vsnrma_normalized.pdf')
par(mar=c(9,5,3,5))
boxplot(exprs(thyroid.vsnrma), 
        col=rainbow(14), 
        cex.axis=0.5, 
        ylab="log intensity", 
        las=2, 
        main="Gene expression in papillary thyroid cancer (GSE35570) \nafter vsnrma normalization (Handkiewicz-Junak et all., 2016)")
abline(h=c(6,8,10,12,14), lty=3)
boxplot(exprs(thyroid.vsnrma), 
        col=rainbow(14), 
        cex.axis=0.5, 
        ylab="log intensity", 
        las=2, 
        add=TRUE,
        main="Gene expression in papillary thyroid cancer (GSE35570) \nafter vsnrma normalization (Handkiewicz-Junak et all., 2016)")
dev.off()

#density function 
#raw data
setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
pdf(file = 'hist_thyroid_rawdata.pdf')
hist(data.thyroid, 
     col=rainbow(14), 
     main="Density function of log Intensity of papillary thyroid cancer (GSE35570) \nbefore normalization (Handkiewicz-Junak et all., 2016)")
abline(h=c(0.1,0.3,0.5), lty=3)
hist(data.thyroid, 
     col=rainbow(14), 
     add=TRUE,
     main="Density function of log Intensity of papillary thyroid cancer (GSE35570) \nbefore normalization (Handkiewicz-Junak et all., 2016)")
dev.off()

#normalized data
setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
pdf(file = 'hist_thyroid_vsnrma_normalised.pdf')
plot(density(exprs(thyroid.vsnrma)[,1]), 
     type="n", 
     xlab="log Intensity", 
     ylim = c(0,1), 
     main="Density function of log Intensity of papillary thyroid cancer (GSE35570) \nafter vsnrma normalization (Handkiewicz-Junak et all., 2016)")
for (i in 1:ncol(exprs(thyroid.vsnrma))) {
  lines(density(exprs(thyroid.vsnrma)[,i]), col=rainbow(14)[i])
}
abline(h=c(0.2,0.4,0.6,0.8), lty=3)
plot(density(exprs(thyroid.vsnrma)[,1]), 
     type="n", 
     xlab="log Intensity", 
     ylim = c(0,1), 
     main="Density function of log Intensity of papillary thyroid cancer (GSE35570) \nafter vsnrma normalization (Handkiewicz-Junak et all., 2016)")
for (i in 1:ncol(exprs(thyroid.vsnrma))) {
  lines(density(exprs(thyroid.vsnrma)[,i]), col=rainbow(14)[i])
}
dev.off()

#RNA degradation plot 
pdf(file = 'rnadeg_human_thyroid_cancer_rawdata.pdf')
rnadeg.raw = AffyRNAdeg(data.thyroid)
plotAffyRNAdeg(rnadeg.raw, col=rainbow(14))
title(sub="human thyroid cancer rawdata (GSE35570) (Handkiewicz-Junak et all., 2016)")
dev.off()

pdf(file = 'rnadeg_human_thyroid_cancer_rawdata_shifted.pdf')
plotAffyRNAdeg(rnadeg.raw, col=rainbow(14), transform="shift.only")
title(sub="human thyroid cancer rawdata (GSE35570) (Handkiewicz-Junak et all., 2016)")
dev.off()

#scatterplots
setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
chipnames.CEL=colnames(exprs(thyroid.vsnrma))
chipnames = substr(chipnames.CEL, 1, nchar(chipnames.CEL)-4) #creating a vector of the chipnames without the .CEL
for (i in 1:13) {
  for (j in (i+1):14) {
    pdf(file = paste('scatterplot_thyroid_vsnrma_normalised', i, "_", j,'.pdf'))
    plot(exprs(thyroid.vsnrma)[,c(i,j)], pch=".", 
         main=paste("plot of probe", chipnames[i], "\nand", chipnames[j], "\n(Handkiewicz-Junak et all., 2016)"))
    abline(0,1,col="red")
    dev.off()
  }
}





### gene annotation

#extract expression values from normalized data 
thyroid.expr.matrix=exprs(thyroid.vsnrma)
head(thyroid.expr.matrix)
ensemble.id = rownames(thyroid.expr.matrix)
expr=substr(ensemble.id,0,15)
rownames(thyroid.expr.matrix)=expr


#read in annotation file 
setwd("/Users/Line/Documents/tra_project/thyroid_annotation_file")
a=read.csv("ensembl_103.txt",sep="\t")
head(a)
transcript.id.ensemble=as.character(a[,2])
symbol.ensembl=as.character(a[,4])
names(symbol.ensembl)=transcript.id.ensemble
head(symbol.ensembl)

#look into data matrix
head(thyroid.expr.matrix)
ensemble.id = rownames(thyroid.expr.matrix)
symbol=symbol.ensembl[ensemble.id]
head(symbol)
dim(thyroid.expr.matrix)
length(symbol)

#re-apply rownames to gene symbols in the data.matrix
rownames(thyroid.expr.matrix)=as.character(symbol)
head(thyroid.expr.matrix)

#save matrix as csv
setwd("/Users/Line/Documents/tra_project/with_30_chips/tables")
write.csv(thyroid.expr.matrix, "thyroid_gene_expression")

##This is the matrix you can calculate with##




###extraction of genes of interest

#extract genes of interest
setwd("/Users/Line/Documents/tra_project/tra_tables")
b=read.csv(file="Human_protein_atlas_TRA_5median_genes_annotated.tsv",sep=",")
tiss=b[,11]
ind=which(tiss=="thyroid")
TRA.symbol=b[,6]
thyroid.TRA1=TRA.symbol[ind]

c=read.csv(file="tra.2014.human.5x.table.tsv",sep="\t")
tiss=c[,11]
ind=which(tiss=="Thyroid")
TRA.symbol=c[,3]
thyroid.TRA2=TRA.symbol[ind]

d=read.csv(file="tra.2014.human.roth.5x.table.tsv",sep="\t")
tiss=d[,11]
ind=which(tiss=="thyroid_gland")
TRA.symbol=d[,3]
thyroid.TRA3=TRA.symbol[ind]

e=read.csv(file="tra.2014.mouse.5x.table.tsv",sep="\t")
tiss=e[,11]
ind=which(tiss=="thyroid")
TRA.symbol=e[,3]
thyroid.TRA4=TRA.symbol[ind]

f=read.csv(file="tra.2014.mouse.4301.5x.table.tsv",sep="\t")
tiss=f[,11]
ind=which(tiss=="thyroid")
TRA.symbol=f[,3]
thyroid.TRA5=TRA.symbol[ind]

g=read.csv(file="tra.2017.human.gtex.5x.table.tsv",sep="\t")
tiss=g[,10]
ind=which(tiss=="Thyroid")
TRA.symbol=g[,3]
thyroid.TRA6=TRA.symbol[ind]

thyroid.TRA = union(thyroid.TRA1, thyroid.TRA2)
thyroid.TRA = union(thyroid.TRA, thyroid.TRA3)
thyroid.TRA = union(thyroid.TRA, thyroid.TRA4)
thyroid.TRA = union(thyroid.TRA, thyroid.TRA6)

#sometimes upper case and lower case and other problems disturb the search (for example levels) 
#which you can eliminate with as.character and toupper

row.ind=which(toupper(thyroid.TRA) %in% as.character(symbol))
thyroid.expr.matrix.sub=thyroid.expr.matrix[row.ind,]

head(thyroid.expr.matrix.sub)
dim(thyroid.expr.matrix.sub)

#deleting the first 46 rows (AFFX...)
thyroid.expr.matrix.sub <- thyroid.expr.matrix.sub[-c(1:46),]
#reorder the matrix
thyroid.expr.matrix.sub <- thyroid.expr.matrix.sub[,c(1,5,6,7,8,10,11,12,13,19,2,3,4,9,14,15,16,17,18,20,21,22,23,24,25,26,27,28,29,30)]

#save matrix as csv
setwd("/Users/Line/Documents/tra_project/with_30_chips/tables")
write.csv(thyroid.expr.matrix.sub, "thyroid_gene_expression_TRA")

#boxplot
par(las=2)
boxplot(t(thyroid.expr.matrix.sub),col=rainbow(length(rownames(thyroid.expr.matrix.sub))),main="gene expression of thyroid-specific genes in thyroid cancer",cex.axis=0.8)
#sort alphabetically
setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
pdf(file="boxplot_thyroid_cancer_genes_of_interest.pdf")
order.vector=sort(colnames(t(thyroid.expr.matrix.sub)),index.return=T)$ix
boxplot(t(thyroid.expr.matrix.sub)[,order.vector],
        col=rainbow(length(rownames(thyroid.expr.matrix.sub))),
        main="gene expression of thyroid-specific genes in papillary thyroid cancer \n(GSE35570) (Handkiewicz-Junak et all., 2016)",
        cex.axis=0.8)
abline(h=c(6,8,10,12,14), lty=3)
boxplot(t(thyroid.expr.matrix.sub)[,order.vector],
        col=rainbow(length(rownames(thyroid.expr.matrix.sub))),
        add=TRUE,
        main="gene expression of thyroid-specific genes in papillary thyroid cancer \n(GSE35570) (Handkiewicz-Junak et all., 2016)",
        cex.axis=0.8)

dev.off()

#create three matrices (radiation_exposed/not_exposed/healthy)
radexp.matrix=thyroid.expr.matrix.sub[,c(1,2,3,4,5,6,7,8,9,10)]
radnonexp.matrix=thyroid.expr.matrix.sub[,c(11,12,13,14,15,16,17,18,19,20)]
healthy.matrix=thyroid.expr.matrix.sub[,c(21,22,23,24,25,26,27,28,29,30)]

###descriptive statistics

#heatmap of all our genes of interest on all 30 chips
#with the pheatmap function
library(pheatmap)
setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
pheatmap(thyroid.expr.matrix.sub, 
         main = "gene expression of thyroid TRA's \nin PTC and healthy thyroid tissue (GSE35570) \n(Handkiewicz-Junak et all., 2016)",
         cluster_rows = F, 
         cluster_cols = F, 
         cellheight = 3, 
         cellwidth = 3,
         legend =TRUE,
         fontsize = 3,
         filename = "pheatmap_thyroid_TRA_expression.pdf",
         gaps_col = c(10,20)
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

###PCA

#PCA of genes
pca.t <- prcomp(t(thyroid.expr.matrix.sub), center = TRUE, scale = TRUE)
summary(pca.t) #with PC1 and PC2 40,87% of the variance are explained
col.pca.t = ifelse(colnames(thyroid.expr.matrix.sub) %in% colnames(radexp.matrix), 
                   "red", 
                   ifelse(colnames(thyroid.expr.matrix.sub) %in% colnames(radnonexp.matrix), 
                          "blue", 
                          "green"))
pdf(file="pca_genes")
plot(pca.t$x[,1], pca.t$x[,2], 
     col= col.pca.t,
     pch=16,
     xlab='PC1 (25,01 %)',ylab='PC2 (15,86 %)',
     main = "PC1 vs. PC2 for radiation-exposed PTC, \nnot-exposed PTC and healthy thyroid tissue \n(Handkiewicz-Junak et all., 2016)") 
legend("topright",c("radiation-exposed PTC","not-exposed PTC", "healthy"),fill=c("red","blue", "green"))

dev.off()

# Percentage of variance explained by dimensions
install.packages('factoextra')
library("factoextra")
eigenvalue <- round(get_eigenvalue(pca.t), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)

#k-means 
km3 = kmeans(t(thyroid.expr.matrix.sub),centers=3,nstart = 100)
km3$cluster
cluster1.1 = c(which(km3$cluster==1))
cluster2.1 = c(which(km3$cluster==2))
cluster3.1 = c(which(km3$cluster==3))

library("factoextra")

setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
pdf(file="pca_genes_kmeans_3")
fviz_cluster(km3, data = t(thyroid.expr.matrix.sub),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
dev.off()

km2 = kmeans(t(thyroid.expr.matrix.sub),centers=2,nstart = 100)
km2$cluster
cluster1.2 = c(which(km2$cluster==1))
cluster2.2 = c(which(km2$cluster==2))
cluster3.2 = c(which(km2$cluster==3))

pdf(file="pca_genes_kmeans_2")
fviz_cluster(km2, data = t(thyroid.expr.matrix.sub),
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
dev.off()

#PCA of the sick patients only
pca.sick <- prcomp(t(thyroid.expr.matrix.sub)[1:20,], center = TRUE, scale = TRUE)
summary(pca.sick) #with PC1 and PC2 40,04% of the variance are explained
col.pca.sick = ifelse(colnames(thyroid.expr.matrix.sub) %in% colnames(radexp.matrix), "red", "blue")
setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
pdf(file="pca_sick")
plot(pca.sick$x[,1], pca.sick$x[,2], 
     col= col.pca.sick,
     pch=16,
     xlab='PC1 (26,06 %)',ylab='PC2 (13,99 %)',
     main = "PC1 vs. PC2 for radiation-exposed PTC \nand not-exposed PTC \n(Handkiewicz-Junak et all., 2016)") 
legend("bottomright",c("radiation-exposed PTC","not-exposed PTC"),fill=c("red","blue"))

dev.off()

#computing the best number of clusters 
#through the elbow method
wss.sick = sapply(1:7,function(k) { 
  kmeans(x=thyroid.expr.matrix.sub[,1:20], centers =k)$tot.withinss
})

plot(1:7,wss.sick,type='b',pch=19,xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares") #3 seems to be the best number of clusters

#k-means of pca with sick patiets only
#with 2 clusters because of 2 groups
km.sick2 = kmeans(t(thyroid.expr.matrix.sub)[1:20,],centers=2,nstart = 100)
km.sick2$cluster

library("factoextra")

setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
pdf(file="pca_sick_kmeans_2")
fviz_cluster(km.sick2, data = t(thyroid.expr.matrix.sub)[1:20,],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
dev.off()

#with 3 clusters because of elbow method
km.sick3 = kmeans(t(thyroid.expr.matrix.sub)[1:20,],centers=3,nstart = 100)
km.sick3$cluster

pdf(file="pca_sick_kmeans_3")
fviz_cluster(km.sick3, data = t(thyroid.expr.matrix.sub)[1:20,],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw(),
             show.clust.cent = FALSE,
)
dev.off()
#3 clusters are a better fit




