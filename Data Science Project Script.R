####Data Science Project####
###Camila Vacas
###07.06.2021

#### 1. Quality Control####

###Preparation
#packages 
library(affy)
library(vsn)
library(AnnotationDbi)
library(hgu133plus2hsenstcdf)
library(hgu133plus2hsenstprobe)

#read in .CEL files
setwd("~/Desktop/Bioinfo/Data Science Project/rawdata/GSE35570 papillary thyroid cancer (PTC)")
data.thyroid = ReadAffy()
data.thyroid@cdfName <- "HGU133Plus2_Hs_ENST"

save.image(file = "rawdata.thyroid.rda")

#quality control

#single chip control 
setwd("~/Desktop/Bioinfo/Data Science Project/plots")
pdf(file="single_chip_control.pdf")
image(data.thyroid, col = rainbow(14, start = 0, end = 0.75)[10:1])
dev.off()
#alle chips scheinen okay 

#normalization
thyroid.vsnrma<-vsnrma(data.thyroid)
setwd("~/Desktop/Bioinfo/Data Science Project/rawdata")
save.image(file = "normalized_data.rda")

#meanSdPlot
setwd("~/Desktop/Bioinfo/Data Science Project/plots")
#install.packages("hexbin")
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
dev.off()

#density function 
#raw data
pdf(file = 'hist_thyroid_rawdata.pdf')
hist(data.thyroid, col=rainbow(14), main="Density function of log Intensity of papillary thyroid cancer (GSE35570) \nbefore normalization (Handkiewicz-Junak et all., 2016)")
dev.off()

#normalized data
pdf(file = 'hist_thyroid_vsnrma_normalised.pdf')
plot(density(exprs(thyroid.vsnrma)[,1]), type="n", xlab="log Intensity", ylim = c(0,1), main="Density function of log Intensity of papillary thyroid cancer (GSE35570) \nafter vsnrma normalization (Handkiewicz-Junak et all., 2016)")
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

#### 2. Gene Annotation####

#extract expression values from normalized data 
thyroid.expr.matrix=exprs(thyroid.vsnrma)
head(thyroid.expr.matrix)
ensemble.id = rownames(thyroid.expr.matrix)
expr=substr(ensemble.id,0,15)
rownames(thyroid.expr.matrix)=expr


#read in annotation file 
setwd("~/Desktop/Bioinfo/Data Science Project/TRA Daten")
a=read.csv("ensembl_103.txt",sep=",")
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
length(symbol)#95721

#re-apply rownames to gene symbols in the data.matrix
rownames(thyroid.expr.matrix)=as.character(symbol)
head(thyroid.expr.matrix)

#save matrix as csv
setwd("~/Desktop/Bioinfo/Data Science Project/tables")
write.csv(thyroid.expr.matrix, "thyroid_gene_expression")

##This is the matrix you can calculate with##



#### 3. Extraction genes of interest####

#extract genes of interest
setwd("~/Desktop/Bioinfo/Data Science Project/TRA Daten")
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
dim(thyroid.expr.matrix.sub)#406 x 30

#deleting the first 46 rows (AFFX...)
thyroid.expr.matrix.sub <- thyroid.expr.matrix.sub[-c(1:46),]
#reorder the matrix
thyroid.expr.matrix.sub <- thyroid.expr.matrix.sub[,c(1,5,6,7,10,11,8,19,13,12,18,2,3,4,9,14,15,16,17,20,21,22,23,24,25,26,27,28,29,30)]

#save matrix as csv
setwd("~/Desktop/Bioinfo/Data Science Project/tables")
write.csv(thyroid.expr.matrix.sub, "thyroid_gene_expression_TRA")

#boxplot
par(las=2)
boxplot(t(thyroid.expr.matrix.sub),col=rainbow(length(rownames(thyroid.expr.matrix.sub))),main="gene expression of thyroid-specific genes in thyroid cancer",cex.axis=0.8)
#sort alphabetically
setwd("~/Desktop/Bioinfo/Data Science Project/plots")
pdf(file="boxplot_thyroid_cancer_genes_of_interest.pdf")

order.vector=sort(colnames(t(thyroid.expr.matrix.sub)),index.return=T)$ix
boxplot(t(thyroid.expr.matrix.sub)[,order.vector],col=rainbow(length(rownames(thyroid.expr.matrix.sub))),main="gene expression of thyroid-specific genes in papillary thyroid cancer \n(GSE35570) (Handkiewicz-Junak et all., 2016)",cex.axis=0.8)

dev.off()

#create three matrices (radiation_exposed/not_exposed/healthy)
radexp.matrix=thyroid.expr.matrix.sub[,c(1,5,6,7,8)]
radnonexp.matrix=thyroid.expr.matrix.sub[,c(2,3,4,9)]
healthy.matrix=thyroid.expr.matrix.sub[,c(10,11,12,13,14)]


#### 4. Descriptive Statistics ####


#heatmap of all our genes of interest on all 14 chips
#with the pheatmap function
library(pheatmap)
setwd("~/Desktop/Bioinfo/Data Science Project/plots")
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



#distribution our genes over chromosomes
#preparation=> I need a matrix with the information "chromosome location" for each gene
tiss1=unique(as.character(b[,11])) # "thyroid"
th.genes1= b [b$Max_tissue %in% tiss1[31],]
dim(th.genes1)

tissues2= unique(as.character(c[,11])) 
length(tissues2) 
th.genes2= c[c$max.tissue %in% tissues2[35],  ]
dim(th.genes2)

tiss3=unique(as.character(d[,11])) #"thyroid_gland"
th.genes3= #distribution our genes over chromosomes (out of human TRA lists=> List e and f can be ignored!)
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

#### 5.Clustering and PCA ####

#computing the best number of clusters 
#through the elbow method
wss = sapply(1:7,function(k) { 
  kmeans(x=thyroid.expr.matrix.sub, centers =k)$tot.withinss
})

plot(1:7,wss,type='b',pch=19,xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

#k-means
km = kmeans(thyroid.expr.matrix.sub,centers=3,nstart = 100)
km$cluster
cluster1 = c(which(km$cluster==1))
cluster2 = c(which(km$cluster==2))
cluster3 = c(which(km$cluster==3))
length(cluster1)
length(cluster2)
length(cluster3)

#scatterplots of all chips sorted by cluster
setwd("~/Desktop/Bioinfo/Data Science Project/plots")
for (i in 1:13) {
  for (j in (i+1):14) {
    pdf(file = paste('scatterplot_thyroid_cluster_', i, "_", j,'.pdf'))
    plot(thyroid.expr.matrix.sub[,c(i,j)], 
         pch=16, 
         col=km$cluster,
         main=paste("plot of probe", chipnames[i], "\nand", chipnames[j], "sorted by cluster\n(Handkiewicz-Junak et all., 2016)"),
         legend=TRUE)
    
    dev.off()
  }
}

#PCA of chips
pca <- prcomp(thyroid.expr.matrix.sub, center = TRUE, scale = TRUE)
summary(pca) #with PC1 and PC2 97.619% auf the variance are explained
pdf(file="pca_chips")
plot(pca$x[,1], pca$x[,2], 
     col= km$cluster, pch=16,
     xlab='PC1',ylab='PC2',
     main = "PC1 vs. PC2 colored according to kmeans-clustering \n(Handkiewicz-Junak et all., 2016)") 
dev.off()

#PCA of chips with ggbiplot (looks not good)

#install.packages("devtools")
library(devtools)
#install_github("vqv/ggbiplot", force=TRUE)
library(ggbiplot)
pdf(file="ggbiplot_pca_chips")
ggbiplot(pca,
         groups=km$cluster,
         varname.size = 3,
         varname.adjust = 1)
dev.off()

#PCA of genes
pca.t <- prcomp(t(thyroid.expr.matrix.sub), center = TRUE, scale = TRUE)
summary(pca.t) #with PC1 and PC2 44.11% auf the variance are explained
col.pca.t = ifelse(colnames(thyroid.expr.matrix.sub) %in% colnames(radexp.matrix), 
                   "red", 
                   ifelse(colnames(thyroid.expr.matrix.sub) %in% colnames(radnonexp.matrix), 
                          "blue", 
                          "green"))
pdf(file="pca_genes")
plot(pca.t$x[,1], pca.t$x[,2], 
     col= col.pca.t,
     pch=16,
     xlab='PC1',ylab='PC2',
     main = "PC1 vs. PC2 for radiation-exposed PTC (red), \nnot-exposed PTC (blue) and healthy thyroid tissue (green) \n(Handkiewicz-Junak et all., 2016)") 
dev.off()

###PCA (with 30 chips)

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
#install.packages('factoextra')
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

#k-means of pca with sich patiets only
#with 2 clusters because of 2 groups
km.sick2 = kmeans(t(thyroid.expr.matrix.sub)[1:20,],centers=2,nstart = 100)
km.sick2$cluster

library("factoextra")

pdf(file="pca_sick_kmeans_2")
fviz_cluster(km.sick2, data = t(thyroid.expr.matrix.sub)[1:20,],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
dev.off()

#with 3 custers because of elbow method
km.sick3 = kmeans(t(thyroid.expr.matrix.sub)[1:20,],centers=3,nstart = 100)
km.sick3$cluster

pdf(file="pca_sick_kmeans_3")
fviz_cluster(km.sick3, data = t(thyroid.expr.matrix.sub)[1:20,],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
dev.off()#3 clusters are a better fit

#### 6. ANOVA and PostHocTest ####


### Creating a new matrix for ANOVA test

groups= c("radexp","radexp","radexp","radexp","radexp","radexp","radexp","radexp","radexp","radexp",
          "radnonexp","radnonexp","radnonexp","radnonexp","radnonexp","radnonexp","radnonexp","radnonexp","radnonexp","radnonexp",
          "healthy","healthy","healthy","healthy","healthy","healthy","healthy","healthy","healthy","healthy")
anovamatrix=t(thyroid.expr.matrix.sub)
anovamatrix= cbind(groups,anovamatrix)

# anovamatrix[,1] = groups(exp, nonexp, healthy)
# anovamatrix[,i] = genes with their expression data


# ANOVA Test for each gene
sig.genes.vector=c()

for (i in 2:dim(anovamatrix)[2]){
  if (summary.aov(aov(anovamatrix[,i]~ anovamatrix[,1]))[[1]][[1,5]] < 0.05){
    sig.genes.vector = rbind(sig.genes.vector, colnames(anovamatrix)[i])
  }
  
}

sig.genes.matrix= thyroid.expr.matrix.sub[sig.genes.vector,]
# matrix with expression data of our significant diffent expressed genes



### post hoc test for a single gene
#anovatry= aov(anovamatrix[,2] ~ anovamatrix[,1]) 
# anova test for testing the following steps

#posthoc=TukeyHSD(anovatry)
#posthoc # results in matrix with differences in means, lower and upper border and p adjustment for each comparison
#plot(posthoc)


#PostHocTest

#install.packages("multcompView")
library(multcompView)

###preparation for boxplot
uniquegroups= unique(groups)

sig.anova.matrix = cbind(groups, t(sig.genes.matrix))  
# [i,] are the chips, [,i] the genes and [,1] are our groups
sig.genes.names= colnames(sig.anova.matrix[,-c(1)])

my_colors <- c( 
  "blue",
  "green", 
  "red"
) #defines colors

###The boxplots
setwd("~/Desktop/Bioinfo/Data Science Project/plots")

for (i in 2:dim(sig.anova.matrix)[2]){
  pdf(file = paste('posthocgene', colnames(sig.anova.matrix)[i], ".pdf"))
  posthoc= TukeyHSD( aov(sig.anova.matrix[,i] ~ sig.anova.matrix[,1])) # makes posthoc for each gene
  
  Tukey.levels = posthoc[[1]][,4]
  Tukey.labels = data.frame(multcompLetters(Tukey.levels)["Letters"])
  Tukey.labels$groups = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels [order(Tukey.labels$groups), ] #defines the 3 groups, and labels, if group is significant different
  
  value= sig.anova.matrix[,i]
  data= data.frame(groups, value)
  data$value= as.numeric(data$value) 
  
  #boxplot coloured by group
  posthoc.plot= boxplot(data$value~data$groups, 
                        col=my_colors[as.factor(Tukey.labels[,1])],
                        xlab="tissue samples",
                        ylab= "expression values",
                        main= paste("posthoc of gene", sig.genes.names[i], "per group"))
  over <- 0.1*max( posthoc.plot$stats[nrow(posthoc.plot$stats),] )
  text( c(1:nlevels(as.factor(data$groups))) +0.5,
        posthoc.plot$stats[nrow(posthoc.plot$stats),], 
        Tukey.labels[,1], pos=1,
        col=my_colors[as.factor(Tukey.labels[,1])] ) #labels the boxes in the plot
  
  dev.off()
}



#### 7. Analysis of cancer subtypes ####


cluster.sick.1 = c(which(km.sick3$cluster==1))
cluster.sick.2 = c(which(km.sick3$cluster==2))
cluster.sick.3 = c(which(km.sick3$cluster==3))

cluster.sick.1.matrix = thyroid.expr.matrix.sub[,cluster.sick.1]
cluster.sick.2.matrix = thyroid.expr.matrix.sub[,cluster.sick.2]
cluster.sick.3.matrix = thyroid.expr.matrix.sub[,cluster.sick.3]

cancer.subtype.matrix = cbind(cluster.sick.1.matrix, cluster.sick.2.matrix, cluster.sick.3.matrix)

cancer.groups= c("subtype 1","subtype 1","subtype 1","subtype 1","subtype 2","subtype 2","subtype 2","subtype 2","subtype 2",
                 "subtype 3","subtype 3","subtype 3","subtype 3","subtype 3","subtype 3","subtype 3","subtype 3","subtype 3",
                 "subtype 3","subtype 3")

cancer.anovamatrix =cbind(cancer.groups,t(cancer.subtype.matrix))

cancer.sig.genes.vector=c()

#reducing the number of genes by finding out which are significantly different expressed in the three subtypes

for (i in 2:dim(cancer.anovamatrix)[2]){
  if (summary.aov(aov(cancer.anovamatrix[,i]~ cancer.anovamatrix[,1]))[[1]][[1,5]] < 0.01){
    cancer.sig.genes.vector = rbind(cancer.sig.genes.vector, colnames(cancer.anovamatrix)[i])
  }

}

cancer.subtype.sig.genes.matrix= cancer.subtype.matrix[cancer.sig.genes.vector,]
# matrix with expression data of our significant different expressed genes between the three cancer subtypes


#### 8. Regression model ####

#posthoc to find out which genes are significantly expressed in all 
#three subtypes to use those for the regression model
subtype.matrix = rbind(cancer.groups, cancer.subtype.sig.genes.matrix)
subtype.matrix = t(subtype.matrix)
modelgenes.sub = c()

for (i in 2:dim(subtype.matrix)[2]){
  posthoc= TukeyHSD( aov(subtype.matrix[,i] ~ subtype.matrix[,1]))
  Tukey.levels = posthoc[[1]][,4]
  Tukey.labels = data.frame(multcompLetters(Tukey.levels)["Letters"])
  if (Tukey.labels[1,1]== "a" & Tukey.labels[2,1]== "b" & Tukey.labels[3,1]== "c") {
    modelgenes.sub= c(modelgenes.sub, colnames(subtype.matrix)[i])}
}

modelgenes.sub

#control
which(colnames(subtype.matrix)=="SYPL1")

posthoc= TukeyHSD( aov(subtype.matrix[,22] ~ subtype.matrix[,1]))
Tukey.levels = posthoc[[1]][,4]
Tukey.labels = data.frame(multcompLetters(Tukey.levels)["Letters"])

Tukey.labels

#intsalling necessary packages for regression model
#install.packages("nnet")
library("nnet")

#defining subtype 1 as reference category
subtype.matrix[,1]= relevel(factor(subtype.matrix[,1]), ref = "subtype 1")

#regeression model (training)
model.matrix.sub=subtype.matrix[,modelgenes.sub]
model.matrix.sub= cbind(cancer.groups,model.matrix.sub)

#Versuch Expressionswerte numerisch zu machen 
model.matrix.sub=as.data.frame(model.matrix.sub)
for (i in 2:12){
  model.matrix.sub[,i]=apply(model.matrix.sub[,i,drop=F],2,as.numeric)
}

subtype.model=multinom(cancer.groups~.,
                       data = as.data.frame(model.matrix.sub[c(1,2,3,5,6,7,8,14,15,16,17,18,19,20),]))

#regressionmodel (prediction)
predict.subtype.1 = predict(subtype.model, model.matrix.sub[c(1,2,3,5,6,7,8,14,15,16,17,18,19,20),])
predict.subtype.2 = predict(subtype.model, model.matrix.sub[c(4,9,10,11,12,13),])    

#dioganal matrix if all subtypes were identified correctly 
tab1=table(predict.subtype.1, model.matrix.sub[c(1,2,3,5,6,7,8,14,15,16,17,18,19,20),1])
tab1 #is diagonal
tab2=table(predict.subtype.2, model.matrix.sub[c(4,9,10,11,12,13),1])
tab2 #is diagonal
#all subtypes are predicted correctly


#### 9 - boxplots of 3 treatment groups.R ####

thyroid.expr.matrix.new=rbind(thyroid.expr.matrix.sub,thyroid.expr.matrix.sub,thyroid.expr.matrix.sub)

thyroid.expr.matrix.new[c(1:360),11:30]=NA
thyroid.expr.matrix.new[c(361:720),c(1:10,21:30)]=NA
thyroid.expr.matrix.new[c(721:1080),1:20]=NA

order.vector=sort(colnames(t(thyroid.expr.matrix.new)),index.return=T)$ix #gives position of genes in alphabetical order

setwd("~/Desktop/Bioinfo/Data Science Project/plots")

#genes starting with A
pdf(file = "boxplot_genes_of_interest_3groups_A")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[1:57]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with B
pdf(file = "boxplot_genes_of_interest_3groups_B")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[58:84]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with C
pdf(file = "boxplot_genes_of_interest_3groups_C1")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[85:138]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

pdf(file = "boxplot_genes_of_interest_3groups_C2")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[139:192]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with D
pdf(file = "boxplot_genes_of_interest_3groups_D")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[193:237]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with E
pdf(file = "boxplot_genes_of_interest_3groups_E")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[238:270]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with F
pdf(file = "boxplot_genes_of_interest_3groups_F")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[271:324]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with G
pdf(file = "boxplot_genes_of_interest_3groups_G")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[325:366]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with H
pdf(file = "boxplot_genes_of_interest_3groups_H")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[367:402]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with I
pdf(file = "boxplot_genes_of_interest_3groups_I")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[403:429]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with J & K
pdf(file = "boxplot_genes_of_interest_3groups_J_and_K")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[430:459]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with L
pdf(file = "boxplot_genes_of_interest_3groups_L")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[460:489]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with M
pdf(file = "boxplot_genes_of_interest_3groups_M1")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[490:528]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

pdf(file = "boxplot_genes_of_interest_3groups_M2")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[529:567]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with N
pdf(file = "boxplot_genes_of_interest_3groups_N")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[568:606]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with O
pdf(file = "boxplot_genes_of_interest_3groups_O")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[607:636]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with P
pdf(file = "boxplot_genes_of_interest_3groups_P1")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[637:684]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

pdf(file = "boxplot_genes_of_interest_3groups_P2")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[685:729]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with Q and R
pdf(file = "boxplot_genes_of_interest_3groups_Q_and_R1")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[730:768]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with S
pdf(file = "boxplot_genes_of_interest_3groups_S1")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[808:855]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

pdf(file = "boxplot_genes_of_interest_3groups_S2")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[856:900]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with T 
pdf(file = "boxplot_genes_of_interest_3groups_T1")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[901:957]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

pdf(file = "boxplot_genes_of_interest_3groups_T2")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[958:1011]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with U and V
pdf(file = "boxplot_genes_of_interest_3groups_U_and_V")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[1012:1041]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()

#genes starting with X, Y and Z
pdf(file = "boxplot_genes_of_interest_3groups_X_Y_and_Z")
par(las=2)
boxplot(t(thyroid.expr.matrix.new)[,order.vector[1042:1080]],
        cex.axis=0.6,
        col=c("red","blue","green"),
        main="Gene expression of thyroid-specific TRAs in thyroid cancer, \nGSE35570 (Handkiewicz-Junak et all., 2016)",
        ylab="log2(gene expression)",
)
abline(h=c(6:12),lty=3)
legend("topright",c("radexp","nonradexp","norm"),col=c("red","blue","green"),pch=15,bg="white")
dev.off()


