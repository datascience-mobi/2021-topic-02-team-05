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
setwd("/Users/Line/Documents/tra_project/thyroid_cancer_sessions/plots")
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
install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot", force=TRUE)
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

















