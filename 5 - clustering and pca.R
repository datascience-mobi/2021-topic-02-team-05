###PCA (with 30 chips)

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

#k-means of pca with sich patiets only
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
dev.off()
#3 clusters are a better fit

















