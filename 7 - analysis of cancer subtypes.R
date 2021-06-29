#seperating the chips according to the cancer subtype
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


