install.packages("nnet")
library("nnet")

subtype.matrix = rbind(cancer.groups, cancer.subtype.sig.genes.matrix)
subtype.matrix = t(subtype.matrix)

#defining subtype 1 as reference category
subtype.matrix[,1] = relevel(factor(subtype.matrix[,1]), ref = "subtype 1")

#regeression model
foldchange.matrix=subtype.matrix[,foldchange.vector]
foldchange.matrix= cbind(cancer.groups,foldchange.matrix)
subtype.model=multinom(cancer.groups~.,
                       data = as.data.frame(foldchange.matrix[c(1,2,5,6,15,16,17,18,19,20),]))

#welche Variablen sind wichtig (nächstes Mal)

#
predict.subtype.1 = predict(subtype.model, foldchange.matrix[c(1,2,5,6,15,16,17,18,19,20),])
predict.subtype.2 = predict(subtype.model, foldchange.matrix[c(3,4,7,8,9,10,11,12,13,14),])     

#dioganal matrix if all subtypes were identified correctly
tab1=table(predict.subtype.1, foldchange.matrix[c(1,2,5,6,15,16,17,18,19,20),1])
tab1 #is diagonal



                  