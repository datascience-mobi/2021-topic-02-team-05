#posthoc to find out which genes are significantly expressed in all 
#three subtypes to use those for the regression model
subtype.matrix = rbind(cancer.groups, cancer.subtype.sig.genes.matrix)
subtype.matrix = t(subtype.matrix)
modelgenes.sub = c()

library(multcompView)
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
install.packages("nnet")
library("nnet")

#defining subtype 1 as reference category
subtype.matrix[,1]= relevel(factor(subtype.matrix[,1]), ref = "subtype 1")

#regeression model (training)
model.matrix.sub=subtype.matrix[,modelgenes.sub]
model.matrix.sub= cbind(cancer.groups,model.matrix.sub)

#Versuch Expressionswerte numerisch zu machen (bis jetzt ohne Erfolg)
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
