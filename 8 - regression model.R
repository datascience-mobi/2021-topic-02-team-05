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
install.packages("nnet")
library("nnet")

#defining subtype 1 as reference category
subtype.matrix[,1]= relevel(factor(subtype.matrix[,1]), ref = "subtype 1")

#regeression model
model.matrix.sub=subtype.matrix[,modelgenes.sub]
model.matrix.sub= cbind(cancer.groups,model.matrix.sub)

#Versuch Expressionswerte numerisch zu machen (bis jetzt ohne Erfolg)
model.matrix.sub =as.data.frame(model.matrix.sub)
model.matrix.sub[,2:12]=numeric(model.matrix.sub[,2:12])
model.matrix.sub[,2:12]=as.matrix(apply(model.matrix.sub[,2:12],2,as.numeric))
model.matrix.sub =cbind(cancer.groups,as.numeric(model.matrix.sub[,2:12]))

subtype.model=multinom(cancer.groups~.,
                       data = as.data.frame(model.matrix.sub[c(1,2,5,6,15,16,17,18,19,20),]))


#welche Variablen sind wichtig (nächstes Mal)

#
predict.subtype.1 = predict(subtype.model, model.matrix.sub[c(1,2,5,6,15,16,17,18,19,20),])
predict.subtype.2 = predict(subtype.model, model.matrix.sub[c(3,4,7,8,9,10,11,12,13,14),])  #Error   

#dioganal matrix if all subtypes were identified correctly
tab1=table(predict.subtype.1, model.matrix.sub[c(1,2,5,6,15,16,17,18,19,20),1])
tab1 #is diagonal



levels(as.factor(model.matrix.sub[c(1,2,5,6,15,16,17,18,19,20),]))

class(model.matrix.sub[,2])
class(model.matrix.sub[3,2])
