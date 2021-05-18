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

##This is the matrix you can calculate with##
