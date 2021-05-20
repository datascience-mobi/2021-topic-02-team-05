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
thyroid.expr.matrix.sub <- thyroid.expr.matrix.sub[,c(1,5,6,7,8,2,3,4,9,10,11,12,13,14)]

#save matrix as csv
setwd("/Users/Line/Documents/tra_project/thyroid_cancer_sessions/tables")
write.csv(thyroid.expr.matrix.sub, "thyroid_gene_expression_TRA")

#boxplot
par(las=2)
boxplot(t(thyroid.expr.matrix.sub),col=rainbow(length(rownames(thyroid.expr.matrix.sub))),main="gene expression of thyroid-specific genes in thyroid cancer",cex.axis=0.8)
#sort alphabetically
setwd("/Users/Line/Documents/tra_project/thyroid_cancer_sessions/rda")
pdf(file="boxplot_thyroid_cancer_genes_of_interest.pdf")

order.vector=sort(colnames(t(thyroid.expr.matrix.sub)),index.return=T)$ix
boxplot(t(thyroid.expr.matrix.sub)[,order.vector],col=rainbow(length(rownames(thyroid.expr.matrix.sub))),main="gene expression of thyroid-specific genes in papillary thyroid cancer \n(GSE35570) (Handkiewicz-Junak et all., 2016)",cex.axis=0.8)

dev.off()

#create three matrices (radiation_exposed/not_exposed/healthy)
radexp.matrix=thyroid.expr.matrix.sub[,c(1,5,6,7,8)]
radnonexp.matrix=thyroid.expr.matrix.sub[,c(2,3,4,9)]
healthy.matrix=thyroid.expr.matrix.sub[,c(10,11,12,13,14)]

