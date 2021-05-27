#### Data Science Project (Breast Cancer) ####
# University of Heidelberg


####Preparation####

## Install Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

## Install packages
BiocManager::install("affy")
BiocManager::install("vsn")
BiocManager::install("AnnotationDbi")

## Set working directory
wd <- "~/Users/macbook/Desktop/Bioinfo/Data Science Project"

# Load packages

library(affy)
library(vsn)
library(AnnotationDbi)

## Brainarray packages are located in .../Bioinfo/Data Science Project/

install.packages("~/Desktop/Bioinfo/Data Science Project/hgu133plus2hsenstprobe_25.0.0.tar.gz", repos = NULL, type = "source")
install.packages("~/Desktop/Bioinfo/Data Science Project/hgu133plus2hsenstcdf_25.0.0.tar.gz", repos = NULL, type = "source")

library(hgu133plus2hsenstcdf)
library(hgu133plus2hsenstprobe)

library(hexbin)

library(ggplot2)

## Read in CEL files, Breast Cancer

setwd(paste0(wd, "/rawdata-breast"))
cels <- list.files(pattern = "CEL")
dataRaw <- ReadAffy(filenames = cels, verbose = TRUE)
dataRaw@cdfName <- "HGU133Plus2_Hs_ENST"

new <- substr(cels, 1, nchar(cels)-4)
colnames(exprs(dataRaw)) <- new
rownames(dataRaw@phenoData) <- new
rownames(dataRaw@protocolData) <- new

save(dataRaw, file = paste0(wd, "/dataRaw.RData"))

####Microarray control and  Matrix fomation ####

## Single chip quality control
setwd(paste0(wd, "/plots"))
for (chipNo in seq_along(new)) { # for each chip
  pdf(file = paste0("QCs/QC_", new[chipNo], ".pdf")) # save in QCs/
  image(dataRaw[, chipNo], col = rainbow(100, start = 0, end = 0.75)[100:1]) # blue -> green -> red
  dev.off()
}

## Vsn normalisation

normdata <- vsnrma(dataRaw)
save(normdata, file = paste0(wd, "/norm.RData"))

## Mean SD plot
meanSdPlot(normdata)
ggsave("Mean_sd_plot.pdf")

## Extract expression values
dataexprs <- exprs(normdata)
head(dataexprs)
dim(dataexprs)
# 95721    10

## Ensemble IDs
dataexprs <- dataexprs[grepl("ENST", rownames(dataexprs)), ]
head(dataexprs)
dim(dataexprs)
# 95659    10

## Remove .x_at suffix
rownames(dataexprs) <- unlist(lapply(strsplit(rownames(dataexprs), split = "\\."), "[", 1))
head(dataexprs)

## Load ensemble IDs and gene symbols
setwd(wd)
ensembltxt <- read.csv("ensembl.103.txt", sep = ",")
head(ensembltxt)

transcriptIDs <- as.character(ensembltxt[, "Transcript.stable.ID"])
geneSymbol <- as.character(ensembltxt[, "HGNC.symbol"])
names(geneSymbol) <- transcriptIDs
head(geneSymbol)

## replace ensemble IDs with HGNC symbol
chipIDs <- rownames(dataexprs)
noMatchIDs <- chipIDs[!chipIDs %in% transcriptIDs] # transcript IDs that are not in ensemble_103.txt
length(noMatchIDs)
# 112 
newChipIDs <- chipIDs[chipIDs %in% transcriptIDs]
dataexprs <- dataexprs[newChipIDs, ]
dim(dataexprs)
# 95547    10

symbol <- geneSymbol[newChipIDs]
rownames(dataexprs) <- as.character(symbol)
head(dataexprs)

###Plots####

## Expression values before and after normalization

#Boxplot
#before normalization
setwd(paste0(wd, "/plots"))
pdf(file = 'boxplot_breast_cancer_rawdata.pdf')
par(mar=c(9,5,3,5))
boxplot(dataRaw,
        col=rainbow(14), 
        cex.axis=0.5, 
        ylab="intensity", 
        las=2, 
        main="Gene expression in Breast Cancer\nbefore normalization ")
dev.off()

#after normalization
setwd(paste0(wd, "/plots"))
pdf(file = 'boxplot_breastcancer_normalized.pdf')
par(mar=c(9,5,3,5))
boxplot(exprs(normdata), 
        col=rainbow(14), 
        cex.axis=0.5, 
        ylab="log intensity", 
        las=2, 
        main="Gene expression in Breast Cancer \nafter vsnrma normalization")
dev.off()

#density function 
#raw data
setwd(paste0(wd, "/plots"))
pdf(file = 'hist_breastcancer_rawdata.pdf')
hist(dataRaw, col=rainbow(14), main="Density function of log Intensity of papillary thyroid cancer (GSE35570) \nbefore normalization (Handkiewicz-Junak et all., 2016)")
dev.off()

#normalized data
setwd(paste0(wd, "/plots"))
pdf(file = 'hist_breastcancer_normalised.pdf')
plot(density(exprs(normdata)[,1]), type="n", xlab="log Intensity", ylim = c(0,1), main="Density function of log Intensity of papillary thyroid cancer (GSE35570) \nafter vsnrma normalization (Handkiewicz-Junak et all., 2016)")
for (i in 1:ncol(exprs(normdata))) {
  lines(density(exprs(normdata)[,i]), col=rainbow(14)[i])
}
dev.off()

#RNA degradation plot 
setwd(paste0(wd, "/plots"))
pdf(file = 'rnadeg_human_breast caner_rawdata.pdf')
rnadeg.raw = AffyRNAdeg(dataRaw)
plotAffyRNAdeg(rnadeg.raw, col=rainbow(14))
title(sub="human breast cancer rawdata ")
dev.off()

pdf(file = 'rnadeg_human_breastcancer_rawdata_shifted.pdf')
plotAffyRNAdeg(rnadeg.raw, col=rainbow(14), transform="shift.only")
title(sub="human breast cancer raw data")
dev.off()


####TRA List####

#set working directory

setwd(wd)

#load data

#data set 1

a=read.csv(file="tra.2014.mouse.5x.table.tsv",sep="\t")

tiss=a[,11]

ind=which(tiss=="thyroid")

TRA.symbol1=a[,3]

thyroid.TRA1=TRA.symbol1[ind]

str(a)

# data set 2 (no thyroid genes)

b=read.csv(file="tra.2014.mouse.4301.5x.table.tsv",sep="\t")

tiss2=b[,11]

ind2=which(tiss2=="Thyroid")

TRA.symbol2=b[,3]

thyroid.TRA2=TRA.symbol2[ind2]

#data set 3

c=read.csv(file="tra.2017.human.gtex.5x.table.tsv",sep="\t")

tiss3=c[,10]

ind3=which(tiss3=="Thyroid")

TRA.symbol3=c[,3]

thyroid.TRA3=TRA.symbol3[ind3]

#data set 4

d=read.csv(file="tra.2014.human.roth.5x.table.tsv",sep="\t")

tiss4=d[,11]

ind4=which(tiss4=="thyroid_gland")

TRA.symbol4=d[,3]

thyroid.TRA4=TRA.symbol4[ind4]

#data set 5

e=read.csv(file="tra.2014.human.5x.table.tsv",sep="\t")

tiss5=e[,11]

ind5=which(tiss5=="Thyroid")

TRA.symbol5=e[,3]

thyroid.TRA5=TRA.symbol5[ind5]


#unite the data 

thyroid.TRA=c(thyroid.TRA1,thyroid.TRA2,thyroid.TRA3,thyroid.TRA4,thyroid.TRA5)

thyroid.Trafinal=na.omit(thyroid.TRA)

y=as.data.frame(thyroid.Trafinal)

final_genes=y[!apply(y == "", 1, all), ]  

####Combining array data with TRA genes####

row.ind=which(final_genes%in%symbol)
dataexprs.sub=dataexprs[row.ind,]

#Box plot

order.vector=sort(colnames(t(dataexprs.sub)),index.return=T)$ix
setwd(paste0(wd, "/plots"))
pdf(file = 'boxplot_breastcancer_TRA.pdf')

boxplot(t(dataexprs.sub)[,order.vector],cex.axis=0.6,col=c("red","blue","yellow"),main="Gene expression of thyroid-specific TRAs in breast cancer",ylab="log2(gene expression)")
abline(h=c(6:12),col="grey",lty=3)

dev.off()
