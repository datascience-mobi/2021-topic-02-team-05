thyroid.expr.matrix.new=rbind(thyroid.expr.matrix.sub,thyroid.expr.matrix.sub,thyroid.expr.matrix.sub)

thyroid.expr.matrix.new[c(1:360),11:30]=NA
thyroid.expr.matrix.new[c(361:720),c(1:10,21:30)]=NA
thyroid.expr.matrix.new[c(721:1080),1:20]=NA

order.vector=sort(colnames(t(thyroid.expr.matrix.new)),index.return=T)$ix #gives position of genes in alphabetical order

setwd("/Users/Line/Documents/tra_project/with_30_chips/plots")
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






