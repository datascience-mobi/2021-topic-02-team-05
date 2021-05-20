#heatmap of all our genes of interest on all 14 chips
#with the pheatmap function
library(pheatmap)
setwd("/Users/Line/Documents/tra_project/thyroid_cancer_sessions/plots")
pheatmap(thyroid.expr.matrix.sub, 
         main = "gene expression of thyroid TRA's \nin PTC and healthy thyroid tissue (GSE35570) \n(Handkiewicz-Junak et all., 2016)",
         cluster_rows = F, 
         cluster_cols = F, 
         cellheight = 3, 
         cellwidth = 3,
         legend =TRUE,
         fontsize = 3,
         filename = "pheatmap_thyroid_TRA_expression.pdf",
         gaps_col = c(5,9)
        )













