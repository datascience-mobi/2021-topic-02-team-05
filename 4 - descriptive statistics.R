#heatmap of all our genes of interest on all 14 chips
#with the pheatmap function
library(pheatmap)
pheatmap(thyroid.expr.matrix.sub, 
         cluster_rows = F, 
         cluster_cols = F, 
         cellheight = 3, 
         cellwidth = 3,
         legend =TRUE,
         fontsize = 3,
         filename = "pheatmap_thyroid_TRA_expression.pdf",
         gaps_col = c(5,9)
        )













