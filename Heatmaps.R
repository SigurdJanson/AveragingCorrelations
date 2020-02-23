
CompareByHeat <- function(d1, d2) {
  require(gplots, quietly	= TRUE, warn.conflicts = FALSE)
  
  DataLabels <- c(levels(d1$group4), levels(d2$group4))
  # compute average deviation for each Rho
  M1 <- aggregate(d1[, 8], list(Descr$group2, Descr$group3), median)
  MH <- aggregate(d2[, 8], list(Descr$group2, Descr$group3), median)
  
  MAll <- merge(M1, MH, by = c(1, 2), suffixes = paste0(".", DataLabels))
  Label1 <- paste0("x.", DataLabels[1])
  Label2 <- paste0("x.", DataLabels[2])
  data   <- matrix(unlist(abs(MAll[Label1]) - abs(MAll[Label2])), nrow=8)
  
  Colors <- c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1','#35978f','#01665e')
  
  # Default Heatmap
  heatmap.2(data, scale="column", Colv = NA, Rowv = NA, dendrogram = "none",
            xlab = "Sample size", labRow = 3:10,
            ylab = "Data sets", labCol = seq(10, 50, 10),
            main = paste("Comparison", DataLabels[1], "-", DataLabels[2]),
            col = Colors, tracecol = "black",
            key.title = "Colors & Histogram", 
            key.xlab = paste0("|", DataLabels[1], "|-|", DataLabels[2], "|"))
}
