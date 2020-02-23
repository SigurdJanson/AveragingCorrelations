# Create a violin plot
PlayViolin <- function(m1, m2, DataStruc = "Matrix") {
  require(ggplot2, quietly	= TRUE, warn.conflicts = FALSE)
  theme_set(theme_gray(base_size = 10))
  
  load(paste0("./data/CoreyDunlapBurke_", DataStruc, "_", m1, "_Avg.Rda"))
  d1 <- Descr
  load(paste0("./data/CoreyDunlapBurke_", DataStruc, "_", m2, "_Avg.Rda"))
  d2 <- Descr
  
  Description <- rbind( d1, d2 )
  ColMethod <- c(None = "#D16103", Fisher = "#FFDB6D", Hotelling = "#C4961A", 
                 MinVar = "#56B4E9", TrueK = "#4E84C4", Precise = "#293352" )
  ColViolin <- c("#999999", ColMethod[levels(Description$group4)])
  ggplot(data = Description, 
         mapping = aes(x = group4, 
                       y = mean,
                       color = ""
         )
  ) + 
    scale_color_manual(values=ColViolin) +
    scale_fill_manual(values=ColViolin) +
    geom_violin(
      trim = TRUE,
      draw_quantiles = c(0.05, 0.5, 0.95),
    ) + 
    geom_boxplot( width=0.05, aes(color = group4) ) + 
    theme(legend.position="none")+
    labs(title="Distribution of Bias", x="Correction Method", y = "Rho")
}
