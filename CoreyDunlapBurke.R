#
#

library(DataExplorer)
library("RColorBrewer")

#' AggregateAvgCorSim
#' Load the complete simulation data and average it for all combinations of 
#' Rho, Sample sizes, and Data sets
#' @param M One or several of these:
#'  None, Fisher, Hotelling, Hotelling2, MinVar, TrueK, Precise
#' @param DataStruc Either 'Matrix' or 'Indie'
#' @note Only needed when data needs to be updated
AggregateAvgCorSim <- function(M, DataStruc) {
  require(psych)
  for(m in M) {
    if(file.exists(file=paste("./data/CoreyDunlapBurke", DataStruc, m, "Raw.Rda", sep = "_"))) {
      load(file=paste("./data/CoreyDunlapBurke", DataStruc, m, "Raw.Rda", sep = "_"))
      summary(Result)
      
      Descr <- describeBy(Result[,"RDelta"], Result[,c("Rho", "SampleSize", "Samples", "Method")], 
                          skew = FALSE, digits=6, omit = TRUE, mat = TRUE )
      save(Descr, file = paste("./data/CoreyDunlapBurke", DataStruc, m, "Avg.Rda", sep = "_"))
    }
  }
  invisible(Result)
}
AggregateAvgCorSim("Hotelling2", "Matrix")
rm(AggregateAvgCorSim)



# PLAUSIBILITY ----
# Check, if the data is plausible
# see https://www.littlemissdata.com/blog/simple-eda
# see https://cran.r-project.org/web/packages/DataExplorer/vignettes/dataexplorer-intro.html

introduce(Result)
describe(Result, omit = TRUE)
plot_bar(Result)

# Compare methods
#bi.bars(Result,"age","gender",ylab="Age",main="Age by males and females")

# Show total distribution of Rho - RObs
old.par <- par(mfrow = c(1, 2))
Method <- "MinVar"
hist(Result$RDelta[Result$Method == Method], 
     main = Method, xlab = "Rho - z", breaks = 24) #plot_histogram
Method <- "Precise"
hist(Result$RDelta[Result$Method == "Hypergeo"], 
     main = Method, xlab = "Rho - z", breaks = 24) #plot_histogram
par(old.par)
rm(old.par, Method)


quantile(Result$RDelta, probs = seq(0, 1, 0.1))




# ANALYSIS ----


#old.par <- par(mfrow = c(2, 2), mar = c(4, 3, 0.5, 1)+0.1)
load("./data/CoreyDunlapBurke_Matrix_None_Avg.Rda")
LinePlot("None", Lines = "Sample Size", Selected = 10)
load("./data/CoreyDunlapBurke_Fisher_Avg.Rda")
MakePlot(50, "Fisher")
load("./data/CoreyDunlapBurke_MinVar_Avg.Rda")
MakePlot(50, "MinVar")
load("./data/CoreyDunlapBurke_Truek_Avg.Rda")
MakePlot(50, "TrueK")
#par(old.par)
#rm(old.par)
