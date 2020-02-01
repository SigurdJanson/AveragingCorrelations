# Analyses for Corey, Dunlap & Burke 

# Only needed when data needs to be updated
if(!exists("Result")) {
  if(file.exists("./data/CoreyDunlapBurke.Raw.Rda"))
    load("data.Rda")
  else
    source("./CoreyDunlapBurke.job.R")
}


library(DataExplorer)
library(psych)
library("RColorBrewer")

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
plot_histogram(Result$RDelta)
quantile(Result$RDelta, probs = seq(0, 1, 0.1))




# ANALYSIS ----

# Get aggregated data
if(!exists("Descr")) {
  Descr <- describeBy(Result[,"RDelta"], Result[,c("Rho", "SampleSize", "Samples", "Method")], 
                      skew = FALSE, digits=3, omit = TRUE, mat = TRUE )
  save(Descr, file = "./data/CoreyDunlapBurke.Avg.Rda")
}


MakePlot <- function( n, m ) {
  # N <- seq(10, 50, 10) # sample size
  # D <- 3:10            # Data sets to correlate
  Palette <- brewer.pal(n = length(D)+1, name = "PuRd")[(length(D)+1):2]
  #for(n in N) {
  for(d in D) {
    nChr <- as.character(n)
    #Color <- Palette[which(n==N)]
    dChr <- as.character(d)
    Color <- Palette[which(d==D)]
    Means <- subset(Descr, group2==nChr & group3==dChr & group4==m, select = mean)
    if(d == D[1])
      plot(R, unlist(Means), ylim = c(-0.01, 0.01), col = Color, type="l", lty=1,
           sub = n, xlab = "r", ylab = "rho - r", oma = c(1, 2, 2, 2))
    else
      lines(R, unlist(Means), col = Color, type="l", lty=1)
  }
  abline(h = 0, col = "gray60")
  abline(h = c(0.01, 0, -0.01), col = "gray60")
  #}
  #legend(0, 0.043, legend=N, col = Palette, cex = 0.8, box.col=Palette[5])
}

old.par <- par(mfrow = c(2, 1), mar = c(4, 3, 0.5, 1)+0.1)
#MakePlot(10, "Default")
#MakePlot(10, "Fisher")
# MakePlot(20, "Default")
# MakePlot(20, "Fisher")
# MakePlot(30, "Default")
# MakePlot(30, "Fisher")
# MakePlot(40, "Default")
# MakePlot(40, "Fisher")
 MakePlot(50, "Default")
 MakePlot(50, "Fisher")
par(old.par)
rm(old.par)
#rep(Palette, each=length(D), length.out=length(N))