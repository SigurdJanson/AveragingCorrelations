#
# Analyses for Corey, Dunlap & Burke 
#

library(DataExplorer)
library(psych)
library("RColorBrewer")

# Only needed when data needs to be updated
M <- c("None", "Fisher", "Hotelling", "MinVar", "TrueK", "Precise")
for(m in M) {
  if(file.exists(file=paste("./data/CoreyDunlapBurke", m, "Raw.Rda", sep = "_"))) {
    load(file=paste("./data/CoreyDunlapBurke", m, "Raw.Rda", sep = "_"))
    summary(Result)
    
    Descr <- describeBy(Result[,"RDelta"], Result[,c("Rho", "SampleSize", "Samples", "Method")], 
                        skew = FALSE, digits=6, omit = TRUE, mat = TRUE )
    save(Descr, file = paste("./data/CoreyDunlapBurke", m, "Avg.Rda", sep = "_"))
  }
}
rm(M, m)


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

# Get aggregated data
# if(!exists("Descr")) {
#   Descr <- describeBy(Result[,"RDelta"], Result[,c("Rho", "SampleSize", "Samples", "Method")], 
#                       skew = FALSE, digits=6, omit = TRUE, mat = TRUE )
#   save(Descr, file = "./data/CoreyDunlapBurke.MinVarTruek.Avg.Rda")
# }



MakePlot <- function( n, m ) {
  # N <- seq(10, 50, 10) # sample size
  # D <- 3:10            # Data sets to correlate
  Palette <- brewer.pal(n = length(D)+1, name = "PuRd")[(length(D)+1):2]
  RangeY <- c(-0.005, 0.005)
  for(d in D) {
    nChr <- as.character(n)
    #Color <- Palette[which(n==N)]
    dChr <- as.character(d)
    Color <- Palette[which(d==D)]
    Means <- subset(Descr, group2==nChr & group3==dChr & group4==m, select = mean)
    if(d == D[1]) {
      plot(R, unlist(Means), ylim = RangeY, col = Color, type="l", lty=1,
           sub = n, xlab = "r", ylab = "rho - r", oma = c(1, 2, 2, 2))
      grid(nx = NA, ny = NULL, col = "lightgray")      
    }
    else
      lines(R, unlist(Means), col = Color, type="l", lty=1)
  }
  abline(h = 0, col = "gray60")
  text(0, y = min(RangeY), m)
}

old.par <- par(mfrow = c(2, 2), mar = c(4, 3, 0.5, 1)+0.1)
load("./data/CoreyDunlapBurke_None_Avg.Rda")
MakePlot(50, "None")
load("./data/CoreyDunlapBurke_Fisher_Avg.Rda")
MakePlot(50, "Fisher")
load("./data/CoreyDunlapBurke_MinVar_Avg.Rda")
MakePlot(50, "MinVar")
load("./data/CoreyDunlapBurke_Truek_Avg.Rda")
MakePlot(50, "TrueK")
par(old.par)
rm(old.par)
#rep(Palette, each=length(D), length.out=length(N))