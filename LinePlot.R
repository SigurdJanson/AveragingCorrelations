
#' SimpleLinePlot
#' ...
SimpleLinePlot <- function( n, m, Height = 0.05 ) {
  require(RColorBrewer, quietly	= TRUE, warn.conflicts = FALSE)
  
  # N <- seq(10, 50, 10) # sample size
  D <- 3:10            # Data sets to correlate
  N <- seq(10, 50, 10)   # sample size
  R <- seq(0.00, 0.95, 0.05)  # Rho
  Palette <- brewer.pal(n = length(D)+1, name = "PuRd")[(length(D)+1):2]
  RangeY <- c(-Height, Height)
  for(d in D) {
    nChr <- as.character(n)
    dChr <- as.character(d)
    Color <- Palette[which(d==D)]
    Means <- subset(Descr, group2==nChr & group3==dChr & group4==m, select = mean)
    ScaleY <- max(abs(Means))
    if(d == D[1]) {
      plot(R, unlist(Means), ylim = RangeY, col = Color, type="l", lty=1,
           sub = n, xlab = "r", ylab = "rho - r", oma = c(1, 2, 2, 2))
      grid(nx = NA, ny = NULL, col = "lightgray")      
    }
    else
      lines(R, unlist(Means), col = Color, type="l", lty=1)
  }
  abline(h = 0, col = "gray60")
  text(0, y = min(RangeY), adj = 0, m)
}



#' LinePlot
#' Draws one or two line plots side by side to compare two correction 
#' methods.
#' @param m1,m2 A string that describes a correction method. m2 is optional.
#' @param Lines Shall lines illustrate different sample size or numbers 
#' of data sets?
#' @param Selected A string or number to select a factor level for
#' which lines will be drawn. If 'Lines' is a sampe size this argument 
#' selects a number of data sets and vice versa.
LinePlot <- function(m1, m2 = NULL, 
                     Lines = c("Sample Size", "Data Sets"), 
                     Selected = 50,
                     DataStruc = "Matrix") {
  require(ggplot2, quietly	= TRUE, warn.conflicts = FALSE)
  if(!is.null(m2)) require(gridExtra, quietly	= TRUE, warn.conflicts = FALSE)
  
  PlotThis <- function(d, m) {
    ThePlot <- ggplot(data = d, 
                      aes_string(x = "group1", y = paste0("mean.", Selected), 
                                 group = LineVar, color = LineVar, 
                                 shape = LineVar, linetype = LineVar)
    ) +
      expand_limits(y=DataRange) +
      scale_x_discrete(breaks=seq(0, 0.95, 0.10)) +
      scale_colour_manual(values = LineCol) + 
      scale_shape_manual(values = c(0:2, 5:6, 15:18)) +
      scale_linetype_manual(values = rep(c("solid", "dashed"), 5)) +
      theme(legend.position="bottom") +
      geom_line(aes_string(x = "group1", y = paste0("mean.", Selected), 
                           group = LineVar, color = LineVar, linetype = LineVar)
      ) +
      geom_point(aes_string(color = LineVar, shape = LineVar)) + 
      labs(title=paste0("Distribution of Bias (", m, ")"), 
           caption = paste0("Correction = ", m, " / ", SelectData, " = ", Selected), 
           x="Rho", y = "Rho - r", color = Lines, linetype = Lines, shape = Lines)
    
    return(ThePlot)
  }
  
  if(Lines == "Sample Size") {
    LineVar  <- "group2"
    ShapeVar <- "group3"
    SelectData <- "Data Sets"
  } else if(Lines == "Data Sets") {
    LineVar  <- "group3"
    ShapeVar <- "group2"
    SelectData <- "Sample Size"
  }
  else stop("Unknown dimension 'dim'")
  
  
  # load data and reshape into wide form
  SelectVar <- paste0("group", 1:4)
  SelectVar <- SelectVar[SelectVar != ShapeVar]
  load(paste0("./data/CoreyDunlapBurke_", DataStruc, "_", m1, "_Avg.Rda"))
  d1 <- reshape(Descr, v.names = "mean", timevar = ShapeVar, direction = "wide", 
                idvar = SelectVar, 
                drop  = c("item", "n", "vars", "sd", "min", "max", "se", "range"))
  DataRange <- range(d1[[paste0("mean.", Selected)]])
  if(Lines == "Data Sets") {
    # Correct the order of the factor levels
    d1[[LineVar]] <- factor(d1[[LineVar]], sort(as.integer(levels(d1$group3))))
  }
  
  if(!is.null(m2)) {
    load(paste0("./data/CoreyDunlapBurke_",  DataStruc, "_", m2, "_Avg.Rda"))
    d2 <- reshape(Descr, v.names = "mean", timevar = ShapeVar, direction = "wide", 
                  idvar = SelectVar, 
                  drop  = c("item", "n", "vars", "sd", "min", "max", "se", "range"))
    D2Range <- range(d2[[paste0("mean.", Selected)]])
    DataRange[1] <- min(DataRange[1], D2Range[1])
    DataRange[2] <- max(DataRange[2], D2Range[2])
    if(Lines == "Data Sets") {
      # Correct the order of the factor levels
      d2[[LineVar]] <- factor(d2[[LineVar]], sort(as.integer(levels(d2$group3))))
    }
  }
  
  # Plot settings
  theme_set(theme_gray(base_size = 10))
  LineCol <- c("#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8",
               "#253494", "#081d58", "#021025", "#000000")
  
  
  
  LP <- PlotThis(d1, m1)
  if(!is.null(m2)) {
    LP2 <- PlotThis(d2, m2)
    grid.arrange(LP, LP2, nrow = 1)
  } else {
    plot(LP)
  }
}

#LinePlot("None", Lines = "Data Sets")
#LinePlot("None", "Fisher", Lines = "Data Sets", Selected = 10)
