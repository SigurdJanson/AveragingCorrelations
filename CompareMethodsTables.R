CompareMethodsTable <- function(d1, d2) {
  require(psych, quietly	= TRUE, warn.conflicts = FALSE)
  
  Select <- c("mean", "sd", "median", "mad", "min", "max", "range")
  Description <- rbind( describe(d1$mean), describe(d2$mean) )
  Description <- Description[,Select] # drop unwanted stats
  # replace "median absolute deviation" with "zero absolute deviation
  Description[1, "mad"] <- mean(abs(d1$mean))
  Description[2, "mad"] <- mean(abs(d2$mean))
  colnames(Description)[colnames(Description) == "mad"] <- "zad"
  # determine number of values closer to zero
  ctz1 <- sum( abs(d1$mean) < abs(d2$mean) )
  ctz2 <- sum( abs(d1$mean) > abs(d2$mean) )
  Description <- cbind(Description, CtZ = c(ctz1, ctz2))
  # print
  rownames(Description) <- c(levels(d1$group4), levels(d2$group4))
  
  knitr::kable(Description,
               caption = paste(rownames(Description)),
               digits = 3)
}

