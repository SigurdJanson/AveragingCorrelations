#' GenData
#' @param Supplied.Data A data set to bootstrap from. Each column represents a 
#' random variable with any desired random distribution. Correlations are 
#' determined between columns.
#' @param N.Factors 
#' @param Max.Trials 
#' @param Initial.Multiplier 
#' @param Seed Seed of the random number generator to replicate results
#' @details This implementation of GenData bootstraps each variable’s score distribution 
#' from a supplied data set. Users should modify this block of the program, as 
#' needed, to generate the desired distribution(s).
#' @note   #' The original comments tells us to replace the line /LINE Q/ above with 
#' one of these expressions:
#' - Distributions[,i] <- sort(rchisq(N, df = 2))
#' - Distributions[,1] <- sort(rnorm(N, 0, 1)) # Standard normal distribution
#' - Distributions[,2] <- sort(runif(N, 0, 1)) # Uniform distribution ranging from 0 - 1
#' - Distributions[,3] <- sort(rlnorm(N, 0, 1)) # Log-normal distribution, log scale M = 0, SD = 1
#' - Distributions[,4] <- sort(rexp(N, rate = 1)) # Exponential distribution with rate = 1
#' - Distributions[,5] <- sort(rpois(N, lambda = 4)) # Poisson distribution with lambda = 4
#' - Distributions[,6] <- sort(rbinom(N, 10, .25) # Binominal distribution, size = 10 and p = .25
#' - Distributions[,7] <- sort(rbinom(N, 2, .25) # Binary distribution with p = .25
#'
#' Because of an error in the code this is not optional but essential. In 
#' other words, the function cannot be used in this version. This version requires that I 
#' already have the data I want to create because it gets the correlation matrix AND the 
#' probability distrbitutions from 'Supplied.Data'.
#' 
#' This implementation of GenData calculates the target correlation matrix 
#' from a supplied data set. 
#' Alternatively, the user can modify the program to generate data with 
#' user-defined sample size, number of variables, and target correlation 
#' matrix by redefining the function as follows:
#'   GenData <- function(N, k, Target.Corr, N.Factors = 0, Max.Trials = 5, Initial.Multiplier = 1, seed = 0)
#' In this case, one would also remove the program lines that 
#' calculate N, k, and Target.Corr. To generate data in which variables 
#' are uncorrelated, one would remove the SsortT function from step 2
#' and terminate the program before step 3 begins by returning the 
#' Distributions object as the data set. 
#' @examples 
#' N <- 50 
#' # Use chi², normal and log-normal distribution
#' sampleData <- cbind(rchisq(N, df = 2), rnorm(N, 0, 1), rlnorm(N, 0, 1))
#' GenData(SampleData)
#' @author John Ruscio and Walter Kaczetow; refactored by Jan Seifert
#' @references Ruscio, J., & Kaczetow, W. (2008). Simulating Multivariate Nonnormal Data Using an Iterative Algorithm. 
#' Multivariate Behavioral Research, 43 (3), 355–381. 
GenData <- function(Supplied.Data, N.Factors = 0, Max.Trials = 5, Initial.Multiplier = 1, 
                    Seed = NA)
{
  # Initialize variables and (if applicable) set random number seed (step 1) ------------
  N <- dim(Supplied.Data)[1]         # Number of cases
  k <- dim(Supplied.Data)[2]         # Number of variables
  Data <- matrix(0, nrow = N, ncol = k) # Matrix to store the simulated data
  Distributions <- matrix(0, nrow = N, ncol = k) # Matrix to store each variable’s score distribution
  Iteration <- 0                     # Iteration counter
  Best.RMSR <- 1                     # Lowest RMSR correlation
  Trials.Without.Improvement <- 0    # Trial counter
  if (!is.na(Seed)) set.seed(Seed)   # If user specified a seed
  
  # Generate distribution for each variable (step 2) ------------------------------------
  for (i in 1:k)
    Distributions[,i] <- sort(sample(Supplied.Data[,i], N, replace = TRUE)) # /LINE Q/

  # All of the commands shown above draw random samples from specified population 
  # distributions. As an alternative, one can reproduce distributions without 
  # sampling error. For example, working with a supplied data set, one can replace 
  # the 2nd line in this block with:
  #    \code{Distributions[,i] <- Supplied.Data[,i]}
  # Alternatively, idealized distributions can be reproduced. For example, uniform 
  # quantiles can be created and used to generate data from common distributions:
  # 
  #    Uniform.Quantiles <- seq(from = 0, to = 1, length = (N + 2))[2:(N + 1)] # quantiles 0, 1 dropped
  #    Distributions[,1] <- qnorm(Uniform.Quantiles, 0, 1) # Standard normal distribution
  #    Distributions[,2] <- qunif(Uniform.Quantiles, 0, 1) # Uniform distribution ranging from 0 to 1
  #    Distributions[,3] <- qchisq(Uniform.Quantiles, df = 2) # Chi-square distribution with 2 df
  #
  # Note that when score distributions are generated from specified populations 
  # rather than bootstrapped from a supplied data set, the user must provide 
  # the target correlation matrix (see the next block). This is true regardless 
  # of whether the distributions incorporate sampling error.

  # Calculate and store a copy of the target correlation matrix (step 3) ----------------
  Target.Corr <- cor(Supplied.Data)
  Intermediate.Corr <- Target.Corr
  
  # If number of latent factors was not specified, determine it through parallel analysis (step 4) ---------
  if (N.Factors == 0) {
    Eigenvalues.Observed <- eigen(Intermediate.Corr)$values
    Eigenvalues.Random <- matrix(0, nrow = 100, ncol = k)
    Random.Data <- matrix(0, nrow = N, ncol = k)
    for (i in 1:100) {
      for (j in 1:k) {
        Random.Data[,j] <- sample(Distributions[,j], size = N, replace = TRUE)
      }
      Eigenvalues.Random[i,] <- eigen(cor(Random.Data))$values
    }
    # calculate mean eigenvalue for each factor
    Eigenvalues.Random <- apply(Eigenvalues.Random, 2, mean) 
    N.Factors <- max(1, sum(Eigenvalues.Observed > Eigenvalues.Random))
  }
  
  # Generate random normal data for shared and unique components, initialize factor loadings (steps 5, 6) --------
  Shared.Comp <- matrix(rnorm(N * N.Factors, 0, 1), nrow = N, ncol = N.Factors)
  Unique.Comp <- matrix(rnorm(N * k, 0, 1), nrow = N, ncol = k)
  Shared.Load <- matrix(0, nrow = k, ncol = N.Factors)
  Unique.Load <- matrix(0, nrow = k, ncol = 1)
  
  # Begin loop that ends when specified number of iterations pass without improvement in RMSR correlation --------
  while (Trials.Without.Improvement < Max.Trials) {
    Iteration <- Iteration + 1
    # Calculate factor loadings and apply to reproduce desired correlations (steps 7, 8) ----
    Fact.Anal <- Factor.Analysis(Intermediate.Corr, Corr.Matrix = TRUE, N.Factors = N.Factors)
    if (N.Factors == 1) Shared.Load[,1] <- Fact.Anal$loadings
    else Shared.Load <- Fact.Anal$loadings
    Shared.Load[Shared.Load > 1] <- 1
    Shared.Load[Shared.Load < -1] <- -1
    if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
    Shared.Load.sq <- Shared.Load * Shared.Load
    for (i in 1:k)
      if (sum(Shared.Load.sq[i,]) < 1) Unique.Load[i,1] <- (1 - sum(Shared.Load.sq[i,]))
    else Unique.Load[i,1] <- 0
    Unique.Load <- sqrt(Unique.Load)
    for (i in 1:k) {
      # %*% operator = matrix multiplication, t() = transpose
      Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
    }

    # Replace normal with nonnormal distributions (step 9) ------------------------------
    for (i in 1:k) {
      Data <- Data[sort.list(Data[,i]),]
      Data[,i] <- Distributions[,i]
    }
    
    # Calculate RMSR correlation, compare to lowest value, take appropriate action (steps 10, 11, 12)
    Reproduced.Corr <- cor(Data)
    Residual.Corr <- Target.Corr - Reproduced.Corr
    RMSR <- sqrt(sum(Residual.Corr[lower.tri(Residual.Corr)] * 
                       Residual.Corr[lower.tri(Residual.Corr)]) /
                   (.5 * (k * k - k)))
    if (RMSR < Best.RMSR) {
      Best.RMSR <- RMSR
      Best.Corr <- Intermediate.Corr
      Best.Res <- Residual.Corr
      Intermediate.Corr <- Intermediate.Corr + Initial.Multiplier * Residual.Corr
      Trials.Without.Improvement <- 0
    } else {
      Trials.Without.Improvement <- Trials.Without.Improvement + 1
      Current.Multiplier <- Initial.Multiplier * .5 ^ Trials.Without.Improvement
      Intermediate.Corr <- Best.Corr + Current.Multiplier * Best.Res
    }
  } # end of the while loop
  
  # Construct the data set with the lowest RMSR correlation (step 13) -------------------
  Fact.Anal <- Factor.Analysis(Best.Corr, Corr.Matrix = TRUE, N.Factors = N.Factors)
  if (N.Factors == 1) Shared.Load[,1] <- Fact.Anal$loadings
  else Shared.Load <- Fact.Anal$loadings
  Shared.Load[Shared.Load > 1] <- 1
  Shared.Load[Shared.Load < -1] <- -1
  if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
  Shared.Load.sq <- Shared.Load * Shared.Load
  for (i in 1:k)
    if (sum(Shared.Load.sq[i,]) < 1) Unique.Load[i,1] <- (1 - sum(Shared.Load.sq[i,]))
  else Unique.Load[i,1] <- 0
  Unique.Load <- sqrt(Unique.Load)
  for (i in 1:k)
    Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
  Data <- apply(Data, 2, scale) # standardizes each variable in the matrix
  for (i in 1:k) {
    Data <- Data[sort.list(Data[,i]),]
    Data[,i] <- Distributions[,i]
  }
  
  # Report the results and return the simulated data set (step 14) ----------------------
  Iteration <- Iteration - Max.Trials
  cat("\nN =", N, ", k =", k, ",", Iteration, "Iterations,",
      N.Factors, "Factors, RMSR r =", round(Best.RMSR,3), "\n")
  return(Data)
}




#' Factor.Analysis
#' Home made factor analysis as support function for 'GenData'.
#' @param Data Raw data or correlation matrix. Will be coerced to a matrix. Will be 
#' interpreted as correlation matrix only if 'Corr.Matrix' is TRUE.
#' @param CorrMatrix Specifies if Data is raw data or a correlation matrix.
#' @param MaxIter Maximum number of iterations for convergence.
#' @param NFactors Number of factors to extract, default is \code{ncol(Data)}
#' @details If 'NFactors' equals the number of variables (i.e. columns of 'Data') 
#' then this function returns the factors with eigenvalues larger than 1.
#' Otherwise it returns the 'NFactors' factors with the highest eigenvalues.
#' @author John Ruscio and Walter Kaczetow; refactored by Jan Seifert
#' @references Ruscio, J., & Kaczetow, W. (2008). Simulating Multivariate Nonnormal Data Using an Iterative Algorithm. 
#' Multivariate Behavioral Research, 43 (3), 355–381. 
Factor.Analysis <- function(Data, CorrMatrix = FALSE, MaxIter = 50, NFactors = NA) {
  Data <- as.matrix(Data)
  k <- ncol(Data) # dim(Data)[2]
  if (is.na(NFactors)) NFactors <- k
  if (!CorrMatrix)
    Cor.Matrix <- cor(Data)
  else 
    Cor.Matrix <- Data

  Criterion <- 0.001
  Old.H2 <- rep(99, k)
  H2 <- rep(0, k)
  Change <- 1
  Iter <- 0
  Factor.Loadings <- matrix(nrow = k, ncol = NFactors)

  while ((Change >= Criterion) & (Iter < MaxIter)) {
    Iter <- Iter + 1
    Eig <- eigen(Cor.Matrix)
    L <- sqrt(Eig$values[1:NFactors])
    for (i in 1:NFactors)
      Factor.Loadings[,i] <- Eig$vectors[,i] * L[i]
    for (i in 1:k)
      H2[i] <- sum(Factor.Loadings[i,] * Factor.Loadings[i,])
    Change <- max(abs(Old.H2 - H2))
    Old.H2 <- H2
    diag(Cor.Matrix) <- H2
  }
  
  if (NFactors == k) NFactors <- sum(Eig$values > 1)
  return(list(loadings = Factor.Loadings[,1:NFactors], factors = NFactors))
}
