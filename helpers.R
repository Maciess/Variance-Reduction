generate_matrix <- function(n) {
  # Create an n x n matrix
  matrix <- matrix(0, nrow = n, ncol = n)
  
  # Fill the matrix using the formula
  for (i in 1:n) {
    for (j in 1:n) {
      matrix[i, j] <- min(i, j) / n
    }
  }
  return(matrix)
}


GetBrownianSample <- function(n)
{
  #return Brownian sample that with lenght n where each point correspond 
  covMatrix <- generate_matrix(n)
  mean <- replicate(n, 0)
  mvrnorm(1, mean, covMatrix)
}

GetGeometricBrownianSample <- function(S, r, sigma, n)
{
  brownianSample <- GetBrownianSample(n)
  time <- (1:n) / n
  mu <- r - (sigma^2) / 2
  S*exp(mu * time + sigma * brownianSample)
}

