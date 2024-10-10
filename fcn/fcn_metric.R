library(dplyr)

metrics <- function(protein_i, 
                    time, 
                    tcc) {
  # remove NAs
  nonNA_idx <- !is.na(protein_i)
  time <- time[nonNA_idx]
  protein_i <- protein_i[nonNA_idx]
  n_pt <- length(time)
  
  # Calculate the slope
  slope_beta <- sum(protein_i * time) / sum(time^2)
  
  # Calculate the half life
  k_dp <- slope_beta - log(2)/tcc
  T_HalfLife <- log(2) / k_dp
  
  # Calculate R^2
  SS_tot <- sum((protein_i - mean(protein_i))^2)
  SS_res <- sum((protein_i - slope_beta * time)^2)
  R2 <- 1 - SS_res / SS_tot
  
  # Calculate the average leave-one-out cross validation sum of squares for the slope
  LeaveOneOut_slope <- c()
  for (i in 1:n_pt) {
    slope_beta_i <- sum(protein_i[-i] * time[-i]) / sum(time[-i]^2)
    LeaveOneOut_slope <- c(LeaveOneOut_slope, slope_beta_i)
  }
  CVSQ <- sum((slope_beta - LeaveOneOut_slope)^2) / n_pt
  
  return(list("slope_beta" = slope_beta,
              "T_HalfLife" = T_HalfLife, 
              "R2" = R2, 
              "CVSQ" = CVSQ,
              "n_pt" = n_pt))
}
