
# function to calculate IFR from Intercept and Slope
IFR_calc <- function(Intercept, Slope, Age_groups, Pop_str){
  sum(exp(Age_groups * Slope + Intercept) * Pop_str/sum(Pop_str))
}

# function to estimate Intercept from Slope and IFR
Int_calc <- function(IFR, Slope, Age_groups, Pop_str){
  # browser()
  # Over a range of intercept values
  Intercept_vals <- seq(-41,-3,by = 0.00001)
  varyInt <- abs(sapply(X = Intercept_vals, FUN = IFR_calc, Slope = Slope, Age_groups = Age_groups, Pop_str = Pop_str) - IFR)
  Intercept_vals[which(min(varyInt)==varyInt)]
}
