# Functions for individual to population niche

# Function to get env values weighted by probability of usage returned by CTMM
# dat = vector of observed environmental values for an individual
# w = vector of weights of the observed environmental values
# x = numeric, the value to look up
getUDVal_wt <- function(x, dat, w){
  d <- density(dat, weights = w,
               from = 15,
               to = 45,
               bw = 0.5)
  
  approx <- approxfun(d$x, d$y)
  
  # Look up the probability associated with some value x.
  prob <- approx(x)
  
  return(prob)
}


# Function to calculate individual contribution of niche breadth and skewness
# x is a data frame that contains three parameters at the individual level: mean, variance and skewness
individual_contribution <- function(x, w = NULL) {
  # x is the individual parameter data, columns are mu and sigma (note that it is the sd not the variance),
  # but the returned population sigma2 is the variance
  n = nrow(x)
  if(!is.null(w)) n = w
  mu_i = as.numeric(x[,"mu"])
  sigma_i = as.numeric(x[,"sigma"])
  if("skew" %in% colnames(x)){
    skew_i = as.numeric(x[,"skew"])
  }
  mu_i = mu_i - mean(mu_i)
  mu = mean(mu_i)
  marginality_sigma2 = (mu_i^2 - mu^2)/n
  specialization_sigma2 = (sigma_i^2)/n
  sigma2 = sum(marginality_sigma2) + sum(specialization_sigma2)
  marginality_skew = (mu_i^3 - mu^3)/(n*sigma2^(3/2))
  specialization_skew = (3*mu_i*(sigma_i^2)-3*mu*sigma2)/(n*sigma2^(3/2))
  if("skew" %in% colnames(x)){
    indiv_skew = (skew_i*(sigma_i^3))/(n*sigma2^(3/2))
    specialization_skew = specialization_skew + indiv_skew
  }
  skew = sum(marginality_skew) + sum(specialization_skew)
  
  output = cbind(mu_pop = mu,marginality_sigma2 = marginality_sigma2,
                 specialization_sigma2 = specialization_sigma2,
                 sigma2_pop = sigma2,
                 marginality_skew = marginality_skew,
                 specialization_skew = specialization_skew,
                 skew_pop = skew)
  if("skew" %in% colnames(x)) output = cbind(output,indiv_skew)
  return(output)
}
