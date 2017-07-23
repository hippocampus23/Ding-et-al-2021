# TODO regularization

all_wls <- function(in_list) {
  # in_list should be list of c(x, y, vars)
  fits <- vapply(in_list,
                 FUN=function(r) wls(r[[1]], r[[2]], r[[3]]),
                 FUN.VALUE=c("b"=0, "stderr"=0, "t"=0, "pval"=0))
  return(fits)
}

wls <- function(x, y, vars) {
  # Returns value of x coefficient, p-val, t-dist, and 
  data <- data.frame(x, y, weights = (1. / vars))
  fit <- lm('y ~ x', data=data, weights=weights)
  out <- summary(fit)
  # Returns u_est, std_err, t, p_val
  return(out$coefficients[2,])
}
