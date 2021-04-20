theta_hat <- function(beta, X1,X2,Phi){
  logit <- beta[1] + beta[2] *X1 + beta[3] *X2 + Phi
  theta <- exp(logit) / (1 + exp(logit))
  return(theta)
}
