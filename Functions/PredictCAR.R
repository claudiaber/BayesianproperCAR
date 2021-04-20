predictCAR <- function(alpha, tau, b_0, b_x1, b_x2, W, x1, x2, seed_mvnorm){
D <- diag(apply(W,1,sum))
set.seed(seed_mvnorm)
system.time({phi <- MASS::mvrnorm(mu=rep(0,nrow(W)), Sigma=solve(tau *(D - alpha *W)))})
logit <- b_0 + b_x1*x1 + b_x2*x2 + phi
prob <- exp(logit) / (1 + exp(logit))
return(c(prob, phi))
}

predictCARSTAN <- function(extfitSTAN,  W, x1, x2, seed_mvnorm){
  predictCAR(b_0=extfitSTAN[1], b_x1=extfitSTAN[2], b_x2=extfitSTAN[3],
    alpha = extfitSTAN[4], tau = extfitSTAN[5], W=W, x1=x1, x2=x2, seed_mvnorm = seed_mvnorm)
  
}
