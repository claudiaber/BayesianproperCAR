rm(list=ls())
library(stringr)
library(rstan)
library(dplyr)
source("./001_STAN_simulated/Functions/tableSummaryTex.R")
# load("./001_STAN_simulated/models_nophi_20191122_1658.rdata")
#### Set up a square lattice region ------
x.easting <- 1:50
x.northing <- 1:50
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)

#### set up distance and neighbourhood (W, based on sharing a common border) matrices ----
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1] <-1 	
isSymmetric(distance)
# is.positive.definite(distance)
# is.positive.definite(0.4 * exp(-0.1 * distance))


#### Generatethe covariates and response data
set.seed(125)
# x1 <- rnorm(K)+1
x1 <- runif(n = K,min = 0,max = 1)
# x1_test <- c(0,0.2, 0.5, 0.8,1)
set.seed(124)
# x2 <- rnorm(K)
x2 <- runif(n = K, min = 0,max = 1)
# x2_test <- c(0,0.2, 0.5, 0.8,1)
# theta <- rnorm(K, sd=0.05)
# phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=0.4 * exp(-0.1 * distance))

logit <- 1 + 4 *x1 - 6 *x2
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(1,K)
set.seed(126)
Y <- rbinom(n=K, size=trials, prob=prob)

# logit_test <- 1 + 4 *x1_test - 6 *x2_test
# prob_test <- exp(logit_test) / (1 + exp(logit_test))
# trials_test <- rep(1,length(logit_test))
# set.seed(126)
# Y_test <- rbinom(n=length(logit_test), size=trials_test, prob=prob_test)

#test base glm logit --------------
set.seed(127)
fitglm <- glm(Y ~ x1+x2, data = cbind.data.frame(Y, x1,x2), family = binomial(link="logit"))
# summary(fitglm)
DT <- cbind.data.frame(
  summary(fitglm)$coefficients[,1:3],
  Pvalue=format(summary(fitglm)$coefficients[,4], digits=2))
tableSummaryTex(summarytab = DT, 
                round = 2,
                rownames = c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$"),
                label = "tab:tableNoPhi_glm",
                caption = "Coefficients estimated by the standard R glm algorithm 
                               on non-autocorrelated data",
                filename = "./001_STAN_simulated/tableNoPhi_glm.tex")
rm(DT)
# 
# DTNoPhi_glm <- data.frame(cbind(
#   round(summary(fitglm)$coefficients[,1:3],2),
#   pvalue=format(summary(fitglm)$coefficients[,4], digits=2)),
#   row.names = c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$"))
# TableNoPhi_glm <- knitr::kable(DTNoPhi_glm,
#                                format = "latex", caption = "\\label{tab:tableNoPhi_glm}
#                                Coefficients estimated by the standard R glm algorithm 
#                                on non-autocorrelated data", escape=FALSE)
# fileTableNoPhi_glm <-  file("./001_STAN_simulated/tableNoPhi_glm.tex")
# writeLines(text = TableNoPhi_glm, fileTableNoPhi_glm)
# close(fileTableNoPhi_glm)
# rm(DTNoPhi_glm, TableNoPhi_glm, fileTableNoPhi_glm)
# write.csv(summary(fitglm)$coefficients, file = "./001_STAN_simulated/NoPhi_glm.csv")

# Call:
#   glm(formula = Y ~ x1 + x2, family = binomial(link = "logit"))
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.9253  -0.6540   0.1604   0.6519   2.8343  
# 
# Coefficients:
#                 Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)   1.0885     0.1289   8.446   <2e-16 ***
#   x1            4.2449     0.2264  18.749   <2e-16 ***
#   x2           -6.2648     0.2543 -24.637   <2e-16 ***
#   ---
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 3462.6  on 2499  degrees of freedom
# Residual deviance: 2133.7  on 2497  degrees of freedom
# AIC: 2139.7
# 
# Number of Fisher Scoring iterations: 5

# pred <- predict(fitglm, newdata = list(x1=x1_test, x2= x2_test), type = "response")
#        1         2         3         4         5 
# 0.7480946 0.6647435 0.5196271 0.3711218 0.2826452 
# table(cbind.data.frame(Pred=pred>0.5,Y_test))

# test stan logit --------------
datalist <- list(N = K, Y = Y, x1 = x1, x2 = x2)
# rstan::stan(file = "./simulated data.stan",data = datalist, iter = 1000, chains = 4,verbose = TRUE )

set.seed(127) #not used by STAN, should be set as an argument of the stan() function
fitstan <- stan(file = "./001_STAN_simulated/Logit_training_2pred.stan",
                data = datalist, iter = 2500, chains = 4,verbose = TRUE,
                seed = 127 )

tableSummaryTex(summarytab = summary(fitstan, pars = fitstan@sim$pars_oi, 
                                     probs = c(0.025,#  0.25, 0.5, 0.75,
                                               0.975))$summary, 
                round = 2,
                rownames = c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$", "lp"),
                label = "tab:tableNoPhi_glmSTAN",
                caption = "Coefficients estimated by a GLM algorithm written in STAN
                                   on non-autocorrelated data",
                filename = "./001_STAN_simulated/tableNoPhi_glmSTAN.tex")
rm(datalist)


pdf(file = "./001_STAN_simulated/NoPhi_stanlogit.pdf", width = 15)
# rstan:::print.stanfit(fitstan, )
plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
n_kept <- fitstan@sim$n_save - fitstan@sim$warmup2
title(paste(
  paste("Inference for Stan model: ", fitstan@model_name, ".\n", sep = ""),
  paste(fitstan@sim$chains, " chains, each with iter=", fitstan@sim$iter,
        "; warmup=", fitstan@sim$warmup, "; thin=", fitstan@sim$thin, "; \n",
        "post-warmup draws per chain=", n_kept[1], ", ", "total post-warmup draws=",
        sum(n_kept), ".\n\n", sep = "")),cex.main=1.5)
# gridExtra::grid.table(round(summary(fitstan, pars = fitstan@sim$pars_oi, probs = c(0.025, 0.25, 0.5, 
#                                    0.75, 0.975))$summary, 4))
# Inference for Stan model: Logit_training_2pred.
# 4 chains, each with iter=1000; warmup=500; thin=1; 
# post-warmup draws per chain=500, total post-warmup draws=2000.
# 
# mean se_mean   sd     2.5%      25%      50%      75%    97.5% n_eff Rhat
# alpha     1.09    0.00 0.13     0.84     1.00     1.09     1.18     1.36   788    1
# b_x1      4.27    0.01 0.23     3.81     4.12     4.27     4.42     4.72  1059    1
# b_x2     -6.28    0.01 0.25    -6.77    -6.45    -6.28    -6.11    -5.80   776    1
# lp__  -1068.38    0.05 1.28 -1071.79 -1068.89 -1068.04 -1067.47 -1066.96   700    1
# 
# Samples were drawn using NUTS(diag_e) at Wed Feb 13 09:30:33 2019.
# For each parameter, n_eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor on split chains (at 
#                                                                   convergence, Rhat=1).
plot(fitstan, pars = c("alpha", "b_x1", "b_x2"))
rstan::traceplot(fitstan, pars = c("alpha", "b_x1", "b_x2"))
acf(fitstan@sim$samples[[1]]$alpha, lag.max = 1000)
dev.off()
rm(n_kept)
# extstan <- extract(fitstan)
# alpha_post <- extstan$alpha
# bx1_post <- extstan$b_x1
# bx2_post <- extstan$b_x2


# # semplificare!!!! -----------
# 
# predstan <- stan(file = "./001_STAN_simulated/Logit_prediction_2pred.stan",
#                  data = list(x1_test = x1_test,x2_test = x2_test, N = length(x1_test),
#                              N_samples = length(alpha_post),
#                              alpha = alpha_post,
#                              b_x1 = bx1_post,
#                              b_x2= bx2_post),
#                  chains = 1, iter = 1,warmup = 0,
#                  algorithm = "Fixed_param", verbose = TRUE)
# 
# # Extract and format output
# ext_pred <- extract(predstan)
# out_mat <- apply(ext_pred[["y_test"]][1,,],2,mean)
# table(cbind.data.frame(Pred=out_mat>0.5,Y_test))

# test CAR BAYES --------- 
formula <- Y ~ x1 + x2
set.seed(127)
model <- CARBayes::S.CARbym(formula=formula, family="binomial", trials=trials,
                            W=W, burnin=2000, n.sample=10000)

tableSummaryTex(summarytab = model$summary.results, 
                round = 2,
                rownames =  c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$", 
                              "$\\tau^{2}$", "$\\sigma^{2}$"),
                label = "tab:tableNoPhi_CARBayes",
                caption = "Coefficients estimated by the CARBayes package in R
                                   on non-autocorrelated data",
                filename = "./001_STAN_simulated/tableNoPhi_CARBayes.tex")

pdf(file = "./001_STAN_simulated/NoPhi_CARBAYES.pdf", width = 15)
# rstan:::print.stanfit(fitstan, )
gridExtra::grid.table(round(model$summary.results, 4))
# Posterior quantities and DIC
# 
#               Median    2.5%   97.5% n.sample % accept n.effective Geweke.diag
# (Intercept)  1.1134  0.8615  1.3784     8000     55.9      2971.7         6.9
# x1           4.3302  3.8680  4.7871     8000     55.9      1223.8         5.9
# x2          -6.3860 -6.9244 -5.8912     8000     55.9       800.4        -8.1
# tau2         0.0159  0.0048  0.1668     8000    100.0         5.5         4.3
# sigma2       0.0787  0.0209  0.1813     8000    100.0         5.7         3.7
# 
# DIC =  2139.71       p.d =  38.17321       LMPL =  -1070.01  
dev.off()

increasesamples <- 50000
set.seed(127)
model2 <- CARBayes::S.CARbym(formula=formula, family="binomial", trials=trials,
                             W=W, burnin=2000, n.sample=increasesamples)
# model$summary.results


tableSummaryTex(summarytab = model2$summary.results, 
                round = 2,
                rownames =  c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$", 
                              "$\\tau^{2}$", "$\\sigma^{2}$"),
                label = "tab:tableNoPhi_CARBayes2",
                caption = paste("Coefficients estimated by the CARBayes package in R
                                on non-autocorrelated data with ",
                                increasesamples ," samples"),
                filename = "./001_STAN_simulated/tableNoPhi_CARBayes2.tex")

#### Toy example for checking
# model <- CARBayes::S.CARbym(formula=formula, family="binomial", trials=trials,
# W=W, burnin=20, n.sample=50)

# test STAN CAR -------------
# shuold extract a 
datalist <- list(N = K, Y = Y, x1 = x1, x2 = x2, W=W, W_n = sum(W)/2)
# rstan::stan(file = "./simulated data.stan",data = datalist, iter = 1000, chains = 4,verbose = TRUE )

# require(coda)
set.seed(127)
fitcarstan <- stan(file = "./001_STAN_simulated/CARlogit_training_2pred.stan",
                   data = datalist, iter = 2500, chains = 4,warmup = floor(1000/3),
                   save_dso=TRUE,
                   verbose = TRUE ,
                   seed = 127)

tableSummaryTex(summarytab = summary(fitcarstan,
                                     pars = fitcarstan@sim$pars_oi,
                                     probs = c(0.025, #0.25, 0.5, 
                                               # 0.75,
                                               0.975))$summary[c(1:4,
                                                                 2502:2506),], 
                round = 2,
                rownames = c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$", 
                             "$\\phi_{1}$", "$\\phi_{2499}$", "$\\phi_{2500}$",
                             "$\\tau$", "$\\alpha$", "$lp$"),
                label = "tab:tableNoPhi_CARSTAN",
                caption = "Coefficients estimated by the sparse CAR implementation in STAN 
                           on non-autocorrelated data. Only few $\\phi$s are shown",
                filename = "./001_STAN_simulated/tableNoPhi_CARSTAN.tex")


pdf(file = "./001_STAN_simulated/NoPhi_CARstan.pdf", width = 15)
# rstan:::print.stanfit(fitstan, )
plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
n_kept <- fitcarstan@sim$n_save - fitcarstan@sim$warmup2
title(paste(
  paste("Inference for Stan model: ", fitcarstan@model_name, ".\n", sep = ""),
  paste(fitcarstan@sim$chains, " chains, each with iter=", fitcarstan@sim$iter,
        "; warmup=", fitcarstan@sim$warmup, "; thin=", fitcarstan@sim$thin, "; \n",
        "post-warmup draws per chain=", n_kept[1], ", ", "total post-warmup draws=",
        sum(n_kept), ".\n\n", sep = "")),cex.main=1.5)
# 
# text(0.4,0.5,adj=c(0,0),lab=
#        paste("Inference for Stan model: ", fitcarstan@model_name, ".\n", sep = ""))
#      # cat(x@sim$chains, " chains, each with iter=", x@sim$iter,
#      #     "; warmup=", x@sim$warmup, "; thin=", x@sim$thin, "; \n",
#      #     "post-warmup draws per chain=", n_kept[1], ", ", "total post-warmup draws=",
#      #     sum(n_kept), ".\n\n", sep = ""))
gridExtra::grid.table(round(summary(fitcarstan,
                                    pars = fitcarstan@sim$pars_oi,
                                    probs = c(0.025, 0.25, 0.5, 
                                              0.75, 
                                              0.975))$summary, 4)[1:10,])
# print(fitcarstan)
plot(fitcarstan, pars = c("alpha", "b_x1", "b_x2"))
traceplot(fitcarstan, pars = c("alpha", "b_x1", "b_x2"))
traceplot(fitcarstan, pars = c("phi[1]", "phi[2]", "phi[3]"))
acf(fitcarstan@sim$samples[[1]]$alpha, lag.max = 1000)
dev.off()

# layout(t(c(1,2)))
Theta_CARbayes <- cbind(prob,
                           apply(model$samples$fitted,2,mean),
                           t(apply(model$samples$fitted,2,quantile, 
                                   probs=c(0.025,0.975))))
matplot(Theta_CARbayes[order(prob),], 
        ylab="theta", xlab="Locations ordered by theta",type = "l", 
        main="Estimated with CARBayes")

Theta_CARbayesinc <- cbind(prob,
                           apply(model2$samples$fitted,2,mean),
                           t(apply(model2$samples$fitted,2,quantile, 
                                   probs=c(0.025,0.975))))
matplot(Theta_CARbayesinc[order(prob),], 
        ylab="theta", xlab="Locations ordered by theta",type = "l", 
        main="Estimated with CARBayes")

extcarstan <- extract(fitcarstan)
Logit_samples_CARStan <-  matrix(extcarstan$alpha,
                                 nrow = length(extcarstan$alpha),
                                 ncol = K, byrow = FALSE) +
  extcarstan$b_x1 %*% t(x1) + 
  extcarstan$b_x2 %*% t(x2)+ 
  extcarstan$phi
Theta_samples_CARStan <-  exp(Logit_samples_CARStan) / (1 + exp(Logit_samples_CARStan))
Theta_CARStan <- cbind(prob,
                       apply(Theta_samples_CARStan,2,mean),
                       t(apply(Theta_samples_CARStan,2,quantile, 
                               probs=c(0.025,0.975))))
matplot(Theta_CARStan[order(prob),], 
        ylab="theta", xlab="Locations ordered by theta",type = "l", 
        main="Estimated with CAR sparse in STAN")

# dev.off()


Coverage_theta_CARbayesinc <- sum(Theta_CARbayesinc[, "prob"]>=Theta_CARbayesinc[, "2.5%"] &
                                    Theta_CARbayesinc[, "prob"]<=Theta_CARbayesinc[, "97.5%"] )/K

Coverage_theta_CARStan <- sum(Theta_CARStan[, "prob"]>=Theta_CARStan[, "2.5%"] &
                                Theta_CARStan[, "prob"]<=Theta_CARStan[, "97.5%"] )/K
Phi_CARbayesinc <- cbind(apply(model2$samples$psi,2,mean),
      t(apply(model2$samples$psi,2,quantile, 
              probs=c(0.025,0.975))))

Coverage_phi_CARbayesinc <- sum(0>=Phi_CARbayesinc[, "2.5%"] &
                                  0<=Phi_CARbayesinc[, "97.5%"] )/K

extcarstan <- extract(fitcarstan)
# # str(extcarstan$phi)
# 
interv_phi <- t(apply(extcarstan$phi, 2,
                       FUN = function(x) quantile(x, probs = c(0.025, 0.975))))
Coverage_phi_CARStan <-  sum(0>=interv_phi[, "2.5%"] &
                               0<=interv_phi[, "97.5%"] )/K
# any(interv_phi[,1]>0)
# any(interv_phi[,2]<0)
# 
# interv_phi <- t(apply(extcarstan$phi, 2,
#                       FUN = function(x) quantile(x, probs = c(0.10, 0.90))))
# any(interv_phi[,1]>0)
# any(interv_phi[,2]<0)
# 
# 
# 
# interv_phi <- t(apply(extcarstan$phi, 2,
#                       FUN = function(x) quantile(x, probs = c(0.15, 0.85))))
# any(interv_phi[,1]>0)
# any(interv_phi[,2]<0)
# 
# 
# interv_phi <- t(apply(extcarstan$phi, 2,
#                       FUN = function(x) quantile(x, probs = c(0.20, 0.80))))
# any(interv_phi[,1]>0)
# any(interv_phi[,2]<0)
# 
# sum(interv_phi[,1]>0)
# sum(interv_phi[,2]<0)
# 
# # asd <- coda::mcmc(extcarstan$phi[,1:10])
# # coda::codamenu()
# 
save(list=ls(), file=paste0("./001_STAN_simulated/models_nophi_",
                            format(Sys.time(),format = "%Y%m%d_%H%M")
                            ,".rdata"))
