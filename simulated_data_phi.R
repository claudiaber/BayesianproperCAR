rm(list=ls())
library(MASS)
library(stringr)
library(rstan)
library(dplyr)

source("./001_STAN_simulated/Functions/tableSummaryTex.R")
source("./001_STAN_simulated/Functions/theta_hat.R")
# load("./001_STAN_simulated/models_phi_20191123_1630.rdata")
# library(shinystan)
# launch_shinystan(fitcarstan)

#### Set up a square lattice region ------
griddim <- 50
x.easting <- 1:griddim
x.northing <- 1:griddim
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)

#### set up distance and neighbourhood (W, based on sharing a common border) matrices ----
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1| distance==(griddim-1)] <-1 	
isSymmetric(distance)
# is.positive.definite(distance)
# is.positive.definite(0.4 * exp(-0.1 * distance))


#### Generate the covariates and response data
set.seed(125)
# x1 <- rnorm(K)+1
x1 <- runif(n = K,min = 0,max = 1)
x1_test <- c(0,0.2, 0.5, 0.8,1)
set.seed(124)
# x2 <- rnorm(K)
x2 <- runif(n = K, min = 0,max = 1)
x2_test <- c(0,0.2, 0.5, 0.8,1)
# theta <- rnorm(K, sd=0.05)
alpha <- 0.8# 1/median(apply(W,1,sum))
tau <- 0.1
D <- diag(apply(W,1,sum))
tau *(D - alpha *W)[1:10, 1:10]
round(solve(tau *(D - alpha *W))[1:10, 1:10],3)
set.seed(123)
# undebug(mvrnorm)
system.time({phi <- mvrnorm(mu=rep(0,K), Sigma=solve(tau *(D - alpha *W)))})
# user  system elapsed 
# 43.79    0.11   44.48 
esempio <- 1946
Grid[c(esempio, which(W[esempio, ]==1)),]
phi[esempio]
mean(phi[which(W[esempio,]==1)])


logit <- 1 + 4 *x1 - 6 *x2 + phi
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
summary(fitglm)
DT <- cbind.data.frame(
  summary(fitglm)$coefficients[,1:3],
  Pvalue=format(summary(fitglm)$coefficients[,4], digits=2))
tableSummaryTex(summarytab = DT, 
                round = 2,
                rownames = c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$"),
                label = "tab:tablePhi_glm",
                caption = "Coefficients estimated by the standard R glm algorithm 
                               on spatially autocorrelated data",
                filename = "./001_STAN_simulated/tablePhi_glm.tex")
rm(DT)

pdf(file = "./001_STAN_simulated/Phi_glm.pdf", width = 15)
plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
gridExtra::grid.table(summary(fitglm)$coefficients)
dev.off()

# pred <- predict(fitglm, newdata = list(x1=x1_test, x2= x2_test), type = "response")
# 
# table(cbind.data.frame(Pred=pred>0.5,Y_test))

# test stan logit --------------
datalist <- list(N = K, Y = Y, x1 = x1, x2 = x2)

set.seed(127)
fitstan <- stan(file = "./001_STAN_simulated/Logit_training_2pred.stan",
                data = datalist, iter = 1000, chains = 4,verbose = TRUE,
                seed = 127)
print(fitstan)
fitstan_sum <- summary(fitstan, 
                       pars = fitstan@sim$pars_oi,
                       probs = c(0.025,# 0.25, 0.5, 
                                 # 0.75, 
                                 0.975))$summary
tableSummaryTex(summarytab = fitstan_sum, 
                round = 2,
                rownames = c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$", "lp"),
                label = "tab:tablePhi_glmSTAN",
                caption = "Coefficients estimated by a GLM algorithm written in STAN 
                             on spatially autocorrelated data",
                filename = "./001_STAN_simulated/tablePhi_glmSTAN.tex")

pdf(file = "./001_STAN_simulated/Phi_stanlogit.pdf", width = 15)
# rstan:::print.stanfit(fitstan, )
plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
n_kept <- fitstan@sim$n_save - fitstan@sim$warmup2
title(paste(
  paste("Inference for Stan model: ", fitstan@model_name, ".\n", sep = ""),
  paste(fitstan@sim$chains, " chains, each with iter=", fitstan@sim$iter,
        "; warmup=", fitstan@sim$warmup, "; thin=", fitstan@sim$thin, "; \n",
        "post-warmup draws per chain=", n_kept[1], ", ", "total post-warmup draws=",
        sum(n_kept), ".\n\n", sep = "")),cex.main=1.5)
gridExtra::grid.table(round(summary(fitstan, pars = fitstan@sim$pars_oi, probs = c(0.025, 0.25, 0.5, 
                                                                                   0.75, 0.975))$summary, 4))

plot(fitstan, pars = c("alpha", "b_x1", "b_x2"))
rstan::traceplot(fitstan, pars = c("alpha", "b_x1", "b_x2"))
acf(fitstan@sim$samples[[1]]$alpha, lag.max = 1000)
dev.off()

extcarstanglm <- extract(fitstan)
png(filename = "./001_STAN_simulated/MatplotTheta_glm_STANOrder.png", width = 880, height = 660)
layout(t(c(1,2)))
Theta_glm <-theta_hat(beta = fitglm$coefficients,
                      X1 = x1, X2=x2,
                      Phi=0)
plot(prob[order(prob)], ylab="theta", xlab="Locations ordered by theta",
     main="Estimated with standard R glm")
lines(Theta_glm[order(prob)], col="red")

Logit_samples_glmStan <-  matrix(extcarstanglm$alpha,
                                 nrow = length(extcarstanglm$alpha),
                                 ncol = K, byrow = FALSE) +
  extcarstanglm$b_x1 %*% t(x1) + 
  extcarstanglm$b_x2 %*% t(x2)#+ 
  # extcarstan$phi
Theta_samples_glmStan <-  exp(Logit_samples_glmStan) / (1 + exp(Logit_samples_glmStan))
Theta_glmStan <- cbind(prob,
                       apply(Theta_samples_glmStan,2,mean),
                       t(apply(Theta_samples_glmStan,2,quantile, 
                               probs=c(0.025,0.975))))
matplot(Theta_glmStan[order(prob),], 
        ylab="theta", xlab="Locations ordered by theta",type = "l", 
        main="Estimated with glm in STAN")
# Theta_glmstan <-theta_hat(beta = fitstan_sum[1:3,"mean"],
#                       X1 = x1, X2=x2,
#                       Phi=0)
# plot(prob[order(prob)], ylab="mean theta", xlab="Locations ordered by theta")
# lines(Theta_glmstan[order(prob)], col="red")

dev.off()
Coverage_theta_glmStan <- sum(Theta_glmStan[, "prob"]>=Theta_glmStan[, "2.5%"] &
  Theta_glmStan[, "prob"]<=Theta_glmStan[, "97.5%"] )/K
rm(Logit_samples_glmStan,Theta_samples_glmStan,Theta_glmStan)




# semplificare!!!! -----------
# alpha_post <- extcarstanglm$alpha
# bx1_post <- extcarstanglm$b_x1
# bx2_post <- extcarstanglm$b_x2
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
set.seed(127)# set.seed(123)# set.seed(127)
model <- CARBayes::S.CARbym(formula=formula, family="binomial", 
                                data=as.data.frame(cbind(Y, x1, x2)),
                                trials=trials,
                            W=W, burnin=2000, n.sample=10000)
print(model)
pdf(file = "./001_STAN_simulated/Phi_CARBAYES.pdf", width = 15)
gridExtra::grid.table(round(model$summary.results, 4))
dev.off()

tableSummaryTex(summarytab = model$summary.results, 
                round = 2,
                rownames =  c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$", 
                              "$\\tau^{2}$", "$\\sigma^{2}$"),
                label = "tab:tablePhi_CARBayes",
                caption = "Coefficients estimated by the CARBayes package in R
                                   on spatially autocorrelated data",
                filename = "./001_STAN_simulated/tablePhi_CARBayes.tex")

# model$modelfit

# prova <- t(apply(model$samples$psi, 2, FUN = function(x) quantile(x, probs = c(0.025, 0.975))))
# matplot(cbind(prova, phi), type="l")
# matplot(cbind(prova, phi)[order(phi),], type="l")
# 
# interv_phi <- t(apply(model$samples$psi, 2,
#                       FUN = function(x) quantile(x, probs = c(0.025, 0.975))))
# sum(interv_phi[,1]>phi)
# sum(interv_phi[,2]<phi)


increasesamples <- 50000
set.seed(127)# set.seed(12323436)#set.seed(123)# set.seed(127) #set.seed(123) #set.seed(12323436)# set.seed(127)
# modelincrease <- CARBayes::S.CARbym(formula=formula, family="binomial", trials=trials,
                            # W=W, burnin=2000, n.sample=increasesamples, thin = 5)
# print(modelincrease)
# modelincrease_old <-modelincrease
modelincrease <- CARBayes::S.CARbym(formula=formula, family="binomial", trials=trials,
                                    W=W, burnin=10000, n.sample=increasesamples, thin = 5)
# modelincrease <- modelincrease2

tableSummaryTex(summarytab = cbind.data.frame(Model=c("First","Increased"),
                                              rbind(model$modelfit,modelincrease$modelfit)), 
                round = 2,
                rownames = NULL,
                label = "tab:modelfitCARBayes",
                caption = "Model Fit Criteria for CARBayes estimations",
                filename = "./001_STAN_simulated/tablePhi_CARBayes_diag.tex")

tableSummaryTex(summarytab = modelincrease$summary.results, 
                round = 2,
                rownames =  c("$\\beta_{0}$", "$\\beta_{1}$", "$\\beta_{2}$", 
                              "$\\tau^{2}$", "$\\sigma^{2}$"),
                label = "tab:tablePhi_CARBayes2",
                caption = paste("Coefficients estimated by the CARBayes package in R
                                on spatially autocorrelated data with ",
                                increasesamples ," samples"),
                filename = "./001_STAN_simulated/tablePhi_CARBayes2.tex")

pdf(file = "./001_STAN_simulated/Phi_CARBAYESincrease.pdf", width = 15)
gridExtra::grid.table(round(modelincrease$summary.results, 4))
dev.off()
# model$samples$psi
prova <- t(apply(model$samples$psi, 2, FUN = function(x) quantile(x, probs = c(0.025, 0.975))))
provaincrease <- t(apply(modelincrease$samples$psi, 2, FUN = function(x) quantile(x, probs = c(0.025, 0.975))))
png(filename = "./001_STAN_simulated/MatplotPhi_CARBayesOrder.png", width = 880, height = 660)
layout(mat = t(c(1,2)))
# matplot(cbind(prova, phi), type="l")
# matplot(cbind(provaincrease, phi), type="l")

modelinterval <- cbind(order=1:2500,cbind(prova, phi)[order(phi),])
colnames(modelinterval)<- c("order","p025", "p975", "phi")
matplot(modelinterval[, 2:4], type="l", ylim = c(-7,7),
        # main="Ordered phi(green) vs estimated 2.5%-97.5% interval", 
        ylab="Phi", xlab="Locations ordered by phi")
abline(h = 0,lty="dashed", col="grey")

modelintervalincrease <- cbind(order=1:2500,cbind(provaincrease, phi)[order(phi),])
colnames(modelintervalincrease)<- c("order","p025", "p975", "phi")
matplot(modelintervalincrease[,2:4], type="l", ylim = c(-7,7), 
        # main="Ordered phi(green) vs new estimated 2.5%-97.5% interval",
        ylab="Phi", xlab="Locations ordered by phi")
abline(h = 0, lty="dashed", col="grey")
dev.off()

# 
# layout(t(c(1,2)))
# Theta_CARbayes <-theta_hat(beta = apply(model$samples$beta, 2, median),
#           X1 = x1, X2=x2,
#           Phi=apply(model$samples$psi, 2, median))
# plot(prob[order(prob)])
# lines(Theta_CARbayes[order(prob)], col="red")
# 
# Theta_CARbayesinc <-theta_hat(beta = apply(modelincrease$samples$beta, 2, median),
#                            X1 = x1, X2=x2,
#                            Phi=apply(modelincrease$samples$psi, 2, median))
# plot(prob[order(prob)])
# lines(Theta_CARbayesinc[order(prob)], col="red")



# NOT used -----------
# TablePhi_CARBayes_phiorderl <- knitr::kable(round(summary(lm(p10 ~ order,
#                                                     data = as.data.frame(modelinterval)))$coefficients,
#                                          4),
#                                    # col.names = colnames(summary(fitglm)$coefficients),
#                                    format = "latex", 
#                                    caption = paste("\\label{tab:tablePhi_CARBayes_phiorderl}
#                                                    Linear regression on the $\\psi$'s 10th percentile
#                                                     as computed by the CARBayes package in R
#                                                    on spatially autocorrelated data")) # with ",
#                                                    # increasesamples ," samples"))
# fileTablePhi_CARBayes_phiorderl <-  file("./001_STAN_simulated/tablePhi_CARBayes_phiorderl.tex")
# writeLines(text = TablePhi_CARBayes_phiorderl, fileTablePhi_CARBayes_phiorderl)
# close(fileTablePhi_CARBayes_phiorderl)
# 
# 
# # summary(lm(p90 ~ order, data = as.data.frame(modelinterval)))$coefficients
# TablePhi_CARBayes_phiorderu <- knitr::kable(round(summary(lm(p90 ~ order,
#                                                              data = as.data.frame(modelinterval)))$coefficients,
#                                                   4),
#                                             # col.names = colnames(summary(fitglm)$coefficients),
#                                             format = "latex", 
#                                             caption = paste("\\label{tab:tablePhi_CARBayes_phiorderu}
#                                                             Linear regression on the $\\psi$'s 90th percentile
#                                                             as computed by the CARBayes package in R
#                                                             on spatially autocorrelated data")) # with ",
# # increasesamples ," samples"))
# fileTablePhi_CARBayes_phiorderu <-  file("./001_STAN_simulated/tablePhi_CARBayes_phiorderu.tex")
# writeLines(text = TablePhi_CARBayes_phiorderu, fileTablePhi_CARBayes_phiorderu)
# close(fileTablePhi_CARBayes_phiorderu)
# 
# 
# # lm(p10 ~ order, data = as.data.frame(modelintervalincrease))
# TablePhi_CARBayesincrease_phiorderl <- knitr::kable(round(summary(lm(p10 ~ order,
#                                                              data = as.data.frame(modelintervalincrease)))$coefficients,
#                                                   4),
#                                             # col.names = colnames(summary(fitglm)$coefficients),
#                                             format = "latex", 
#                                             caption = paste("\\label{tab:tablePhi_CARBayesincrease_phiorderl}
#                                                             Linear regression on the $\\psi$'s 10th percentile
#                                                             as computed by the CARBayes package in R
#                                                             on spatially autocorrelated data with ",
#                                                             increasesamples ," samples"))
# fileTablePhi_CARBayesincrease_phiorderl <-  file("./001_STAN_simulated/tablePhi_CARBayesincrease_phiorderl.tex")
# writeLines(text = TablePhi_CARBayesincrease_phiorderl, fileTablePhi_CARBayesincrease_phiorderl)
# close(fileTablePhi_CARBayesincrease_phiorderl)
# 
# # lm(p90 ~ order, data = as.data.frame(modelintervalincrease))
# TablePhi_CARBayesincrease_phiorderu <- knitr::kable(round(summary(lm(p90 ~ order,
#                                                                      data = as.data.frame(modelintervalincrease)))$coefficients,
#                                                           4),
#                                                     # col.names = colnames(summary(fitglm)$coefficients),
#                                                     format = "latex", 
#                                                     caption = paste("\\label{tab:tablePhi_CARBayesincrease_phiorderu}
#                                                                     Linear regression on the $\\psi$'s 90th percentile
#                                                                     as computed by the CARBayes package in R
#                                                                     on spatially autocorrelated data with ",
#                                                                     increasesamples ," samples"))
# fileTablePhi_CARBayesincrease_phiorderu <-  file("./001_STAN_simulated/tablePhi_CARBayesincrease_phiorderu.tex")
# writeLines(text = TablePhi_CARBayesincrease_phiorderu, fileTablePhi_CARBayesincrease_phiorderu)
# close(fileTablePhi_CARBayesincrease_phiorderu)

## acf ------
acf(modelincrease$samples$psi[,220], lag.max =2000)
acf(modelincrease$samples$beta[,1], lag.max=2000)
acf(modelincrease$samples$beta[,2], lag.max=2000)
acf(modelincrease$samples$beta[,3], lag.max=2000)
acf(modelincrease$samples$tau2, lag.max=5000)
acf(model$samples$tau2, lag.max=5000)
acf(modelincrease$samples$sigma2, lag.max=2000)

## traceplot -----
plot(modelincrease$samples$psi[,220])
plot(modelincrease$samples$beta[,1])
plot(modelincrease$samples$beta[,2])
plot(modelincrease$samples$beta[,3])
plot(modelincrease$samples$sigma2)
plot(modelincrease$samples$tau2)
plot(model$samples$tau2)

# check if the chains are the same untill a certain point ----
which(modelincrease$samples$sigma2 == model$samples$sigma2[5000])
# # test CARBayes 200 000 -----------
# increasesamples2 <- 200000
# set.seed(127)
# system.time({modelincrease2 <- CARBayes::S.CARbym(formula=formula, family="binomial", trials=trials,
#                                                   W=W, burnin=10000, n.sample=increasesamples2, thin = 5)})
# acf(modelincrease2$samples$psi[,150], lag.max =2000)
# acf(modelincrease2$samples$beta[,1], lag.max=2000)
# acf(modelincrease2$samples$beta[,2], lag.max=2000)
# acf(modelincrease2$samples$beta[,3], lag.max=2000)
# acf(modelincrease2$samples$tau2, lag.max=2000)
# acf(modelincrease2$samples$sigma2, lag.max=2000)
# # Finished in  1368.5 seconds.
# print(modelincrease2)
# plot(modelincrease2$samples$sigma2)

# test STAN CAR -------------
# shuold extract a 
datalist <- list(N = K, Y = Y, x1 = x1, x2 = x2, W=W, W_n = sum(W)/2)
# rstan::stan(file = "./simulated data.stan",data = datalist, iter = 1000, chains = 4,verbose = TRUE )
require(rstan)
require(dplyr)
# require(coda)
set.seed(127)
fitcarstan <- stan(file = "./001_STAN_simulated/CARlogit_training_2pred.stan",
                   data = datalist, iter = 2500, chains = 4,warmup = 1000,# floor(1000/3),
                   save_dso=TRUE,
                   verbose = TRUE,
                   seed = 127 )
# Warning messages:
#   1: There were 24 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 2: There were 4 chains where the estimated Bayesian Fraction of Missing Information was low. See
# http://mc-stan.org/misc/warnings.html#bfmi-low 
# 3: Examine the pairs() plot to diagnose sampling problems
# print(fitcarstan, pars = c('beta', 'tau', 'alpha'))
# plot(fitcarstan, pars = c("alpha", "b_x1", "b_x2"))
# print(fitcarstan, pars = c("alpha", "b_x1", "b_x2"))
# Inference for Stan model: CARlogit_training_2pred.
# 4 chains, each with iter=1000; warmup=333; thin=1; 
# post-warmup draws per chain=667, total post-warmup draws=2668.
# 
# mean se_mean     sd     2.5%      25%      50%      75%    97.5% n_eff Rhat
# alpha         0.99    0.08   0.25     0.63     0.82     0.95     1.11     1.62    10 1.35
# b_x1          3.42    0.27   0.68     2.63     2.98     3.21     3.63     5.37     6 1.81
# b_x2         -5.40    0.43   1.05    -8.52    -5.70    -5.04    -4.72    -4.30     6 1.94
# phi[1]       -0.18    0.02   1.35    -2.72    -0.95    -0.20     0.58     2.48  3011 1.00
# phi[2]       -0.65    0.16   1.32    -3.81    -1.32    -0.53     0.16     1.53    68 1.05
# Samples were drawn using NUTS(diag_e) at Tue Feb 12 17:38:09 2019.
# For each parameter, n_eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor on split chains (at 
#                                                                   convergence, Rhat=1).
pdf(file = "./001_STAN_simulated/Phi_CARstan.pdf", width = 15)
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
                                    probs = c(0.025, # 0.25, 0.5, 
                                              #0.75,
                                              0.975))$summary, 4)[1:10,])
# print(fitcarstan)
# plot(fitcarstan, pars = c("alpha", "b_x1", "b_x2"))
traceplot(fitcarstan, pars = c("alpha", "b_x1", "b_x2"))
# traceplot(fitcarstan, pars = c("phi[1]", "phi[2]", "phi[3]"))
acf(fitcarstan@sim$samples[[1]]$alpha, lag.max = 2500)
acf(fitcarstan@sim$samples[[1]]$b_x1, lag.max = 2500)
acf(fitcarstan@sim$samples[[1]]$b_x2, lag.max = 2500)
acf(fitcarstan@sim$samples[[1]]$`phi[220]`, lag.max = 2500)
dev.off()

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
                label = "tab:tablePhi_CARSTAN",
                caption = "Coefficients estimated by the sparse CAR implementation in STAN 
                           on spatially autocorrelated data. Only few $\\phi$s are shown",
                filename = "./001_STAN_simulated/tablePhi_CARSTAN.tex")


names(fitcarstan)[1:3] <- c("beta_0", "beta_1", "beta_2")
names(fitcarstan)[4:2503] <-paste0("phi_",1:2500)
names(fitcarstan)[2505] <-"alpha"

png(filename = "./001_STAN_simulated/TraceplotPhi_CARSTAN.png", width = 880, height = 660)
traceplot(fitcarstan, pars = c("beta_0", "beta_1", "beta_2", "phi_220", "tau", "alpha"))
dev.off()


png(filename = "./001_STAN_simulated/ACFPhi_CARSTAN.png", width = 880, height = 660)
layout(mat = rbind(c(1,2),
                   c(3,4),
                   c(5,6)))
chain1 <- fitcarstan@sim$samples[[1]]
names(chain1)[1:3] <- c("beta_0", "beta_1", "beta_2")
names(chain1)[4:2503] <-paste0("phi_",1:2500)
names(chain1)[2505] <-"alpha"
# traceplot(fitcarstan, pars = c("phi[1]", "phi[2]", "phi[3]"))
acf(chain1$beta_0, lag.max = 2500)
acf(chain1$beta_1, lag.max = 2500)
acf(chain1$beta_2, lag.max = 2500)
acf(chain1$phi_220, lag.max = 2500)
acf(chain1$tau, lag.max = 2500)
acf(chain1$alpha, lag.max = 2500)
dev.off()

# 
extcarstan <- extract(fitcarstan)
# # str(extcarstan$phi)
interv_phi <- t(apply(extcarstan$phi, 2, FUN = function(x) quantile(x, probs = c(0.025, 0.975))))

# matplot(cbind(interv_phi, phi), type="l")


png(filename = "./001_STAN_simulated/MatplotPhi_CARSTANOrder.png", width = 880, height = 660)
# layout(mat = t(c(1,2)))
# matplot(cbind(prova, phi), type="l")
# matplot(cbind(provaincrease, phi), type="l")

modelintervalSTAN <- cbind(order=1:2500,cbind(interv_phi, phi)[order(phi),])
colnames(modelintervalSTAN)<- c("order","p025", "p975", "phi")
matplot(modelintervalSTAN[, 2:4], type="l", ylim = c(-7,7), 
        # main="Ordered phi(green) vs estimated 2.5%-97.5% interval",
        ylab="Phi", xlab="Locations ordered by phi")
abline(h = 0,lty="dashed", col="grey")
# 
# modelintervalincrease <- cbind(order=1:2500,cbind(provaincrease, phi)[order(phi),])
# colnames(modelintervalincrease)<- c("order","p10", "p90", "phi")
# matplot(modelintervalincrease[,2:4], type="l", ylim = c(-7,7), main="Ordered phi(green) vs new estimated 10%-90% interval")
# abline(h = 0, lty="dashed", col="grey")
dev.off()

png(filename = "./001_STAN_simulated/MatplotTheta_CARBayes_CARSTANOrder.png", width = 880, height = 660)
layout(t(c(1,2)))
# Theta_CARbayesinc <- theta_hat(beta = apply(modelincrease$samples$beta, 2, median),
#                               X1 = x1, X2=x2,
#                               Phi=apply(modelincrease$samples$psi, 2, median))
Theta_CARbayesinc <- cbind(prob,
                           apply(modelincrease$samples$fitted,2,mean),
                           t(apply(modelincrease$samples$fitted,2,quantile, 
                                   probs=c(0.025,0.975))))
matplot(Theta_CARbayesinc[order(prob),], 
        ylab="theta", xlab="Locations ordered by theta",type = "l", 
        main="Estimated with CARBayes")
# plot(prob[order(prob)], ylab="theta", xlab="Locations ordered by theta")
# lines(Theta_CARbayesinc[order(prob)], col="red")


# Theta_CARStan <- theta_hat(beta = summary(fitcarstan,
#                                          pars = fitcarstan@sim$pars_oi,
#                                          probs = c(0.025, 0.975))$summary[1:3,"mean"],
#                               X1 = x1, X2=x2,
#                               Phi=summary(fitcarstan,
#                                           pars = fitcarstan@sim$pars_oi,
#                                           probs = c(0.025, 0.975))$summary[4:2503,"mean"])
# plot(prob[order(prob)], ylab="theta", xlab="Locations ordered by theta")
# lines(Theta_CARStan[order(prob)], col="red")
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

dev.off()

Theta_CARbayes <- cbind(prob,
                           apply(model$samples$fitted,2,mean),
                           t(apply(model$samples$fitted,2,quantile, 
                                   probs=c(0.025,0.975))))

Coverage_theta_CARbayes <- sum(Theta_CARbayes[, "prob"]>=Theta_CARbayes[, "2.5%"] &
                                 Theta_CARbayes[, "prob"]<=Theta_CARbayes[, "97.5%"] )/K


Coverage_theta_CARbayesinc <- sum(Theta_CARbayesinc[, "prob"]>=Theta_CARbayesinc[, "2.5%"] &
                                    Theta_CARbayesinc[, "prob"]<=Theta_CARbayesinc[, "97.5%"] )/K


Coverage_theta_CARStan <- sum(Theta_CARStan[, "prob"]>=Theta_CARStan[, "2.5%"] &
                                Theta_CARStan[, "prob"]<=Theta_CARStan[, "97.5%"] )/K

rm(Logit_samples_CARStan,Theta_samples_CARStan,Theta_CARStan)

Phi_CARbayes <- cbind(apply(model$samples$psi,2,mean),
                         t(apply(model$samples$psi,2,quantile, 
                                 probs=c(0.025,0.975))))

Coverage_phi_CARbayes <- sum(phi>=Phi_CARbayes[, "2.5%"] &
                                  phi<=Phi_CARbayes[, "97.5%"] )/K

Phi_CARbayesinc <- cbind(apply(modelincrease$samples$psi,2,mean),
                         t(apply(modelincrease$samples$psi,2,quantile, 
                                 probs=c(0.025,0.975))))

Coverage_phi_CARbayesinc <- sum(phi>=Phi_CARbayesinc[, "2.5%"] &
                                  phi<=Phi_CARbayesinc[, "97.5%"] )/K

# extcarstan <- extract(fitcarstan)
# # str(extcarstan$phi)
# 
# interv_phi <- t(apply(extcarstan$phi, 2,
                      # FUN = function(x) quantile(x, probs = c(0.025, 0.975))))
Coverage_phi_CARStan <-  sum(phi>=interv_phi[, "2.5%"] &
                               phi<=interv_phi[, "97.5%"] )/K


# 
# any(interv_phi[,1]>0)
# any(interv_phi[,2]<0)
# interv_phi <- t(apply(extcarstan$phi, 2, FUN = function(x) quantile(x, probs = c(0.10, 0.90))))
# any(interv_phi[,1]>0)
# any(interv_phi[,2]<0)
# sum(interv_phi[,1]>0)
# sum(interv_phi[,2]<0)
# 
# round(phi[which(interv_phi[,1]>0)],3)
# round(apply(extcarstan$phi[,which(interv_phi[,1]>0)],2,median),3)
# round(apply(extcarstan$phi[,which(interv_phi[,1]>0)],2,mean),3)
# traceplot(fitcarstan, pars = eval(parse(text=paste("c(", paste(paste('"phi[',which(interv_phi[,1]>0), ']"', sep=""), collapse=","), ")"))))
# 
# round(phi[which(interv_phi[,2]<0)],3)
# round(apply(extcarstan$phi[,which(interv_phi[,2]<0)],2,median),3)
# round(apply(extcarstan$phi[,which(interv_phi[,2]<0)],2,mean),3)
# traceplot(fitcarstan, pars = eval(parse(text=paste("c(", paste(paste('"phi[',which(interv_phi[,2]<0), ']"', sep=""), collapse=","), ")"))))
# 
# median(extcarstan$phi[,esempio])
# phi[esempio]

# asd <- coda::mcmc(extcarstan$phi[,1:10])
# coda::gelman.diag(asd)

# # test STAN CAR 200 000 ------------
# set.seed(127)
# fitcarstanincrease <- stan(file = "./001_STAN_simulated/CARlogit_training_2pred.stan",
#                    data = datalist, iter = 50000, chains = 4,warmup = 2500,# floor(1000/3),
#                    save_dso=TRUE,
#                    verbose = TRUE )
# test STAN CAR 200 000 ------------
# 
# set.seed(127)
# increasesamplesStan <- 50000
# fitcarstanincrease <- stan(file = "./001_STAN_simulated/CARlogit_training_2pred.stan",
#                            data = datalist, thin = 5, iter = increasesamplesStan,
#                            chains = 4,warmup = 2500,# floor(1000/3),
#                            save_dso=TRUE,
#                            verbose = TRUE )
# 
# 
# TablePhi_CARSTANincrease <- knitr::kable(round(summary(fitcarstanincrease,
#                                                pars = fitcarstanincrease@sim$pars_oi,
#                                                probs = c(0.025, #0.25, 0.5, 
#                                                          # 0.75,
#                                                          0.975))$summary, 4)[c(1:10, 
#                                                                                2500:2506),],
#                                  # col.names = colnames(summary(fitglm)$coefficients),
#                                  format = "latex", 
#                                  caption = paste("\\label{tab:tablePhi_CARSTANincrease}
#                                                  Coefficients estimated by the sparse CAR implementation in STAN 
#                                                  on spatially autocorrelated data with ", increasesamplesStan," samples. Only the first $\\phi$s are shown"))
# fileTablePhi_CARSTANincrease <-  file("./001_STAN_simulated/tablePhi_CARSTANincrease.tex")
# writeLines(text = TablePhi_CARSTANincrease, fileTablePhi_CARSTANincrease)
# close(fileTablePhi_CARSTANincrease)
# 
# png(filename = "./001_STAN_simulated/TraceplotPhi_CARSTANincrease.png", width = 880, height = 660)
# traceplot(fitcarstanincrease, pars = c("alpha", "b_x1", "b_x2", "phi[220]", "alpha_car", "tau"))
# dev.off()
# 
# png(filename = "./001_STAN_simulated/ACFPhi_CARSTANincrease.png", width = 880, height = 660)
# layout(mat = rbind(c(1,2),
#                    c(3,4),
#                    c(5,6)))
# chain1 <- fitcarstanincrease@sim$samples[[1]]
# # traceplot(fitcarstan, pars = c("phi[1]", "phi[2]", "phi[3]"))
# acf(chain1$alpha, lag.max = 2500)
# acf(chain1$b_x1, lag.max = 2500)
# acf(chain1$b_x2, lag.max = 2500)
# acf(chain1$`phi[220]`, lag.max = 2500)
# acf(chain1$tau, lag.max = 2500)
# acf(chain1$alpha_car, lag.max = 2500)
# dev.off()
# 
# 
# extcarstanincrease <- extract(fitcarstanincrease)
# # str(extcarstan$phi)
# interv_phiincrease <- t(apply(extcarstanincrease$phi, 2, FUN = function(x) quantile(x, probs = c(0.1, 0.9))))
# 
# # matplot(cbind(interv_phi, phi), type="l")
# 
# 
# png(filename = "./001_STAN_simulated/MatplotPhi_CARSTANOrderincrease.png", width = 880, height = 660)
# # layout(mat = t(c(1,2)))
# # matplot(cbind(prova, phi), type="l")
# # matplot(cbind(provaincrease, phi), type="l")
# 
# modelintervalSTANincrease <- cbind(order=1:2500,cbind(interv_phiincrease, phi)[order(phi),])
# colnames(modelintervalSTANincrease)<- c("order","p10", "p90", "phi")
# matplot(modelintervalSTANincrease[, 2:4], type="l", ylim = c(-7,7), main="Ordered phi(green) vs estimated 10%-90% interval")
# abline(h = 0,lty="dashed", col="grey")
# # 
# # modelintervalincrease <- cbind(order=1:2500,cbind(provaincrease, phi)[order(phi),])
# # colnames(modelintervalincrease)<- c("order","p10", "p90", "phi")
# # matplot(modelintervalincrease[,2:4], type="l", ylim = c(-7,7), main="Ordered phi(green) vs new estimated 10%-90% interval")
# # abline(h = 0, lty="dashed", col="grey")
# dev.off()

# save models ------------
# save(fitcarstanincrease, fitcarstan, model, modelincrease, modelincrease2, file="./001_STAN_simulated/models_phi.rdata")

save(list=ls(), file=paste0("./001_STAN_simulated/models_phi_20191123_1630_bugfix_",
                            format(Sys.time(),format = "%Y%m%d_%H%M"),".rdata"))
# THE BUGFIX WAS ONLY ON THE COVERAGE OF PHI, the rest is identical to the previous version
# load("./001_STAN_simulated/models_phi_20191123_1630_.rdata")


# save(list=ls(), file=paste0("./001_STAN_simulated/models_phi_",
#                       format(Sys.time(),format = "%Y%m%d_%H%M")
#                              ,".rdata"))
# load("./001_STAN_simulated/models_phi_20191123_1630.rdata")
# layout(mat = rbind(c(1,2),
#                    c(3,4)))
# chain1 <- fitcarstan@sim$samples[[1]]
# # traceplot(fitcarstan, pars = c("phi[1]", "phi[2]", "phi[3]"))
# acf(chain1$alpha, lag.max = 2500)
# acf(chain1$b_x1, lag.max = 2500)
# acf(chain1$b_x2, lag.max = 2500)
# acf(chain1$`phi[220]`, lag.max = 2500)
