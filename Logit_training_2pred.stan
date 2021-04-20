data {                          
int<lower=0> N;                // number of observations
int<lower=0,upper=1> Y[N];  // setting the dependent variable (vote) as binary
vector[N] x1;             // independent variable 1
vector[N] x2;              // independent variable 2
// vector[N] age;                 // independent variable 3
}
// transformed data {
//vector[N] age_sq;              // create new variable (4), age squared (no dots in the variable name)
//age_sq = age .* age;          // formula for the variable, do not forget the . before multiplication
//}
parameters {
real alpha;                    // intercept
real b_x1;                // beta for x1, etc
real b_x2; 
//real b_age;
//real b_age_sq; 
}
model {
alpha ~ normal(0,100);         // you can set priors for all betas
b_x1 ~ normal(0,100);     // if you prefer not to, uniform priors will be used
b_x2 ~ normal(0,100);
//b_age ~ normal(0,100);
//b_age_sq ~ normal(0,100);
Y ~ bernoulli_logit(alpha + b_x1 * x1 + b_x2 * x2); # model
}
//generated quantities {         // simulate quantities of interest
//real y_hat;                    // create a new variable for the predicted values
//y_hat <- inv_logit(alpha + b_x1 * 0.2 + b_x2 * 0.8 ); // model
//}
