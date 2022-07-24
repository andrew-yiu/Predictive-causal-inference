## Treatment effects with the cattaneo dataset
## Marginal prediction with inverse probability weighting
options(java.parameters = "-Xmx14000m")
library(haven)
library(bartMachine)
library(DirichletReg)

set.seed(1)

dat <- read_dta("cattaneo2.dta") # Download dataset from http://www.stata-press.com/data/r13/cattaneo2.dta

# restrict to white and non-hispanic only
ind <- intersect(which(dat[,17]==1), which(dat[,3]==0))
dat <- dat[ind,]

y <- as.vector(t(dat[,1])) # outcome (birthweight)
t <- 1-as.vector(t(dat[,16])) # maternal smoking cessation
x <- dat[, c(2,5,7,8,9,10,11,22)] # pre-treatment covariates
x <- data.frame(x)
names(x) <- c("Marital status", "Foreign-born status", "Previous newborn death", "Mother's age", "Mother's education", "Father's age", "Father's education", "First baby")
age <- as.vector(t(dat[,8]))

Mskip <- 1000 # number of burn-in MCMC samples
M <- 2000 # number of actual MCMC samples used
n <- length(y)

# Run probit BART regression to obtain estimate of the propensity score
output_bart <- bartMachine(x, as.factor(t),verbose = FALSE, seed = 1, num_burn_in = Mskip, num_iterations_after_burn_in = M)
pidraws <- 1-bart_machine_get_posterior(output_bart,x)$y_hat_posterior_samples
pimean <- rowMeans(pidraws)

n <- length(y)

hwts1 <- t/pimean
hwts1 <- hwts1/sum(hwts1)
hwts0 <- (1-t)/(1-pimean)
hwts0 <- hwts0/sum(hwts0)
ind_0 <- which(t==0)
ind_1 <- which(t==1)

# Compute effective sample sizes for the observed method
n_til_1 <- length(ind_1)
n_til_0 <- length(ind_0)

# Generate posterior draws for the observed ESS method
ate_obs_post <- numeric(M)
for (m in 1:M) {
  ome_1 <- numeric(n)
  ome_1[ind_1] <- DirichletReg::rdirichlet(1, alpha=n_til_1*hwts1[ind_1])
  ome_0 <- numeric(n)
  ome_0[ind_0] <- DirichletReg::rdirichlet(1, alpha=n_til_0*hwts0[ind_0])
  ate_obs_post[m] <- sum(ome_1*y)-sum(ome_0*y)
}

# Compute effective sample sizes for the importance sampling method
n_til_1 <- 1/sum(hwts1^2)
n_til_0 <- 1/sum(hwts0^2)

# Generate posterior draws for the importance sampling ESS method
ate_is_post <- numeric(M)
for (m in 1:M) {
  ome_1 <- numeric(n)
  ome_1[ind_1] <- DirichletReg::rdirichlet(1, alpha=n_til_1*hwts1[ind_1])
  ome_0 <- numeric(n)
  ome_0[ind_0] <- DirichletReg::rdirichlet(1, alpha=n_til_0*hwts0[ind_0])
  ate_is_post[m] <- sum(ome_1*y)-sum(ome_0*y)
}


# Plot the posterior density for the ATE
dat_ate <-  data.frame(ate_obs_post=ate_obs_post, ate_is_post, ate_is_post)
library(ggplot2)
g <- ggplot()
g <- g + geom_density(data = dat_ate, aes(ate_obs_post, colour= "Marg-Obs"))
g <- g + geom_density(data = dat_ate, aes(ate_is_post, colour= "Marg-IS"))
g <- g + labs(colour = "Posterior") 
g <- g + xlab("Average treatment effect") + ylab("Posterior density")
g
