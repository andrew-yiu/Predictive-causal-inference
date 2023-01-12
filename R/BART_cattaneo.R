## Treatment effects with the cattaneo dataset
## Apply BART without the clever covariate
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

w <- data.frame(cbind(x,t))
x1 <- cbind(x, rep(1,length(t)))
x0 <- cbind(x, rep(0,length(t)))
names(x1)[dim(w)[2]] <- "Treatment"
names(x0)[dim(w)[2]] <- "Treatment"
names(w) <- names(x1)
x_comb <- rbind(x1, x0)

# run bart for y on w
output_bart <- bartMachine(w, y,verbose = FALSE, seed = 1, num_burn_in = Mskip, num_iterations_after_burn_in = M)

# posterior draws for t = 1
fdraws1 <- t(bart_machine_get_posterior(output_bart,x1)$y_hat_posterior_samples)
# posterior draws for t = 0
fdraws0 <- t(bart_machine_get_posterior(output_bart,x0)$y_hat_posterior_samples)

# Generate uniform Dirichlet weights for the Bayesian bootstrap
wts_mx <- matrix(,nrow = M, ncol = n)
for (m in 1:M) {
  wts_mx[m,] <- DirichletReg::rdirichlet(1, alpha=rep(1,n))
}

post_draws_bart <- matrix(,nrow =M, ncol = 35-15+1)
post_mean_bart <- numeric(35-15+1)
upp_bart <- numeric(35-15+1)
low_bart <- numeric(35-15+1)

# Compute posterior draws and intervals for the CATE for ages 15-35
for (j in 15:35) {
  ind <- which(age == j)
  sub_samp1 <- fdraws1[,ind]
  sub_samp0 <- fdraws0[,ind]
  post_mean_bart[j-15+1] <- mean(rowMeans(sub_samp1))-mean(rowMeans(sub_samp0))
  wts_sub <- as.vector(wts_mx[m,ind])
  for (m in 1:M) {
    post_draws_bart[m,j-15+1] <- sum(wts_sub*as.vector(sub_samp1[m,]))-sum(wts_sub*as.vector(sub_samp0[m,]))
    post_draws_bart[m,j-15+1] <- post_draws_bart[m,j-15+1]/sum(wts_sub)
  }
  upp_bart[j-15+1] <- quantile(as.vector(post_draws_bart[,j-15+1]), probs = 0.975)
  low_bart[j-15+1] <- quantile(as.vector(post_draws_bart[,j-15+1]), probs = 0.025)
}

# Plot cate curves and intervals
plot(15:35, post_mean_bart, type = "l", ylim = c(0,400),  xlab = "Mother's age", ylab = "CATE")
polygon(c(15:35, rev(15:35)), c(low_bart, rev(upp_bart)),
        col = "#a8d5e5", lty = 0)
lines(15:35, post_mean_bart)

# Compute posterior draws for the ATE
ate_post <- numeric(M)
for (m in 1:M) {
  wts <- as.vector(wts_mx[m,])
  ate_post[m] <- sum(wts*as.vector(fdraws1[m,]))-sum(wts*as.vector(fdraws0[m,]))
}

# Plot the posterior density for the ATE
dat_ate <-  data.frame(ate_post=ate_post)
library(ggplot2)
g <- ggplot()
g <- g + geom_density(data = dat_ate, aes(ate_post, colour= "BART"))
g <- g + labs(colour = "Posterior") 
g <- g + xlab("Average treatment effect") + ylab("Posterior density")
g
