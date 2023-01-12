## Treatment effects with the cattaneo dataset
## All methods
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

# Compute posterior draws for the ATE
ate_post <- numeric(M)
for (m in 1:M) {
  wts <- as.vector(wts_mx[m,])
  ate_post[m] <- sum(wts*as.vector(fdraws1[m,]))-sum(wts*as.vector(fdraws0[m,]))
}

set.seed(1)
# Run probit BART regression to obtain estimate of the propensity score
output_bart <- bartMachine(x, as.factor(t),verbose = FALSE, seed = 1, num_burn_in = Mskip, num_iterations_after_burn_in = M)
pidraws <- 1-bart_machine_get_posterior(output_bart,x)$y_hat_posterior_samples
pimean <- rowMeans(pidraws)

# Compute clever covariates
cc <- t/pimean - (1-t)/(1-pimean)

n <- length(y)

w <- data.frame(cbind(z,cc, t))
x1 <- cbind(x,1/pimean, rep(1,length(t)))
x0 <- cbind(x,-1/(1-pimean), rep(0,length(t)))
names(x1)[dim(w)[2]] <- "Treatment"
names(x0)[dim(w)[2]] <- "Treatment"
names(x1)[dim(w)[2] - 1] <- "Clever covariate"
names(x0)[dim(w)[2] - 1] <- "Clever covariate"
names(w) <- names(x1)
x_comb <- rbind(x1, x0)

# run bart for y on w
output_bart <- bartMachine(w, y,verbose = FALSE, seed = 1, num_burn_in = Mskip, num_iterations_after_burn_in = M)

# posterior draws for t = 1
fdraws1 <- t(bart_machine_get_posterior(output_bart,x1)$y_hat_posterior_samples)
# posterior draws for t = 0
fdraws0 <- t(bart_machine_get_posterior(output_bart,x0)$y_hat_posterior_samples)

## find inclusion proportions for the different covariates (not run by default)
#investigate_var_importance(output_bart, num_replicates_for_avg = 100)

# Generate uniform Dirichlet weights for the Bayesian bootstrap
wts_mx <- matrix(,nrow = M, ncol = n)
for (m in 1:M) {
  wts_mx[m,] <- DirichletReg::rdirichlet(1, alpha=rep(1,n))
}

post_draws_bart_cc <- matrix(,nrow =M, ncol = 35-15+1)
post_mean_bart_cc <- numeric(35-15+1)
upp_bart_cc <- numeric(35-15+1)
low_bart_cc <- numeric(35-15+1)
for (j in 15:35) {
  ind <- which(age == j)
  sub_samp1 <- fdraws1[,ind]
  sub_samp0 <- fdraws0[,ind]
  post_mean_bart_cc[j-15+1] <- mean(rowMeans(sub_samp1))-mean(rowMeans(sub_samp0))
  # post_draws <- numeric(M)
  wts_sub <- as.vector(wts_mx[m,ind])
  for (m in 1:M) {
    post_draws_bart_cc[m,j-15+1] <- sum(wts_sub*as.vector(sub_samp1[m,]))-sum(wts_sub*as.vector(sub_samp0[m,]))
    post_draws_bart_cc[m,j-15+1] <- post_draws_bart_cc[m,j-15+1]/sum(wts_sub)
  }
  upp_bart_cc[j-15+1] <- quantile(as.vector(post_draws_bart_cc[,j-15+1]), probs = 0.975)
  low_bart_cc[j-15+1] <- quantile(as.vector(post_draws_bart_cc[,j-15+1]), probs = 0.025)
}

# Compute posterior draws for the ATE
ate_post_cc <- numeric(M)
for (m in 1:M) {
  wts <- as.vector(wts_mx[m,])
  ate_post_cc[m] <- sum(wts*as.vector(fdraws1[m,]))-sum(wts*as.vector(fdraws0[m,]))
}

# Plot cate curves and intervals
par(mfrow = c(1,2))
plot(15:35, post_mean_bart, type = "l", ylim = c(0,400),  xlab = "Mother's age", ylab = "CATE")
polygon(c(15:35, rev(15:35)), c(low_bart, rev(upp_bart)),
        col = "#a8d5e5", lty = 0)
lines(15:35, post_mean_bart)

plot(15:35, post_mean_bart_cc, type = "l", ylim = c(0,400),  xlab = "Mother's age", ylab = "CATE")
polygon(c(15:35, rev(15:35)), c(low_bart_cc, rev(upp_bart_cc)),
        col = "#a8d5e5", lty = 0)
lines(15:35, post_mean_bart_cc)
par(mfrow = c(1,1))

## Marginal prediction with inverse probability weighting
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

## Create the density strips plot in Figure 7
library(denstrip)

col1 <- "#FE9929"; col2="#4292C6"; col3="#41AB5D"; col4="#6A51A3"
colmin <- "gray92"
wd <- 0.3
ci <- c(0.025, 0.975)
base <- 1
dy <- 0.5
dg <- 0.7
base2 <- base + 4*dy + dg
base3 <- base2 + 3*dy + dg

prev.plot <- function(sam1, sam2, sam3, sam4) { 
  
  par(mar=c(3, 0, 0, 0), mgp=c(2,1,0))
  
  xmax <- 370
  plot(0, type="n", xlim=c(110, xmax), ylim=c(1,3), axes=FALSE, xlab="Average treatment effect", ylab="")
  lim <- par("usr")
  rect(0.6, lim[3], 1.7, lim[4], col="gray92", border="gray92")
  axis(1, at=seq(150,370,by=20))
  
  denstrip(sam1, at=base+0.2 + 3*dy, width=wd, colmax=col1, colmin=colmin, from=150, to=370, ticks=quantile(sam1, ci))
  denstrip(sam2,        at=base+0.2 + 2*dy, width=wd, colmax=col2, colmin=colmin, from=150, to=370, ticks=quantile(sam2, ci))
  denstrip(sam3,   at=base+0.2 + dy,  width=wd, colmax=col3, colmin=colmin, from=150, to=370, ticks=quantile(sam3, ci))
  denstrip(sam4, at=base + 0.2,  width=wd, colmax=col4, colmin=colmin, from=150, to=370, ticks=quantile(sam4, ci))
  
  text(110, base+0.2 + c(3*dy,2*dy,dy,0), c("BART", "BART-CC", "Marg-Obs", "Marg-IS"), pos=4, col=c(col1,col2,col3, col4))
  
}

prev.plot(ate_post, ate_post_cc, ate_obs_post, ate_is_post)

