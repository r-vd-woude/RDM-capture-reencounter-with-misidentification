######################################
# This file contains function that we sue to simulate the data and then solve it with naive model and also CJSm
# Some of these functions are the same as in Kery and Schaub (2011) book
# while other were developed specifically for CJSm models

withAutoprint({
cat('\ninitializing run_CJSm_c_c_c()') 


run_CJSm_c_c_c<-function(CH, run_naive=TRUE, run_CJSm=TRUE, ni=3000, nb=1500) {
  CH_input<-CH
  if (run_naive) {
    cat('running naive CJS-c-c model')
    # flatten the data
    CH<-apply(CH, c(1,2), sign)
    # Create vector with occasion of marking
    f <- apply(CH, 1, get.first)
    Order<-(order(f))
    CH<-CH[Order,]
    f<-f[Order]
# Specify model in BUGS language
sink("cjs-c-c.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- mean.phi
      p[i,t] <- mean.p
      } #t
   } #i

mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean recapture

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2])

inits <- function(){list(mean.phi=runif(1,0,1), mean.p = runif(1,0,1), z = known.state.cjs(CH))}

# Parameters monitored
parameters <- c("mean.p", "mean.phi")

# MCMC settings
nt <- 6
nc <- 5

# Call JAGS from R (BRT 1 min)
cjs.c.c <- jags.parallel(jags.data, inits, parameters, "cjs-c-c.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, , export_obj_names=c('cjs.init.z', 'known.state.cjs', 'nb', 'ni', 'nt', 'CH', 'f'), envir=environment())

cat('   Done!\n')
  } # end of run_naive

  if (run_CJSm) {
  cat('running CJS-c-c-c model')
  CH<-CH_input
  # Create vector with occasion of marking
  f <- apply(CH, 1, get.first)
  Order<-(order(f))
  CH<-CH[Order,]
  f<-f[Order]

# Specify model in BUGS language
sink("cjs-mult-c-c-c.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- mean.phi
      p_r[i,t] <- mean.p
	  theta[i,t] <- theta.t[t]

      } #t
   } #i
mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean recapture
mean.theta ~ dunif(0,1)        # Prior for mean correct identification
for (t in 1:(n.occasions-1)) {
   theta.t[t]<-mean.theta
}

for (t in 1:(n.occasions-1)) {
   mu_wrong_assign[t] <-  N_reads[t] *(1-theta.t[t]) / sum_N_marks_used[t] 
}

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
	  lambda_obs[i,t]<- -log(1-p_r[i,t-1])*z[i,t]*theta[i,t-1] + mu_wrong_assign[t-1]
	  y[i,t] ~ dpois(lambda_obs[i,t])
      } #t
   } #i
}
",fill = TRUE)
sink()

sightings_and_fresh<-colSums(CH) # count per year

N_marked<-rle(apply(CH, 1, get.first))$lengths # marked per year

if (length(sightings_and_fresh)>length(N_marked)) N_marked<-c(N_marked, 0)

N_reads<-sightings_and_fresh-N_marked
N_reads<-N_reads[-1]

sum_N_marks_used<-cbind(cumsum(rle(apply(CH, 1, get.first))$lengths), rle(apply(CH, 1, get.first))$values)

sum_N_marks_used<-sum_N_marks_used[,1]

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], sum_N_marks_used=sum_N_marks_used, N_reads=N_reads)

inits <- function(){list(z = known.state.cjs.mult(CH), mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), mean.theta=runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.phi" , "mean.p", "mean.theta")

cjs.c.c.c <- jags.parallel(jags.data, inits, parameters, "cjs-mult-c-c-c.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,  export_obj_names=c('cjs.init.z', 'known.state.cjs.mult', 'nb', 'ni', 'nt', 'CH', 'f', 'p'), envir=environment())
cat('   Done!\n')
}

Res<-list(cjs.c.c=cjs.c.c, cjs.c.c.c=cjs.c.c.c, CH=CH)
return(Res)
}

cat('   Done\n')
}, echo=FALSE)