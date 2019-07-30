# Supplementary code to 'Biased survival estimates by Cormack-Jolly-Seber models due to misidentification: evidence and the solution' by E. Rakhimberdiev et al., Methods in Ecology and Evolution, submitted

_code by Eldar Rakhimberdiev_

You might enjoy [html version of this supplementary](http://htmlpreview.github.io/?https://raw.githubusercontent.com/eldarrak/CJS-with-misidentification/master/code/All_code.html).

## 0. Set up output and run options
```{r}
set.seed(123)
```
### 0.1 To run or to download?
Parts of this code take time to run, so if you want to go through examples fast and do not want to install all the software
```{r}
Run_everything=FALSE # Change to TRUE if you want run all the models on your side
```

### 0.2 Online or local data?
This script can operate with online or local data. If `Use_local_data` set to TRUE working directory expected to have folder `data`. If `Run_everything` set to TRUE `results` directory will be created and all results will be saved there.
```{r}
Use_local_data=FALSE # Change to TRUE if you want to use local data and to not download them from the GitHub page.
```

### 0.3 Checks and downloads 
```{r}
if (Use_local_data) {
   if (!'results' %in% list.dirs(full.names=FALSE, recursive=FALSE)) {
       if (Run_everything) {
          dir.create('results')
	   } else {
	      stop('results directory not found!')
	   }
   } 
   if (!'data' %in% list.dirs(full.names=FALSE, recursive=FALSE)) stop('data directory not found!')
} else {
   if (!'results' %in% list.dirs(full.names=FALSE, recursive=FALSE)) dir.create('results')
   # we will now download all the result files
   Parameter_combos<-expand.grid(Phi=c(0.9, 0.5), p=c(0.9, 0.5, 0.3), Theta=c(1, 0.99, 0.95, 0.9))
   Results_dir<-'https://github.com/eldarrak/CJS-with-misidentification/blob/master/results/'
   for (i in 1:nrow(Parameter_combos)) {
      Filename<-paste0('Sim_study_100_simulations_Phi_dot_', 
        Parameter_combos$Phi[i], '_P_dot_', Parameter_combos$p[i], 
        '_theta_dot_', Parameter_combos$Theta[i],
	    '_n_occasions_25_marked_25_RHat1.01_.RData')
      download.file(paste0(Results_dir, Filename, '?raw=true'),
        destfile=paste0('./results/', Filename), mode='wb')
   }
   download.file('https://git.io/fjeU2', destfile='./results/cjs.T.c.no_misr.RDS', mode='wb')
   download.file('https://git.io/fjeUV',
          destfile='./results/cjs.T.c.with.misidentification.RDS', mode='wb')
   download.file('https://git.io/fjeUw',
          destfile='./results/cjsm.T.c.c.RDS', mode='wb')
   download.file('https://git.io/fjeU1', destfile='./results/cjs.naive.godwits.RDS', mode='wb')
   download.file('https://git.io/fjeUy',
          destfile='./results/cjsm.godwits.RDS', mode='wb')
}
if (Use_local_data & !'data' %in% list.dirs(full.names=FALSE, recursive=FALSE)) {
    stop('data directory not found') 
}

```

### 0.4 load libraries and source CJSm functions
We need `jagsUI` and/or `R2jags`to run jags code. Note that to run the models you will also have to install jags software (Plummer 2011) from [here]( http://mcmc-jags.sourceforge.net).

```{r}
library('jagsUI')
```

Now we source helper functions
```{r}
# general functions
source('https://git.io/fhhtZ')

# some locally used functions
source('https://git.io/fjPgo')
```

## 1 Simulate and solve Phi_dot_P_dot_Theta_dot models
For this exersice we simulate data with the function `simul.cjs.multiple.sightings` and then with the function `run_CJSm_c_c_c` estimate two cjs models over these data Phi_dot_P_dot model (CJS-c-c) and Phi_dot_P_dot_Theta_dot (cjsm-c-c-c)

### 1.1 Define wrapper function to run simulations and run the models output
```{r}
run_CJSm_c_c_c<-function(CH, run_naive=TRUE, run_CJSm=TRUE, 
                         n.adapt=5000, Rhat.limit=1.01, 
						 max.iter=150000, iter.increment2=10000) {
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
    if (TRUE) {
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
   }
   # Bundle data
   jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2])

   inits <- function(){list(mean.phi=runif(1,0,1), mean.p = runif(1,0,1), z = known.state.cjs(CH))}

   # Parameters monitored
   parameters <- c("mean.p", "mean.phi")

   # Call JAGS from R 
   cjs.c.c <- autojags(jags.data, inits, parameters, "cjs-c-c.jags",
                       n.chains = 6, parallel=TRUE, 
					   n.cores=6,Rhat.limit=Rhat.limit, n.adapt=n.adapt,
					   max.iter=max.iter, iter.increment=4000)
   
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
  if (TRUE) {
     # Specify model in BUGS language
    sink("cjs-mult-c-c-c.jags")
    cat("
    model {

    # Priors and constraints
    for (i in 1:nind){
       for (t in f[i]:(n.occasions-1)){
          phi[i,t] <- mean.phi
          p_r[i,t] <- mean.p
       } #t
    } #i
    mean.phi ~ dunif(0, 1)         # Prior for mean survival
    mean.p ~ dunif(0, 1)           # Prior for mean recapture
    mean.theta ~ dunif(0,1)        # Prior for mean correct identification
    for (t in 1:(n.occasions-1)) {
    theta.t[t]<-mean.theta
  }

  for (t in 1:(n.occasions-1)) {
      mu_wrong_assign[t] <-  N_reads[t] *(1-theta.t[t]) /
                             (sum_N_marks_used[t]-1)
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
    	  lambda_obs[i,t]<- -log(1-p_r[i,t-1])*z[i,t]*theta.t[t-1]
                             + mu_wrong_assign[t-1]
    	  y[i,t] ~ dpois(lambda_obs[i,t])
        } #t
    } #i
  }
  ",fill = TRUE)
  sink()
  }
   sightings_and_fresh<-colSums(CH) # count per year
   N_marked<-rle(apply(CH, 1, get.first))$lengths # marked per year
   if (length(sightings_and_fresh)>length(N_marked)) N_marked<-c(N_marked, 0)

   N_reads<-sightings_and_fresh-N_marked
   N_reads<-N_reads[-1]

   sum_N_marks_used<-cbind(cumsum(rle(apply(CH, 1, get.first))$lengths),
                           rle(apply(CH, 1, get.first))$values)

   sum_N_marks_used<-sum_N_marks_used[,1]

   # Bundle data
   jags.data <- list(y = CH, f = f, nind = dim(CH)[1],
                     n.occasions = dim(CH)[2],
					 sum_N_marks_used=sum_N_marks_used,
					 N_reads=N_reads)

   inits <- function(){list(z = known.state.cjs.mult(CH), 
                            mean.phi = runif(1, 0, 1),
                            mean.p = runif(1, 0, 1),
                            mean.theta=runif(1, 0, 1))}

   # Parameters monitored
   parameters <- c("mean.phi" , "mean.p", "mean.theta")

   # Call JAGS from R 
   cjs.c.c.c <- autojags(jags.data, inits, parameters,
                       "cjs-mult-c-c-c.jags", n.chains = 6, 
					   parallel=TRUE, n.cores=6, 
					   Rhat.limit=Rhat.limit, n.adapt=NULL, 
					   max.iter=max.iter, iter.increment= iter.increment2)
   cat('   Done!\n')
   }

   Res<-list(cjs.c.c=cjs.c.c, cjs.c.c.c=cjs.c.c.c, CH=CH)
   return(Res)
}
```
### 1.2 Define parameters for a sample simulation (Phi=0.9, P=0.9, theta=0.95)
```{r}
n.occasions <- 25                   # Number of capture occasions
marked <- rep(25, n.occasions-1) # Annual number of newly marked individuals

Phi_dot=0.9
P_dot=0.3
theta_dot=0.95
RHat_limit=1.01
max.iter=160000
n.adapt=5000
iter.increment2=20000
nRuns=100
n.occasions <-25  

phi <- rep(Phi_dot, n.occasions-1) # survival
p <- rep(P_dot, n.occasions-1)  # resighting
# Define matrices with survival recapture and correct identification
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

theta<- rep(theta_dot, n.occasions-1) # correctly identified
THETA<-matrix(theta, ncol = n.occasions-1, nrow = sum(marked))
```

### 1.3 Run loop for 100 simulations
Note, if you turn `Run_everything` to `TRUE`  the loop below will try to run 200 models! and it will take time (days). 

```{r}
if (Run_everything) {
for (i in 1:nRuns) {
   cat('i = ', i, '\n')
   CH <- simul.cjs.multiple.sightings(PHI, P, THETA,  marked)
   Res_list<-run_CJSm_c_c_c(CH$CH, Rhat.limit=RHat_limit, 
                            max.iter=max.iter, 
							iter.increment2=iter.increment2, 
							n.adapt=n.adapt)
   Mean.c.c.c<-rbind(Mean.c.c.c, 
                   c('mean.phi'=Res_list$cjs.c.c.c$mean$mean.phi, 
				   'mean.p'=Res_list$cjs.c.c.c$mean$mean.p, 
				   'mean.theta'=Res_list$cjs.c.c.c$mean$mean.theta,
				   'deviance'=Res_list$cjs.c.c.c$mean$deviance)) 
   Mean.c.c<-rbind(Mean.c.c,
                   c('mean.phi'=Res_list$cjs.c.c$mean$mean.phi, 
				   'mean.p'=Res_list$cjs.c.c$mean$mean.p, 
				   'deviance'=Res_list$cjs.c.c$mean$deviance))

   ccc_CI_cur<-cbind(c('mean.phi'=phi[1], 'mean.p'=p[1], 
                    'mean.theta'=theta[1]), 
                  c(Res_list$cjs.c.c.c$q2.5$mean.phi, 
				    Res_list$cjs.c.c.c$q2.5$mean.p,
					Res_list$cjs.c.c.c$q2.5$mean.theta), 
				  c(Res_list$cjs.c.c.c$q97.5$mean.phi, 
				    Res_list$cjs.c.c.c$q97.5$mean.p,
					Res_list$cjs.c.c.c$q97.5$mean.theta))

   cc_CI_cur<-cbind(c('mean.phi'=phi[1], 'mean.p'=p[1]), 
                  c(Res_list$cjs.c.c$q2.5$mean.phi, 
				    Res_list$cjs.c.c$q2.5$mean.p), 
				  c(Res_list$cjs.c.c$q97.5$mean.phi, 
				    Res_list$cjs.c.c$q97.5$mean.p))

   Cov.c.c.c<-rbind(Cov.c.c.c,apply(ccc_CI_cur, 1, FUN=function(x) x[[1]]>=x[[2]] & x[[1]]<=x[[3]]))
   Cov.c.c<-rbind(Cov.c.c,apply(cc_CI_cur, 1, FUN=function(x) x[[1]]>=x[[2]] & x[[1]]<=x[[3]]))

   ccc_CI_width<-rbind(ccc_CI_width, apply(ccc_CI_cur, 1, FUN=function(x) x[[3]]- x[[2]]))
   cc_CI_width<-rbind(cc_CI_width, apply(cc_CI_cur, 1, FUN=function(x) x[[3]]- x[[2]]))

   N_eff_ccc<-rbind(N_eff_ccc, unlist(Res_list$cjs.c.c.c$n.eff))
   Rhat_ccc<-rbind(Rhat_ccc, unlist(Res_list$cjs.c.c.c$Rhat))

   N_eff_cc<-rbind(N_eff_cc, unlist(Res_list$cjs.c.c$n.eff))
   Rhat_cc<-rbind(Rhat_cc, unlist(Res_list$cjs.c.c$Rhat))

   MCMCInfo_ccc<-rbind(MCMCInfo_ccc, 
      c('n.chains'=Res_list$cjs.c.c.c$mcmc.info$n.chains,
        'n.burnin'=Res_list$cjs.c.c.c$mcmc.info$n.burnin,
        'n.iter'=Res_list$cjs.c.c.c$mcmc.info$n.iter,
        'elapsed.mins'=Res_list$cjs.c.c.c$mcmc.info$elapsed.mins))

   MCMCInfo_cc<-rbind(MCMCInfo_cc, 
       c('n.chains'=Res_list$cjs.c.c$mcmc.info$n.chains,
         'n.burnin'=Res_list$cjs.c.c$mcmc.info$n.burnin,
         'n.iter'=Res_list$cjs.c.c$mcmc.info$n.iter,
         'elapsed.mins'=Res_list$cjs.c.c$mcmc.info$elapsed.mins))
   #######
   # once in e.g. 10 runs we will save the result..
   if (i%%10 ==0) { 
     save(Mean.c.c.c, Mean.c.c, Cov.c.c.c, Cov.c.c, 
        ccc_CI_width, cc_CI_width, 
	    N_eff_ccc, N_eff_cc, 
		Rhat_ccc, Rhat_cc, 
		MCMCInfo_ccc, MCMCInfo_cc,
		file=paste0('LBurnin_Phi_dot_', Phi_dot, '_P_dot_', P_dot, 
            '_theta_dot_', theta_dot, '_n_occasions_', n.occasions, 
	        '_marked_', marked[1], '_iteration_', i, '_RHat_', RHat_limit, '.Rdata'))
    }
   }
			 
   N_ccc<-  length(which(apply(Rhat_ccc, 1,max)<=RHat_limit))

   cur_parameters<-list(Phi=Phi, p=p, theta=theta, RHat=RHat_limit, N_ccc=N_ccc)

   Filename<-paste0('Sim_study_100_simulations_Phi_dot_', 
      Phi, '_P_dot_', p, 
      '_theta_dot_', theta,
	  '_n_occasions_25_marked_25_RHat1.01_.RData')	
	  
   save(cur_parameters,
     cc_CI_width,  ccc_CI_width, 
     Cov.c.c, Cov.c.c.c, 
     MCMCInfo_cc, MCMCInfo_ccc, 
	 Mean.c.c, Mean.c.c.c,
	 N_eff_cc, N_eff_ccc,
	 Rhat_cc, Rhat_ccc,
	 file=Filename)
 
}
```

### 1.4 Combine simulation results
```{r}
C_Files<-list.files(path='./results/', 
                    pattern='Sim_study_100_simulations', 
					recursive=TRUE,full.names=TRUE)
Res<-c()
for (i in 1:length(C_Files)) {
   Res<-rbind(Res, make_two_cells(C_Files[i]))
}

print(Res, digits=3)
```
### 1.5 Make Figure 3

First we prepare indices for colors, color intencity and jitter.
```{r}
Res$pchtype<-1
Res$pchtype[Res$Phi==0.5]<-2

Res$colintensity<-NA
Res$colintensity[Res$p==0.3]<-1
Res$colintensity[Res$p==0.5]<-2
Res$colintensity[Res$p==0.9]<-3

Res$X<-NA
Res$X[Res$theta==1]<-1
Res$X[Res$theta==0.99]<-2
Res$X[Res$theta==0.95]<-3
Res$X[Res$theta==0.9]<-4

# this order is needed to shift points to the left and to the right..
Res$X.order<-NA
Res$X.order[Res$Phi==0.9 & Res$p==0.9 & Res$type=='naive']<-1
Res$X.order[Res$Phi==0.9 & Res$p==0.9 & Res$type=='CJSm']<-2

Res$X.order[Res$Phi==0.9 & Res$p==0.5 & Res$type=='naive']<-3
Res$X.order[Res$Phi==0.9 & Res$p==0.5 & Res$type=='CJSm']<-4

Res$X.order[Res$Phi==0.9 & Res$p==0.3 & Res$type=='naive']<-5
Res$X.order[Res$Phi==0.9 & Res$p==0.3 & Res$type=='CJSm']<-6

Res$X.order[Res$Phi==0.5 & Res$p==0.9 & Res$type=='naive']<-7
Res$X.order[Res$Phi==0.5 & Res$p==0.9 & Res$type=='CJSm']<-8

Res$X.order[Res$Phi==0.5 & Res$p==0.5 & Res$type=='naive']<-9
Res$X.order[Res$Phi==0.5 & Res$p==0.5 & Res$type=='CJSm']<-10

Res$X.order[Res$Phi==0.5 & Res$p==0.3 & Res$type=='naive']<-11
Res$X.order[Res$Phi==0.5 & Res$p==0.3 & Res$type=='CJSm']<-12
```
Now plotting the figure:
```{r}
line_col=grey(0.5)
linelty=2
par(mar=c(3,3,2,1))
par(mfrow=c(2,3))
orange.colors<-colorRampPalette(c('white', '#D95F02') )

##################################
### Panel A
##### Phi bias
plot(x=c(0.5,4.5), y=c(-0.1, 0.5),
     xaxs='i', type='n', axes=F,
	 xlab='Theta', ylab='Phi bias')

axis(1, labels=c(1, 0.99, 0.95, 0.9), at=c(1:4))
axis(2, las=1)
abline(h=seq(-0.4, 0.4, by=0.2), col=line_col, lty=linelty)
abline(v=c(1.5, 2.5, 3.5), col=line_col, lty=linelty)
abline(h=0, lwd=1, col='red')
box()

sets4lines<-unique(subset(Res, select=c('Phi', 'p', 'type')))
for (i in 1:nrow(sets4lines)) {
   Cur_data<-Res[Res$Phi==sets4lines$Phi[i] 
                 & Res$p==sets4lines$p[i]
				 & Res$type==sets4lines$type[i],]
   Cur_data<-Cur_data[order(Cur_data$theta),]
   lines(Phi.bias ~ I(X + (X.order-6.5)/14), data=Cur_data,
         pch=ifelse(pchtype==1, 15, 17), col=ifelse(type=='naive',
         orange.colors(4)[colintensity+1],
		 rev(grey.colors(4))[colintensity+1]))
}

points(Phi.bias~ I(X + (X.order-6.5)/14), data=Res, 
       pch=ifelse(pchtype==1, 15, 17), 
	   col=ifelse(type=='naive', orange.colors(4)[colintensity+1],
	   rev(grey.colors(4))[colintensity+1]))

##################################
### Panel B
##### p bias
plot(x=c(0.5,4.5), y=range(Res$p.bias),
     xaxs='i', type='n', axes=F,
	 xlab='Theta', ylab='p bias')

axis(1, labels=c(1, 0.99, 0.95, 0.9), at=c(1:4))
axis(2, las=1)
abline(h=seq(-0.4, 0.4, by=0.2), col=line_col, lty=linelty)
abline(v=c(1.5, 2.5, 3.5), col=line_col, lty=linelty)
abline(h=0, lwd=1, col='red')
box()
sets4lines<-unique(subset(Res, select=c('Phi', 'p', 'type')))

for (i in 1:nrow(sets4lines)) {
   Cur_data<-Res[Res$Phi==sets4lines$Phi[i]
                & Res$p==sets4lines$p[i] 
				& Res$type==sets4lines$type[i],]
   Cur_data<-Cur_data[order(Cur_data$theta),]
   lines(p.bias ~ I(X + (X.order-6.5)/14), data=Cur_data,
         pch=ifelse(pchtype==1, 15, 17), 
		 col=ifelse(type=='naive', orange.colors(4)[colintensity+1],
		 rev(grey.colors(4))[colintensity+1]))
}

points(p.bias~ I(X + (X.order-6.5)/14), data=Res, 
       pch=ifelse(pchtype==1, 15, 17), 
	   col=ifelse(type=='naive', orange.colors(4)[colintensity+1],
	   rev(grey.colors(4))[colintensity+1]))

##################################
### Panel C
##### theta bias
plot(x=c(0.5,4.5), y=range(-0.3, 0.3), 
     xaxs='i', type='n', axes=F, 
	 xlab='Theta', ylab='theta bias')

axis(1, labels=c(1, 0.99, 0.95, 0.9), at=c(1:4))
axis(2, las=1)
abline(v=c(1.5, 2.5, 3.5),, col=line_col, lty=linelty)
abline(h=seq(-0.4, 0.4, by=0.2), col=line_col, lty=linelty)
abline(h=0, lwd=1, col='red')
box()

sets4lines<-unique(subset(Res, select=c('Phi', 'p', 'type')))
for (i in 1:nrow(sets4lines)) {
    Cur_data<-Res[Res$Phi==sets4lines$Phi[i]
               	  & Res$p==sets4lines$p[i]
	              & Res$type==sets4lines$type[i],]
    Cur_data<-Cur_data[order(Cur_data$theta),]

    lines(theta.bias ~ I(X + (X.order-6.5)/14), data=Cur_data,
          pch=ifelse(pchtype==1, 15, 17), col=ifelse(type=='naive', 
		  orange.colors(4)[colintensity+1],
		  rev(grey.colors(4))[colintensity+1]))
}

points(theta.bias~ I(X + (X.order-6.5)/14), data=Res, 
       pch=ifelse(pchtype==1, 15, 17), col=ifelse(type=='naive',
	   orange.colors(4)[colintensity+1],
	   rev(grey.colors(4))[colintensity+1]))

##########################################
##########################################
# the second set of panels - coverage

##################################
### Panel D
##### Phi coverage
plot(x=c(0.5,4.5), y=c(0,1),  xaxs='i',
     type='n', axes=F, xlab='Theta', 
	 ylab='Phi 95 CI coverage')

axis(1, labels=c(1, 0.99, 0.95, 0.9), at=c(1:4))
axis(2, las=1)
abline(h=0.95, lwd=1, col='red')
abline(v=c(1.5, 2.5, 3.5), col=line_col, lty=linelty)
abline(h=c(0, 0.25, 0.5, 0.75, 1), col=line_col, lty=linelty)
box()

sets4lines<-unique(subset(Res, select=c('Phi', 'p', 'type')))
for (i in 1:nrow(sets4lines)) {
   Cur_data<-Res[Res$Phi==sets4lines$Phi[i] 
                & Res$p==sets4lines$p[i] 
                & Res$type==sets4lines$type[i],]
   Cur_data<-Cur_data[order(Cur_data$theta),]

   lines(Phi.CI.coverage~ I(X + (X.order-6.5)/14), data=Cur_data,
         pch=ifelse(pchtype==1, 15, 17), col=ifelse(type=='naive', 
		 orange.colors(4)[colintensity+1],
		 rev(grey.colors(4))[colintensity+1]))
}

points(Phi.CI.coverage~ I(X + (X.order-6.5)/14), data=Res, 
       pch=ifelse(pchtype==1, 15, 17), col=ifelse(type=='naive', 
	   orange.colors(4)[colintensity+1],
	   rev(grey.colors(4))[colintensity+1]))

##################################
### Panel E
##### p coverage
plot(x=c(0.5,4.5), y=c(0,1), type='n', axes=F,
     xlab='Theta', ylab='p 95 CI coverage')

axis(1, labels=c(1, 0.99, 0.95, 0.9), at=c(1:4))
axis(2, las=1)
abline(h=0.95, lwd=1, col='red')
abline(v=c(1.5, 2.5, 3.5), col=line_col, lty=linelty)
abline(h=c(0, 0.25, 0.5, 0.75, 1), col=line_col, lty=linelty)
box()

sets4lines<-unique(subset(Res, select=c('Phi', 'p', 'type')))
for (i in 1:nrow(sets4lines)) {
    Cur_data<-Res[Res$Phi==sets4lines$Phi[i]
   	& Res$p==sets4lines$p[i] 
	& Res$type==sets4lines$type[i],]
    Cur_data<-Cur_data[order(Cur_data$theta),]
    lines(p.CI.coverage~ I(X + (X.order-6.5)/14), data=Cur_data,
          pch=ifelse(pchtype==1, 15, 17), col=ifelse(type=='naive', 
		  orange.colors(4)[colintensity+1],
		  rev(grey.colors(4))[colintensity+1]))
}

points(p.CI.coverage~ I(X + (X.order-6.5)/14), data=Res, 
       pch=ifelse(pchtype==1, 15, 17), col=ifelse(type=='naive',
	   orange.colors(4)[colintensity+1],
	   rev(grey.colors(4))[colintensity+1]))

##################################
### Panel F
##### theta coverage
plot(x=c(0.5,4.5), y=c(0,1),, type='n', 
     axes=F, xlab='Theta', ylab='theta 95 CI coverage')

axis(1, labels=c(1, 0.99, 0.95, 0.9), at=c(1:4))
axis(2, las=1)
abline(h=0.95, lwd=1, col='red')
abline(v=c(1.5, 2.5, 3.5), col=line_col, lty=linelty)
abline(h=c(0, 0.25, 0.5, 0.75, 1), col=line_col, lty=linelty)
box()
Res_c<-Res[Res$theta<1,]

sets4lines<-unique(subset(Res_c, select=c('Phi', 'p', 'type')))
for (i in 1:nrow(sets4lines)) {
    Cur_data<-Res_c[Res_c$Phi==sets4lines$Phi[i] 
	& Res_c$p==sets4lines$p[i] 
	& Res_c$type==sets4lines$type[i],]
    Cur_data<-Cur_data[order(Cur_data$theta),]
    lines(theta.CI.coverage~ I(X + (X.order-6.5)/14), data=Cur_data,
          pch=ifelse(pchtype==1, 15, 17), 
		  col=ifelse(type=='naive', orange.colors(4)[colintensity+1],
		  rev(grey.colors(4))[colintensity+1]))
}

points(theta.CI.coverage~ I(X + (X.order-6.5)/14), data=Res_c, 
       pch=ifelse(pchtype==1, 15, 17), 
	   col=ifelse(type=='naive', orange.colors(4)[colintensity+1],
	   rev(grey.colors(4))[colintensity+1]))
```

## 2. Estimation of trend in survival probability over time
In this part we simulate data with no trend in survival over time and then run naive CJS and CJSm with Time trend

### 2.1 simulate data
#### Define parameter values
```{r}
n.occasions <- 25                   # Number of capture occasions
marked <- rep(100, n.occasions-1)   # Annual number of newly marked individuals
phi <- rep(0.7, n.occasions-1) # survival
p <- rep(0.9, n.occasions-1)  # resighting
theta<- rep(0.95, n.occasions-1) # correctly identified
```

#### Define matrices with survival recapture and correct identification
```{r}
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))
THETA<-matrix(theta, ncol = n.occasions-1, nrow = sum(marked))
```

#### Simulate data
```{r}
CH_misr <- simul.cjs.multiple.sightings(PHI, P, THETA,  marked)
```

### 2.2 run naive CJS model on data _with_ misidentificaion

```{r}
CH<-CH_misr$CH
# flatten the data
CH<-apply(CH_misr$CH, c(1,2), sign) # all positive values become 1s

# Create vector with occasion of marking
f <- apply(CH, 1, get.first)
Order<-order(f)
CH<-CH[Order,]
f<-f[Order]

# Specify model in BUGS language
sink("cjs-T-c.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- phi.mu+phi.beta*t #intercept and slope for survival
      p[i,t] <- mean.p
      } #t
   } #i

phi.mu ~ dnorm(0,0.1)          # Prior for intercept in survival
phi.beta ~ dnorm(0,1)          # Prior for slope in survival

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
jags.data <- list(y = CH, f = f, nind = dim(CH)[1],
                  n.occasions = dim(CH)[2])

inits <- function(){list(phi.mu = 1, phi.beta = 0, mean.p = 0.7,
                         z = known.state.cjs(CH))}

# Parameters monitored
parameters <- c("mean.p", "phi.mu", "phi.beta")

# MCMC settings
ni <- 7000
nt <- 1
nb <- 2500
nc <- 6

detach('package:jagsUI', unload=TRUE) 
# we used this package for the first set of simulations
# It is safer to detach it before uploading R2jags
library(R2jags)

if (Run_everything) {
# Call JAGS from R (BRT ~ 30 min)
cjs.T.c.with.misidentification <- jags.parallel(jags.data, inits,
               parameters, "cjs-T-c.jags",
			   n.chains = nc, n.thin = nt,
			   n.iter = ni, n.burnin = nb,
			   export_obj_names=c('cjs.init.z', 'known.state.cjs', 
			                      'nb', 'ni', 'nt', 'CH', 'f'))

saveRDS(cjs.T.c.with.misidentification, 
        file='./results/cjs.T.c.with.misidentification.RDS')
}

cjs.T.c.with.misidentification<-readRDS('./results/cjs.T.c.with.misidentification.RDS')
print(cjs.T.c.with.misidentification, digits=3)
```

### 2.3 run naive CJS model on data _without_ misidentificaion
Here we will run a naive CJS model over the 'good' data.

```{r}
# flatten the data
CH_no_misr<-apply(CH_misr$CH_no_misr, c(1,2), sign)

# Create vector with occasion of marking
f_no_misr <- apply(CH_no_misr, 1, get.first)
Order_nm<-(order(f_no_misr))
CH_no_misr<-CH_no_misr[Order_nm,]
f_no_misr<-f_no_misr[Order_nm]

# Bundle data
jags.data.no_misr <- list(y = CH_no_misr, f = f_no_misr,
                          nind = dim(CH_no_misr)[1], n.occasions = dim(CH_no_misr)[2])

inits.no_misr <- function(){list(phi.mu = 1, phi.beta = 0,
                                 mean.p = 0.7, z = known.state.cjs(CH_no_misr))}

# Parameters monitored
parameters <- c("mean.p", "phi.mu", "phi.beta")

# MCMC settings
ni <- 5000
nt <- 1
nb <- 2500
nc <- 6

if (Run_everything) {
# Call JAGS from R (BRT ~ 30 min)

cjs.T.c.no_misr <- jags.parallel(jags.data.no_misr,
                           inits.no_misr, parameters, "cjs-T-c.jags",
						   n.chains = nc, n.thin = nt,
						   n.iter = ni, n.burnin = nb,
						   export_obj_names=c('cjs.init.z', 'known.state.cjs',
						   'nb', 'ni', 'nt', 
						   'CH_no_misr', 'f_no_misr'))

saveRDS(cjs.T.c.no_misr, 
        file='./results/cjs.T.c.no_misr.RDS')
}

cjs.T.c.no_misr<-readRDS('./results/cjs.T.c.no_misr.RDS')
print(cjs.T.c.no_misr, digits=3)
```

### 2.4 run CJSm model on data _with_ misidentificaion
We now run Phi_T_P_dot_Theta_dot model over the data with misidentification

```{r}
CH<-CH_misr$CH
# Create vector with occasion of marking
f <- apply(CH, 1, get.first)
Order<-(order(f))
CH<-CH[Order,]
f<-f[Order]

# Specify model in BUGS language
sink("cjsm-T-c-c.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- phi.mu+phi.beta*t
      p_r[i,t] <- mean.p
	  theta[i,t] <- theta.t[t]
      } #t
   } #i

phi.mu ~ dnorm(0,0.1)          # Prior for intercept in survival
phi.beta ~ dnorm(0,1)          # Prior for slope in survival
mean.p ~ dunif(0, 1)           # Prior for mean recapture
mean.theta ~ dunif(0.5,1)        # Prior for mean correct identification
for (t in 1:(n.occasions-1)) {
   theta.t[t]<-mean.theta
}

for (t in 1:(n.occasions-1)) {
   mu_wrong_assign[t] <-  N_reads[t] *(1-theta.t[t]) / (sum_N_marks_used[t]-1)
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
	  lambda_obs[i,t]<- -log(1-p_r[i,t-1])*z[i,t]*theta[i,t-1] +
	                         mu_wrong_assign[t-1]
	  y[i,t] ~ dpois(lambda_obs[i,t])
      } #t
   } #i
}
",fill = TRUE)
sink()

sightings_and_fresh<-colSums(CH_misr$CH) # count per year

N_marked<-rle(apply(CH, 1, get.first))$lengths # marked per year

if (length(sightings_and_fresh)>length(N_marked)) N_marked<-c(N_marked, 0)

N_reads<-sightings_and_fresh-N_marked
N_reads<-N_reads[-1]

sum_N_marks_used<-cbind(cumsum(rle(apply(CH, 1, get.first))$lengths),
                           rle(apply(CH, 1, get.first))$values)

sum_N_marks_used<-sum_N_marks_used[,1]

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],
                  sum_N_marks_used=sum_N_marks_used, N_reads=N_reads)

inits <- function(){list(phi.mu = 1, phi.beta = 0, mean.p = unique(p),
                         z = known.state.cjs.mult(CH),mean.theta=runif(1, 0.5, 1))}

# Parameters monitored
parameters <- c("phi.mu", "phi.beta" , "mean.p", "mean.theta")

# MCMC settings
ni <- 5000
nt <- 1
nb <- 2500
nc <- 6

if (Run_everything) {
   # (BRT ~ 29 min)
   cjsm.T.c.c <- jags.parallel(jags.data, inits, parameters, "cjsm-T-c-c.jags",
                     n.chains = nc, n.thin = nt,
					 n.iter = ni, n.burnin = nb,
					 export_obj_names=c('cjs.init.z',
                       'known.state.cjs.mult', 'nb',
					   'ni', 'nt', 'CH', 'f', 'p'))
   saveRDS(cjsm.T.c.c, file='./results/cjsm.T.c.c.RDS')
}

cjsm.T.c.c<-readRDS('./results/cjsm.T.c.c.RDS')
print(cjsm.T.c.c, digits=3)
```

### 2.5 Plot the results

#### Extract model outputs
```{r}
XX<-seq(1, n.occasions, by=1)

CJS_no_mr_Median<-sapply(XX, FUN=function(x) 
       plogis(quantile(cjs.T.c.no_misr$BUGSoutput$sims.list$phi.mu +
       x*cjs.T.c.no_misr$BUGSoutput$sims.list$phi.beta, 0.5)))
CJS_no_mr_LCI<-sapply(XX, FUN=function(x) 
       plogis(quantile(cjs.T.c.no_misr$BUGSoutput$sims.list$phi.mu +
       x*cjs.T.c.no_misr$BUGSoutput$sims.list$phi.beta, 0.025)))
CJS_no_mr_UCI<-sapply(XX, FUN=function(x)
       plogis(quantile(cjs.T.c.no_misr$BUGSoutput$sims.list$phi.mu +
	   x*cjs.T.c.no_misr$BUGSoutput$sims.list$phi.beta, 0.925)))

CJSm_Median<-sapply(XX, FUN=function(x) 
       plogis(quantile(cjsm.T.c.c$BUGSoutput$sims.list$phi.mu + 
	   x*cjsm.T.c.c$BUGSoutput$sims.list$phi.beta, 0.5)))
CJSm_LCI<-sapply(XX, FUN=function(x) 
       plogis(quantile(cjsm.T.c.c$BUGSoutput$sims.list$phi.mu + 
	   x*cjsm.T.c.c$BUGSoutput$sims.list$phi.beta, 0.025)))
CJSm_UCI<-sapply(XX, FUN=function(x) 
       plogis(quantile(cjsm.T.c.c$BUGSoutput$sims.list$phi.mu +
       x*cjsm.T.c.c$BUGSoutput$sims.list$phi.beta, 0.925)))

Naive_Median<-sapply(XX, FUN=function(x) 
      plogis(quantile(
	  cjs.T.c.with.misidentification$BUGSoutput$sims.list$phi.mu +
	  x*cjs.T.c.with.misidentification$BUGSoutput$sims.list$phi.beta, 
	  0.5)))
Naive_LCI<-sapply(XX, FUN=function(x) 
      plogis(quantile(
	  cjs.T.c.with.misidentification$BUGSoutput$sims.list$phi.mu + 
	  x*cjs.T.c.with.misidentification$BUGSoutput$sims.list$phi.beta,
	  0.025)))
Naive_UCI<-sapply(XX, FUN=function(x) 
      plogis(quantile(
	  cjs.T.c.with.misidentification$BUGSoutput$sims.list$phi.mu +
	  x*cjs.T.c.with.misidentification$BUGSoutput$sims.list$phi.beta,
	  0.925)))
```
#### Plot
```{r}
plot(Naive_Median~XX, ylim=range(c(Naive_LCI, Naive_UCI, CJS_no_mr_LCI, CJS_no_mr_UCI)),
     las=1, col='#ef8a62',  ylab='apparent survial probability',
	 xlab='Time', type='n')

abline(h=seq(0,1, by=0.05), lty=2)
abline(v=seq(0,25, by=5), lty=2)

lines(Naive_Median~XX, lwd=2, col='#d95f02')
lines(Naive_LCI~XX, lty=2, col='#d95f02', lwd=2)
lines(Naive_UCI~XX, lty=2, col='#d95f02', lwd=2)
abline(h=unique(phi), lwd=2, col=grey(0.5))

lines(CJS_no_mr_Median~XX, lty=1, col='#7570b3', lwd=2)
lines(CJS_no_mr_LCI~XX, lty=2, col='#7570b3',lwd=2)
lines(CJS_no_mr_UCI~XX, lty=2, col='#7570b3',lwd=2)

points(CJSm_Median~XX,  bg='#1b9e77', col=NULL, pch=23)
points(CJS_no_mr_LCI~XX,  bg='#1b9e77', col=NULL, pch=24)
points(CJS_no_mr_UCI~XX, bg='#1b9e77', col=NULL,pch=25)
```

## 3 Black-tailed godwits annual survival
Here we analyse capture-resight dataset of black-tailed godwits captured in Friesland, The Netherlands. The raw data consists of three files:
- the main dataset `black-tailed_godwits_CH.csv`;
- `black-tailed_x.Phi.csv` - the matrix for Phi that indicates whether individual [i,] is juvenile `1` or adult `2` at the time [,t]
- `black-tailed_x.P.csv` - the matrix with information on whether individual [i,] at time [,t] is 1st year juvenile `1`, 2nd year juvenile `2`, adult marked as juvenile `3` or adult marked as breeding adult `4`.

### 3.1 CJSm model with black-tailed godwit data

#### Read and prepare the data

```{r}
if (Use_local_data) {
   CH<-read.csv('./data/black-tailed_godwits_CH.csv', check.names=FALSE)
   x.Phi<-read.csv('./data/black-tailed_x.Phi.csv',
                   check.names=FALSE, header=FALSE)
   x.P<-read.csv('./data/black-tailed_x.P.csv',
                   check.names=FALSE, header=FALSE)
} else {
   CH<-as.matrix(read.csv('https://git.io/fhhtJ', check.names=FALSE))
   x.Phi<-as.matrix(read.csv('https://git.io/fje0H', check.names=FALSE, header=FALSE))
   x.P<-as.matrix(read.csv('https://git.io/fje0Q', check.names=FALSE, header=FALSE))
}

# make vector with years of marking
f <- apply(CH, 1, get.first)

# make vector of age at marking (1 for juveniles and 2 for adults)
markedas<-apply(x.Phi, 1, min, na.rm=TRUE)

# Order data by year of marking
Order<-(order(f))
CH<-CH[Order,]
f<-f[Order]
x.P<-x.P[Order,]
x.Phi<-x.Phi[Order,]
markedas<-markedas[Order]

# Make sure that at the year of marking the value is 1 (not more)!
# this is very important for CJSm model
for (i in 1:nrow(CH)) {
   CH[i, f[i]]<-1
}
```
Now, for the CJSm model we need to know how many reads we had each year and how many tags were already used
```{r}
sightings_and_fresh<-colSums(CH) # count per year
N_marked<-rle(apply(CH, 1, get.first))$lengths # marked per year

# add zero in the end if there was no birds marked in the last year
# Note that this part should be corrected if you have missing years with no birds marked
if (length(sightings_and_fresh)>length(N_marked)) N_marked<-c(N_marked, 0)

N_reads<-sightings_and_fresh-N_marked

# delete first year with zero reads.
N_reads<-N_reads[-1]

# estimate how many tags were already used
sum_N_marks_used<-cbind(cumsum(rle(apply(CH, 1, get.first))$lengths),
                           rle(apply(CH, 1, get.first))$values)
sum_N_marks_used<-sum_N_marks_used[,1]
```
Data exploration how many juveniles and adults was marked overall and every year
```{r}
print(table(markedas))
print(table(markedas, f))
```
#### Specify model in BUGS language
```{r}
sink("cjsm-Phi_t_plus_age_plus_year_P_tsm_plus_year_plus_raneff_Theta_c.jags")
cat("
model {
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- phi.t[t, x.Phi[i,t]]
      logit(p_r[i,t]) <- p.age[x.P[i,t]] + p.y[t] +  p.epsilon[i] 
	  theta[i,t] <- theta.t[t]
      } #t
   } #i

for (t in 1:(n.occasions-1)){
   phi.t[t,1] ~ dunif(0, 1)          # Priors for time-spec. survival
   phi.t[t,2] ~ dunif(0, 1)          # Priors for time-spec. survival
}

p.age[1]~ dnorm(0, 0.001)
p.age[2]~ dnorm(0, 0.001)
p.age[3]~ dnorm(0, 0.001)
p.age[4]<-0
   
for (t in 1:(n.occasions-2)){
	  p.y[t] ~ dnorm(0, 0.001)
}
# for the last occasion we fix resighting rate to the mean of previous 5 years
p.y[n.occasions-1]<-mean(p.y[(n.occasions-6):(n.occasions-2)])
	  
for (i in 1:nind){
      p.epsilon[i] ~ dnorm(0, tau)
}

sigma ~ dunif(0, 10) # Prior for standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)
      
mean.theta ~ dunif(0,1)        # Prior for mean correct identification

for (t in 1:(n.occasions-1)) {
   theta.t[t]<-mean.theta
}

for (t in 1:(n.occasions-1)) {
   mu_wrong_assign[t] <-  N_reads[t] *(1-theta.t[t]) / (sum_N_marks_used[t]-1)
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
	  lambda_obs[i,t]<- -log(1-p_r[i,t-1])*z[i,t]*theta[i,t-1] +
	                         mu_wrong_assign[t-1]
	  y[i,t] ~ dpois(lambda_obs[i,t])
      } #t
   } #i
}
",fill = TRUE)
sink()
```
#### Bundle data, generate intitials and prepare model run
```{r}
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],
                  N_reads=N_reads, sum_N_marks_used=sum_N_marks_used, x.Phi=x.Phi,x.P=x.P)

# Initial values
inits <- function(){list(z = known.state.cjs.mult(CH), phi.t = 
     cbind(runif((dim(CH)[2]-1), 0, 1),
	      runif((dim(CH)[2]-1), 0, 1)),
		  p.age =c(rnorm(3, 0, 5), NA),
		  p.y= c(rnorm((dim(CH)[2]-2), 0, 5), NA),
		  mean.theta=runif(1, 0.5, 1),
		  sigma = runif(1, 0, 5))}

# Parameters monitored
parameters <- c("phi.t", "p.age", "mean.theta", "sigma2", "p.y")

# MCMC settings
ni <- 50000
nt <- 1
nb <- 10000
nc <- 6
```
#### Run the model
```{r}
if (Run_everything) {
a<-Sys.time()
   # (BRT ~ 2.5 hours)
   cjsm.godwits <- jags.parallel(jags.data, inits, parameters,
        "cjsm-Phi_t_plus_age_plus_year_P_tsm_plus_year_plus_raneff_Theta_c.jags",
		n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,  
		export_obj_names=c('cjs.init.z', 'known.state.cjs.mult', 'nb', 'ni', 
		                   'nt', 'CH', 'f', 'x.P', 'x.Phi', 'markedas'))
Sys.time()-a
   saveRDS(cjsm.godwits, file='./results/cjsm.godwits.RDS')
}

cjsm.godwits<-readRDS(file='./results/cjsm.godwits.RDS')
print(cjsm.godwits)
```

### 3.2 Naive CJS model with black-tailed godwit data
#### Prepare the data
The only difference from CJSm model here is that we flatten the data replacing all positive values with 1.
```{r}
# flatten data
CH_flat<-apply(CH, c(1,2), sign)
```

#### Specify model in BUGS language
```{r}
sink("cjs-naive-Phi_t_plus_age_plus_year_P_tsm_plus_year_plus_raneff.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- phi.t[t, x.Phi[i,t]]
      logit(p_r[i,t]) <- p.age[x.P[i,t]] + p.y[t] +  p.epsilon[i] 
      } #t
   } #i

for (t in 1:(n.occasions-1)){
   phi.t[t,1] ~ dunif(0, 1)          # Priors for time-spec. survival
   phi.t[t,2] ~ dunif(0, 1)          # Priors for time-spec. survival
}

p.age[1]~ dnorm(0, 0.001)
p.age[2]~ dnorm(0, 0.001)
p.age[3]~ dnorm(0, 0.001)
p.age[4]<-0
   
for (t in 1:(n.occasions-2)){
	  p.y[t] ~ dnorm(0, 0.001)
}
	  p.y[n.occasions-1]<-mean(p.y[(n.occasions-6):(n.occasions-2)])
	  
for (i in 1:nind){
      p.epsilon[i] ~ dnorm(0, tau)
}


sigma ~ dunif(0, 10) # Prior for standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)
  
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
      mu2[i,t] <- p_r[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()
```

#### Bundle data, generate intitials and prepare model run
```{r}
# Bundle data
jags.data <- list(y = CH_flat, f = f, nind = dim(CH_flat)[1],
                  n.occasions = dim(CH_flat)[2], x.Phi=x.Phi,x.P=x.P)


# Initial values
inits <- function(){list(z = known.state.cjs(CH_flat), phi.t = 
     cbind(runif((dim(CH_flat)[2]-1), 0, 1),
	      runif((dim(CH_flat)[2]-1), 0, 1)),
		  p.age =c(rnorm(3, 0, 5), NA),
		  p.y= c(rnorm((dim(CH_flat)[2]-2), 0, 5), NA),
		  mean.theta=runif(1, 0.5, 1),
		  sigma = runif(1, 0, 5))}

# Parameters monitored
parameters <- c("phi.t", "p.age", "sigma2", "p.y")

# MCMC settings
ni <- 25000
nt <- 1
nb <- 5000
nc <- 6
```

#### Run the model
```{r}
if (Run_everything) {
   # (BRT ~ 2.5 hours)
  cjs.naive.godwits <- jags.parallel(jags.data, inits, parameters, 
   "cjs-naive-Phi_t_plus_age_plus_year_P_tsm_plus_year_plus_raneff.jags",
   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,  
   export_obj_names=
   c('known.state.cjs', 'nb', 'ni', 'nt', 'CH_flat', 'f', 'x.Phi', 'x.P'))
  saveRDS(cjs.naive.godwits, file='./results/cjs.naive.godwits.RDS')
}

cjs.naive.godwits<-readRDS(file='./results/cjs.naive.godwits.RDS')
print(cjs.naive.godwits)
```

### 3.3 Plot godwits survival estimates

#### Extract output
```{r}
Res_adults<-data.frame(
   Year.Phi=as.numeric(dimnames(CH)[[2]])[-length(dimnames(CH)[[2]])],
   Phi.cjsm.median= apply(cjsm.godwits$BUGSoutput$sims.list$phi.t[,,2],
                          2,quantile, probs=0.5),
   Phi.cjsm.lci= apply(cjsm.godwits$BUGSoutput$sims.list$phi.t[,,2],
                          2,quantile, probs=0.025),
   Phi.cjsm.uci= apply(cjsm.godwits$BUGSoutput$sims.list$phi.t[,,2],
                          2,quantile, probs=0.975),
   Phi.cjsn.median= apply(cjs.naive.godwits$BUGSoutput$sims.list$phi.t[,,2],
                          2,quantile, probs=0.5),
   Phi.cjsn.lci= apply(cjs.naive.godwits$BUGSoutput$sims.list$phi.t[,,2],
                          2,quantile, probs=0.025),
   Phi.cjsn.uci= apply(cjs.naive.godwits$BUGSoutput$sims.list$phi.t[,,2],
                          2,quantile, probs=0.975),
   Year.p=as.numeric(dimnames(CH)[[2]])[-1],
   p.cjsm.median= plogis(apply(cjsm.godwits$BUGSoutput$sims.list$p.y,
                          2,quantile, probs=0.5)),
   p.cjsm.lci= plogis(apply(cjsm.godwits$BUGSoutput$sims.list$p.y, 
                          2,quantile, probs=0.025)),
   p.cjsm.uci= plogis(apply(cjsm.godwits$BUGSoutput$sims.list$p.y, 
                          2,quantile, probs=0.975)),
   p.cjsn.median= plogis(apply(cjs.naive.godwits$BUGSoutput$sims.list$p.y,
                          2,quantile, probs=0.5)),
   p.cjsn.lci= plogis(apply(cjs.naive.godwits$BUGSoutput$sims.list$p.y, 
                          2,quantile, probs=0.025)),
   p.cjsn.uci= plogis(apply(cjs.naive.godwits$BUGSoutput$sims.list$p.y,
                          2,quantile, probs=0.975)))

Res_juveniles<-data.frame(
   Year.Phi=as.numeric(dimnames(CH)[[2]])[-length(dimnames(CH)[[2]])],
   Phi.cjsm.median= apply(cjsm.godwits$BUGSoutput$sims.list$phi.t[,,1],
                          2,quantile, probs=0.5),
   Phi.cjsm.lci= apply(cjsm.godwits$BUGSoutput$sims.list$phi.t[,,1],
                          2,quantile, probs=0.025),
   Phi.cjsm.uci= apply(cjsm.godwits$BUGSoutput$sims.list$phi.t[,,1], 
                          2,quantile, probs=0.975),
   Phi.cjsn.median= apply(cjs.naive.godwits$BUGSoutput$sims.list$phi.t[,,1],
                          2,quantile, probs=0.5),
   Phi.cjsn.lci= apply(cjs.naive.godwits$BUGSoutput$sims.list$phi.t[,,1],
                          2,quantile, probs=0.025),
   Phi.cjsn.uci= apply(cjs.naive.godwits$BUGSoutput$sims.list$phi.t[,,1],
                          2,quantile, probs=0.975))
```
#### Plot

```{r}
#pdf('godwit_figure_2.pdf', 10,5, useDingbats=FALSE)
Shift=0.1  
par(mfrow=c(1,2), mar=c(5.1, 4.1, 0.5, 0.5))
  
#Panel A juveniles

plot(Res_juveniles$Phi.cjsm.median~Res_juveniles$Year.Phi,
     ylim=c(0, 1), las=1, pch=19, ylab='apparent survival',
	 xlab='Year', type='n')  

abline(v=2004:2020, lty=3, col=grey(0.5))
abline(v=seq(2005, 2020, by=5), lty=2, col=grey(0.5))

abline(h=seq(0, 1, by=0.1), lty=3, col=grey(0.5))
abline(h=seq(0, 1, by=0.5), lty=2, col=grey(0.5))

points(Res_juveniles$Phi.cjsm.median~I(Res_juveniles$Year.Phi+Shift), 
       ylim=c(0.6, 1), las=1, pch=23, bg='#1b9e77', col='#1b9e77')   
segments(y0=Res_juveniles$Phi.cjsm.lci,
         y1=Res_juveniles$Phi.cjsm.uci,
		 x0=Res_juveniles$Year.Phi+Shift, lwd=2, col='#1b9e77')

points(Res_juveniles$Phi.cjsn.median~I(Res_juveniles$Year.Phi-Shift),
       ylim=c(0.6, 1), las=1, pch=19, col='#d95f02')   
segments(y0=Res_juveniles$Phi.cjsn.lci,
         y1=Res_juveniles$Phi.cjsn.uci,
		 x0=Res_juveniles$Year.Phi-Shift, lwd=2, col='#d95f02')

#Panel B adults
plot(Res_adults$Phi.cjsm.median~Res_adults$Year.Phi,
     ylim=c(0.75, 1), las=1, ylab='',
	 xlab='Year', type='n')  

abline(v=2004:2020, lty=3, col=grey(0.5))
abline(v=seq(2005, 2020, by=5), lty=2, col=grey(0.5))

abline(h=seq(0.5, 1, by=0.02), lty=3, col=grey(0.5))
abline(h=seq(0.5, 1, by=0.1), lty=2, col=grey(0.5))

points(Res_adults$Phi.cjsm.median~I(Res_adults$Year.Phi+Shift), 
       ylim=c(0.6, 1), las=1, pch=23, bg='#1b9e77', col='#1b9e77')   

segments(y0=Res_adults$Phi.cjsm.lci,
         y1=Res_adults$Phi.cjsm.uci,
		 x0=Res_adults$Year.Phi+Shift, lwd=2, col='#1b9e77')

points(Res_adults$Phi.cjsn.median~I(Res_adults$Year.Phi-Shift),
       ylim=c(0.6, 1), las=1, pch=19 , col='#d95f02')   
segments(y0=Res_adults$Phi.cjsn.lci,
         y1=Res_adults$Phi.cjsn.uci,
		 x0=Res_adults$Year.Phi-Shift, lwd=2, col='#d95f02')
```

### 3.4 Summary statistics used in the paper
#### CJSm model results

```{r}
# average survival probability of adults
print(quantile(cjsm.godwits$BUGSoutput$sims.list$phi.t[,,2],
         c(0.025, 0.5, 0.975)), digits=3)

# average survival probability of juveniles
print(quantile(cjsm.godwits$BUGSoutput$sims.list$phi.t[,,1],
         c(0.025, 0.5, 0.975)), digits=3)

# average adult resighting probability
print(plogis(quantile(cjsm.godwits$BUGSoutput$sims.list$p.y,
       c(0.025, 0.5, 0.975))), digits=5)

# average expected number of sightings	   
print(-log(1-plogis(quantile(cjsm.godwits$BUGSoutput$sims.list$p.y,
       c(0.025, 0.5, 0.975)))), digits=3)

# Average 1st year juvenile resighting rate
print(plogis(quantile(cjsm.godwits$BUGSoutput$sims.list$p.y+
             cjsm.godwits$BUGSoutput$sims.list$p.age[,1],
			 c(0.025, 0.5, 0.975))), digits=3)

# Average 2nd year juvenile resighting rate
print(plogis(quantile(cjsm.godwits$BUGSoutput$sims.list$p.y+
             cjsm.godwits$BUGSoutput$sims.list$p.age[,2],
			 c(0.025, 0.5, 0.975))), digits=3)

# Average 2nd year juvenile resighting rate
print(plogis(quantile(cjsm.godwits$BUGSoutput$sims.list$p.y+
             cjsm.godwits$BUGSoutput$sims.list$p.age[,3],
			 c(0.025, 0.5, 0.975))), digits=3)

# sigma - RE on resighting prob.
print(sqrt(quantile(cjsm.godwits$BUGSoutput$sims.list$sigma2,
           c(0.025, 0.5, 0.975))), digits=3)

# and theta
print(quantile(cjsm.godwits$BUGSoutput$sims.list$mean.theta,
           c(0.025, 0.5, 0.975)), digits=3)

print(sd(cjsm.godwits$BUGSoutput$sims.list$mean.theta))
```

#### naive CJS model results

```{r}
# average survival probability of adults
print(quantile(cjs.naive.godwits$BUGSoutput$sims.list$phi.t[,,2],
         c(0.025, 0.5, 0.975)), digits=3)

# average survival probability of juveniles
print(quantile(cjs.naive.godwits$BUGSoutput$sims.list$phi.t[,,1],
         c(0.025, 0.5, 0.975)), digits=3)

# average adult resighting probability
print(plogis(quantile(cjs.naive.godwits$BUGSoutput$sims.list$p.y,
       c(0.025, 0.5, 0.975))), digits=5)

# average expected number of sightings	   
print(-log(1-plogis(quantile(cjs.naive.godwits$BUGSoutput$sims.list$p.y,
       c(0.025, 0.5, 0.975)))), digits=3)

# Average 1st year juvenile resighting rate
print(plogis(quantile(cjs.naive.godwits$BUGSoutput$sims.list$p.y+
             cjs.naive.godwits$BUGSoutput$sims.list$p.age[,1],
			 c(0.025, 0.5, 0.975))), digits=3)

# Average 2nd year juvenile resighting rate
print(plogis(quantile(cjs.naive.godwits$BUGSoutput$sims.list$p.y+
             cjs.naive.godwits$BUGSoutput$sims.list$p.age[,2],
			 c(0.025, 0.5, 0.975))), digits=3)

# Average 2nd year juvenile resighting rate
print(plogis(quantile(cjs.naive.godwits$BUGSoutput$sims.list$p.y+
             cjs.naive.godwits$BUGSoutput$sims.list$p.age[,3],
			 c(0.025, 0.5, 0.975))), digits=3)

# sigma - RE on resighting prob.
print(sqrt(quantile(cjs.naive.godwits$BUGSoutput$sims.list$sigma2,
           c(0.025, 0.5, 0.975))), digits=3)
```

P.S. The html file from this markdown can be recereated with folloing code
```{r eval=FALSE}
download.file(
   'https://raw.githubusercontent.com/eldarrak/CJS-with-misidentification/master/code/All_code.md',
   'tmp.rmd', cacheOK = FALSE)
rmarkdown::render('tmp.rmd', output_format = 'html_document',
        output_options=list(toc=TRUE, toc_float=list(collapsed=FALSE)), 
        encoding='utf-8')
```

