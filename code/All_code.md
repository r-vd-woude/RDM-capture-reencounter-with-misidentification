# Supplementary code to 'Misidentification errors in reencounters result in biased estimates of survival from CJS models: evidence and a possible solution using the robust design' by E. Rakhimberdiev et al., Methods in Ecology and Evolution, submitted

_code by Eldar Rakhimberdiev_

You might enjoy [html version of this supplementary](http://htmlpreview.github.io/?https://raw.githubusercontent.com/eldarrak/RDM-capture-reencounter-with-misidentification/master/code/All_code.html).

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
library(httr)
library(tidyverse)
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
   Result_req <- GET("https://api.github.com/repos/eldarrak/RDM-capture-reencounter-with-misidentification/contents/results")
   Results_files <- unlist(lapply(content(Result_req), "[", "name"), use.names = F)
   Results_dir<-'https://github.com/eldarrak/CJS-with-misidentification/blob/master/results/'
   for (i in 1:length(Results_files)) {
      Filename<-Results_files[i]
      download.file(paste0(Results_dir, Filename, '?raw=true'),
        destfile=paste0('./results/', Filename), mode='wb')
   }
  }
if (Use_local_data & !'data' %in% list.dirs(full.names=FALSE, recursive=FALSE)) {
    stop('data directory not found') 
}
```

### 0.4 load libraries and source CJSm functions
We need `jagsUI` to run jags code. Note that to run the models you will also have to install jags software (Plummer 2011) from [here]( http://mcmc-jags.sourceforge.net).

```{r}
library('jagsUI')
```

Now we source helper functions
```{r}
# general functions
source('https://git.io/fhhtZ')
```

## 1. Simulate and solve Phi_dot_P_dot_Theta_dot models
For this exersice we simulate data with the function `simul.cjs.multiple.sightings` and then with the function `run_CJSm_c_c_c` estimate two cjs models over these data Phi_dot_P_dot model (CJS-c-c) and Phi_dot_P_dot_Theta_dot (cjsm-c-c-c)

### 1.1 Define wrapper function to run simulations and run the models
```{r}
run_c_c_c_models<-function(CH, n.chains=6,  n.adapt=NULL, RHat_limit=1.1, max.iter=200000, 
                           iter.increment=5000, models.to.run=c('CJS', 'RDM', 'RDMa')) {
  # this function runs Phi_dot_p_dot_Theta_dot models
  Res<-list(CH=CH)
  
  if ('CJS' %in% models.to.run) {
     cat('running CJS-c-c model')
     # flatten the data
     CH_flat<-apply(CH, c(1,2), sign)
     # Create vector with occasion of marking
     f <- apply(CH_flat, 1, get.first)
     Order<-(order(f))
     CH_flat<-CH_flat[Order,]
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
    jags.data <- list(y = CH_flat, f = f, nind = dim(CH_flat)[1], n.occasions = dim(CH_flat)[2])

    inits <- function(){list(mean.phi=runif(1,0,1), mean.p = runif(1,0,1), z = known.state.cjs(CH_flat))}

    # Parameters monitored
    parameters <- c("mean.p", "mean.phi")

    # Call JAGS from R 
    cjs.c.c <- autojags(jags.data, inits, 
                       parameters, "cjs-c-c.jags", n.chains = n.chains, 
                       parallel=TRUE, 
                       Rhat.limit=RHat_limit, n.adapt=NULL, 
                       max.iter=max.iter, iter.increment= iter.increment, save.all.iter=FALSE )  
    cat('   Done!\n')
    Res$CJS<-cjs.c.c
  }
  if ('RDM' %in% models.to.run | 'RDMa' %in% models.to.run) {
    # Create vector with occasion of marking
    f <- apply(CH, 1, get.first)
    Order<-(order(f))
    CH<-CH[Order,]
    f<-f[Order]
	
	sightings_and_fresh<-colSums(CH) # count per year
    All.Years=data.frame(Years=1:dim(CH)[2], N_marked=0)
    Rle<-rle(apply(CH, 1, get.first))
    for (year in 1:length(Rle$lengths)) {
      All.Years$N_marked[Rle$values[year]]<-Rle$lengths[year]
    }

    N_marked<-  All.Years$N_marked
    N_reads<-sightings_and_fresh-N_marked
    N_reads<-N_reads[-1]
    sum_N_marks_used<-cumsum(All.Years$N_marked)

    # Bundle data
    jags.data <- list(y = CH, f = f, nind = dim(CH)[1],
                     n.occasions = dim(CH)[2],
                     sum_N_marks_used=sum_N_marks_used,
                     N_reads=N_reads, N_reads2=N_reads)

    inits <- function(){list(z = known.state.cjs.mult(CH), 
                            mean.phi = runif(1, 0, 1),
                            mean.p = runif(1, 0, 1),
                            mean.theta=runif(1, 0, 1))}

    # Parameters monitored
    parameters <- c("mean.phi" , "mean.p", "mean.theta")
  } 
  
   if ('RDM' %in% models.to.run) {
    # Specify model in BUGS language
    sink("RDM-c-c-c.jags")
    cat("
      model {
        # Priors and constraints
        for (i in 1:nind){
          for (t in f[i]:(n.occasions-1)){
             phi[i,t] <- mean.phi
             lambda_true[i,t]<-mean.lambda.true
          } #t
        } #i
        mean.phi ~ dunif(0, 1)         # Prior for mean survival
        mean.p ~ dunif(0, 1)           # Prior for mean recapture
        mean.lambda.true<- -log(1-mean.p)
        mean.theta ~ dunif(0,1)        # Prior for mean correct identification
        for (t in 1:(n.occasions-1)) {
          theta.t[t]<-mean.theta
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
            lambda_obs[i,t-1]<- lambda_true[i,t-1]*z[i,t]*theta.t[t-1]+
		                    Lambda_z_sum[t-1] -
  						    lambda_true[i,t-1]*z[i,t]*(1-theta.t[t-1])/(sum_N_marks_used[t-1]-1)
          } #t
        } #i
	
        for (t in 1:(n.occasions-1)) {
	    # here I add the total sum over all lambda_z      
	    Lambda_z_sum[t]<-sum(lambda_true[1:(sum_N_marks_used[t]),t]*
	               z[1:(sum_N_marks_used[t]),t+1])*
	    		   (1-theta.t[t])/(sum_N_marks_used[t]-1)
        # Below is the same but in RDMa model
		#Lambda_z_sum[t]<-N_reads2[t]*
	   	#	  (1-theta.t[t])/(sum_N_marks_used[t]-1)		
        Lambda_obs_sum[t]<-sum(lambda_obs[1:(sum_N_marks_used[t]),t])
        P_multin[1:(sum_N_marks_used[t]),t] <- 
		      lambda_obs[1:(sum_N_marks_used[t]),t]/
			  ifelse(Lambda_obs_sum[t]==0, 1, Lambda_obs_sum[t])
        N_reads[t]~ dpois(Lambda_obs_sum[t])
        y[1:sum_N_marks_used[t],t+1] ~ 
		      dmulti(P_multin[1:sum_N_marks_used[t],t],
	          N_reads2[t]) 
        }
    }
    ",fill = TRUE)
    sink()
   
    # Call JAGS from R 
    RDM.c.c.c <- autojags(jags.data, inits, parameters,
                       "RDM-c-c-c.jags", n.chains = n.chains, 
                       parallel=TRUE, 
                       Rhat.limit=RHat_limit, n.adapt=NULL, 
                       max.iter=max.iter, iter.increment= iter.increment, save.all.iter=FALSE)
					   cat('   Done!\n')
    Res$RDM<-RDM.c.c.c
  }

  if ('RDMa' %in% models.to.run) {
    cat('running RDMa-c-c-c model')
    # Specify model in BUGS language
    sink("RDMa-c-c-c.jags")
    cat("
      model {
        # Priors and constraints
        for (i in 1:nind){
          for (t in f[i]:(n.occasions-1)){
             phi[i,t] <- mean.phi
             lambda_true[i,t]<-mean.lambda.true
          } #t
        } #i
        mean.phi ~ dunif(0, 1)         # Prior for mean survival
        mean.p ~ dunif(0, 1)           # Prior for mean recapture
        mean.lambda.true<- -log(1-mean.p)
        mean.theta ~ dunif(0,1)        # Prior for mean correct identification
        for (t in 1:(n.occasions-1)) {
          theta.t[t]<-mean.theta
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
            lambda_obs[i,t-1]<- lambda_true[i,t-1]*z[i,t]*theta.t[t-1]+
		                   Lambda_z_sum[t-1] -	
   						   lambda_true[i,t-1]*z[i,t]*(1-theta.t[t-1])/(sum_N_marks_used[t-1]-1)
          } #t
        } #i
	
        for (t in 1:(n.occasions-1)) {
	    # here I add the total sum over all lambda_z      
	    #Lambda_z_sum[t]<-sum(lambda_true[1:(sum_N_marks_used[t]),t]*
	    #           z[1:(sum_N_marks_used[t]),t+1])*
	    #		   (1-theta.t[t])/(sum_N_marks_used[t]-1)
        # above is the RDM model
		Lambda_z_sum[t]<-N_reads2[t]*
	   		  (1-theta.t[t])/(sum_N_marks_used[t]-1)		
        Lambda_obs_sum[t]<-sum(lambda_obs[1:(sum_N_marks_used[t]),t])
        P_multin[1:(sum_N_marks_used[t]),t] <- 
		      lambda_obs[1:(sum_N_marks_used[t]),t]/
			  ifelse(Lambda_obs_sum[t]==0, 1, Lambda_obs_sum[t])
        N_reads[t]~ dpois(Lambda_obs_sum[t])
        y[1:sum_N_marks_used[t],t+1] ~ 
		      dmulti(P_multin[1:sum_N_marks_used[t],t],
	          N_reads2[t]) 
        }
    }
    ",fill = TRUE)
    sink()

    # Call JAGS from R 
    RDMa.c.c.c <- autojags(jags.data, inits, parameters,
                       "RDMa-c-c-c.jags", n.chains = n.chains, 
                       parallel=TRUE, 
                       Rhat.limit=RHat_limit, n.adapt=NULL, 
                       max.iter=max.iter, iter.increment= iter.increment, save.all.iter=FALSE)
					   cat('   Done!\n')
    Res$RDMa<-RDMa.c.c.c
   }
   return(Res)
}
```

### 1.2 Define parameters for a sample simulation (Phi=0.9, P=0.9, theta=0.95)

```{r}
n.occasions <- 5                       # Number of capture occasions
n.marked <- 15                         # Number marked each year
marked <- rep(n.marked, n.occasions-1) # Annual number of newly marked individuals

Phi_dot=0.9
P_dot=0.9
Theta_dot=0.95
RHat_limit=1.01
max.iter=20000
n.adapt=5000
iter.increment=4000
nRuns=5

phi <- rep(Phi_dot, n.occasions-1) # survival
p <- rep(P_dot, n.occasions-1)  # resighting
# Define matrices with survival recapture and correct identification
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

theta<- rep(Theta_dot, n.occasions-1) # correctly identified
THETA<-matrix(theta, ncol = n.occasions-1, nrow = sum(marked))
```

### 1.3 Run some simulations
Note, if you turn `Run_everything` to `TRUE`  the loop below will try to run 200 models! and it will take time. 

```{r}
if (Run_everything) {

Output_df<-data.frame(Model=character(), 
            Phi.true=numeric(), P.true=numeric(), N.marked = integer(), Theta.true=numeric(),
            Phi.mean=numeric(), P.mean =numeric(), Theta.mean=numeric(), 
			Deviance =numeric(), 
			Phi.CI.width =numeric(),P.CI.width=numeric(), Theta.CI.width=numeric(),
			Phi.CI.coverage = logical(), P.CI.coverage = logical(), Theta.CI.coverage = logical(),
			Phi.n.eff=integer(), P.n.eff=integer(), Theta.n.eff=integer(), Deviance.n.eff=integer(),
			Phi.RHat =numeric(), P.RHat =numeric(), Theta.RHat =numeric(), Deviance.RHat =numeric(),
			n.chains =integer(), n.burnin=integer(),  n.iter=integer(), n.elapsed.mins  = integer(), 
			Converged = logical(), Remark = character(),
			Phi.bias = numeric(),P.bias = numeric(), Theta.bias = numeric())

						   
for (i in 1:nRuns) {
   cat('i = ', i, '\n')
   # simulate data
   CH_list <- simul.cjs.multiple.sightings(PHI, P, THETA,  marked)
   # estimate models
   Res_list<-run_c_c_c_models(CH_list$CH, n.chains=6,  
                              n.adapt=n.adapt, RHat_limit=RHat_limit,
							  max.iter=max.iter, iter.increment=iter.increment,
							models.to.run=c('CJS', 'RDM', 'RDMa'))

   # extract results					
    if (!is.null(Res_list$CJS)) {
	    CI_cur.CJS<-cbind(c('mean.phi'=phi[1], 'mean.p'=p[1], 
                    'mean.theta'=theta[1]), 
                  c(Res_list$CJS$q2.5$mean.phi, 
				    Res_list$CJS$q2.5$mean.p,
					NA), 
				  c(Res_list$CJS$q97.5$mean.phi, 
				    Res_list$CJS$q97.5$mean.p,
					NA)) 
    	Output_line_CJS<-data.frame(Model='CJS', 
            Phi.true=Phi_dot, P.true=P_dot, N.marked = n.marked, Theta.true=Theta_dot,
            Phi.mean=Res_list$CJS$mean$mean.phi, P.mean =Res_list$CJS$mean$mean.p, Theta.mean=NA, 
			Deviance =Res_list$CJS$mean$deviance, 
			Phi.CI.width =CI_cur.CJS[1,3]-CI_cur.CJS[1,2],
			P.CI.width=CI_cur.CJS[2,3]-CI_cur.CJS[2,2], Theta.CI.width=NA,
			Phi.CI.coverage =  CI_cur.CJS[1,1]>CI_cur.CJS[1,2] & CI_cur.CJS[1,1]<=CI_cur.CJS[1,3], 
			P.CI.coverage =    CI_cur.CJS[2,1]>CI_cur.CJS[2,2] & CI_cur.CJS[2,1]<=CI_cur.CJS[2,3],
			Theta.CI.coverage = NA,
			Phi.n.eff=unlist(Res_list$CJS$n.eff)[1], 
			P.n.eff=unlist(Res_list$CJS$n.eff)[2], 
			Theta.n.eff=NA, 
			Deviance.n.eff=unlist(Res_list$CJS$n.eff)[3],
			Phi.RHat = unlist(Res_list$CJS$Rhat)[1], P.RHat =unlist(Res_list$CJS$Rhat)[2],
			Theta.RHat =NA, Deviance.RHat =unlist(Res_list$CJS$Rhat)[3],
			n.chains =Res_list$CJS$mcmc.info$n.chains, 
			n.burnin=Res_list$CJS$mcmc.info$n.burnin,
			n.iter=Res_list$CJS$mcmc.info$n.iter, 
			n.elapsed.mins  = Res_list$CJS$mcmc.info$elapsed.mins, 
			Converged = max(unlist(Res_list$CJS$Rhat), na.rm=TRUE)<1.01,
			Remark = '',
			Phi.bias = Res_list$CJS$mean$mean.phi-Phi_dot,
			P.bias = Res_list$CJS$mean$mean.p-P_dot,
			Theta.bias = NA)
	Output_df<-rbind(Output_df, Output_line_CJS)		
    }
	
    if (!is.null(Res_list$RDM)) {	
    	CI_cur.RDM<-cbind(c('mean.phi'=phi[1], 'mean.p'=p[1], 
                    'mean.theta'=theta[1]), 
                  c(Res_list$RDM$q2.5$mean.phi, 
				    Res_list$RDM$q2.5$mean.p,
					Res_list$RDM$q2.5$mean.theta), 
				  c(Res_list$RDM$q97.5$mean.phi, 
				    Res_list$RDM$q97.5$mean.p,
					Res_list$RDM$q97.5$mean.theta))
    	Output_line_RDM<-data.frame(Model='RDM', 
            Phi.true=Phi_dot, P.true=P_dot, N.marked = n.marked, Theta.true=Theta_dot,
            Phi.mean=Res_list$RDM$mean$mean.phi, P.mean =Res_list$RDM$mean$mean.p, Theta.mean=Res_list$RDM$mean$mean.theta, 
			Deviance =Res_list$RDM$mean$deviance, 
			Phi.CI.width =CI_cur.RDM[1,3]-CI_cur.RDM[1,2],
			P.CI.width=CI_cur.RDM[2,3]-CI_cur.RDM[2,2], 
			Theta.CI.width=CI_cur.RDM[3,3]-CI_cur.RDM[3,2],
			Phi.CI.coverage =  CI_cur.RDM[1,1]>CI_cur.RDM[1,2] & CI_cur.RDM[1,1]<=CI_cur.RDM[1,3], 
			P.CI.coverage =    CI_cur.RDM[2,1]>CI_cur.RDM[2,2] & CI_cur.RDM[2,1]<=CI_cur.RDM[2,3],
			Theta.CI.coverage =  CI_cur.RDM[3,1]>CI_cur.RDM[3,2] & CI_cur.RDM[2,1]<=CI_cur.RDM[3,3],
			Phi.n.eff=unlist(Res_list$RDM$n.eff)[1], 
			P.n.eff=unlist(Res_list$RDM$n.eff)[2], 
			Theta.n.eff=unlist(Res_list$RDM$n.eff)[3], 
			Deviance.n.eff=unlist(Res_list$RDM$n.eff)[4],
			Phi.RHat = unlist(Res_list$RDM$Rhat)[1], P.RHat =unlist(Res_list$RDM$Rhat)[2],
			Theta.RHat =unlist(Res_list$RDM$Rhat)[3], Deviance.RHat =unlist(Res_list$RDM$Rhat)[4],
			n.chains =Res_list$RDM$mcmc.info$n.chains, 
			n.burnin=Res_list$RDM$mcmc.info$n.burnin,
			n.iter=Res_list$RDM$mcmc.info$n.iter, 
			n.elapsed.mins  = Res_list$RDM$mcmc.info$elapsed.mins, 
			Converged = max(unlist(Res_list$RDM$Rhat), na.rm=TRUE)<1.01,
			Remark = '',
			Phi.bias = Res_list$RDM$mean$mean.phi-Phi_dot,
			P.bias = Res_list$RDM$mean$mean.p-P_dot,
			Theta.bias =  Res_list$RDM$mean$mean.theta-Theta_dot)
		
		Output_df<-rbind(Output_df, Output_line_RDM)		
		}
		
    if (!is.null(Res_list$RDMa)) {
     CI_cur.RDMa<-cbind(c('mean.phi'=phi[1], 'mean.p'=p[1], 
                    'mean.theta'=theta[1]), 
                  c(Res_list$RDMa$q2.5$mean.phi, 
				    Res_list$RDMa$q2.5$mean.p,
					Res_list$RDMa$q2.5$mean.theta), 
				  c(Res_list$RDMa$q97.5$mean.phi, 
				    Res_list$RDMa$q97.5$mean.p,
					Res_list$RDMa$q97.5$mean.theta))

  	Output_line_RDMa<-data.frame(Model='RDMa', 
            Phi.true=Phi_dot, P.true=P_dot, N.marked = n.marked, Theta.true=Theta_dot,
            Phi.mean=Res_list$RDMa$mean$mean.phi, P.mean =Res_list$RDMa$mean$mean.p, Theta.mean=Res_list$RDMa$mean$mean.theta, 
			Deviance =Res_list$RDMa$mean$deviance, 
			Phi.CI.width =CI_cur.RDMa[1,3]-CI_cur.RDMa[1,2],
			P.CI.width=CI_cur.RDMa[2,3]-CI_cur.RDMa[2,2], 
			Theta.CI.width=CI_cur.RDMa[3,3]-CI_cur.RDMa[3,2],
			Phi.CI.coverage =  CI_cur.RDMa[1,1]>CI_cur.RDMa[1,2] & CI_cur.RDMa[1,1]<=CI_cur.RDMa[1,3], 
			P.CI.coverage =    CI_cur.RDMa[2,1]>CI_cur.RDMa[2,2] & CI_cur.RDMa[2,1]<=CI_cur.RDMa[2,3],
			Theta.CI.coverage =  CI_cur.RDMa[3,1]>CI_cur.RDMa[3,2] & CI_cur.RDMa[2,1]<=CI_cur.RDMa[3,3],
			Phi.n.eff=unlist(Res_list$RDMa$n.eff)[1], 
			P.n.eff=unlist(Res_list$RDMa$n.eff)[2], 
			Theta.n.eff=unlist(Res_list$RDMa$n.eff)[3], 
			Deviance.n.eff=unlist(Res_list$RDMa$n.eff)[4],
			Phi.RHat = unlist(Res_list$RDMa$Rhat)[1], P.RHat =unlist(Res_list$RDMa$Rhat)[2],
			Theta.RHat =unlist(Res_list$RDMa$Rhat)[3], Deviance.RHat =unlist(Res_list$RDMa$Rhat)[4],
			n.chains =Res_list$RDMa$mcmc.info$n.chains, 
			n.burnin=Res_list$RDMa$mcmc.info$n.burnin,
			n.iter=Res_list$RDMa$mcmc.info$n.iter, 
			n.elapsed.mins  = Res_list$RDMa$mcmc.info$elapsed.mins, 
			Converged = max(unlist(Res_list$RDMa$Rhat), na.rm=TRUE)<1.01,
			Remark = '',
			Phi.bias = Res_list$RDMa$mean$mean.phi-Phi_dot,
			P.bias = Res_list$RDMa$mean$mean.p-P_dot,
			Theta.bias =  Res_list$RDMa$mean$mean.theta-Theta_dot)
	
    	Output_df<-rbind(Output_df, Output_line_RDM)		
	}	
	
   ###
   # once in e.g. 10 runs we will save the result..
   if (i%%10 ==0) { 
     save(Output_df,
		file=paste0('Output_df_Phi_dot_', Phi_dot, '_P_dot_', P_dot, 
            '_theta_dot_', Theta_dot, '_n_occasions_', n.occasions, 
	        '_marked_', n.marked, '_iteration_', i, '_RHat_', RHat_limit, '.Rdata'))
    }
   }
			 
   Filename<-paste0('Output_df_fin_', nRuns, '_simulations_Phi_dot_', 
      Phi_dot, '_p_dot_', P_dot, 
      '_Theta_dot_', Theta_dot,
	  '_n_occasions_',n.occasions, '_marked_', n.marked, '_RHat_', RHat_limit, '.RData')	
	  
   save(Output_df, file=Filename)
}

```

### 1.4 Make Figure 1

load packages and data
```{r}
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(ggsci)
library(stringr)

out_all<-readRDS('./results/Sim_results_25_years.RDS')
```

prepare data for plotting

```{r}

Res <- out_all %>% 
       rename(Phi.estimate = Phi.mean, P.estimate = P.mean, Theta.estimate = Theta.mean) %>%
       filter(Converged & N.marked !=100 & !(Model=='CJS_naive' & is.na(Theta.true))) %>% 
	   group_by(Model, Phi.true, P.true,  N.marked, Theta.true) %>%
	   group_modify(~ head(.x, 100)) %>%  
	   summarize(Phi.mean=mean(Phi.estimate),
          Phi.lci =	 quantile(Phi.estimate, 0.025),
          Phi.uci =	 quantile(Phi.estimate, 0.975),
		  Phi.median = quantile(Phi.estimate, 0.5),
          Phi.bias.mean = mean(Phi.estimate) - unique(Phi.true),
          Phi.bias.median = median(Phi.estimate) - unique(Phi.true),
          Phi.bias.uci = quantile(Phi.estimate, 0.975) - unique(Phi.true),
          Phi.bias.lci = quantile(Phi.estimate, 0.025) - unique(Phi.true),
	  	  Phi.sc.MSE =  sum((Phi.estimate-unique(Phi.true))^2)/((n())-1)*unique(Phi.true)^2,
		  Phi.CI.width= quantile(Phi.estimate, 0.975)-quantile(Phi.estimate, 0.025),
		  Phi.CI.coverage= sum(Phi.CI.coverage)/n(),
          Phi.mean.n.eff=round(mean(Phi.n.eff)),
		  P.mean=mean(P.estimate), 
          P.lci = quantile(P.estimate, 0.025),
          P.uci = quantile(P.estimate, 0.975),	
		  P.median = quantile(P.estimate, 0.5),
          P.bias.mean = mean(P.estimate) - unique(P.true),
          P.bias.median = median(P.estimate) - unique(P.true),
          P.bias.uci = quantile(P.estimate, 0.975) - unique(P.true),
          P.bias.lci = quantile(P.estimate, 0.025) - unique(P.true),
	  	  P.sc.MSE =  sum((P.estimate-unique(P.true))^2)/((n())-1)*unique(P.true)^2,
		  P.CI.width= quantile(P.estimate, 0.975)-quantile(P.estimate, 0.025),
		  P.CI.coverage= sum(P.CI.coverage)/n(),
          P.mean.n.eff=round(mean(P.n.eff)),
		  Theta.mean=mean(Theta.estimate),
          Theta.lci = quantile(Theta.estimate, 0.025, na.rm=TRUE),
          Theta.uci = quantile(Theta.estimate, 0.975, na.rm=TRUE),
          Theta.median = quantile(Theta.estimate, 0.5, na.rm=TRUE),
          Theta.bias.mean = mean(Theta.estimate) - unique(Theta.true),
          Theta.bias.median = median(Theta.estimate) - unique(Theta.true),
          Theta.bias.uci = quantile(Theta.estimate, 0.975, na.rm=TRUE) - unique(Theta.true),
          Theta.bias.lci = quantile(Theta.estimate, 0.025, na.rm=TRUE) - unique(Theta.true),
	  	  Theta.sc.MSE =  sum((Theta.mean-unique(Theta.true))^2)/((n())-1)*unique(Theta.true)^2,
		  Theta.CI.width= quantile(Theta.estimate, 0.975, na.rm=TRUE)-quantile(Theta.estimate, 0.025, na.rm=TRUE),
		  Theta.CI.coverage= sum(Theta.CI.coverage)/n(),
          Theta.mean.n.eff=round(mean(Theta.n.eff)),
          N_models=n()) %>%
		  as.data.frame()

		  
Res_Phi_long<-Res %>% select(starts_with('Phi.'), Model, N_models, P.true, Theta.true)%>% add_column(Parameter='Phi')
names(Res_Phi_long)[1:13]<-substr(names(Res_Phi_long)[1:13], 5, 100)
Res_Phi_long$Phi.true<-Res_Phi_long$true

Res_P_long<-Res %>% select(starts_with('P.'), Model, N_models,Phi.true, Theta.true) %>% add_column(Parameter='P') 
names(Res_P_long)[1:13]<-substr(names(Res_P_long)[1:13], 3, 100)
Res_P_long$P.true<-Res_P_long$true
Res_P_long$Parameter<-'P'
Res_Theta_long<-Res %>% select(starts_with('Theta.'), Model, N_models, Phi.true, P.true)%>% add_column(Parameter='Theta') 
names(Res_Theta_long)[1:13]<-substr(names(Res_Theta_long)[1:13], 7, 100)
Res_Theta_long$Theta.true<-Res_Theta_long$true

Res_long<-bind_rows(Res_Phi_long, Res_P_long, Res_Theta_long) %>% mutate(Parameter=factor(Parameter, levels=c('Phi', 'P', 'Theta')))

```{r}

plot

```{r}
p.bias <- ggplot(data = Res_long, aes(x=factor(Theta.true), 
               color=as.factor(Model),  group=1)) +
		 xlab('Theta values') +
		 ylab('Phi bias') +
         geom_linerange(aes(x= factor(Theta.true), ymin=bias.lci,  ymax=bias.uci), 
		                position=position_dodge2(width=0.8)) + 
	     geom_point(aes(x= factor(Theta.true), y=bias.median,  shape = factor(Phi.true)), 
		                position=position_dodge2(width=0.8), size=1.5) +
	     facet_grid( P.true ~ Parameter ,labeller = labeller(.rows = label_both, .cols = label_both)) +
		 geom_hline( aes(yintercept = 0),  colour = grey(0.5), linetype='dashed') +
		 theme_bw()+
		 scale_shape(name= 'Phi value') +
		 scale_color_npg(name='Model',
            breaks=c('CJS_naive', 'CJSm_ccc_multinom', 'CJSm_ccc_multinom_approx'),
            labels=c('Naive CJS', 'CJSm', 'CJSm approx')) +
		 theme(legend.position = "bottom") 

p.bias

ggsave("bias_25_year.pdf", width = 20, height = 20, units = "cm", dpi=600)		  
```

### 1.5 Make Figure 2

This figure is somewhat similat to the previous one, but uses simulations where animals were marked only at the first occasion.

First we load result files and format data

```{r} 
out_all<-readRDS('./results/Sim_results_one_cohort.RDS')

# ok, now I want to combine these  
library(tidyr)
library(dplyr)

Res <- out_all %>% 
       rename(Phi.estimate = Phi.mean, P.estimate = P.mean, Theta.estimate = Theta.mean) %>%
       filter(Converged & N.marked !=100 & !(Model %in% c('CJS','CJS_no_misr') & is.na(Theta.true))) %>% 
	   group_by(Model, Phi.true, P.true,  N.marked, Theta.true) %>%
	   group_modify(~ head(.x, 100)) %>%  
	   summarize(Phi.mean=mean(Phi.estimate),
          Phi.lci =	 quantile(Phi.estimate, 0.025),
          Phi.uci =	 quantile(Phi.estimate, 0.975),
          Phi.median =	 quantile(Phi.estimate, 0.5),
          Phi.bias.mean = mean(Phi.estimate) - unique(Phi.true),
          Phi.bias.median = median(Phi.estimate) - unique(Phi.true),
          Phi.bias.uci = quantile(Phi.estimate, 0.975) - unique(Phi.true),
          Phi.bias.lci = quantile(Phi.estimate, 0.025) - unique(Phi.true),
	  	  Phi.sc.MSE =  sum((Phi.estimate-unique(Phi.true))^2)/((n())-1)*unique(Phi.true)^2,
		  Phi.CI.width= quantile(Phi.estimate, 0.975)-quantile(Phi.estimate, 0.025),
		  Phi.CI.coverage= sum(Phi.CI.coverage)/n(),
          Phi.mean.n.eff=mean(Phi.n.eff),
		  P.mean=mean(P.estimate),
          P.lci = quantile(P.estimate, 0.025),
          P.uci = quantile(P.estimate, 0.975),
          P.median = quantile(P.estimate, 0.5),
          P.bias.mean = mean(P.estimate) - unique(P.true),
          P.bias.median = median(P.estimate) - unique(P.true),
          P.bias.uci = quantile(P.estimate, 0.975) - unique(P.true),
          P.bias.lci = quantile(P.estimate, 0.025) - unique(P.true),
	  	  P.sc.MSE =  sum((P.estimate-unique(P.true))^2)/((n())-1)*unique(P.true)^2,
		  P.CI.width= quantile(P.estimate, 0.975)-quantile(P.estimate, 0.025),
		  P.CI.coverage= sum(P.CI.coverage)/n(),
          P.mean.n.eff=mean(P.n.eff),
		  Theta.mean=mean(Theta.estimate) ,
          Theta.lci = quantile(Theta.estimate, 0.025, na.rm=TRUE),
          Theta.uci = quantile(Theta.estimate, 0.975, na.rm=TRUE),
          Theta.median = quantile(Theta.estimate, 0.5, na.rm=TRUE),
          Theta.bias.mean = mean(Theta.estimate) - unique(Theta.true),
          Theta.bias.median = median(Theta.estimate) - unique(Theta.true),
          Theta.bias.uci = quantile(Theta.estimate, 0.975, na.rm=TRUE) - unique(Theta.true),
          Theta.bias.lci = quantile(Theta.estimate, 0.025, na.rm=TRUE) - unique(Theta.true),
	  	  Theta.sc.MSE =  sum((Theta.mean-unique(Theta.true))^2)/((n())-1)*unique(Theta.true)^2,
		  Theta.CI.width= quantile(Theta.estimate, 0.975, na.rm=TRUE)-quantile(Theta.estimate, 0.025, na.rm=TRUE),
		  Theta.CI.coverage= sum(Theta.CI.coverage)/n(),
          Theta.mean.n.eff=mean(Theta.n.eff),
          N_models=n()) %>%
		  as.data.frame()
```

and now plot the results

```{r}
library(ggplot2)
Ylim<-range(c(Res$Phi.bias.median, Res$P.bias.median, Res$Theta.bias.median), na.rm=TRUE)

p.Phi.bias <- ggplot(data = Res, aes(x=factor(Theta.true), 
               color=as.factor(Model),  group=1)) +
		 xlab('Theta values') +
		 ylab('Phi bias') +
		 ylim(Ylim) +
		 geom_hline( aes(yintercept = 0),  colour = grey(0.5), linetype='dashed') +
         geom_linerange(aes(x= factor(Theta.true), ymin=Phi.bias.lci,  ymax=Phi.bias.uci), position=position_dodge2(width=0.5)) + 
	     geom_point(aes(x= factor(Theta.true), y=Phi.bias.median), position=position_dodge2(width=0.5)) +
	     facet_grid( P.true ~ N.marked ,labeller = labeller(.rows = label_both, .cols = label_both)) +
		 theme_bw() +
		 scale_color_npg(name='Model',
                         breaks=c('CJS', 'RDM', 'RDMa'),
                         labels=c('CJS', 'RDM', 'RDMa')) +
		 theme(legend.position = "bottom") 
				 		 
p.Phi.bias

ggsave("Phi_10_years.pdf", width = 20, height = 10, units = "cm", dpi=600)		  

p.P.bias <- ggplot(data = Res, aes(x=factor(Theta.true), 
               color=as.factor(Model),  group=1)) +
		 xlab('Theta values') +
		 ylab('P bias') +
 		 geom_hline( aes(yintercept = 0),  colour = grey(0.5), linetype='dashed') +
         geom_linerange(aes(x= factor(Theta.true), ymin=P.bias.lci,  ymax=P.bias.uci), position=position_dodge2(width=0.5)) + 
	     geom_point(aes(x= factor(Theta.true), y=P.bias.median), position=position_dodge2(width=0.5)) +
	     facet_grid( P.true ~ N.marked ,labeller = labeller(.rows = label_both, .cols = label_both)) +
		 theme_bw() +
		 scale_color_npg(name='Model',
                         breaks=c('CJS', 'RDM', 'RDMa'),
                         labels=c('CJS', 'RDM', 'RDMa')) +
		 theme(legend.position = "bottom") 


p.P.bias
ggsave("P_10_years.pdf", width = 20, height = 10, units = "cm", dpi=600)		  

p.Theta.bias <- ggplot(data = Res, aes(x=factor(Theta.true), 
               color=as.factor(Model),  group=1)) +
		 xlab('Theta values') +
		 ylab('Theta') +
		 geom_hline( aes(yintercept = 0),  colour = grey(0.5), linetype='dashed') +
         geom_linerange(aes(x= factor(Theta.true), ymin=Theta.bias.lci,  ymax=Theta.bias.uci),
              		 position=position_dodge2(width=0.5)) + 
	     geom_point(aes(x= factor(Theta.true), y=Theta.bias.median), position=position_dodge2(width=0.5)) +
	     facet_grid( P.true ~ N.marked ,labeller = labeller(.rows = label_both, .cols = label_both)) +
		 theme_bw() +
		 scale_color_npg(name='Model',
                         breaks=c('CJS', 'RDM', 'RDMa'),
                         labels=c('CJS', 'RDM', 'RDMa')) +
		 theme(legend.position = "bottom") 

p.Theta.bias
ggsave("Theta_10_years.pdf", width = 20, height = 10, units = "cm", dpi=600)		  

```

## 2. Estimation of trend in survival probability over time
In this part we simulate data with no trend in survival over time and then run CJS, RDM and RDMa with Time trend

## 2.1 wrapper function to run the models

```{r}

run_T_c_c_models<-function(CH, n.chains=6,  n.adapt=NULL, RHat_limit=1.1, m
        ax.iter=200000, iter.increment=5000, models.to.run=c('CJS', 'RDM', 'RDMa')) {
  # this function runs Phi_dot_p_dot_Theta_dot models
  Res<-list(CH=CH)
  
  if ('CJS' %in% models.to.run) {
     cat('running CJS-T-c model')
     # flatten the data
     CH_flat<-apply(CH, c(1,2), sign)
     # Create vector with occasion of marking
     f <- apply(CH_flat, 1, get.first)
     Order<-(order(f))
     CH_flat<-CH_flat[Order,]
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
   phi.beta ~ dnorm(0,10)          # Prior for slope in survival
   mean.p ~ dunif(0, 1)           # Prior for mean recapture
   
   # estimate mean phi
   logit(mean.phi)<-phi.mu+phi.beta*(n.occasions/2)

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
    jags.data <- list(y = CH_flat, f = f, nind = dim(CH_flat)[1], 
	                  n.occasions = dim(CH_flat)[2])

    inits <- function(){list(phi.mu = rnorm(1, 0,1), 
                         phi.beta = rnorm(1, 0,0.5),
						 mean.p = runif(1, 0.5,1),
                         z = known.state.cjs(CH_flat))}

   # Parameters monitored
    parameters <- c("mean.p", "mean.phi", "phi.mu", "phi.beta")

    # Call JAGS from R 
    cjs.T.c <- autojags(jags.data, inits, 
                       parameters, "cjs-T-c.jags", n.chains = n.chains, 
                       parallel=TRUE, 
                       Rhat.limit=RHat_limit, n.adapt=NULL, 
                       max.iter=max.iter, iter.increment= iter.increment, save.all.iter=FALSE)  
    cat('   Done!\n')
    Res$CJS<-cjs.T.c
  }
  if ('RDM' %in% models.to.run | 'RDMa' %in% models.to.run) {
    # Create vector with occasion of marking
    f <- apply(CH, 1, get.first)
    Order<-(order(f))
    CH<-CH[Order,]
    f<-f[Order]
	
	sightings_and_fresh<-colSums(CH) # count per year
    All.Years=data.frame(Years=1:dim(CH)[2], N_marked=0)
    Rle<-rle(apply(CH, 1, get.first))
    for (year in 1:length(Rle$lengths)) {
      All.Years$N_marked[Rle$values[year]]<-Rle$lengths[year]
    }

    N_marked<-  All.Years$N_marked
    N_reads<-sightings_and_fresh-N_marked
    N_reads<-N_reads[-1]
    sum_N_marks_used<-cumsum(All.Years$N_marked)

    # Bundle data
    jags.data <- list(y = CH, f = f, nind = dim(CH)[1],
                     n.occasions = dim(CH)[2],
                     sum_N_marks_used=sum_N_marks_used,
                     N_reads=N_reads, N_reads2=N_reads)

    inits <- function(){list(phi.mu = rnorm(1, 0,1), 
                         phi.beta = rnorm(1, 0,0.5),
						 mean.p = runif(1, 0.5,1),
                         z = known.state.cjs.mult(CH),
						 mean.theta=runif(1, 0.5, 1))}

    # Parameters monitored
    parameters <- c("mean.phi", "phi.mu", "phi.beta" , "mean.p", "mean.theta")
  } 
  
   if ('RDM' %in% models.to.run) {
    # Specify model in BUGS language
    sink("RDM-T-c-c.jags")
    cat("
    model {
    # Priors and constraints
    for (i in 1:nind){
       for (t in f[i]:(n.occasions-1)){
          logit(phi[i,t]) <- phi.mu+phi.beta*t
          lambda_true[i,t]<-mean.lambda.true
       } #t
    } #i
    phi.mu ~ dnorm(0,1)          # Prior for intercept in survival
    phi.beta ~ dnorm(0,10)         # Prior for slope in survival
    mean.p ~ dunif(0.5, 1)           # Prior for mean recapture
	mean.lambda.true<- -log(1-mean.p)
    mean.theta ~ dunif(0.5,1)        # Prior for mean correct identification
    for (t in 1:(n.occasions-1)) {
    theta.t[t]<-mean.theta
  }

    # estimate mean phi
    logit(mean.phi)<-phi.mu+phi.beta*(n.occasions/2)

  
        # Likelihood 
        for (i in 1:nind){
          # Define latent state at first capture
          z[i,f[i]] <- 1
          for (t in (f[i]+1):n.occasions){
            # State process
            z[i,t] ~ dbern(mu1[i,t])
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]
            # Observation process
            lambda_obs[i,t-1]<- lambda_true[i,t-1]*z[i,t]*theta.t[t-1]+
		                Lambda_z_sum[t-1] -
                        lambda_true[i,t-1]*z[i,t]*(1-theta.t[t-1])/(sum_N_marks_used[t-1]-1)
          } #t
        } #i
	
        for (t in 1:(n.occasions-1)) {
	    # here I add the total sum over all lambda_z      
	    Lambda_z_sum[t]<-sum(lambda_true[1:(sum_N_marks_used[t]),t]*
	               z[1:(sum_N_marks_used[t]),t+1])*
	    		   (1-theta.t[t])/(sum_N_marks_used[t]-1)
        # Below is the same but in RDMa model
		#Lambda_z_sum[t]<-N_reads2[t]*
	   	#	  (1-theta.t[t])/(sum_N_marks_used[t]-1)		
        Lambda_obs_sum[t]<-sum(lambda_obs[1:(sum_N_marks_used[t]),t])
        P_multin[1:(sum_N_marks_used[t]),t] <- 
		      lambda_obs[1:(sum_N_marks_used[t]),t]/
			  ifelse(Lambda_obs_sum[t]==0, 1, Lambda_obs_sum[t])
        N_reads[t]~ dpois(Lambda_obs_sum[t])
        y[1:sum_N_marks_used[t],t+1] ~ 
		      dmulti(P_multin[1:sum_N_marks_used[t],t],
	          N_reads2[t]) 
        }
    }
    ",fill = TRUE)
    sink()
   
    # Call JAGS from R 
    RDM.T.c.c <- autojags(jags.data, inits, parameters,
                       "RDM-T-c-c.jags", n.chains = n.chains, 
                       parallel=TRUE, 
                       Rhat.limit=RHat_limit, n.adapt=NULL, 
                       max.iter=max.iter, iter.increment= iter.increment, save.all.iter=FALSE)
					   cat('   Done!\n')
    Res$RDM<-RDM.T.c.c
  }

  if ('RDMa' %in% models.to.run) {
    cat('running RDMa-T-c-c model')
    # Specify model in BUGS language
    sink("RDMa-T-c-c.jags")
    cat("
    model {
    # Priors and constraints
    for (i in 1:nind){
       for (t in f[i]:(n.occasions-1)){
          logit(phi[i,t]) <- phi.mu+phi.beta*t
          lambda_true[i,t]<-mean.lambda.true
       } #t
    } #i
    phi.mu ~ dnorm(0,1)          # Prior for intercept in survival
    phi.beta ~ dnorm(0,10)         # Prior for slope in survival
    mean.p ~ dunif(0.5, 1)           # Prior for mean recapture
	mean.lambda.true<- -log(1-mean.p)
    mean.theta ~ dunif(0.5,1)        # Prior for mean correct identification
    for (t in 1:(n.occasions-1)) {
    theta.t[t]<-mean.theta
    }

	# estimate mean phi
    logit(mean.phi)<-phi.mu+phi.beta*(n.occasions/2)

        # Likelihood 
        for (i in 1:nind){
          # Define latent state at first capture
          z[i,f[i]] <- 1
          for (t in (f[i]+1):n.occasions){
            # State process
            z[i,t] ~ dbern(mu1[i,t])
            mu1[i,t] <- phi[i,t-1] * z[i,t-1]
            # Observation process
            lambda_obs[i,t-1]<- lambda_true[i,t-1]*z[i,t]*theta.t[t-1]+
		                  Lambda_z_sum[t-1] -
  						  lambda_true[i,t-1]*z[i,t]*(1-theta.t[t-1])/(sum_N_marks_used[t-1]-1)
          } #t
        } #i
	
        for (t in 1:(n.occasions-1)) {
	    # here I add the total sum over all lambda_z      
	    #Lambda_z_sum[t]<-sum(lambda_true[1:(sum_N_marks_used[t]),t]*
	    #           z[1:(sum_N_marks_used[t]),t+1])*
	    #		   (1-theta.t[t])/(sum_N_marks_used[t]-1)
        # above is the RDM model
		Lambda_z_sum[t]<-N_reads2[t]*
	   		  (1-theta.t[t])/(sum_N_marks_used[t]-1)		
        Lambda_obs_sum[t]<-sum(lambda_obs[1:(sum_N_marks_used[t]),t])
        P_multin[1:(sum_N_marks_used[t]),t] <- 
		      lambda_obs[1:(sum_N_marks_used[t]),t]/
			  ifelse(Lambda_obs_sum[t]==0, 1, Lambda_obs_sum[t])
        N_reads[t]~ dpois(Lambda_obs_sum[t])
        y[1:sum_N_marks_used[t],t+1] ~ 
		      dmulti(P_multin[1:sum_N_marks_used[t],t],
	          N_reads2[t]) 
        }
    }
    ",fill = TRUE)
    sink()

    # Call JAGS from R 
    RDMa.T.c.c <- autojags(jags.data, inits, parameters,
                       "RDMa-T-c-c.jags", n.chains = n.chains, 
                       parallel=TRUE, 
                       Rhat.limit=RHat_limit, n.adapt=NULL, 
                       max.iter=max.iter, iter.increment= iter.increment, save.all.iter=FALSE)
					   cat('   Done!\n')
    Res$RDMa<-RDMa.T.c.c
   }
   return(Res)
}
```
### 2.1 define simulation parameters
#### Define parameter values
```{r}
n.occasions <- 25                   # Number of capture occasions
marked <- rep(25, n.occasions-1)   # Annual number of newly marked individuals
phi <- rep(0.7, n.occasions-1) # survival
p <- rep(0.9, n.occasions-1)  # resighting
theta<- rep(0.95, n.occasions-1) # correctly identified
```

#### Define matrices with survival recapture and correct identification
```{r}
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))
THETA<-matrix(theta, ncol = n.occasions-1, nrow = sum(marked))

RHat_limit=1.01
max.iter=8000
n.adapt=4000
iter.increment=4000
nRuns=5
```
### 2.2 run models
```{r}
if (Run_everything) {

Output_df<-data.frame(Model=character(), 
            Phi.true=numeric(), P.true=numeric(), N.marked = integer(), Theta.true=numeric(),
            Phi.mu=numeric(), Phi.beta=c(), P.mean =numeric(), Theta.mean=numeric(), 
			Deviance =numeric(), 
			Phi.CI.width =numeric(),P.CI.width=numeric(), Theta.CI.width=numeric(),
			Phi.CI.coverage = logical(), P.CI.coverage = logical(), Theta.CI.coverage = logical(),
			Phi.n.eff=integer(), P.n.eff=integer(), Theta.n.eff=integer(), Deviance.n.eff=integer(),
			Phi.RHat =numeric(), P.RHat =numeric(), Theta.RHat =numeric(), Deviance.RHat =numeric(),
			n.chains =integer(), n.burnin=integer(),  n.iter=integer(), n.elapsed.mins  = integer(), 
			Converged = logical(), Remark = character(),
			Phi.bias = numeric(),P.bias = numeric(), Theta.bias = numeric())

						   
for (i in 1:nRuns) {
   cat('i = ', i, '\n')
   # simulate data
   CH_list <- simul.cjs.multiple.sightings(PHI, P, THETA,  marked)
   # estimate models
   Res_list<-run_T_c_c_models(CH_list$CH, n.chains=6,  
                              n.adapt=n.adapt, RHat_limit=RHat_limit,
							  max.iter=max.iter, iter.increment=iter.increment,
							models.to.run=c('CJS', 'RDM', 'RDMa'))

   # extract results					
    if (!is.null(Res_list$CJS)) {
	    CI_cur.CJS<-cbind(c('mean.phi' = phi[1], 'phi.mu'=0, 
		            'phi.beta'=0, 'mean.p'=p[1], 
                    'mean.theta'=theta[1]), 
                  c(Res_list$CJS$q2.5$mean.phi,
				    Res_list$CJS$q2.5$phi.mu,
				    Res_list$CJS$q2.5$phi.beta, 
				    Res_list$CJS$q2.5$mean.p,
					NA), 
				  c(Res_list$CJS$q97.5$mean.phi,
				    Res_list$CJS$q97.5$phi.mu,
				    Res_list$CJS$q97.5$phi.beta, 
				    Res_list$CJS$q97.5$mean.p,
					NA)) 
    	Output_line_CJS<-data.frame(Model='CJS', 
            Phi.true=phi[1], P.true=p[1], N.marked = marked[1], Theta.true=theta[1],
            Phi.mu=Res_list$CJS$mean$phi.mu, Phi.beta=Res_list$CJS$mean$phi.mu,
			Phi.mean =Res_list$CJS$mean$mean.phi, P.mean =Res_list$CJS$mean$mean.p, Theta.mean=NA, 
			Deviance =Res_list$CJS$mean$deviance, 
			Phi.CI.width =CI_cur.CJS[1,3]-CI_cur.CJS[1,2],
			P.CI.width=CI_cur.CJS[2,3]-CI_cur.CJS[2,2], Theta.CI.width=NA,
			Phi.CI.coverage =  CI_cur.CJS[1,1]>CI_cur.CJS[1,2] & CI_cur.CJS[1,1]<=CI_cur.CJS[1,3], 
			P.CI.coverage =    CI_cur.CJS[2,1]>CI_cur.CJS[2,2] & CI_cur.CJS[2,1]<=CI_cur.CJS[2,3],
			Theta.CI.coverage = NA,
			Phi.n.eff=unlist(Res_list$CJS$n.eff)[1], 
			P.n.eff=unlist(Res_list$CJS$n.eff)[2], 
			Theta.n.eff=NA, 
			Deviance.n.eff=unlist(Res_list$CJS$n.eff)[3],
			Phi.RHat = unlist(Res_list$CJS$Rhat)[1], P.RHat =unlist(Res_list$CJS$Rhat)[2],
			Theta.RHat =NA, Deviance.RHat =unlist(Res_list$CJS$Rhat)[3],
			n.chains =Res_list$CJS$mcmc.info$n.chains, 
			n.burnin=Res_list$CJS$mcmc.info$n.burnin,
			n.iter=Res_list$CJS$mcmc.info$n.iter, 
			n.elapsed.mins  = Res_list$CJS$mcmc.info$elapsed.mins, 
			Converged = max(unlist(Res_list$CJS$Rhat), na.rm=TRUE)<1.01,
			Remark = '',
			Phi.bias = Res_list$CJS$mean$mean.phi-phi[1],
			P.bias = Res_list$CJS$mean$mean.p-p[1],
			Theta.bias = NA)
	Output_df<-rbind(Output_df, Output_line_CJS)		
    }
	
    if (!is.null(Res_list$RDM)) {	
	    CI_cur.RDM<-cbind(c('mean.phi' = phi[1], 'phi.mu'=0, 
		            'phi.beta'=0, 'mean.p'=p[1], 
                    'mean.theta'=theta[1]), 
                  c(Res_list$RDM$q2.5$mean.phi,
				    Res_list$RDM$q2.5$phi.mu,
				    Res_list$RDM$q2.5$phi.beta, 
				    Res_list$RDM$q2.5$mean.p,
					Res_list$RDM$q2.5$mean.theta), 
				  c(Res_list$RDM$q97.5$mean.phi,
				    Res_list$RDM$q97.5$phi.mu,
				    Res_list$RDM$q97.5$phi.beta, 
				    Res_list$RDM$q97.5$mean.p,
					Res_list$RDM$q97.5$mean.theta)) 
    	Output_line_RDM<-data.frame(Model='RDM', 
            Phi.true=phi[1], P.true=p[1], N.marked = marked[1], Theta.true=theta[1],
            Phi.mu=Res_list$RDM$mean$phi.mu, Phi.beta=Res_list$RDM$mean$phi.mu,
			Phi.mean =	Res_list$RDM$mean$mean.phi, P.mean =Res_list$RDM$mean$mean.p
			           , Theta.mean=Res_list$RDM$mean$mean.theta, 
			Deviance =Res_list$RDM$mean$deviance, 
			Phi.CI.width =CI_cur.RDM[1,3]-CI_cur.RDM[1,2],
			P.CI.width=CI_cur.RDM[2,3]-CI_cur.RDM[2,2], Theta.CI.width=CI_cur.RDM[3,3]-CI_cur.RDM[3,2],
			Phi.CI.coverage =  CI_cur.RDM[1,1]>CI_cur.RDM[1,2] & CI_cur.RDM[1,1]<=CI_cur.RDM[1,3], 
			P.CI.coverage =    CI_cur.RDM[2,1]>CI_cur.RDM[2,2] & CI_cur.RDM[2,1]<=CI_cur.RDM[2,3],
			Theta.CI.coverage =  CI_cur.RDM[3,1]>CI_cur.RDM[3,2] & CI_cur.RDM[3,1]<=CI_cur.RDM[3,3],
			Phi.n.eff=unlist(Res_list$RDM$n.eff)[1], 
			P.n.eff=unlist(Res_list$RDM$n.eff)[2], 
			Theta.n.eff=unlist(Res_list$RDM$n.eff)[3], 
			Deviance.n.eff=unlist(Res_list$RDM$n.eff)[3],
			Phi.RHat = unlist(Res_list$RDM$Rhat)[1], P.RHat =unlist(Res_list$RDM$Rhat)[2],
			Theta.RHat =unlist(Res_list$RDM$Rhat)[3], Deviance.RHat =unlist(Res_list$RDM$Rhat)[3],
			n.chains =Res_list$RDM$mcmc.info$n.chains, 
			n.burnin=Res_list$RDM$mcmc.info$n.burnin,
			n.iter=Res_list$RDM$mcmc.info$n.iter, 
			n.elapsed.mins  = Res_list$RDM$mcmc.info$elapsed.mins, 
			Converged = max(unlist(Res_list$RDM$Rhat), na.rm=TRUE)<1.01,
			Remark = '',
			Phi.bias = Res_list$RDM$mean$mean.phi-phi[1],
			P.bias = Res_list$RDM$mean$mean.p-p[1],
			Theta.bias = Res_list$RDM$mean$mean.theta-theta[1])
		
		Output_df<-rbind(Output_df, Output_line_RDM)		
		}
		
    if (!is.null(Res_list$RDMa)) {
	    CI_cur.RDMa<-cbind(c('mean.phi' = phi[1], 'phi.mu'=0, 
		            'phi.beta'=0, 'mean.p'=p[1], 
                    'mean.theta'=theta[1]), 
                  c(Res_list$RDMa$q2.5$mean.phi,
				    Res_list$RDMa$q2.5$phi.mu,
				    Res_list$RDMa$q2.5$phi.beta, 
				    Res_list$RDMa$q2.5$mean.p,
					Res_list$RDMa$q2.5$mean.theta), 
				  c(Res_list$RDMa$q97.5$mean.phi,
				    Res_list$RDMa$q97.5$phi.mu,
				    Res_list$RDMa$q97.5$phi.beta, 
				    Res_list$RDMa$q97.5$mean.p,
					Res_list$RDMa$q97.5$mean.theta)) 
    	Output_line_RDMa<-data.frame(Model='RDMa', 
            Phi.true=phi[1], P.true=p[1], N.marked = marked[1], Theta.true=theta[1],
            Phi.mu=Res_list$RDMa$mean$phi.mu, Phi.beta=Res_list$RDMa$mean$phi.mu,
			Phi.mean =Res_list$RDMa$mean$mean.phi, P.mean =Res_list$RDMa$mean$mean.p,
			Theta.mean=Res_list$RDMa$mean$mean.theta, 
			Deviance =Res_list$RDMa$mean$deviance, 
			Phi.CI.width =CI_cur.RDMa[1,3]-CI_cur.RDMa[1,2],
			P.CI.width=CI_cur.RDMa[2,3]-CI_cur.RDMa[2,2], Theta.CI.width=CI_cur.RDMa[3,3]-CI_cur.RDMa[3,2],
			Phi.CI.coverage =  CI_cur.RDMa[1,1]>CI_cur.RDMa[1,2] & CI_cur.RDMa[1,1]<=CI_cur.RDMa[1,3], 
			P.CI.coverage =    CI_cur.RDMa[2,1]>CI_cur.RDMa[2,2] & CI_cur.RDMa[2,1]<=CI_cur.RDMa[2,3],
			Theta.CI.coverage =  CI_cur.RDMa[3,1]>CI_cur.RDMa[3,2] & CI_cur.RDMa[3,1]<=CI_cur.RDMa[3,3],
			Phi.n.eff=unlist(Res_list$RDMa$n.eff)[1], 
			P.n.eff=unlist(Res_list$RDMa$n.eff)[2], 
			Theta.n.eff=unlist(Res_list$RDMa$n.eff)[3], 
			Deviance.n.eff=unlist(Res_list$RDMa$n.eff)[3],
			Phi.RHat = unlist(Res_list$RDMa$Rhat)[1], P.RHat =unlist(Res_list$RDMa$Rhat)[2],
			Theta.RHat =unlist(Res_list$RDMa$Rhat)[3], Deviance.RHat =unlist(Res_list$RDMa$Rhat)[3],
			n.chains =Res_list$RDMa$mcmc.info$n.chains, 
			n.burnin=Res_list$RDMa$mcmc.info$n.burnin,
			n.iter=Res_list$RDMa$mcmc.info$n.iter, 
			n.elapsed.mins  = Res_list$RDMa$mcmc.info$elapsed.mins, 
			Converged = max(unlist(Res_list$RDMa$Rhat), na.rm=TRUE)<1.01,
			Remark = '',
			Phi.bias = Res_list$RDMa$mean$mean.phi-phi[1],
			P.bias = Res_list$RDMa$mean$mean.p-p[1],
			Theta.bias = Res_list$RDMa$mean$mean.theta-theta[1])
	
    	Output_df<-rbind(Output_df, Output_line_RDMa)		
	}	

   ###
   # once in e.g. 10 runs we will save the result..
   if (i%%10 ==0) { 
     save(Output_df,
		file=paste0('Output_df_slope_Phi_dot_', Phi_dot, '_P_dot_', P_dot, 
            '_theta_dot_', Theta_dot, '_n_occasions_', n.occasions, 
	        '_marked_', n.marked, '_iteration_', i, '_RHat_', RHat_limit, '.Rdata'))
    }
   }
			 
   Filename<-paste0('Output_df_slope_fin_', nRuns, '_simulations_Phi_dot_', 
      Phi_dot, '_p_dot_', P_dot, 
      '_Theta_dot_', Theta_dot,
	  '_n_occasions_',n.occasions, '_marked_', n.marked, '_RHat_', RHat_limit, '.RData')	
	  
   save(Output_df, file=Filename)
}
```



### 2.6 Plot the results

#### load model outputs
We have repeated the runds above 100 times, and combined the results into a data.frame.
  
saveRDS(Res_sim_trend, file='Res_sim_trend.RDS')

```{r}
# look at the bias in the results from CJS model.

summary(Res_sim_trend$Phi.beta[Res_sim_trend$Model=='CJS'])


Res_sim_trend<-readRDS('Res_sim_trend.RDS') 
  
XX<-seq(1, 25, by=1)  
Res_sim_trend_CJS_no_misr<-Res_sim_trend %>% filter(Model=='CJS_no_misr')

CJS_no_mr_Median<-sapply(XX, FUN=function(x) 
       plogis(quantile(Res_sim_trend_CJS_no_misr$Phi.mu +
       x*Res_sim_trend_CJS_no_misr$Phi.beta, 0.5)))
	   
CJS_no_mr_LCI<-sapply(XX, FUN=function(x) 
       plogis(quantile(Res_sim_trend_CJS_no_misr$Phi.mu +
       x*Res_sim_trend_CJS_no_misr$Phi.beta, 0.025)))
CJS_no_mr_UCI<-sapply(XX, FUN=function(x)
       plogis(quantile(Res_sim_trend_CJS_no_misr$Phi.mu +
	   x*Res_sim_trend_CJS_no_misr$Phi.beta, 0.975)))
	   
Res_sim_trend_CJS<-Res_sim_trend %>% filter(Model=='CJS')

CJS_Median<-sapply(XX, FUN=function(x) 
       plogis(quantile(Res_sim_trend_CJS$Phi.mu +
       x*Res_sim_trend_CJS$Phi.beta, 0.5)))
CJS_LCI<-sapply(XX, FUN=function(x) 
       plogis(quantile(Res_sim_trend_CJS$Phi.mu +
       x*Res_sim_trend_CJS$Phi.beta, 0.025)))
CJS_UCI<-sapply(XX, FUN=function(x)
       plogis(quantile(Res_sim_trend_CJS$Phi.mu +
	   x*Res_sim_trend_CJS$Phi.beta, 0.975)))

Res_sim_trend_RDM<-Res_sim_trend %>% filter(Model=='RDM')
	   
RDM_Median<-sapply(XX, FUN=function(x) 
       plogis(quantile(Res_sim_trend_RDM$Phi.mu +
       x*Res_sim_trend_RDM$Phi.beta, 0.5)))
RDM_LCI<-sapply(XX, FUN=function(x) 
       plogis(quantile(Res_sim_trend_RDM$Phi.mu +
       x*Res_sim_trend_RDM$Phi.beta, 0.025)))
RDM_UCI<-sapply(XX, FUN=function(x)
       plogis(quantile(Res_sim_trend_RDM$Phi.mu +
	   x*Res_sim_trend_RDM$Phi.beta, 0.975)))

Res_sim_trend_RDMa<-Res_sim_trend %>% filter(Model=='RDMa')
	   
RDMa_Median<-sapply(XX, FUN=function(x) 
       plogis(quantile(Res_sim_trend_RDMa$Phi.mu +
       x*Res_sim_trend_RDMa$Phi.beta, 0.5)))
RDMa_LCI<-sapply(XX, FUN=function(x) 
       plogis(quantile(Res_sim_trend_RDMa$Phi.mu +
       x*Res_sim_trend_RDMa$Phi.beta, 0.025)))
RDMa_UCI<-sapply(XX, FUN=function(x)
       plogis(quantile(Res_sim_trend_RDMa$Phi.mu +
	   x*Res_sim_trend_RDMa$Phi.beta, 0.975)))
```
#### Plot
```{r}

library(dplyr)
library(ggsci)

CJS_no_mr_T<-bind_cols(XX=XX, Median=CJS_no_mr_Median, UCI=CJS_no_mr_UCI, LCI=CJS_no_mr_LCI, Model='CJS_no_mr')
CJS_T<-bind_cols(XX=XX, Median=CJS_Median, UCI=CJS_UCI, LCI=CJS_LCI, Model='CJS')
RDM_T<-bind_cols(XX=XX, Median=RDM_Median, UCI=RDM_UCI, LCI=RDM_LCI, Model='RDM')
RDMa_T<-bind_cols(XX=XX, Median=RDMa_Median, UCI=RDMa_UCI, LCI=RDMa_LCI, Model='RDMa')

all_models_T<-bind_rows(CJS_no_mr_T, CJS_T, RDM_T, RDMa_T)


ggplot(data = all_models_T, aes(x=XX, 
               color=as.factor(Model),  group=Model)) +
		 ylab('apparent survial probability') +
		 xlab('Time') +
         geom_linerange(aes(x= XX, ymin=LCI,  ymax=UCI),
	     	 position=position_dodge2(width=0.8)) + 
	     geom_line(aes(y=Median), position=position_dodge2(width=0.8)) 	+	
		 geom_hline(aes(yintercept=0.7), colour = grey(0.5), linetype='dashed') +
		 #scale_color_npg(name='Model',
         #   breaks=c('CJS', 'RDM', 'RDMa', 'CJS_no_mr'),
         #   labels=c('Naive CJS', 'RDM', 'RDMa', 'CJS_naive no misidentification')) +
		 scale_color_manual(name='Model',
            breaks=c('CJS', 'CJS_no_mr','RDM', 'RDMa'),
            labels=c('CJS', 'CJS over data without misidentification' , 'RDM', 'RDMa' ),
			values = c("#E64B35", "black", "#4DBBD5", "#00A087")) +
		 theme_bw() +
		 theme(legend.position = "bottom") 


ggsave("Slopes_figure_22_06_2020.pdf", width = 10, height = 10, units = "cm", dpi=600)		  
		 
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
                   check.names=FALSE, header=TRUE)
   x.P<-read.csv('./data/black-tailed_x.P.csv',
                   check.names=FALSE, header=TRUE)
} else {
   CH<-as.matrix(read.csv('https://git.io/fhhtJ', check.names=FALSE))
   x.Phi<-as.matrix(read.csv('https://git.io/fje0H', check.names=FALSE, header=TRUE))
   x.P<-as.matrix(read.csv('https://git.io/fje0Q', check.names=FALSE, header=TRUE))
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

  sightings_and_fresh<-colSums(CH) # count per year
  All.Years=data.frame(Years=1:dim(CH)[2], N_marked=0)
  Rle<-rle(apply(CH, 1, get.first))
  for (year in 1:length(Rle$lengths)) {
  All.Years$N_marked[Rle$values[year]]<-Rle$lengths[year]
  }

   N_marked<-  All.Years$N_marked
   N_reads<-sightings_and_fresh-N_marked
   N_reads<-N_reads[-1]
   sum_N_marks_used<-cumsum( All.Years$N_marked)

print(table(markedas))
print(table(markedas, f))

# third model
sink("Theta_t.jags")
cat("
model {
# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- phi.t[t, x.Phi[i,t]]
      #logit(p_r[i,t]) <- p.age[x.P[i,t]] + p.y[t]  
	  #lambda_true[i,t]<- -log(1-p_r[i,t]) + lambda.epsilon[i]
      log(lambda_true[i,t])<- lambda.age[x.P[i,t]] + lambda.y[t]   + lambda.epsilon[i]

	  theta[i,t] <- theta.t[t]
      } #t
   } #i

lambda.age[1]~ dnorm(0, 0.001)
lambda.age[2]~ dnorm(0, 0.001)
lambda.age[3]~ dnorm(0, 0.001)
lambda.age[4]<-0
   
for (t in 1:(n.occasions-1)){
	  lambda.y[t] ~ dnorm(0, 0.001)
}

for (i in 1:nind){
      lambda.epsilon[i] ~ dnorm(0, lambda.tau)
}

lambda.sigma ~ dunif(0, 1) # Prior for standard deviation
lambda.tau <- pow(lambda.sigma, -2)
lambda.sigma2 <- pow(lambda.sigma, 2)

####
# trend and RE for adults
# I need on intercept

#  survival
for (g in 1:2){
    phi.mu[g] ~ dnorm(0,1)          # Prior for intercept in survival
    phi.beta[g] ~ dnorm(0,10)         # Prior for slope in survival
	sigma[g] ~ dunif(0, 2) # Prior for standard deviation
    tau[g] <- pow(sigma[g], -2)
    sigma2[g] <- pow(sigma[g], 2)
}
	  
for (t in 1:(n.occasions-1)){
      for (g in 1:2){
      phi.epsilon[g,t] ~ dnorm(phi.mu[g]+phi.beta[g]*t, tau[g])
	  logit(phi.t[t,g])<- phi.epsilon[g,t]
	  }
}


for (g in 1:2){
logit(phi.average[g])<-phi.mu[g]+phi.beta[g]*(n.occasions/2)
}
      
mean.theta ~ dnorm(3, 0.01)        # Prior for mean correct identification
theta.sigma ~ dunif(0, 2) # Prior for standard deviation
theta.tau <- pow(theta.sigma, -2)
theta.sigma2 <- pow(theta.sigma, 2)


for (t in 1:(n.occasions-1)) {
     cur_theta[t]~ dnorm(mean.theta, theta.tau)
     logit(theta.t[t]) <- cur_theta[t] 
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
          lambda_obs[i,t-1]<- lambda_true[i,t-1]*z[i,t]*theta.t[t-1]+
		                  Lambda_z_sum[t-1] -
 						  lambda_true[i,t-1]*z[i,t]*(1-theta.t[t-1])/(sum_N_marks_used[t-1]-1)
        } #t
    } #i
	
	for (t in 1:(n.occasions-1)) {		
       Lambda_z_sum[t]<-N_reads2[t]*
	   		   (1-theta.t[t])/(sum_N_marks_used[t]-1)		

	   # here I add the total sum over all lambda_z      
	    #Lambda_z_sum[t]<-sum(lambda_true[1:(sum_N_marks_used[t]),t]*
	    #          z[1:(sum_N_marks_used[t]),t+1])*
	    #		   (1-theta.t[t])/(sum_N_marks_used[t]-1)
        # above is the RDM model
		
		
		
       Lambda_obs_sum[t]<-sum(lambda_obs[1:(sum_N_marks_used[t]),t])
	   
       P_multin[1:(sum_N_marks_used[t]),t] <- lambda_obs[1:(sum_N_marks_used[t]),t]/ifelse(Lambda_obs_sum[t]==0, 1, Lambda_obs_sum[t])
	   
	   N_reads[t]~ dpois(Lambda_obs_sum[t])
	   
       y[1:sum_N_marks_used[t],t+1] ~ dmulti(P_multin[1:sum_N_marks_used[t],t], N_reads[t]) 
  }
}
",fill = TRUE)
sink()

#### Bundle data, generate intitials and prepare model run
```{r}
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2],
                  N_reads=N_reads, N_reads2=N_reads, sum_N_marks_used=sum_N_marks_used, x.Phi=x.Phi,x.P=x.P)

# Initial values
inits <- function(){list(z = known.state.cjs.mult(CH),
          phi.mu = c(rnorm(1, 2,1), rnorm(1, 1, 1)),
          phi.beta = rnorm(2, 0,0.5),
		  lambda.age =c(rnorm(3, -1, 1), NA),
		  lambda.y= c(rnorm((dim(CH)[2]-2), 2, 1), NA),
		  mean.theta=rnorm(1, 3, 1),
		  sigma =runif(2, 0, 2),
		  lambda.sigma = runif(1, 0, 1),
		  theta.sigma = runif(1, 0, 1))
	  }

# Parameters monitored
parameters <- c("phi.t", "lambda.age", "mean.theta", "sigma2", "lambda.sigma2", "lambda.y", 'phi.beta', 'phi.mu', 'phi.average', "theta.sigma2", "theta.t")

# MCMC settings
max.iter<-50#4000
iter.increment <- 50
RHat_limit=1.01
n.adapt=10 #500
nt <- 1
nc <- 6
```
#### Run the model
```{r}
if (Run_everything) {
   # (brt - days)
   rdma_theta_t  <-autojags(jags.data, inits, parameters,
                       "Theta_t.jags", n.chains = nc, 
					   n.thin = nt, parallel=TRUE, n.cores=6, 
                       Rhat.limit=RHat_limit, n.adapt=n.adapt, 
                       max.iter=max.iter, iter.increment= iter.increment)

 rdma_theta_t<-update( rdma_theta_t, n.iter=500)
 # the real run rquires way longer runs.
 
 rdma.godwits<-rdma_theta_t_WI_10_25
 # here I do more updates to reach convergence
 
 # the final model outcome exceeds the 100 Mb, allowed for free storage on GitHub, so I will not cut it in pieces are glue together after reading.
 
 rdma.godwits.small<-rdma.godwits
 rdma.godwits.chunk.1<-rdma.godwits$sims.list
 rdma.godwits.chunk.2<-rdma.godwits$samples

 rdma.godwits.small$sims.list<-NA
 rdma.godwits.small$samples<-NA

 saveRDS(rdma.godwits.small, file='./results/rdma.godwits.small.RDS')
 saveRDS(rdma.godwits.chunk.1, file='./results/rdma.godwits.chunk.1.RDS')
 saveRDS(rdma.godwits.chunk.2, file='./results/rdma.godwits.chunk.2.RDS')
 }
 ############
 # now I read the files back

 rdma.godwits.small<-readRDS('./results/rdma.godwits.small.RDS')
 rdma.godwits.chunk.1<-readRDS('./results/rdma.godwits.chunk.1.RDS')
 rdma.godwits.chunk.2<-readRDS('./results/rdma.godwits.chunk.2.RDS')


 rdma.godwits.small$sims.list<-rdma.godwits.chunk.1
 rdma.godwits.small$samples<-rdma.godwits.chunk.2
 rdma.godwits<-rdma.godwits.small
print(rdma.godwits)
```

### 3.2 Naive CJS model with black-tailed godwit data
#### Prepare the data
The only difference from RDM model here is that we flatten the data replacing all positive values with 1.
```{r}
CH_flat<-apply(CH, c(1,2), sign)


sink("cjs_with_slope.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- phi.t[t, x.Phi[i,t]]
      logit(p_r[i,t]) <- p.age[x.P[i,t]] + p.y[t] +  p.epsilon[i] 
      } #t
   } #i


p.age[1]~ dnorm(0, 0.001)
p.age[2]~ dnorm(0, 0.001)
p.age[3]~ dnorm(0, 0.001)
p.age[4]<-0
   
for (t in 1:(n.occasions-1)){
	  p.y[t] ~ dnorm(0, 0.001)
}

for (i in 1:nind){
      p.epsilon[i] ~ dnorm(0, p.tau)
}

p.sigma ~ dunif(0, 10) # Prior for standard deviation
p.tau <- pow(p.sigma, -2)
p.sigma2 <- pow(p.sigma, 2)


####
# trend and RE for adults
# I need on intercept

#  survival
for (g in 1:2){
    phi.mu[g] ~ dnorm(0,1)          # Prior for intercept in survival
    phi.beta[g] ~ dnorm(0,10)         # Prior for slope in survival
	sigma[g] ~ dunif(0, 3) # Prior for standard deviation
    tau[g] <- pow(sigma[g], -2)
    sigma2[g] <- pow(sigma[g], 2)
}
	  
for (t in 1:(n.occasions-1)){
      for (g in 1:2){
      phi.epsilon[g,t] ~ dnorm(phi.mu[g]+phi.beta[g]*t, tau[g])
	  logit(phi.t[t,g])<- phi.epsilon[g,t]
	  }
}

##########
# and I want to output average survival
##########
# this simply means predisction for the average time.. 
for (g in 1:2){
   logit(phi.average[g])<-phi.mu[g]+phi.beta[g]*(n.occasions/2)
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
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p_r[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()


# Bundle data
jags.data <- list(y = CH_flat, f = f, nind = dim(CH_flat)[1],
                  n.occasions = dim(CH_flat)[2], x.Phi=x.Phi,x.P=x.P)


# Initial values
inits <- function(){list(z = known.state.cjs(CH_flat),
		  p.age =c(rnorm(3, 0, 5), NA),
		  p.y= c(rnorm((dim(CH_flat)[2]-2), 0, 5), NA),
          phi.mu = rnorm(2, 0,1), 
          phi.beta = rnorm(2, 0,0.5),
		  sigma = runif(2, 0, 3),
		  p.sigma = runif(1, 0, 5))}

# Parameters monitored
parameters <- c("phi.t", "p.age", "sigma2", "p.sigma2", "p.y", 'phi.mu', 'phi.beta', 'phi.average')

# MCMC settings
max.iter<-5000
iter.increment <- 5000
RHat_limit=1.01
n.adapt=5000 #500
nt <- 1
nc <- 6
```

#### Run the model
```{r}
if (Run_everything) {

a<-Sys.time()
   # (brt - days)
   CJS.godwits_T_t  <-autojags(jags.data, inits, parameters,
                       "cjs_with_slope.jags", n.chains = nc, 
					   n.thin = nt, parallel=TRUE, n.cores=6, 
                       Rhat.limit=RHat_limit, n.adapt=n.adapt, 
                       max.iter=max.iter, iter.increment= iter.increment)
Sys.time()-a


cjs.godwits<-update(CJS.godwits_T_t, n.iter=20000)

# again split and save approach for GitHub..

cjs.godwits.small<-cjs.godwits
cjs.godwits.chunk.1<-cjs.godwits$sims.list[1:4]
cjs.godwits.chunk.2<-cjs.godwits$sims.list[5:9]
cjs.godwits.small$sims.list<-NA

cjs.godwits.chunk.3<-cjs.godwits$samples[1:3]
cjs.godwits.chunk.4<-cjs.godwits$samples[4:6]
cjs.godwits.small$samples<-NA

saveRDS(cjs.godwits.small, file='./results/cjs.godwits.small.RDS')
saveRDS(cjs.godwits.chunk.1, file='./results/cjs.godwits.chunk.1.RDS')
saveRDS(cjs.godwits.chunk.2, file='./results/cjs.godwits.chunk.2.RDS')
saveRDS(cjs.godwits.chunk.3, file='./results/cjs.godwits.chunk.3.RDS')
saveRDS(cjs.godwits.chunk.4, file='./results/cjs.godwits.chunk.4.RDS')
}

cjs.godwits.small<-readRDS('./results/cjs.godwits.small.RDS')
cjs.godwits.chunk.1<-readRDS('./results/cjs.godwits.chunk.1.RDS')
cjs.godwits.chunk.2<-readRDS('./results/cjs.godwits.chunk.2.RDS')
cjs.godwits.chunk.3<-readRDS('./results/cjs.godwits.chunk.3.RDS')
cjs.godwits.chunk.4<-readRDS('./results/cjs.godwits.chunk.4.RDS')


cjs.godwits.small$sims.list<-c(cjs.godwits.chunk.1, cjs.godwits.chunk.2)
cjs.godwits.small$samples<-c(cjs.godwits.chunk.3, cjs.godwits.chunk.4)
cjs.godwits<-cjs.godwits.small

print(cjs.godwits)
```

### 3.3 Plot godwits survival estimates

#### Extract output
```{r}
Res_adults<-data.frame(
   Year.Phi=as.numeric(dimnames(CH)[[2]])[-length(dimnames(CH)[[2]])],
   Phi.rdma.median= apply(rdma.godwits$sims.list$phi.t[,,2],
                          2,quantile, probs=0.5),
   Phi.rdma.lci= apply(rdma.godwits$sims.list$phi.t[,,2],
                          2,quantile, probs=0.025),
   Phi.rdma.uci= apply(rdma.godwits$sims.list$phi.t[,,2],
                          2,quantile, probs=0.975),
   Phi.cjs.median= apply(cjs.godwits$sims.list$phi.t[,,2],
                          2,quantile, probs=0.5),
   Phi.cjs.lci= apply(cjs.godwits$sims.list$phi.t[,,2],
                          2,quantile, probs=0.025),
   Phi.cjs.uci= apply(cjs.godwits$sims.list$phi.t[,,2],
                          2,quantile, probs=0.975))

Res_juveniles<-data.frame(
   Year.Phi=as.numeric(dimnames(CH)[[2]])[-length(dimnames(CH)[[2]])],
   Phi.rdma.median= apply(rdma.godwits$sims.list$phi.t[,,1],
                          2,quantile, probs=0.5),
   Phi.rdma.lci= apply(rdma.godwits$sims.list$phi.t[,,1],
                          2,quantile, probs=0.025),
   Phi.rdma.uci= apply(rdma.godwits$sims.list$phi.t[,,1], 
                          2,quantile, probs=0.975),
   Phi.cjs.median= apply(cjs.godwits$sims.list$phi.t[,,1],
                          2,quantile, probs=0.5),
   Phi.cjs.lci= apply(cjs.godwits$sims.list$phi.t[,,1],
                          2,quantile, probs=0.025),
   Phi.cjs.uci= apply(cjs.godwits$sims.list$phi.t[,,1],
                          2,quantile, probs=0.975))
```
#### Plot

```{r}
pdf('godwit_figure_2.pdf', 10,5, useDingbats=FALSE)
Shift=0.1  
par(mfrow=c(1,2), mar=c(5.1, 4.1, 0.5, 0.5))
  
#Panel A juveniles

plot(Res_juveniles$Phi.rdma.median~Res_juveniles$Year.Phi,
     ylim=c(0, 1), las=1, pch=19, ylab='apparent survival',
	 xlab='Year', type='n')  

abline(v=2004:2020, lty=3, col=grey(0.5))
abline(v=seq(2005, 2020, by=5), lty=2, col=grey(0.5))

abline(h=seq(0, 1, by=0.1), lty=3, col=grey(0.5))
abline(h=seq(0, 1, by=0.5), lty=2, col=grey(0.5))

points(Res_juveniles$Phi.rdma.median~I(Res_juveniles$Year.Phi+Shift),
         las=1, pch=23, bg='#00A087', col='#00A087')   
segments(y0=Res_juveniles$Phi.rdma.lci,
         y1=Res_juveniles$Phi.rdma.uci,
		 x0=Res_juveniles$Year.Phi+Shift, lwd=2, col='#00A087')

points(Res_juveniles$Phi.cjs.median~I(Res_juveniles$Year.Phi-Shift),
       las=1, pch=19, col='#E64B35')   
segments(y0=Res_juveniles$Phi.cjs.lci,
         y1=Res_juveniles$Phi.cjs.uci,
		 x0=Res_juveniles$Year.Phi-Shift, lwd=2, col='#E64B35')

#Panel B adults
plot(Res_adults$Phi.rdma.median~Res_adults$Year.Phi,
     ylim=c(0, 1), las=1, ylab='',
	 xlab='Year', type='n')  

abline(v=2004:2020, lty=3, col=grey(0.5))
abline(v=seq(2005, 2020, by=5), lty=2, col=grey(0.5))

abline(h=seq(0, 1, by=0.1), lty=3, col=grey(0.5))
abline(h=seq(0, 1, by=0.5), lty=2, col=grey(0.5))

points(Res_adults$Phi.rdma.median~I(Res_adults$Year.Phi+Shift), 
       las=1, pch=23, bg='#00A087', col='#00A087')   

segments(y0=Res_adults$Phi.rdma.lci,
         y1=Res_adults$Phi.rdma.uci,
		 x0=Res_adults$Year.Phi+Shift, lwd=2, col='#00A087')

points(Res_adults$Phi.cjs.median~I(Res_adults$Year.Phi-Shift),
       las=1, pch=19 , col='#E64B35')   
segments(y0=Res_adults$Phi.cjs.lci,
         y1=Res_adults$Phi.cjs.uci,
		 x0=Res_adults$Year.Phi-Shift, lwd=2, col='#E64B35')
dev.off()
```

### 3.4 Summary statistics used in the paper
#### RDMa model results

```{r}
# average survival probability of adults
print(quantile(rdma.godwits$sims.list$phi.average[,2],
         c(0.025, 0.5, 0.975)), digits=3)

# average survival probability of juveniles
print(quantile(rdma.godwits$sims.list$phi.average[,1],
         c(0.025, 0.5, 0.975)), digits=3)

# slope over time in adult survival probability 		 
print(quantile(rdma.godwits$sims.list$phi.beta[,2],
         c(0.025, 0.5, 0.975)), digits=3)

# slope over time in juvenile survival probability 		 
print(quantile(rdma.godwits$sims.list$phi.beta[,1],
         c(0.025, 0.5, 0.975)), digits=3)
	   		 
# average adult resighting probability adults
print(1-exp(c(-quantile(rdma.godwits$sims.list$lambda.y,
       c(0.025, 0.5, 0.975)))), digits=3)

# average expected number of sightings adults   
print(quantile(rdma.godwits$sims.list$lambda.y,
       c(0.025, 0.5, 0.975)), digits=3)

# Average 1st year juvenile resighting rate
print(1-exp(-plogis(quantile(rdma.godwits$sims.list$lambda.y+
             rdma.godwits$sims.list$lambda.age[,1],
			 c(0.025, 0.5, 0.975)))), digits=3)

# Average 2nd year juvenile resighting rate
print(1-exp(-plogis(quantile(rdma.godwits$sims.list$lambda.y+
             rdma.godwits$sims.list$lambda.age[,2],
			 c(0.025, 0.5, 0.975)))), digits=3)

# Average adults marked as chicks resighting rate
print(1-exp(-plogis(quantile(rdma.godwits$sims.list$lambda.y+
             rdma.godwits$sims.list$lambda.age[,3],
			 c(0.025, 0.5, 0.975)))), digits=3)


# and theta
print(plogis(quantile(rdma.godwits$sims.list$mean.theta,
           c(0.025, 0.5, 0.975))), digits=3)
# sigma theta
print(sd(rdma.godwits$sims.list$mean.theta))
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

