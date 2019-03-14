######################################
# This file contains functions to simulate CJS data with misidentification
# and also some helper functions for formatting these data and generating initial values for JAGS
# Some of these functions are the same as in Kery and Schaub (2011) book
# while other were developed specifically for CJSm models

withAutoprint({
cat('sourcing functions...') 
#####################################################
# 1.The function generally assumes that there are several sightings possible, so it converts P to Lambda and does the simulation
  cat('\nsimul.cjs.multiple.sightings()' )
simul.cjs.multiple.sightings <- function(PHI, P, THETA, marked, CH=NULL){
   # in this function I introduce multiple resightings.
   
   # Note, that if CH is not null function will not look at PHI and P and will take supplied CH and just add misreadings..
   
   if (is.null(CH)) {
   # first we convert P values to lambda
   Lambda_trials<-apply(P, c(1,2), FUN=function(x) -log(1-x))
   
   n.occasions <- dim(PHI)[2] + 1
   CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- rep(1:length(marked), marked[1:length(marked)])
   # Fill the CH matrix
   for (i in 1:sum(marked)){
      CH[i, mark.occ[i]] <- -1       # Write an - 1 at the release occasion and correct it to one later.
      if (mark.occ[i]==n.occasions) next
      for (t in (mark.occ[i]+1):n.occasions){
         # Bernoulli trial: does individual survive occasion?
         sur <- rbinom(1, 1, PHI[i,t-1])
         if (sur==0) break		# If dead, move to next individual 
         # Poisson trial: how many times was individual resighted 
         rp <- rpois(1, Lambda_trials[i,t-1])
         if (rp>0) CH[i,t] <- rp		 
         } #t
      } #i
   }
   
   if (min(CH)==0) { # here I make the freshly marked -1
      First<-apply(CH, 1, FUN=function(x) min(which(x!=0)))
   for (ii in 1:nrow(CH)) CH[ii, First[ii]]<--1
   }

   CH_no_misr<-CH
   ##### 
   # now we introduce misreadings.
   #    
   for (t in 2:n.occasions){
     # select the resigthted individuals
     Resighted<-which(CH[,t]>0)
	 # Generate misreadings
     Misread<-sapply(Resighted ,  FUN=function(x)  rbinom(1, CH[x,t], 1-THETA[x,t-1])) #
     if (sum(Misread>0)) {
        CH[Resighted,t]<-CH[Resighted,t] - Misread
        # assign misreading to any other number that was already in use
        # 1. which individuals are marked   
        Already_marked<-sum(marked[1:(t-1)])
        # 2. randomly add some ones there..
        Assinged<-sample.int(Already_marked, sum(Misread), replace=TRUE)
        # give them one
        Rle<-rle(sort(Assinged))
        CH[Rle$values,t]<-CH[Rle$values,t]+Rle$lengths
     }
   }
   CH<-apply(CH, c(1,2), FUN=function(x) ifelse(x==-1,1, x))
   CH_no_misr<-apply(CH_no_misr, c(1,2), FUN=function(x) ifelse(x==-1,1, x))
   return(list(CH=CH, CH_no_misr=CH_no_misr))
}
  cat('   Done')

#####################################################
# 2. This function from Schaub and Kery 2011 helps to set up initital under single readings
  cat('\nknown.state.cjs()' )
known.state.cjs <- function(ch){
   state <- ch
   for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }
   cat('   Done')

#####################################################
# 3. This function from Schaub and Kery 2011 helps to create a matrix of initial values for latent state z
  cat('\ncjs.init.z()' )
cjs.init.z <- function(ch,f){
   for (i in 1:dim(ch)[1]){
      if (sum(ch[i,])==1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,f[i]:n2] <- NA
      }
   for (i in 1:dim(ch)[1]){
   ch[i,1:f[i]] <- NA
   }
   return(ch)
   }
   cat('   Done')

#####################################################
# 4. This function sets up inititals under multiple readings values
  cat('\nknown.state.cjs.mult()' )
known.state.cjs.mult<- function(ch) {
  state <- sign(ch)
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]>0))
    n2 <- max(which(ch[i,]>0))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}
  cat('   Done')

  cat('\nget.first()' )

#####################################################
# 5. This function from from Schaub and Kery (2011) siply finds the first positive values in the vector.
get.first <- function(x) min(which(x!=0))
  cat('   Done')

cat('\nyou are all set \n')
}, echo=FALSE)
