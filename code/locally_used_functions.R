######################################
# This file contains a simple function to combine simulation results.

withAutoprint({
cat('sourcing make_two_cells()') 

make_two_cells<-function(filename=NULL, N_Models_to_use=100) {
# the idea is that this function takes output from the function that did combine and wll just make a results table
load(filename)
  Phi<-cur_parameters$Phi
  p<-cur_parameters$p
  theta<-cur_parameters$theta
  RHat<-cur_parameters$RHat
  #N_cc<-  length(which(apply(Rhat_cc, 1,max)<=RHat))
  #N_ccc<-  length(which(apply(Rhat_ccc, 1,max)<=RHat))
  Converged_cc<-na.omit(which(apply(Rhat_cc, 1,max)<=RHat)[1:N_Models_to_use])
  Converged_ccc<-na.omit(which(apply(Rhat_ccc, 1,max)<=RHat))[1:N_Models_to_use]
  N_cc<-length(Converged_cc)
  N_ccc<-length(Converged_ccc)
  Res<-c()

  if (length(Converged_cc>1)) {
  Mean.c.c<-as.data.frame(Mean.c.c)
  cc_CI_width<-as.data.frame(cc_CI_width)
  Means_cc<-colMeans(Mean.c.c[Converged_cc,])
  Res_cc<-data.frame(
     'type'='naive',
     'Phi'=Phi,
	 'p'=p,
	 'theta'=theta,
	 'RHat'=RHat,
	 'Phi.mean'=Means_cc[1],
	 'Phi.lci'=	 quantile(Mean.c.c$mean.phi[Converged_cc], 0.025),
	 'Phi.uci'=	 quantile(Mean.c.c$mean.phi[Converged_cc], 0.975),
	 'Phi.bias'=  (Means_cc[1]-Phi),
	 'Phi.sc.MSE'=  sum((Mean.c.c$mean.phi[Converged_cc]-Phi)^2)/((N_cc-1)*Phi^2),
	 'Phi.CI.width'= mean(cc_CI_width$mean.phi[Converged_cc]),
	 'Phi.CI.coverage'= sum(Cov.c.c[Converged_cc,1])/N_cc,
     'Phi.mean.n.eff'=mean(N_eff_cc[Converged_cc,1]),
	 'p.mean'=Means_cc[2],
	 'p.lci'=	 quantile(Mean.c.c$mean.p[Converged_cc], 0.025),
	 'p.uci'=	 quantile(Mean.c.c$mean.p[Converged_cc], 0.975),
	 'p.bias'=  (Means_cc[2]-p),
	 'p.sc.MSE'=  sum((Mean.c.c$mean.p[Converged_cc]-p)^2)/((N_cc-1)*p^2),
	 'p.CI.width'= mean(cc_CI_width$mean.p[Converged_cc]),
	 'p.CI.coverage'= sum(Cov.c.c[Converged_cc,2])/N_cc,
     'p.mean.n.eff'=mean(N_eff_cc[Converged_cc,2]),
	 'theta.mean'=NA,
	 'theta.lci'=	 NA,
	 'theta.uci'=	 NA,
	 'theta.bias'=  NA,
	 'theta.sc.MSE'=  NA,
	 'theta.CI.width'= NA,
	 'theta.CI.coverage'= NA,
     'theta.mean.n.eff'=NA,
	 'N_models'=N_cc

	)
  Res<-rbind(Res, Res_cc)
  }
  if (length(Converged_ccc>1)) {
  Mean.c.c.c<-as.data.frame(Mean.c.c.c)
  ccc_CI_width<-as.data.frame(ccc_CI_width)
  Means_ccc<-colMeans(Mean.c.c.c[Converged_ccc,])

  Res_ccc<-data.frame(
     'type'='CJSm',
     'Phi'=Phi,
	 'p'=p,
	 'theta'=theta,
	 'RHat'=RHat,
	 'Phi.mean'=Means_ccc[1],
	 'Phi.lci'=	 quantile(Mean.c.c.c$mean.phi[Converged_ccc], 0.025),
	 'Phi.uci'=	 quantile(Mean.c.c.c$mean.phi[Converged_ccc], 0.975),
	 'Phi.bias'=  (Means_ccc[1]-Phi),
	 'Phi.sc.MSE'=  sum((Mean.c.c.c$mean.phi[Converged_ccc]-Phi)^2)/((N_ccc-1)*Phi^2),
	 'Phi.CI.width'= mean(ccc_CI_width$mean.phi[Converged_ccc]),
	 'Phi.CI.coverage'= sum(Cov.c.c.c[Converged_ccc,1])/N_ccc,
     'Phi.mean.n.eff'=mean(N_eff_ccc[Converged_ccc,1]),
	 'p.mean'=Means_ccc[2],
	 'p.lci'=	 quantile(Mean.c.c.c$mean.p[Converged_ccc], 0.025),
	 'p.uci'=	 quantile(Mean.c.c.c$mean.p[Converged_ccc], 0.975),
	 'p.bias'=  (Means_ccc[2]-p),
	 'p.sc.MSE'=  sum((Mean.c.c.c$mean.p[Converged_ccc]-p)^2)/((N_ccc-1)*p^2),
	 'p.CI.width'= mean(ccc_CI_width$mean.p[Converged_ccc]),
	 'p.CI.coverage'= sum(Cov.c.c.c[Converged_ccc,2])/N_ccc,
     'p.mean.n.eff'=mean(N_eff_ccc[Converged_ccc,2]),
     'theta.mean'=Means_ccc[3],
	 'theta.lci'= quantile(Mean.c.c.c$mean.theta[Converged_ccc], 0.025),
	 'theta.uci'= quantile(Mean.c.c.c$mean.theta[Converged_ccc], 0.975),
	 'theta.bias'= (Means_ccc[3]-theta),
	 'theta.sc.MSE'=  sum((Mean.c.c.c$mean.theta[Converged_ccc]-theta)^2)/((N_ccc-1)*theta^2),
	 'theta.CI.width'= mean(ccc_CI_width$mean.theta[Converged_ccc]),
	 'theta.CI.coverage'= sum(Cov.c.c.c[Converged_ccc,3])/N_ccc,
     'theta.mean.n.eff'=mean(N_eff_ccc[Converged_ccc,3]),
	 'N_models'=N_ccc
	)
   Res<-rbind(Res, Res_ccc)
  }
  return(Res)  
  }

  cat('   Done\n')
}, echo=FALSE)
