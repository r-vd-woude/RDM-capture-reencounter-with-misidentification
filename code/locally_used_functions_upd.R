######################################
# This file contains a simple function to combine simulation results.

withAutoprint({
cat('sourcing make_two_cells()') 

make_two_cells<-function(out_all, Grid_line, N_Models_to_use=100) {
# now the function will take out_all
# and the grid_line - named vector with 
# Phi, P, Theta
cur_out_all<-out_all[out_all$Phi.true==Grid_line$Phi & out_all$P.true==Grid_line$P & out_all$Theta.true==Grid_line$Theta,]
  Phi<-Grid_line$Phi
  p<-Grid_line$P
  theta<-Grid_line$Theta
  RHat<-1.01
  
  Converged_cc<-na.omit(which(cur_out_all$Model=='CJS_naive' &  pmax(cur_out_all$Phi.RHat, cur_out_all$P.RHat, cur_out_all$Deviance.RHat)<=RHat)[1:N_Models_to_use])
  
  Converged_ccc<-na.omit(which(pmax(cur_out_all$Phi.RHat, cur_out_all$P.RHat, cur_out_all$Theta.RHat, cur_out_all$Deviance.RHat)<=RHat)[1:N_Models_to_use])
  N_cc<-length(Converged_cc)
  N_ccc<-length(Converged_ccc)
  Res<-c()

  if (N_cc>1) {
  Res_cc<-data.frame(
     'type'='naive',
     'Phi'=Phi,
	 'p'=p,
	 'theta'=theta,
	 'RHat'=RHat,
	 'Phi.mean'=  mean(cur_out_all$Phi.mean[Converged_cc]),
	 'Phi.lci'=	 quantile(cur_out_all$Phi.mean[Converged_cc], 0.025),
	 'Phi.uci'=	 quantile(cur_out_all$Phi.mean[Converged_cc], 0.975),
	 'Phi.bias'=  (mean(cur_out_all$Phi.mean[Converged_cc])-Phi),
	 'Phi.sc.MSE'=  sum((cur_out_all$Phi.mean[Converged_cc]-Phi)^2)/((N_cc-1)*Phi^2),
	 'Phi.CI.width'= mean(cur_out_all$Phi.CI.width[Converged_cc]),
	 'Phi.CI.coverage'= sum(cur_out_all$Phi.CI.coverage[Converged_cc])/N_cc,
     'Phi.mean.n.eff'= mean(cur_out_all$Phi.n.eff[Converged_cc]),
	 'p.mean'=mean(cur_out_all$P.mean[Converged_cc]),
	 'p.lci'=	 quantile(cur_out_all$P.mean[Converged_cc], 0.025),
	 'p.uci'=	 quantile(cur_out_all$P.mean[Converged_cc], 0.975),
	 'p.bias'=  (mean(cur_out_all$P.mean[Converged_cc])-p),
	 'p.sc.MSE'=  sum((cur_out_all$P.mean[Converged_cc]-p)^2)/((N_cc-1)*p^2),
	 'p.CI.width'= mean(cur_out_all$P.CI.width[Converged_cc]),
	 'p.CI.coverage'= sum(cur_out_all$P.CI.coverage[Converged_cc])/N_cc,
     'p.mean.n.eff'= mean(cur_out_all$P.n.eff[Converged_cc]),
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
  if (N_ccc>1) {

  Res_ccc<-data.frame(
     'type'='CJSm',
     'Phi'=Phi,
	 'p'=p,
	 'theta'=theta,
	 'RHat'=RHat,
	 'Phi.mean'=  mean(cur_out_all$Phi.mean[Converged_ccc]),
	 'Phi.lci'=	 quantile(cur_out_all$Phi.mean[Converged_ccc], 0.025),
	 'Phi.uci'=	 quantile(cur_out_all$Phi.mean[Converged_ccc], 0.975),
	 'Phi.bias'=  (mean(cur_out_all$Phi.mean[Converged_ccc])-Phi),
	 'Phi.sc.MSE'=  sum((cur_out_all$Phi.mean[Converged_ccc]-Phi)^2)/((N_ccc-1)*Phi^2),
	 'Phi.CI.width'= mean(cur_out_all$Phi.CI.width[Converged_ccc]),
	 'Phi.CI.coverage'= sum(cur_out_all$Phi.CI.coverage[Converged_ccc])/N_ccc,
     'Phi.mean.n.eff'= mean(cur_out_all$Phi.n.eff[Converged_ccc]),
	 'p.mean'=mean(cur_out_all$P.mean[Converged_ccc]),
	 'p.lci'=	 quantile(cur_out_all$P.mean[Converged_ccc], 0.025),
	 'p.uci'=	 quantile(cur_out_all$P.mean[Converged_ccc], 0.975),
	 'p.bias'=  (mean(cur_out_all$P.mean[Converged_ccc])-p),
	 'p.sc.MSE'=  sum((cur_out_all$P.mean[Converged_ccc]-p)^2)/((N_ccc-1)*p^2),
	 'p.CI.width'= mean(cur_out_all$P.CI.width[Converged_ccc]),
	 'p.CI.coverage'= sum(cur_out_all$P.CI.coverage[Converged_ccc])/N_ccc,
     'p.mean.n.eff'= mean(cur_out_all$P.n.eff[Converged_ccc]),
     'theta.mean'= mean(cur_out_all$Theta.mean[Converged_ccc]),
	 'theta.lci'= quantile(cur_out_all$Theta.mean[Converged_ccc], 0.025),
	 'theta.uci'= quantile(cur_out_all$Theta.mean[Converged_ccc], 0.975),
	 'theta.bias'= (mean(cur_out_all$Theta.mean[Converged_ccc])-theta),
	 'theta.sc.MSE'=  sum((cur_out_all$Theta.mean[Converged_ccc]-theta)^2)/((N_ccc-1)*theta^2),
	 'theta.CI.width'= mean(cur_out_all$Theta.CI.width[Converged_ccc]),
	 'theta.CI.coverage'= sum(cur_out_all$Theta.CI.coverage[Converged_ccc])/N_ccc,
     'theta.mean.n.eff'=mean(cur_out_all$Theta.n.eff[Converged_ccc]),
	 'N_models'=N_ccc
	)
   Res<-rbind(Res, Res_ccc)
  }
  return(Res)  
  }

  cat('   Done\n')
}, echo=FALSE)
#make_two_cells(out_all, Grid_line)