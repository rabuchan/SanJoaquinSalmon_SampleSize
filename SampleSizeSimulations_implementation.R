######## Example of implementation

## Analysis steps
# Step 1. Define parameter values
# Step 2. Define matrix of parameter sets - assumes the possibility of simulating sample size for more than one parameter set.
# Step 3. Define candidate release sizes and number of simulations
# Step 4. Define vectors of parameter names used in likelihood components (optional)
# Step 5. Define vectors of derived parameter names used in likelihood components (optional)
# Step 6. Define expressions for detection history probability parameterizations for multinomial likelihood
# Step 7. Simulate data and estimate parameters using maximum likelihood (use analytical solutions for method-of-moment estimators = MLE)
# Step 8. Compile maximum likelihood estimates into matrix
# Step 9. Compute criteria performance summaries for criteria C1-C4
# Step 10. Examine results


In order to implement this, I need to find pars.mat.2 from the 2016 sample size simulations, but it appears to have come from on high
Maybe it is defined in the 2015 sample size simulations?
  see "MoM simulations.r" from 2015



### Step 1. Define parameter values





### Step 2. Define matrix of parameter sets - assumes the possibility of simulating sample size for more than one parameter set.
pars_mat<-....





### Step 3. Define candidate release sizes and number of simulations
R_DF<-450
R_supp<-450
N_sim<-1000





### Step 4. Define vectors of parameter names used in likelihood components (optional)
# Note: not required because built into sim.supp.seeds.pari.R1.R2()
{
  pars.full.DF=c("sR0","sA1","sA2","sB","sF","psiA1","psiA2","pA1","pA2","pB1","pF1","pG2")
  pars.full.supp=c("sA1.supp","sA2","sF","psiA2","pA2","pF1","pG2")
  pars.dual.A1=c("pA1a","pA1b")
  pars.dual.A2=c("pA2a","pA2b")
  pars.dual.B1=c("pB1a","pB1b")
  pars.dual.F1=c("pF1a","pF1b")
  pars.dual.G2=c("pG2a","pG2b")  
}





### Step 5. Define vectors of derived parameter names used in likelihood components (optional)
# Note: not required because built into sim.supp.seeds.pari.R1.R2()
{
  vars.full.DF=c()
  vars.full.DF["qA1"] = "1 - pA1"
  vars.full.DF["qA2"] = "1 - pA2"
  vars.full.DF["qB1"] = "1 - pB1"
  vars.full.DF["qF1"] = "1 - pF1"
  vars.full.DF["qG2"] = "1 - pG2"
  vars.full.DF["psiB1"] = "1 - psiA1"
  vars.full.DF["psiF2"] = "1 - psiA2"
  vars.full.DF["xF1"] = "1 - sF + sF*qG2"
  vars.full.DF["xA2"] = "1 - sA2 + sA2*qG2"
  vars.full.DF["xB1"] = "1 - sB + sB*qG2"
  vars.full.DF["xA1"] = "1 - sA1 + sA1*(psiA2*qA2*xA2 + psiF2*qF1*xF1)"
  vars.full.DF["xR"] = "1 - sR0 + sR0*(psiA1*qA1*xA1 + psiB1*qB1*xB1)"
  vars.full.DF["gammaA1G2"] = "sA1*(psiA2*qA2*sA2 + psiF2*qF1*sF)"
  vars.full.DF["gammaRG2"] = "sR0*(psiA1*qA1*gammaA1G2 + psiB1*qB1*sB)"
  
  vars.full.supp=c()
  vars.full.supp["qA2"] = "1 - pA2"
  vars.full.supp["qF1"] = "1 - pF1"
  vars.full.supp["qG2"] = "1 - pG2"
  vars.full.supp["psiF2"] = "1 - psiA2"
  vars.full.supp["xF1"] = "1 - sF + sF*qG2"
  vars.full.supp["xA2"] = "1 - sA2 + sA2*qG2"
  vars.full.supp["xSupp"] = "1 - sA1.supp + sA1.supp*(psiA2*qA2*xA2 + psiF2*qF1*xF1)"
  
  vars.dual.A1=c()
  vars.dual.A1["pA1"] = "1 - (1 - pA1a)*(1 - pA1b)"
  
  vars.dual.A2=c()
  vars.dual.A2["pA2"] = "1 - (1 - pA2a)*(1 - pA2b)"
  
  vars.dual.B1=c()
  vars.dual.B1["pB1"] = "1 - (1 - pB1a)*(1 - pB1b)"
  
  vars.dual.F1=c()
  vars.dual.F1["pF1"] = "1 - (1 - pF1a)*(1 - pF1b)"
  
  vars.dual.G2=c()
  vars.dual.G2["pG2"] = "1 - (1 - pG2a)*(1 - pG2b)"  
}





### Step 6. Define expressions for detection history probability parameterizations for multinomial likelihood
# Note: not required because built into sim.supp.seeds.pari.R1.R2()
{
  like.full.DF=c()
  like.full.DF["DF A1 A2 G2"] = "sR0*psiA1*pA1*sA1*psiA2*pA2*sA2*pG2"
  like.full.DF["DF A1 A2 x"] = "sR0*psiA1*pA1*sA1*psiA2*pA2*xA2"
  like.full.DF["DF A1 F1 G2"] = "sR0*psiA1*pA1*sA1*psiF2*pF1*sF*pG2"
  like.full.DF["DF A1 F1 x"] = "sR0*psiA1*pA1*sA1*psiF2*pF1*xF1"
  like.full.DF["DF A1 G2"] = "sR0*psiA1*pA1*gammaA1G2*pG2"
  like.full.DF["DF A1 x"] = "sR0*psiA1*pA1*xA1"
  like.full.DF["DF B1 G2"] = "sR0*psiB1*pB1*sB*pG2"
  like.full.DF["DF B1 x"] = "sR0*psiB1*pB1*xB1"
  like.full.DF["DF A2 G2"] = "sR0*psiA1*qA1*sA1*psiA2*pA2*sA2*pG2"
  like.full.DF["DF A2 x"] = "sR0*psiA1*qA1*sA1*psiA2*pA2*xA2"
  like.full.DF["DF F1 G2"] = "sR0*psiA1*qA1*sA1*psiF2*pF1*sF*pG2"
  like.full.DF["DF F1 x"] = "sR0*psiA1*qA1*sA1*psiF2*pF1*xF1"
  like.full.DF["DF G2"] = "gammaRG2*pG2"
  like.full.DF["DF x"] = "xR"
  
  like.full.supp["supp A2 G2"] = "sA1.supp*psiA2*pA2*sA2*pG2"
  like.full.supp["supp A2 x"] = "sA1.supp*psiA2*pA2*xA2"
  like.full.supp["supp F1 G2"] = "sA1.supp*psiF2*pF1*sF*pG2"
  like.full.supp["supp F1 x"] = "sA1.supp*psiF2*pF1*xF1"
  like.full.supp["supp G2"] = "sA1.supp*(psiA2*qA2*sA2 + psiF2*qF1*sF)*pG2"
  like.full.supp["supp x"] = "xSupp"
  
  like.dual.A1=c()
  like.dual.A1["A1 AB"] = "pA1a*pA1b/pA1"
  like.dual.A1["A1 A0"] = "pA1a*(1-pA1b)/pA1"
  like.dual.A1["A1 B0"] = "(1-pA1a)*pA1b/pA1"
  
  like.dual.A2=c()
  like.dual.A2["A2 AB"] = "pA2a*pA2b/pA2"
  like.dual.A2["A2 A0"] = "pA2a*(1-pA2b)/pA2"
  like.dual.A2["A2 B0"] = "(1-pA2a)*pA2b/pA2"
  
  like.dual.F1=c()
  like.dual.F1["F1 AB"] = "pF1a*pF1b/pF1"
  like.dual.F1["F1 A0"] = "pF1a*(1-pF1b)/pF1"
  like.dual.F1["F1 B0"] = "(1-pF1a)*pF1b/pF1"
  
  like.dual.B1=c()
  like.dual.B1["B1 AB"] = "pB1a*pB1b/pB1"
  like.dual.B1["B1 A0"] = "pB1a*(1-pB1b)/pB1"
  like.dual.B1["B1 B0"] = "(1-pB1a)*pB1b/pB1"
  
  like.dual.G2=c()
  like.dual.G2["G2 AB"] = "pG2a*pG2b/pG2"
  like.dual.G2["G2 A0"] = "pG2a*(1-pG2b)/pG2"
  like.dual.G2["G2 B0"] = "(1-pG2a)*pG2b/pG2"  
}




### Step 7. Simulate data and estimate parameters using maximum likelihood (use analytical solutions for method-of-moment estimators = MLE)
dat<-apply(pars_mat[,-1],1,sim.supp.seeds.pari.R1.R2,R1=R_DF,R2=R_supp,nsim=N_sim)

# full function call naming arguments for DF.pars, supp.pars, etc.:
 # dat<-apply(pars_mat[,-1],1,sim.supp.seeds.pari.R1.R2,R1=R_DF,R2=R_supp,nsim=N_sim,pars.full.DF,pars.full.supp,pars.dual.A1,pars.dual.A2,pars.dual.B1,pars.dual.F1,pars.dual.G2, vars.full.DF,vars.full.supp,vars.dual.A1,vars.dual.A2,vars.dual.B1,vars.dual.F1,vars.dual.G2, like.full.DF,like.full.supp,like.dual.A1,like.dual.A2,like.dual.B1,like.dual.F1,like.dual.G2)

# Note: 
 # Can also estimate MLE using optim via sim.supp.mle.pari.R1.R2() - but code does not reliably generate estimates. 
 # Estimating via "seeds" (actually = MoM=MLE) is more robust - this is the approach used here via sim.supp.seeds.pari.R1.R2





### Step 8. Compile maximum likelihood estimates into matrix
est.par.mat<-lapply(c(1:ncol(dat)),get.est.mat.wrap,xx=dat,Nsim=N_sim)




### Step 9. Compute criteria performance summaries for criteria C1-C4
C1.out<-do.call(rbind,lapply(est.par.mat,function(xx) apply(xx,2,get.C1)))
C2.out<-do.call(rbind,lapply(est.par.mat,function(xx) apply(xx,2,get.C2)))
C3.out<-do.call(rbind,lapply(est.par.mat,function(xx) apply(xx,2,get.C3)))
C4.out<-do.call(rbind,lapply(est.par.mat,function(xx) apply(xx,2,get.C4)))
mean.out<-do.call(rbind,lapply(est.par.mat,function(xx) apply(xx,2,get.mean)))
max.out<-do.call(rbind,lapply(est.par.mat,function(xx) apply(xx,2,get.max)))




### Step 10. Examine results




