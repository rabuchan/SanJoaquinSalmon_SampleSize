######## Example of implementation: explore expected performance for single sample size across varying parameter sets

## Analysis steps
# Step 0. Load in functions file.
# Step 1. Define parameter values
# Step 2. Define matrix of parameter sets - assumes the possibility of simulating sample size for more than one parameter set.
# Step 3. Define candidate release sizes and number of simulations
# Step 4. Define vectors of parameter names used in likelihood components (optional)
# Step 5. Define vectors of derived parameter names used in likelihood components (optional)
# Step 6. Define expressions for detection history probability parameterizations for multinomial likelihood (optional)
# Step 7. Simulate data and estimate parameters using maximum likelihood (use analytical solutions for method-of-moment estimators = MLE)
# Step 8. Compile maximum likelihood estimates into matrix
# Step 9. Compute criteria performance summaries for criteria C1-C4
# Step 10. Examine results

################################################################################################################

### Step 0. Load in functions file.

# Define directory where "SampleSizeSimulations_functions.r" file is stored
functions_dir<-"D:/Documents/CVPIA/LongtermStudyDesign/SanJoaquinSalmon_SampleSize"

# Source functions file.
source(file.path(functions_dir,"SampleSizeSimulations_functions.r"))





### Step 1. Define parameter values

# Parameters = sR, sR0, sA1, sA2, sB, sF, psiA1, psiA2, sA1.supp, pA1a, pA1b, pA2a, pA2b, pB1a, pB1b, pF1a, pF1b, pG2a, pG2b
#   where sA1 = survival from site A1 to A2 for fish from primary release (i.e., Durham Ferry release)
#     and sA1.supp = survival from supplemental release location to A2 for fish from supplemental release

# derived parameters = sR, sA, pA1, pA2, pB1, pF1, pG2

# Note: each parameter is assigned a vector of one or more candidate values.
sR<-c(0.01, 0.05, 0.1) # this defines colum for sR - will overwrite below
sR0<-c(0.2, 0.7, 0.9)
sA1<-c(0.1, 0.4, 0.5)
sA2<-c(0.02, 0.1, 0.16)
sB<-c(0.01, 0.1, 0.16)
sF<-c(0.01, 0.05)
psiA1<-c(0.4, 0.6, 0.98)
psiA2<-c(0.75, 0.95)
sA1.supp<-sA1 # will equate sA1.supp to sA1 in each parameter set
pA1a<-0.9
pA1b<-0.9
pA2a<-c(0.85)
pA2b<-c(0.85)
pB1a<-c(0.9)
pB1b<-c(0.9)
pF1a<-c(0.85)
pF1b<-c(0.85)
pG2a<-c(0.5,0.9)
pG2b<-pG2a # pG2b will be set equal to pG2a in each parameter set

# will vary sR0, sA1, sA2, sB together to define three survival scenarios: low, medium, high



### Step 2. Define matrix of parameter sets - assumes the possibility of simulating sample size for more than one parameter set.

# define character vector of all parameters, including sR
all_par_names<-c("sR","sR0","sA1","sA2","sB","sF","psiA1","psiA2","sA1.supp","pA1a","pA1b","pA2a","pA2b","pB1a","pB1b","pF1a","pF1b","pG2a","pG2b")

# create matrix of all combinations, leaving spots empty for the parameters that vary together
pars_mat<-expand.grid(sR,NA,NA,NA,NA,sF,psiA1,psiA2,NA,pA1a,pA1b,pA2a,pA2b,pB1a,pB1b,pF1a,pF1b,pG2a,NA)
names(pars_mat)<-all_par_names 

# vary sR0, sA1, sA2, sB together to define four survival scenarios: very low, low, medium, high
for(i in 1:length(sR)) pars_mat[which(pars_mat$sR==sR[i]),c("sR0","sA1","sA2","sB")]=rep(c(sR0[i],sA1[i],sA2[i],sB[i]),each=length(which(pars_mat$sR==sR[i])))

# equate sA1.supp and sA1
pars_mat$sA1.supp<-pars_mat$sA1

# equate pG2a and pG2b
pars_mat$pG2b<-pars_mat$pG2a

# update value of sR: is derived from other parameters
pars_mat[,"sR"]<-apply(pars_mat,1,get.sR)

# define overall conditional detection probabilities at dual arrays
pars_mat$pA1 <- 1 - (1-pars_mat$pA1a)*(1-pars_mat$pA1b)
pars_mat$pA2 <- 1 - (1-pars_mat$pA2a)*(1-pars_mat$pA2b)
pars_mat$pB1 <- 1 - (1-pars_mat$pB1a)*(1-pars_mat$pB1b)
pars_mat$pF1 <- 1 - (1-pars_mat$pF1a)*(1-pars_mat$pF1b)
pars_mat$pG2 <- 1 - (1-pars_mat$pG2a)*(1-pars_mat$pG2b)




### Step 3. Define candidate release sizes and number of simulations
R_DF<-324
R_supp<-324
N_sim<-20  # 1000 or 5000 would be better but takes longer



### Step 4. Define vectors of parameter names used in likelihood components (optional)
# Note: not required because built into sim.supp.seeds.pari.R1.R2()
{
  pars.full.DF<-c("sR0","sA1","sA2","sB","sF","psiA1","psiA2","pA1","pA2","pB1","pF1","pG2")
  pars.full.supp<-c("sA1.supp","sA2","sF","psiA2","pA2","pF1","pG2")
  pars.dual.A1<-c("pA1a","pA1b")
  pars.dual.A2<-c("pA2a","pA2b")
  pars.dual.B1<-c("pB1a","pB1b")
  pars.dual.F1<-c("pF1a","pF1b")
  pars.dual.G2<-c("pG2a","pG2b")  
}





### Step 5. Define vectors of derived parameter names used in likelihood components (optional)
# Note: not required because built into sim.supp.seeds.pari.R1.R2()
{
  vars.full.DF<-c()
  vars.full.DF["qA1"] <- "1 - pA1"
  vars.full.DF["qA2"] <- "1 - pA2"
  vars.full.DF["qB1"] <- "1 - pB1"
  vars.full.DF["qF1"] <- "1 - pF1"
  vars.full.DF["qG2"] <- "1 - pG2"
  vars.full.DF["psiB1"] <- "1 - psiA1"
  vars.full.DF["psiF2"] <- "1 - psiA2"
  vars.full.DF["xF1"] <- "1 - sF + sF*qG2"
  vars.full.DF["xA2"] <- "1 - sA2 + sA2*qG2"
  vars.full.DF["xB1"] <- "1 - sB + sB*qG2"
  vars.full.DF["xA1"] <- "1 - sA1 + sA1*(psiA2*qA2*xA2 + psiF2*qF1*xF1)"
  vars.full.DF["xR"] <- "1 - sR0 + sR0*(psiA1*qA1*xA1 + psiB1*qB1*xB1)"
  vars.full.DF["gammaA1G2"] <- "sA1*(psiA2*qA2*sA2 + psiF2*qF1*sF)"
  vars.full.DF["gammaRG2"] <- "sR0*(psiA1*qA1*gammaA1G2 + psiB1*qB1*sB)"
  
  vars.full.supp<-c()
  vars.full.supp["qA2"] <- "1 - pA2"
  vars.full.supp["qF1"] <- "1 - pF1"
  vars.full.supp["qG2"] <- "1 - pG2"
  vars.full.supp["psiF2"] <- "1 - psiA2"
  vars.full.supp["xF1"] <- "1 - sF + sF*qG2"
  vars.full.supp["xA2"] <- "1 - sA2 + sA2*qG2"
  vars.full.supp["xSupp"] <- "1 - sA1.supp + sA1.supp*(psiA2*qA2*xA2 + psiF2*qF1*xF1)"
  
  vars.dual.A1<-c()
  vars.dual.A1["pA1"] <- "1 - (1 - pA1a)*(1 - pA1b)"
  
  vars.dual.A2<-c()
  vars.dual.A2["pA2"] <- "1 - (1 - pA2a)*(1 - pA2b)"
  
  vars.dual.B1<-c()
  vars.dual.B1["pB1"] <- "1 - (1 - pB1a)*(1 - pB1b)"
  
  vars.dual.F1<-c()
  vars.dual.F1["pF1"] <- "1 - (1 - pF1a)*(1 - pF1b)"
  
  vars.dual.G2<-c()
  vars.dual.G2["pG2"] <- "1 - (1 - pG2a)*(1 - pG2b)"  
}





### Step 6. Define expressions for detection history probability parameterizations for multinomial likelihood (optional)
# Note: not required because built into sim.supp.seeds.pari.R1.R2()
{
  like.full.DF<-c()
  like.full.DF["DF A1 A2 G2"] <- "sR0*psiA1*pA1*sA1*psiA2*pA2*sA2*pG2"
  like.full.DF["DF A1 A2 x"] <- "sR0*psiA1*pA1*sA1*psiA2*pA2*xA2"
  like.full.DF["DF A1 F1 G2"] <- "sR0*psiA1*pA1*sA1*psiF2*pF1*sF*pG2"
  like.full.DF["DF A1 F1 x"] <- "sR0*psiA1*pA1*sA1*psiF2*pF1*xF1"
  like.full.DF["DF A1 G2"] <- "sR0*psiA1*pA1*gammaA1G2*pG2"
  like.full.DF["DF A1 x"] <- "sR0*psiA1*pA1*xA1"
  like.full.DF["DF B1 G2"] <- "sR0*psiB1*pB1*sB*pG2"
  like.full.DF["DF B1 x"] <- "sR0*psiB1*pB1*xB1"
  like.full.DF["DF A2 G2"] <- "sR0*psiA1*qA1*sA1*psiA2*pA2*sA2*pG2"
  like.full.DF["DF A2 x"] <- "sR0*psiA1*qA1*sA1*psiA2*pA2*xA2"
  like.full.DF["DF F1 G2"] <- "sR0*psiA1*qA1*sA1*psiF2*pF1*sF*pG2"
  like.full.DF["DF F1 x"] <- "sR0*psiA1*qA1*sA1*psiF2*pF1*xF1"
  like.full.DF["DF G2"] <- "gammaRG2*pG2"
  like.full.DF["DF x"] <- "xR"
  
  like.full.supp<-c()
  like.full.supp["supp A2 G2"] <- "sA1.supp*psiA2*pA2*sA2*pG2"
  like.full.supp["supp A2 x"] <- "sA1.supp*psiA2*pA2*xA2"
  like.full.supp["supp F1 G2"] <- "sA1.supp*psiF2*pF1*sF*pG2"
  like.full.supp["supp F1 x"] <- "sA1.supp*psiF2*pF1*xF1"
  like.full.supp["supp G2"] <- "sA1.supp*(psiA2*qA2*sA2 + psiF2*qF1*sF)*pG2"
  like.full.supp["supp x"] <- "xSupp"
  
  like.dual.A1<-c()
  like.dual.A1["A1 AB"] <- "pA1a*pA1b/pA1"
  like.dual.A1["A1 A0"] <- "pA1a*(1-pA1b)/pA1"
  like.dual.A1["A1 B0"] <- "(1-pA1a)*pA1b/pA1"
  
  like.dual.A2<-c()
  like.dual.A2["A2 AB"] <- "pA2a*pA2b/pA2"
  like.dual.A2["A2 A0"] <- "pA2a*(1-pA2b)/pA2"
  like.dual.A2["A2 B0"] <- "(1-pA2a)*pA2b/pA2"
  
  like.dual.F1<-c()
  like.dual.F1["F1 AB"] <- "pF1a*pF1b/pF1"
  like.dual.F1["F1 A0"] <- "pF1a*(1-pF1b)/pF1"
  like.dual.F1["F1 B0"] <- "(1-pF1a)*pF1b/pF1"
  
  like.dual.B1<-c()
  like.dual.B1["B1 AB"] <- "pB1a*pB1b/pB1"
  like.dual.B1["B1 A0"] <- "pB1a*(1-pB1b)/pB1"
  like.dual.B1["B1 B0"] <- "(1-pB1a)*pB1b/pB1"
  
  like.dual.G2<-c()
  like.dual.G2["G2 AB"] <- "pG2a*pG2b/pG2"
  like.dual.G2["G2 A0"] <- "pG2a*(1-pG2b)/pG2"
  like.dual.G2["G2 B0"] <- "(1-pG2a)*pG2b/pG2"  
}




### Step 7. Simulate data and estimate parameters using maximum likelihood (use analytical solutions for method-of-moment estimators = MLE)
dat<-apply(pars_mat[,-1],1,sim.supp.seeds.pari.R1.R2, R1=R_DF, R2=R_supp, nSim=N_sim)

# full function call naming arguments for DF.pars, supp.pars, etc.:
 # dat<-apply(pars_mat[,-1],1,sim.supp.seeds.pari.R1.R2, R1=R_DF, R2=R_supp, nSim=N_sim, pars.full.DF,pars.full.supp,pars.dual.A1,pars.dual.A2,pars.dual.B1,pars.dual.F1,pars.dual.G2, vars.full.DF,vars.full.supp,vars.dual.A1,vars.dual.A2,vars.dual.B1,vars.dual.F1,vars.dual.G2, like.full.DF,like.full.supp,like.dual.A1,like.dual.A2,like.dual.B1,like.dual.F1,like.dual.G2)

# Note: 
 # Can also estimate MLE using optim via sim.supp.mle.pari.R1.R2() - but code does not reliably generate estimates. 
 # Estimating via "seeds" (actually = MoM=MLE) is more robust - this is the approach used here via sim.supp.seeds.pari.R1.R2





### Step 8. Compile maximum likelihood estimates into matrix

# define the parameters to be estimate
pars_est<-all_par_names[-1]
pars_est<-c(pars_est,c("sA","sR"))

est_pars_mat<-lapply(c(1:ncol(dat)),get.est.mat.wrap,xx=dat,Nsim=N_sim, Pnames=pars_est)
# est_pars_mat is a list whose length is equal to the number of parameter sets (=nrow(pars_mat))
# each element of est_pars_mat is a matrix with N_sim+1 rows and columns = parameters estimated (ncol = length of pars_est)
#   -- first row of est_pars_mat defines the true value for each parameter
#   -- remaining rows contain the MLEs for the 1:N_sim simulations (row = simulation)



### Step 9. Compute criteria performance summaries for criteria C1-C4
C1_out<-do.call(rbind,lapply(est_pars_mat,function(xx) apply(xx,2,get.C1)))
C2_out<-do.call(rbind,lapply(est_pars_mat,function(xx) apply(xx,2,get.C2)))
C3_out<-do.call(rbind,lapply(est_pars_mat,function(xx) apply(xx,2,get.C3)))
C4_out<-do.call(rbind,lapply(est_pars_mat,function(xx) apply(xx,2,get.C4)))
mean_out<-do.call(rbind,lapply(est_pars_mat,function(xx) apply(xx,2,get.mean)))
max_out<-do.call(rbind,lapply(est_pars_mat,function(xx) apply(xx,2,get.max)))

# Note: Each of C1_out, C2_out, ..., mean_out, max_out is a matrix with nrows = number of parameter sets and ncol = number of parameters estimated
 # -- row = parameter set
 # -- column = parameter estimated
 # -- entry:
 # ---- C1: proportion of simulations that generate estimates for parameter
 # ---- C2: proportion of simulations that generate estimates > cutoff value
 # ---- C3: standard deviation of simulation distribution of parameter estimates
 # ---- C4: bias of parameter estimates
 # ---- mean: mean of simulation distribution of parameter estimates
 # ---- max: maximum of simulation distribution of parameter estimates




### Step 10. Examine results
head(est_pars_mat[[1]]) # the true values of parameters from parameter set 1 (row 1) and the parameter MLEs for the first 5 simulations for parameter set 1 (rows 2-6)
C1_out[1,] # the proportion of all simulations for which the parameter was estimable for parameter set 1
C2_out[1,] # the proportion of all simulations for which the parameter MLE was not greater than 1.1 for parameter set 1
C3_out[1,] # the proportion of all simulations for which the parameter MLE standard error was not greater than 0.5 for parameter set 1
C4_out[1,] # the proportion of all simulations for which the bias in parameter MLE was not greater than 0.5 for parameter set 1





