######## Example of implementation: explore expected performance for varying sample size across for single parameter set

## Analysis steps
# Step 0. Load in functions file.
# Step 1. Define parameter values and derived parameters
# Step 2. Define candidate release sizes and number of simulations
# Step 3. Define matrix of release sizes
# Step 4. Define vectors of parameter names used in likelihood components (optional)
# Step 5. Define vectors of derived parameter names used in likelihood components (optional)
# Step 6. Define expressions for detection history probability parameterizations for multinomial likelihood (optional)
# Step 7. Simulate data and estimate parameters using maximum likelihood (use analytical solutions for method-of-moment estimators = MLE)
# Step 8. Compute criteria performance summaries for criteria C1-C4 (optional)
# Step 9. Summarize results: caluculate minimum release size required to meet the criteria

################################################################################################################


### Step 0. Load in functions file.

# Define directory where "SampleSizeSimulations_functions.r" file is stored
functions_dir<-"D:/Documents/CVPIA/LongtermStudyDesign/SanJoaquinSalmon_SampleSize"

# Source functions file.
source(file.path(functions_dir,"SampleSizeSimulations_functions.r"))



### Step 1. Define parameter values and derived parameters

# Parameters = sR, sR0, sA1, sA2, sB, sF, psiA1, psiA2, sA1.supp, pA1a, pA1b, pA2a, pA2b, pB1a, pB1b, pF1a, pF1b, pG2a, pG2b
#   where sA1 = survival from site A1 to A2 for fish from primary release (i.e., Durham Ferry release)
#     and sA1.supp = survival from supplemental release location to A2 for fish from supplemental release

# derived parameters = sR, sA, pA1, pA2, pB1, pF1, pG2

# Note: each parameter is assigned a single value.
pars_vec<-c()
pars_vec[c("sR","sR0","sA1","sA2","sB","sF")]<-c(0.01,0.2,0.1,0.1,0.1,0.05) # will overwrite sR below
pars_vec[c("psiA1","psiA2")]<-c(0.4,0.75)
pars_vec[c("sA1.supp")]<-pars_vec[c("sA1")]
pars_vec[c("pA1a","pA1b")]<-0.9
pars_vec[c("pA2a","pA2b")]<-0.85
pars_vec[c("pB1a","pB1b")]<-0.9
pars_vec[c("pF1a","pF1b")]<-0.85
pars_vec[c("pG2a","pG2b")]<-0.7


# define character vector of all parameters, including sR
all_par_names<-c("sR","sR0","sA1","sA2","sB","sF","psiA1","psiA2","sA1.supp","pA1a","pA1b","pA2a","pA2b","pB1a","pB1b","pF1a","pF1b","pG2a","pG2b")

# create data.frame of parameter values (nrow=1)
pars_mat<-data.frame(matrix(pars_vec,nrow=1,byrow=T))
names(pars_mat)<-all_par_names

# update value of sR: is derived from other parameters
pars_mat[,"sR"]<-apply(pars_mat,1,get.sR)

# define overall conditional detection probabilities at dual arrays
pars_mat$pA1 <- 1 - (1-pars_mat$pA1a)*(1-pars_mat$pA1b)
pars_mat$pA2 <- 1 - (1-pars_mat$pA2a)*(1-pars_mat$pA2b)
pars_mat$pB1 <- 1 - (1-pars_mat$pB1a)*(1-pars_mat$pB1b)
pars_mat$pF1 <- 1 - (1-pars_mat$pF1a)*(1-pars_mat$pF1b)
pars_mat$pG2 <- 1 - (1-pars_mat$pG2a)*(1-pars_mat$pG2b)



### Step 2. Define candidate release sizes and number of simulations
R_DF<-c(100,324,600)
R_supp<-c(100,324,600)
N_sim<-20  # 1000 or 5000 would be better but takes longer



### Step 3. Define matrix of release sizes
R_mat<-expand.grid(R_DF,R_supp)



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
dat<-list()
for(i in 1:nrow(R_mat))
{
  R_1<-R_mat[i,1]
  R_2<-R_mat[i,2]
  dat[[i]]<-list(R_DF=R_1,R_supp=R_2,mle_mat=t(sim.supp.seeds.pari.R1.R2(parvec<-pars_mat[,-1], R1=R_1, R2=R_2, nSim=N_sim)))
}

# dat is a length equal to the number of sample size combinations (=nrow(R_mat))
# - each element of dat is a list of length 3
# -- element 1 = R_DF = sample size of primary release group at Durham Ferry
# -- element 2 = R_supp = sample size of supplement release group
# -- element 3 = matrix of simulated MLEs for parameters (row 1 = true values); nrows = N_sim + 1; ncol = number of estimated parameters



### Step 8. Compute criteria performance summaries for criteria C1-C4 (optional)
 # Note: alternative to this step is Step 9.
C1_out_byR<-do.call(rbind,lapply(dat,calculate.Ci,"C1"))
C2_out_byR<-do.call(rbind,lapply(dat,calculate.Ci,"C2"))
C3_out_byR<-do.call(rbind,lapply(dat,calculate.Ci,"C3"))
C4_out_byR<-do.call(rbind,lapply(dat,calculate.Ci,"C4"))



### Step 9. Summarize results: caluculate minimum release size required to meet the criteria

# summarize minimum release size at Durham Ferry required to meet the 4 criteria
R_DF.min<-sapply(all_par_names,get.minR.1param,dat,R.type="DF")
row.names(R_DF.min)<-c("minR.overall","minR.c1","minR.c2","minR.c3","minR.c4",paste("N.ests_Rsupp",R_supp,sep="_"))
# the "N.ests_..." rows report the number of simulations that generated estimates for R_DF=minR.overall for the various levels of R_supp (sample size of supplemental release)

# summarize minimum supplemental release size required to meet the 4 criteria
R_supp.min<-sapply(all_par_names,get.minR.1param,dat,R.type="supp")
row.names(R_supp.min)<-c("minR.overall","minR.c1","minR.c2","minR.c3","minR.c4",paste("N.ests_RDF",R_DF,sep="_"))
# the "N.ests_..." rows report the number of simulations that generated estimates for R_supp=minR.overall for the various levels of R_DF (sample size of Durham Ferry release)


