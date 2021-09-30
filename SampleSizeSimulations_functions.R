############ Source file: functions

##### List of functions in this file
# sim.data.f         # simulate data using parameters named in par.names using values in par.values from multinomial likelihood with cell probabilities in like.vec and release size N.
# get.sA              # compute route-specific survival from Lathrop (site A1) to Chipps Island (site G2) via San Joaquin River route (= route A)
# get.sR              # compute total survival from release at Durham Ferry to Chipps Island (site G2)
# get.seeds.supp.wrap   # wrapper for get.seeds.supp
# get.seeds.supp        # compute seeds for parameters from supplemental release likelihood
# get.stat.supp         # get detailed summary statistics (a and b) from dat.vec (named, output = list)
# sim.supp.seeds.pari.R1.R2  # Simulate data, calculate summary statistics and seeds for each simulation, return matrix of seeds
# sim.supp.mle.pari.R1.R2    # Simulate data, calculate MLE of parameters using optim, return matrix of MLE estimates (NOT IMPLEMENTED)
# get.est.mat           # reshape output of apply(...,sim.supp.seeds.pari.R1.R2) 
# get.est.mat.wrap      # wrapper for get.est.mat so can use with sapply
# log.like.f            # returns the negative log-likelihood value for multinomial model specified in like.vec for data in dat.vec and values in pars
# get.C1                # compute proportion of parameter estimates that are neither NA nor infinite
# get.C2                # compute proportion of parameter estimates greater than specified cutoff value
# get.C3                # compute standard deviation of vector of parameter estimates; remove NA and infinite values
# get.C4                # compute absolute value of bias of vector of parameter estimates; remove NA and infinite values
# get.mean              # compute mean of vector of parameter estimates; remove NA and infinite values
# get.max               # compute maximum of vector of parameter estimates; remove NA and infinite values
# extract.astat         # extract "a" statistic at specified site from single simulation
# extract.bstat         # extract "b" statistic at specified site from single simulation
# calculate.Ci          # for given release size combination, evaluate summarization pertinent to specified criterion 
# get.minR.1param       # for param, gets minimum R for specified release type according to criteria C1, C2, and C3 and C4 from output of sim.supp.seeds.pari.R1.R2 (applied as list)


####################################################################################################################################
require("compiler")
enableJIT(3)


##### The functions:

sim.data.f<-function(N,par.values,par.names,var.vec,like.vec,nsim)
{
  # simulate data using parameters named in par.names using values in par.values from multinomial likelihood with cell probabilities in like.vec and release size N.
  # var.vec = named character vector of "variables" (as in USER model)
  # nsim = number of data sets to simulate
  # N = number of trials (release size)
  
  # return named vector of cell counts.
  
  for(i in 1:length(par.names)) assign(par.names[i],par.values[i])
  f1<-function(x) eval(parse(text=x))
  for(i in 1:length(var.vec)) assign(names(var.vec)[i],f1(var.vec[i]))
  p.vec<-sapply(like.vec,f1)
  sim.d<-rmultinom(n=nsim, size=N, prob=p.vec) # column = cell counts for single simulation; row = different simulations
  return(sim.d)
}


get.sA<-function(pars)
{
  # compute route-specific survival from Lathrop (site A1) to Chipps Island (site G2) via San Joaquin River route (= route A)

  if(any(is.na(pars[c("psiA2","sA1","sA2")])) | any(!is.finite(pars[c("psiA2","sA1","sA2")]))) return(NA)
  if(pars["psiA2"]!=1) return(pars["sA1"]*(pars["sA2"]*pars["psiA2"] + pars["sF"]*(1-pars["psiA2"])))
  return(prod(pars[c("sA1","sA2")]))
}

get.sR<-function(pars)
{
  # compute total survival from release at Durham Ferry to Chipps Island (site G2)

  sA<-get.sA(pars)
  if(is.na(sA) | !is.finite(sA)) return(NA)
  if(any(is.na(pars[c("psiA1","sR0")])) | any(!is.finite(pars[c("psiA1","sR0")]))) return(NA)
  if(pars["psiA1"]!=1) return(pars["sR0"]*(pars["psiA1"]*sA + (1-pars["psiA1"])*pars["sB"]))
  return(pars["sR0"]*sA)
}

get.seeds.supp.wrap<-function(sim,DF.list,supp.like,A1.mat,A2.mat,B1.mat,F1.mat,G2.mat,p.full=unique(c(pars.full.DF,pars.full.supp,pars.dual.A1,pars.dual.A2,pars.dual.B1,pars.dual.F1,pars.dual.G2)))
{
  # wrapper for get.seeds.supp
  d.stats<-DF.list[[sim]]
  s.stats<-supp.like[[sim]]
  a1.dat<-A1.mat[,sim]
  a2.dat<-A2.mat[,sim]
  b1.dat<-B1.mat[,sim]
  f1.dat<-F1.mat[,sim]
  g2.dat<-G2.mat[,sim]
  return(get.seeds.supp(d.stats,s.stats,a1.dat,a2.dat,b1.dat,f1.dat,g2.dat,pars.full=p.full))
}


get.seeds.supp<-function(DF.stats,supp.stats,A1.dat,A2.dat,B1.dat,F1.dat,G2.dat,pars.full=unique(c(pars.full.DF,pars.full.supp,pars.dual.A1,pars.dual.A2,pars.dual.B1,pars.dual.F1,pars.dual.G2)))
{
  # compute seeds for parameters from supplemental release likelihood
  # DF.stats and supp.stats are lists of length 2 with names "a" and "b"
  # A1.dat, etc., are vectors of data of length 3, form = "AB","A0","B0"
  
  seeds<-rep(NA,length(pars.full)); names(seeds)<-pars.full
  pars0<-c(); pars1<-c(); parsNA<-c()
  
  # actually do not need pA1, etc. for dual arrays
  seeds<-seeds[-match(c("pA1","pA2","pB1","pF1","pG2"),names(seeds))]
  
  seeds["pA1a"]<-A1.dat[grep("AB",names(A1.dat))]/sum(A1.dat[grep("B",names(A1.dat))])
  seeds["pA1b"]<-A1.dat[grep("AB",names(A1.dat))]/sum(A1.dat[grep("A",names(A1.dat))])
  
  seeds["pA2a"]<-A2.dat[grep("AB",names(A2.dat))]/sum(A2.dat[grep("B",names(A2.dat))])
  seeds["pA2b"]<-A2.dat[grep("AB",names(A2.dat))]/sum(A2.dat[grep("A",names(A2.dat))])
  
  seeds["pB1a"]<-B1.dat[grep("AB",names(B1.dat))]/sum(B1.dat[grep("B",names(B1.dat))])
  seeds["pB1b"]<-B1.dat[grep("AB",names(B1.dat))]/sum(B1.dat[grep("A",names(B1.dat))])
  
  seeds["pF1a"]<-F1.dat[grep("AB",names(F1.dat))]/sum(F1.dat[grep("B",names(F1.dat))])
  seeds["pF1b"]<-F1.dat[grep("AB",names(F1.dat))]/sum(F1.dat[grep("A",names(F1.dat))])
  
  seeds["pG2a"]<-G2.dat[grep("AB",names(G2.dat))]/sum(G2.dat[grep("B",names(G2.dat))])
  seeds["pG2b"]<-G2.dat[grep("AB",names(G2.dat))]/sum(G2.dat[grep("A",names(G2.dat))])
  
  pA1<-1-(1-seeds["pA1a"])*(1-seeds["pA1b"])
  pA2<-1-(1-seeds["pA2a"])*(1-seeds["pA2b"])
  pB1<-1-(1-seeds["pB1a"])*(1-seeds["pB1b"])
  pF1<-1-(1-seeds["pF1a"])*(1-seeds["pF1b"])
  pG2<-1-(1-seeds["pG2a"])*(1-seeds["pG2b"])
  
  seeds["sR0"]<-(1/DF.stats[["a"]]["N"])*(DF.stats[["a"]]["A1"]/pA1 + DF.stats[["a"]]["B1"]/pB1)
  seeds["psiA1"]<-(DF.stats[["a"]]["A1"]/pA1)/(DF.stats[["a"]]["N"]*seeds["sR0"])
  seeds["sB"]<-(DF.stats[["b"]]["B1"]/DF.stats[["a"]]["B1"])*(1/pG2)
  
  seeds["sA1"]<-(pA1/DF.stats[["a"]]["A1"])*(DF.stats[["a"]]["A2"]/pA2 + DF.stats[["a"]]["F1"]/pF1)
  seeds["sA1.supp"]<-(1/supp.stats[["a"]]["N"])*(supp.stats[["a"]]["A2"]/pA2 + supp.stats[["a"]]["F1"]/pF1)
  seeds["psiA2"]<-(DF.stats[["a"]]["A2"] + supp.stats[["a"]]["A2"])/((DF.stats[["a"]]["A1"]*seeds["sA1"] + supp.stats[["a"]]["N"]*seeds["sA1.supp"])*pA2)
  seeds["sA2"]<-(DF.stats[["b"]]["A2"] + supp.stats[["b"]]["A2"])/((DF.stats[["a"]]["A2"] + supp.stats[["a"]]["A2"])*pG2)
  seeds["sF"]<-(DF.stats[["b"]]["F1"] + supp.stats[["b"]]["F1"])/((DF.stats[["a"]]["F1"] + supp.stats[["a"]]["F1"])*pG2)
  
  parsNA<-unique(c(parsNA,names(seeds[is.na(seeds) & !is.finite(seeds)])))
  pars1<-unique(c(pars1,names(seeds[!is.na(seeds) & is.finite(seeds) & seeds==1])))
  pars0<-unique(c(pars0,names(seeds[!is.na(seeds) & is.finite(seeds) & seeds==0])))
  
  return(seeds)   # so not actually returning parsNA, pars1, pars0 - they are all in seeds
}



get.stat.supp<-function(dat.vec,release.site)
{
  # get detailed summary statistics (a and b) from dat.vec (named, output = list)
  # release.site = "DF" or "supp"
  
  a<-c()
  b<-c()
  x<-list()  # indices of a and b
  
  x[["A2"]]<-grep("A2",names(dat.vec)); x[["F1"]]<-grep("F1",names(dat.vec)); x[["G2"]]<-grep("G2",names(dat.vec))
  
  if(release.site=="DF")
  {
    x[["A1"]]<-grep("A1",names(dat.vec)); x[["B1"]]<-grep("B1",names(dat.vec))
    a["A1"]<-sum(dat.vec[x[["A1"]]]); a["B1"]<-sum(dat.vec[x[["B1"]]])
    b["A1"]<-sum(dat.vec[x[["A1"]][x[["A1"]] %in% c(x[["A2"]],x[["F1"]],x[["G2"]])]])
    b["B1"]<-sum(dat.vec[x[["B1"]][x[["B1"]] %in% x[["G2"]]]])
  }
  
  a["N"]<-sum(dat.vec)
  a["A2"]<-sum(dat.vec[x[["A2"]]]); a["F1"]<-sum(dat.vec[x[["F1"]]]); a["G2"]<-sum(dat.vec[x[["G2"]]])
  b["A2"]<-sum(dat.vec[x[["A2"]][x[["A2"]] %in% x[["G2"]]]])
  b["F1"]<-sum(dat.vec[x[["F1"]][x[["F1"]] %in% x[["G2"]]]])
  
  out<-list(a,b)
  names(out)<-c("a","b")
  return(out)
}

sim.supp.seeds.pari.R1.R2<-function(parvec,R1,R2,nSim,DF.pars,supp.pars,A1.pars,A2.pars,B1.pars,F1.pars,G2.pars,DF.vars,supp.vars,A1.vars,A2.vars,B1.vars,F1.vars,G2.vars,DF.like,supp.like,A1.like,A2.like,B1.like,F1.like,G2.like)
{
  # Simulate data, calculate summary statistics and seeds for each simulation, return matrix of seeds
  
  # For release size R1 at Durham Ferry and release size R2 downstream of A1, and single set of parameter values in parvec (names in xx.pars)
  #   Simulate data using release-recapture model defined in xx.like and xx.vars
  #   Calculate summary statistics for each simulation
  #   Calculate seeds - for each simulation
  #   Store results: 1 row per data set, just the estimates and release size
  # Return matrix of results

  
    
  ########################################################################################
  
  #################### Arguments
  
  # parvec = named numeric vector of parameter values
  
  # R1 = sample size of Durham Ferry release group (primary release)
  # R2 = sample size of supplemental release group
  # nSim = number of simulations to run
  
  ###### Names of model parameters
  ## Basic model parameters (these will be estimated directly using maximum likelihood)
  
  # DF.pars = parameters used in Durham Ferry likelihood: default value = c("sR0","sA1","sA2","sB","sF","psiA1","psiA2","pA1","pA2","pB1","pF1","pG2")
  # supp.pars = parameters used in supplemental release likelihood: default value = c("sA1.supp","sA2","sF","psiA2","pA2","pF1","pG2")
  # A1.pars = parameters used in auxiliary likelihood for dual array at A1: default value = c("pA1a","pA1b")
  # A2.pars = parameters used in auxiliary likelihood for dual array at A2: default value = c("pA2a","pA2b")
  # B1.pars = parameters used in auxiliary likelihood for dual array at B1: default value = c("pB1a","pB1b")
  # F1.pars = parameters used in auxiliary likelihood for dual array at F1: default value = c("pF1a","pF1b")
  # G2.pars = parameters used in auxiliary likelihood for dual array at G2: default value = c("pG2a","pG2b")
  
  
  ## Derived parameters that are either used in the likelihood expression or are estimation targets
  # - each of these items is a named character vector containing the analytic expression that defines the relevant derived parameters
  # - default values are provided 
  
  # DF.vars = analytic definitions of parameters used in likelihood for Durham Ferry release; default = definitions of qA1, qA2, qB1, qF1, qG2, psiB1, psiF2, xF1, xA2, xB1, xA1, xR, gammaA1G2, gammaRG2
  # supp.vars = analytic definitions of parameters used in likelihood for supplemental release; default = definitions of qA2, qF1, qG2, psiF2, xF1, xA2, xSupp
  # A1.vars = analytic definition of pA1 in terms of pA1a, pA1b
  # A2.vars = analytic definition of pA2 in terms of pA2a, pA2b
  # B1.vars = analytic definition of pB1 in terms of pB1a, pB1b
  # F1.vars = analytic definition of pF1 in terms of pF1a, pF1b
  # G2.vars = analytic definition of pG2 in terms of pG2a, pG2b
  

  
  ##### Analytic expressions of parameterizations of likelihood components
  # - each of these items is a named character vector containing the analytic expression that defines the parameterizations associated with detection histories for given likelihood
  # - default values are provided

  # DF.like = analytic parameterizations that define likelihood for Durham Ferry release; default = 
   # "DF A1 A2 G2", "DF A1 A2 x", "DF A1 F1 G2", "DF A1 F1 x", "DF A1 G2", "DF A1 x", 
   # "DF B1 G2", "DF B1 x", "DF A2 G2", "DF A2 x", "DF F1 G2", "DF F1 x", "DF G2", "DF x"
  # supp.like = analytic parameterizations that define likelihood for supplemental release; default = "supp A2 G2", "supp A2 x", "supp F2 G2", "supp F2 x", "supp G2", "supp x"
  # A1.like = analytic parameterizations that define auxiliary likelihood for A1 dual array
  # A2.like = analytic parameterizations that define auxiliary likelihood for A2 dual array
  # B1.like = analytic parameterizations that define auxiliary likelihood for B1 dual array
  # F1.like = analytic parameterizations that define auxiliary likelihood for F1 dual array
  # G2.like = analytic parameterizations that define auxiliary likelihood for G2 dual array
  
  ########################################################################################
  
  parvec<-unlist(parvec)
  
  # compile vector of likelihood detection history parameterizations
  if(missing(DF.like))
  {
    DF.like<-c()
    DF.like["DF A1 A2 G2"] <- "sR0*psiA1*pA1*sA1*psiA2*pA2*sA2*pG2"
    DF.like["DF A1 A2 x"] <- "sR0*psiA1*pA1*sA1*psiA2*pA2*xA2"
    DF.like["DF A1 F1 G2"] <- "sR0*psiA1*pA1*sA1*psiF2*pF1*sF*pG2"
    DF.like["DF A1 F1 x"] <- "sR0*psiA1*pA1*sA1*psiF2*pF1*xF1"
    DF.like["DF A1 G2"] <- "sR0*psiA1*pA1*gammaA1G2*pG2"
    DF.like["DF A1 x"] <- "sR0*psiA1*pA1*xA1"
    DF.like["DF B1 G2"] <- "sR0*psiB1*pB1*sB*pG2"
    DF.like["DF B1 x"] <- "sR0*psiB1*pB1*xB1"
    DF.like["DF A2 G2"] <- "sR0*psiA1*qA1*sA1*psiA2*pA2*sA2*pG2"
    DF.like["DF A2 x"] <- "sR0*psiA1*qA1*sA1*psiA2*pA2*xA2"
    DF.like["DF F1 G2"] <- "sR0*psiA1*qA1*sA1*psiF2*pF1*sF*pG2"
    DF.like["DF F1 x"] <- "sR0*psiA1*qA1*sA1*psiF2*pF1*xF1"
    DF.like["DF G2"] <- "gammaRG2*pG2"
    DF.like["DF x"] <- "xR"
  }
  
  if(missing(supp.like))
  {
    supp.like<-c()
    supp.like["supp A2 G2"] <- "sA1.supp*psiA2*pA2*sA2*pG2"
    supp.like["supp A2 x"] <- "sA1.supp*psiA2*pA2*xA2"
    supp.like["supp F1 G2"] <- "sA1.supp*psiF2*pF1*sF*pG2"
    supp.like["supp F1 x"] <- "sA1.supp*psiF2*pF1*xF1"
    supp.like["supp G2"] <- "sA1.supp*(psiA2*qA2*sA2 + psiF2*qF1*sF)*pG2"
    supp.like["supp x"] <- "xSupp"    
  }

  if(missing(A1.like)) 
  {
    A1.like<-c()
    A1.like["A1 AB"] <- "pA1a*pA1b/pA1"
    A1.like["A1 A0"] <- "pA1a*(1-pA1b)/pA1"
    A1.like["A1 B0"] <- "(1-pA1a)*pA1b/pA1"
  }
  
  if(missing(A2.like)) 
  {
    A2.like<-c()
    A2.like["A2 AB"] <- "pA2a*pA2b/pA2"
    A2.like["A2 A0"] <- "pA2a*(1-pA2b)/pA2"
    A2.like["A2 B0"] <- "(1-pA2a)*pA2b/pA2"
  }
  
  if(missing(F1.like)) 
  {
    F1.like<-c()
    F1.like["F1 AB"] <- "pF1a*pF1b/pF1"
    F1.like["F1 A0"] <- "pF1a*(1-pF1b)/pF1"
    F1.like["F1 B0"] <- "(1-pF1a)*pF1b/pF1"
  }
  
  if(missing(B1.like)) 
  {
    B1.like<-c()
    B1.like["B1 AB"] <- "pB1a*pB1b/pB1"
    B1.like["B1 A0"] <- "pB1a*(1-pB1b)/pB1"
    B1.like["B1 B0"] <- "(1-pB1a)*pB1b/pB1"
  }
  
  if(missing(G2.like)) 
  {
    G2.like<-c()
    G2.like["G2 AB"] <- "pG2a*pG2b/pG2"
    G2.like["G2 A0"] <- "pG2a*(1-pG2b)/pG2"
    G2.like["G2 B0"] <- "(1-pG2a)*pG2b/pG2"
  }
  
  like.all<-c(DF.like,supp.like,A1.like,A2.like,B1.like,F1.like,G2.like)
  

  # compile vector of definitions of derived parameters
  if(missing(DF.vars))
  {
    DF.vars<-c()
    DF.vars["qA1"] <- "1 - pA1"
    DF.vars["qA2"] <- "1 - pA2"
    DF.vars["qB1"] <- "1 - pB1"
    DF.vars["qF1"] <- "1 - pF1"
    DF.vars["qG2"] <- "1 - pG2"
    DF.vars["psiB1"] <- "1 - psiA1"
    DF.vars["psiF2"] <- "1 - psiA2"
    DF.vars["xF1"] <- "1 - sF + sF*qG2"
    DF.vars["xA2"] <- "1 - sA2 + sA2*qG2"
    DF.vars["xB1"] <- "1 - sB + sB*qG2"
    DF.vars["xA1"] <- "1 - sA1 + sA1*(psiA2*qA2*xA2 + psiF2*qF1*xF1)"
    DF.vars["xR"] <- "1 - sR0 + sR0*(psiA1*qA1*xA1 + psiB1*qB1*xB1)"
    DF.vars["gammaA1G2"] <- "sA1*(psiA2*qA2*sA2 + psiF2*qF1*sF)"
    DF.vars["gammaRG2"] <- "sR0*(psiA1*qA1*gammaA1G2 + psiB1*qB1*sB)"
  }
  
  if(missing(supp.vars))
  {
    supp.vars<-c()
    supp.vars["qA2"] <- "1 - pA2"
    supp.vars["qF1"] <- "1 - pF1"
    supp.vars["qG2"] <- "1 - pG2"
    supp.vars["psiF2"] <- "1 - psiA2"
    supp.vars["xF1"] <- "1 - sF + sF*qG2"
    supp.vars["xA2"] <- "1 - sA2 + sA2*qG2"
    supp.vars["xSupp"] <- "1 - sA1.supp + sA1.supp*(psiA2*qA2*xA2 + psiF2*qF1*xF1)"
  }
  
  if(missing(A1.vars)) {A1.vars<-c(); A1.vars["pA1"] <- "1 - (1 - pA1a)*(1 - pA1b)"}
  if(missing(A2.vars)) {A2.vars<-c(); A2.vars["pA2"] <- "1 - (1 - pA2a)*(1 - pA2b)"}
  if(missing(B1.vars)) {B1.vars<-c(); B1.vars["pB1"] <- "1 - (1 - pB1a)*(1 - pB1b)"}
  if(missing(F1.vars)) {F1.vars<-c(); F1.vars["pF1"] <- "1 - (1 - pF1a)*(1 - pF1b)"}
  if(missing(G2.vars)) {G2.vars<-c(); G2.vars["pG2"] <- "1 - (1 - pG2a)*(1 - pG2b)"}
  
  vars.all<-c(A1.vars,A2.vars,B1.vars,F1.vars,G2.vars,DF.vars,supp.vars)
  vars.all<-vars.all[!duplicated(names(vars.all))]
  
  # compile vector of parameter names from primary and supplemental releases
  {
    if(missing(DF.pars)) DF.pars<-c("sR0","sA1","sA2","sB","sF","psiA1","psiA2","pA1","pA2","pB1","pF1","pG2")
    if(missing(supp.pars)) supp.pars<-c("sA1.supp","sA2","sF","psiA2","pA2","pF1","pG2")
    if(missing(A1.pars)) A1.pars<-c("pA1a","pA1b")
    if(missing(A2.pars)) A2.pars<-c("pA2a","pA2b")
    if(missing(B1.pars)) B1.pars<-c("pB1a","pB1b")
    if(missing(F1.pars)) F1.pars<-c("pF1a","pF1b")
    if(missing(G2.pars)) G2.pars<-c("pG2a","pG2b")    
  }

  pars.supp.all<-unique(c(DF.pars,supp.pars,A1.pars,A2.pars,B1.pars,F1.pars,G2.pars))
  pars.supp.fit<-pars.supp.all[-match(c("pA1","pA2","pB1","pF1","pG2"),pars.supp.all)]
  
  # set up matrix with results (col = simulation, row = parameter (first col = parameter values))
  mat.out<-matrix(NA,nrow=length(pars.supp.fit)+2,ncol=nSim+1)
  row.names(mat.out)<-c(pars.supp.fit,"sA","sR")
  mat.out[,1]<-unlist(c(parvec[pars.supp.fit],get.sA(parvec),get.sR(parvec)))
  
  # simulate release-recapture data from Durham Ferry release
  data.DF<-sim.data.f(par.values=parvec[DF.pars], par.names=DF.pars, var.vec=DF.vars, like.vec=DF.like, N=R1, nsim=nSim)
  
  # simulate release-recapture data from supplemental release
  data.supp<-sim.data.f(par.values=parvec[supp.pars], par.names=supp.pars, var.vec=supp.vars, like.vec=supp.like, N=R2, nsim=nSim)
  
  # compute sufficient statistics for Durham Ferry release
  stats.DF<-apply(data.DF,2,get.stat.supp,"DF")
  
  # compute sufficient statistics for supplemental release
  stats.supp<-apply(data.supp,2,get.stat.supp,"supp")
  
  # simulate detections at the two lines of each dual array, conditional on total number of tags detected at each dual array
  data.A1<-sapply(unlist(lapply(stats.DF,extract.astat,"A1")),sim.data.f,par.values=parvec[A2.pars], par.names=A1.pars, var.vec=A1.vars, like.vec=A1.like, nsim=1)
  data.A2<-sapply(unlist(lapply(stats.DF,extract.astat,"A2"))+unlist(lapply(stats.supp,extract.astat,"A2")),sim.data.f,par.values=parvec[A2.pars], par.names=A2.pars, var.vec=A2.vars, like.vec=A2.like, nsim=1)
  data.B1<-sapply(unlist(lapply(stats.DF,extract.astat,"B1")),sim.data.f,par.values=parvec[B1.pars], par.names=B1.pars, var.vec=B1.vars, like.vec=B1.like, nsim=1)
  data.F1<-sapply(unlist(lapply(stats.DF,extract.astat,"F1"))+unlist(lapply(stats.supp,extract.astat,"F1")),sim.data.f,par.values=parvec[F1.pars], par.names=F1.pars, var.vec=F1.vars, like.vec=F1.like, nsim=1)
  data.G2<-sapply(unlist(lapply(stats.DF,extract.astat,"G2"))+unlist(lapply(stats.supp,extract.astat,"G2")),sim.data.f,par.values=parvec[G2.pars], par.names=G2.pars, var.vec=G2.vars, like.vec=G2.like, nsim=1)
  
  row.names(data.A1)<-paste("A1",c("AB","A0","B0"),sep=" ")
  row.names(data.A2)<-paste("A2",c("AB","A0","B0"),sep=" ")
  row.names(data.B1)<-paste("B1",c("AB","A0","B0"),sep=" ")
  row.names(data.F1)<-paste("F1",c("AB","A0","B0"),sep=" ")
  row.names(data.G2)<-paste("G2",c("AB","A0","B0"),sep=" ")
  
  # combine all data
  data.sim<-rbind(data.DF,data.supp,data.A1,data.A2,data.B1,data.F1,data.G2)
  
  # get seeds, store, and return (haven't dealt with reducing model if necessary)
  seeds<-sapply(1:nSim,get.seeds.supp.wrap,stats.DF,stats.supp,data.A1,data.A2,data.B1,data.F1,data.G2,p.full=pars.supp.all)
  seeds<-rbind(seeds,matrix(NA,nrow=2,ncol=nSim))
  row.names(seeds)[nrow(seeds)-c(1,0)]<-c("sA","sR")
  seeds["sA",]<-apply(seeds[1:18,],2,get.sA)
  seeds["sR",]<-apply(seeds[1:18,],2,get.sR)
  mat.out[,-1]<-seeds
  return(mat.out)
}




sim.supp.mle.pari.R1.R2<-function(parvec,R1,R2,nSim,DF.pars,supp.pars,A1.pars,A2.pars,B1.pars,F1.pars,G2.pars,DF.vars,supp.vars,A1.vars,A2.vars,B1.vars,F1.vars,G2.vars,DF.like,supp.like,A1.like,A2.like,B1.like,F1.like,G2.like)
{
  # Simulate data, calculate MLE of parameters using optim, return matrix of MLE estimates (NOT IMPLEMENTED)

  # For release size R1 at Durham Ferry and release size R2 downstream of A1, and single set of parameter values in parvec (names in xx.pars)
  #   Simulate data using release-recapture model defined in xx.like and xx.vars
  #   Calculate summary statistics for each simulation
  #   Calculate seeds - for each simulation
  #   Estimate parameters using optim (MLEs) - without hessian
  #   Store results: 1 row per data set, just the estimates and release size
  # Return matrix of results
  
  ########################################################################################
  
  #################### Arguments
  
  # parvec = named numeric vector of parameter values
  
  # R1 = sample size of Durham Ferry release group (primary release)
  # R2 = sample size of supplemental release group
  # nSim = number of simulations to run
  
  ###### Names of model parameters
  ## Basic model parameters (these will be estimated directly using maximum likelihood)
  
  # DF.pars = parameters used in Durham Ferry likelihood: default value = c("sR0","sA1","sA2","sB","sF","psiA1","psiA2","pA1","pA2","pB1","pF1","pG2")
  # supp.pars = parameters used in supplemental release likelihood: default value = c("sA1.supp","sA2","sF","psiA2","pA2","pF1","pG2")
  # A1.pars = parameters used in auxiliary likelihood for dual array at A1: default value = c("pA1a","pA1b")
  # A2.pars = parameters used in auxiliary likelihood for dual array at A2: default value = c("pA2a","pA2b")
  # B1.pars = parameters used in auxiliary likelihood for dual array at B1: default value = c("pB1a","pB1b")
  # F1.pars = parameters used in auxiliary likelihood for dual array at F1: default value = c("pF1a","pF1b")
  # G2.pars = parameters used in auxiliary likelihood for dual array at G2: default value = c("pG2a","pG2b")
  
  
  ## Derived parameters that are either used in the likelihood expression or are estimation targets
  # - each of these items is a named character vector containing the analytic expression that defines the relevant derived parameters
  # - default values are provided 
  
  # DF.vars = analytic definitions of parameters used in likelihood for Durham Ferry release; default = definitions of qA1, qA2, qB1, qF1, qG2, psiB1, psiF2, xF1, xA2, xB1, xA1, xR, gammaA1G2, gammaRG2
  # supp.vars = analytic definitions of parameters used in likelihood for supplemental release; default = definitions of qA2, qF1, qG2, psiF2, xF1, xA2, xSupp
  # A1.vars = analytic definition of pA1 in terms of pA1a, pA1b
  # A2.vars = analytic definition of pA2 in terms of pA2a, pA2b
  # B1.vars = analytic definition of pB1 in terms of pB1a, pB1b
  # F1.vars = analytic definition of pF1 in terms of pF1a, pF1b
  # G2.vars = analytic definition of pG2 in terms of pG2a, pG2b
  
  
  
  ##### Analytic expressions of parameterizations of likelihood components
  # - each of these items is a named character vector containing the analytic expression that defines the parameterizations associated with detection histories for given likelihood
  # - default values are provided
  
  # DF.like = analytic parameterizations that define likelihood for Durham Ferry release; default = 
  # "DF A1 A2 G2", "DF A1 A2 x", "DF A1 F1 G2", "DF A1 F1 x", "DF A1 G2", "DF A1 x", 
  # "DF B1 G2", "DF B1 x", "DF A2 G2", "DF A2 x", "DF F1 G2", "DF F1 x", "DF G2", "DF x"
  # supp.like = analytic parameterizations that define likelihood for supplemental release; default = "supp A2 G2", "supp A2 x", "supp F2 G2", "supp F2 x", "supp G2", "supp x"
  # A1.like = analytic parameterizations that define auxiliary likelihood for A1 dual array
  # A2.like = analytic parameterizations that define auxiliary likelihood for A2 dual array
  # B1.like = analytic parameterizations that define auxiliary likelihood for B1 dual array
  # F1.like = analytic parameterizations that define auxiliary likelihood for F1 dual array
  # G2.like = analytic parameterizations that define auxiliary likelihood for G2 dual array
  
  ########################################################################################
  
  parvec<-unlist(parvec)

  
  # compile vector of likelihood detection history parameterizations
  if(missing(DF.like))
  {
    DF.like<-c()
    DF.like["DF A1 A2 G2"] <- "sR0*psiA1*pA1*sA1*psiA2*pA2*sA2*pG2"
    DF.like["DF A1 A2 x"] <- "sR0*psiA1*pA1*sA1*psiA2*pA2*xA2"
    DF.like["DF A1 F1 G2"] <- "sR0*psiA1*pA1*sA1*psiF2*pF1*sF*pG2"
    DF.like["DF A1 F1 x"] <- "sR0*psiA1*pA1*sA1*psiF2*pF1*xF1"
    DF.like["DF A1 G2"] <- "sR0*psiA1*pA1*gammaA1G2*pG2"
    DF.like["DF A1 x"] <- "sR0*psiA1*pA1*xA1"
    DF.like["DF B1 G2"] <- "sR0*psiB1*pB1*sB*pG2"
    DF.like["DF B1 x"] <- "sR0*psiB1*pB1*xB1"
    DF.like["DF A2 G2"] <- "sR0*psiA1*qA1*sA1*psiA2*pA2*sA2*pG2"
    DF.like["DF A2 x"] <- "sR0*psiA1*qA1*sA1*psiA2*pA2*xA2"
    DF.like["DF F1 G2"] <- "sR0*psiA1*qA1*sA1*psiF2*pF1*sF*pG2"
    DF.like["DF F1 x"] <- "sR0*psiA1*qA1*sA1*psiF2*pF1*xF1"
    DF.like["DF G2"] <- "gammaRG2*pG2"
    DF.like["DF x"] <- "xR"
  }
  
  if(missing(supp.like))
  {
    supp.like<-c()
    supp.like["supp A2 G2"] <- "sA1.supp*psiA2*pA2*sA2*pG2"
    supp.like["supp A2 x"] <- "sA1.supp*psiA2*pA2*xA2"
    supp.like["supp F1 G2"] <- "sA1.supp*psiF2*pF1*sF*pG2"
    supp.like["supp F1 x"] <- "sA1.supp*psiF2*pF1*xF1"
    supp.like["supp G2"] <- "sA1.supp*(psiA2*qA2*sA2 + psiF2*qF1*sF)*pG2"
    supp.like["supp x"] <- "xSupp"    
  }
  
  if(missing(A1.like)) 
  {
    A1.like<-c()
    A1.like["A1 AB"] <- "pA1a*pA1b/pA1"
    A1.like["A1 A0"] <- "pA1a*(1-pA1b)/pA1"
    A1.like["A1 B0"] <- "(1-pA1a)*pA1b/pA1"
  }
  
  if(missing(A2.like)) 
  {
    A2.like<-c()
    A2.like["A2 AB"] <- "pA2a*pA2b/pA2"
    A2.like["A2 A0"] <- "pA2a*(1-pA2b)/pA2"
    A2.like["A2 B0"] <- "(1-pA2a)*pA2b/pA2"
  }
  
  if(missing(F1.like)) 
  {
    F1.like<-c()
    F1.like["F1 AB"] <- "pF1a*pF1b/pF1"
    F1.like["F1 A0"] <- "pF1a*(1-pF1b)/pF1"
    F1.like["F1 B0"] <- "(1-pF1a)*pF1b/pF1"
  }
  
  if(missing(B1.like)) 
  {
    B1.like<-c()
    B1.like["B1 AB"] <- "pB1a*pB1b/pB1"
    B1.like["B1 A0"] <- "pB1a*(1-pB1b)/pB1"
    B1.like["B1 B0"] <- "(1-pB1a)*pB1b/pB1"
  }
  
  if(missing(G2.like)) 
  {
    G2.like<-c()
    G2.like["G2 AB"] <- "pG2a*pG2b/pG2"
    G2.like["G2 A0"] <- "pG2a*(1-pG2b)/pG2"
    G2.like["G2 B0"] <- "(1-pG2a)*pG2b/pG2"
  }
  
  like.all<-c(DF.like,supp.like,A1.like,A2.like,B1.like,F1.like,G2.like)
  
  
  # compile vector of definitions of derived parameters
  if(missing(DF.vars))
  {
    DF.vars<-c()
    DF.vars["qA1"] <- "1 - pA1"
    DF.vars["qA2"] <- "1 - pA2"
    DF.vars["qB1"] <- "1 - pB1"
    DF.vars["qF1"] <- "1 - pF1"
    DF.vars["qG2"] <- "1 - pG2"
    DF.vars["psiB1"] <- "1 - psiA1"
    DF.vars["psiF2"] <- "1 - psiA2"
    DF.vars["xF1"] <- "1 - sF + sF*qG2"
    DF.vars["xA2"] <- "1 - sA2 + sA2*qG2"
    DF.vars["xB1"] <- "1 - sB + sB*qG2"
    DF.vars["xA1"] <- "1 - sA1 + sA1*(psiA2*qA2*xA2 + psiF2*qF1*xF1)"
    DF.vars["xR"] <- "1 - sR0 + sR0*(psiA1*qA1*xA1 + psiB1*qB1*xB1)"
    DF.vars["gammaA1G2"] <- "sA1*(psiA2*qA2*sA2 + psiF2*qF1*sF)"
    DF.vars["gammaRG2"] <- "sR0*(psiA1*qA1*gammaA1G2 + psiB1*qB1*sB)"
  }
  
  if(missing(supp.vars))
  {
    supp.vars<-c()
    supp.vars["qA2"] <- "1 - pA2"
    supp.vars["qF1"] <- "1 - pF1"
    supp.vars["qG2"] <- "1 - pG2"
    supp.vars["psiF2"] <- "1 - psiA2"
    supp.vars["xF1"] <- "1 - sF + sF*qG2"
    supp.vars["xA2"] <- "1 - sA2 + sA2*qG2"
    supp.vars["xSupp"] <- "1 - sA1.supp + sA1.supp*(psiA2*qA2*xA2 + psiF2*qF1*xF1)"
  }
  
  if(missing(A1.vars)) {A1.vars<-c(); A1.vars["pA1"] <- "1 - (1 - pA1a)*(1 - pA1b)"}
  if(missing(A2.vars)) {A2.vars<-c(); A2.vars["pA2"] <- "1 - (1 - pA2a)*(1 - pA2b)"}
  if(missing(B1.vars)) {B1.vars<-c(); B1.vars["pB1"] <- "1 - (1 - pB1a)*(1 - pB1b)"}
  if(missing(F1.vars)) {F1.vars<-c(); F1.vars["pF1"] <- "1 - (1 - pF1a)*(1 - pF1b)"}
  if(missing(G2.vars)) {G2.vars<-c(); G2.vars["pG2"] <- "1 - (1 - pG2a)*(1 - pG2b)"}
  
  vars.all<-c(A1.vars,A2.vars,B1.vars,F1.vars,G2.vars,DF.vars,supp.vars)
  vars.all<-vars.all[!duplicated(names(vars.all))]
  
  # compile vector of parameter names from primary and supplemental releases
  if(missing(DF.pars)) DF.pars<-c("sR0","sA1","sA2","sB","sF","psiA1","psiA2","pA1","pA2","pB1","pF1","pG2")
  if(missing(supp.pars)) supp.pars<-c("sA1.supp","sA2","sF","psiA2","pA2","pF1","pG2")
  if(missing(A1.pars)) A1.pars<-c("pA1a","pA1b")
  if(missing(A2.pars)) A2.pars<-c("pA2a","pA2b")
  if(missing(B1.pars)) B1.pars<-c("pB1a","pB1b")
  if(missing(F1.pars)) F1.pars<-c("pF1a","pF1b")
  if(missing(G2.pars)) G2.pars<-c("pG2a","pG2b")
  
  pars.supp.all<-unique(c(DF.pars,supp.pars,A1.pars,A2.pars,B1.pars,F1.pars,G2.pars))
  pars.supp.fit<-pars.supp.all[-match(c("pA1","pA2","pB1","pF1","pG2"),pars.supp.all)]
  
  # set up matrix with results (col = simulation, row = parameter (first col = parameter values))
  mat.out<-matrix(NA,nrow=length(pars.supp.fit)+2,ncol=nSim+1)
  row.names(mat.out)<-c(pars.supp.fit,"sA","sR")
  mat.out[,1]<-unlist(c(parvec[pars.supp.fit],get.sA(parvec),get.sR(parvec)))
  
  # simulate release-recapture data from Durham Ferry release
  data.DF<-sim.data.f(par.values=parvec[DF.pars], par.names=DF.pars, var.vec=DF.vars, like.vec=DF.like, N=R1, nsim=nSim)
  
  # simulate release-recapture data from supplemental release
  data.supp<-sim.data.f(par.values=parvec[supp.pars], par.names=supp.pars, var.vec=supp.vars, like.vec=supp.like, N=R2, nsim=nSim)
  
  # compute sufficient statistics for Durham Ferry release
  stats.DF<-apply(data.DF,2,get.stat.supp,"DF")
  
  # compute sufficient statistics for supplemental release
  stats.supp<-apply(data.supp,2,get.stat.supp,"supp")
  
  # simulate detections at the two lines of each dual array, conditional on total number of tags detected at each dual array
  data.A1<-sapply(unlist(lapply(stats.DF,extract.astat,"A1")),sim.data.f,par.values=parvec[A2.pars], par.names=A1.pars, var.vec=A1.vars, like.vec=A1.like, nsim=1)
  data.A2<-sapply(unlist(lapply(stats.DF,extract.astat,"A2"))+unlist(lapply(stats.supp,extract.astat,"A2")),sim.data.f,par.values=parvec[A2.pars], par.names=A2.pars, var.vec=A2.vars, like.vec=A2.like, nsim=1)
  data.B1<-sapply(unlist(lapply(stats.DF,extract.astat,"B1")),sim.data.f,par.values=parvec[B1.pars], par.names=B1.pars, var.vec=B1.vars, like.vec=B1.like, nsim=1)
  data.F1<-sapply(unlist(lapply(stats.DF,extract.astat,"F1"))+unlist(lapply(stats.supp,extract.astat,"F1")),sim.data.f,par.values=parvec[F1.pars], par.names=F1.pars, var.vec=F1.vars, like.vec=F1.like, nsim=1)
  data.G2<-sapply(unlist(lapply(stats.DF,extract.astat,"G2"))+unlist(lapply(stats.supp,extract.astat,"G2")),sim.data.f,par.values=parvec[G2.pars], par.names=G2.pars, var.vec=G2.vars, like.vec=G2.like, nsim=1)
  
  row.names(data.A1)<-paste("A1",c("AB","A0","B0"),sep=" ")
  row.names(data.A2)<-paste("A2",c("AB","A0","B0"),sep=" ")
  row.names(data.B1)<-paste("B1",c("AB","A0","B0"),sep=" ")
  row.names(data.F1)<-paste("F1",c("AB","A0","B0"),sep=" ")
  row.names(data.G2)<-paste("G2",c("AB","A0","B0"),sep=" ")
  
  # combine all data
  data.sim<-rbind(data.DF,data.supp,data.A1,data.A2,data.B1,data.F1,data.G2)
  
  # get seeds (haven't dealt with reducing model if necessary)
  
  seeds<-sapply(1:nSim,get.seeds.supp.wrap,stats.DF,stats.supp,data.A1,data.A2,data.B1,data.F1,data.G2,p.full=pars.supp.all)
  print("have seeds")
  
  for(i in 1:nSim)
  {
    seeds.tmp<-seeds[,i]
    
    # remove parameters whose seeds are on boundary of [0,1] or are NA - will not be using optim to estimate these parameters
    pars.0<-names(seeds.tmp[!is.na(seeds.tmp) & is.finite(seeds.tmp) & seeds.tmp==0])
    pars.1<-names(seeds.tmp[!is.na(seeds.tmp) & is.finite(seeds.tmp) & seeds.tmp==1])
    pars.na<-names(seeds.tmp[is.na(seeds.tmp) | !is.finite(seeds.tmp)])
    seeds.tmp<-seeds.tmp[-match(unique(c(pars.0,pars.1,pars.na)),names(seeds.tmp))]
    
    if(length(pars.na)>0) next
    
    # initial fit of model
    j<-1
    x.init<-optim(seeds.tmp,log.like.f.cmp,dat.vec=data.sim[,i], par.names=pars.supp.fit, var.vec=vars.all, like.vec=like.all, pars0=pars.0, pars1=pars.1,parsNA=pars.na) #, control=list(maxit=10000))

    # refit using current estimates as seeds    
    pars.old<-x.init$par; val.old<-x.init$value
    change<-10000
    while(change>0.001)
    {
      j<-j+1
      x.step<-optim(pars.old,log.like.f.cmp,dat.vec=data.sim[,i],par.names=pars.supp.fit, var.vec=vars.all, like.vec=like.all, pars0=pars.0,pars1=pars.1,parsNA=pars.na) #,control=list(maxit=10000))
      pars.new<-x.step$par; val.new<-x.step$value
      change<-abs(val.old-val.new)
      val.old<-val.new; pars.old<-pars.new
    }
    print(paste("iteration = ",j, sep=""))
    
    pars.new<-x.step$par # the MLEs from optim
    
    # re-insert the parameters whose seeds were 0, 1, or NA
    pars.new<-c(pars.new,rep(0,length(pars.0)),rep(1,length(pars.1)),rep(NA,length(pars.na)))
    names(pars.new)<-c(pars.supp.fit[-match(unique(c(pars.0,pars.1,pars.na)),pars.supp.fit)],pars.0,pars.1,pars.na)
    pars.new<-pars.new[pars.supp.fit]
    sA<-get.sA(pars.new)
    sR<-get.sR(pars.new)
    pars.new<-c(pars.new,sA,sR)
    names(pars.new)[length(pars.new)+(-1:0)]<-c("sA","sR")
    
    mat.out[,i+1]<-unlist(pars.new)
  } # end of simulation i
  
  return(mat.out)
}

get.est.mat<-function(x,nSim,parnames=pars_est)
{
  # reshape output of apply(...,sim.supp.seeds.pari.R1.R2)
  out.mat<-matrix(x,nrow=nSim+1,byrow=T)
  dimnames(out.mat)<-list(c("true",paste("Sim",1:nSim,sep=" ")),parnames)
  return(out.mat)
}

get.est.mat.wrap<-function(i,xx,Nsim,Pnames=pars_est)
{
  # wrapper for get.est.mat so can use with sapply
  return(get.est.mat(xx[,i],Nsim,Pnames))
}

log.like.f<-function(pars,dat.vec,par.names,var.vec,like.vec,pars1=c(),pars0=c(),parsNA=c())
{
  # returns the negative log-likelihood value for multinomial model specified in like.vec for data in dat.vec and values in pars
  
  #### Arguments:
  
  # pars=numeric vector of parameter values
  # par.names=character vector of parameter names in same order as pars
  # var.vec=character vector of variables (as in USER model)
  # pars1, pars0 are character vectors of parameters that are fixed to 1 and 0, respectively - will not try to estiamte them numerically
  # parsNA = character vector of parameters that are set to NA - will not try to estimate them

  #  Note: parameters named in pars1, pars0, parsNA should not be included in the vector "pars" passed to the function.
  
  #### 
  
  # Put back in values for parameters that were removed because seeds were either 0 or 1, and fix up par.names accordingly
  if(length(pars1)>0) {pars<-c(pars,rep(1,length(pars1))); par.name<-par.names[!(par.names %in% pars1)]; par.names<-c(par.names,pars1)}
  if(length(pars0)>0) {pars<-c(pars,rep(0,length(pars0))); par.names<-par.names[!(par.names %in% pars0)]; par.names<-c(par.names,pars0)}
  # for NA parameters, give value of 0.5; they do not contribute to value of likelihood at all, regardless of their value
  if(length(parsNA)>0) {pars<-c(pars,rep(0.5,length(parsNA))); par.names<-par.names[!(par.names %in% parsNA)]; par.names<-c(par.names,parsNA)}
  
  # correct pars vector - none of the parameters can be outside [0,1]
  pars<-abs(pars)
  pars[pars>0.9999999999]<-0.9999999999
  pars[pars<0]<-0.000001
  
  for(i in 1:length(par.names)) assign(par.names[i],pars[i])
  f1<-function(x) eval(parse(text=x))
  f1.cmp<-cmpfun(f1)
  for(i in 1:length(var.vec)) assign(names(var.vec[i]),f1.cmp(var.vec[i]))
  l.vec<-sapply(like.vec,f1.cmp)
  if(sum(is.na(l.vec))>0) print("warning: probability vector in likelihood has one or more NA components")
  
  # if l.vec=0 (i.e., if parameters=0 or 1), then log(l.vec)=Inf.
  # remove cells with counts = 0 (they do not contribute to likelihood)
  neglogl<-(-1)*sum(dat.vec[dat.vec>0]*log(l.vec[dat.vec>0]))
  return(neglogl)
}

log.like.f.cmp<-cmpfun(log.like.f)


get.C1<-function(x) 
  {
    # compute proportion of parameter estimates that are neither NA nor infinite
    # assumes first entry is true value
    x<-x[-1]; return(sum(!is.na(x) & is.finite(x))/length(x))
  } 

get.C2<-function(x,cutoff=1.1) 
  {
    # compute proportion of parameter estimates greater than specified cutoff value
    # assume first entry of x is true value
    # remove NA And infinite values
    x<-x[-1]; return(sum(!is.na(x) & is.finite(x) & x>cutoff)/length(x[!is.na(x) & is.finite(x)]))
  } 

get.C3<-function(x) 
  {
    # compute standard deviation of vector of parameter estimates; remove NA and infinite values
    # assumes first entry is true value
    x<-x[-1]; return(sd(x[!is.na(x) & is.finite(x)]))
  } 

get.C4<-function(x)
{
  # compute absolute value of bias of vector of parameter estimates; remove NA and infinite values
  # assumes first entry is true value
  value<-x[1]
  est<-x[-1]
  dif<-abs(value-mean(est[!is.na(est) & is.finite(est)]))
  return(dif)
}

get.mean<-function(x) 
  {
    # compute mean of vector of parameter estimates; remove NA and infinite values
    # assumes first entry is true value
    x<-x[-1]; return(mean(x[!is.na(x) & is.finite(x)]))
  } 

get.max<-function(x) 
  {
    # compute maximum of vector of parameter estimates; remove NA and infinite values
    # assumes first entry is true value
    x<-x[-1]; return(max(x[!is.na(x) & is.finite(x)]))
  } 

extract.astat=function(x,site) 
  {
    # extract "a" statistic at specified site from single simulation 
    # where a = number of tags detected at site (pooled over dual array, if applicable)
    x$a[site]
  }

extract.bstat<-function(x,site)
{
  # extract "b" statistic at specific site from single simulation
  # where b = number of tags detected both at site and downstream 
  x$b[site]
}

calculate.Ci<-function(dat.object,criterion="C1")
{
  # for given release size combination, evaluate summarization pertinent to specified criterion
  
  ### Arguments
  # dat.object = list of length 3 (e.g., dat[[1]])
  # -- element 1 = R_DF = sample size of primary release group at Durham Ferry
  # -- element 2 = R_supp = sample size of supplement release group
  # -- element 3 = matrix of simulated MLEs for parameters (row 1 = true values); nrows = N_sim + 1; ncol = number of estimated parameters
  
  # criterion = one among c("C1","C2","C3","C4")
  # ---- C1: proportion of simulations that generate estimates for parameter
  # ---- C2: proportion of simulations that generate estimates > cutoff value
  # ---- C3: standard deviation of simulation distribution of parameter estimates
  # ---- C4: bias of parameter estimates
  
  # note that default cutoff value for get.C2 is 1.1. If need to compute for different cutoff value, change in get.C2 definition
  
  fun.tmp<-get(paste("get",criterion,sep="."))
  C_tmp<-apply(dat.object[["mle_mat"]],2,fun.tmp)
  C_tmp<-c(dat.object[["R_DF"]],dat.object[["R_supp"]],C_tmp)
  names(C_tmp)[1:2]<-c("R_DF","R_supp")
  
  return(C_tmp)
}

get.minR.1param=function(param,dat.list,R.type=c("DF","supp"),c1=0.95,c2=0.95,Max.se=0.05,Nsims=N_sim,Max.bias=0.05)
{
  # for param, gets minimum R for specified release type according to criteria C1, C2, and C3 and C4 from output of sim.supp.seeds.pari.R1.R2 (applied as list)
  # dat.list = dat
  # R.type = release type (DF = primary, supp = supplemental)
  # c1 and c2 = cutoff percentage for C1 and C2
  # Max.se = maximum se allowed for C3
  # Max.bias = maximum bias allowed for C4
  # Nsims = number of simulations
  
  # this function is used in function min.samp.size
  
  # calculate C1, C2, C3, and C4 matrices and convert to data.frames
  c1.full<-data.frame(do.call(rbind,lapply(dat.list,calculate.Ci,"C1")))
  c2.full<-data.frame(do.call(rbind,lapply(dat.list,calculate.Ci,"C2")))
  c3.full<-data.frame(do.call(rbind,lapply(dat.list,calculate.Ci,"C3")))
  c4.full<-data.frame(do.call(rbind,lapply(dat.list,calculate.Ci,"C4")))
  
  # identify column index of desired release size
  R.index<-ifelse(R.type=="DF",1,ifelse(R.type=="supp",2,NA))
  if(is.na(R.index)) return(NA)
  
  # limit each ci.mat to the release size and parameter of interest
  c1.mat<-c1.full[,c(R.index,which(names(c1.full)==param))]
  c2.mat<-c2.full[,c(R.index,which(names(c2.full)==param))]
  c3.mat<-c3.full[,c(R.index,which(names(c3.full)==param))]
  c4.mat<-c4.full[,c(R.index,which(names(c4.full)==param))]
  
  # calculate minimum release size necessary to achieve criteria
  min.c1<-min(c1.mat[c1.mat[,2]>=c1,1]); if(!is.finite(min.c1)) min.c1<-NA
  min.c2<-min(c2.mat[c2.mat[,2]<=(1-c2),1]); if(!is.finite(min.c2)) min.c2<-NA
  min.c3<-min(c3.mat[c3.mat[,2]<=Max.se,1]); if(!is.finite(min.c3)) min.c3<-NA
  min.c4<-min(c4.mat[c4.mat[,2]<=Max.bias,1]); if(!is.finite(min.c4)) min.c4<-NA
  if(!is.na(min.c1)) N.ests<-round(c1.mat[c1.mat[,1]==min.c1,2]*Nsims,0) else N.ests<-NA
  if(sum(!is.na(N.ests))>0) names(N.ests)<-c1.full[c1.mat[,1]==min.c1,setdiff(c(1,2),R.index)]
  
  minR<-max(c(min.c1,min.c2,min.c3,min.c4))
  return(c(minR,min.c1,min.c2,min.c3,min.c4,N.ests))
  
  # named elements are the number of successful estimates for the possible sample sizes of the other release
}


