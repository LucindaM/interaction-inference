
source('params.R')
load('init.rda')

source('msi_den_mod.R')
dyn.load(solib)

# set non initial cond params
params <- c(
  beta.sd=beta.sd,
  mu=mu,
  beta=beta,
  gamma=gamma,
  delta=delta_o,
  theta=theta,
  chi=chi_o,
  omega=omega,
  rho.1=rho.1,
  rho.2=rho.2,
  N.0=N)

# init cond params
icompnames <- array(NA, 80)
for(j in 1:80) icompnames[j] <- paste("X",j,".0",sep='')
icompnames <- as.character(icompnames)

lp <- length(params)
params <- c(params, rep(0,80))
names(params)[-(1:lp)] <- icompnames
params[-(1:lp)] <- as.numeric(init)

msi.po <- make.pomp(dt=0.005,lt=40,smp=smp,rmeasure=rmeasure, dmeasure=dmeasure)
