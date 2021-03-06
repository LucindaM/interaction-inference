
source('params.R')
load('init.rda')

# dir with the model files
dir.mod <- '~/DenSerIntInf/model/'

# current dir
dir.current <- getwd()

# change dir to the dir with model
setwd(dir.mod)

datapath <- paste(dir.current,'/data',sep='')

source('msi_den_mod.R')
dyn.load(solib)

# change dir to current
setwd(dir.current)
# set non initial cond params
params <- c(
            beta.sd=beta.sd,
            mu=mu,
            beta=beta,
            gamma=gamma,
            delta=delta,
            theta=theta,
            chi=chi,
            omega=omega,
            rho.1=rho.1,
            rho.2=rho.2,
            N.0=N
)

# init cond params
icompnames <- array(NA, 80)
for(j in 1:80) icompnames[j] <- paste("X",j,".0",sep='')
icompnames <- as.character(icompnames)

lp <- length(params)
params <- c(params, rep(0,80))
names(params)[-(1:lp)] <- icompnames
params[-(1:lp)] <- as.numeric(init)

msi.po <- make.pomp(dt=0.005,lt=40,smp=smp,rmeasure=rmeasure, dmeasure=dmeasure)

null.params <- params 
