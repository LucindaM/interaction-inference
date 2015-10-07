library(pomp)
library(foreach)
library(doParallel)
registerDoParallel(cores=2)  

setwd("/Users/Xi/Desktop/New/")
load('init.rda')

#set parameters
rmeasure <- "poisson_rmeasure"
dmeasure <- "poisson_dmeasure"
N <- 10000000
mu <- 0.02
beta <- 70
gamma <- 26
delta_o <- 365
theta <- 0
chi_o <- 1
omega <- 5e-7
rho.1 <- 1
rho.2 <- 1
# time steps per year
ts <- 1000
dt <- 1/ts
# samples per year
smp <- 365
# number of years
nt <- 1000
tmax <- nt
beta.sd <- 0.005

modelStem <- "msi_den_mod"
modelCfile <- paste(modelStem,".c",sep="")
solib <- paste(modelStem,.Platform$dynlib.ext,sep="")

if (!file.exists(modelCfile))
  stop("file ",sQuote(modelCfile)," not found")

rexe <- paste(R.home(),"bin/R",sep="/")
pomp.inc <- system.file("include/pomp.h",package="pomp")
pomp.lib <- system.file("libs/pomp.so",package="pomp")
system(paste("cp",pomp.inc,".",sep=" "))
system(paste("cp",pomp.lib,".",sep=" "))
system(paste(rexe,"CMD SHLIB -o",solib,modelCfile,pomp.lib,sep=" "))
system("rm pomp.h")

casedata <- read.csv("casedata.csv", row=1)

if (!file.exists(solib))
  stop("file ",sQuote(solib)," not found")

compnames <- array(NA, 80)
for(j in 1:80) compnames[j] <- paste("X",j,sep='')
compnames <- as.character(compnames)

make.pomp <- function(dt=0.001, lt=1000, smp=12, rmeasure, dmeasure,...) {
  pomp(
    data = casedata[c('time','y1','y2','y3','y4','yp1','yp2','yp3','yp4')],
    times= 'time',                         
    t0 = 0,
    rmeasure=rmeasure,
    dmeasure=dmeasure,
    zeronames=c('cases1','cases2','cases3','cases4','casesp1','casesp2','casesp3','casesp4'),
    PACKAGE=modelStem,
    rprocess = euler.sim(delta.t = dt, step.fun = "msi_euler_simulate", PACKAGE=modelStem),		
    obsnames = c('y1','y2','y3','y4','yp1','yp2','yp3','yp4'),
    comp.names = compnames,
    statenames = c(compnames,'N','W', 'cases1', 'cases2', 'cases3','cases4','casesp1','casesp2','casesp3','casesp4'),
    paramnames = c('beta.sd','mu','beta','gamma','delta','theta','chi','omega','rho.1','rho.2'),
    log.params = c('beta.sd','mu','beta','gamma','delta','theta','chi','omega', paste(compnames,".0",sep=''),'N.0'),
    logit.params = c('rho.1','rho.2'),
    initializer = function (params, t0, ...) {
      comp.ic.names <- paste(compnames,".0",sep='')
      states <- numeric(length(c(compnames,'N','W', 'cases1', 'cases2', 'cases3','cases4','casesp1','casesp2','casesp3','casesp4')))
      names(states) <- c(compnames,'N','W', 'cases1', 'cases2', 'cases3','cases4','casesp1','casesp2','casesp3','casesp4')
      frac <- exp(params[comp.ic.names])
      states['N'] <- round(exp(params['N.0']))
      states[compnames] <- round(states['N']*frac/sum(frac))
      states
    },
    transform.fn = function (params, log.params, logit.params, ...) {
      x <- as.array(params)
      nx <- dimnames(x)
      dx <- dim(x)
      dim(x) <- c(dx[1],prod(dx)/dx[1])
      log.ind <- match(log.params,nx[[1]])
      logit.ind <- match(logit.params,nx[[1]])
      x[log.ind,] <- log(x[log.ind,])
      logit <- function(p){log(p/(1-p))}
      x[logit.ind,] <- logit(x[logit.ind,])
      if (length(dx)>1) {
        dim(x) <- dx
        dimnames(x) <- nx
      } else {
        dim(x) <- NULL
        names(x) <- nx[[1]]
      }
      x
    },
    untransform.fn = function (params, log.params, logit.params, ...) {
      x <- as.array(params)
      nx <- dimnames(x)
      dx <- dim(x)
      dim(x) <- c(dx[1],prod(dx)/dx[1])
      log.ind <- match(log.params,nx[[1]])
      logit.ind <- match(logit.params,nx[[1]])
      x[log.ind,] <- exp(x[log.ind,])
      expit <- function(r){1/(1+exp(-r))}
      x[logit.ind,] <- expit(x[logit.ind,])
      if (length(dx)>1) {
        dim(x) <- dx
        dimnames(x) <- nx
      } else {
        dim(x) <- NULL
        names(x) <- nx[[1]]
      }
      x
    } 
  )
}

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

###### particle filtering

npart <- 10
nfil <- 2
delta.set <- c(100, 243, 365)
chi.set <- c(0.95,1.1, 1.3)

ind <- 1

joblist <- vector(mode='list',length=length(delta.set)*length(chi.set)*nfil)

for( jd in 1:length(delta.set)){
  for (jc in 1:length(chi.set)){
    for( jt in 1:nfil){
      params['delta'] <- delta.set[jd]
      params['chi'] <- chi.set[jc]
      joblist[[ind]] <- list(params=params, npart=npart)
      ind <- ind+1
    }
  }
}


tic <- Sys.time()
result <- foreach(i=1:length(joblist)) %dopar%
{
  params<-joblist[[i]]$params
  pf <- pfilter(msi.po,
                params=msi.po@userdata$transform.fn(params, msi.po@userdata$log.params, msi.po@userdata$logit.params),
                pred.mean=F, max.fail=500, tol=1e-68, Np=joblist[[i]]$npart)
  out <- list(params=params,pf=pf)
  out
}

toc <- Sys.time()
print(toc-tic)

dat <- result

delta <- sapply(dat, function(x)x$params['delta'])
chi <- sapply(dat, function(x)x$params['chi'])
loglik <- sapply(dat, function(x)x$pf$loglik)

dat1<-cbind(delta,chi,loglik)
rownames(dat1)<-c(1:nrow(dat1))

dat_paramnames <- c("delta","chi")
dat_mle <- dat1[which.max(dat1[,"loglik"]),][dat_paramnames]

params_mle <- c(
  beta.sd=beta.sd,
  mu=mu,
  beta=beta,
  gamma=gamma,
  delta=dat_mle[1],
  theta=theta,
  chi=dat_mle[2],
  omega=omega,
  rho.1=rho.1,
  rho.2=rho.2,
  N.0=N)

############# iterated filtering

my_rw.sd <- 0.00001
my_cooling.fraction.50 <- 0.00005


t1 <- Sys.time()
mifs <- foreach(i=1:20,.packages='pomp', .combine=c) %dopar%
{
  mif2(
    msi.po,
    start=params_mle,
    Np=20,
    Nmif=5,
    cooling.type="geometric",
    cooling.fraction.50=my_cooling.fraction.50,
    transform=TRUE,
    rw.sd=rw.sd(
      delta=my_rw.sd,
      chi=my_rw.sd
    )
  )
}

t2 <- Sys.time()
print(t2-t1)


