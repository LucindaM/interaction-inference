library(pomp)
source("params.R")

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

dyn.load(solib)

load('init.rda')

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
  N.0=N)

compnames <- array(NA, 80)
for(j in 1:80) compnames[j] <- paste("X",j,sep='')
compnames <- as.character(compnames)


sim <- pomp(
  data = data.frame(time=seq(0,40,length=481),
                    y1=NA, y2=NA,y3=NA,y4=NA,yp1=NA,yp2=NA,yp3=NA,yp4=NA),
  times= 'time',                         
  t0 = 0,
  rmeasure=rmeasure,
  #dmeasure=dmeasure,
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

sim <- simulate(sim,params=params)


