
require(pomp)

modelStem <- "msi_den_mod"
modelCfile <- paste(modelStem,".c",sep="")
solib <- paste(modelStem,.Platform$dynlib.ext,sep="")

if (!file.exists(modelCfile))
  stop("file ",sQuote(modelCfile)," not found")

rexe <- paste(R.home(),"bin/R",sep="/")
pomp.inc <- system.file("include/pomp.h",package="pomp")
pomp.lib <- system.file("libs/pomp.so",package="pomp")
system(paste("cp",pomp.inc,".",sep=" "))
system(paste(rexe,"CMD SHLIB -o",solib,modelCfile,pomp.lib,sep=" "))
system("rm pomp.h")

datapath <- "../data"
casedatafile <- "casedata.csv"
readcasedata <- file.path(datapath,casedatafile)
casedata <- read.csv(readcasedata, row=1)

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
    initializer = function (params, t0, statenames, comp.names, ...) {
      comp.ic.names <- paste(comp.names,".0",sep='')
      states <- numeric(length(statenames))
      names(states) <- statenames
      frac <- exp(params[comp.ic.names])
      states['N'] <- round(exp(params['N.0']))
      states[comp.names] <- round(states['N']*frac/sum(frac))
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



