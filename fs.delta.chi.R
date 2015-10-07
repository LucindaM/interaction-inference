library(foreach)
library(doParallel)
registerDoParallel(cores=20)

#setwd("") 
source('filter.study.setup.R')

fname <- 'fs.delta.chi.rda'

# parameters
npart <- 5000
nfil <- 5
delta.set <- c(0.1, 0.5, 1, 2, 6, 12, 52, 365)
chi.set <- c(0, 0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 1, 1.01, 1.025, 1.05, 1.1, 1.2, 1.5, 2)

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

result <- mapply(c,result,SIMPLIFY=F)
save(result,file=fname)
