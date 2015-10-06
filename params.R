## paramaters file
## observation model
rmeasure <- "poisson_rmeasure"
dmeasure <- "poisson_dmeasure"
####
N <- 10000000
mu <- 0.02
beta <- 70
gamma <- 26
delta <- 365
theta <- 0
chi <- 1
omega <- 5e-7
rho.1 <- 1
rho.2 <- 1

## 
# time steps per year
ts <- 1000
dt <- 1/ts
# samples per year
smp <- 365
# number of years
nt <- 1000

tmax <- nt

##
beta.sd <- 0.005

