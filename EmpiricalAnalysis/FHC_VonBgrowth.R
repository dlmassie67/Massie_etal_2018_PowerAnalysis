#### Load libraries
library(MCMCpack) # rwish function
library(R2jags)


#read in data
dat<-read.csv("FHCLengthAtAgeData.csv")

#remove prob pop's
dat<-dat[dat$River!="North Raccoon", ] #remove population due to unrealistic Linf values (see Massie et al.2018 for details)
length(unique(dat$River)) #check that there are 25 populations

#must drop levels before ordering river in next step
dat <- droplevels(dat) 

#river name column
dat <- dat[order(dat$River), ]
dat$River_name <- dat$River
dat$River <- as.numeric(dat$River)
dat <- droplevels(dat)

#final check:
length(unique(dat$River)) #should be 25
unique(dat$River) #should be in a consecutive numeric order
unique(dat$River_name) #does not include removed pops above
#################################################################
########## BUGS CODE ############################################
#################################################################
sink("model.txt")
cat("
    model{
    for(i in 1:n){
    y[i] ~ dnorm(y.hat[i], tau.y)
    y.hat[i] <- Linf[g[i]] * (1-exp(-k[g[i]] * (age[i] - t0[g[i]] )))
    }

    tau.y <- pow(sigma.y,-2)
    sigma.y ~ dunif(0,100)

    # Parameters modeled on log-scale
    for(j in 1:J){
    Linf[j] <- exp(BB[j,1])
    k[j] <- exp(BB[j,2])
    t0[j] ~ dnorm(mu.t0, tau.t0)
    BB[j,1:K] ~ dmnorm (BB.hat[j,], Tau.B[,])
    BB.hat[j,1] <- mu.Linf
    BB.hat[j,2] <- mu.k
    }

    # Priors for population-average parameters
    tau.t0 <- pow(sigmat0,-2)
    sigmat0 ~ dunif(0.001,1)

    mu.Linf ~ dnorm(0,.0001)
    mu.k ~ dnorm(0,.0001)
    mu.t0 ~ dnorm(0,.0001)



    # Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }


     } # end model
    ",fill=TRUE)
sink()
################
J <- length(unique(dat$River))


# Create identity matrix for Wishart dist'n
#!!Number of varying parameters to estimate (K)
K <- 2

# Create identity matrix for Wishart dist'n
W <- diag(K)

# load data
data <- list(y = dat$TL, age = dat$Age, g = dat$River, n = dim(dat)[1],
             J = J, W=W, K=K)

# Initial values

inits <- function(){list(mu.Linf = rnorm(1,3,0.001), mu.k = rnorm(1,1,0.001), mu.t0 = rnorm(1,0.7,0.001),
                         sigma.y = runif(1,1,10), sigmat0 = runif(1), 
                         BB=array(c(rep(log(950) +rnorm(1,0.01,0.01),J),rep(log(0.04)+rnorm(1,0.001,0.1),J)),
                                  c(J,K)), Tau.B=rwish(K+1,diag(K)) ) }

# Parameters monitored
params1 <- c("mu.Linf", "mu.k", "mu.t0", "sigma.y","BB","Sigma.B","rho.B", "sigmat0")
rep(log(-2+10)+rnorm(1,0.01,0.1),J)



# MCMC settings
ni <- 100000
nt <- 3
nb <- 50000
nc <- 3
###########################
start.time = Sys.time()  

out <- jags(data = data, inits = inits, parameters.to.save = params1, 
                     model.file = "model.txt", n.chains = nc, 
                     n.thin = nt, n.iter = ni, n.burnin = nb)

#Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
##############################
##check convergence
# Find which parameters, if any, have Rhat > 1.1
which(out$BUGSoutput$summary[, c("Rhat")] > 1.1) #three increase iterations

traceplot(out)

# Save output as RDS file
saveRDS(out, file=paste0("FHCgrowth",1,'.rds'))


