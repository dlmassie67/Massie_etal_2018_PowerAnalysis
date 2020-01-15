# Load required packages
library(R2jags)
library(MASS)
library(dplyr)
library(MCMCpack)
library(data.table)
# # Use ggmcmc package for assessing convergence
library(ggmcmc)
library(gridExtra)
library(ggthemes)
library(coda)


# Read in rds file that contains all mcmc output from Von Bertalanffy growth model
params <- readRDS('FHCgrowth1.rds')
str(params)

### input parameters
# b.true is a vector of effect sizes that describe the % change in growth paramater with a 1 unit increase in effect
#we illistrate by assessing the ability to detect a water temperature effect on the growth coefficient (K) 
b.true=(0.05) 

# Number of sites to sample (rivers, lakes, etc.)
n.river=15
# Simply a vector of made up of a range of covariate values that the growth cofficient, K, will be a function of
water.temp=seq(10,32,length=n.river)
# Standardize water temperature - helps with interpretation of parameters and convergence
water.temp=(water.temp-mean(water.temp))/sd(water.temp)

# "True" Population number of fish for each site
n.popn=5000 
# Number of fish sampled from each site (from n.popn) to use for fitting the growth model
n.sample=500
# Proportion at age information used to simulate length at age data
age.dist=read.table('age_structure.txt',header=T,sep='\t',row.names=NULL,as.is=c(1:2)) 

# Number of covarying parameters (ie K, Linf, t0)
K <- 3
# Create identity matrix for Wishart dist'n
W <- diag(K)
# Number of sites indicator for jags model
J <- n.river

######## Number of simulations (iterations) to run
num.sim <-100

# Chain length from analysis
chainLength <- params$BUGSoutput$n.sims
# Select thinned steps in chain for posterior predictions to ensure we take values from length of posterior
#accounts for parameter uncertainty from the empirical growth analysis
ID = round(seq( 1 , chainLength , length=num.sim ))


# The grand-mean L-infinity on the raw (natural) scale, across all sites
mean.Linf.true=exp(params$BUGSoutput$sims.list$mu.Linf[ID]) # raw scale (mm)
# The grand mean K on the log-scale, across all sites
grand.mean.k.true=params$BUGSoutput$sims.list$mu.k[ID]
# The grand-mean t0 (across all sites) on raw scale
mean.t0.true=params$BUGSoutput$sims.list$mu.t0[ID]
# Residual standard deviation on raw scale #sigmaY
sigma.age.true=params$BUGSoutput$sims.list$sigma.y[ID]


# create empty results object to bind jags results into. This is a list, with each element of the list corresponding
# to a different effects size, i.e., corresponding to a different value in b.true
results <- list()
for(m in 1:length(b.true)){
  results[[m]] <- list(NULL)
}


#Create container to to hold posterior estimates of b from each simulation run
bPosterior <- list()

################################################start simulation
# Loop over effect sizes in b
start.time = Sys.time()  #start timer
for(j in 1:length(b.true)){

  
  # parameters to save during simulations
  ### matrix to save result
  pars.to.save <- c("b","Linfbar","muK","t0bar")
  resultCols<-c('parameter','Mean','SD','2.5%','97.5%','rHat','iterNumber','effect_size')
  res = array(NA,dim=c(length(pars.to.save),length(resultCols))) %>%
    data.frame()
  names(res)<-resultCols
  res$parameter<-pars.to.save
  
  # Loop over the number of simulations
  for(sim in 1:num.sim){
    
    # # Generate true K that is a function of water temperature
    # mean.k.true=exp(grand.mean.k.true[sim]+b.true[j]*water.temp)
    
    ###################### Generate true popln
    Linf.true=rep(NA,n.river) #placeholder for mu.Linf for each pop/n.river
    k.true=rep(NA,n.river)
    t0.true=rep(NA,n.river)
    # Variance covariance matrix for mvnorm
    vacov.true=params$BUGSoutput$sims.list$Sigma.B[ID[sim],,]
    
    # Generate L-infinity, K, and t0 for each site, from a multivariate normal (i.e., we want to account for dependency in these parameters)
    for(r in 1:n.river){
      para.true=mvrnorm(n=1,mu=c(log(mean.Linf.true[sim]),grand.mean.k.true[sim]+b.true[j]*water.temp[r],mean.t0.true[sim]),Sigma=vacov.true)
      Linf.true[r]=round(exp(para.true[1]),2)
      k.true[r]=round(exp(para.true[2]),2)
      t0.true[r]=round(para.true[3],2)
    } #end r

    
    dd.popn=matrix(NA,nrow=1,ncol=3) # initialize data matrix, will delete later
    colnames(dd.popn)=c('RiverCode','Age','TL')
    # Generate length-at-age data for each river
    for(r in 1:n.river){
      n.byAge=round(n.popn*age.dist[,'Portion'])
      Age=rep(age.dist[,'Age'],n.byAge)
      TL=round(Linf.true[r]*(1-exp(-k.true[r]*(Age-t0.true[r])))+(rnorm(length(Age),mean=0,sd=sigma.age.true[sim])),2)
      RiverCode=rep(r,length(Age))
      dd.popn.tem=cbind(RiverCode,Age,TL)
      dd.popn=rbind(dd.popn,dd.popn.tem)
    }
    dd.popn=dd.popn[-1,] #delete 1st row of NAs	
    n.age.popn=nrow(dd.popn)
    
    
    
    #################### Sample true popln (sample n.sample from n.popn)
    dd.age=matrix(NA,nrow=1,ncol=3) #initialize matrix, delete later
    colnames(dd.age)=colnames(dd.popn)
    for(r in 1:n.river){
      dd.popn.tem=subset(dd.popn,dd.popn[,'RiverCode']==r)
      dd.age.tem=dd.popn.tem[sample(1:nrow(dd.popn.tem),size=n.sample,replace=F),]
      
      dd.age=rbind(dd.age,dd.age.tem)		
    }#end r
    dd.age=dd.age[-1,] #delete first row of NAs
    dd.age=dd.age[which(dd.age[,'TL']>0),] #remove those data with nonpositive TL value due to t0 setup and random number generating
    n.age=nrow(dd.age)
    
   
    
#####################fit growth model to sampled data
    
    datafit=list(y=dd.age[,3],age=dd.age[,2],g=dd.age[,1],n=nrow(dd.age),water.temp=water.temp,J = J, W=W, K=K)
    para=c('b','mu.Linf','mu.t0','sigma.y','mu.k',"Sigma.B","Linf","k","t0")
    inits <- function(){list(mu.Linf = rnorm(1,params$BUGSoutput$mean$mu.Linf,0.001), mu.k = rnorm(1,params$BUGSoutput$mean$k,0.001), 
                             mu.t0 = rnorm(1,params$BUGSoutput$mean$mu.t0,0.001),b=rnorm(1), 
                             sigma.y = runif(1,50,150), 
                             BB=array(c(rep(log(950) +rnorm(1,0.01,0.01),J),rep(log(0.04)+rnorm(1,0.001,0.1),J),rep(-0.03+rnorm(1,0.001,0.1),J)),
                                      c(J,K)), Tau.B=rwish(K+1,diag(K)) ) }
    
    # Jags could crash due to random number generation, so using try to avoid nullifying the whole simulation
    jagsWorked<-FALSE
    try(expr={
      fit=jags(data=datafit,inits=inits,parameters.to.save=para,model.file='model.txt',n.chains=3,n.iter=100000,n.burnin=50000,n.thin=3,DIC=T)
      jagsWorked<-TRUE
    })
    
    if(!jagsWorked){
      #if jags throws an error, then just save NAs as the result to indicate which models failed
      res[,c('Mean','SD','2.5%','97.5%','rHat')] <- as.numeric(NA)
      res$iterNumber<- sim
      res$effect_size <- b.true[j]
      results[[j]]<-rbind(results[[j]],res)
      next
    }
    
    #save the results summary and simulation info
    res[,c('Mean','SD','2.5%','97.5%','rHat')] <- fit$BUGSoutput$summary[c("b","mu.Linf","mu.k","mu.t0"),c(1,2,3,7,8)]
    res$iterNumber <- sim
    res$effect_size <- b.true[j]
    
    #bind with previous simulations
    results[[j]]<-rbind(results[[j]],res)
    
    print(sim)
    
    #!!!!!!!! Retain posterior samples for b - effect of temp
    name <- paste("b:",j,sim,sep='')
    bPosterior[[name]] <- fit$BUGSoutput$sims.list$b
    
  }# end sim loop
} # end b.true loop


end.time = Sys.time() #end timer
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

####################################################
#save last fit iteration to access MCMC chains/test convergence
saveRDS(fit, file=paste0("fit",1,'.rds'))

# Save output as RDS file
saveRDS(results, file=paste0("sim",1,'.rds'))

# Read simulation output back into R for further processing
simOutput <- readRDS(file='sim1.rds')


#check convergence
#read in MCMC chains for last iteration
fit <- readRDS(file='fit1.rds')
# Rhat > 1.1
which(fit$BUGSoutput$summary[, c("Rhat")] > 1.1) #three increase iterations
which(fit$BUGSoutput$summary[, c("n.eff")] < 100)

#traceplots
out.mcmc <- as.mcmc(fit)
# # Creates full report
S <- ggs(out.mcmc)
ggmcmc(S)


############Calculate detection probability
saveRDS(bPosterior, file="bPosterior1.rds") #save bPosterior list
bPosterior <- readRDS(file='bPosterior1.rds')  
bPosterior.matrix <- do.call(rbind, bPosterior)#convert bPosterior from list into matrix
sum(bPosterior.matrix > 0)/nrow(bPosterior.matrix)#prob that posterior estimates of b are greater than 0


############Calculate statistical power
# # Converve simOutput from a list into a matrix
out_matrix <- do.call(rbind, simOutput)
# Convert to dataframe
out_df <- data.frame(out_matrix)
# Select out only the effect of water temperature paramter (b)
out_b <- out_df[out_df$parameter=='b',]
# Create index showing if b is different than zero (i.e., significant): 1 = significant; 0 = not significant
out_b$Sig <- ifelse(out_b$X2.5. * out_b$X97.5. > 0, 1, 0)
# Calculate power as the proportion of interations (by effect size) that we detected an effect of water temperature,
# i.e., that the effect of b is different than zero
# Convert out_b to data.table for power summaries
out_b <- data.table(out_b)
Power <- out_b[, (sum(Sig, na.rm = TRUE)/num.sim),by = effect_size]
# Rename newly calculated power value to "power"
setnames(Power, old = "V1", new = "power")
Power


