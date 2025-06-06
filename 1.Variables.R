setwd("C:\\Users\\david\\Documents\\DAVID\\CS\\DD\\Cours\\Master_Thesis\\ENVM7130_Master_Thesis\\Experiment")


#############################################

### 1. Variables

#############################################

# Set the variables of the experiment. Mostly similar to the classic MetaPopGen ones, with some add-ons.

library(MetaPopGen)

l <- 2              # alleles
m <- l*(l+1) / 2    # genotypes
n <- 1              # demes
z <- 10             # age classes (life expectancy expressed in number of generations)

# age-classes grouped in same characteristics (a generation might not correspond to a real transition from an age-class to the next one)
seedling <- 1:1
juvenile <- 2:3
mature <- 4:10

T_max <- 100         # simulation time (in number of generations)




# Initial population in each genotype, deme and age class

N0 <- array(10,dim=c(m,n,z))
dimnames(N0) <- list(genotype=c(1:m),deme=c(1:n),age=c(1:z))


# Initial population in each genotype, deme and age class

SB0 <- array(100,dim=c(m,n))
dimnames(SB0) <- list(genotype=c(1:m),deme=c(1:n))






#############################################

### Demographic data

#############################################


return_to_pre_fire <- 11 # time for the parameters to recover their pre-fire values (linear modelling back to the stable-state)


# Survival probabilities "sigma"

sigma <- array(0.8,dim=c(m,n,z,T_max),
               dimnames = list(genotype=c(1:m), deme=c(1:n), age=c(1:z), time=c(1:T_max)))

sigma[,,1,] <- 0.3

for (t in 1:return_to_pre_fire){
  sigma[,,1,t] <- 0.4 + (0.3-0.4)/(return_to_pre_fire-1)*(t-1) #linear decrease from 0.4 to 0.3
}



# Fecundities "phi"

phi_F <- array(0,dim=c(m,n,z,T_max),
               dimnames = list(genotype=c(1:m), deme=c(1:n), age=c(1:z), time=c(1:T_max)))

phi_F[,,mature,] <- 260    #mature are the only fertile ones


phi_M <- array(0,dim=c(m,n,z,T_max),
               dimnames = list(genotype=c(1:m), deme=c(1:n), age=c(1:z), time=c(1:T_max)))

phi_M[,,mature,] <- 1000000   #mature are the only fertile ones




# Carrying capacity "kappa0"

recr.dd <- "adults"
kappa0 <- array(1600,dim=c(n,T_max),
                dimnames = list(deme=c(1:n), time=c(1:T_max)))

for (t in 1:return_to_pre_fire){
  kappa0[,t] <- 6400 + (1600-6400)/(return_to_pre_fire-1)*(t-1) #linear decrease from 6400 to 1600
}




# Mutation probabilities "mu" (each row has to sum to 1)

mu <- array(1e-6 / 19, dim=c(l,l))
diag(mu) <- 1 - 1e-6   



# Dispersal probabilities "delta" between demes (each row has to sum to 1)

if (n==1){
  delta <- matrix(1, nrow=n, ncol=n) #if only one deme, no dispersal out of this one deme
  } else {
    delta <- matrix(0.1/(n-1), nrow=n, ncol=n)
    diag(delta) <- 1 - 0.1
  }
delta <- array(rep(delta,T_max),c(n,n,T_max))




##### NEW demographic data: DISTURBANCE

# Fire fatality rates "psi" (changes according to age classes z)

psi <- array(0.08, dim = c(1, z),
             dimnames = list("Probabilities of death by fire", paste0("age-class ", 1:z)))

#future add-on: control these probabilities according to fire intensity


# Generation time of the species in years

generation_time <- 1 # in years



##### Added demographic data: SEEDBANK

# Ideal probability of germination "gamma": gives the probability for a seed to germinate in ideal conditions (new seed in perfect environmental conditions)
  #Depends on time since fire

gamma <- array(0.25, dim = c(1, T_max),
               dimnames = list("Ideal probability of germination"=c(1:1), time=c(1:T_max)))

for (t in 1:return_to_pre_fire){
  gamma[,t] <- 0.30 + (0.25-0.30)/(return_to_pre_fire-1)*(t-1) #linear decrease from 0.30 to 0.25
}



# Seed life expectancy: gives the characteristic time after which seeds usually do not germinate or die 

seed_life_expectancy <- 2 #has to be expressed in number of generation time (i.e. if the species typical generation time is 2yrs, then 5 would mean the seed can survive for 10yrs)



# Degeneration type: select the model of degeneration that the seed will follow (linear or exponential)

#degeneration_type <- "linear"
degeneration_type <- "exponential"



# Death model: choose if the algorithm should kill the seeds after a certain age to reduce calculation or not ("True" / "False")

is_death <- TRUE
#is_death <- FALSE





#############################################

### Control parameters and data saving

#############################################

input.type = "array"
verbose <- F
save.res <- F
save.res.T <- seq(1,T_max)
save(list=ls(),file="Data.Experiment.RData")


