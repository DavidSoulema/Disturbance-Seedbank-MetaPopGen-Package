#############################################

### 6. INITIAL POPULATION GENERATION

#############################################


### To build a random but realistic initial population & seedbank to uniform all simulations.
  # Simulating with a consistent (N1,SB1) for all generations enables to focus on the impact of the disturbance.


# Initial variables
Tries <- 100

N0 <- array(10,dim=c(m,n,z))
dimnames(N0) <- list(genotype=c(1:m),deme=c(1:n),age=c(1:z))

SB0 <- array(100,dim=c(m,n))
dimnames(SB0) <- list(genotype=c(1:m),deme=c(1:n))

res.starting.population <- get_starting_population(input.type, demographic.data,
                                                   N0,
                                                   SB0,
                                                   sigma,
                                                   phi_F, phi_M,
                                                   fec.distr_F = "poisson", fec.distr_M = "poisson",
                                                   mu,
                                                   migration="forward", migr=NULL,
                                                   delta, 
                                                   recr.dd="settlers",
                                                   kappa0,
                                                   psi,
                                                   gamma,
                                                   seed_life_expectancy,
                                                   degeneration_type = "linear",
                                                   T_max,
                                                   Tries)

N1 <- res.starting.population$N1
SB1 <- res.starting.population$SB1

save(N1, file="N1.Rdata")
save(SB1, file="SB1.Rdata")