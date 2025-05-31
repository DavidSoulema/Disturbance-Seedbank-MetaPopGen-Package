#############################################

### 5. PRE SIMULATION FUNCTION

#############################################


# Function to generate a consistent and random population & seedbank to start with.
  # Use the DISTURBANCE SIMULATION FUNCTION with no disturbance for a specified period of time.
  # The "Tries" argument is to repeat the simulation several time and then average the results.

get_starting_population <- function(input.type, demographic.data,
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
                                    Tries) {
  
  ########
  # Initialisation
  ########
  
  # Set the simulation vector to all 0 to have undisturbed simulations 
  simulation_vector <- rep(0, T_max)
  
  # Retrieve some parameters
  m <- dim(N0)[1]
  n <- dim(N0)[2]
  z <- dim(N0)[3]
  
  # Set the arrays that will store the results of each simulation before averaging
  pool_N <- array(NA, dim=c(m,n,z,Tries))
  dimnames(pool_N) <- list(genotype=c(1:m),deme=c(1:n),age=c(1:z), Try=c(1:Tries))
  
  pool_SB <- array(NA, dim=c(m,n,Tries))
  dimnames(pool_SB) <- list(genotype=c(1:m),deme=c(1:n), Try=c(1:Tries))
  
  ########
  # Run many simulations then average the results
  ########
  
  for (try in 1:Tries){
    print(paste0("Try number ", try))
    res <- sim.DISTURBANCE.metapopgen.monoecious(input.type, demographic.data,
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
                                                 simulation_vector = simulation_vector,
                                                 save.res=F, save.res.T = seq(1,T_max),
                                                 verbose=F)
                                          
    N_res <- res$N[,,,T_max, drop = FALSE]
    SB_res <- res$SB[T_max,,, drop = FALSE]

    pool_N[,,,try] <- N_res
    pool_SB[,,try] <- SB_res
  }
  
  # Get N1 and SB1 as the average of all tries
  
  N1 <- array(NA,dim=c(m,n,z))
  dimnames(N1) <- list(genotype=c(1:m),deme=c(1:n),age=c(1:z))
  
  N1 <- apply(pool_N, c(1, 2, 3), function(x) floor(mean(x)))
  
  SB1 <- array(NA,dim=c(T_max,m,n))
  dimnames(SB1) <- list(generation=c(1:T_max),genotype=c(1:m),deme=c(1:n))
  
  SB1 <- apply(pool_SB, c(1, 2), function(x) floor(mean(x)))
  
  return(list(N1=N1,SB1=SB1))
}

