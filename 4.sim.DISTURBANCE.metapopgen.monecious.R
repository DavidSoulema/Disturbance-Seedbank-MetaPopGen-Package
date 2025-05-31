#############################################

### 4. DISTURBANCE SIMULATION FUNCTION

#############################################

# This function has been redesigned from the MetaPopGen 'sim.metapopgen.monoecious' for 2 purposes:
  # 1. refine the recruitment model for plants, incorporating a seedbank.
  # 2. allowing disturbances such as fire to return periodically.

# It necessitates new arguments: 
  # generation_time: number of years per generation of the species
  # disturbance_return_interval: periodicity, in years, of the disturbance - if you want to study a disturbance that is not periodic, feel free to modify the 'get_disturbance_vector' function
  # psi: a new demographic parameter defining the disturbance fatality rate, according to age classes z)
  # gamma: another new demographic parameter defining the ideal germination rate, according to time since disturbance
  # seed_life_expectancy: time, in generations, after which the seeds are dead / not fertile anymore
  # degeneration_type = "linear" or "exponential": describes the pace of the decrease in fertility over generations


# The idea of the simulation is the following:
## 1. Let the normal MetaPopGen run the simulation for a certain time until reaching a stable state ('stable_state_time')
## 2. Run a following number of generations T_max with disturbance returning on a specific pattern ('fire_return_interval)
## 3. Run a final undisturbed simulation of a specified number of generation after the last disturbance

# For this project, the chosen disturbance is fire and will return periodically every X years (fire_return_interval). Modify this pattern for the disturbance you want to study.
# The code uses a vector full of 0, except for the generations when the disturbance comes. A 1 will be put for the first disturbance, then a 2, etc...
# The simulation then checks at each generation if the vector's position corresponding to the generation is a 0 (then runs undisturbed), or not (then simulate the disturbance before running).



# MAIN FUNCTION: this has been adapted from the MetaPopGen function sim.metapopgen.monoecious

sim.DISTURBANCE.metapopgen.monoecious <- function(input.type, demographic.data,
                                               N1,
                                               SB1,
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
                                               T_max, #time of the disturbed simulation (it would be the one of the original function)
                                               simulation_vector = disturbance_vector,
                                               save.res=F, save.res.T = seq(1,T_max),
                                               verbose=F) {
  
  
  ##########################################################################
  
  # Initial definitions
  
  ##########################################################################
  
  T_total <- length(simulation_vector)
  save.res.T <- seq(1,T_total)
  
  
  if (input.type=="data.frame") {
    
    #print("Input type = data.frame")
    
    a <- metapopgen.input.convert.monoecious(demographic.data)
    N1    <- a[[1]]
    sigma <- a[[2]]
    phi_M <- a[[3]]
    phi_F <- a[[4]]
    rm(a)
    
    
  } else {
    if (input.type == "array") {
      
      #print("Input type = array")
      
    } else {
      stop("Unknown value for argument input.type. It must be either data.frame, array or txt")
      
    }
  }
  
  
  
  # Define basic variables
  
  m <- dim(N1)[1]                 # Number of genotypes
  l <- (sqrt(1+8*m)-1)/2          # Nuber of alleles
  n <- dim(N1)[2]                 # Number of demes
  z <- dim(N1)[3]                 # Number of age-classes
  
  # Check fec.distr_F and fec.distr_M
  if(!fec.distr_F %in% c("fixed","poisson")) stop(paste("Unknown parameter value for fec.distr_F:",fec.distr_F))
  if(!fec.distr_M %in% c("fixed","poisson")) stop(paste("Unknown parameter value for fec.distr_M:",fec.distr_M))
  
  
  # If only one age-class and recr.dd=="adults", gives an error
  if (z == 1 & recr.dd == "adults") {
    stop("Detected only one age class (z=1) and recruitment probability dependent on adult density (recr.dd == 'adults'). This combination is not supported. Use recr.dd == 'settlers' instead.")
  }
  
  ##########################################################################
  
  # Check if input data are time-dependent or not
  
  ##########################################################################
  
  # Survival
  if (is.na(dim(sigma)[4])) {
    sigma <- array(rep(sigma,T_total),c(m,n,z,T_total))
  }
  
  # Female fecundity
  if (is.na(dim(phi_F)[4])) {
    phi_F <- array(rep(phi_F,T_total),c(m,n,z,T_total))
  }
  
  # Male fecundity
  if (is.na(dim(phi_M)[4])) {
    phi_M <- array(rep(phi_M,T_total),c(m,n,z,T_total))
  }
  
  # Dispersal
  if (is.na(dim(delta)[3])) {
    delta <- array(rep(delta,T_total),c(n,n,T_total))
  }
  
  # Carrying capacity
  if (is.vector(kappa0)) {
    kappa0 <- array(rep(kappa0,T_total),c(n,T_total))
  }
  
  
  
  ##########################################################################
  
  # Initialize state variables
  
  ##########################################################################
  
  #print("Initializing variables...")
  if (save.res){
    N <- N1
    rm(N1)
    
    SB <- SB1
    rm(SB1)
    
  } else {
    N       <- array(NA,dim=c(m,n,z,T_total))
    dimnames(N) <- list(genotype=c(1:m),deme=c(1:n), age_class=c(1:z), generation=c(1:T_total))
    
    N[,,,1] <- N1                             
    
    L       <- array(NA,dim=c(m,n,T_total))
    S       <- array(0,dim=c(m,n,T_total))
    
    
    SB <- array(0,dim=c(T_total,m,n))
    dimnames(SB) <- list(generation=c(1:T_total),genotype=c(1:m),deme=c(1:n))
    
    SB[1,,] <- SB1
  }

  
  
  ##########################################################################
  
  # Simulate metapopulation genetics
  
  ##########################################################################
  if (save.res){
    dir.res.name <- paste(getwd(),format(Sys.time(), "%Y-%b-%d-%H.%M.%S"),sep="/")
    dir.create(dir.res.name)
    
    if (1 %in% save.res.T) {
      file.name <- "N1.RData"
      save(N,file=paste(dir.res.name,file.name,sep="/"))
    }
  }
  
  
  #print(paste0("Disturbances to run: ", max(disturbance_vector)))
  
  #print(paste0("Total number of generations to run: ", T_total))
  
  
  # Keeping track of the time since disturbance (TSD), which will be used for post-disturbance parameters
  TSD <- 0
  past_disturbance <- FALSE
  
  for (t in 1 : T_total) {

    TSD <- TSD + 1
    
    
    # If we are still in the stable-state transition, use the undisturbed demographic parameter (their maximum temporal value)
    # this section can be optimized (very heavy to reupdate the ..._sim parameters at each step)
    if (past_disturbance == FALSE) {

      sigma_sim <- array( rep(sigma[,,,dim(sigma)[4]], times = dim(sigma)[4]),  dim = dim(sigma))
      phi_F_sim <- array( rep(phi_F[,,,dim(phi_F)[4]], times = dim(phi_F)[4]),  dim = dim(phi_F))
      phi_M_sim <- array( rep(phi_M[,,,dim(phi_M)[4]], times = dim(phi_M)[4]),  dim = dim(phi_M))
      delta_sim <- array( rep(delta[,,dim(delta)[3]], times = dim(delta)[3]),  dim = dim(delta))
      kappa0_sim <- array( rep(kappa0[,dim(kappa0)[2]], times = dim(kappa0)[2]),  dim = dim(kappa0)) 
    }
    
    
    # If we have entered the disturbed pattern, use the disturbed demographic parameters (the ones provided by the user as the function arguments)
    
    #print(paste0("t = ", t))
    #if (t %% 10 == 0) print(paste0("t = ", t))
    
    # At each time-step, redefine variable Nprime
    # If save.res, redefine also larval and settlers numbers
    if (save.res) {
      Nprime  <- array(NA,dim=c(m,n,z))
      L       <- array(NA,dim=c(m,n))
      S       <- array(0,dim=c(m,n))
      
    } else {
      Nprime  <- array(NA,dim=c(m,n,z))
    }
    
    
    ### Survival
    
    # If there is only one age-class, we must force the third dimension. What if only one year?
    if (length(dim(sigma))==2){
      dim(sigma)[3] <- 1
      dim(sigma_sim)[3] <- 1
    } 
    
    if (verbose) print("Apply survival function")
    for (i in 1 : n) {
      for (x in 1 : z) {
        for (k in 1 : m) {
          
          if (save.res){
            Nprime[k,i,x] = surv(sigma_sim[k,i,x,TSD],N[k,i,x])
            
          } else {
            
            if (is.na(sigma_sim[k,i,x,TSD])) warning(paste("sigma: NA detected for k =", k, ", i =", i, ", x =", x, ", t =", TSD))
            
            if (is.na(N[k,i,x,t])) warning(paste("N: NA detected for k =", k, ", i =", i, ", x =", x, ", t =", t))
            
            Nprime[k,i,x] = surv(sigma_sim[k,i,x,TSD],N[k,i,x,t])
          }
        }
      }
    }
    
    
    ## Reproduction
    
    if (verbose) print("Apply reproduction function")
    
    # If there is only one age-class, we must force the third dimension
    if (length(dim(phi_F))==2){
      dim(phi_F)[3] <- 1
      dim(phi_F_sim)[3] <- 1
    } 
    if (length(dim(phi_M))==2){
      dim(phi_M)[3] <- 1
      dim(phi_M_sim)[3] <- 1
    } 
    
    if (migration == "backward") {
      G_F <- array(0,c(l,n))
      G_M <- array(0,c(l,n))
    }
    
    for (i in 1 : n) {
      
      if (save.res) {
        
        if (sum(Nprime[,i,])==0) { # To save computing time
          L[,i] = 0
          next
        } else {
          res.repr <- repr.monoecious.onelocus(Nprime, phi_F_sim[,,,TSD], phi_M_sim[,,,TSD], mu, i, l, m, n, z,
                                               fec.distr_F = fec.distr_F, fec.distr_M = fec.distr_M, migration=migration, verbose=verbose)
          if (migration == "forward") {
            L[,i] <- res.repr
          } else {
            G_F[,i] <- res.repr$G_F
            G_M[,i] <- res.repr$G_M
          }
        }
        
      } else {
        
        if (sum(Nprime[,i,])==0) { # To save computing time
          L[,i,t] = 0
          next
        } else {
          res.repr <- repr.monoecious.onelocus(Nprime,phi_F_sim[,,,TSD], phi_M_sim[,,,TSD],mu, i, l, m, n, z,
                                               fec.distr_F=fec.distr_F, fec.distr_M=fec.distr_M, migration=migration, verbose=verbose)
          if (migration == "forward") {
            L[,i,t] <- res.repr
          } else {
            G_F[,i] <- res.repr$G_F
            G_M[,i] <- res.repr$G_M
            
          }
        }
        
      }
    }
    
    
    ## Propagule dispersal
    if (migration == "forward") {
      if (verbose) print("Apply dispersal function")
      for (i in 1 : n) {
        for (k in 1 : m) {
          if (save.res) {
            y = disp(L[k,i],delta_sim[,i,TSD])
            S[k,] <- S[k,] + y[1:n]       
          } else {
            y = disp(L[k,i,t],delta_sim[,i,TSD])
            S[k,,t] <- S[k,,t] + y[1:n]
          }
        }
      }
    }
    
    
    
    # Save results if save.res=T
    if (save.res){
      if ((t+1) %in% save.res.T) {
        file.name <- paste("N",(t+1),".RData",sep="")
        save(N,file=paste(dir.res.name,file.name,sep="/"))
      }
    }
    
    
    # Recruitment
    if (t == length(save.res.T)) break # otherwise attempts to write on T_total + 1
    
    ###The main difference between the original function and this new one is that the recruitment is now done after the seeds are stored in the seedbank, and are then given a chance to germinate.
    # The recruitment cohort is here the seeds that successfully germinated.
    # I added two steps: storing the propagules in the seedbank, and then giving a germination try.
    
    # Seedbank storage
    # store the propagules in the seedbank for actual generation 't'.
    if (save.res) {
      SB[t+1,,] <- S[,]       
    } else {
      SB[t+1,,] <- S[,,t]
    }
    
        # substep: Germination: returns a vector (updated_seedbank, germination_list)
    
    updated_seedbank <- array(NA, dim = c(T_max,m,n))
    germination_list <- array(NA, dim = c(m,n))
    
    for (i in 1:n){
      germination_result <-germination(seedbank = SB,
                                       deme = i,
                                       gamma,
                                       time_since_disturbance = TSD,
                                       seed_life_expectancy,
                                       t+1,
                                       degeneration_type)
      
      updated_seedbank <- germination_result$seedbank
      germination_list[,i] <- germination_result$germination_list
    }
    
    
    SB <- updated_seedbank
    
    if (save.res){
      S[,] <- germination_list
    } else {
      S[,,t] <- germination_list
    }
    
    
    
    
    # Recruitment with backward migration
    if (migration == "backward") {
      for (j in 1 : n) {
        res.recr <- recr.backward.migration.onelocus(G_F, G_M, migr, l, m, n, z, j.local.deme=j, kappa0_sim[j,TSD+1],
                                                     sexuality="monoecious", TSD=TSD+1)
        if (save.res) {
          N[,j,1] <- res.recr # 1 because only 1 age class is supported now
        } else {
          N[,j,1,t+1] <- res.recr # idem
        }
        
      }
      next
    }
    
    
    if (verbose) print("Apply recruitment function")
    for (i in 1 : n) {
      if (save.res) {
        N[,i,] <- recr(N = array(Nprime[,i,], dim=c(m,z)),
                       S = array(S[,i], dim=c(m,1)),
                       m = m,
                       z = z,
                       kappa0 = kappa0_sim[i,TSD+1], # Pay attention here: using carrying capacity of next generation
                       recr.dd = recr.dd,
                       sexuality="monoecious")
      } else {
        N[,i,,t+1] <- recr(N = array(Nprime[,i,], dim=c(m,z)),
                           S = array(S[,i,t], dim=c(m,1)),
                           m = m,
                           z = z,
                           kappa0 = kappa0_sim[i,TSD+1],
                           recr.dd = recr.dd,
                           sexuality = "monoecious")
      }
    }
    
    
    ## Disturbance
    
      # Check if the generation has to be disturbed or not. If yes (i.e. the value in the vector is not 0), apply the disturbance.
    
    if (simulation_vector[t]!=0){
      
      #fire_index <- simulation_vector[t]
      #print(paste0("Fire number ", fire_index))

      
      if (save.res){
        #charge the last created file and retrieve N
        latest_dir <- get_latest_folder(parent)
        file_to_load <- paste0(latest_dir, "/N", t, ".RData")
        
        #only print an error if the file does not exist
        if (!file.exists(file_to_load)) {
          stop("The specified file does not exist: ", file_to_load)
        }
        load(file_to_load)
        
        #apply fire
        N_prime <- N
        #SB_prime <- SB
        
        N_burnt <- fire_save_res(N_prime, psi)
        #SB_prime <- fire_seedbank(SB_prime, m, n, z, psi_seedbank) #if effect of fire on seedbank was to be implemented
        
        N <- N_burnt
        #SB <- SB_prime
        
        #save the data
        file.name <- paste("N",(t),".RData",sep="")
        save(N,file=paste(dir.res.name,file.name,sep="/"))
        
        
        
      } else { #if save.res==F
        
        N_prime <- N
        #SB_prime <- SB
        
        #apply fire
        N_burnt <- fire(N_prime, t, psi)
        #SB_prime <- fire_seedbank(SB_prime, m, n, z, psi_seedbank) #if effect of fire on seedbank was to be implemented
        
        N <- N_burnt
        #SB <- SB_prime
      }
      
      #Resetting the time since disturbance
      TSD <- 0
      
      #and if it is the first disturbance occurring, sets the demographic parameters to the disturbed ones
      if (past_disturbance == FALSE){
        sigma_sim <- sigma
        phi_F_sim <- phi_F
        phi_M_sim <- phi_M
        delta_sim <- delta
        kappa0_sim <- kappa0_sim
        
        past_disturbance <- TRUE
      }
      
    }
    
    
  }
  #End of temporal loop.
  #print("...done")
  
  if (save.res){
    file.name <- paste("seedbank",".RData",sep="")
    save(SB,file=paste(dir.res.name,file.name,sep="/"))
    
  } else { #if save.res==F
    return(list(N=N,SB=SB))
  }
  
}