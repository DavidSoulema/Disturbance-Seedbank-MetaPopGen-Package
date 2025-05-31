#############################################

### 3. DISTURBANCE PACKAGE

#############################################

# MetaPopGen does not allow to simulate disturbance during the simulation.
# Disturbances kill individuals and modify environmental characteristics.
# This subpackage adds a function fire that kills individuals with a certain probability.
# This probability vary across age-classes but not across demes (same effect).
# To model the change in environmental characteristics post-fire, I chose to modify the parameters (phi, sigma, delta, kappa0, gamma).


## The disturbance here is fire. Adapt the fire function below to your disturbance if necessary. 




# Function get_latest_folder: fetch the latest folder created. Used when save.res==T

get_latest_folder <- function(path = ".") {
  folders <- list.dirs(path, full.names = TRUE, recursive = FALSE)
  if (length(folders) == 0) return(NULL) # No folder found
  
  latest_folder <- folders[which.max(file.info(folders)$mtime)]
  return(latest_folder)
}

#------------------



# Function fire: takes a population N[m,n,z,t] and returns the population after having killed some individuals according to the parameter psi
  # Does not affect the seedbank. Could copy/paste this function to create a function fire_seedbank acting on SB instead of N.

fire <- function(N, t, psi){
  
  m <- dim(N)[1]
  n <- dim(N)[2]
  z <- dim(N)[3]
  
  N_burnt <- N

  for (i in 1:n){ #won't change anything : works the same way on all demes.
    for (x in 1:z){
      for (k in 1:m){
        individuals <- N_burnt[k,i,x,t]
        
        if (is.na(individuals)) {
          warning(paste("NA detected for k =", k, ", i =", i, ", x =", x, ", t =", t))
        }
        
        if (individuals!=0){
          deaths <- sum(rbinom(individuals, 1, psi[1,x]))
          survivals <- individuals - deaths
          #optional: print the number of casualties
          #message(deaths, "individuals have been killed by fire")
          
          #now update the surviving individuals, leave the loop and update the population after the fire
          N_burnt[k,i,x,t] <- survivals
        }
      }
    }
  }
  return(N_burnt)
}

#------------------




# Function fire_save_res
 # same function if save.res==T (will not manipulate arrays of the same size because no time dimension)

fire_save_res <- function(N, psi){
  
  m <- dim(N)[1]
  n <- dim(N)[2]
  z <- dim(N)[3]
  
  N_burnt <- N
  
  for (i in n){ #won't change anything: works the same way on all demes.
    for (x in z){
      for (k in m){
        individuals <- N_burnt[k,i,x]
        
        if (is.na(individuals)) {
          warning(paste("NA detected for k =", k, ", i =", i, ", x =", x, ", t =", t))
        }
        
        if (individuals!=0){
          deaths <- sum(rbinom(individuals, 1, psi[1,x]))
          survivals <- individuals - deaths
          #optional: print the number of casualties
          #message(deaths, "individuals have been killed by fire")
          
          #now update the surviving individuals, leave the loop and update the population after the fire
          N_burnt[k,i,x] <- survivals
        }
      }
    }
  }
  return(N_burnt)
}

#------------------





# Function to create a vector indicating the generations when the disturbance should return

get_disturbance_vector <- function(T_max,
                                   generation_time,
                                   stable_state_time = 0,
                                   disturbance_return_interval,
                                   final_time_since_disturbance = NA) {
  
  SST <- stable_state_time
  FTSD <- final_time_since_disturbance
  
  if (disturbance_return_interval==0){ # in case we don't want any disturbance ("disturbance_return_interval = 0")
    T_total <- SST + T_max
    return(rep(0, T_total))
  }
  
  # compute the disturbance_return_interval in generations instead of years (FRI)
  DRI <- floor(disturbance_return_interval / generation_time)
  
  
  if (DRI==0){
    stop("Warning: Fire return interval shorter than generation time. ")
  }
  
  
  if (DRI >= T_max){ #if DRI is bigger than the simulation time, then only one disturbance occurs.
    
    T_total <- SST + T_max
    disturbance_vector <- rep(0, T_total)
    disturbance_vector[SST+1] <- 1
    
    } else {
      
      # compute the total number of generations separated by fires
      disturbance_number <- T_max %/% DRI            #quotient of Euclidian division
      remaining_generations <- T_max %% DRI   #remainder of Euclidian division
      
      # compute the total length of the simulation
      
      if (is.na(FTSD)) {
        T_total <- SST + T_max
        
      } else {
        # to get the same number of generations since last disturbance at the end of the simulation - and make appropriate comparisons - we need to take out the remainder and then add FTSF
        T_total <- SST + T_max - remaining_generations + FTSD + 1
      }
      
      disturbance_vector <- rep(0, T_total)
      
      # Replace the 0s when the disturbance returns by the index of the disturbance (1 if it is the 1st time, 2 if it is the 2nd, etc...)
      for (index in 1:(disturbance_number+1)){
        disturbed_generation <- SST + (index-1) * DRI + 1
        disturbance_vector[disturbed_generation] <- index 
      }
  
    }
  
  return(disturbance_vector)
}
