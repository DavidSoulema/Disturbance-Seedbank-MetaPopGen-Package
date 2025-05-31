#############################################

### 7. OUTPUT FUNCTIONS

#############################################

# To retrieve allele number (Na) and Heterozygosity (He) from population (N)



##########################################

# General functions

##########################################

#from MetaPopGen package

freq.all = function(N){ #N must be a vector.
  m <- length(N)              # Number of genotypes
  l <- (sqrt(1+8*m)-1)/2      # Number of alleles
  p <- array(0,dim=l)
  ntot=sum(N)
  pos.gen <- genotype.index(l)
  for (i in 1 : l) {
    p[i] <- 0.5 * ( sum(N[pos.gen[i,]]) + N[pos.gen[i,i]] ) / ntot
  }
  return(p)
}

#-----------



#new function for Allele Frequency

get_Fa <- function(N_deme_ageclass){ #must be a 1-dimension array.
  
  N_vector <- as.vector(N_deme_ageclass)
  p <- freq.all(N_vector)
  return(p)
}

#-----------



##########################################

# Allele frequencies Fa

##########################################

##To apply on N[m,n,z] (no t!)


#for each deme:

get_Fa_deme <- function(N){
  m <- dim(N)[1]
  n <- dim(N)[2]
  
  l <- (sqrt(1+8*m)-1)/2
  
    #summing N on the age-classes to focus on the variation between demes
  N_all_age <- apply(N, MARGIN = c(1,2), FUN = sum)
  
    #creating a matrix of alleles frequencies among demes
  Fa_deme <- apply(N_all_age, MARGIN = 2, FUN = get_Fa)
  dimnames(Fa_deme)=list(alleles=c(1:l), deme=c(1:n))
  
  return(Fa_deme)
}

#-----------




#for each age-class:

get_Fa_age <- function(N){
  m <- dim(N)[1]
  z <- dim(N)[3]
  
  l <- (sqrt(1+8*m)-1)/2
  
    #summing N on the demes to focus on the variation between age-classes
  N_all_deme <- apply(N, MARGIN = c(1,3), FUN = sum)

    #creating a matrix of alleles frequencies among age-classes
  Fa_age <- apply(N_all_deme, MARGIN = 2, FUN = get_Fa)
  dimnames(Fa_age)=list(alleles=c(1:l), age=c(1:z))
  
  return(Fa_age)
}

#-----------




#for the global population:

get_Fa_global <- function(N){
  m <- dim(N)[1]

  l <- (sqrt(1+8*m)-1)/2
  
    #summing N on the demes and the age-classes to focus on the variation among the whole population
  N_all_deme_all_age <- apply(N, MARGIN = 1, FUN = sum)

    #creating a matrix of global alleles frequencies
  Fa_global <- get_Fa(N_all_deme_all_age)
  names(Fa_global) <- paste0("allele_", seq_along(Fa_global))
  
  return(Fa_global)
}





##########################################

# Allele richness Na and heterozygosity He

##########################################

# Allele richness

allele_richness <- function(Fa){#Fa must be a vector / a 1-dimension array
  
  Fa_vector <- as.vector(Fa)
  Na <- length( which( Fa_vector>0 ) )
  
  return(Na) #a scalar
}

#-----------



get_Na_deme <- function(N){ #must be N[m,n,z]
  
  Fa_deme <- get_Fa_deme(N)
  Na_deme <- apply(Fa_deme, MARGIN = 2, FUN = allele_richness)
  names(Na_deme) <- paste0("deme_", seq_along(Na_deme))
  
  return(Na_deme) #a vector
}

#-----------



get_Na_age <- function(N){ #must be N[m,n,z]
  
  Fa_age <- get_Fa_age(N)
  Na_age <- apply(Fa_age, MARGIN = 2, FUN = allele_richness)
  names(Na_age) <- paste0("age class_", seq_along(Na_age))
  
  return(Na_age) #a vector
}

#-----------



get_Na_global <- function(N){ #must be N[m,n,z]
  
  Fa_global <- get_Fa_global(N)
  Na_global <- allele_richness(Fa_global)
  
  return(Na_global) #a scalar
}
#-----------



# Heterozygosity

heterozygosity <- function(Fa){ #Fa must be a vector / a 1-dimension array
  
  Fa_vector <- as.vector(Fa)
  He <- 1 - sum( Fa_vector^2 )
  
  return(He) #a scalar
}

#-----------



get_He_deme <- function(N){ #must be N[m,n,z]
  
  Fa_deme <- get_Fa_deme(N)
  He_deme <- apply(Fa_deme, MARGIN = 2, FUN = heterozygosity)
  names(He_deme) <- paste0("deme_", seq_along(He_deme))
  
  return(He_deme) #a vector
}

#-----------



get_He_age <- function(N){ #must be N[m,n,z]
  
  Fa_age <- get_Fa_age(N)
  He_age <- apply(Fa_age, MARGIN = 2, FUN = heterozygosity)
  names(He_age) <- paste0("age class_", seq_along(He_age))
  
  return(He_age) #a vector
}

#-----------



get_He_global <- function(N){ #must be N[m,n,z]
  
  Fa_global <- get_Fa_global(N)
  He_global <- heterozygosity(Fa_global)
  
  return(He_global) #a scalar
}







##########################################

# Output over time or different simulations

##########################################

#to study the evolution of Na and He over time (or over different simulation results)


# Allele diversity

  #Grouped by deme

get_Na_deme_over_time <- function(N) { #must be N[m,n,z,T_max]
  T_max <- dim(N)[4]
  
  Na_over_time <- sapply(1:T_max, function(t) { #extraction of all the N (dimensions m,n,z)
    get_Na_deme(N[,,,t, drop = FALSE])
  })
  
  return(Na_over_time) #a list of vectors (length = T_max)
}

#-----------



  #Grouped by age class

get_Na_age_over_time <- function(N) { #must be N[m,n,z,T_max]
  T_max <- dim(N)[4]
  
  Na_over_time <- sapply(1:T_max, function(t) { #extraction of all the N (dimensions m,n,z)
    get_Na_age(N[,,,t, drop = FALSE])
  })
  
  return(Na_over_time) #a list of vectors (length = T_max)
}

#-----------



  #Among the whole population

get_Na_global_over_time <- function(N) { #must be N[m,n,z,T_max]
  T_max <- dim(N)[4]
  
  Na_over_time <- sapply(1:T_max, function(t) { #extraction of all the N (dimensions m,n,z)
    get_Na_global(N[,,,t, drop = FALSE])
  })
  
  return(Na_over_time) #a list of scalar (length = T_max)
}

#-----------



# Heterozygosity

  #Grouped by deme

get_He_deme_over_time <- function(N) { #must be N[m,n,z,T_max]
  T_max <- dim(N)[4]
  
  He_over_time <- sapply(1:T_max, function(t) { #extraction of all the N (dimensions m,n,z)
    get_He_deme(N[,,,t, drop = FALSE])
  })
  
  return(He_over_time) #a list of vectors (length = T_max)
}

#-----------




  #Grouped by age class

get_He_age_over_time <- function(N) { #must be N[m,n,z,T_max]
  T_max <- dim(N)[4]
  
  He_over_time <- sapply(1:T_max, function(t) { #extraction of all the N (dimensions m,n,z)
    get_He_age(N[,,,t, drop = FALSE])
  })
  
  return(He_over_time) #a list of vectors (length = T_max)
}

#-----------




  #Among whole population

get_He_global_over_time <- function(N) { #must be N[m,n,z,T_max]
  T_max <- dim(N)[4]
  
  He_over_time <- sapply(1:T_max, function(t) { #extraction of all the N (dimensions m,n,z)
    get_He_global(N[,,,t, drop = FALSE])
  })
  
  return(He_over_time) #a list of scalar (length = T_max)
}

