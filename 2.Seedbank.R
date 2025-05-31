#############################################

### 2. SEEDBANK PACKAGE

#############################################

# MetaPopGen uses an age-class system that does not allow to model the seedbank.
# The seedbank should store the offsprings issued by mature adults through sexual reproduction until the seeds germinate or die.
# I represent the seedbank as a matrix: one line for each generation time, and one column for each genotype. Values of the cells represent the number of individuals for each generation and genotype.
# There is one seedbank per deme.
# Looking at Marco Andrello's graph to explain the steps of his functions, the seedbank should intervene between 'Propagule dispersal' and 'Recruitment'.






#Degeneration function: generates a probability of germination according to the seed age and the model of degeneration

degeneration <- function(best_proba_germination, seed_life_expectancy, seed_age, degeneration_type = "linear"){
  
  #Choice of degeneration model: linear or exponential
  if (degeneration_type=="linear"){
    proba <- max( 0, best_proba_germination * (1 - (seed_age-1) / seed_life_expectancy)) #is equal to best_proba_germination when the seed is fresh new
    
  } else {
    if (degeneration_type=="exponential"){
      proba <- best_proba_germination * exp( -2 * ln(10) / seed_life_expectancy * (seed_age-1)) #is equal to best_proba_germination when the seed is fresh new. When the seed will reach its life expectancy, the probability will be of 1% of the ideal germination probability
    } else {
      print("Specify degeneration model")
    }
  }
  return(proba)
}

#------------------




#Death function: use it to save calculation time. Will kill all seeds which germination probability is below 1%

death <- function(seeds, proba){
  if (proba <= 0.01){
    seeds <- 0
  }
  return (seeds)
}

#------------------




#Germination function: look at the seedbank matrix and gives all seeds a try to germinate. Extracts the germinated ones into another matrix that can be turned into the first age class z=1, and updates the seed bank matrix.
 #Returns 2 elements in a list: (seedbank, germination_list)


germination <- function(seedbank,
                        deme,
                        gamma,
                        time_since_disturbance,
                        seed_life_expectancy,
                        generation,
                        degeneration_type="linear"){
  
  #generation in (1:T_max)
  #seedbank dimensions: (T_max, m, n)
  
  matrix <- seedbank
  i <- deme
  germination_list <- array(0, dim=c(m)) #list of the seeds that will germinate  and become the z=1
  
  #Let's look at the seeds of each age and genotype:
  
  for (t in 1:generation){ #We don't want to look at the whole matrix to save calculation time
    for (k in 1:m){
      seeds <- matrix[t,k,i]
      
      if (seeds!=0){
        
        #Now calculate the age of each seed, the associated degeneration and germination probability ("proba"). If "proba" is too low, kill the seeds
        age <- generation - t #as we go forward in time (generation), the seeds of the first rows (small t) grow older
        
        best_proba_germination <- gamma[1,time_since_disturbance]
        
        proba <- degeneration( best_proba_germination, seed_life_expectancy, age, degeneration_type )
        
        if (is_death==TRUE){
          seeds <- death( seeds, proba ) #to save calculation time in the next iteration
        }
        #Now give the surviving seeds a chance to germinate (using rbinom)
        germinated <- sum(rbinom(seeds, 1, proba))
        non_germinated <- seeds - germinated
        
        germination_list[k] <- germination_list[k] + germinated
        matrix[t,k,i] <- non_germinated
       }
       
    }
  }
  #Update the seedbank with the dead seeds and the germinated ones
  seedbank <- matrix
  return(list(seedbank=seedbank,germination_list=germination_list))
}
