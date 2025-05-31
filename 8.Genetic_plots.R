##########################################

# 8. Genetic Plots

##########################################


# To plot He, Na, N and SB over time.
  # Possibility to plot them for each deme, age-class, or global (whole population).
  # You might want to complete this code to plot for each deme and age-class with a simple copy-paste of the last lines.


result_sim <- N_res #(m,n,z,T_max)
SB_sim <- SB_res #(T_max,m,n)

# Fetch the outputs from the simulation 
He_global_over_time <- get_He_global_over_time(result_sim)
He_deme_over_time <- get_He_deme_over_time(result_sim)
He_age_over_time <- get_He_age_over_time(result_sim)

Na_global_over_time <- get_Na_global_over_time(result_sim)
Na_deme_over_time <- get_Na_deme_over_time(result_sim)
Na_age_over_time <- get_Na_age_over_time(result_sim)

N_global_over_time <- apply(result_sim, MARGIN = 4, FUN = sum)
SB_global_over_time <- apply(SB_sim, MARGIN = 1, FUN = sum)



# Plot them over time
generations <- 1:T_total

#create the plots
# He global
plot(generations, He_global_over_time, type = "n", xlab = "Generations", ylab = "Heterozygosity",
     ylim = c(0, 0.5), main = paste0("Heterozygosity over Generations (l = ", l, ")"), lwd = 2)
abline(v = generations[disturbance_vector != 0], col = "red", lty = 2)
lines(generations, He_global_over_time, type = "l", lwd = 2)

# Na global
plot(generations, Na_global_over_time, type = "n", xlab = "Generations", ylab = "Allele diversity", 
     main = paste0("Allele diveristy over Generations (l = ",l,")"), lwd = 2) 
abline(v = generations[disturbance_vector != 0], col = "red", lty = 2)
lines(generations, Na_global_over_time, type = "l", lwd = 2)

# N global
plot(generations, N_global_over_time, type = "n", xlab = "Generations", ylab = "Total population", ylim = c(0, 1000),
     main = paste0("Total population over Generations (l = ",l,")"), lwd = 2)
abline(v = generations[disturbance_vector != 0], col = "red", lty = 2)
lines(generations, N_global_over_time, type = "l", lwd = 2)

# SB global
plot(generations, SB_global_over_time, type = "n", xlab = "Generations", ylab = "Total seedbank", ylim = c(0, 1000),
     main = paste0("Total seedbank over Generations (l = ",l,")"), lwd = 2) 
abline(v = generations[disturbance_vector != 0], col = "red", lty = 2)
lines(generations, SB_global_over_time, type = "l", lwd = 2)



