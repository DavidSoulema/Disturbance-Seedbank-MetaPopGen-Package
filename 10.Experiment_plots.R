#############################################

### 9. RETRIEVE THE EXPERIMENT PLOTS

#############################################


# If you've run the experiment already and just want to retrieve different plots without having to re-run the simulation, just run this code. It follows the same structure as  Experiment.R




##########################################

# Required packages

##########################################

# Required packages
library(data.table)  # For efficient data manipulation
library(ggplot2)     # For advanced plotting
library(scales)      # For formatting axes in ggplot
library(gridExtra)
library(grid)



##########################################

# Initialization

##########################################

# Path to the main folder
base_dir <- "C:\\Users\\david\\Documents\\DAVID\\CS\\DD\\Cours\\Master Thesis\\Coding\\Codes\\05-26\\Experiment\\Experiment"

# DRI values
DRI_values <- c(1:10, 20, 100, 0)


# .csv file to get the numerical He values
table_matrix <- matrix(NA, nrow = 3, ncol = length(DRI_values))
rownames(table_matrix) <- c("mean He", "upper He", "lower He")
colnames(table_matrix) <- as.character(DRI_values)

df_table <- data.frame(DRI = rownames(table_matrix), table_matrix, check.names = FALSE)


# Change to TRUE to save
save.plots.separated = FALSE    #curves of He over time, in separated images
save.plots.summary = FALSE      #curves of He over time, in a single image
save.individual.box = FALSE     #individual boxplot for each simulation
save.global.boxes = FALSE       #all boxplots combined
save.global.curves = FALSE      #all curves combined in a mosaic


#To plot the global.curves in a mosaic (setup)
plots <- list()
n_cols <- 4             #number of column in the mosaic
y_lim = c(0.40, 0.51)   #limit of y-axis




##########################################

# Retrieve the He_all.csv for each Disturbance Return Interval (DRI)

##########################################


# Loop over each DRI value
for (DRI in DRI_values) {
  
  # Generate the disturbance vector
  disturbance_vector <- get_disturbance_vector(T_max,
                                               generation_time,
                                               stable_state_time=0,
                                               disturbance_return_interval=DRI,
                                               final_time_since_disturbance = 1)
  T_total <- length(disturbance_vector)
  
  
  
  DRI_dir <- file.path(base_dir, paste0("Experiment_DRI_", DRI))
  file_path <- file.path(DRI_dir, "He_all.csv")
  
  He_all <- fread(file_path)
  He_all_matrix <- as.matrix(He_all)
  
  
  
  ###-----------------------
  # Generate and save Plots
  ###-----------------------
  
  # Plot He_all.jpg
  # Prepare data: melt matrix to long format for ggplot
  df_plot <- melt(as.data.table(He_all_matrix), id.vars = NULL,
                  variable.name = "Simulation", value.name = "He")
  df_plot[, Generation := rep(1:T_total, Sim)]
  
  # Compute the mean He per generation
  df_plot[, Mean := ave(He, Generation, FUN = mean)]
  
  # Create data frame for disturbance lines
  disturbance_indices <- which(disturbance_vector != 0)
  
  disturbance_df <- if (length(disturbance_indices) > 0) {
    data.frame(x = disturbance_indices,
               type = "Disturbance occurrence")
  } else {
    data.frame(x = numeric(0),
               type = character(0))
  }
  
  # Compute the statistical data
  He_mean <- rowMeans(He_all_matrix)
  He_sd <- apply(He_all_matrix, 1, sd)
  error_margin <- qnorm(0.975) * He_sd / sqrt(Sim)
  lower <- He_mean - error_margin
  upper <- He_mean + error_margin
  
  # Save it in the .csv file
  DRI_col <- as.character(DRI)
  df_table[df_table$DRI == "mean He", DRI_col] <- tail(He_mean, 1)
  df_table[df_table$DRI == "lower He", DRI_col] <- tail(lower, 1)
  df_table[df_table$DRI == "upper He", DRI_col] <- tail(upper, 1)
  
  
  if (save.plots.separated==TRUE){
    
    # Plot  all individual He curves + mean
    
    p_all <- ggplot(df_plot, aes(x = Generation, y = He, group = Simulation)) +
      geom_line(color = "lightgrey", alpha = 0.1, size = 0.3) +  # individual curves
      geom_line(aes(y = Mean), color = "blue", size = 1) +     # mean curve
      geom_vline(data = disturbance_df,
                 aes(xintercept = x, linetype = type, color = type),
                 alpha = 0.4,
                 show.legend = TRUE) +
      scale_color_manual(values = c("Disturbance occurrence" = "indianred2")) +
      scale_linetype_manual(values = c("Disturbance occurrence" = "dotted")) +
      labs(x = "Generations", y = "He", title = paste("All He curves for DRI =", DRI),
           color = "", linetype = "") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, vjust = 1),
        legend.position = "top"
      )
    
    ggsave(filename = file.path(DRI_dir, "He_all.jpg"), plot = p_all, width = 8, height = 6)
    
    
    
    # Plot mean He + 95% CI
    
    df_summary <- data.table(Generation = 1:T_total, Mean = He_mean, Lower = lower, Upper = upper)
    
    
    p_mean <- ggplot(df_summary, aes(x = Generation, y = Mean)) +
      geom_line(color = "blue", linewidth = 1.2) +
      geom_ribbon(data = df_summary, aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "blue") +
      geom_vline(data = disturbance_df,
                 aes(xintercept = x, linetype = type, color = type),
                 alpha = 0.4,
                 show.legend = TRUE) +
      scale_color_manual(values = c("Disturbance occurrence" = "indianred2")) +
      scale_linetype_manual(values = c("Disturbance occurrence" = "dotted")) +
      labs(x = "Generations", y = "He", title = paste("Mean He with 95% CI for DRI =", DRI),
           color = "", linetype = "") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, vjust = 1),
        legend.position = "top"
      )
    
    ggsave(filename = file.path(DRI_dir, "He_mean.jpg"), plot = p_mean, width = 8, height = 6)
  }
  
  
  
  # Plot all He curves + mean + 95% CI
  
  if (save.plots.summary==TRUE){
    
    df_summary <- data.table(Generation = 1:T_total, Mean = He_mean, Lower = lower, Upper = upper)
    
    
    p_curves <- ggplot(df_plot, aes(x = Generation, y = He, group = Simulation)) +
      geom_line(color = "lightgrey", alpha = 0.1, size = 0.3) +
      geom_line(aes(y = Mean), color = "blue", size = 1) +
      geom_ribbon(data = df_summary, aes(x = Generation, ymin = Lower, ymax = Upper), 
                  inherit.aes = FALSE, alpha = 0.3, fill = "blue") +
      geom_vline(data = disturbance_df,
                 aes(xintercept = x, linetype = type, color = type),
                 alpha = 0.4,
                 show.legend = TRUE) +
      scale_color_manual(values = c("Disturbance occurrence" = "indianred2")) +
      scale_linetype_manual(values = c("Disturbance occurrence" = "dotted")) +
      labs(x = "Generations", y = "He", title = paste("DRI =", DRI),
           color = "", linetype = "") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, vjust = 1),
        legend.position = "top"
      )
    
    ggsave(filename = file.path(DRI_dir, "He_curves.jpg"), plot = p_curves, width = 8, height = 6)
    
  }
  
  
  
  # Plot all the curves for each DRI in a mosaic
  
  if(save.global.curves==TRUE){
    
    p <- ggplot(df_plot, aes(x = Generation)) +
      geom_line(aes(y = He, group = Simulation), color = "lightgrey", alpha = 0.1, linewidth = 0.3) +
      geom_ribbon(data = df_summary, aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.3) +
      geom_line(aes(y = Mean), color = "blue", linewidth = 1) +
      geom_vline(data = disturbance_df,
                 aes(xintercept = x, linetype = type, color = type),
                 alpha = 0.4, show.legend = FALSE) +
      scale_color_manual(values = c("Disturbance occurrence" = "red")) +
      scale_linetype_manual(values = c("Disturbance occurrence" = "dotted")) +
      labs(title = paste("DRI =", DRI), x = "Generations", y = "He") +
      coord_cartesian(ylim = y_lim) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, vjust = 1),
            legend.position = "none")
    
    
    plots <- append(plots,list(p))
  }
  
  
  
  # Plot He_box.jpg (boxplot of last generation)
  
  if (save.individual.box==TRUE){
    
    last_gen_values <- He_all_matrix[T_total, ]
    df_box <- data.table(DRI = as.factor(DRI), He = last_gen_values)
    p_box <- ggplot(df_box, aes(x = DRI, y = He)) +
      geom_boxplot(fill = "skyblue") +
      labs(x = "", y = "He", title = paste("He Boxplot at last generation for DRI =", DRI)) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(filename = file.path(DRI_dir, "He_box.jpg"), plot = p_box, width = 6, height = 4)
    
  }
  
  # END OF DRI LOOP : back to Global Experiment repository
}



##########################################

# Save the global plots (curves and box plot)

##########################################

# He statistics for each DRI

csv_file <- file.path(base_dir, "He_statistics_final.csv")
write.csv(df_table, csv_file, row.names = FALSE, na = "")



# global He_box.jpg for all DRI

if (save.global.boxes==TRUE){
  
  box_data <- rbindlist(lapply(DRI_values, function(DRI) {
    df <- fread(file.path(base_dir, paste0("Experiment_DRI_", DRI), "He_all.csv"))
    T_total_local <- nrow(df)
    He_vals <- df[T_total_local, ]
    melt(data.table(DRI = as.factor(DRI), t(He_vals)), 
         id.vars = "DRI", 
         measure.vars = patterns("^V"), 
         value.name = "He")
  }))
  
  p_final_box <- ggplot(box_data, aes(x = DRI, y = He)) +
    geom_boxplot(fill = "lightgreen") +
    labs(x = "Disturbance Return Interval (DRI)", y = "He") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = file.path(base_dir, "He_boxes.jpg"), plot = p_final_box, width = 10, height = 6)
}


# global He_curves.jpg for all DRI

if (save.global.curves==TRUE){
  
  plots_grid <- arrangeGrob(grobs = plots, ncol = n_cols)
  
  # Create legend plot for curve descriptions
  legend_grid <- ggplot() +
    geom_line(aes(x = 1:4, y = 1:4, color = "Individual He")) +
    geom_line(aes(x = 1:4, y = 1:4, color = "Mean He"), size = 1) +
    geom_ribbon(aes(x = 1:4, ymin = 0.8, ymax = 1.2, fill = "95% CI"), alpha = 0.3) +
    geom_vline(aes(xintercept = 2, linetype = "Disturbance occurrence"), color = "red", alpha = 1) +
    scale_color_manual(values = c("Individual He" = "lightgrey", "Mean He" = "blue")) +
    scale_fill_manual(values = c("95% CI" = "blue")) +
    scale_linetype_manual(values = c("Disturbance occurrence" = "dotted")) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 2)),
      fill = guide_legend(order = 2),
      linetype = guide_legend(order = 3)
    ) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  
  
  final_plot <- arrangeGrob(
    plots_grid,
    legend_grid,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - unit(1, "lines"), unit(1, "lines"))
  )
  
  ggsave(
    filename = file.path(base_dir, "He_grid3.jpg"),
    plot = final_plot,
    width = 12,  
    height = 10, 
    units = "in",
    dpi = 300
  )
  
  
}
