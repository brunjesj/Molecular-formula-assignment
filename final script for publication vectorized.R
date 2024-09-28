library(data.table)
library(mgcv)
library(dplyr)
library(quantreg)
library(igraph)
library(plotly)
# Define the directory containing the .csv files
directory_path_blanks <- "C:/Users/jonas/OneDrive/Data/FTMS/All measured samples/Instrument Blanks"
directory_path_samples <- "C:/Users/jonas/OneDrive/Data/FTMS/All measured samples/Water samples"

# Get a list of all .csv files in the directories
csv_files_blanks <- list.files(path = directory_path_blanks, pattern = "\\.csv$", full.names = TRUE)
csv_files_samples <- list.files(path = directory_path_samples, pattern = "\\.csv$", full.names = TRUE)

# Read each .csv file into a list, using the file name (without extension) as the list name
blank_list <- lapply(csv_files_blanks, fread)
sample_list <- lapply(csv_files_samples, fread)

names(blank_list) <- tools::file_path_sans_ext(basename(csv_files_blanks))
names(sample_list) <- tools::file_path_sans_ext(basename(csv_files_samples))

#remove outliers of resolution vs mass. From Merder et al, 2020

ResPow_outlier <- function(dataset) {
  # Ensure dataset is a data.table
  setDT(dataset)
  
  # Fit quantile regression model
  d <- rq(log(ResPow) ~ log(m.z), data = dataset, method = "fn")
  
  # Compute density of residuals
  j <- density(d$residuals)
  
  # Use data.table to process the residuals' density
  ex <- data.table(x = j$x[-1], y = diff(j$y))
  ex <- ex[x > 0]
  ex <- ex[which.min(y):.N]
  e <- ex[y > 0, x[1]]
  
  # Handle case where no threshold is found
  if (is.na(e)) e <- Inf
  
  # Filter dataset based on residuals
  dataset <- dataset[d$residuals <= abs(e)]
  
  return(dataset)
}

MDL_master <- lapply(blank_list, function(blank){
  colnames(blank)[1] <- "m.z"
  
  
  #following Riedel 2014
  blank[, nominal_mass := floor(m.z + 0.5)]
  
  # Compute MDL for each nominal mass and return only m.z and MDL
  MDL_master <- blank[, .(
    MDL = qt(0.99, df = .N - 1) * sd(I) + mean(I)
  ), by = nominal_mass]
  
  setnames(MDL_master, "nominal_mass", "m.z")
  MDL_master
})

# Combine all data.tables into one
MDL_master <- rbindlist(MDL_master)

# Compute mean and sd of MDL for each m.z
MDL_master <- MDL_master[, .(
  MDL = mean(MDL),
  sd = sd(MDL)
), by = m.z]



#########
process_samples_MDL <- function(sample, sample_name, MDL_master) {
  # Rename column efficiently
  setnames(sample, "m/z", "m.z", skip_absent = TRUE)
  
  # Calculate nominal_mz without modifying original sample
  sample[, nominal_mz := floor(m.z + 0.5)]
  
  # Merge using `on` to avoid setting keys
  merged_data <- merge(sample, MDL_master, by.x = "nominal_mz", by.y = "m.z", all.x = TRUE)
  
  # Filter and process in one step
  data <- merged_data[I >= MDL, .(
    m.z = m.z, 
    I = I, 
    MDL = MDL, 
    ResPow = Res., 
    index = sample_name
  )]
  
  # Remove outliers in mass resolution
  data <- ResPow_outlier(data)
  
  return(data)
}



# Apply function across sample list
MDL_water <- lapply(seq_along(sample_list), function(i) {
  print(i)
  sample <- sample_list[[i]]
  process_samples_MDL(sample = sample, sample_name = names(sample_list)[i], MDL_master = MDL_master)
  
})

names(MDL_water)<- names(sample_list)

#clean up
rm(list = setdiff(ls(), c( "MDL_water")))
gc()

#################################
#formula processing for calibration
generate_formulas <- function(C_range = 1:50, H_range = 2:120, O_range = 0:50,
                              N_range = 0:4, S_range = 0:2, P_range = 0:1,
                              mass_range = c(50, 1000), HC = c( 0.2, 3), OC = c(0, 1.2), insert_isotopes = T) {
  
  # Define atomic masses with increased precision
  mass_C  <- 12.00000000
  mass_H  <- 1.00782503
  mass_O  <- 15.99491463
  mass_N  <- 14.00307400
  mass_S  <- 31.97207117
  mass_P  <- 30.97376203
  
  # Isotopic masses with increased precision
  mass_C13 <- 13.00335540
  mass_O18 <- 17.99915960
  mass_S34 <- 33.96786690
  
  # Generate combinations of C, N, S, P
  dt <- CJ(C = C_range, N = N_range, S = S_range, P = P_range)
  
  # Apply initial constraints using vectorized filtering
  dt <- dt[C > 0 & (N / C) <= 1.3 & (S / C) <= 0.8 & (P / C) <= 0.3]
  dt <- dt[!(N > 1 & (S + P > 0)) & !(N != 0 & S != 0 & P != 0) & !(S + P > 2) & !(N == 1 & (S > 1 | P > 1))]
  
  # Precompute H and O valid ranges based on C
  dt[, `:=`(H_min = pmax(ceiling(HC[1] * C), H_range[1]),
            H_max = pmin(floor(HC[2] * C), max(H_range)),
            O_max = pmin(floor(OC[2] * C), max(O_range)))]
  
  # Expand for H and O combinations for each formula
  results <- dt[, {
    if (H_min <= H_max) {
      H_values <- H_min:H_max
      O_values <- 0:O_max
      H_O_dt <- CJ(H = H_values, O = O_values)
      H_O_dt <- H_O_dt[(H / C) >= HC[1] & (H / C) <= HC[2] & (O / C) <= OC[2] & (O / C) >= OC[1]]
    } else {
      H_O_dt <- NULL
    }
    H_O_dt
  }, by = .(C, N, S, P)]
  
  # Vectorized DBE and valence calculations
  results[, `:=`(DBE = 1 + (2 * C - H + N + P) / 2,
                 valence = C * 4 + H * 1 + O * 2 + N * 3 + S * 2 + P * 5)]
  
  # Filter valid DBE and valence
  results <- results[DBE >= 0 & DBE %% 1 == 0 & valence %% 2 == 0]
  
  # Remove valence and DBE columns
  results[, `:=`(valence = NULL, DBE = NULL)]
  
  # Vectorized mass calculation
  results[, mz := C * mass_C + H * mass_H + O * mass_O + N * mass_N + S * mass_S + P * mass_P]
  results <- results[mz >= mass_range[1] & mz <= mass_range[2]]
  
  
  
  print("monoisotopic formulas generated")
  if(insert_isotopes){
    
    # Precompute maximum isotope counts for all rows
    results[, `:=`(max_O18 = pmin(1, O),
                   max_C13 = pmin(2, C),
                   max_S34 = pmin(1, S))]
    
    # Generate all isotope combinations without unnecessary replication
    isotope_combinations <- function(formulas_dt) {
      # Generate combinations of isotopes in a more efficient way
      isotope_grid <- CJ(O_18 = 0:max(formulas_dt$max_O18),
                         C_13 = 0:max(formulas_dt$max_C13),
                         S_34 = 0:max(formulas_dt$max_S34))
      
      isotope_grid <- isotope_grid[
        (O_18 + S_34 + C_13== 1) | (C_13 ==2 & O_18 ==0 & S_34==0)
      ]
      
      # Expand the original formulas for each isotope combination
      expanded_formulas <- formulas_dt[rep(1:.N, each = nrow(isotope_grid))]
      
      # Now add the isotope combinations to the expanded rows
      expanded_isotopes <- isotope_grid[rep(1:.N, times = formulas_dt[, .N])]
      
      # Combine the two tables
      combined_formulas <- cbind(expanded_formulas, expanded_isotopes)
      
      # Filter by isotope presence and constraints
      combined_formulas <- combined_formulas[
        (C_13 >0 & O_18 ==0 & S_34==0) | 
          ((S_34 > 0 & S > 0) | (O_18 > 0 & O > 0))
      ]
      
      # Adjust for isotopes and compute mz
      combined_formulas[, mz :=  C * mass_C + H * mass_H + O * mass_O + N * mass_N + S * mass_S + P * mass_P +
                          C_13 * (mass_C13-mass_C) + O_18 * (mass_O18-mass_O) + S_34 * (mass_S34-mass_S)]
      
      # Filter by mass range
      combined_formulas <- combined_formulas[mz >= mass_range[1] & mz <= mass_range[2]]
      
      
      return(combined_formulas)
    }
    
    # Generate isotope formulas
    results_isotopes <- isotope_combinations(results)
    
    # Bind isotope formulas with monoisotopic formulas
    results <- rbind(results, results_isotopes, use.names = TRUE, fill = TRUE)
  }else{
    new_columns <- c("O_18", "C_13", "S_34")
    results[, (new_columns) := 0]
  }
  
  # Final cleanup and ordering
  results <- results[, .(C, H, O, N, S, P, C_13, O_18, S_34, calculated_m.z=mz)]
  results[is.na(results)] <- 0
  setorder(results, calculated_m.z)
  
  return(results)
}
#maybe it makes sense to keep to full dataset and then only later after assignments keep CHO formulas?
formulas <- generate_formulas(C_range = 1:100, H_range = 2:150, O_range = 0:50, HC= c(0.2,3), OC = c(0,1.2), mass_range = c(90, 1000), insert_isotopes = T)
#formulas_calibration <- generate_formulas(C_range = 1:50, H_range = 2:100, O_range = 1:25, HC= c(1,2), OC = c(0.2,0.8), mass_range = c(100, 700), insert_isotopes = F)

#test <- generate_formulas(C_range = 1:50, H_range = 2:100, O_range = 0:50, HC= c(1,2), OC = c(0.1,1), mass_range = c(90, 800), insert_isotopes = T)
#test 

formula_assignment <- function(formulas, dataset, ppm_tolerance = 0.5, ion = 1.007825032 + -1 * 0.00054857990907) {
  
  # Adjust measured m/z values
  dataset$m.z_ion <- dataset$m.z
  dataset$m.z <- dataset$m.z + ion

  # Create data.table for measured m/z values
  
  
  # Calculate ppm tolerance limits
  dataset[, `:=`(
    lower_limit = m.z * (1 - ppm_tolerance / 1e6),
    upper_limit = m.z * (1 + ppm_tolerance / 1e6)
  )]
  
  # Create interval columns in results (assuming 'results' is available in the environment)
  formulas[, `:=`(
    start = calculated_m.z,
    end = calculated_m.z
  )]
  
  # Perform the overlap join using foverlaps
  setkey(formulas, start, end)
  setkey(dataset, lower_limit, upper_limit)
  final_matches <- foverlaps(formulas, dataset, 
                             by.x = c("start", "end"), 
                             by.y = c("lower_limit", "upper_limit"), 
                             type = "within", nomatch = 0)
  
  # Calculate actual ppm difference and absolute error difference
  # Calculate the actual ppm difference and absolute error difference
  final_matches[, absolute_error := calculated_m.z - m.z]
  final_matches[, ppm := (absolute_error / m.z) * 1e6]
  
  # Select columns of interest
  filtered_matches <- final_matches[, .(m.z_ion,m.z, calculated_m.z, absolute_error, ppm,
                                        C, H, O, N, S, P , C_13, O_18, S_34, I, MDL
  )]

  #filter out multi assignments
  # Step 1: Identify m.z values that are unique (appear exactly once)
  unique_mz <- filtered_matches[, .N, by = m.z][N == 1, m.z]
  
  # Step 2: Filter the original data.table to keep only rows with unique m.z values
  filtered_matches <- filtered_matches[m.z %in% unique_mz]
  
  filtered_matches
}

calculate_connections <- function(dt,   diffs = list(
  CH2 = c(C = 1, H = 2, O = 0, N = 0, S = 0, P = 0),
  #CH2O = c(C = 1, H = 2, O = 1, N = 0, S = 0, P = 0),
  #C2HO = c(C = 2, H = 1, O = 1, N = 0, S = 0, P = 0),
  #CO2 = c(C = 1, H = 0, O = 2, N = 0, S = 0, P = 0),
  #H2O = c(C = 0, H = 2, O = 1, N = 0, S = 0, P = 0),
  H = c(C = 0, H = 1, O = 0, N = 0, S = 0, P = 0)#,
  #O = c(C = 0, H = 0, O = 1, N = 0, S = 0, P = 0),
  #NH3 = c(C = 0, H = 3, O = 0, N = 1, S = 0, P = 0)
)) {
  
  
  
  setorder(dt, C, H, O)
  dt[, index := .I]
  dt[, group := fifelse(N > 0 & P > 0 , "NP", 
                        fifelse(N > 0 & S > 0 |N > 0 & S_34 > 0, "NS", 
                                fifelse(S > 0 & P >0|S_34 > 0 & P > 0 , "SP", 
                                        fifelse(N > 0 & S == 0 & P ==0 & S_34 ==0 , "N", 
                                                fifelse(S > 0 & N == 0 & P ==0 |S_34 > 0 & N == 0 & P ==0 , "S", 
                                                        fifelse(P > 0 & S == 0 & N ==0 & S_34 ==0 , "P", 
                                                                fifelse(O > 0 & P == 0 & S == 0 & N ==0 & S_34  == 0, "CHO", 
                                                                        fifelse(O == 0 & P == 0 & S == 0 & N ==0 & S_34 ==0 , "CH",
                                                                                "None"))))))))]
  
  
  results <- list()
  
  # Function to calculate connections for a specific group
  calculate_group_connections <- function(dt_group) {
    dt_group_unique <- dt_group#[!duplicated(MF)]
    group_results <- list()
    
    for (diff_name in names(diffs)) {
      diff_vector <- diffs[[diff_name]]
      
      # Calculate the new data.table after applying the difference
      diff_dt <- dt_group_unique[, .(
        C = C + diff_vector["C"],
        H = H + diff_vector["H"],
        O = O + diff_vector["O"],
        N = N + diff_vector["N"],
        S = S + diff_vector["S"],
        P = P + diff_vector["P"],
        from_index = index
      )]
      
      # Match the diff_dt with the original dt using a self-join
      setkeyv(diff_dt, c("C", "H", "O", "N", "S", "P"))
      setkeyv(dt_group_unique, c("C", "H", "O", "N", "S", "P"))
      
      matched <- dt_group_unique[diff_dt, nomatch = 0L]
      
      if (nrow(matched) > 0) {
        connections <- matched[, .(
          from_index = from_index,
          to_index = index,
          connection_type = diff_name
        )]
        group_results[[diff_name]] <- connections
      }
    }
    
    # Combine all results into one data.table
    all_group_connections <- rbindlist(group_results, use.names = TRUE, fill = TRUE)
    
    # Add entries with no connections
    no_connections <- dt_group_unique[!index %in% unique(c(all_group_connections$from_index, all_group_connections$to_index)), .(
      from_index = index,
      to_index = NA,
      connection_type = "No Connection"
    )]
    
    # Combine connections with no connections
    all_group_connections <- rbind(all_group_connections, no_connections, fill = TRUE)
    
    return(all_group_connections)
  }
  
  # Calculate connections for each group
  for (grp in unique(dt$group)) {
    dt_group <- dt[group == grp]
    group_connections <- calculate_group_connections(dt_group)
    results[[grp]] <- group_connections
  }
  
  # Combine all group results
  all_connections <- rbindlist(results, use.names = TRUE, fill = TRUE)
  
  # Create a graph and find components
  Ak <- graph_from_data_frame(all_connections[!is.na(to_index)], directed = FALSE) %>% components()
  
  # Create a data.table for membership and merge with dt by index
  membership_dt <- data.table(index = as.integer(names(Ak$membership)), homologues_membership = Ak$membership)
  
  
  # Count the number of members in each component
  # x <- membership_dt$homologues_membership
  # for (i in unique(x)) {
  #   count <- x[x == i] %>% length()
  #   membership_dt$homologues_membership[membership_dt$homologues_membership == i] <- count
  # }
  
  # Efficiently count the number of members in each homologues_membership
  membership_dt[, homologues_membership := .N, by = homologues_membership]
  
  out <- merge(dt, membership_dt, by = "index", all.x = TRUE)
  
  out[is.na(homologues_membership), homologues_membership := 1]
  
  return(dt_with_membership = out)
}


mz_recalibration <- function(formulas, dataset, ppm_tolerance = 0.5, ion =  1.007825032 + -1 * 0.00054857990907, S_MDL_limit = 1.5, hom_member_min = 1, showplot = T) {
  # Formula assignment and initial filtering
  filtered_matches <- formula_assignment(formulas, dataset, ppm_tolerance = ppm_tolerance, ion = ion)[P==0]#[N == 0 & S == 0 & P == 0 & C_13 == 0 & O_18 == 0 & S_34 == 0]
  out <- filtered_matches[I / MDL > S_MDL_limit]
  
  if (nrow(out) > 20) {
    out <- calculate_connections(out)
    out <- out[homologues_membership >= hom_member_min]
  }
  
  if (nrow(out) > 20) {
    # Apply quantile filtering
    Q1 <- quantile(out$ppm, 0.25)
    Q3 <- quantile(out$ppm, 0.75)
    IQR <- Q3 - Q1
    bounds <- c(Q1 - 1.5 * IQR, Q3 + 1.5 * IQR)
    out <- out[ppm > bounds[1] & ppm < bounds[2]]
    
    # Apply GAM model
    gam_model <- gam(ppm ~ s(m.z, bs = "cs", k = 4), data = out)
    # Function for applying drift correction
    correct_new_data <- function(m.z, gam_model) {
      new_data <- data.table(m.z)
      new_data[, predicted_drift_ppm := predict(gam_model, newdata = new_data)]
      new_data[, corrected_mz := m.z * (1 + predicted_drift_ppm / 1e6)]
      return(new_data)
    }
    
    pred <- correct_new_data(dataset$m.z, gam_model)
    dataset[, m.z := pred$corrected_mz]
    
    if(showplot){
    pred <- subset(pred, m.z <max(filtered_matches$m.z)) #for plotting
    
    # Recalculate and filter again after correction
    filtered_matches_2 <- formula_assignment(formulas, dataset, ppm_tolerance = ppm_tolerance, ion = ion)[P==0]#[N == 0 & S == 0 & P == 0 & C_13 == 0 & O_18 == 0 & S_34 == 0]
    out2 <- filtered_matches_2[I / MDL > S_MDL_limit]
    out2 <- calculate_connections(out2)[homologues_membership >= hom_member_min]
    
    # Apply quantile filter to the second output
    Q1_2 <- quantile(out2$ppm, 0.25)
    Q3_2 <- quantile(out2$ppm, 0.75)
    IQR_2 <- Q3_2 - Q1_2
    bounds_2 <- c(Q1_2 - 1.5 * IQR_2, Q3_2 + 1.5 * IQR_2)
    out2 <- out2[ppm > bounds_2[1] & ppm < bounds_2[2]]
    
    
    # Create the plot using plotly
    plot1 <- plot_ly() %>%
      add_trace(data = filtered_matches, x = ~m.z, y = ~ppm, type = 'scatter', mode = 'markers',
                marker = list(color = 'grey', size = 5), name = 'All assignments') %>%
      add_trace(data = out, x = ~m.z, y = ~ppm, type = 'scatter', mode = 'markers',
                marker = list(color = 'black', size = 5), name = 'Filtered assignments') %>%
      add_trace(data = out, x = pred$corrected_mz, y = pred$predicted_drift_ppm, type = 'scatter', mode = 'lines',
                line = list(color = 'red', width = 2), name = 'GAM Model') %>%
      layout(title = dataset$index[1], xaxis = list(title = "m/z"), yaxis = list(title = "ppm"))
    
    plot2 <- plot_ly() %>%
      add_trace(data = filtered_matches_2, x = ~m.z, y = ~ppm, type = 'scatter', mode = 'markers',
                marker = list(color = 'grey', size = 5), name = 'All assignments', showlegend = F) %>%
      add_trace(data = out2, x = ~m.z, y = ~ppm, type = 'scatter', mode = 'markers',
                marker = list(color = 'black', size = 5), name = 'Filtered assignments', showlegend = F) %>%
      layout(xaxis = list(title = "m/z"), yaxis = list(title = "ppm"))
    
    plot <- subplot(plot1, plot2, shareY = T, titleX = T, titleY = T)
    
    print(paste0(dataset$index[1], " SSE before: ",round(sum(out$ppm^2),2),
                 "| SSE after: ",round(sum(out2$ppm^2),2)))
    
    }else{
      plot <- plot_ly()
    }
    dataset <- dataset[, .(m.z, I, MDL, ResPow, index)]
  } else {
    print(paste0(dataset$index[1]," Not enough points"))
    plot <- plot_ly()
  }
  
  print(paste0("Finished ",dataset$index[1]))
  return(list(data = dataset, plot = plot))
}


MDL_water_recal <- lapply(MDL_water, function(dataset){
  res <- mz_recalibration(formulas, dataset, showplot = T,ppm_tolerance = 0.5)
})
names(MDL_water_recal) <- names(MDL_water)

MDL_water_recal_data<- lapply(MDL_water_recal, function(x){x$data})
MDL_water_recal_plots<- lapply(MDL_water_recal, function(x){x$plot})


##########all samples recalibrated. now merge?
#merge samples and increase mass accuracy
merge_mz_values <- function(dt, ppm_tolerance) {

  dt$I <- as.numeric(dt$I)
  
  # Rename the third column to "MDL" if necessary
  colnames(dt)[3] <- "MDL"
  colnames(dt)[1] <- "mz"
  
  # Sort the data by m/z (ascending), I (descending), and index (ascending)
  setorder(dt, mz)
  
  # Calculate the difference in ppm between consecutive m/z values
  dt[, diff_ppm := c(Inf, diff(mz)) / mz * 1e6]
  
  # Create a group identifier for values within the tolerance
  dt[, group := cumsum(diff_ppm > ppm_tolerance)]
  
  # Calculate the weighted mean m/z value, mean MDL, and standard error within each group
  merged <- dt[, {
    sqrt_I <- sqrt(I)                     # Compute sqrt(I) once
    sum_sqrt_I <- sum(sqrt_I)             # Sum of sqrt(I) only once
    sum_mz_sqrt_I <- sum(mz * sqrt_I)     # Weighted sum of mz
    
    weighted_mz <- sum_mz_sqrt_I / sum_sqrt_I  # Weighted mean of mz
    
    # Standard deviation can be optimized if needed, but we leave it as is for now
    .(
      merged_mz = weighted_mz,
      mean_MDL = MDL[1],
      SE = sd(mz) / weighted_mz * 1e6  # Standard error (in ppm)
    )
  }, by = .(group)]
  
  # Merge the weighted mean, mean MDL, and SE back to the original data
  dt <- merge(dt, merged, by = "group", all.x = TRUE)
  
  # Remove duplicates by selecting one row per group and sample
  dt2 <- dt[, .SD[1], by = .(group, index)]
  
  # Reshape data to wide format where each index gets a column
  dt_wide <- dcast(dt2, merged_mz + mean_MDL + SE ~ index, value.var = "I", fill = NA)
  
  # Calculate the number of samples (columns) with non-NA values for each merged m/z
  #dt_wide[, sample_count := rowSums(!is.na(.SD)), .SDcols = !c("merged_mz", "mean_MDL", "SE", "mean_mz")]
  
  # Remove rows where sample_count is 1
  #dt_wide <- dt_wide[sample_count > 1]
  
  # Optionally, remove the sample_count column if no longer needed
  #dt_wide[, sample_count := NULL]
  
  # Return the final data.table
  
  colnames(dt_wide)[c(1,2)]<- c("m.z", "MDL")
  
  return(dt_wide)
}

#combine all samples back together for increased mass precision?

results_water_dt <-  merge_mz_values(rbindlist(MDL_water_recal_data), ppm_tolerance = 0.6)

#function for formula assignments
calculate_isotopes <- function(filtered_matches, dataset){
  
  # Define the isotope_deviance function
  isotope_deviance <- function(intensity_parent, intensity_child, number_of_atoms, q, p, number_of_isotopes = 1) {
    expected_intensity_ratio <- dbinom(number_of_isotopes, number_of_atoms, q) / dbinom(number_of_atoms, number_of_atoms, p)
    ratio <- (intensity_child / intensity_parent)
    
    # Vectorized weighted mean with na.rm
    rdeviance_weighted <- ((weighted.mean(ratio, intensity_parent, na.rm = TRUE) / expected_intensity_ratio) - 1) * 1000
    return(rdeviance_weighted)
  }
  
  
  
  
  # Select intensity columns
  intensity_cols <- setdiff(names(dataset), c("m.z", "MDL", "SE", "ResPow", "m.z_ion","index"))
  
  # Check if there's more than one column, and apply rowMeans conditionally
  if (length(intensity_cols) > 1) {
    # If more than one intensity column, calculate rowMeans and keep m.z
    intensities <- dataset[, .(m.z, intensities = rowMeans(.SD, na.rm = TRUE)), .SDcols = intensity_cols]
  } else {
    # If only one intensity column, keep m.z and the single intensity column
    intensities <- dataset[, .(m.z, intensities = get(intensity_cols))]
  }
  
  # Merge intensities into filtered_matches
  filtered_matches <- merge(filtered_matches, intensities, by.x = "m.z", by.y = "m.z", all.x = TRUE)
  
  # Subset parent and child isotopes using data.table filtering
  parent_all <- filtered_matches[C_13 == 0 & O_18 == 0 & S_34 == 0]
  child_C_1_all <- filtered_matches[C_13 == 1 & O_18 == 0 & S_34 == 0]
  child_C_2_all <- filtered_matches[C_13 == 2 & O_18 == 0 & S_34 == 0]
  child_O_all <- filtered_matches[C_13 == 0 & O_18 == 1 & S_34 == 0]
  child_S_all <- filtered_matches[C_13 == 0 & O_18 == 0 & S_34 == 1]
  
  # Perform the merges and filtering
  
  
  result_C1 <- merge(parent_all, child_C_1_all, 
                     by = c("C", "H", "O", "N", "S", "P"),
                     all = FALSE, suffixes = c("_parent", "_child"))
  result_C2 <- merge(parent_all, child_C_2_all, 
                     by = c("C", "H", "O", "N", "S", "P"),
                     all = FALSE, suffixes = c("_parent", "_child"))
  
  result_O <- merge(parent_all, child_O_all, 
                    by = c("C", "H", "O", "N", "S", "P"),
                    all = FALSE, suffixes = c("_parent", "_child"))
  
  result_S <- merge(parent_all, child_S_all, 
                    by = c("C", "H", "O", "N", "S", "P"),
                    all = FALSE, suffixes = c("_parent", "_child"))
  
  result_C1 <- result_C1[intensities_parent > intensities_child]
  
  result_C2 <- result_C2[intensities_parent > intensities_child]
  
  result_O <- result_O[intensities_parent > intensities_child]
  
  result_S <- result_S[intensities_parent > intensities_child]
  
  
  #remove all non-matches?
  colnames_parent <- paste0("intensities" ,"_parent")
  colnames_child <- paste0("intensities", "_child")
  
  #can actually pre-compute the expected ratio for each C/O/S already
  #does not seem to be necessary though, already extremely fast
  
  
  # Calculate isotope deviance
  result_C1[, isotope_deviance_C_13_1 := isotope_deviance(
    intensity_parent = as.numeric(get(colnames_parent)),    # Column from parent
    intensity_child = as.numeric(get(colnames_child)),      # Column from child
    number_of_atoms = C,                   
    q = 0.0108,                             
    p = 1-0.0108,
    number_of_isotopes =1
  ), by = seq_len(nrow(result_C1))]
  
  result_C2[, isotope_deviance_C_13_2 := isotope_deviance(
    intensity_parent = as.numeric(get(colnames_parent)),    # Column from parent
    intensity_child = as.numeric(get(colnames_child)),      # Column from child
    number_of_atoms = C,                   
    q = 0.0108,                             
    p = 1-0.0108,
    number_of_isotopes =2
  ), by = seq_len(nrow(result_C2))]
  
  result_O[, isotope_deviance_O_18 := isotope_deviance(
    intensity_parent = as.numeric(get(colnames_parent)),    # Column from parent
    intensity_child = as.numeric(get(colnames_child)),      # Column from child
    number_of_atoms = O,                   
    q = 0.00205,                             
    p = 1-0.00205,
    number_of_isotopes =1
  ), by = seq_len(nrow(result_O))]
  
  result_S[, isotope_deviance_S_34 := isotope_deviance(
    intensity_parent = as.numeric(get(colnames_parent)),    # Column from parent
    intensity_child = as.numeric(get(colnames_child)),      # Column from child
    number_of_atoms = S,                   
    q = 0.0437,                             
    p = 1-0.0437,
    number_of_isotopes =1
  ), by = seq_len(nrow(result_S))]
  
  result <-rbind(parent_all, result_C1, result_C2, result_O, result_S, fill = T)
  #child and parent columns need to be un-merged.
  
  
  # Identify columns with '_parent' and '_child'
  parent_cols <- grep("_parent$", colnames(result_C1), value = TRUE)
  child_cols <- grep("_child$", colnames(result_C1), value = TRUE)
  
  # Combine with C H O N S P columns (which are common)
  common_cols <- c("C", "H", "O", "N", "S", "P")
  
  # Create separate data.tables for parent and child
  parent_C1 <- result_C1[, c(common_cols, parent_cols,"isotope_deviance_C_13_1"), with = FALSE]
  result_C1 <- result_C1[, c(common_cols, child_cols,"isotope_deviance_C_13_1"), with = FALSE]
  
  parent_C2 <- result_C2[, c(common_cols, parent_cols,"isotope_deviance_C_13_2"), with = FALSE]
  result_C2 <- result_C2[, c(common_cols, child_cols,"isotope_deviance_C_13_2"), with = FALSE]
  
  parent_O <- result_O[, c(common_cols, parent_cols,"isotope_deviance_O_18"), with = FALSE]
  result_O <- result_O[, c(common_cols, child_cols,"isotope_deviance_O_18"), with = FALSE]
  
  parent_S <- result_S[, c(common_cols, parent_cols,"isotope_deviance_S_34"), with = FALSE]
  result_S <- result_S[, c(common_cols, child_cols,"isotope_deviance_S_34"), with = FALSE]
  
  # Optionally, remove the '_parent' and '_child' suffixes for better readability
  setnames(parent_C1, gsub("_parent$", "", colnames(parent_C1)))
  setnames(result_C1, gsub("_child$", "", colnames(result_C1)))
  setnames(parent_C2, gsub("_parent$", "", colnames(parent_C2)))
  setnames(result_C2, gsub("_child$", "", colnames(result_C2)))
  setnames(parent_O, gsub("_parent$", "", colnames(parent_O)))
  setnames(result_O, gsub("_child$", "", colnames(result_O)))
  setnames(parent_S, gsub("_parent$", "", colnames(parent_S)))
  setnames(result_S, gsub("_child$", "", colnames(result_S)))
  
  parent_C1[, intensities := NULL]
  parent_C2[, intensities := NULL]
  parent_O[, intensities := NULL]
  parent_S[, intensities := NULL]
  result_C1[, intensities := NULL]
  result_C2[, intensities := NULL]
  result_O[, intensities := NULL]
  result_S[, intensities := NULL]
  
  parent_all[, intensities := NULL]
  
  
  # parents <- merge(parent_C1, parent_C2, 
  #                  by = colnames(parent_C1)[-ncol(parent_C1)], all = T) #in this case drop non-matches for C2, exclude isotope peaks for C2 where no C1 is present
  # idx <-which(is.na(parents$isotope_deviance_C_13_1))
  # parents <- parents [-idx,]
  # parents <- merge(parents, parent_O, 
  #                  by = colnames(parent_C1)[-ncol(parent_C1)], all = T)
  # parents <- merge(parents, parent_S, 
  #                  by = colnames(parent_C1)[-ncol(parent_C1)], all = T)
  isotopes <-merge(result_C1, result_C2, 
                   by = colnames(result_C1)[-ncol(result_C1)], all = T) 
  #idx <-which(is.na(isotopes$isotope_deviance_C_13_1))
  #isotopes <- isotopes [-idx,]
  isotopes <- merge(isotopes, result_O, 
                    by = colnames(parent_C1)[-ncol(parent_C1)], all = T)
  isotopes <- merge(isotopes, result_S, 
                    by = colnames(parent_C1)[-ncol(parent_C1)], all = T)
  #final_data <- merge(filtered_matches, parents,by = colnames(parent_C1)[-ncol(parent_C1)], all = T)
  #...and so on for other isotopes. Insert for parent as well?
  
  # final <- merge(isotopes, parents, 
  #                            by = colnames(isotopes), all = T)
  out <- rbind(parent_all, isotopes, fill =T)
}

formula_assignment_multi_intensity <- function(formulas, dataset, ppm_tolerance = 0.5, ion = 1.007825032 + -1 * 0.00054857990907, threshold = 1000, return_likeliest = T,
                                               homologues_network = list(
                                                 CH2 = c(C = 1, H = 2, O = 0, N = 0, S = 0, P = 0),
                                                 #CH2O = c(C = 1, H = 2, O = 1, N = 0, S = 0, P = 0),
                                                 #C2HO = c(C = 2, H = 1, O = 1, N = 0, S = 0, P = 0),
                                                 #CO2 = c(C = 1, H = 0, O = 2, N = 0, S = 0, P = 0),
                                                 #H2O = c(C = 0, H = 2, O = 1, N = 0, S = 0, P = 0),
                                                 H = c(C = 0, H = 1, O = 0, N = 0, S = 0, P = 0)#,
                                                 #O = c(C = 0, H = 0, O = 1, N = 0, S = 0, P = 0),
                                                 #NH3 = c(C = 0, H = 3, O = 0, N = 1, S = 0, P = 0)
                                               )) {
  
  # Adjust measured m/z values
  dataset$m.z_ion <- dataset$m.z
  dataset$m.z <- dataset$m.z + ion
  
  # Create data.table for measured m/z values
  
  
  # Calculate ppm tolerance limits
  dataset[, `:=`(
    lower_limit = m.z * (1 - ppm_tolerance / 1e6),
    upper_limit = m.z * (1 + ppm_tolerance / 1e6)
  )]
  
  # Create interval columns in results (assuming 'results' is available in the environment)
  formulas[, `:=`(
    start = calculated_m.z,
    end = calculated_m.z
  )]
  
  # Perform the overlap join using foverlaps
  setkey(formulas, start, end)
  setkey(dataset, lower_limit, upper_limit)
  all_matches <- foverlaps(formulas, dataset, 
                             by.x = c("start", "end"), 
                             by.y = c("lower_limit", "upper_limit"), 
                             type = "within", nomatch = 0)
  
  dataset[, `:=`(lower_limit = NULL, upper_limit= NULL)]
  
  # Calculate actual ppm difference and absolute error difference
  # Calculate the actual ppm difference and absolute error difference
  all_matches[, absolute_error := calculated_m.z - m.z]
  all_matches[, ppm := (absolute_error / m.z) * 1e6]
  
  # Select columns of interest
  filtered_matches <- all_matches[, .(m.z_ion,m.z, calculated_m.z, absolute_error, ppm,
                                        C, H, O, N, S, P , C_13, O_18, S_34
  )]
  all_matches <- NULL
  #isotope peak verification
  #this function automatically removes isotopes where no parent is present in the dataset and isotopes where intensity of parent is below intensity of child.
  filtered_matches <- calculate_isotopes(filtered_matches, dataset = dataset)
  #setorder(filtered_matches, m.z)
  
  #only drop isotope formulas with high deviance if not multi assignment?
  #homologues series calculation
  filtered_matches <- calculate_connections ( filtered_matches,  diffs = homologues_network)
  
  filtered_matches[, `:=`(index = NULL)]

  # Create a boolean column that indicates if any isotope deviance is below threshold
  filtered_matches[, isotope_below_threshold := (isotope_deviance_C_13_1 < threshold |
                                                   isotope_deviance_C_13_2 < threshold |
                                                   isotope_deviance_O_18 < threshold |
                                                   isotope_deviance_S_34 < threshold)]
  
  if(return_likeliest){

  # Create a new column for the absolute value of ppm
  filtered_matches[, abs_ppm := abs(ppm)]
  
  # Now sort the data.table by m.z, isotope_below_threshold, homologues_membership, and abs_ppm
  setorder(filtered_matches, m.z, -isotope_below_threshold, -homologues_membership, abs_ppm)
  
  # Keep only the first row for each unique m.z after applying the sorting
  filtered_matches <- filtered_matches[, .SD[1], by = m.z]
  
  # Remove the helper columns 'isotope_below_threshold' and 'abs_ppm' if they're not needed
  filtered_matches[, `:=`(abs_ppm = NULL)]
  
  #if formulas were removed this way, recalculate the homologues network
  filtered_matches[, `:=`(homologues_membership = NULL)]
  
  filtered_matches <- calculate_connections ( filtered_matches,  diffs = homologues_network)
  }
  #to do as next step: insert molecular formula, calculate indices and assign compound groups.
  #merge with intensity.
  #formula
  # Construct molecular formula
  construct_formula <- function(elements) {
    formula_parts <- paste0(names(elements), ifelse(elements == 1, "", elements))
    formula <- paste(formula_parts[elements > 0], collapse = "")
    return(formula)
  }
  
  filtered_matches[, MF := construct_formula(
    list(
      C = C ,
      H = H ,
      O = O ,
      N = N ,
      S = S ,
      P = P
    )
  ), by = 1:nrow(filtered_matches)]
  
  
  #
  #calculate further indices
  # Calculate AI.mod and add it as a new column
  filtered_matches[, AI_mod := (1 + C - 0.5 * O - S - 0.5 * (H + N + P )) / 
                        ((C) - 0.5 * (O) - (S) - (N) - P)]
  # Set AI_mod values below 0 and infinite values to 0
  filtered_matches[, AI_mod := ifelse(AI_mod < 0 | is.infinite(AI_mod), 0, AI_mod)]
  
  # Remove rows where AI_mod is above 1, unlikely formulas
  filtered_matches <- filtered_matches[AI_mod <= 1]
  
  filtered_matches[, NOSC := (-((4* (C) + (H) -3*(N) -2*(O) +5*P  -2* (S))/(C))+4)]
  
  filtered_matches[, O.C := O/C]
  filtered_matches[, H.C := H/C]
  filtered_matches[, N.C := N/C]
  filtered_matches[, S.C := S/C]
  filtered_matches[, P.C := P/C]
  # Create the 'class' column 
  filtered_matches[, class := fifelse(O.C > 0 & O.C <= 0.3 & H.C > 1.5 & H.C <= 2.5, "Lipid",
                                         fifelse(O.C <= 0.125 & H.C >= 1 & H.C <= 1.5 & AI_mod > 0.5, "Unsaturated hydrocarbon",
                                                 fifelse(O.C > 0.3 & O.C <= 0.55 & H.C > 1.5 & H.C <= 2.3, "Protein",
                                                         fifelse(O.C > 0.55 & O.C <= 0.7 & H.C > 1.5 & H.C <= 2.2, "Aminosugar",
                                                                 fifelse(O.C > 0.7 & O.C <= 1.05 & H.C > 1.5 & H.C <= 2.2, "Carbohydrate",
                                                                         fifelse(AI_mod >= 0.67, "Condensed hydrocarbon",
                                                                                 fifelse(AI_mod < 0.67 & AI_mod > 0.5 & O.C > 0.125, "Aromatic",
                                                                                         fifelse(O.C > 0.125 & AI_mod <= 0.5 & H.C <= 1.5, "Highly unsaturated", 
                                                                                                 "Unassigned"))))))))]

  
    # Select columns of interest
  filter <- c("m.z_ion", "m.z","calculated_m.z", "absolute_error", "ppm","MF",  "group","class",
              "C", "H", "O", "N", "S", "P" , "C_13", "O_18", "S_34", "H.C", "O.C", "N.C", "S.C", "P.C","AI_mod", "NOSC","homologues_membership",
              "isotope_deviance_C_13_1","isotope_deviance_C_13_2","isotope_deviance_O_18","isotope_deviance_S_34","isotope_below_threshold")
  filtered_matches <-  filtered_matches[, ..filter]
  # Merge intensities into filtered_matches
  dataset[, `:=`(m.z = NULL)]
  filtered_matches <- merge(filtered_matches, dataset, by.x = "m.z_ion", by.y = "m.z_ion", all.x = TRUE)
  
  setorder(filtered_matches, m.z)
  
  return(filtered_matches)
}
diffs <- list(
  CH2 = c(C = 1, H = 2, O = 0, N = 0, S = 0, P = 0),
  CH2O = c(C = 1, H = 2, O = 1, N = 0, S = 0, P = 0),
  C2HO = c(C = 2, H = 1, O = 1, N = 0, S = 0, P = 0),
  CO2 = c(C = 1, H = 0, O = 2, N = 0, S = 0, P = 0),
  H2O = c(C = 0, H = 2, O = 1, N = 0, S = 0, P = 0),
  H = c(C = 0, H = 1, O = 0, N = 0, S = 0, P = 0),
  O = c(C = 0, H = 0, O = 1, N = 0, S = 0, P = 0),
  NH3 = c(C = 0, H = 3, O = 0, N = 1, S = 0, P = 0)
)
crosstab_water_unambigous <- formula_assignment_multi_intensity(formulas, results_water_dt,ppm_tolerance = 0.3, return_likeliest = T, homologues_network = diffs)

#test for Na and Cl additions?
formula_assignment_check_adducts <- function(formulas, dataset, ppm_tolerance = 0.3,  homologues_network = list(
  CH2 = c(C = 1, H = 2, O = 0, N = 0, S = 0, P = 0),
  CH2O = c(C = 1, H = 2, O = 1, N = 0, S = 0, P = 0),
  C2HO = c(C = 2, H = 1, O = 1, N = 0, S = 0, P = 0),
  CO2 = c(C = 1, H = 0, O = 2, N = 0, S = 0, P = 0),
  H2O = c(C = 0, H = 2, O = 1, N = 0, S = 0, P = 0),
  H = c(C = 0, H = 1, O = 0, N = 0, S = 0, P = 0),
  O = c(C = 0, H = 0, O = 1, N = 0, S = 0, P = 0),
  NH3 = c(C = 0, H = 3, O = 0, N = 1, S = 0, P = 0)
)){
test_H <- formula_assignment_multi_intensity(formulas, dataset,ppm_tolerance = ppm_tolerance, homologues_network = homologues_network)
test_Na <- formula_assignment_multi_intensity(formulas, dataset,ppm_tolerance = ppm_tolerance, homologues_network = homologues_network, ion = -20.974666)
test_Cl <- formula_assignment_multi_intensity(formulas, dataset,ppm_tolerance = ppm_tolerance, homologues_network = homologues_network, ion = -34.969402)

test_H[, type := "H"]
test_Na[, type := "Na"]
test_Cl[, type := "Cl"]

merge <- rbind(test_H, test_Na,test_Cl)
# Create a new column for the absolute value of ppm
merge[, abs_ppm := abs(ppm)]

# Now sort the data.table by m.z, isotope_below_threshold, homologues_membership, and abs_ppm
setorder(merge, m.z_ion, -isotope_below_threshold, -homologues_membership, abs_ppm)

# Keep only the first row for each unique m.z after applying the sorting
merge <- merge[, .SD[1], by = m.z_ion]

# Remove the helper columns 'isotope_below_threshold' and 'abs_ppm' if they're not needed
merge[, `:=`(abs_ppm = NULL)]

#this checks for sodium adducts. seems to remove likely wrongly assigned P formulas.
test1 <- subset(merge, type == "H")

#recalculate homologues series
test1[, `:=`(homologues_membership = NULL)]

test1 <- calculate_connections ( test1,  diffs = homologues_network)
}

#unite data first in duplicates before assignment?
data_single <- lapply(seq(1,24,2), function(i){
  dat <- rbind(MDL_water_recal_data[[i]], MDL_water_recal_data[[i+1]])
  data <- merge_mz_values(dat, ppm_tolerance = 0.5)
  b <- rowMeans(data[,c(4,5)])
  idx <- which(!is.na(b))
  data <- data.table(data[idx,c(1:3),],b[idx])
})

crosstab_single <- lapply(data_single, function(sample){
  out <- formula_assignment_check_adducts(formulas,sample, ppm_tolerance = 0.3)[group != "CH" &homologues_membership >1 & C_13 ==0 & O_18 ==0 & S_34==0]
  #apply the first filter already
  
})





