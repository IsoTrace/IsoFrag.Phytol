####### Modelling isotopic structure of Phytols #############
#------------------
# created by Merve Öztoprak
# m.oeztoprak@gmail.com
#------------------

# --- Constants ---
R_13C.VPDB <- 0.01118  # VPDB 13C/12C ratio
R_2H.SMOW <- 0.00015575
R_17O.SMOW <- 0.0003799
R_18O.SMOW <- 0.0020672

# --- Define δ13C values for each carbon position in phytol (modifiable) ---
configure_d13C <- function(mode = "homogeneous", 
                           homogeneous_value = -20.3,
                           head_value = NULL, 
                           tail_value = NULL,
                           custom_values = NULL) {
  carbon_positions <- paste0("C", 1:20)
  # Define fragment groups for easy configuration
  fragment_groups <- list(
    head = c("C1","C2","C3","C4","C5","C6","C7","C17","C18"),
    tail = c("C8","C9","C10","C11","C12","C13","C14","C15","C16","C19","C20"),
    whole = carbon_positions
  )
  # Initialize all values to NA to catch configuration errors
  d13C <- setNames(rep(NA, 20), carbon_positions)
  if (mode == "homogeneous") {
    # All carbons same value
    d13C[] <- homogeneous_value
  } else if (mode == "head_tail") {
    # Head and tail different values
    if (is.null(head_value) || is.null(tail_value)) {
      stop("In head_tail mode, both head_value and tail_value must be provided")
    }
    d13C[fragment_groups$head] <- head_value
    d13C[fragment_groups$tail] <- tail_value
  } else if (mode == "custom") {
    # Fully custom values for each position
    if (is.null(custom_values) || length(custom_values) != 20) {
      stop("In custom mode, you must provide exactly 20 values (one for each carbon)")
    }
    d13C[] <- custom_values
  } else {
    stop("Invalid mode. Choose 'homogeneous', 'head_tail', or 'custom'")
  }
  # Check for any NA values that weren't properly configured
  if (any(is.na(d13C))) {
    stop("Some carbon positions were not properly configured. Missing values for: ",
         paste(names(d13C)[is.na(d13C)], collapse = ", "))
  }
  return(d13C)
}

# Hydrogen and oxygen isotope values (modifiable)
d2H_Phyt.H <- -200   # Phytol hydrogen δ2H
d2H_H.OH <- 0           # Hydroxyl hydrogen δ2H
d2H_H.APCI <- 0         # APCI hydrogen δ2H
d17O_Phyt.O <- 0        # Phytol oxygen δ17O

# --- Conversion functions ----
# Convert δ13C to fractional abundance (f)
delta_to_fraction <- function(delta, R_std) {
  R_sample <- R_std * (1 + delta / 1000)
  return(R_sample / (1 + R_sample))
}
# Convert δ13C to R (ratio)
delta_to_R <- function(delta, R_std) {
  R_std * (1 + delta / 1000)
}
# --- Calculate fractional abundances for all positions ---
calculate_fractional_abundances <- function(d13C, d2H_Phyt.H, d2H_H.OH, d2H_H.APCI, d17O_Phyt.O) {
  # Carbon positions (1-20)
  f_13C <- sapply(d13C, function(x) delta_to_fraction(x, R_13C.VPDB))
  # Phytol hydrogens (37 positions)
  f_2H_Phyt <- delta_to_fraction(d2H_Phyt.H, R_2H.SMOW)
  # Hydroxyl hydrogen (1 position)
  f_2H_OH <- delta_to_fraction(d2H_H.OH, R_2H.SMOW)
  # APCI hydrogen (1 position)
  f_2H_APCI <- delta_to_fraction(d2H_H.APCI, R_2H.SMOW)
  # Phytol oxygen (1 position)
  f_17O_Phyt <- delta_to_fraction(d17O_Phyt.O, R_17O.SMOW)
  return(list(
    f_13C = f_13C,
    f_2H_Phyt = f_2H_Phyt,
    f_2H_OH = f_2H_OH,
    f_2H_APCI = f_2H_APCI,
    f_17O_Phyt = f_17O_Phyt
  ))
}
# --- Calculate molecular averages ---
calculate_molecular_averages <- function(fractions) {
  # Calculate unsubstituted molecular abundance 
  # This is the product of all (1-f) terms for 12C, plus other isotopes
  unsubstituted <- prod(1 - fractions$f_13C) *       # All carbons as 12C
    (1 - fractions$f_17O_Phyt) *                     # 16O
    (1 - fractions$f_2H_OH)^1 *                      # 1H in OH
    (1 - fractions$f_2H_APCI)^1 *                    # 1H in APCI
    (1 - fractions$f_2H_Phyt)^37                     # 37 1H in phytol chain
  return(unsubstituted)
}
# --- Calculate M+1 relative abundances---
calculate_M1_abundances <- function(fractions, unsubstituted) {
  # Initialize vector for M+1 abundances
  M1_abund <- numeric(32)  # For positions C1-C20, Phyt.H, H.OH, H.APCI, Phyt.O
  names(M1_abund) <- c(paste0("C", 1:20), "Phyt.H", "H.OH", "H.APCI", "Phyt.O")
  # Calculate for each carbon position
  for (i in 1:20) {
    carbon <- paste0("C", i)
    M1_abund[carbon] <- (unsubstituted / (1 - fractions$f_13C[i])) * fractions$f_13C[i]
  }
  # Calculate for phytol hydrogens
  M1_abund["Phyt.H"] <- 37 * (unsubstituted / (1 - fractions$f_2H_Phyt)) * fractions$f_2H_Phyt
  # Calculate for hydroxyl hydrogen
  M1_abund["H.OH"] <- (unsubstituted / (1 - fractions$f_2H_OH)) * fractions$f_2H_OH
  # Calculate for APCI hydrogen (Excel AG6)
  M1_abund["H.APCI"] <- (unsubstituted / (1 - fractions$f_2H_APCI)) * fractions$f_2H_APCI
  # Calculate for phytol oxygen (Excel AH6)
  M1_abund["Phyt.O"] <- (unsubstituted / (1 - fractions$f_17O_Phyt)) * fractions$f_17O_Phyt
  return(M1_abund)
}
# --- Calculate relative abundances ---
calculate_relative_abundances <- function(M1_abund) {
  total_M1 <- sum(M1_abund)
  rel_abund <- M1_abund / total_M1
  return(rel_abund)
}
# --- Define which carbon positions belong to each fragment ----
fragment_definitions <- list(
  C4H7 = c("C16", "C15", "C14", "C20"),
  C4H9 = c("C16", "C15", "C14", "C20"),
  C4H7O = c("C1", "C2", "C3", "C17"),
  C5H7 = c("C1", "C2", "C3", "C17", "C4"),
  C5H9 = c("C16", "C15", "C14", "C20","C13"),
  C5H11 = c("C16", "C15", "C14", "C20","C13"),
  C5H9O = c("C1", "C2", "C3", "C17", "C4"),
  C6H9 = c("C1", "C2", "C3", "C17", "C4", "C5"),
  C6H11 = c("C16", "C15", "C14", "C20", "C13", "C12"),
  C6H13 = c("C16", "C15", "C14", "C20", "C13", "C12"),
  C6H11O = c("C1", "C2", "C3", "C17", "C4", "C5"),
  C7H11 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6"),
  C7H13 = c("C16", "C15", "C14", "C13", "C12", "C11", "C20"),
  C7H15 = c("C16", "C15", "C14", "C13", "C12", "C11", "C20"),
  C7H13O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6"),
  C8H13 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7"),
  C8H15 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20"),
  C8H17 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20"),
  C8H15O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7"),
  C9H15 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18"),
  C9H17 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10"),
  C9H19 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10"),
  C9H17O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18"),
  C10H17 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8"),
  C10H19 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9"),
  C10H21 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9"),
  C10H19O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8"),
  C11H19 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9"),
  C11H21 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8"),
  C11H23 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8"),
  C11H21O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9"),
  C12H21 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9", "C10"),
  C12H23 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8", "C7"),
  C12H25 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8", "C7"),
  C12H23O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9", "C10"),
  C13H23 = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9", "C10", "C11"),
  C13H25 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8", "C7", "C18"),
  C13H27 = c("C16", "C15", "C14", "C13", "C12", "C11", "C19", "C20", "C10", "C9", "C8", "C7", "C18"),
  C13H25O = c("C1", "C2", "C3", "C17", "C4", "C5", "C6", "C7", "C18", "C8", "C9", "C10", "C11")
) 
# --- Calculate r13C theory for each fragment ----
calculate_r13C_theory <- function(rel_abund, fragments) {
  r13C_theory <- sapply(fragments, function(carbons) {
    sum(rel_abund[carbons])
  })
  return(r13C_theory)
}
# --- run the model ----
run_isotope_model <- function(d13C, d2H_Phyt.H, d2H_H.OH, d2H_H.APCI, d17O_Phyt.O) {
  # Calculate all fractional abundances
  fractions <- calculate_fractional_abundances(d13C, d2H_Phyt.H, d2H_H.OH, d2H_H.APCI, d17O_Phyt.O)
  # Calculate molecular averages
  unsubstituted <- calculate_molecular_averages(fractions)
  # Calculate M+1 abundances
  M1_abund <- calculate_M1_abundances(fractions, unsubstituted)
  # Calculate relative abundances
  rel_abund <- calculate_relative_abundances(M1_abund)
  # Calculate r13C theory for each fragment
  r13C_theory <- calculate_r13C_theory(rel_abund, fragment_definitions)
  return(list(
    fractional_abundances = fractions,
    unsubstituted = unsubstituted,
    M1_abundances = M1_abund,
    relative_abundances = rel_abund,
    r13C_theory = r13C_theory
  ))
}

# ---Main wrapper function for each filename -----
calculate_theoretical_fragments <- function(sample_data, d13C_config) {
  # Validate input
  required_cols <- c("filename", "deltaC", "deltaH", "deltaOH", "deltaAPCI", "deltaO", "deltaO18")
  missing_cols <- setdiff(required_cols, names(sample_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Process each sample
  results <- lapply(1:nrow(sample_data), function(i) {
    sample <- sample_data[i, ]
    
    # Configure d13C
    d13C_values <- configure_d13C(
      mode = d13C_config$mode,
      homogeneous_value = if (d13C_config$mode == "homogeneous") sample$deltaC else NULL,
      head_value = if (d13C_config$mode == "head_tail") d13C_config$head_value else NULL,
      tail_value = if (d13C_config$mode == "head_tail") d13C_config$tail_value else NULL,
      custom_values = if (d13C_config$mode == "custom") d13C_config$custom_values else NULL
    )
    
    # Run model with sample-specific parameters
    model_results <- run_isotope_model(
      d13C = d13C_values,
      d2H_Phyt.H = sample$deltaH,
      d2H_H.OH = sample$deltaOH,
      d2H_H.APCI = sample$deltaAPCI,
      d17O_Phyt.O = sample$deltaO
    )
    
    # Convert to dataframe
    data.frame(
      filename = sample$filename,
      ID = sample$ID,
      date = sample$date,
      replicate = sample$replicate,
      Fragment = names(model_results$r13C_theory),
      r13C_Theory = unname(model_results$r13C_theory),
      stringsAsFactors = FALSE
    )
  })
  
  # Combine all results
  bind_rows(results)
}

#---- Optimization function to minimize difference between observed and theoretical ----
optimize_d13C <- function(observed, fragment_defs, lower_bound, upper_bound, d2H, initial_guess = NULL) {
    # Carbon positions in phytol
  carbon_positions <- paste0("C", 1:20)
    # Set initial guess if not provided
  if(is.null(initial_guess)) {
    initial_guess <- setNames(rep(-20, 20), carbon_positions)
  }
    # Error function to minimize
  error_fn <- function(params) {
    # Extract current d13C values
    current_d13C <- setNames(params, carbon_positions)
        # Run model with current parameters
    results <- run_isotope_model(
      d13C = current_d13C,
      d2H_Phyt.H = d2H,
      d2H_H.OH = 0,
      d2H_H.APCI = 0,
      d17O_Phyt.O = 0
    )
        # Calculate squared errors for observed fragments
    errors <- sapply(observed$Fragment, function(frag) {
      if(frag %in% names(results$r13C_theory)) {
        (results$r13C_theory[frag] - observed$r13C_Observed[observed$Fragment == frag])^2
      } else {
        NA # Fragment not modeled
      }
    })
        # Return sum of squared errors (excluding NAs)
    sum(errors, na.rm = TRUE)
  }
    # Run optimization
  opt_result <- optim(
    par = initial_guess,
    fn = error_fn,
    method = "L-BFGS-B",
    lower = lower_bound, # Reasonable lower bound
    upper = upper_bound     # Reasonable upper bound
  )
    # Get optimized d13C values
  optimized_d13C <- setNames(opt_result$par, carbon_positions)
    # Run final model with optimized values
  final_results <- run_isotope_model(
    d13C = optimized_d13C,
    d2H_Phyt.H = d2H,
    d2H_H.OH = 0,
    d2H_H.APCI = 0,
    d17O_Phyt.O = 0
  )
    # Prepare output
  list(
    optimized_d13C = optimized_d13C,
    theoretical_values = data.frame(
      Fragment = names(final_results$r13C_theory),
      r13C_Theory = unname(final_results$r13C_theory)
    ),
    observed_values = observed,
    convergence = opt_result$convergence,
    message = opt_result$message
  )
}
# Modified optimization function with positional constraints for 85m/z to be 10 permil lower then 71 m/z
optimize_d13C_constrained <- function(observed, fragment_defs_raw, d2H, initial_guess = NULL, 
                                      use_analytical_gradient = FALSE) {
  
  carbon_positions <- paste0("C", 1:20)
  
  fragment_defs <- lapply(fragment_defs_raw, function(carb_vec) {
    list(carbons = carb_vec)
  })
  
  if(is.null(initial_guess)) {
    initial_guess <- setNames(rep(-20, 20), carbon_positions)
  }
  
  # Define constraint relationships
  constraints <- list(
    list(parent = "C4", child = "C5", offset = -10),
    list(parent = "C8", child = "C9", offset = -10),
    list(parent = "C12", child = "C13", offset = -10)
  )
  
  # Determine which parameters to optimize (exclude constrained children)
  optim_params <- setdiff(carbon_positions, sapply(constraints, function(x) x$child))
  
  # Adjust initial guess to satisfy constraints
  for(constraint in constraints) {
    initial_guess[constraint$child] <- initial_guess[constraint$parent] + constraint$offset
  }
  
  # Error function with constraints
  error_fn <- function(params) {
    current_d13C <- setNames(params, carbon_positions)
    
    # Apply constraints
    for(constraint in constraints) {
      current_d13C[constraint$child] <- current_d13C[constraint$parent] + constraint$offset
    }
    
    results <- run_isotope_model(
      d13C = current_d13C,
      d2H_Phyt.H = d2H,
      d2H_H.OH = 0,
      d2H_H.APCI = 0,
      d17O_Phyt.O = 0
    )
    
    errors <- sapply(observed$Fragment, function(frag) {
      if(frag %in% names(results$r13C_theory)) {
        (results$r13C_theory[frag] - observed$r13C_Observed[observed$Fragment == frag])^2
      } else {
        NA
      }
    })
    
    sum(errors, na.rm = TRUE)
  }
  
  # Analytical gradient function (if available)
  if (use_analytical_gradient) {
    gradient_fn <- function(params) {
      current_d13C <- setNames(params, carbon_positions)
      
      # Apply constraints
      for(constraint in constraints) {
        current_d13C[constraint$child] <- current_d13C[constraint$parent] + constraint$offset
      }
      
      # Run model to get theoretical values
      results <- run_isotope_model(
        d13C = current_d13C,
        d2H_Phyt.H = d2H,
        d2H_H.OH = 0,
        d2H_H.APCI = 0,
        d17O_Phyt.O = 0
      )
      
      # Initialize gradient vector
      grad <- numeric(length(params))
      names(grad) <- names(params)
      
      # Compute gradient for each fragment
      for(frag in observed$Fragment) {
        if(frag %in% names(results$r13C_theory)) {
          # Get the carbon atoms contributing to this fragment
          contributing_carbons <- fragment_defs[[frag]]$carbons
          
          # Compute derivative for each contributing carbon
          for(carbon in contributing_carbons) {
            # Skip constrained parameters (they're not independent)
            if(carbon %in% optim_params) {
              error <- results$r13C_theory[frag] - observed$r13C_Observed[observed$Fragment == frag]
              grad[carbon] <- grad[carbon] + 2 * error
            }
          }
        }
      }
      
      # Return only the gradient for optimizable parameters
      grad[optim_params]
    }
  } else {
    gradient_fn <- NULL
  }
  
  # Wrapper function for optim that reconstructs full parameter set
  optim_wrapper <- function(opt_params) {
    params <- initial_guess
    params[optim_params] <- opt_params
    
    # Reapply constraints in case optim suggests values that would violate them
    for(constraint in constraints) {
      params[constraint$child] <- params[constraint$parent] + constraint$offset
    }
    
    error_fn(params)
  }
  
  # Gradient wrapper (if using analytical gradient)
  if (use_analytical_gradient) {
    gradient_wrapper <- function(opt_params) {
      params <- initial_guess
      params[optim_params] <- opt_params
      
      # Reapply constraints
      for(constraint in constraints) {
        params[constraint$child] <- params[constraint$parent] + constraint$offset
      }
      
      gradient_fn(params)
    }
  } else {
    gradient_wrapper <- NULL
  }
  
  # Run optimization
  opt_result <- optim(
    par = initial_guess[optim_params],
    fn = optim_wrapper,
    gr = gradient_wrapper,
    method = "L-BFGS-B",
    lower = rep(-50, length(optim_params)),
    upper = rep(20, length(optim_params))
  )
  
  # Reconstruct full solution
  optimized_d13C <- initial_guess
  optimized_d13C[optim_params] <- opt_result$par
  for(constraint in constraints) {
    optimized_d13C[constraint$child] <- optimized_d13C[constraint$parent] + constraint$offset
  }
  
  # Run final model
  final_results <- run_isotope_model(
    d13C = optimized_d13C,
    d2H_Phyt.H = d2H,
    d2H_H.OH = 0,
    d2H_H.APCI = 0,
    d17O_Phyt.O = 0
  )
  
  list(
    optimized_d13C = optimized_d13C,
    theoretical_values = data.frame(
      Fragment = names(final_results$r13C_theory),
      r13C_Theory = unname(final_results$r13C_theory)
    ),
    observed_values = observed,
    constraints = constraints,
    convergence = opt_result$convergence,
    message = opt_result$message,
    gradient_used = use_analytical_gradient
  )
}


# ---- Prediction model based on regression ----
predict_intramolecular_ro <- function(df, fragment_definitions) {
  # Validate input
  required_cols <- c("Fragment","Group","m.ro.13C.corr","homogeneous.theory")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  # Step 1: Fit regression models for each Group
  group_models <- df %>%
    group_by(Group) %>%
    do(model = lm(m.ro.13C.corr ~ homogeneous.theory, data = .))
  
  # Step 2: Predict for all fragments in df
  fragments_with_predictions <- df %>%
    rowwise() %>%
    mutate(
      PredictedValue = {
        model <- group_models$model[group_models$Group == Group][[1]]
        predict(model, newdata = data.frame(homogeneous.theory = homogeneous.theory))
      }
    ) %>%
    ungroup()
  
  # Step 3: Map to carbon positions
  carbon_contributions <- list()
  
  for (i in 1:nrow(fragments_with_predictions)) {
    frag <- fragments_with_predictions$Fragment[i]
    if (frag %in% names(fragment_definitions)) {
      carbons <- fragment_definitions[[frag]]
      pred_value <- fragments_with_predictions$PredictedValue[i]
      
      for (carbon in carbons) {
        if (is.null(carbon_contributions[[carbon]])) {
          carbon_contributions[[carbon]] <- numeric(0)
        }
        carbon_contributions[[carbon]] <- c(carbon_contributions[[carbon]], pred_value)
      }
    }
  }
  
  # Step 4: Aggregate by carbon position
  carbon_isotopes <- data.frame(
    Carbon = names(carbon_contributions),
    MeanIsotope = sapply(carbon_contributions, mean),
    SDIsotope = sapply(carbon_contributions, sd),
    N_Fragments = sapply(carbon_contributions, length),
    stringsAsFactors = FALSE
  ) %>%
    arrange(Carbon)
  
  return(list(
    fragment_predictions = fragments_with_predictions,
    carbon_isotopes = carbon_isotopes
  ))
}

#---- stat model ----
perform_statistical_tests_by_ID <- function(data) {
  # Initialize list to store all results
  results_by_ID <- list()
  
  # Get unique ID groups
  id_groups <- unique(data$ID)
  
  for (id in id_groups) {
    # Filter data for current ID
    id_data <- data %>% filter(ID == id)
    
    # Skip if insufficient data
    if (nrow(id_data) < 5) {
      warning(paste("Insufficient data for ID:", id, "- Skipping tests"))
      next
    }
    
    # Initialize storage for this ID's results
    id_results <- list(
      ID = id,
      n_samples = nrow(id_data),
      models = list(),
      diagnostics = list(),
      multivariate = list()
    )
    
    # 1. Linear Model ----
    lm_formula <- m.delta.corr.t ~ m.ro.13C.corr
    lm_fit <- try(lm(lm_formula, data = id_data), silent = TRUE)
    
    if (!inherits(lm_fit, "try-error")) {
      # Store model summary
      id_results$models$linear <- summary(lm_fit)
      
      # Model diagnostics
      id_results$diagnostics <- list(
        autocorrelation = list(
          dwtest = lmtest::dwtest(lm_fit)
        ),
        heteroskedasticity = list(
          ncvTest = car::ncvTest(lm_fit),
          bptest = lmtest::bptest(lm_fit)
        ),
        normality = shapiro.test(residuals(lm_fit)),
        vif = tryCatch(car::vif(lm_fit), error = function(e) NA)
      )
    }
    
    # 2. Multivariate Analyses ----
    
    # Prepare numeric data matrix
    num_data <- id_data %>% 
      select(m.delta.corr.t, m.ro.13C.corr) %>% 
      scale() %>% 
      na.omit()
    
    # a) PCA
    pca_fit <- try(prcomp(num_data), silent = TRUE)
    if (!inherits(pca_fit, "try-error")) {
      id_results$multivariate$pca <- list(
        summary = summary(pca_fit),
        components = pca_fit$x[, 1:2],  # First two PCs
        importance = pca_fit$sdev^2 / sum(pca_fit$sdev^2)
      )
    }
    
    # b) K-means clustering (k=3)
    set.seed(123)
    km_fit <- try(kmeans(num_data, centers = 3), silent = TRUE)
    if (!inherits(km_fit, "try-error")) {
      id_results$multivariate$kmeans <- list(
        clusters = km_fit$cluster,
        centers = km_fit$centers,
        withinss = km_fit$withinss,
        tot.withinss = km_fit$tot.withinss
      )
    }
    
    # 3. Additional Models (if applicable) ----
    
    # Only run if multiple files exist for this ID
    if (length(unique(id_data$Group)) > 1) {
      # a) Multinomial model
      mn_fit <- try(
        nnet::multinom(m.ro.13C.corr ~ n.C * Group, data = id_data),
        silent = TRUE
      )
      
      if (!inherits(mn_fit, "try-error")) {
        id_results$models$multinomial <- list(
          coefficients = coef(mn_fit),
          summary = summary(mn_fit),
          AIC = AIC(mn_fit)
        )
      }
    }
    
    # Store this ID's results
    results_by_ID[[id]] <- id_results
  }
  
  # Add metadata about the analysis
  attr(results_by_ID, "analysis_date") <- Sys.Date()
  attr(results_by_ID, "variables_used") <- c("m.delta.corr.t", "m.ro.13C.corr")
  attr(results_by_ID, "n_IDs_analyzed") <- length(results_by_ID)
  
  return(results_by_ID)
}

# Usage example:
# statistical_results <- perform_statistical_tests_by_ID(results)

# To access results for a specific ID:
# statistical_results[["STD"]]  # Replace "STD" with your ID of interest

# To extract all linear model summaries:
# lapply(statistical_results, function(x) x$models$linear)