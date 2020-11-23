# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# Run source code to import/format data
source("./Code/Format data.R")

# Load additional required packages
library(lspline); library(forecast)

# Set storage directory for outputs
out <- paste0("./Results/")

## Functions -------------------------------------------------------------------

# Function to summarise the number of measures that exceed a cutoff value
# (each measure has its own specific cutoff value)
# Arguments: 
# DATA = full policies dataframe (i.e. contains all measures) for ONE country;
# MEASURES = a list of the measures to evaluate;
# CUTOFFS = a list of the cutoff points corresponding to each individual measure;
# Returns: 
# Dataframe with 2 columns: (1) Date; (2) N_measures
Summarise_N_Measures <- function(DATA, MEASURES, CUTOFFS) {
  
  # Print warning if each measure doesn't have its own cutoff point
  if (length(MEASURES) != length(CUTOFFS)) {
    print("Error: Each measure must have its own cutoff point.")
  }
  
  # Subset dataset to include only selected Measures
  data <- DATA %>% ungroup() %>% select(Date, unlist(MEASURES))
  
  # For each measure, determine whether it is greater than or equal to its cutoff
  # If yes, assign 1; if no, assign 0
  for (i in 1:length(MEASURES)) {
    measure <- MEASURES[[i]]
    data[, measure] <- ifelse(data[, measure] >= CUTOFFS[[i]], 1, 0)
  }
  
  # Calculate number of measures above cutoff (i.e. row sums)
  data$N_measures <- apply(data[, unlist(MEASURES)], 1, sum, na.rm = TRUE)
  
  # Retain Date, N_measures variables
  data <- data %>% select(Date, N_measures)
  
  # Return dataframe with Date, N_measures variables
  return(data)
  
}

# Function to calculate Poisson deviance between two vectors
# from: https://en.wikipedia.org/wiki/Deviance_(statistics)
# Arguments: obs = vector of observed values, sim = vector of simulated/predicted values
Calc_Pois_Dev <- function(obs, sim) {
  
  #D <- 2 * sum(obs * log(obs / sim) - (obs - sim))
  D <- 2 * sum(obs * log(obs / sim) - (obs - sim), na.rm = TRUE)
  return(D)
  
}

# ------------------------------------------------------------------------------
# Summarise countries 
# ------------------------------------------------------------------------------

# Create summary table to store info on:
# Max number of restrictions,
# Important dates: Date_0, Date_pop_pct, Date_max, Date_first_restriction, 
# Date_restrictions_eased, Date_lockdown, Date_lockdown_eased, Date_lockdown_end
summary_eur <- data_eur %>% select(Country, Date_0, Date_pop_pct, Date_max) %>% unique %>%
  mutate(Date_first_restriction = as.Date(NA),
         Date_restrictions_eased = as.Date(NA),
         Date_lockdown = as.Date(NA),
         Date_lockdown_eased = as.Date(NA),
         Date_lockdown_end = as.Date(NA),
         Max_number_restrictions = as.numeric(NA))

# Define requirements of interest and cutoff points for whether measures have been implemented
# [1 corresponds to measure recommended, 2-3 corresponds to measure required]
## (1) Lockdown measure
measures_lockdown <- list("C6_Stay_at_home_requirements")
cutoffs_lockdown <- list(2)
#length(measures_lockdown) == length(cutoffs_lockdown)
## (2) Alternate group of measures which together are broadly equivalent to lockdown
measures_lockdown_alt <- list("C1_School_closing", "C2_Workplace_closing", "C3_Cancel_public_events",
                              "C5_Close_public_transport", "C7_Restrictions_on_internal_movement")
cutoffs_lockdown_alt <- list(2, 2, 2, 2, 2)
#length(measures_lockdown_alt) == length(cutoffs_lockdown_alt)
## (3) Full group of measures (i.e. any restrictions)
measures_any_restriction <- as.list(sort(unlist(c(measures_lockdown, measures_lockdown_alt))))
cutoffs_any_restriction <- list(1, 1, 1, 1, 1, 1)
#length(measures_any_restriction) == length(cutoffs_any_restriction)

# Define threshold values, i.e. how many measures together constitute a given restriction
threshold_lockdown <- 1
threshold_lockdown_alt <- 3
threshold_any_restriction <- 1

# Record country summary data
for (i in 1:nrow(summary_eur)) {
  
  # Define country
  country <- summary_eur[[i, "Country"]] %>% as.character()
  
  # Filter cases/deaths and policies datasets by country
  policies_eur_i <- policies_eur %>% filter(Country == country)
  data_eur_i <- data_eur %>% filter(Country == country)
  
  # Define date_max (i.e. last date for which data can be reasonably assumed complete)
  date_max <- data_eur_i %>% pull(Date_max) %>% head(1)
  
  # Filter policies data by date_max
  policies_eur_i <- policies_eur_i %>% filter(Date <= date_max)
  
  # Restrictions  --------
  
  # (either recommended or required)
  
  # Determine number of recommended or required measures values on each date
  data_any_restriction <- Summarise_N_Measures(DATA = policies_eur_i,
                                               MEASURES = measures_any_restriction,
                                               CUTOFFS = cutoffs_any_restriction)
  
  # Determine change in number of restrictions from previous day
  # (negative value indicates restriction(s) lifted, ...
  # ...zero indicates no change, positive value indicates restriction(s) increased)
  data_any_restriction <- data_any_restriction %>%
    mutate(N_measures_diff = N_measures - lag(N_measures, n = 1, default = NA))
  
  # Determine whether any measures were recommended/required (i.e. threshold of 1 exceeded) for each date,
  # and whether number of measures decreases (i.e. diff is negative) for each date
  data_any_restriction <- data_any_restriction %>% 
    mutate(Threshold_exceeded = ifelse(N_measures >= threshold_any_restriction, TRUE, FALSE),
           Diff_neg = ifelse(N_measures_diff < 0, TRUE, FALSE)) %>%
    relocate(-contains("diff"))
  
  # Create indicator for whether any restriction was recommended/required
  any_restriction_tf <- any(data_any_restriction$Threshold_exceeded == TRUE, na.rm = TRUE)
  
  # Analyse countries which recommended/required any restriction
  if (any_restriction_tf == TRUE) {
    
    ## Max number of restrictions ##
    
    # Determine max number of recommended or required measures
    max_number_restrictions <- data_any_restriction %>% summarise(max(N_measures)) %>% pull()
    # Add max number of measures to summary table
    summary_eur[[i, "Max_number_restrictions"]] <- max_number_restrictions
    
    ## First restriction ##
    ## (first date where any restrictions were either recommended or required) ##
    
    # Filter restrictions data by dates where any measures were recommended
    data_any_restriction_filt <- data_any_restriction %>% filter(Threshold_exceeded == TRUE)
    # Define date of first restriction as first date where threshold was exceeded
    date_first_restriction <- data_any_restriction_filt %>% pull(Date) %>% min(na.rm = TRUE)
    # Add to summary table
    summary_eur[[i, "Date_first_restriction"]] <- date_first_restriction
    
    ## Restrictions eased ##
    ## (first date after date of first restriction where number of measures decreased) ##
    
    # Filter restrictions data by dates greater than date of first restriction
    data_any_restriction_filt <- data_any_restriction %>% filter(Date >= date_first_restriction)
    
    # Create indicator for whether restrictions were eased
    # (TRUE if difference is negative on any date, FALSE if difference always >= 0)
    restrictions_eased_tf <- any(data_any_restriction_filt$Diff_neg == TRUE, na.rm = TRUE)
    
    # Analyse countries whose restrictions were eased
    if (restrictions_eased_tf == TRUE) {
      # Re-filter data by dates where difference is negative
      data_any_restriction_filt <- data_any_restriction_filt %>% filter(Diff_neg == TRUE)
      # Define date of restrictions eased as first date with negative difference
      date_restrictions_eased <- data_any_restriction_filt %>% pull(Date) %>% min(na.rm = TRUE)
      # Add to summary table
      summary_eur[[i, "Date_restrictions_eased"]] <- date_restrictions_eased
    }  
    
  }  # (close inner loop - countries which which recommended/required any restriction)
  
  # Summarise lockdown --------
  
  # Determine number of the following required measures on each date
  # (1) lockdown (general or targeted)
  data_lockdown <- Summarise_N_Measures(DATA = policies_eur_i,
                                        MEASURES = measures_lockdown,
                                        CUTOFFS = cutoffs_lockdown)
  # (2) alternate lockdown measures
  data_lockdown_alt <- Summarise_N_Measures(DATA = policies_eur_i, 
                                            MEASURES = measures_lockdown_alt, 
                                            CUTOFFS = cutoffs_lockdown_alt)
  
  # Determine change in number of restrictions from previous day
  # (negative value indicates restriction(s) lifted, ...
  # ...zero indicates no change, positive value indicates restriction(s) increased)
  data_lockdown <- data_lockdown %>% 
    mutate(N_measures_diff = N_measures - lag(N_measures, n = 1, default = NA))
  data_lockdown_alt <- data_lockdown_alt %>% 
    mutate(N_measures_diff = N_measures - lag(N_measures, n = 1, default = NA))
  
  # Merge two dataframes together
  data_lockdown_all <- full_join(data_lockdown, data_lockdown_alt, by = "Date", suffix = c("_lockdown", "_lockdown_alt"))
  
  # Determine whether lockdown measures or alternate lockdown threshold is exceeded for each date,
  # and whether number of lockdown measures or alternate lockdown measures decreases (i.e. diff is negative) for each date
  data_lockdown_all <- data_lockdown_all %>% 
    mutate(Threshold_exceeded = ifelse(N_measures_lockdown >= threshold_lockdown |
                                         N_measures_lockdown_alt >= threshold_lockdown_alt,
                                       TRUE, FALSE),
           Diff_neg = ifelse(N_measures_diff_lockdown < 0 | N_measures_diff_lockdown_alt < 0, 
                             TRUE, FALSE)) %>%
    relocate(-contains("diff"))
  
  # Create indicator for whether lockdown was entered
  # (TRUE if threshold exceeded on any date, FALSE if threshold never exceeded)
  lockdown_tf <- any(data_lockdown_all$Threshold_exceeded == TRUE, na.rm = TRUE)
  
  # Analyse countries which entered lockdown
  if (lockdown_tf == TRUE) {
    
    ## Beginning of lockdown ##
    ## (first date where either lockdown or alternate lockdown measures were required) ##
    
    # Filter lockdown data by dates where lockdown threshold was exceeded
    data_lockdown_all_filt <- data_lockdown_all %>% filter(Threshold_exceeded == TRUE)
    # Define lockdown date as first date where threshold was exceeded
    date_lockdown <- data_lockdown_all_filt %>% pull(Date) %>% min(na.rm = TRUE)
    # Add to summary table
    summary_eur[[i, "Date_lockdown"]] <- date_lockdown
    
    ## Lockdown eased ##
    ## (first date after lockdown date where number of measures decreased) ##
    
    # Filter lockdown data by dates greater than lockdown date
    data_lockdown_all_filt <- data_lockdown_all %>% filter(Date >= date_lockdown)
    
    # Create indicator for whether lockdown was eased
    # (TRUE if difference is negative on any date, FALSE if difference always >= 0)
    lockdown_eased_tf <- any(data_lockdown_all_filt$Diff_neg == TRUE, na.rm = TRUE)
    
    # Analyse countries which eased lockdown
    if (lockdown_eased_tf == TRUE) {
      # Re-filter data by dates where difference is negative
      data_lockdown_all_filt <- data_lockdown_all_filt %>% filter(Diff_neg == TRUE)
      # Define lockdown easing date as first date with negative difference
      date_lockdown_eased <- data_lockdown_all_filt %>% pull(Date) %>% min(na.rm = TRUE)
      # Add to summary table
      summary_eur[[i, "Date_lockdown_eased"]] <- date_lockdown_eased
    }
    
    ## End of lockdown ##
    ## (first date AFTER lockdown date where neither lockdown measures no alt lockdown measures were required) ##
    
    # Filter lockdown data by dates greater than lockdown date
    data_lockdown_all_filt <- data_lockdown_all %>% filter(Date >= date_lockdown)
    
    # Create indicator for whether lockdown was ended
    lockdown_end_tf <- any(data_lockdown_all_filt$Threshold_exceeded == FALSE, na.rm = TRUE)
    
    # Analyse countries which ended lockdown
    if (lockdown_end_tf == TRUE) {
      # Re-filter data by dates where lockdown threshold wasn't exceeded
      data_lockdown_all_filt <- data_lockdown_all_filt %>% filter(Threshold_exceeded == FALSE) 
      # Define lockdown end date
      date_lockdown_end <- data_lockdown_all_filt %>% pull(Date) %>% min(na.rm = TRUE)
      # Add to summary table
      summary_eur[[i, "Date_lockdown_end"]] <- date_lockdown_end
    }
    
  }  # (close inner loop - countries which entered lockdown)
  
}  # (close outer loop - all countries)

# Remove measures/cutoffs, loop objects
rm(measures_lockdown, measures_lockdown_alt, measures_any_restriction,
   cutoffs_lockdown, cutoffs_lockdown_alt, cutoffs_any_restriction,
   threshold_lockdown, threshold_lockdown_alt, threshold_any_restriction)
rm(i, country, data_eur_i, policies_eur_i, date_max,
   any_restriction_tf, restrictions_eased_tf, 
   lockdown_tf, lockdown_eased_tf, lockdown_end_tf,
   data_any_restriction, data_any_restriction_filt, 
   data_lockdown, data_lockdown_alt, data_lockdown_all, data_lockdown_all_filt,
   max_number_restrictions, date_first_restriction, date_restrictions_eased,
   date_lockdown, date_lockdown_eased, date_lockdown_end)

# Calculate date_T (last date to include data from) as either...
# Date_max, Date_restrictions_eased + 28, or Date_lockdown_eased + 28, whichever comes first
summary_eur <- summary_eur %>% group_by(Country) %>%
  mutate(Date_T = min(Date_max, Date_restrictions_eased + 28, Date_lockdown_eased + 28, na.rm = TRUE)) %>%
  relocate(Date_T, .before = Max_number_restrictions) %>% ungroup

# Determine European countries which entered lockdown
countries_eur_lockdown <- summary_eur %>% filter(!is.na(Date_lockdown)) %>% 
  pull(Country) %>% as.character %>% as.list
## print note about any countries which did not enter lockdown:
if (length(countries_eur_lockdown) != length(countries_eur)) {
  unavail <- setdiff(unlist(countries_eur), unlist(countries_eur_lockdown))
  cat(paste0("Note that the following countries did not enter lockdown:\n", 
             paste0(unavail, collapse = ", ")))
  rm(unavail)
} 

# Export summary table
write_csv(summary_eur, file = paste0(out, "Country summaries.csv"))

# ------------------------------------------------------------------------------
# Estimate knot points
# ------------------------------------------------------------------------------

# Create list to store best knots for each country
knots_best <- list()

# (1) Iterate through countries
start <- Sys.time()
for (i in countries_eur) {
  
  # Define country
  country <- i
  
  # Filter cases/deaths, summary dataframes by country
  data_eur_i <- data_eur %>% filter(Country == country)
  summary_eur_i <- summary_eur %>% filter(Country == country)
  
  # Record important dates
  date_pop_pct <- summary_eur_i %>% pull(Date_pop_pct)
  date_first_restriction <- summary_eur_i %>% pull(Date_first_restriction)
  date_lockdown <- summary_eur_i %>% pull(Date_lockdown)
  date_T <- summary_eur_i %>% pull(Date_T)
  
  # Create copy of cases/deaths dataframe where 
  # cumulative cases >= population percentage threshold and up to date_T
  data_eur_pop_pct_i <- data_eur_i %>% filter(Date >= date_pop_pct & Date <= date_T)
  
  # Create indicator to skip to next iteration (country)
  skip_to_next_i <- FALSE
  
  # Define potential knot dates (from dates of first restriction and lockdown to 28 days subsequent),
  # And create grid of all possible combinations of knot dates, with restrictions that...
  # (a) first knot date must be before or at the same time as second knot date, and
  # (b) knot dates must fall within modelling period 
  # (i.e. after the first date at which cumulative cases >= pop pct threshold and less than date_T)
  if (is.na(date_lockdown) | date_first_restriction == date_lockdown) {
    possible_knot_dates_1 <- seq(from = date_first_restriction, to = date_first_restriction + 28, by = 1)
    grid <- tibble("Knot_date_1" = possible_knot_dates_1) %>% 
      filter(Knot_date_1 >= date_pop_pct, Knot_date_1 < date_T)
  } else {
    possible_knot_dates_1 <- seq(from = date_first_restriction, to = date_first_restriction + 28, by = 1)
    possible_knot_dates_2 <- seq(from = date_lockdown, to = date_lockdown + 28, by = 1)
    grid <- tibble(expand.grid(possible_knot_dates_2, possible_knot_dates_1))
    names(grid) <- c("Knot_date_2", "Knot_date_1")
    grid <- grid %>% select("Knot_date_1", "Knot_date_2") %>% 
      filter(Knot_date_1 <= Knot_date_2, Knot_date_1 >= date_pop_pct, Knot_date_2 < date_T)  
    # If first knot date equals second knot date, replace second with NA
    for (g in 1:nrow(grid)) {
      k_1 <- grid[[g, "Knot_date_1"]]
      k_2 <- grid[[g, "Knot_date_2"]]
      if (k_1 == k_2) {grid[[g, "Knot_date_2"]] <- NA}
    }
  }
  
  # If no rows in grid, skip to next country 
  # (modelling is not possible for current country, 
  # likely because restrictions happened long before threshold for modelling is reached)
  if (nrow(grid) == 0) { skip_to_next_i <- TRUE }
  if (skip_to_next_i) { next }
  
  # Create dataframe to store summary statistics for all possible combinations of knot dates
  knots <- bind_rows(tibble(Knot_date_1 = as.Date(character()),
                            Knot_date_2 = as.Date(character()),
                            N_knots = as.numeric(),
                            Growth_factor_1 = as.numeric(),
                            Growth_factor_2 = as.numeric(),
                            Growth_factor_3 = as.numeric(),
                            Growth_factor_1_sd = as.numeric(),
                            Growth_factor_2_sd = as.numeric(),
                            Growth_factor_3_sd = as.numeric(),
                            Intercept_1 = as.numeric(),
                            Intercept_2 = as.numeric(),
                            Intercept_3 = as.numeric(),
                            Pois_dev_inc = as.numeric(),
                            Pois_dev_cum = as.numeric()),
                     grid)
  rm(grid)
  
  # Set dates over which to simulate growth
  dates <- seq.Date(from = date_pop_pct, to = date_T, by = 1)
  
  # Create matrices for simulated data (daily and cumulative cases)
  # (1 row per simulation run, 1 col per date)
  daily_cases_sim <- cumulative_cases_end_sim <- 
    matrix(nrow = 1, ncol = length(dates) + 1,
           dimnames = list(1, as.character(seq.Date(from = date_pop_pct - 1, to = date_T, by = 1))))
  # Initialise matrices with data at date_pop_pct - 1
  daily_cases_sim[, 1] <- data_eur_i %>% 
    filter(Date == (date_pop_pct - 1)) %>% pull(Daily_cases_MA7)
  cumulative_cases_end_sim[, 1] <- data_eur_i %>% 
    filter(Date == (date_pop_pct - 1)) %>% pull(Cumulative_cases_end_MA7)
  
  # (2) Iterate through pairs of candidate knot points
  for (j in 1:nrow(knots)) {
    
    # Set knot dates
    knot_date_1 <- knots[[j, "Knot_date_1"]]
    knot_date_2 <- knots[[j, "Knot_date_2"]]
    
    # Create indicator to skip to next iteration (pair of knot dates)
    skip_to_next_j <- FALSE
    
    # Estimate growth parameters
    ## If first knot occurs at first date for which cases exceeded pop pct threshold (i.e. when we begin modelling),
    ## there may be either no knots (i.e. knot occured before or at date_pop_pct)
    ## OR 1 knot (occurring at knot_date_2).
    ## Otherwise, there may be either 1 knot (occurring at knot_date_1 (= knot_date_2, if it exists))
    ## OR 2 knots (occurring at knot_date_1 and knot_date_2)
    if (knot_date_1 == date_pop_pct) {
      
      if (is.na(knot_date_2)) {  # NO knot points
        
        # Set number of knot points and record
        n_knots <- 0
        knots[[j, "N_knots"]] <- n_knots
        
        # Fit regular Arima model (with intercept, since this is not technically first segment)
        model <- tryCatch(Arima(data_eur_pop_pct_i$Daily_cases, order = c(2, 0, 0), 
                                seasonal = list(order = c(1, 0, 0), period = 7),
                                xreg = as.matrix(data_eur_pop_pct_i[, "Cumulative_cases_beg"]), 
                                include.constant = TRUE, method = "ML"), 
                          error = function(e) { skip_to_next_j <<- TRUE } )
        if (skip_to_next_j) { next }
        
        # Record model parameters (intercept, slope, and SD of slope)
        intercept_1 <- as.numeric(coef(model)["intercept"])
        slope_1 <- as.numeric(coef(model)["Cumulative_cases_beg"])
        slope_1_sd <- sqrt(diag(model$var.coef))[["Cumulative_cases_beg"]]
        
        # Calculate and record growth factor, record model parameters
        knots[[j, "Growth_factor_1"]] <- growth_factor_1 <- slope_1 + 1
        knots[[j, "Growth_factor_1_sd"]] <- slope_1_sd
        knots[[j, "Intercept_1"]] <- intercept_1
        
      } else {  # ONE knot point (at knot_date_2)
        
        # Set number of knot points and record
        n_knots <- 1
        knots[[j, "N_knots"]] <- n_knots
        
        # Set knot point
        knot_1 <- data_eur_pop_pct_i %>% filter(Date == knot_date_2) %>% pull(Cumulative_cases_beg)
        
        # Create dataframe for fitting manual splines
        data_j <- data.frame(lspline(data_eur_pop_pct_i$Cumulative_cases_beg, knots = c(knot_1)))
        names(data_j) <- names <- paste0("Cumulative_cases_beg_", 1:2)
        data_j <- bind_cols(Daily_cases = data_eur_pop_pct_i$Daily_cases, data_j)
        
        # Fit ARIMA spline model w/ specified knot point (with intercept)
        model <- tryCatch(Arima(data_j$Daily_cases, order = c(2, 0, 0), 
                                seasonal = list(order = c(1, 0, 0), period = 7),
                                xreg = as.matrix(data_j[, names]), 
                                include.constant = TRUE, method = "ML"), 
                          error = function(e) { skip_to_next_j <<- TRUE } )
        if (skip_to_next_j) { next }
        
        # Record model parameters (intercept, slope, and SD of slope)
        intercept_1 <- as.numeric(coef(model)["intercept"])
        slope_1 <- as.numeric(coef(model)["Cumulative_cases_beg_1"])
        slope_2 <- as.numeric(coef(model)["Cumulative_cases_beg_2"])
        slope_1_sd <- sqrt(diag(model$var.coef))[["Cumulative_cases_beg_1"]]
        slope_2_sd <- sqrt(diag(model$var.coef))[["Cumulative_cases_beg_2"]]
        intercept_2 <- (intercept_1 + slope_1*knot_1) - slope_2*knot_1
        
        # Calculate and record growth factor, record model parameters
        knots[[j, "Growth_factor_1"]] <- growth_factor_1 <- slope_1 + 1
        knots[[j, "Growth_factor_2"]] <- growth_factor_2 <- slope_2 + 1
        knots[[j, "Growth_factor_1_sd"]] <- slope_1_sd
        knots[[j, "Growth_factor_2_sd"]] <- slope_2_sd
        knots[[j, "Intercept_1"]] <- intercept_1
        knots[[j, "Intercept_2"]] <- intercept_2
        
      }
      
    } else {
      
      if (is.na(knot_date_2)) {  # ONE knot point (at knot_date_1)
        
        # Set number of knot points and record
        n_knots <- 1
        knots[[j, "N_knots"]] <- n_knots
        
        # Set knot point
        knot_1 <- data_eur_pop_pct_i %>% filter(Date == knot_date_1) %>% pull(Cumulative_cases_beg)
        
        # Create dataframe for fitting manual splines
        data_j <- data.frame(lspline(data_eur_pop_pct_i$Cumulative_cases_beg, knots = c(knot_1)))
        names(data_j) <- names <- paste0("Cumulative_cases_beg_", 1:2)
        data_j <- bind_cols(Daily_cases = data_eur_pop_pct_i$Daily_cases, data_j)
        
        # Fit ARIMA spline model w/ specified knot point (with intercept)
        model <- tryCatch(Arima(data_j$Daily_cases, order = c(2, 0, 0), 
                                seasonal = list(order = c(1, 0, 0), period = 7),
                                xreg = as.matrix(data_j[, names]), 
                                include.constant = TRUE, method = "ML"), 
                          error = function(e) { skip_to_next_j <<- TRUE } )
        if (skip_to_next_j) { next }
        
        # Record model parameters (intercept, slope, and SD of slope)
        intercept_1 <- as.numeric(coef(model)["intercept"])
        slope_1 <- as.numeric(coef(model)["Cumulative_cases_beg_1"])
        slope_2 <- as.numeric(coef(model)["Cumulative_cases_beg_2"])
        slope_1_sd <- sqrt(diag(model$var.coef))[["Cumulative_cases_beg_1"]]
        slope_2_sd <- sqrt(diag(model$var.coef))[["Cumulative_cases_beg_2"]]
        intercept_2 <- (intercept_1 + slope_1*knot_1) - slope_2*knot_1
        
        # Calculate and record growth factor, record model parameters
        knots[[j, "Growth_factor_1"]] <- growth_factor_1 <- slope_1 + 1
        knots[[j, "Growth_factor_2"]] <- growth_factor_2 <- slope_2 + 1
        knots[[j, "Growth_factor_1_sd"]] <- slope_1_sd
        knots[[j, "Growth_factor_2_sd"]] <- slope_2_sd
        knots[[j, "Intercept_1"]] <- intercept_1
        knots[[j, "Intercept_2"]] <- intercept_2
        
      } else {  # TWO knot points (at knot_date_1 and knot_date_2)
        
        # Set number of knot points and record
        n_knots <- 2
        knots[[j, "N_knots"]] <- n_knots
        
        # Set knot points
        knot_1 <- data_eur_pop_pct_i %>% filter(Date == knot_date_1) %>% pull(Cumulative_cases_beg)
        knot_2 <- data_eur_pop_pct_i %>% filter(Date == knot_date_2) %>% pull(Cumulative_cases_beg)
        
        # Create dataframe for fitting manual splines
        data_j <- data.frame(lspline(data_eur_pop_pct_i$Cumulative_cases_beg, knots = c(knot_1, knot_2)))
        names(data_j) <- names <- paste0("Cumulative_cases_beg_", 1:3)
        data_j <- bind_cols(Daily_cases = data_eur_pop_pct_i$Daily_cases, data_j)
        
        # Fit ARIMA spline model w/ specified knot points (no intercept)
        model <- tryCatch(Arima(data_j$Daily_cases, order = c(2, 0, 0), 
                                seasonal = list(order = c(1, 0, 0), period = 7),
                                xreg = as.matrix(data_j[, names]), 
                                include.constant = FALSE, method = "ML"), 
                          error = function(e) { skip_to_next_j <<- TRUE } )
        if (skip_to_next_j) { next }
        
        # Record model parameters (intercept, slope, and SD of slope)
        intercept_1 <- 0
        slope_1 <- as.numeric(coef(model)["Cumulative_cases_beg_1"])
        slope_2 <- as.numeric(coef(model)["Cumulative_cases_beg_2"])
        slope_3 <- as.numeric(coef(model)["Cumulative_cases_beg_3"])
        slope_1_sd <- sqrt(diag(model$var.coef))[["Cumulative_cases_beg_1"]]
        slope_2_sd <- sqrt(diag(model$var.coef))[["Cumulative_cases_beg_2"]]
        slope_3_sd <- sqrt(diag(model$var.coef))[["Cumulative_cases_beg_3"]]
        intercept_2 <- (intercept_1 + slope_1*knot_1) - slope_2*knot_1
        intercept_3 <- (intercept_2 + slope_2*knot_2) - slope_3*knot_2
        
        # Calculate and record growth factor, record model parameters
        knots[[j, "Growth_factor_1"]] <- growth_factor_1 <- slope_1 + 1
        knots[[j, "Growth_factor_2"]] <- growth_factor_2 <- slope_2 + 1
        knots[[j, "Growth_factor_3"]] <- growth_factor_3 <- slope_3 + 1
        knots[[j, "Growth_factor_1_sd"]] <- slope_1_sd
        knots[[j, "Growth_factor_2_sd"]] <- slope_2_sd
        knots[[j, "Growth_factor_3_sd"]] <- slope_3_sd
        knots[[j, "Intercept_1"]] <- intercept_1
        knots[[j, "Intercept_2"]] <- intercept_2
        knots[[j, "Intercept_3"]] <- intercept_3
        
      }
      
    }  # (close if-else section)
    
    # Record model summaries
    knots[[j, "BIC"]] <- BIC(model)
    
    # (3) Estimate growth of cases using knot point(s)
    for (t in as.list(dates)) {
      
      # Get daily and cumulative cases from time t-1
      inc_tminus1 <- daily_cases_sim[, as.character(t-1)]
      cum_tminus1 <- cumulative_cases_end_sim[, as.character(t-1)]
      
      # Define growth parameters
      if (n_knots == 0) {  # NO knot points
        growth <- growth_factor_1  
      } else if (n_knots == 1) {  # ONE knot point
        if (t <= knot_date_1) {
          growth <- growth_factor_1
        } else {
          growth <- growth_factor_2
        }
      } else {  # TWO knot points
        if (t <= knot_date_1) {
          growth <- growth_factor_1
        } else if (t <= knot_date_2) {
          growth <- growth_factor_2
        } else {
          growth <- growth_factor_3
        }
      }
      
      # Calculate daily cases at time t and record
      inc_t <- growth*inc_tminus1
      daily_cases_sim[, as.character(t)] <- inc_t
      
      # Calculate cumulative cases at end of time t and record
      cum_t <- cum_tminus1 + inc_t
      cumulative_cases_end_sim[, as.character(t)] <- cum_t
      
    }  # (close loop 3)
    
    # Calculate and record Poisson deviance
    ## (1) For predicted vs true (7-day moving average) incident cases
    true_inc <- data_eur_pop_pct_i$Daily_cases_MA7
    pred_inc <- daily_cases_sim[1, -1]
    knots[[j, "Pois_dev_inc"]] <- Calc_Pois_Dev(obs = true_inc, sim = pred_inc)
    ## (2) For true vs predicted cumulative cases
    true_cum <- data_eur_pop_pct_i$Cumulative_cases_end_MA7
    pred_cum <- cumulative_cases_end_sim[1, -1]
    knots[[j, "Pois_dev_cum"]] <- Calc_Pois_Dev(obs = true_cum, sim = pred_cum)
    
    # Display progress 
    cat('\r', paste(round((j / nrow(knots) * 100), 0), 
                    "% done of country", grep(country, unlist(countries_eur)), "of", 
                    length(countries_eur), "          ", sep = " "))
    
  }  # (close loop 2)
  
  # Remove from consideration knot date pairs for which ...
  # (1) growth factor 1 is less than 1 AND growth factor 3 exists 
  # (because where there are 3 segments, the first scenario must represent initial uncontrolled growth)
  # (2) any of growth factors are negative
  remove_1 <- knots %>% filter(Growth_factor_1 < 1 & !is.na(Growth_factor_3)) 
  remove_2 <- knots %>% filter(Growth_factor_1 < 0 | Growth_factor_2 < 0 | Growth_factor_3 < 0)
  knots <- anti_join(knots, remove_1, by = names(knots)) %>% 
    anti_join(., remove_2, by = names(knots))
  
  # Find best knot points (by lowest Pois_dev_inc) for each country and label 
  knots_best_i <- knots %>% arrange(Pois_dev_inc) %>% head(10) %>% 
    mutate(Country = country) %>% relocate(Country)
  
  # Add best knots for country i to list of best knots
  knots_best[[i]] <- knots_best_i
  
}  # (close loop 1)
end <- Sys.time()
end - start  # ~9 mins

# Remove loop variables
rm(i, j, t, g, country, data_eur_i, summary_eur_i, data_eur_pop_pct_i, 
   date_pop_pct, date_first_restriction, date_lockdown, date_T,
   possible_knot_dates_1, possible_knot_dates_2, k_1, k_2, knots, dates, 
   daily_cases_sim, cumulative_cases_end_sim, knot_date_1, knot_date_2,
   skip_to_next_i, skip_to_next_j, names, n_knots, knot_1, knot_2, data_j, model, 
   intercept_1, intercept_2, intercept_3,
   slope_1, slope_2, slope_3, slope_1_sd, slope_2_sd, slope_3_sd,
   growth_factor_1, growth_factor_2, growth_factor_3,
   inc_tminus1, cum_tminus1, inc_t, cum_t, growth,
   true_inc, pred_inc, true_cum, pred_cum, 
   knots_best_i, remove_1, remove_2, start, end)

# Combine best knots from all countries into single dataframe
knots_best <- bind_rows(knots_best)

# Group dataframe of best knots by country, arrange by knot dates
knots_best <- knots_best %>% group_by(Country) %>% arrange(Country, Knot_date_1, Knot_date_2)

# Construct (unequal) probability of each pair of knot points according to Pois_dev_cum
### (create inverse of Pois_dev_cum values so that lower values are ranked higher,
### calculate normaliser for rescaling Pois_dev_cum inverse values, and
### calculate probability by multiplying Pois_dev_cum inverse values by normaliser)
knots_best <- knots_best %>% mutate(Pois_dev_cum_inv = 1 / Pois_dev_cum,
                                    Norm = 1 / sum(Pois_dev_cum_inv), 
                                    Prob_unequal = Pois_dev_cum_inv * Norm) %>%
  select(-c(Pois_dev_cum_inv, Norm))

# Calculate minimum number of simulation runs per knot date required to retain relative probability of each
### (find the probability of the least likely knot date,
### calculate the multiplier required for this knot date to be simulated once, and
### multiply the probability of each knot date by multiplier)
knots_best <- knots_best %>% mutate(Prob_min = min(Prob_unequal),
                           Mult = ceiling(1 / Prob_min),
                           N_unequal = round(Prob_unequal * Mult)) %>%
  select(-c(Prob_min, Mult))

# Export knots_best dataframe
write_csv(knots_best, file = paste0(out, "Best knot points.csv"))
