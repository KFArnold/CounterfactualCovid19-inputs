# ------------------------------------------------------------------------------
# Set up
# ------------------------------------------------------------------------------

# Restore package library to last snapshot
packrat::restore()

# Load required packages
library(tidyverse); library(scales); library(ggpubr)

# Define storage directory for formatted data
data_directory_f <- paste0("./Data/Formatted/")

# Define storage directory for results
results_directory <- paste0("./Results/")

# Load formatted data
data_eur <- read_csv(paste0(data_directory_f, "Cases_deaths_data_europe.csv"))

# Import files containing best knot point pairs and country summaries
knots_best <- read_csv(paste0(results_directory, "Best knot points.csv"))
summary_eur <- read_csv(paste0(results_directory, "Country summaries.csv"))

# Load list of European countries for which we have both cases/deaths data and policy data
load(paste0(results_directory, "countries_eur.RData"))

# ------------------------------------------------------------------------------
# Simulation
# ------------------------------------------------------------------------------

# Countries tested so far: UK, Germany, Austria, Croatia, Italy
# Albania only works for natural history because one of 'best' knot points has zero knots - knot occurred before or at time we started modelling (don't know which)
# Not yet adapted code for countries where all best pairs have only 1 knot (e.g. Andorra)

# Define country
country <- "Germany"

# Filter datasets by country
knots_best_country <- knots_best %>% filter(Country == country)  # best knots
data_eur_country <- data_eur %>% filter(Country == country)  # cases/deaths data 
summary_eur_country <- summary_eur %>% filter(Country == country)

# Record important dates
date_pop_pct <- summary_eur_country %>% pull(Date_pop_pct)
date_T <- summary_eur_country %>% pull(Date_T)
date_first_restriction <- summary_eur_country %>% pull(Date_first_restriction)
date_lockdown <- summary_eur_country %>% pull(Date_lockdown)

# Calculate total number of simulation runs
n_runs <- knots_best_country %>% pull(N_unequal) %>% sum

# Needs to be adapted for only one knot point!!


# Calculate maximum number of days earlier we can estimate FIRST RESTRICTION
# (minimum value of knot_date_1 greater than or equal to date_pop_pct)
max_days_counterfactual_first_restriction <- knots_best_country %>% 
  mutate(Diff = Knot_date_1 - date_pop_pct) %>% pull(Diff) %>% min %>% as.numeric
# Calculate maximum number of days earlier we can estimate LOCKDOWN
# (minimum value of knot_date_2 must be greater than or equal to date_pop_pct, and greater than than knot_date_1)
max_days_counterfactual_lockdown <- knots_best_country %>% 
  mutate(Diff = Knot_date_2 - (date_pop_pct + 1)) %>% pull(Diff) %>% min(na.rm = TRUE) %>% as.numeric

# Determine all possible combinations of counterfactual days
# (filter so that minimum value of knot_date_2 is always greater than knot_date_1)
possible_days_counterfactual <- expand_grid(Poss_days_counterfactual_first_restriction = seq(0, max_days_counterfactual_first_restriction, 1),
                                            Poss_days_counterfactual_lockdown = seq(0, max_days_counterfactual_lockdown, 1)) %>% 
  filter((Poss_days_counterfactual_lockdown - Poss_days_counterfactual_first_restriction) <=
           (max_days_counterfactual_lockdown - max_days_counterfactual_first_restriction - 1))

# Set dates to simulate
dates <- seq.Date(from = date_pop_pct, to = date_T, by = 1)

# Define number of best knot point pairs
n_knots_best <- nrow(knots_best_country)

# Define number of days earlier to implement FIRST RESTRICTION and LOCKDOWN
# (0 for both indicates natural history, i.e. both first restrictins and lockdown as implemented)
n_days_counterfactual_first_restriction <- 0
n_days_counterfactual_lockdown <- 0
n_days_counterfactual <- tibble(Poss_days_counterfactual_first_restriction = n_days_counterfactual_first_restriction,
                                Poss_days_counterfactual_lockdown = n_days_counterfactual_lockdown)

# Determine if specified counterfactual conditions are possible,
# and print warning if not
match <- plyr::match_df(possible_days_counterfactual, n_days_counterfactual,
                  on = c("Poss_days_counterfactual_first_restriction", "Poss_days_counterfactual_lockdown"))
if (nrow(match) == 0) {
  print("Stop - cannot estimate counterfactual")
}

# Create best knots dataframe for counterfactual scenario
knots_best_country_counterfactual <- knots_best_country %>%
  mutate(Knot_date_1 = Knot_date_1 - n_days_counterfactual_first_restriction,
         Knot_date_2 = Knot_date_2 - n_days_counterfactual_lockdown)

# Create empty matrices for simulated data
# (1 row per simulation run, 1 col per date)
daily_cases_sim <- matrix(nrow = 0, ncol = length(dates) + 1, 
                          dimnames = list(NULL, as.character(seq.Date(from = date_pop_pct - 1, to = date_T, by = 1))))

## Create matrix for estimated deaths
## (1 row per simulation run, 1 col per CFR)
#deaths_sim <- matrix(nrow = n_runs, ncol = 2,
#                       dimnames = list(NULL, c("cfr_hosp_dod", "cfr_all_dod")))

start <- Sys.time()

# (1) Iterate through knot date pairs
for (i in 1:n_knots_best) {
  
  # Refilter best knots dataframe by row i
  knots_best_country_counterfactual_i <- knots_best_country_counterfactual %>% filter(row_number() == i)
  
  # Record number of knots
  n_knots <- knots_best_country_counterfactual_i %>% pull(N_knots)
  
  # Define number of simulation runs for specified knot dates
  n_runs_i <- knots_best_country_counterfactual_i %>% pull(N_unequal)
  
  # Set knot dates
  knot_date_1_i <- knots_best_country_counterfactual_i %>% pull(Knot_date_1)
  knot_date_2_i <- knots_best_country_counterfactual_i %>% pull(Knot_date_2) 
  
  # Define mean growth parameters
  growth_factor_1_i <- knots_best_country_counterfactual_i %>% pull(Growth_factor_1)
  growth_factor_2_i <- knots_best_country_counterfactual_i %>% pull(Growth_factor_2)
  growth_factor_3_i <- knots_best_country_counterfactual_i %>% pull(Growth_factor_3)
  
  # Define number of simulation runs for specified knot dates
  n_runs_i <- knots_best_country_counterfactual_i %>% pull(N_unequal)  
  
  # Create matrices for simulated data for given knot dates
  # (1 row per simulation run, 1 col per date)
  daily_cases_sim_i <- matrix(nrow = n_runs_i, ncol = length(dates) + 1, 
                              dimnames = list(NULL, as.character(c(min(dates) - 1, dates))))
  # Initialise matrices with data at date_pop_pct - 1
  daily_cases_sim_i[, 1] <- data_eur_country %>% 
    filter(Date == (date_pop_pct - 1)) %>% pull(Daily_cases_MA7)
  
  # (2) Iterate through dates
  for (t in as.list(dates)) {
    
    # Get daily cases from time t-1
    inc_tminus1 <- daily_cases_sim_i[, as.character(t-1)]

    # Define growth parameters
    if (n_knots == 0) {  # NO knot points
      growth <- growth_factor_1_i  
    } else if (n_knots == 1) {  # ONE knot point
      if (t <= knot_date_1_i) {
        growth <- growth_factor_1_i
      } else {
        growth <- growth_factor_2_i
      }
    } else {  # TWO knot points
      if (t <= knot_date_1_i) {
        growth <- growth_factor_1_i
      } else if (t <= knot_date_2_i) {
        growth <- growth_factor_2_i
      } else {
        growth <- growth_factor_3_i
      }
    }
    
    # Calculate daily cases at time t and record
    inc_t <- growth*inc_tminus1
    daily_cases_sim_i[, as.character(t)] <- inc_t
    
  }  # (close loop 2 - t)
  
  # Bind knot-specific dataframes to full scenario dataframe
  daily_cases_sim <- rbind(daily_cases_sim, daily_cases_sim_i)

}  # (close loop 1 - i)

# Compute cumulative cases 
a <- daily_cases_sim
a[, 1] <- data_eur_country %>% filter(Date == (date_pop_pct - 1)) %>% pull(Cumulative_cases_end_MA7)
cumulative_cases_end_sim <- apply(X = a, MARGIN = 1, FUN = cumsum) %>% t

# Calculate mean incident/cumulative cases
mean_daily_cases_sim <- apply(X = daily_cases_sim, MARGIN = 2, FUN = mean) %>% 
  as_tibble(rownames = "Date") %>% mutate(Date = as.Date(Date))
mean_cumulative_cases_end_sim <- apply(X = cumulative_cases_end_sim, MARGIN = 2, FUN = mean) %>%
  as_tibble(rownames = "Date") %>% mutate(Date = as.Date(Date))

end <- Sys.time()
end - start  # ~ 0.2 seconds

mean_cumulative_cases_end_sim %>% filter(Date == date_T) %>% pull(value)

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

# Create figures of incident and cumulative cases under natural vs counterfactual histories

# Incident cases
plot_inc <- ggplot(data = mean_daily_cases_sim, aes(x = Date, y = value)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_col(data = filter(data_eur_country, Date <= date_T), 
           aes(x = Date, y = Daily_cases), alpha = 0.4) +
  #geom_vline(xintercept = date_first_restriction - n_days_counterfactual_first_restriction, col = "red4") +
  #geom_vline(xintercept = date_lockdown - n_days_counterfactual_lockdown, col = "red4") +
  geom_line(color = "navyblue", size = 1) +
  scale_x_date(name = "Date", 
               limits = c(as.Date(NA), date_T + 7), 
               date_breaks = "1 week", 
               date_labels = "%d %b %C") +
  scale_y_continuous(name = "Daily number of lab-confirmed cases",
                     labels = comma_format(accuracy = 1))

# Cumulative cases
plot_cum <- ggplot(data = mean_cumulative_cases_end_sim, aes(x = Date, y = value)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #labs(title = simulation_i,
  #     subtitle = description_i) +
  geom_col(data = filter(data_eur_country, Date <= date_T), 
           aes(x = Date, y = Cumulative_cases_end), alpha = 0.4) +
  geom_line(color = "navyblue", size = 1) +
  geom_text(data = mean_cumulative_cases_end_sim, 
            aes(x = Date, y = value,
                label = ifelse(Date == date_T, formatC(value, 
                                                       format = "f", big.mark = ",", digits = 0), "")),
            vjust = -1, size = 3, color = "navyblue", inherit.aes = FALSE) +
  scale_x_date(name = "Date", 
               limits = c(as.Date(NA), date_T + 7), 
               date_breaks = "1 week", 
               date_labels = "%d %b %C") +
  scale_y_continuous(name = "Cumulative number of lab-confirmed cases",
                     labels = comma_format(accuracy = 1))
#plot_cum

# Plot incident and cumulative cases
text <- paste(country, 
              paste("First restriction", n_days_counterfactual_first_restriction, "days earlier", sep = " "),
              paste("Lockdown", n_days_counterfactual_lockdown, "days earlier", sep = " "),
              sep = "\n")
ggarrange(plotlist = list(plot_inc, plot_cum), nrow = 2) %>%
  annotate_figure(., top = text_grob(text))

