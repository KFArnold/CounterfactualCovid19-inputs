# ------------------------------------------------------------------------------
# Notes
# ------------------------------------------------------------------------------

# This script downloads most recent data from:

# (1) COVID-19 Data Repository by the Center for Systems Science and Engineering (CSSE) at Johns Hopkins University
# (2) Oxford Covid-19 Government Response Tracker
# (3) World Bank (from wbstats package)

# Data is then saved to project repository.

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# Restore package library to last snapshot
packrat::restore()

# Load required packages
library(tidyverse); library(wbstats); library(sjlabelled)

# Define project directory where unformatted data is located
data_directory_u <- paste0("./Data/Unformatted/")

# ------------------------------------------------------------------------------
# Import data
# ------------------------------------------------------------------------------

## Cases and deaths data -------------------------------------------------------

# Set storage directory for saved data
out <- paste0(data_directory_u, "CSSE data/")

# Define web repository where CSSE data is located
source <- paste0("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/")

# Define data files to download
filenames <- list("time_series_covid19_confirmed_global.csv",
                  "time_series_covid19_deaths_global.csv")  # also data on recovered cases, if necessary

# Download data files and save to repository Data folder
for (i in filenames) {
  filenames_i <- i  # define filename
  data <- read_csv(url(paste0(source, filenames_i)))  # download file
  names(data) <- str_replace_all(names(data), c("/" = "_", " " = "_"))  # replace slashes and spaces with underscores in var names
  write_csv(x = data,
            file = paste0(out, filenames_i))  # write to repo folder
}

## Policy data -----------------------------------------------------------------

# Set storage directory for saved data
out <- paste0(data_directory_u, "OxCGRT data/")

# Define web repository where OxCGRT data is located
source <- paste0("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/")

# Define data files to download
filenames <- list("OxCGRT_latest.csv")  

# Download data files and save to repository Data folder
for (i in filenames) {
  filenames_i <- i  # define filenames
  data <- read_csv(url(paste0(source, filenames_i)),
                   col_types = cols(RegionName = col_character(),
                                    RegionCode = col_character()))  # download file
  names(data) <- str_replace_all(names(data), c("/" = "_", " " = "_"))  # replace slashes and spaces with underscores in var names
  write_csv(x = data,
            file = paste0(out, filenames_i))  # write to repo folder
}

## World Bank data -------------------------------------------------------------

# Set storage directory for saved data
out <- paste0(data_directory_u, "World Bank data/")

# Download data - population size
data <- wb_data(indicator = "SP.POP.TOTL", mrv = 1) %>%
  rename(Population = SP.POP.TOTL, Year = date) %>% remove_all_labels()

# Capitalise first letter of every variable name
names(data) <- str_to_title(names(data))

# Select relevant variables
data <- data %>% select(contains("iso"), Country, Year, Population)

# Save to repository Data folder
write_csv(x = data, file = paste0(out, "Worldbank_data.csv"))

# Remove variables from environment
rm(data_directory_u, out, source, filenames, i, filenames_i, data)
