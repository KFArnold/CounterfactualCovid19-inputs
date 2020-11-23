# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# Restore package library to last snapshot
packrat::restore()

# Load required packages
library(tidyverse); library(caTools)

# Run source code to update external data
#source("./Code/Update data.R")

# Define all European countries
# (from https://www.worldometers.info/geography/how-many-countries-in-europe/)
countries <- list("Albania", "Andorra", "Austria", "Belarus", "Belgium", 
                  "Bosnia and Herzegovina", "Bulgaria", "Croatia", "Czech Republic",
                  "Denmark", "Estonia", "Finland", "France", "Germany", 
                  "Greece", "Holy See", "Hungary", "Iceland", "Ireland", 
                  "Italy", "Latvia", "Liechtenstein", "Lithuania", 
                  "Luxembourg", "Malta", "Moldova", "Monaco", "Montenegro",
                  "Netherlands", "North Macedonia", "Norway", "Poland",
                  "Portugal", "Romania", "Russia", "San Marino",
                  "Serbia", "Slovak Republic", "Slovenia", "Spain", "Sweden", 
                  "Switzerland", "Ukraine", "United Kingdom")

# ------------------------------------------------------------------------------
# Import and format data
# ------------------------------------------------------------------------------

# Load data
## (1) CSSE data (cumulative cases and deaths)
cases <- read_csv("./Data/CSSE data/time_series_covid19_confirmed_global.csv")
deaths <- read_csv("./Data/CSSE data/time_series_covid19_deaths_global.csv")
## (2) OxCGRT data (government policies)
policies <- read_csv("./Data/OxCGRT data/OxCGRT_latest.csv",
                     col_types = cols(RegionName = col_character(),
                                      RegionCode = col_character())) %>% 
  mutate(Date = as.Date(as.character(Date), format = "%Y%m%d"))
## World Bank data (population size)
worldbank <- read_csv("./Data/World Bank data/Worldbank_data.csv")

# Convert datasets to long form, select relevant variables, 
# rename country variable, convert characters to factors,
# and order / group by Province and Country
cases <- cases %>% gather(Date, Cumulative_cases_end, -(1:4)) %>% 
  mutate(Date = as.Date(Date, format = "%m_%d_%y")) %>% 
  select(-c(Lat, Long)) %>% rename(c(Country = Country_Region)) %>% 
  mutate_if(is.character, as.factor) %>%
  group_by(Province_State, Country) %>% arrange(Country, Date)
deaths <- deaths %>% gather(Date, Cumulative_deaths_end, -(1:4)) %>% 
  mutate(Date = as.Date(Date, format = "%m_%d_%y")) %>% 
  select(-c(Lat, Long)) %>% rename(c(Country = Country_Region)) %>% 
  mutate_if(is.character, as.factor) %>%
  group_by(Province_State, Country) %>% arrange(Country, Date)
policies <- policies %>% rename(c(Country = CountryName,
                                  Province_State = RegionName)) %>% 
  select(-contains("Code"), c(Date, C1_School_closing:C8_International_travel_controls)) %>%
  mutate_if(is.character, as.factor) %>%
  group_by(Province_State, Country) %>% arrange(Country, Date)
worldbank <- worldbank %>% mutate_if(is.character, as.factor) %>% 
  group_by(Country) %>% arrange(Country)

# Calculate daily cases and deaths
cases <- cases %>% mutate(Cumulative_cases_beg = lag(Cumulative_cases_end, n = 1, default = 0),
                          Daily_cases = Cumulative_cases_end - Cumulative_cases_beg) %>%
  relocate(c(Cumulative_cases_beg, Daily_cases), .before = Cumulative_cases_end)
deaths <- deaths %>% mutate(Cumulative_deaths_beg = lag(Cumulative_deaths_end, n = 1, default = 0), 
                            Daily_deaths = Cumulative_deaths_end - Cumulative_deaths_beg) %>%
  relocate(c(Cumulative_deaths_beg, Daily_deaths), .before = Cumulative_deaths_end)

# Calculate 7-day moving averages of daily and cumulative cases and deaths
cases <- cases %>% mutate(Daily_cases_MA7 = round(runmean(Daily_cases, k = 7, alg = "C", endrule = "mean"), 3),
                          Cumulative_cases_end_MA7 = cumsum(Daily_cases_MA7),
                          Cumulative_cases_beg_MA7 = lag(Cumulative_cases_end_MA7, n = 1, default = 0)) %>%
  relocate(Cumulative_cases_beg_MA7, .before = Daily_cases_MA7)
deaths <- deaths %>% mutate(Daily_deaths_MA7 = round(runmean(Daily_deaths, k = 7, alg = "C", endrule = "mean"), 3),
                            Cumulative_deaths_end_MA7 = cumsum(Daily_deaths_MA7),
                            Cumulative_deaths_beg_MA7 = lag(Cumulative_deaths_end_MA7, n = 1, default = 0)) %>%
  relocate(Cumulative_deaths_beg_MA7, .before = Daily_deaths_MA7)

# Merge cases and deaths datasets into single dataframe
data_all <- full_join(cases, deaths, by = c("Province_State", "Country", "Date")) 
rm(cases, deaths)  # (remove separate datasets)

# Create variables for: date of first case (Date_0), 
# date for which data can be reasonably assumed complete (Date_max)
data_all <- data_all %>% mutate(Date_0 = Date[which(Daily_cases >= 1)[1]],
                                Date_max = max(Date) - 7)

# Remove data after Date_max, since this is likely incomplete
data_all <- data_all %>% filter(Date <= Date_max)

# Rename Czechia and Slovakia in cases/deaths dataset,
# and Russian Federation in World Bank data
data_all <- data_all %>% mutate(Country = recode(Country, "Czechia" = "Czech Republic"),
                                Country = recode(Country, "Slovakia" = "Slovak Republic"))
worldbank <- worldbank %>% mutate(Country = recode(Country, "Russian Federation" = "Russia"))

# Retain data for countries in Europe only,
# remove dependencies (of the Netherlands, UK, France, and Denmark),
# and drop unused levels
data_eur <- data_all %>% filter(Country %in% countries) %>% ungroup(Province_State) %>%
  filter(is.na(Province_State)) %>% select(-Province_State) %>% droplevels
policies_eur <- policies %>% filter(Country %in% countries) %>% ungroup(Province_State) %>%
  filter(is.na(Province_State)) %>% select(-Province_State) %>% droplevels

# Create list of European countries for which we have both cases/deaths data and policy data
countries_eur <- as.list(intersect(levels(data_eur$Country), levels(policies_eur$Country)))
#### print note about any countries which both are not avaiable:
##if (length(countries_eur) != length(countries)) {
##  unavail <- setdiff(levels(data_eur$Country), levels(policies_eur$Country))
##  cat(paste0("Note that the following countries do not have both cases/deaths
##  data and policy data available:\n", 
##             paste0(unavail, collapse = ", ")))
##  rm(unavail)
##} 

# Retain data for European countries for which we have both cases/deaths data and policy data
data_eur <- data_eur %>% filter(Country %in% countries_eur) %>% droplevels
policies_eur <- policies_eur %>% filter(Country %in% countries_eur) %>% droplevels
worldbank_eur <- worldbank %>% filter(Country %in% countries_eur) %>% droplevels

# Calculate 0.0001% of population for each country
pct <- worldbank_eur %>% filter(Year == 2019) %>% 
  mutate(Pop_pct = 0.000001 * Population) %>% select(Country, Pop_pct)
# Determine date for which cumulative cases first exceeded this percent (Date_pop_pct)
# and add to data_eur dataframe
data_eur <- full_join(data_eur, pct, by = "Country") %>%
  mutate(Date_pop_pct = Date[which(Cumulative_cases_beg >= Pop_pct)[1]]) %>%
  select(-Pop_pct) %>% relocate(Date_pop_pct, .before = Date_max)
rm(pct)

# Remove non-European dataframes and lists
rm(countries, data_all, policies, worldbank)

# Ungroup dataframes
data_eur <- data_eur %>% ungroup
policies_eur <- policies_eur %>% ungroup
worldbank_eur <- worldbank_eur %>% ungroup
