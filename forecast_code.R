install.packages('fpp3') # package for applying simple forecasting methods
install.packages('tsibble') # package for dealing with time series data sets and tsibble objects

library(tidyverse)
library(lubridate)

source('R/ignore_sigpipe.R')

#read in the targets data
targets <- readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/terrestrial_daily/terrestrial_daily-targets.csv.gz", guess_max = 1e6) |> 
  na.omit()

# read in the sites data
terrestrial_sites <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(terrestrial == 1)

dbl_sites <- terrestrial_sites %>%
  filter(phenocam_vegetation == 'deciduous broadleaf')

targets <- targets %>%
  filter(site_id %in% dbl_sites$field_site_id)

# Only use nee
targets <- targets %>%
  filter(variable == 'nee')


# get NOAA forecasts
# Create the connections to data products
df_past <- neon4cast::noaa_stage3()

forecast_date <- Sys.Date() 
noaa_date <- forecast_date - days(1)
df_future <- neon4cast::noaa_stage2()

# specify the covariates
variables <- c("air_temperature", "surface_downwelling_shortwave_flux_in_air")

lm_forecast <- NULL
model_fit <- NULL


## Stash NOAA locally. This will speed up the access to NOAA by writing a local copy that we query
df_past |> 
  dplyr::filter(datetime >= ymd('2017-01-01'),
                variable %in% variables,
                site_id  %in% dbl_sites$field_site_id) |> 
  arrow::write_dataset("noaa_past", partitioning="site_id")

df_future |> 
  dplyr::filter(reference_datetime == noaa_date,
                datetime >= forecast_date,
                variable %in% variables,
                site_id %in% dbl_sites$field_site_id,
                horizon < 744,
                parameter < 31) |> 
  arrow::write_dataset("noaa_future", partitioning="site_id")


for(i in 1:length(dbl_sites$field_site_id)) {  
  
  site <- dbl_sites$field_site_id[i]
  
  # 1. Download historic NOAA data
  noaa_past <- arrow::open_dataset("noaa_past") |> 
    dplyr::filter(site_id %in% site,
                  datetime >= ymd('2017-01-01'),
                  variable %in% variables) |> 
    dplyr::collect()
  
  # calculate a daily mean to fit the model
  noaa_past_daily <- noaa_past |> 
    mutate(datetime = as_date(datetime)) |> 
    group_by(datetime, site_id, variable) |> 
    summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
    pivot_wider(names_from = variable, values_from = prediction) |> 
    # convert air temp to C
    rename(shortwave = surface_downwelling_shortwave_flux_in_air) |> 
    mutate(air_temperature = air_temperature - 273.15)
  
  # message('Stage 3 for ', site)
  
  
  #2. Download future NOAA data
  # Download the stage2 data
  noaa_future <- arrow::open_dataset("noaa_future") |> 
    dplyr::filter(reference_datetime == noaa_date,
                  datetime >= forecast_date,
                  site_id %in% site,
                  variable %in% variables,
                  horizon < 744,
                  parameter < 31) |> 
    dplyr::collect()
  
  # generate a mean daily forecast to use in the forecast
  noaa_future_daily <- noaa_future |> 
    mutate(datetime = as_date(datetime)) |> 
    # mean daily forecasts at each site per ensemble
    group_by(datetime, site_id, parameter, variable) |> 
    summarize(prediction = mean(prediction)) |>
    pivot_wider(names_from = variable, values_from = prediction) |>
    # convert to Celsius
    mutate(air_temperature = air_temperature - 273.15) |> 
    rename(shortwave = surface_downwelling_shortwave_flux_in_air) |> 
    select(datetime, site_id, air_temperature, shortwave, parameter)
  
  # message('Stage 2 for ', site)
  
  #3. Fit the model
  # targets data reformatted to aid model fitting
  
  targets_lm <- targets |> 
    filter(site_id == site) |> 
    pivot_wider(names_from = 'variable', values_from = 'observation') |> 
    left_join(noaa_past_daily, 
              by = c("datetime","site_id"))
  
  #Fit linear model based on past data
  fit <- lm(targets_lm$nee ~ targets_lm$air_temperature + targets_lm$shortwave)
  
  # use linear regression to forecast NEE for each ensemble member
  forecasted_nee <- fit$coefficients[1] + 
    (fit$coefficients[2] * noaa_future_daily$air_temperature) + 
    (fit$coefficients[3] * noaa_future_daily$shortwave)
  
  # put all the relevant information into a tibble that we can bind together
  NEE <- tibble(datetime = noaa_future_daily$datetime,
                site_id = site,
                parameter = noaa_future_daily$parameter,
                prediction = forecasted_nee,
                variable = "nee")
  
  lm_forecast <- dplyr::bind_rows(lm_forecast, NEE)
  message(site, ' NEE forecast run')
  
  # extract the model fit
  # you can comment/uncomment this out to extract the R-squared from the model summary
  
  model_fit <- dplyr::bind_rows(model_fit, data.frame(site_id = site,
                                                      r_squared = summary(fit)$r.squared))
  
  
}

# Remember to change the model_id when you make changes to the model structure!
my_model_id <- 'test_model'

lm_forecast_EFI <- lm_forecast %>%
  mutate(model_id = my_model_id,
         reference_datetime = as_date(min(datetime)) - days(1),
         family = 'ensemble',
         parameter = as.character(parameter)) %>%
  select(model_id, datetime, reference_datetime, site_id, family, parameter, variable, prediction)


# Start by writing the forecast to file
theme <- 'terrestrial_daily'
date <- lm_forecast_EFI$reference_datetime[1]
forecast_name <- paste0(lm_forecast_EFI$model_id[1], ".csv")

# Write the file locally
forecast_file <- paste(theme, date, forecast_name, sep = '-')
forecast_file

if (!dir.exists('Forecasts')) {
  dir.create('Forecasts')
}


write_csv(lm_forecast_EFI, file.path('Forecasts', forecast_file))


# neon4cast::submit(forecast_file = file.path('Forecasts', forecast_file),
#                   ask = FALSE) # if ask = T (default), it will produce a pop-up box asking if you want to submit

