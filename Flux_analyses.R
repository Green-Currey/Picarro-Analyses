# Picaro flux analyses


# Packages
library(readr)
library(dplyr)
library(lubridate)
library(stringr)


# Directory where PF monthly files are stored
dir <- 'c:/users/bryce/Dropbox/Bryce and Justin Folder/Picarro/PF_Flux_Data_Picarro/'
setwd(dir)

# Data where flux data is accessed
df <- read_csv('c:/users/bryce/Dropbox/Bryce and Justin Folder/Picarro/PF_Flux_2019.csv')

# vector of months for finding and extracting the picarro data
months <- c('January',
            'February',
            'March',
            'April',
            'May',
            'June',
            'July',
            'August',
            'September',
            'October',
            'November',
            'December')

# Dimensions of chamber ---------------------------------------------------

d <- 0.2032 # diameter (m)
r <- d/2 # radius (m)
h <- 0.15 # height (m)

A <- pi*r^2 # area (m^2)
V <- A*h * 1000 # volume (liters)
R <- 8.314 # gas constant (J/mol K)




# Flux calculation --------------------------------------------------------

for (i in seq(nrow(df))) {
  
  # date of measurement (e.g. May 3rd, 2019)
  date <- ymd(paste(df$Year[i],df$Month[i],df$Sample_Day[i]))
  
  # date in formate for extracting the data (eg 'May2019/03')
  date.str <- paste0(months[month(date)],
                     year(date), 
                     '/', 
                     str_pad(df$Sample_Day[i], width = 2, side = 'left', '0')
  )
  
  # time on is when the button on the picarro was clicked to start recording of flux data
  # time off is 10 minues after the button was clicked
  time.on <- df$ON[i]
  time.off <- time.on+10
  
  # extracts the data file that corresponds with 'time.on' from the pathway created from 'dir' variable and 'date.str' variable
  flux.file <- list.files(path = paste0(dir, '/', date.str), pattern = paste(time.on), full.names = T)
  
  # checks if flux.file doesnt exist (i.e. there was no file corresponding to 'time.on' variable)
  # also checks if flux.file is more than one file, indicating that the button was clicked twice in the same
  # minute. This script chooses the second file created, corresponding to the second time the button was clicked.
  if (flux.file %>% length == 0) { next }
  if (flux.file %>% length > 1) { flux.file <- flux.file[2] }
    
  # formating the data.
  # 1) reads the data from read_table
  # 2) from the DATE and TIME columns from the picarro data file, creates one variable in mountain time (e.g. '2019-05-03 12:00:00 MDT')
  # 3) puts the time variable into a numeric format (e.g. 12:34:51 turns into 123451)
  # 4) chooses the variables of interest
  # 5) filters for times that occur after 'time.on' and before 10 minutes following 'time.on'
  data <- read_table(flux.file) %>%
    mutate(time = ymd_hms(paste(DATE, TIME), tz = 'America/Denver')) %>%
    mutate(time = strftime(time, format = "%H%M%S") %>% as.numeric) %>%
    mutate(time = time) %>%
    select(time, N2O_dry, CO2, CH4, NH3) %>%
    filter(time >= paste0(time.on,'00') %>% as.numeric) %>%
    filter(time <= paste0(time.off,'00') %>% as.numeric)
  
  
  # obtaining the slope coefficints from linear models between time and flux
  N2O.slope <- coefficients(lm(N2O_dry ~ time, data))[2] %>% as.numeric
  CO2.slope <- coefficients(lm(CO2 ~ time, data))[2] %>% as.numeric
  CH4.slope <- coefficients(lm(CH4 ~ time, data))[2] %>% as.numeric
  NH3.slope <- coefficients(lm(NH3 ~ time, data))[2] %>% as.numeric
  
  
  # flux calculations
  N2O.flux = N2O.slope * (V / (A * R * (273.15 + df$Temp[i]))) * 3600
  CO2.flux = CO2.slope * (V / (A * R * (273.15 + df$Temp[i]))) * 3600 * 0.001
  CH4.flux = CH4.slope * (V / (A * R * (273.15 + df$Temp[i]))) * 3600
  NH3.flux = NH3.slope * (V / (A * R * (273.15 + df$Temp[i]))) * 3600
  
  # appending the data frame
  df$Flux_N2O[i] <- N2O.flux
  df$Mean_N2O[i] <- mean(data$N2O_dry, na.rm = T)
  df$SD_N2O[i] <- sd(data$N2O_dry, na.rm = T)
  
  df$Flux_CO2[i] <- CO2.flux
  df$Mean_CO2[i] <- mean(data$CO2, na.rm = T)
  df$SD_CO2[i] <- sd(data$CO2, na.rm = T)
  
  df$Flux_CH4[i] <- CH4.flux
  df$Mean_CH4[i] <- mean(data$CH4, na.rm = T)
  df$SD_CH4[i] <- sd(data$CH4, na.rm = T)
  
  df$Flux_NH3[i] <- NH3.flux
  df$Mean_NH3[i] <- mean(data$NH3, na.rm = T)
  df$SD_NH3[i] <- sd(data$NH3, na.rm = T)
}

# rewrite data
write_csv(df, 'c:/users/bryce/Dropbox/Bryce and Justin Folder/Picarro/PF_Flux_2019.csv')
