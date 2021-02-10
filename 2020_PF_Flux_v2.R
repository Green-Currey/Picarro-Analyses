# Picaro flux analyses

rm(list = ls())

# Packages
library(readr)
library(dplyr)
library(lubridate)
library(stringr)


# Directory where PF monthly files are stored
dir <- '~/Other Projects/Justin Picarro/'
setwd(dir)

# Data where flux data is accessed
df <- read.csv('PF_Flux_2020.csv')

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

alpha <- 0.05

# Dimensions of chamber ---------------------------------------------------

d <- 0.2032 # diameter (m)
r <- d/2 # radius (m)
h <- 0.15 # height (m)

A <- pi*r^2 # area (m^2)
V <- A*h * 1000 # volume (liters)
R <- .08314 # gas constant (LK mol -1) #updated from Christiansen




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
  #     this is important for locating the flux.file in the the directory
  time.on <- df$ON[i] 
  minutes <- time.on %>% as.character() %>% substr(3,4) %>% as.numeric()

  # For filtering the data, we use one minute past time on (time.on + 1)
  #     and 5 minutes after that time.
  #     This is complicated when the time goes over the hour (e.g. 10:57-11:02)
  #     because 1057 + 5 = 1062 or 10:62, which is a fake time.
  
  # So, if time.on is 5 minutes from the next hour (e.g. 10:55-10:59), 
  #     time.off (time.on+5) would be a nonexistent time (e.g., 10:60-10:64). 
  #     To fix this, we add 100 to the nonexistent time.off to move up to the next hour (e.g., 10:62 -> 11:62)
  #     then subtract off the excess 60 minutes to go back to the hour (e.g. 11:62 -> 11:02).
  #     This is equivalent to adding 45 minutes to time.on (time.on + 5) + 100 - 60
  # 
  #     Thus, time.on = 10:57 now has a time.off that is 5 minutes later (11:02)
  
  if (minutes >= 55 & minutes < 59) {
    
    new.time.on <- time.on + 1
    new.time.off <- (time.on + 45) + 1
    
  }
  
  # if time.on is :59 (e.g. 10:59), time.on+1 would make it a fake time (e.g., 10:60)
  #     by adding 100, we move up to the next hour (e.g., 10:59 -> 11:59)
  #     then subtract off the excess 59 minutes (11:59 -> 11:00).
  #     This is the equivalent of adding 41 to time on (time.on + 100 - 59) 
  # if time.on isn't :59, simply adds one (above).
  if (minutes == 59) {
    new.time.on <- time.on + 41
    new.time.off <- new.time.on + 5
    
  }
  
  # Else this
  
  if (minutes < 55) {
    new.time.on <- time.on + 1
    new.time.off <- new.time.on + 5
    
  }
  
  # extracts the data file that corresponds with 'time.on' from the pathway created from 'dir' variable and 'date.str' variable
  flux.file <- list.files(path = paste0(dir,date.str), pattern = paste(time.on), full.names = T)
  
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
  # 5) filters the data between the time.on+1 and time.off+1 points.
  
  data <- read_table(flux.file) %>%
    mutate(time = ymd_hms(paste(DATE, TIME), tz = 'America/Denver')) %>%
    mutate(time = strftime(time, format = "%H%M%S") %>% as.numeric) %>%
    select(time, N2O_dry, CO2, CH4_dry, NH3)%>% 
    filter(time >= paste0(new.time.on,'00') %>% as.numeric) %>%
    filter(time <= paste0(new.time.off,'00') %>% as.numeric) %>%
    mutate(time2 = seq(n()))
  
  
  # Significance filter
  N2O_p <- if_else(summary(lm(N2O_dry ~ time2, data = data))$coefficients[2,4] < alpha, 1, 0)
  CO2_p <- if_else(summary(lm(CO2 ~ time2, data = data))$coefficients[2,4] < alpha, 1, 0)
  CH4_p <- if_else(summary(lm(CH4_dry ~ time2, data = data))$coefficients[2,4] < alpha, 1, 0)
  NH3_p <- if_else(summary(lm(NH3 ~ time2, data = data))$coefficients[2,4] < alpha, 1, 0)
  
  
  # obtaining the slope coefficints from linear models between time and flux
  N2O.slope <- coefficients(lm(N2O_dry ~ time2, data))[2] %>% as.numeric
  CO2.slope <- coefficients(lm(CO2 ~ time2, data))[2] %>% as.numeric
  CH4.slope <- coefficients(lm(CH4_dry ~ time2, data))[2] %>% as.numeric
  NH3.slope <- coefficients(lm(NH3 ~ time2, data))[2] %>% as.numeric
  
  
  # flux calculations
  N2O.flux = N2O.slope * (V / (A * R * (273.15 + df$Temp[i]))) * 3600 
  CO2.flux = CO2.slope * (V / (A * R * (273.15 + df$Temp[i]))) # multiply my 0.001 to get mmol per hour
  CH4.flux = CH4.slope * (V / (A * R * (273.15 + df$Temp[i]))) * 3600
  NH3.flux = NH3.slope * (V / (A * R * (273.15 + df$Temp[i]))) * 3600
  #####* 24 is to get to mmol per m^2 per day (not hour)
  
  # appending the data frame. If P-value > alpha, Flux = NA
  df$Flux_N2O[i] <- ifelse(N2O_p==1, N2O.flux, NA)
  df$Mean_N2O[i] <- ifelse(N2O_p==1, mean(data$N2O_dry, na.rm = T), NA)
  df$SD_N2O[i] <- ifelse(N2O_p==1, sd(data$N2O_dry, na.rm = T), NA)
  
  df$Flux_CO2[i] <- ifelse(CO2_p==1, CO2.flux, NA)
  df$Mean_CO2[i] <- ifelse(CO2_p==1, mean(data$CO2, na.rm = T), NA)
  df$SD_CO2[i] <- ifelse(CO2_p==1, sd(data$CO2, na.rm = T), NA)
  
  df$Flux_CH4[i] <- ifelse(CO2_p==1, CH4.flux, NA)
  df$Mean_CH4[i] <- ifelse(CO2_p==1, mean(data$CH4_dry, na.rm = T), NA)
  df$SD_CH4[i] <- ifelse(CO2_p==1, sd(data$CH4_dry, na.rm = T), NA)
  
  df$Flux_NH3[i] <- ifelse(CO2_p==1, NH3.flux, NA)
  df$Mean_NH3[i] <- ifelse(CO2_p==1, mean(data$NH3, na.rm = T), NA)
  df$SD_NH3[i] <- ifelse(CO2_p==1, sd(data$NH3, na.rm = T), NA)
}

# rewrite data
write.csv(df, '~/Dropbox/Bryce_Justin_Folder/Picarro/PF_Flux_Data_Picarro/PF_Flux_2020.csv')
