# Picaro flux analyses by single date + plots

###General notes: Make sure that raw data file in dropbox has two digits. i.e the 3rd of May should be 03
##make sure to clear global environment before every run

rm(list = ls())

# Packages
library(readr)
library(dplyr)
library(lubridate)
library(stringr)
library(ggpubr)

# Directory where PF monthly files are stored
dir <- 'c:/users/bryce/onedrive/Documents/Other Projects/Justin Picarro/'
setwd(dir)


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
R <- .08314 # gas constant (LK mol-1) 


paste(test)
# Plotting by sampling day and time ------------------------------------------------

# Enter Year
Year <- 2020
df.main <- read_csv(paste0('PF_Flux_', Year, '.csv'))
# searches for all days
print(unique(df.main$Month))


# Enter Month numerically
Month <- 09
m = months[Month]

df <- df.main %>% filter(Month == m)# Filter by day
Month <- str_pad(Month, width = 2, side = 'left', pad = '0')

# searches for all days
print(unique(df$Sample_Day))


# Enter day
Day <- 14
# Filter by day
df <- df %>% filter(Sample_Day == Day)
Day <- str_pad(Day, width = 2, side = 'left', pad = '0')


paste('Running script for:',m,Day,Year)

for (i in seq(nrow(df))) {
  
  time.on <- df$ON[i]
  minutes <- time.on %>% as.character() %>% substr(3,4) %>% as.numeric()
  
  # See main script (e.g., 2020_PF_Flux_v2.R) for more information on these if statements.
  if (minutes >= 55 & minutes < 59) {
    
    new.time.on <- time.on + 1
    new.time.off <- (time.on + 45) + 1
    
  }
  
  if (minutes == 59) {
    new.time.on <- time.on + 41
    new.time.off <- new.time.on + 5
    
  }
  
  if (minutes < 55) {
    new.time.on <- time.on + 1
    new.time.off <- new.time.on + 5
    
  }
  
  Tx <- df$Treatment[i]
  Plot <- df$Plot[i]
  # Finds specified file
  
  flux.file <- list.files(path = paste0(m,Year,'/',Day,'/'), pattern = paste0(Year, Month, Day, '-', time.on), full.names = T)
  
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
    dplyr::select(time, N2O_dry, CO2, CH4, NH3) %>%
    filter(time >= paste0(new.time.on,'00') %>% as.numeric) %>%
    filter(time <= paste0(new.time.off,'00') %>% as.numeric) %>%
    mutate(time2 = seq(n()),
           Treatment = rep(Tx),
           Plot = rep(Plot),
           Plot_Tx = paste0(rep(Plot),'. ', Tx),
           On = rep(new.time.on))
  
  if (i == 1) {to.plot <- data} else {to.plot <- to.plot %>% rbind(data)}
}

cols <- c('grey47', 'darkolivegreen3', 'darkorange')

N2O <- ggplot(data = to.plot) +
  aes(x = time2, y = N2O_dry, group = Treatment, color = Treatment) +
  geom_point() + geom_smooth() +
  facet_grid(.~Plot_Tx, scales = 'free') +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none') +
  labs(x = '')



CO2 <- ggplot(data = to.plot) +
  aes(x = time2, y = CO2, group = Treatment, color = Treatment) +
  geom_point() + geom_smooth() +
  facet_grid(.~Plot_Tx, scales = 'free') +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none') +
  labs(x = '')



CH3 <- ggplot(data = to.plot) +
  aes(x = time2, y = CH4, group = Treatment, color = Treatment) +
  geom_point() + geom_smooth() +
  facet_grid(.~Plot_Tx, scales = 'free') +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none') +
  labs(x = '')



NH3 <- ggplot(data = to.plot) +
  aes(x = time2, y = NH3, group = Treatment, color = Treatment) +
  geom_point() + geom_smooth() +
  facet_grid(.~Plot_Tx, scales = 'free') +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none') +
  labs(x = '5 minute Interval')



tiff(filename = paste0('Raw flux measurements for ',m,Day,'_',Year, '.tif'), height = 16, width = 16, res = 150, units = 'in')
ggarrange(N2O, CO2, CH3, NH3,
          nrow = 4, 
          align = 'hv') %>%
  annotate_figure(top = text_grob(paste0(m,' ',Day,', ',Year), size = 25))
dev.off()
