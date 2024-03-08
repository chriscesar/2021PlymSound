# inf.data.import.01.R #### 
# Title: Import infaunal data and save for analysis

# Set up ####
source("R/datfol.R")

## Load required packages ####
ld_pkgs <- c("tidyverse") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

# Import data ####
#dir(paste0(datfol,"Data"))
df0_2011  <- read.csv(paste0(datfol,"Data/2021 data/Infauna/2011_inf.csv"))
df0_2015  <- read.csv(paste0(datfol,"Data/2021 data/Infauna/2015_inf.csv"))
df0_2021  <- read.csv(paste0(datfol,"Data/2021 data/Infauna/2021_inf.csv"))

dfl_2011 <- df0_2011 %>% 
  pivot_longer(cols = "Astrorhiza":"Branchiostoma.lanceolatum",
             names_to = "Taxon", values_to = "Abund")

dfl_2015 <- df0_2015 %>% 
  pivot_longer(cols = "Suberitidae":"Insecta.Larva",
               names_to = "Taxon", values_to = "Abund")

dfl_2021 <- df0_2021 %>% 
  pivot_longer(cols = "ANIMALIA":"SEEDS.TOMATO",
               names_to = "Taxon", values_to = "Abund")

dfl_2011$Year <- 2011
dfl_2015$Year <- 2015
dfl_2021$Year <- 2021

dfl <- dplyr::bind_rows(dfl_2011,dfl_2015,dfl_2021)
rm(df0_2011,df0_2015,df0_2021)
rm(dfl_2011,dfl_2015,dfl_2021)

tail(dfl)

dfw <- dfl %>% 
  group_by(Site, Mesh, Easting, Northing, Taxon, Year) %>% 
  summarise(Abund = sum(Abund),.groups = "drop") %>% 
  pivot_wider(names_from = Taxon, values_from = Abund,values_fill=list(Abund = 0))

write.csv(dfw, file = paste0(datfol,"Data/2021 data/Infauna/AllInfWide.csv"), row.names = FALSE)
write.csv(dfl, file = paste0(datfol,"Data/2021 data/Infauna/AllInfLong.csv"), row.names = FALSE)

## tidy up
rm(list=ls(pattern="^df"))
rm(datfol)

detach("package:tidyverse", unload=TRUE)
