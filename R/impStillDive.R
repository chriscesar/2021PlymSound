# impStillDive.R #### 
# Import Still & Dive data, save in Long format

# Set up ####
source("R/datfol.R")
ld_pkgs <- c("tidyverse") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

# set data folder:
df_dive <- openxlsx::read.xlsx(paste0(datfol,"Data/Dive surveys/Dive Survey Data 2017-2021.xlsx"),
                               sheet = "Standardised 2017-2021")

df_still <- openxlsx::read.xlsx(paste0(datfol,"Data/2021 data/DDV/2021 DDV Stills .xlsx"),
                               sheet = "Standardised ALL")

df_diveL <- df_dive %>% 
  pivot_longer(cols = "Acrosorium.ciliolatum":"Tunicata",
               names_to = "Taxon", values_to = "Abund")

df_stillL <- df_still %>% 
  pivot_longer(cols = "Alcyonidium.diaphanum":"Terebellidae",
               names_to = "Taxon", values_to = "Abund")

df <- rbind(df_diveL,df_stillL)

write.csv(df, file = paste0(datfol,"Data/2021 data/","StillDiveLong.csv"),row.names = FALSE)
write.csv(unique(df$Taxon),file = paste0(datfol,"Data/2021 data/","StillDivetaxa.csv"),row.names = FALSE)

## tidy up
rm(list = ls(pattern = "^df"))
rm(datfol)

detach("package:tidyverse", unload=TRUE)