# an.inf_v2.R #### 
# Import and analyse new version of infauna data
# Only comparing sites within Plymouth Sound

# Set up ####
source("R/metadata.R")

## load required packages ####
ld_pkgs <- c("tidyverse","vegan","lmerTest","rstatix", "mvabund",
             "MASS","ggtext","ggpmisc", "gllvm") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

# Copy data:
file_name <- "Inf_Sed_all_years_Long_updated.xlsx"

# Set the source and destination file paths
source_file <- paste0(datfol,"Data/2021 data/Infauna/Updated/",file_name)
destination_file <- paste0("data/in/",file_name)

# Check if the file exists in the data folder
if (!file.exists(destination_file)) {
  # If not, copy the file from the source folder to the data folder
  file.copy(source_file, destination_file)
  cat("File copied successfully.\n")
} else {
  # If the file already exists in the data folder, do nothing
  cat("File already exists in the data folder.\n")
}

# load data
df0 <- as_tibble(openxlsx::read.xlsx("data/in/Inf_Sed_all_years_Long_updated.xlsx",
                                     sheet = "Matched sites"))

df0 %>% rename("MatchedSite"=`2021.matched.site`) -> df0

# 1. filter out taxa to be 'discarded' and tidy up MatchedSite variable
df0 %>% 
  # change taxa flagged as -999 to values of 1. Consider removal of these taxa if required
  # mutate(Abund = ifelse(Abund == -999, 1, Abund)) %>% 
  # filter(.,Abund != 0) %>% #remove 'zero' abundances
  mutate(MatchedSite = paste0("PL_", substr(MatchedSite, 6, 7))) %>% #rename MatchedSite contents
  mutate(MatchedSiteYr = paste0(MatchedSite, "_", Year)) %>% #append Year value
  relocate(MatchedSiteYr, .after = MatchedSite) %>% #move new variable somewhere sensible
  #remove taxa flagged for removal
  filter(.,!str_detect(Flag, "^Remove")) -> dfl0

# next step:
# 2. generate mean values for taxa in replicate samples
# (i.e., those with non-numeric values at end of dfl$Site variable)

## extract Replicates: rows with numeric value at end of Site name
dfl_rep <- dfl0 %>% 
  filter(str_sub(Site,-1) %in% letters) %>% # extract rows with numeric value at end of Site name
  mutate(Site = substr(Site, 1, nchar(Site) - 1)) %>% #remove final character from Site
  dplyr::select(.,-Taxon) %>% 
  group_by(across(-Abund)) %>% 
  summarise(Abund=mean(Abund),.groups = "drop")

### set aside non-Replicate rows:
dfl_nonrep <- dfl0 %>% 
  filter(!str_sub(Site,-1) %in% letters) %>% 
  dplyr::select(.,-Taxon)

dfl <- rbind(dfl_nonrep,dfl_rep)

# 3. 