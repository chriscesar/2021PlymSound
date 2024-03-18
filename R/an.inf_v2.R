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
                 