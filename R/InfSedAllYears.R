# InfSedAllYears.R #### 
# Generate infaunal sediment data for analysis

# Set up ####
source("R/datfol.R")

## load required packages ####
ld_pkgs <- c("tidyverse") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

# load data (long format) ####
df0 <- openxlsx::read.xlsx(paste0(datfol,"Data/2021 data/Infauna/Inf_Sed_all_years.xlsx"),sheet = "out")

dfw <- as_tibble(df0) %>% 
  dplyr::select(., -c(Site_Yr:Mesh, Taxon)) %>% 
  filter(!grepl("^Rem",Flag)) %>% 
  dplyr::select(.,-Flag)

dfw <- dfw %>% 
  group_by(Site_Yr_mesh,Easting, Northing, Year, Taxon.USE, Station.Code,
           GEAR,Latitude,Longitude,OSGB36.Easting,OSGB36.Northing,Sample.use,GRAVELPCT,
           SANDPCT,MUDPCT,FOLK,BSH,BSH_CODE) %>% 
  summarise(Abund = sum(Abund)) %>% #sum values that share a taxon name AND a site site.year.mesh code
  ungroup()

dfw <- dfw %>% 
  pivot_wider(names_from = Taxon.USE, values_from = Abund,values_fill = list(Abund = 0))

write.csv(x=dfw,file=paste0(datfol,"Data/2021 data/Infauna/Inf_Sed_all_yearsWIDE.csv"))
saveRDS(dfw, file =paste0(datfol,"Data/2021 data/Infauna/Inf_Sed_all_yearsWIDE.Rdat"))
dfw <- readRDS(paste0(datfol,"Data/2021 data/Infauna/Inf_Sed_all_yearsWIDE.Rdat"))

# Tidy up ####
rm(list=ls(pattern="^df"))
rm(datfol,libfolder)

detach(package:tidyverse, unload = TRUE)
