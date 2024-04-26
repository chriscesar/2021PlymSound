# an.inf_v2.R #### 
# Import and analyse new version of infauna data
# Only comparing sites within Plymouth Sound

# Set up ####
source("R/metadata.R")

## load required packages ####
ld_pkgs <- c("tidyverse","vegan","lmerTest","rstatix", "mvabund","tictoc",
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
tic()
df0 <- as_tibble(openxlsx::read.xlsx("data/in/Inf_Sed_all_years_Long_updated.xlsx",
                                     sheet = "Matched sites"))
toc()

df0 %>% rename("MatchedSite"=`2021.matched.site`) -> df0

# 1. filter out taxa to be 'discarded' and tidy up MatchedSite variable
tic();df0 %>% 
  # change taxa flagged as -999 to values of 1. Consider removal of these taxa if required
  mutate(Abund = ifelse(Abund == -999, 1, Abund)) %>%
  # filter(.,Abund != 0) %>% #remove 'zero' abundances
  mutate(MatchedSite = paste0("PL", substr(MatchedSite, 6, 7))) %>% #rename MatchedSite contents
  mutate(MatchedSiteYr = paste0(MatchedSite, "_", Year)) %>% #append Year value
  relocate(MatchedSiteYr, .after = MatchedSite) %>% #move new variable somewhere sensible
  #remove taxa flagged for removal
  filter(.,!str_detect(Flag, "^Remove")) %>% 
  # remove potentially problematic variables
  dplyr::select(., -c(
    Flag, # "Flag" (now empty)
    Site_Yr_mesh, #superseded by MatchedSite
    Site_Yr, #superseded by MatchedSite
    Site, #superseded by MatchedSite
    Station.Code,#superseded by MatchedSite
    Mesh, #not needed
    Notes, #not needed
    Latitude, #potential issues with condensing data
    Longitude, #potential issues with condensing data
    OSGB36.Easting,#potential issues with condensing data
    OSGB36.Northing,#potential issues with condensing data
    Easting,#potential issues with condensing data
    Northing,#potential issues with condensing data
    Taxon, # not needed
    Area
    )) -> dfl0;toc()

#2. Calculate mean by sampling event and widen

# Function to identify numeric columns containing only zero values
contains_only_zero <- function(x) {
  is_numeric <- is.numeric(x)
  all_zero <- all(x == 0)
  is_numeric && all_zero
}

tic();dfl0 %>%
  group_by(across(-Abund)) %>% 
  summarise(Abund=mean(Abund),.groups = "drop") %>% #summarise abundances by station/year
  ungroup() %>% 
  filter(., Abund != 0) %>% #remove zero values
  filter(., !is.na(BSH)) %>% ## remove empty BSH values
  pivot_wider(names_from = Taxon.USE,
              values_from = Abund,values_fill = 0) %>% 
  dplyr::select(where(~!contains_only_zero(.))) -> dfw;toc() # remove columns summing to zero

# 3. Run ordination
## create taxon-only data and feed into ordination
set.seed(22); dfw %>% 
  dplyr::select(.,-c(MatchedSite:BSH_CODE)) %>% metaMDS(., trymax = 200) -> ord
plot(ord)## this version is across all years

## create taxon-only data and feed into ordination
set.seed(22); dfw %>% 
  filter(.,Year==2021) %>% # retain current year only
  dplyr::select(where(~!contains_only_zero(.))) %>% 
  dplyr::select(.,-c(MatchedSite:BSH_CODE)) %>% 
  metaMDS(., trymax = 200) -> ord
plot(ord)## this version is for 2021 only

### plot ordination through ggplot
## extract Site scores
mds_scores <- as_tibble(as.data.frame(scores(ord,"site")))
mds_scores$BSH <- dfw %>% 
  filter(.,Year==2021) %>% # retain current year only
  dplyr::select(where(~!contains_only_zero(.))) %>% 
  dplyr::select(.,BSH_CODE)
mds_scores$station <- dfw %>% 
  filter(.,Year==2021) %>% # retain current year only
  dplyr::select(where(~!contains_only_zero(.))) %>% 
  dplyr::select(.,MatchedSite)

### extract species scores and groups
spp_scores <- as.data.frame(scores(ord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)

mds_scores %>% 
  ggplot(.,aes(x=NMDS1, y=NMDS2))+
  geom_text(data=spp_scores,
            aes(x=NMDS1, y= NMDS2,label=species_sh),
            col=2,
            size=2)+
  geom_text(aes(
    # label=station$MatchedSite
    label=BSH$BSH_CODE
    ),
            fontface=2)+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(ord$stress, 3))))
