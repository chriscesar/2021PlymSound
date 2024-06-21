# an.inf_post_model_digging.R #### 
# digging into the multivariate abundance models produced in the an.inf_v2.R script 

# Set up ####
source("R/metadata.R")
cbPalettetxt <- c("#994F00", "#0C7BDC", # colour palette for plots
                           "#d41159", "#009E73",
                           "#F0E442", "#0072B2",
                           "#D55E00", "#CC79A7")
                           ## load required packages ####
ld_pkgs <- c("tidyverse","vegan","lmerTest","rstatix", "mvabund","tictoc",
             "MASS","ggtext","ggpmisc", "gllvm") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)
tictoc::tic.clearlog()

# load & prep raw data ####
tic("Load & prep raw data")
# load data
df0 <- as_tibble(openxlsx::read.xlsx("data/in/Inf_Sed_all_years_Long_updated.xlsx",
                                     sheet = "Matched sites"))
df0 %>% rename("MatchedSite"=`2021.matched.site`) -> df0

# 1. filter out taxa to be 'discarded' and tidy up MatchedSite variable
df0 %>% 
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
  )) -> dfl0;toc(log=TRUE)

#2. Calculate mean by sampling event and widen
tic("DATA IMPORT: Calc means across reps, widen for analysis, export")
print("Calc means across reps, widen for analysis, export")
# Function to identify numeric columns containing only zero values
contains_only_zero <- function(x) {
  is_numeric <- is.numeric(x)
  all_zero <- all(x == 0)
  is_numeric && all_zero
}

dfl0 %>%
  group_by(across(-Abund)) %>% 
  summarise(Abund=mean(Abund),.groups = "drop") %>% #summarise abundances by station/year
  ungroup() %>% 
  filter(., Abund != 0) %>% #remove zero values
  filter(., !is.na(BSH)) %>% ## remove empty BSH values
  pivot_wider(names_from = Taxon.USE,
              values_from = Abund,values_fill = 0) %>% 
  dplyr::select(where(~!contains_only_zero(.))) %>% 
  ungroup() -> dfw # remove columns summing to zero

## remove 5.1 (only 2 samples)
dfw_trim <- dfw %>% filter(., BSH_CODE != "A5.1") %>% 
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0))

# bsh_data <- subset(dfw_trim, BSH_CODE == "A5.2")

## remove 'empty' columns
bsh_data <- dfw_trim %>%
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

# create data for ordination
bsh_data %>% 
  dplyr::select(-c(1:12)) -> bsh_dataord

# load model outputs ####
## A5.2 ####
dfw %>% 
  filter(.,BSH_CODE == "A5.2") %>% 
  dplyr::select(-c(1,2,4:12)) %>% 
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) -> dfA5.2## remove numeric variables that sum to 0


fit.glm_A5.2 <- readRDS("outputs/mvabund.inf.A5.2.rdat")
fit.glm.summary_A5.2 <- readRDS("outputs/mvabund.inf.A5.2.summary.rdat")
fit.glm.out_A5.2 <- readRDS("outputs/mvabund.inf.A5.2.pw.rdat")
fit.anosim_A5.2 <- readRDS("outputs/anosim.inf.A5.2.rdat")
fit.adonis2 <- readRDS("outputs/adonis2.inf.A5.2.rdat")
fit.simper_A5.2 <- readRDS("outputs/simper.inf.A5.2.rdat")
yr_A5.2 <- dfA5.2$Year
x_A5.2 <- mvabund(dfA5.2[,-c(1)])

## taxa per year:
print(paste0(ncol(x_A5.2)," unique taxa recorded in the A5.2 BSH across all sampling years"))

## by year:
dfA5.2 %>%
  pivot_longer(
    cols = -Year,
    names_to = "Species",
    values_to = "Abundance"
  ) %>% 
  filter(Abundance > 0) %>%
  group_by(Year) %>%
  summarise(SpeciesCount = n_distinct(Species))

(best_rsq_A5.2 <- mvabund::best.r.sq(x_A5.2 ~ yr_A5.2))

## A5.3 ####
dfw %>% 
  filter(.,BSH_CODE == "A5.3") %>% 
  dplyr::select(-c(1,2,4:12)) %>% 
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) -> dfA5.3## remove numeric variables that sum to 0


fit.glm_A5.3 <- readRDS("outputs/mvabund.inf.A5.3.rdat")
fit.glm.summary_A5.3 <- readRDS("outputs/mvabund.inf.A5.3.summary.rdat")
fit.glm.out_A5.3 <- readRDS("outputs/mvabund.inf.A5.3.pw.rdat")

yr_A5.3 <- dfA5.3$Year
x_A5.3 <- mvabund(dfA5.3[,-c(1)])

## taxa per year:
print(paste0(ncol(x_A5.3)," unique taxa recorded in the A5.3 BSH across all sampling years"))

## by year:
dfA5.3 %>%
  pivot_longer(
    cols = -Year,
    names_to = "Species",
    values_to = "Abundance"
  ) %>% 
  filter(Abundance > 0) %>%
  group_by(Year) %>%
  summarise(SpeciesCount = n_distinct(Species))


(best_rsq_A5.3 <- mvabund::best.r.sq(x_A5.3 ~ yr_A5.3))

## A5.4 ####
dfw %>% 
  filter(.,BSH_CODE == "A5.4") %>% 
  dplyr::select(-c(1,2,4:12)) %>% 
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) -> dfA5.4## remove numeric variables that sum to 0

fit.glm_A5.4 <- readRDS("outputs/mvabund.inf.A5.4.rdat")
fit.glm.summary_A5.4 <- readRDS("outputs/mvabund.inf.A5.4.summary.rdat")
fit.glm.out_A5.4 <- readRDS("outputs/mvabund.inf.A5.4.pw.rdat")

yr_A5.4 <- dfA5.4$Year
x_A5.4 <- mvabund(dfA5.4[,-c(1)])

## taxa per year:
print(paste0(ncol(x_A5.4)," unique taxa recorded in the A5.4 BSH across all sampling years"))

## by year:
dfA5.4 %>%
  pivot_longer(
    cols = -Year,
    names_to = "Species",
    values_to = "Abundance"
  ) %>% 
  filter(Abundance > 0) %>%
  group_by(Year) %>%
  summarise(SpeciesCount = n_distinct(Species))

(best_rsq_A5.4 <- mvabund::best.r.sq(x_A5.4 ~ yr_A5.4))

