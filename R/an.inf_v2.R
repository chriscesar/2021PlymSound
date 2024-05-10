# an.inf_v2.R #### 
# Import and analyse new version of infauna data
# Only comparing sites within Plymouth Sound

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

# Copy data:
tic("DATA IMPORT: Copying & loading data");print("Copying & loading data")
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
toc(log=TRUE)

tic("DATA IMPORT: Initial data tidy");print("Initial data tidy")
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

write.csv(dfw,file="outputs/infauna_wide.csv",row.names = FALSE)
toc(log=TRUE)

# 3. Run ordination
## create taxon-only data and feed into ordination
tic("DATA IMPORT: Run ordination across all years")
print("Run ordination across all years")
set.seed(22); dfw %>% 
  dplyr::select(.,-c(MatchedSite:BSH_CODE)) %>% metaMDS(., trymax = 200) -> ord
# plot(ord)## this version is across all years
toc(log=TRUE)

## create taxon-only data and feed into ordination
tic("DATA IMPORT: Run ordination for 2021")
print("Run ordination for 2021")
set.seed(22); dfw %>% 
  filter(.,Year==2021) %>% # retain current year only
  dplyr::select(where(~!contains_only_zero(.))) %>% 
  dplyr::select(.,-c(MatchedSite:BSH_CODE)) %>% 
  metaMDS(., trymax = 200) -> ord
plot(ord)## this version is for 2021 only
toc(log=TRUE)

tic("DATA IMPORT: Extract ordination data for plotting");print("Extract ordination data for plotting")
### plot ordination through ggplot
## extract Site scores
mds_scores <- as_tibble(as.data.frame(scores(ord,"site")))
# mds_scores$BSH <- dfw %>% 
#   filter(.,Year==2021) %>% # retain current year only
#   dplyr::select(where(~!contains_only_zero(.))) %>% 
#   dplyr::select(.,BSH_CODE)
mds_scores$BSH <- dfw$BSH_CODE[dfw$Year==2021]
mds_scores$Station <- dfw$MatchedSite[dfw$Year==2021]
# mds_scores$station <- dfw %>% 
#   filter(.,Year==2021) %>% # retain current year only
#   dplyr::select(where(~!contains_only_zero(.))) %>% 
#   dplyr::select(.,MatchedSite)

### extract species scores and groups
spp_scores <- as.data.frame(scores(ord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)

pdf("figs/inf_mds_all_2021.pdf", width=14,height = 14)
mds_scores %>% 
  ggplot(.,aes(x=NMDS1, y=NMDS2))+
  geom_text(data=spp_scores,
            aes(
              x=NMDS1,
              y= NMDS2,
              label=species_sh
              ),
            col="darkgrey",
            size=4, show.legend = FALSE, alpha = 0.5)+
  geom_text(aes(
    label=Station,
    colour=BSH
    ),
    size=8,
    fontface=2)+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(ord$stress, 3))))+
  scale_colour_manual(values = cbPalettetxt)+
  theme(
    legend.title = element_text(face=2),
    legend.text = element_text(face=2))
dev.off()

toc(log=TRUE)

# tidy
rm(mds_scores,ord,spp_scores)

# Temporal ordinations by BSH (FOR LOOP) ####
## remove 5.1 (only 2 samples)
tic("ANALYSES by BSH: Ordinations by BSH")
dfw_trim <- dfw %>% filter(., BSH_CODE != "A5.1") %>% 
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0))

# Loop through unique WB_Names
for(bshcode in unique(dfw_trim$BSH_CODE)) {
  # Subset data for current BSH
  bsh_data <- subset(dfw_trim, BSH_CODE == bshcode)
  
  ## remove 'empty' columns
  bsh_data <- bsh_data %>%
    dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0
  
  # create data for ordination
  bsh_data %>% 
    dplyr::select(-c(1:12)) -> bsh_dataord
  
  #Run ordination
  set.seed(80);bsh_dataordout <- metaMDS(bsh_dataord, try = 30, trymax = 200,k=2)
  
  ### ordination plot ####
  ### extract site scores and groups
  mds_scores <- as_tibble(as.data.frame(scores(bsh_dataordout,"site")))
  mds_scores$Year <- bsh_data$Year
  mds_scores$MatchedSite <- bsh_data$MatchedSite
  mds_scores$MatchedSiteYR <- bsh_data$MatchedSiteYr
  mds_scores$GEAR <- bsh_data$GEAR
  mds_scores$GRAVELPCT <- bsh_data$GRAVELPCT
  mds_scores$SANDPCT <- bsh_data$SANDPCT
  mds_scores$MUDPCT <- bsh_data$MUDPCT
  mds_scores$BSH <- bsh_data$BSH
  mds_scores$BSH_CODE <- bsh_data$BSH_CODE
  
  ### extract species scores and groups
  spp_scores <- as.data.frame(scores(bsh_dataordout, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
  spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
  spp_scores$species_sh <- make.cepnames(spp_scores$species)
  
  # generate plot ####
  ### amend to provide station names
  pl <- ggplot(mds_scores,
               aes(
                 x=NMDS1, y = NMDS2
               ))+
    geom_hline(yintercept = 0,lty=2, col="grey")+
    geom_vline(xintercept = 0,lty=2, col="grey")+
    geom_text(data=spp_scores, aes(x=NMDS1, y= NMDS2,label=species_sh),
              col="grey",
              alpha=0.5,
              size=6,
              inherit.aes = FALSE)+
    geom_text(aes(label=MatchedSite,
                  col=factor(Year)),
              size=8,
              fontface="bold")+
    scale_colour_manual(name = "Year",values=cbPalettetxt)+
    coord_fixed()+
    geom_text_npc(aes(npcx = .99, npcy = .99,
                      label=paste("Stress = ",
                                  round(bsh_dataordout$stress, 3))))+
    labs(x="nmMDS Axis 1",
         y="nmMDS Axis 2",
         caption=paste0(unique(bsh_data$BSH)[1]," infaunal data"))+
    theme(plot.caption = element_text(size=16,face="bold"),
          axis.title = element_text(size=14, face="bold"),
          legend.text = element_text(face=2),
          legend.title = element_text(face=2))
  
  ggsave(filename = paste0("figs/infauna_",unique(bsh_data$BSH)[1],"_mds.png"),
         width = 14, height = 14, units="in",plot=pl)
  rm(bsh_data,bsh_dataord,bsh_dataordout,mds_scores,pl,spp_scores,bshcode)
}
toc(log=TRUE)

## MVABUND ####
tic("ANALYSES by BSH: Run gllvm models by BSH")
dfw_trim <- dfw %>% filter(., BSH_CODE != "A5.1") %>% 
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0))

for (bshcode in unique(dfw_trim$BSH_CODE)) {
  
  # subset data by BSH
  bsh_data <- subset(dfw_trim, BSH_CODE == bshcode)
  
  ## remove 'empty' columns
  bsh_data <- bsh_data %>%
    dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0
  
  # create data for ordination
  bsh_data %>% 
    dplyr::select(-c(1:12)) -> bsh_dataord
    # mvabund::mvabund(.)-> bsh_dataord
  
  ### produce meanvar plot
  png(file = paste0("figs/infMeanVar.",unique(bsh_data$BSH)[1],".png"),
      width=12*ppi, height=6*ppi, res=ppi)
  mvpl <- mvabund::meanvar.plot(mvabund(bsh_dataord),
                                xlab="Mean",
                                ylab="Variance",
                                table=TRUE)
  
  # Step 1: Find the minimum and maximum values
  min_value <- min(mvpl[,2])
  max_value <- max(mvpl[,2])
  
  min_order <- floor(log10(min_value))
  max_order <- floor(log10(max_value))
  orders_of_magnitude_covered <- max_order - min_order
  
  ttl <- paste0("Very strong mean-variance relationship in invertebrate abundances in the ", unique(bsh_data$BSH)[1]," broadscale habitat")
  sbtt <- paste0("Variance within the dataset covers *",orders_of_magnitude_covered," orders of magnitude*.")
  
  mtext(side=3, line = 1, at =-0.07, adj=0, cex = 1, ttl, font=1)
  mtext(side=3, line = 0.25, at =-0.07, adj=0, cex = 0.7, sbtt)
  
  dev.off()
  rm(min_order,max_order,mvpl,min_value,max_value,orders_of_magnitude_covered,ttl,sbtt)
  
  ## run model
  fit.glm <- manyglm(mvabund::as.mvabund(bsh_dataord) ~ bsh_data$Year, family = "negative.binomial")
  fit.glm.summary <- summary(fit.glm)
  saveRDS(fit.glm, file = paste0("outputs/mvabund.inf.",unique(bsh_data$BSH_CODE)[1],".rdat"))
  saveRDS(fit.glm.summary,file=paste0("outputs/mvabund.inf.",unique(bsh_data$BSH_CODE)[1],".summary.rdat"))
  fit.glm.out <- mvabund::anova.manyglm(fit.glm,p.uni = "adjusted")
  saveRDS(fit.glm.out, file = paste0("outputs/mvabund.inf.",unique(bsh_data$BSH_CODE)[1],".pw.rdat"))
  
  m2tmp1 <- t(as.data.frame(fit.glm.out$uni.p))[,2]
  names(m2tmp1) <- names(bsh_dataord)
  
  # m2tmp1[m2tmp1<0.056]; range(m2tmp1[m2tmp1<0.051])
  print(paste0(length(names(m2tmp1[m2tmp1<0.056]))," 'significant' taxa"))
  
  m2tx <- names(m2tmp1[m2tmp1<0.056])#which taxa are 'significantly' different?
  kptx <- names(bsh_dataord) %in% m2tx
  m2tx <- bsh_dataord[, kptx]
  m2tx$Year <- bsh_data$Year
  ##make long
  m2txl <- m2tx %>% 
    relocate(Year) %>% #move Year to start
    pivot_longer(cols = c(2:ncol(.)))
  
  ggplot(m2txl)+
    geom_boxplot(aes(x=name,y=value),varwidth = TRUE)+
    facet_wrap(.~Year)+
    coord_flip()+
    labs(y="Taxon abundance",
         caption=paste0(unique(bsh_data$BSH)[1]," BSH"))+
    theme(axis.title.y = element_blank(),
          # axis.text.x = element_blank(),
          strip.text = element_text(face="bold")) -> pl2
  
  m2txl$Year <- as.factor(m2txl$Year)
  (ggplot(m2txl,aes(x=log(value+1),y=name,
                    fill=Year,
                    # colour=Year,stroke=1.5,
                    shape = Year))+
      geom_jitter(data=m2txl[,c(2:3)], inherit.aes = FALSE,
                  aes(x=log(value+1),y=name,),
                  height = 0.05,size=1, alpha = 0.5, colour = "grey") +
      geom_jitter(height = 0.05,size=3, alpha = 0.9) +
      scale_shape_manual(values = c(21:24))+
      # scale_shape_manual(values = c(1:3))+
      labs(title = paste0(unique(bsh_data$BSH)[1]," BSH"),
           x="log(Taxon abundance (n+1))",
           caption="Each facet displays the abundance data for each of the taxa which showed significant differences by year.
           Taxon abundances for a given year are displayed by larger, coloured icons.")+
      scale_fill_manual(values = cbPalette)+
      scale_colour_manual(values = cbPalette)+
      facet_wrap(.~Year)+
      theme(
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12,face="italic"),
        axis.title.x = element_text(face="bold"),
        strip.text.x = element_text(face="bold",size=12),
        # plot.caption = element_text(face="bold",size=12),
        plot.title.position = "plot",
        plot.title = element_text(face="bold",size=14)
        ) -> pl3)
  
  ggsave(filename = paste0("figs/infauna_",unique(bsh_data$BSH)[1],"_relabund.pdf"),
         width = 14, height = 6, units="in",plot=pl2)
  ggsave(filename = paste0("figs/infauna_",unique(bsh_data$BSH)[1],"_relabund_ver2.pdf"),
         width = 14, height = 6, units="in",plot=pl3)
}

toc(log=TRUE)

unlist(tic.log())

##############################################################


rm(wtmp,wtmpord,mv_wtmpord,m2,m2out,m2tx,m2txl,pl2,pl3)

