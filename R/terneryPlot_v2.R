### terneryPlot_v2.R ####
# produce ternery plot for sediment data ####

# Set up ####
source("R/metadata.R")
cbPalettetxt <- c("#994F00", "#0C7BDC", # colour palette for plots
                           "#d41159", "#009E73",
                           "#F0E442", "#0072B2",
                           "#D55E00", "#CC79A7")
                           ## load required packages ####
ld_pkgs <- c("grid","tidyverse","tictoc","ggtern") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)
tictoc::tic.clearlog()

tic("Import and prep data")
### load data
df0 <- as.data.frame(read_csv("outputs/infauna_wide.csv"))

## throw away species data
df0 %>% 
  dplyr::select(c(MatchedSite:BSH_CODE)) -> df

df %>% 
  mutate(Gravel = GRAVELPCT/100,
         Sand=SANDPCT/100,
         Mud=MUDPCT/100) %>% 
  filter(., !is.na(BSH_CODE))-> df
toc(log=TRUE)

# Plot by BSH
tic("Plot by BSH")

for (bshcode in unique(df$BSH_CODE)) {
  # subset data by BSH
  bsh_data <- subset(df, BSH_CODE == bshcode)
  
  tic(paste0(unique(bsh_data$BSH_CODE)[1], " prep data"))

  # png(file=paste0("figs/PSA_BSH_.",unique(bsh_data$BSH)[1],".png"),
  #     width=12*ppi, height=12*ppi, res=ppi)
  ggplot(data=bsh_data,
       aes(Mud,
           Gravel,
           Sand,
           colour=as.factor(Year))) +
    coord_tern()+
    theme_bw() +
    theme_showarrows() + custom_percent('%') + 
    geom_mask() +
    geom_text(aes(label=MatchedSiteYr),show.legend = F,
              fontface=2) + 
    labs(Tarrow = "Gravel",
         Larrow = "Mud",
         Rarrow = "Sand")+
    labs(title = paste0(unique(bsh_data$BSH)[1]," BSH"),
         caption = "Sample identities coloured by sample year")+
    theme(axis.title = element_text(face=2),
          plot.title = element_text(hjust=0.5,vjust=0,face="bold"))
  ggsave(filename = paste0("figs/PSA_BSH_.",unique(bsh_data$BSH)[1],".png"),
         device = "png",
         width=10, height=10, units = "in")
  # dev.off()
  toc(log=TRUE)
}

toc(log=TRUE)
unlist(tic.log())
