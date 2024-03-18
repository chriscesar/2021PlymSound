# an.divddv.01.R #### 
# Import and analyse DDV & dive data

# Set up ####
source("R/metadata.R")

## load required packages ####
ld_pkgs <- c("tidyverse","vegan","lmerTest","rstatix", "mvabund",
             "MASS","ggtext","ggpmisc") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

# set data folder:
df0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,"Data/2021 data/DDV/StillDiveLong.xlsx"), sheet = "StillDiveLong"))

# throw away 'Remove'-flagged rows
df <- df0 %>% 
  filter(Flag == "0") %>% # keep flagged as 0
  dplyr::select(-Flag) %>% #remove Flag column
  filter(Abund != 0) #remove Abundance = 0 rows

minval <- min(df$Abund)

### summarise to remove 'duplicate' taxon names
#'hi res' taxa
dfw <- df %>% 
  # mutate(Abund = ifelse(Abund != 0,
  #                       as.integer(Abund/min(Abund[Abund != 0])),0)) %>% #OPTIONAL: convert 'continuous' data into integers
  dplyr::select(-c(Taxon:Species,MTG:MTG.level)) %>% #remove taxon info
  group_by(Data.type,Year,Dive.Site,Station.Code,Sample.Ref,Easting,Northing,
           `Depth.(m)`,Quality,Kelp,MNCR.Code,BSH, TaxonUSE) %>% 
  summarise(Abund = sum(Abund)) %>% 
  pivot_wider(names_from = TaxonUSE, values_from = Abund,
              values_fill = 0) %>% ungroup()

#'Maj Taxonomic group' taxa
dfwMTG <- df %>% 
  # dplyr::mutate(Abund = ifelse(Abund != 0,
  #                       as.integer(Abund/min(Abund[Abund != 0])),0)) %>% #OPTIONAL: convert 'continuous' data into integers
  dplyr::select(-c(Taxon:TaxonUSE,MTG.level)) %>% #remove taxon info
  group_by(Data.type,Year,Dive.Site,Station.Code,Sample.Ref,Easting,Northing,
           `Depth.(m)`,Quality,Kelp,MNCR.Code,BSH, MTG) %>% 
  summarise(Abund = sum(Abund)) %>% 
  pivot_wider(names_from = MTG, values_from = Abund,
              values_fill = 0) %>% ungroup()

  
### tweak Dive Site values
dfw$Dive.Site <- ifelse(dfw$Dive.Site == "Duke Rock S ","Duke Rock S",dfw$Dive.Site)
dfwMTG$Dive.Site <- ifelse(dfwMTG$Dive.Site == "Duke Rock S ","Duke Rock S",dfwMTG$Dive.Site)

## export data
#write.csv(dfw, file="data_in/diveDDVWide.csv")

##########################FROM HERE#################################
# look into removing samples containing few taxa
##############################################################################
### keep only 2021 data
dfw <- dfw[dfw$Year==2021,]
dfwMTG <- dfwMTG[dfwMTG$Year==2021,]

# 'hi res' ####
## calculate taxon richnesses and remove samples with low S from ordination data
dfw$S <- vegan::specnumber(dfw[,-c(1:12)])
min_sp <- 3 ### set minimum taxon richness for MDS plotting
keeps <- dfw$S>min_sp
dfw_ord <- dfw[keeps,]

### replace NA values with zero
dfwtmp <- dfw_ord %>% 
  dplyr::select(-c(Data.type:BSH)) %>% 
  replace(is.na(.), 0)

### run ordinations
set.seed(80);ord <- metaMDS(dfwtmp, try = 30, trymax = 100)
#plot(ord,type="t")

## make neater version using ggplot ####
### extract site scores and groups
mds_scores <- as_tibble(as.data.frame(scores(ord,"site")))
mds_scores$Type <- dfw[keeps,]$Data.type
mds_scores$Year <- dfw[keeps,]$Year
mds_scores$Dive.Site <- dfw[keeps,]$Dive.Site
mds_scores$Station.Code <- dfw[keeps,]$Station.Code
mds_scores$Sample.Ref <- dfw[keeps,]$Sample.Ref
mds_scores$`Depth.(m)` <- dfw[keeps,]$`Depth.(m)`
mds_scores$Quality <- dfw[keeps,]$Quality
mds_scores$Kelp <- dfw[keeps,]$Kelp
mds_scores$MNCR.Code <- dfw[keeps,]$MNCR.Code
mds_scores$BSH <- dfw[keeps,]$BSH

### extract species scores and groups
spp_scores <- as.data.frame(scores(ord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)
head(spp_scores)  #look at the data

### create centroids by BSH
centroid_DiveSite <- mds_scores %>% 
  group_by(Dive.Site) %>%
  summarise(NMDS1=mean(NMDS1),NMDS2=mean(NMDS2))## 2d version

### create link to BSH centroids
mds_scores <- mds_scores %>% group_by(Dive.Site) %>% 
  mutate(BSH_centroid1=mean(NMDS1),
         BSH_centroid2=mean(NMDS2)) %>% ungroup ##2d version

# generate plot ####
pl <- ggplot(mds_scores,
             aes(x=NMDS1,y=NMDS2,colour=Dive.Site,shape=Type,
                 xend=BSH_centroid1,
                 yend=BSH_centroid2))+  
  geom_hline(yintercept = 0,lty=2, col="grey")+
  geom_vline(xintercept = 0,lty=2, col="grey")+
  geom_text(data=spp_scores, aes(x=NMDS1, y= NMDS2,label=species_sh),
            # col=Phylum),
            col="grey",
            alpha=0.5,
            inherit.aes = FALSE)+
  geom_segment(show.legend = FALSE,alpha=0.6)+## add 'spider legs' joining samples-centroids
  geom_point(show.legend = TRUE,colour=1,aes(fill=Dive.Site),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_DiveSite,aes(x=NMDS1,y=NMDS2,label=Dive.Site,fill=Dive.Site),
    width = unit(0.05, "npc"),
    inherit.aes = FALSE,show.legend = FALSE,hjust=0.5,
    fontface="bold")+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  coord_fixed()+
  labs(x="nmMDS Axis 1",
       y="nmMDS Axis 2")+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(ord$stress, 3))))

pdf("figs/inf_mds.pdf", width=14,height = 14)
print(pl)
dev.off()

# does S differ by survey type? ####
ggplot(dfw, aes(x = Data.type, y = S)) +
  geom_boxplot(varwidth = TRUE, outlier.colour = NA)+
  geom_jitter(height = 0, width = 0.25, alpha = 0.25,show.legend = TRUE,
              size = 3,
              aes(colour = Dive.Site, fill=Dive.Site))+
  facet_wrap(.~Dive.Site)

# 'MTG' ####
## calculate taxon richnesses and remove samples with low S from ordination data
dfwMTG$S <- vegan::specnumber(dfwMTG[,-c(1:12)])
min_sp <- 3 ### set minimum taxon richness for MDS plotting
keeps <- dfwMTG$S>min_sp
dfwMTG_ord <- dfwMTG[keeps,]

### replace NA values with zero
dfwtmp <- dfwMTG_ord %>% 
  dplyr::select(-c(Data.type:BSH,S)) %>% 
  replace(is.na(.), 0)

### run ordinations
set.seed(80);ord <- metaMDS(dfwtmp, try = 30, trymax = 200)
#plot(ord,type="t")

## make neater version using ggplot ####
### extract site scores and groups
mds_scores <- as_tibble(as.data.frame(scores(ord,"site")))
mds_scores$Type <- dfwMTG[keeps,]$Data.type
mds_scores$Year <- dfwMTG[keeps,]$Year
mds_scores$Dive.Site <- dfwMTG[keeps,]$Dive.Site
mds_scores$Station.Code <- dfwMTG[keeps,]$Station.Code
mds_scores$Sample.Ref <- dfwMTG[keeps,]$Sample.Ref
mds_scores$`Depth.(m)` <- dfwMTG[keeps,]$`Depth.(m)`
mds_scores$Quality <- dfwMTG[keeps,]$Quality
mds_scores$Kelp <- dfwMTG[keeps,]$Kelp
mds_scores$MNCR.Code <- dfwMTG[keeps,]$MNCR.Code
mds_scores$BSH <- dfwMTG[keeps,]$BSH

### extract species scores and groups
spp_scores <- as.data.frame(scores(ord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)
head(spp_scores)  #look at the data

### create centroids by BSH
centroid_DiveSite <- mds_scores %>% 
  group_by(Dive.Site) %>%
  summarise(NMDS1=mean(NMDS1),NMDS2=mean(NMDS2))## 2d version

### create link to BSH centroids
mds_scores <- mds_scores %>% group_by(Dive.Site) %>% 
  mutate(BSH_centroid1=mean(NMDS1),
         BSH_centroid2=mean(NMDS2)) %>% ungroup ##2d version

# generate plot ####
pl <- ggplot(mds_scores,
             aes(x=NMDS1,y=NMDS2,colour=Dive.Site,shape=Type,
                 xend=BSH_centroid1,
                 yend=BSH_centroid2))+  
  geom_hline(yintercept = 0,lty=2, col="grey")+
  geom_vline(xintercept = 0,lty=2, col="grey")+
  geom_text(data=spp_scores, aes(x=NMDS1, y= NMDS2,label=species_sh),
            # col=Phylum),
            col="grey",
            alpha=0.5,
            inherit.aes = FALSE)+
  geom_segment(show.legend = FALSE,alpha=0.6)+## add 'spider legs' joining samples-centroids
  geom_point(show.legend = TRUE,colour=1,aes(fill=Dive.Site),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_DiveSite,aes(x=NMDS1,y=NMDS2,label=Dive.Site,fill=Dive.Site),
    width = unit(0.05, "npc"),
    inherit.aes = FALSE,show.legend = FALSE,hjust=0.5,
    fontface="bold")+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  coord_fixed()+
  labs(x="nmMDS Axis 1",
       y="nmMDS Axis 2")+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(ord$stress, 3))))

pdf("figs/inf_mds.pdf", width=14,height = 14)
print(pl)
dev.off()

# does S differ by survey type? ####
ggplot(dfwMTG, aes(x = Data.type, y = S)) +
  geom_boxplot(varwidth = TRUE, outlier.colour = NA)+
  geom_jitter(height = 0, width = 0.25, alpha = 0.25,show.legend = FALSE,
              size = 3,
              aes(colour = Dive.Site, fill=Dive.Site))+
  facet_wrap(.~Dive.Site)

# ANOSIM: statistical comparison of MTG assemblages ####
(anoMTG2021 <- adonis2(dfwtmp ~ dfwMTG_ord$Data.type,
                            permutations = perm))

# MVABUND: statistical comparison of MTG assemblages ####
mv_dfwtmp <- mvabund(dfwtmp)#create mvabund object
mvabund::meanvar.plot(mv_dfwtmp)# mean-variance plot

## poisson:
system.time(mod1 <- manyglm(mv_dfwtmp ~ dfwMTG_ord$Data.type,
                                    family="poisson"))
summary(mod1)
plot(mod1)

## negative binomial
system.time(mod2 <- manyglm(mv_dfwtmp ~ dfwMTG_ord$Data.type*dfwMTG_ord$Dive.Site,
                            family="negative_binomial"))
summary(mod2);plot(mod2)
# system.time(anova_mod2 <- mvabund::anova.manyglm(mod2,p.uni = "adjusted"))
# saveRDS(anova_mod2,file="outputs/mvabund.MTG.2021.rds")#proportion values
saveRDS(anova_mod2,file="outputs/mvabund.MTG.2021.integer.rds")#integer version
anova_mod2 <- readRDS("outputs/mvabund.MTG.2021.rds")

###try different distributions
# system.time(mod3 <- manyany)