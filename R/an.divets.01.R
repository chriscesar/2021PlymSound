# an.divets.01.R #### 
# Import and analyse dive time series data

# Set up ####
source("R/datfol.R")
source("R/metadata.R")

## load packages ####
ld_pkgs <- c("tidyverse","vegan","lmerTest","rstatix", "mvabund",
             "MASS","ggtext","ggpmisc") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

# import data ####
df0 <- as_tibble(openxlsx::read.xlsx(paste0(datfol,"Data/2021 data/DDV/StillDiveLong.xlsx"), sheet = "StillDiveLong"))

# throw away 'Remove'-flagged rows
df <- df0 %>% 
  filter(Flag == "0") %>% # keep flagged as 0
  dplyr::select(-Flag) %>% #remove Flag column
  filter(Abund != 0) #remove Abundance = 0 rows

### retain only 'Dive' Data
df <- df %>% 
  filter(Data.type == "Dive")

### summarise to remove 'duplicate' taxon names
#'hi res' taxa
dfw <- df %>% 
  dplyr::select(-c(Taxon:Species,MTG:MTG.level)) %>% #remove taxon info
  group_by(Data.type,Year,Dive.Site,Station.Code,Sample.Ref,Easting,Northing,
           `Depth.(m)`,Quality,Kelp,MNCR.Code,BSH, TaxonUSE) %>% 
  summarise(Abund = sum(Abund),.groups = "drop") %>% #sum duplicate taxon names in each sample
  pivot_wider(names_from = TaxonUSE, values_from = Abund, #create 'wide' version
              values_fill = 0) %>% ungroup()

### tweak Dive Site values
dfw$Dive.Site <- ifelse(dfw$Dive.Site == "Duke Rock S ","Duke Rock S",dfw$Dive.Site)

# ### export wide data
# write.csv(dfw, file="data_in/diveWide.csv")

# initial comparisons ####

## calculate taxon richnesses compare between years
dfw$S <- vegan::specnumber(dfw[,-c(1:12)])
dfw$Year <- as.factor(dfw$Year)
summary(m01 <- lm(S ~ Year, data=dfw))
visreg::visreg(m01)

dfw %>% ggplot(., aes(S, fill = Year))+
  geom_histogram(alpha=0.5, col=1, show.legend = FALSE)+
  # geom_density(alpha=0.2, show.legend = FALSE)+
  # facet_wrap(. ~ Year, ncol = 1)+
  # facet_grid(Year~.)+
  facet_grid(Year~Dive.Site)+
  scale_fill_manual(values=cbPalette)+
  labs(x = "Taxon richness", y = "Taxon richness",
  # labs(x = "Taxon richness", y = "Density",
       title = "Density of taxon richness values",
       subtitle = "Reported in dive surveys over four years")+
  theme(strip.text.y = element_text(face="bold"))

## remove samples with low S from ordination data
min_sp <- 0 ### set minimum taxon richness for MDS plotting (1-3 causes disparity)
keeps <- dfw$S>(min_sp-1)
dfw$S <- NULL
dfw_ord <- dfw[keeps,]

### replace NA values with zero
dfwtmp <- dfw_ord %>% 
  dplyr::select(-c(Data.type:BSH)) %>% 
  replace(is.na(.), 0)

### run ordinations
set.seed(80);ord <- metaMDS(dfwtmp, try = 30, trymax = 200)
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

### create centroids by Year
centroid_DiveSite <- mds_scores %>% 
  group_by(Year) %>%
  summarise(NMDS1=mean(NMDS1),NMDS2=mean(NMDS2))## 2d version

### create link to Year centroids
mds_scores <- mds_scores %>% group_by(Year) %>% 
  mutate(BSH_centroid1=mean(NMDS1),
         BSH_centroid2=mean(NMDS2)) %>% ungroup ##2d version

# generate plot ####
pl <- ggplot(mds_scores,
             aes(x=NMDS1,y=NMDS2,colour=Year,shape=Year,
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
  geom_point(show.legend = FALSE,colour=1,aes(fill=Year),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_DiveSite,aes(x=NMDS1,y=NMDS2,label=Year,fill=Year),
    width = unit(0.05, "npc"),
    inherit.aes = FALSE,show.legend = FALSE,hjust=0.5,
    fontface="bold")+
  scale_shape_manual(values = c(21:24))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  coord_fixed()+
  labs(x="nmMDS Axis 1",
       y="nmMDS Axis 2")+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(ord$stress, 3))))
pl
rm(pl, centroid_DiveSite,mds_scores,spp_scores, ord)

# ANOSIM: statistical comparison of Dive years ####
(anoDive <- adonis2(dfwtmp ~ dfw_ord$Year,
                    permutations = perm))

# MVABUND: statistical comparison of MTG assemblages ####
mv_dfwtmp <- mvabund(dfwtmp)#create mvabund object
mvabund::meanvar.plot(mv_dfwtmp)# mean-variance plot shows STRONG mean-var relationship

# split data into chunks by dive site ####
dfwDvDvPnt <- dfw[dfw$Dive.Site == "Devils Point",]
dfwDvDkRck <- dfw[dfw$Dive.Site == "Duke Rock S",]
dfwDvENDGB <- dfw[dfw$Dive.Site == "E N DG Buoy",]
dfwDvEKng <- dfw[dfw$Dive.Site == "Eastern Kings",]
dfwDvMSW <- dfw[dfw$Dive.Site == "Mew Stone W",]

## Devils Point ####
wtmp <- dfwDvDvPnt

## remove 'empty' columns
wtmp <- wtmp %>% 
  # dplyr::select(where(is.character)|across(wher(is.numeric),-sum(.) !=0))
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

# create data for ordination
wtmp %>% 
  dplyr::select(-c(1:12)) -> wtmpord

## how many unique taxa?
print(paste0(ncol(wtmpord)," unique taxa"))

## by year(?)
vegan::specnumber(wtmpord, groups=wtmp$Year)

## run ordination ##
set.seed(80);tmpord <- metaMDS(wtmpord, try = 30, trymax = 200)
plot(tmpord)

### ordination plot using ggplot ####
### extract site scores and groups
mds_scores <- as_tibble(as.data.frame(scores(tmpord,"site")))
mds_scores$Year <- wtmp$Year
mds_scores$Station.Code <- wtmp$Station.Code
mds_scores$Sample.Ref <- wtmp$Sample.Ref
mds_scores$`Depth.(m)` <- wtmp$`Depth.(m)`
mds_scores$Quality <- wtmp$Quality
mds_scores$Kelp <- wtmp$Kelp
mds_scores$MNCR.Code <- wtmp$MNCR.Code
mds_scores$BSH <- wtmp$BSH

### extract species scores and groups
spp_scores <- as.data.frame(scores(tmpord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)
head(spp_scores)  #look at the data

### create centroids by Year
centroid_Site <- mds_scores %>% 
  group_by(Year) %>%
  summarise(NMDS1=mean(NMDS1),NMDS2=mean(NMDS2))## 2d version

### create link to Year centroids
mds_scores <- mds_scores %>% group_by(Year) %>% 
  mutate(BSH_centroid1=mean(NMDS1),
         BSH_centroid2=mean(NMDS2)) %>% ungroup ##2d version

### generate plot ####
pl <- ggplot(mds_scores,
             aes(x=NMDS1,y=NMDS2,colour=Year,shape=Year,
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
  geom_point(show.legend = FALSE,colour=1,aes(fill=Year),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_Site,aes(x=NMDS1,y=NMDS2,label=Year,fill=Year),
    width = unit(0.05, "npc"),
    inherit.aes = FALSE,show.legend = FALSE,
    hjust=0.5,vjust=0.5,halign=0.5,valign=0.5,
    fontface="bold")+
  scale_shape_manual(values = c(21:24))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  coord_fixed()+
  labs(x="nmMDS Axis 1",
       y="nmMDS Axis 2")+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(tmpord$stress, 3))))+
  labs(caption=paste0(unique(wtmp$Dive.Site)," dive site"))
pl
# pdf(paste0("figs/dive_",unique(wtmp$Dive.Site),"_mds.pdf"), width=14,height = 14)
ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_mds.pdf"),
       width = 14, height = 14, units="in",plot=pl)
rm(pl,mds_scores,centroid_Site,spp_scores)

## analyse for sig diffs between years ####
# (m1 <- adonis2(wtmpord~wtmp$Year))
mv_wtmpord <- mvabund::as.mvabund(wtmpord)
mvabund::meanvar.plot(mv_wtmpord)
m2 <- manylm(mv_wtmpord ~ wtmp$Year);summary(m2)
m2out <- mvabund::anova.manylm(m2,p.uni = "adjusted")
m2tmp1 <- t(as.data.frame(m2out$uni.p))[,2]
names(m2tmp1) <- names(wtmpord)
m2tmp1[m2tmp1<0.056]

m2tx <- names(m2tmp1[m2tmp1<0.056])#which taxa are 'significantly' different?
kptx <- names(wtmpord) %in% m2tx
m2tx <- wtmpord[, kptx]
m2tx$Year <- wtmp$Year
##make long
m2txl <- m2tx %>% 
  relocate(Year) %>% #move Year to start
  pivot_longer(cols = c(2:ncol(.)))

ggplot(m2txl)+
  geom_boxplot(aes(x=name,y=value),varwidth = TRUE)+
  facet_wrap(.~Year)+
  coord_flip()+
  labs(y="Relative abundance",
       caption=paste0(unique(wtmp$Dive.Site)," dive site\nAbundances for each year on the same relative scale"))+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) ->pl2

ggplot(m2txl,aes(x=value,y=name, fill=Year, shape = Year))+
  geom_jitter(height = 0.15,size=3, alpha = 0.4) +
  scale_shape_manual(values = c(21:24))+
  labs(x="Relative abundance",
       caption=paste0(unique(wtmp$Dive.Site)," dive site\nAbundances for each year on the same relative scale"))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl3

ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_relabund.pdf"),
       width = 14, height = 6, units="in",plot=pl2)
ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_relabund_ver2.pdf"),
       width = 14, height = 6, units="in",plot=pl3)
rm(wtmp,wtmpord,m1,mv_wtmpord,m2,m2out,m2tx,m2txl,pl2,pl3)

##########
## Duke Rock ####
wtmp <- dfwDvDkRck

## remove 'empty' columns
wtmp <- wtmp %>% 
  # dplyr::select(where(is.character)|across(wher(is.numeric),-sum(.) !=0))
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

# create data for ordination
wtmp %>% 
  dplyr::select(-c(1:12)) -> wtmpord

## how many unique taxa?
print(paste0(ncol(wtmpord)," unique taxa"))

## by year(?)
vegan::specnumber(wtmpord, groups=wtmp$Year)

## run ordination ##
set.seed(80);tmpord <- metaMDS(wtmpord, try = 30, trymax = 200)
plot(tmpord)

### ordination plot using ggplot ####
### extract site scores and groups
mds_scores <- as_tibble(as.data.frame(scores(tmpord,"site")))
mds_scores$Year <- wtmp$Year
mds_scores$Station.Code <- wtmp$Station.Code
mds_scores$Sample.Ref <- wtmp$Sample.Ref
mds_scores$`Depth.(m)` <- wtmp$`Depth.(m)`
mds_scores$Quality <- wtmp$Quality
mds_scores$Kelp <- wtmp$Kelp
mds_scores$MNCR.Code <- wtmp$MNCR.Code
mds_scores$BSH <- wtmp$BSH

### extract species scores and groups
spp_scores <- as.data.frame(scores(tmpord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)
head(spp_scores)  #look at the data

### create centroids by Year
centroid_Site <- mds_scores %>% 
  group_by(Year) %>%
  summarise(NMDS1=mean(NMDS1),NMDS2=mean(NMDS2))## 2d version

### create link to Year centroids
mds_scores <- mds_scores %>% group_by(Year) %>% 
  mutate(BSH_centroid1=mean(NMDS1),
         BSH_centroid2=mean(NMDS2)) %>% ungroup ##2d version

### generate plot ####
pl <- ggplot(mds_scores,
             aes(x=NMDS1,y=NMDS2,colour=Year,shape=Year,
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
  geom_point(show.legend = FALSE,colour=1,aes(fill=Year),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_Site,aes(x=NMDS1,y=NMDS2,label=Year,fill=Year),
    width = unit(0.05, "npc"),
    inherit.aes = FALSE,show.legend = FALSE,
    hjust=0.5,vjust=0.5,halign=0.5,valign=0.5,
    fontface="bold")+
  scale_shape_manual(values = c(21:24))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  coord_fixed()+
  labs(x="nmMDS Axis 1",
       y="nmMDS Axis 2")+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(tmpord$stress, 3))))+
  labs(caption=paste0(unique(wtmp$Dive.Site)," dive site"))
pl
# pdf(paste0("figs/dive_",unique(wtmp$Dive.Site),"_mds.pdf"), width=14,height = 14)
ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_mds.pdf"),
       width = 14, height = 14, units="in",plot=pl)
rm(pl,mds_scores,centroid_Site,spp_scores)

## analyse for sig diffs between years ####
(m1 <- adonis2(wtmpord~wtmp$Year))
mv_wtmpord <- mvabund::as.mvabund(wtmpord)
mvabund::meanvar.plot(mv_wtmpord)
m2 <- manylm(mv_wtmpord ~ wtmp$Year);summary(m2)
m2out <- mvabund::anova.manylm(m2,p.uni = "adjusted")
m2tmp1 <- t(as.data.frame(m2out$uni.p))[,2]
names(m2tmp1) <- names(wtmpord)
m2tmp1[m2tmp1<0.056]

m2tx <- names(m2tmp1[m2tmp1<0.056])#which taxa are 'significantly' different?
kptx <- names(wtmpord) %in% m2tx
m2tx <- wtmpord[, kptx]
m2tx$Year <- wtmp$Year

##make long
m2txl <- m2tx %>% 
  relocate(Year) %>% #move Year to start
  pivot_longer(cols = c(2:ncol(.)))

ggplot(m2txl)+
  geom_boxplot(aes(x=name,y=value),varwidth = TRUE)+
  facet_wrap(.~Year)+
  coord_flip()+
  labs(y="Relative abundance",
       caption=paste0(unique(wtmp$Dive.Site)," dive site\nAbundances for each year on the same relative scale"))+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl2

ggplot(m2txl,aes(x=value,y=name, fill=Year, shape = Year))+
  geom_jitter(height = 0.15,size=3, alpha = 0.4) +
  scale_shape_manual(values = c(21:24))+
  labs(x="Relative abundance",
       caption=paste0(unique(wtmp$Dive.Site)," dive site\nAbundances for each year on the same relative scale"))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl3

ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_relabund.pdf"),
       width = 14, height = 6, units="in",plot=pl2)
ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_relabund_ver2.pdf"),
       width = 14, height = 6, units="in",plot=pl3)
rm(wtmp,wtmpord,m1,mv_wtmpord,m2,m2out,m2tx,m2txl,pl2,pl3)

##########
## ENDG Buoy ####
wtmp <- dfwDvENDGB

## remove 'empty' columns
wtmp <- wtmp %>% 
  # dplyr::select(where(is.character)|across(wher(is.numeric),-sum(.) !=0))
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

# create data for ordination
wtmp %>% 
  dplyr::select(-c(1:12)) -> wtmpord

## how many unique taxa?
print(paste0(ncol(wtmpord)," unique taxa"))

## by year(?)
vegan::specnumber(wtmpord, groups=wtmp$Year)

## run ordination ##
set.seed(80);tmpord <- metaMDS(wtmpord, try = 30, trymax = 200)
plot(tmpord)

### ordination plot using ggplot ####
### extract site scores and groups
mds_scores <- as_tibble(as.data.frame(scores(tmpord,"site")))
mds_scores$Year <- wtmp$Year
mds_scores$Station.Code <- wtmp$Station.Code
mds_scores$Sample.Ref <- wtmp$Sample.Ref
mds_scores$`Depth.(m)` <- wtmp$`Depth.(m)`
mds_scores$Quality <- wtmp$Quality
mds_scores$Kelp <- wtmp$Kelp
mds_scores$MNCR.Code <- wtmp$MNCR.Code
mds_scores$BSH <- wtmp$BSH

### extract species scores and groups
spp_scores <- as.data.frame(scores(tmpord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)
head(spp_scores)  #look at the data

### create centroids by Year
centroid_Site <- mds_scores %>% 
  group_by(Year) %>%
  summarise(NMDS1=mean(NMDS1),NMDS2=mean(NMDS2))## 2d version

### create link to Year centroids
mds_scores <- mds_scores %>% group_by(Year) %>% 
  mutate(BSH_centroid1=mean(NMDS1),
         BSH_centroid2=mean(NMDS2)) %>% ungroup ##2d version

### generate plot ####
pl <- ggplot(mds_scores,
             aes(x=NMDS1,y=NMDS2,colour=Year,shape=Year,
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
  geom_point(show.legend = FALSE,colour=1,aes(fill=Year),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_Site,aes(x=NMDS1,y=NMDS2,label=Year,fill=Year),
    width = unit(0.05, "npc"),
    inherit.aes = FALSE,show.legend = FALSE,
    hjust=0.5,vjust=0.5,halign=0.5,valign=0.5,
    fontface="bold")+
  scale_shape_manual(values = c(21:24))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  coord_fixed()+
  labs(x="nmMDS Axis 1",
       y="nmMDS Axis 2")+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(tmpord$stress, 3))))+
  labs(caption=paste0(unique(wtmp$Dive.Site)," dive site"))
pl
# pdf(paste0("figs/dive_",unique(wtmp$Dive.Site),"_mds.pdf"), width=14,height = 14)
ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_mds.pdf"),
       width = 14, height = 14, units="in",plot=pl)
rm(pl,mds_scores,centroid_Site,spp_scores)

## analyse for sig diffs between years ####
(m1 <- adonis2(wtmpord~wtmp$Year))
mv_wtmpord <- mvabund::as.mvabund(wtmpord)
mvabund::meanvar.plot(mv_wtmpord)
m2 <- manylm(mv_wtmpord ~ wtmp$Year);summary(m2)
m2out <- mvabund::anova.manylm(m2,p.uni = "adjusted")
m2tmp1 <- t(as.data.frame(m2out$uni.p))[,2]
names(m2tmp1) <- names(wtmpord)
m2tmp1[m2tmp1<0.056]
range(m2tmp1[m2tmp1<0.051])

m2tx <- names(m2tmp1[m2tmp1<0.051])#which taxa are 'significantly' different?
kptx <- names(wtmpord) %in% m2tx
m2tx <- wtmpord[, kptx]
m2tx$Year <- wtmp$Year

##make long
m2txl <- m2tx %>% 
  relocate(Year) %>% #move Year to start
  pivot_longer(cols = c(2:ncol(.)))

ggplot(m2txl)+
  geom_boxplot(aes(x=name,y=value),varwidth = TRUE)+
  facet_wrap(.~Year)+
  coord_flip()+
  labs(y="Relative abundance",
       caption=paste0(unique(wtmp$Dive.Site)," dive site\nAbundances for each year on the same relative scale"))+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl2

ggplot(m2txl,aes(x=value,y=name, fill=Year, shape = Year))+
  geom_jitter(height = 0.15,size=3, alpha = 0.4) +
  scale_shape_manual(values = c(21:24))+
  labs(x="Relative abundance",
       caption=paste0(unique(wtmp$Dive.Site)," dive site\nAbundances for each year on the same relative scale"))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl3

ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_relabund.pdf"),
       width = 14, height = 6, units="in",plot=pl2)
ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_relabund_ver2.pdf"),
       width = 14, height = 6, units="in",plot=pl3)
rm(wtmp,wtmpord,m1,mv_wtmpord,m2,m2out,m2tx,m2txl,pl2,pl3)

##########
## Eastern Kings ####
wtmp <- dfwDvEKng

## remove 'empty' columns
wtmp <- wtmp %>% 
  # dplyr::select(where(is.character)|across(wher(is.numeric),-sum(.) !=0))
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

# create data for ordination
wtmp %>% 
  dplyr::select(-c(1:12)) -> wtmpord

## how many unique taxa?
print(paste0(ncol(wtmpord)," unique taxa"))

## by year(?)
vegan::specnumber(wtmpord, groups=wtmp$Year)

## run ordination ##
set.seed(80);tmpord <- metaMDS(wtmpord, try = 30, trymax = 200)
plot(tmpord)

### ordination plot using ggplot ####
### extract site scores and groups
mds_scores <- as_tibble(as.data.frame(scores(tmpord,"site")))
mds_scores$Year <- wtmp$Year
mds_scores$Station.Code <- wtmp$Station.Code
mds_scores$Sample.Ref <- wtmp$Sample.Ref
mds_scores$`Depth.(m)` <- wtmp$`Depth.(m)`
mds_scores$Quality <- wtmp$Quality
mds_scores$Kelp <- wtmp$Kelp
mds_scores$MNCR.Code <- wtmp$MNCR.Code
mds_scores$BSH <- wtmp$BSH

### extract species scores and groups
spp_scores <- as.data.frame(scores(tmpord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)
head(spp_scores)  #look at the data

### create centroids by Year
centroid_Site <- mds_scores %>% 
  group_by(Year) %>%
  summarise(NMDS1=mean(NMDS1),NMDS2=mean(NMDS2))## 2d version

### create link to Year centroids
mds_scores <- mds_scores %>% group_by(Year) %>% 
  mutate(BSH_centroid1=mean(NMDS1),
         BSH_centroid2=mean(NMDS2)) %>% ungroup ##2d version

### generate plot ####
pl <- ggplot(mds_scores,
             aes(x=NMDS1,y=NMDS2,colour=Year,shape=Year,
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
  geom_point(show.legend = FALSE,colour=1,aes(fill=Year),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_Site,aes(x=NMDS1,y=NMDS2,label=Year,fill=Year),
    width = unit(0.05, "npc"),
    inherit.aes = FALSE,show.legend = FALSE,
    hjust=0.5,vjust=0.5,halign=0.5,valign=0.5,
    fontface="bold")+
  scale_shape_manual(values = c(21:24))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  coord_fixed()+
  labs(x="nmMDS Axis 1",
       y="nmMDS Axis 2")+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(tmpord$stress, 3))))+
  labs(caption=paste0(unique(wtmp$Dive.Site)," dive site"))
pl
# pdf(paste0("figs/dive_",unique(wtmp$Dive.Site),"_mds.pdf"), width=14,height = 14)
ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_mds.pdf"),
       width = 14, height = 14, units="in",plot=pl)
rm(pl,mds_scores,centroid_Site,spp_scores)

## analyse for sig diffs between years ####
(m1 <- adonis2(wtmpord~wtmp$Year))
mv_wtmpord <- mvabund::as.mvabund(wtmpord)
mvabund::meanvar.plot(mv_wtmpord)
m2 <- manylm(mv_wtmpord ~ wtmp$Year);summary(m2)
m2out <- mvabund::anova.manylm(m2,p.uni = "adjusted")
m2tmp1 <- t(as.data.frame(m2out$uni.p))[,2]
names(m2tmp1) <- names(wtmpord)
m2tmp1[m2tmp1<0.056]
range(m2tmp1[m2tmp1<0.051])

m2tx <- names(m2tmp1[m2tmp1<0.051])#which taxa are 'significantly' different?
kptx <- names(wtmpord) %in% m2tx
m2tx <- wtmpord[, kptx]
m2tx$Year <- wtmp$Year

##make long
m2txl <- m2tx %>% 
  relocate(Year) %>% #move Year to start
  pivot_longer(cols = c(2:ncol(.)))

ggplot(m2txl)+
  geom_boxplot(aes(x=name,y=value),varwidth = TRUE)+
  facet_wrap(.~Year)+
  coord_flip()+
  labs(y="Relative abundance",
       caption=paste0(unique(wtmp$Dive.Site)," dive site\nAbundances for each year on the same relative scale"))+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl2

ggplot(m2txl,aes(x=value,y=name, fill=Year, shape = Year))+
  geom_jitter(height = 0.15,size=3, alpha = 0.4) +
  scale_shape_manual(values = c(21:24))+
  labs(x="Relative abundance",
       caption=paste0(unique(wtmp$Dive.Site)," dive site\nAbundances for each year on the same relative scale"))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl3

ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_relabund.pdf"),
       width = 14, height = 6, units="in",plot=pl2)
ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_relabund_ver2.pdf"),
       width = 14, height = 6, units="in",plot=pl3)
rm(wtmp,wtmpord,m1,mv_wtmpord,m2,m2out,m2tx,m2txl,pl2,pl3)
##########

## Mew Stone West####
wtmp <- dfwDvMSW

## remove 'empty' columns
wtmp <- wtmp %>% 
  # dplyr::select(where(is.character)|across(wher(is.numeric),-sum(.) !=0))
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

# create data for ordination
wtmp %>% 
  dplyr::select(-c(1:12)) -> wtmpord

## how many unique taxa?
print(paste0(ncol(wtmpord)," unique taxa"))

## by year(?)
vegan::specnumber(wtmpord, groups=wtmp$Year)

## run ordination ##
set.seed(80);tmpord <- metaMDS(wtmpord, try = 30, trymax = 200)
plot(tmpord)

### ordination plot using ggplot ####
### extract site scores and groups
mds_scores <- as_tibble(as.data.frame(scores(tmpord,"site")))
mds_scores$Year <- wtmp$Year
mds_scores$Station.Code <- wtmp$Station.Code
mds_scores$Sample.Ref <- wtmp$Sample.Ref
mds_scores$`Depth.(m)` <- wtmp$`Depth.(m)`
mds_scores$Quality <- wtmp$Quality
mds_scores$Kelp <- wtmp$Kelp
mds_scores$MNCR.Code <- wtmp$MNCR.Code
mds_scores$BSH <- wtmp$BSH

### extract species scores and groups
spp_scores <- as.data.frame(scores(tmpord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)
head(spp_scores)  #look at the data

### create centroids by Year
centroid_Site <- mds_scores %>% 
  group_by(Year) %>%
  summarise(NMDS1=mean(NMDS1),NMDS2=mean(NMDS2))## 2d version

### create link to Year centroids
mds_scores <- mds_scores %>% group_by(Year) %>% 
  mutate(BSH_centroid1=mean(NMDS1),
         BSH_centroid2=mean(NMDS2)) %>% ungroup ##2d version

### generate plot ####
pl <- ggplot(mds_scores,
             aes(x=NMDS1,y=NMDS2,colour=Year,shape=Year,
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
  geom_point(show.legend = FALSE,colour=1,aes(fill=Year),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_Site,aes(x=NMDS1,y=NMDS2,label=Year,fill=Year),
    width = unit(0.05, "npc"),
    inherit.aes = FALSE,show.legend = FALSE,
    hjust=0.5,vjust=0.5,halign=0.5,valign=0.5,
    fontface="bold")+
  scale_shape_manual(values = c(21:24))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  coord_fixed()+
  labs(x="nmMDS Axis 1",
       y="nmMDS Axis 2")+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(tmpord$stress, 3))))+
  labs(caption=paste0(unique(wtmp$Dive.Site)," dive site"))
pl
# pdf(paste0("figs/dive_",unique(wtmp$Dive.Site),"_mds.pdf"), width=14,height = 14)
ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_mds.pdf"),
       width = 14, height = 14, units="in",plot=pl)
rm(pl,mds_scores,centroid_Site,spp_scores)

## analyse for sig diffs between years ####
(m1 <- adonis2(wtmpord~wtmp$Year))
mv_wtmpord <- mvabund::as.mvabund(wtmpord)
mvabund::meanvar.plot(mv_wtmpord)
m2 <- manylm(mv_wtmpord ~ wtmp$Year);summary(m2)
m2out <- mvabund::anova.manylm(m2,p.uni = "adjusted")
m2tmp1 <- t(as.data.frame(m2out$uni.p))[,2]
names(m2tmp1) <- names(wtmpord)
m2tmp1[m2tmp1<0.056]
range(m2tmp1[m2tmp1<0.051])

m2tx <- names(m2tmp1[m2tmp1<0.051])#which taxa are 'significantly' different?
kptx <- names(wtmpord) %in% m2tx
m2tx <- wtmpord[, kptx]
m2tx$Year <- wtmp$Year

##make long
m2txl <- m2tx %>% 
  relocate(Year) %>% #move Year to start
  pivot_longer(cols = c(2:ncol(.)))

ggplot(m2txl)+
  geom_boxplot(aes(x=name,y=value),varwidth = TRUE)+
  facet_wrap(.~Year)+
  coord_flip()+
  labs(y="Relative abundance",
       caption=paste0(unique(wtmp$Dive.Site)," dive site\nAbundances for each year on the same relative scale"))+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl2

ggplot(m2txl,aes(x=value,y=name, fill=Year, shape = Year))+
  geom_jitter(height = 0.15,size=3, alpha = 0.4) +
  scale_shape_manual(values = c(21:24))+
  labs(x="Relative abundance",
       caption=paste0(unique(wtmp$Dive.Site)," dive site\nAbundances for each year on the same relative scale"))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl3

ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_relabund.pdf"),
       width = 14, height = 6, units="in",plot=pl2)
ggsave(filename = paste0("figs/dive_",unique(wtmp$Dive.Site),"_relabund_ver2.pdf"),
       width = 14, height = 6, units="in",plot=pl3)
rm(wtmp,wtmpord,m1,mv_wtmpord,m2,m2out,m2tx,m2txl,pl2,pl3)

##########

