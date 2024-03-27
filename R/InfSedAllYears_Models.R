# InfSedAllYears_Models.R (file name) #### 
# Title: 

# Set up ####
# source("R/datfol.R")
source("R/metadata.R")

## load required packages ####
ld_pkgs <- c("tidyverse","vegan","lmerTest","rstatix","mvabund","gllvm",
             "MASS","ggtext","ggpmisc") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

# load data (long format) ####
dfw <- readRDS(paste0(datfol,"Data/2021 data/Infauna/Inf_Sed_all_yearsWIDE.Rdat"))
dfw$S <- vegan::specnumber(dfw[,-c(1:17)])
dfw$YrFactor <- as.factor(dfw$Year)
dfw$YrStandard <- dfw$Year-2011 #create standardised 'year' variable so that 2011=0
#(therefore intercept represents 2011, rather than the year 0)

#reorder columns
dfw <- dfw %>% 
  relocate(YrFactor, .after = Year) %>% 
  relocate(YrStandard, .after = YrFactor)

# initial explorations ####
hist(dfw$S, breaks = 20); range(dfw$S)

ggplot(dfw[!is.na(dfw$BSH_CODE),], aes(x = Year, y = S,group=BSH_CODE))+
  geom_jitter(width=0.2)+
  geom_smooth()+
  facet_wrap(.~BSH_CODE)

ggplot(dfw[!is.na(dfw$BSH_CODE) & dfw$Year == 2021,], aes(x = BSH_CODE, y = S,group=BSH_CODE))+
  geom_jitter(width=0.2) +
  ggtitle("2021 taxon richness by BSH")

# model taxon richness vs BSH ####
summary(m01 <- lm(S ~ BSH_CODE, data = dfw[dfw$Year == 2021,]))
anova(m01)
sjPlot::plot_model(m01)
report::report(m01)

# models of tax richness over time ####
## basic model
summary(m1 <- lm(S ~ YrStandard, data = dfw))
table(dfw$Year, dfw$BSH_CODE)
performance::check_model(m1)
summary(mixedmod1 <- lmer(S ~ YrStandard + (1|BSH_CODE), data = dfw))
sjPlot::plot_model(m1)

## model with interaction
summary(m2 <- lm(S ~ YrStandard*BSH_CODE, data = dfw))
sjPlot::plot_model(m2)
performance::check_model(m2)
report::report(m2)

#### run as Poisson as it';'s count data
summary(m_pois <- glm(S ~ YrStandard*BSH_CODE, family = poisson(),data=dfw))
sjPlot::plot_model(m_pois)
performance::check_model(m_pois)

#### is negative binomial better?
summary(m_nb <- glm.nb(S ~ YrStandard*BSH_CODE, data=dfw))
sjPlot::plot_model(m_nb)
performance::check_model(m_nb)
report::report(m_nb)

AIC(m1,m2,m_pois,m_nb) #lowest AIC reflects 'best' fitting model

# select most recent data and run ordination ####
# most recent year is 2021
dfwcur <- dfw %>% 
  filter(.,Year == 2021) %>% ## keep current year
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

### look into presence-absence transform AND/OR account for -999 values

tmp <- dfwcur %>% 
  dplyr::select(-c(S,Site_Yr_mesh:BSH_CODE))

tmp <- tmp %>% 
  mutate_all(~ ifelse(. < 0, 1, .))

set.seed(pi);ord <- vegan::metaMDS(tmp, trymax = 500)
plot(ord)

### run unconstrained ordination using gllvm. See: https://jenniniku.github.io/gllvm/articles/vignette3.html
ptm <- Sys.time() # start timer
# xps <- gllvm::gllvm(tmp,family = poisson(), num.lv = 2)
# xnb <- gllvm::gllvm(tmp,family = "negative.binomial", num.lv = 2)
# saveRDS(xnb, file = "outputs/gllvm.inf.2021.rdat")
xnb <- readRDS("outputs/gllvm.inf.2021.rdat")
Sys.time() - ptm; rm(ptm)

ordiplot.gllvm(xnb, biplot = TRUE)

## make neater version using ggplot ####
### extract site scores and groups
mds_scores <- as_tibble(as.data.frame(scores(ord,"site")))
mds_scores$BSH <- dfwcur$BSH_CODE
mds_scores$station <- dfwcur$Station.Code
mds_scores$code <- dfwcur$Site_Yr_mesh

### extract species scores and groups
spp_scores <- as.data.frame(scores(ord, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
spp_scores$species <- rownames(spp_scores)  # create a column of species, from the rownames of species.scores
spp_scores$species_sh <- make.cepnames(spp_scores$species)
head(spp_scores)  #look at the data

### create centroids by BSH
centroid_BSH <- mds_scores %>% 
  group_by(BSH) %>%
  summarise(NMDS1=mean(NMDS1),NMDS2=mean(NMDS2))## 2d version

### create link to BSH centroids
mds_scores <- mds_scores %>% group_by(BSH) %>% 
  mutate(BSH_centroid1=mean(NMDS1),
         BSH_centroid2=mean(NMDS2)) %>% ungroup ##2d version

## append taxonomy to species data ####
sp_tx <- openxlsx::read.xlsx(paste0(datfol_Rproj,"data_in/plyinftax_matched.xlsx"),
                             sheet = "tx")
# sp_tx$species <- sp_tx$ScientificName_accepted
# sp_tx <- sp_tx %>% 
#   dplyr::select(.,-c(ScientificName.Concat:Authority_accepted))

spp_scores <- left_join(spp_scores,sp_tx, by="species")

pl <- ggplot(mds_scores,
       aes(x=NMDS1,y=NMDS2,colour=BSH,shape=BSH,
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
  geom_point(show.legend = FALSE,colour=1,aes(fill=BSH),
             size = 2.5)+
  # stat_ellipse(geom="polygon",
  #              level = 0.9,aes(fill=BSH),
  #              alpha=0.25,show.legend = FALSE)+
  geom_textbox(#size=6,
               data=centroid_BSH,aes(x=NMDS1,y=NMDS2,label=BSH,fill=BSH),
               width = unit(0.05, "npc"),
               inherit.aes = FALSE,show.legend = FALSE,hjust=0.5,
               fontface="bold")+
  scale_shape_manual(values = c(21:24))+
  scale_fill_manual(name=NULL,
                    breaks = c("A5.1","A5.2","A5.3","A5.4"),
                    labels = c("A5.1","A5.2","A5.3","A5.4"),
                    values = cbPalette)+
  scale_colour_manual(name=NULL,
                      breaks = c("A5.1","A5.2","A5.3","A5.4"),
                      labels = c("A5.1","A5.2","A5.3","A5.4"),
                      values = cbPalette)+
  coord_fixed()+
  labs(x="nmMDS Axis 1",
       y="nmMDS Axis 2")+
  geom_text_npc(aes(npcx = .99, npcy = .99, label=paste("Stress = ",
                                                        round(ord$stress, 3))))

pdf("figs/inf_mds.pdf", width=14,height = 14)
print(pl)
dev.off()
########

# compare temporal trends within each BSH ####
dfwtmp <- dfw %>% 
  filter(!is.na(BSH_CODE)) %>% 
  dplyr::select(-S)
dfwA51 <- dfwtmp[dfwtmp$BSH_CODE == "A5.1",]
dfwA52 <- dfwtmp[dfwtmp$BSH_CODE == "A5.2",]
dfwA53 <- dfwtmp[dfwtmp$BSH_CODE == "A5.3",]
dfwA54 <- dfwtmp[dfwtmp$BSH_CODE == "A5.4",]

## A5.1 ####
wtmp <- dfwA51

## remove 'empty' columns
wtmp <- wtmp %>%
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

# create data for ordination
wtmp %>% 
  dplyr::select(-c(1:19)) -> wtmpord

wtmpord <- wtmpord %>% # convert -999 to 1
  mutate_all(~ ifelse(. < 0, 1, .))

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
mds_scores$Easting <- wtmp$Easting
mds_scores$Northing <- wtmp$Northing
mds_scores$Station.Code <- wtmp$Station.Code
mds_scores$GEAR <- wtmp$GEAR
mds_scores$GRAVELPCT <- wtmp$GRAVELPCT
mds_scores$SANDPCT <- wtmp$SANDPCT
mds_scores$MUDPCT <- wtmp$MUDPCT
mds_scores$FOLK <- wtmp$FOLK
mds_scores$BSH <- wtmp$BSH
mds_scores$BSH_CODE <- wtmp$BSH_CODE

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
             aes(x=NMDS1,y=NMDS2,colour=as.factor(Year),shape=as.factor(Year),
                 xend=BSH_centroid1,
                 yend=BSH_centroid2))+  
  geom_hline(yintercept = 0,lty=2, col="grey")+
  geom_vline(xintercept = 0,lty=2, col="grey")+
  geom_text(data=spp_scores, aes(x=NMDS1, y= NMDS2,label=species_sh),
            col="grey",
            alpha=0.5,
            size=6,
            inherit.aes = FALSE)+
  geom_segment(show.legend = FALSE,alpha=0.6)+## add 'spider legs' joining samples-centroids
  geom_point(show.legend = FALSE,colour=1,aes(fill=as.factor(Year)),
             size = 2.5)+
  geom_textbox(size=5,
    data=centroid_Site,aes(x=NMDS1,y=NMDS2,label=Year,fill=as.factor(Year)),
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
  labs(caption=paste0(unique(wtmp$BSH)[1]," infaunal data"))+
  theme(plot.caption = element_text(size=16,face="bold"),
        axis.title = element_text(size=14, face="bold"))
pl

ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_mds.pdf"),
       width = 14, height = 14, units="in",plot=pl)
rm(pl,mds_scores,centroid_Site,spp_scores)

## analyse for sig diffs between years ####
mv_wtmpord <- mvabund::as.mvabund(wtmpord)

##################
##mean-variance plot ####
png(file = "figs/infMeanVarA5.1.png",
    width=12*ppi, height=6*ppi, res=ppi)

mvpl <- mvabund::meanvar.plot(mv_wtmpord,
                              xlab="Mean",
                              ylab="Variance",
                              table=TRUE)

# Step 1: Find the minimum and maximum values
min_value <- min(mvpl[,2])
max_value <- max(mvpl[,2])

min_order <- floor(log10(min_value))
max_order <- floor(log10(max_value))
orders_of_magnitude_covered <- max_order - min_order

ttl <- "Very strong mean-variance relationship in invertebrate abundances in the A5.1 Subtidal coarse sediment broadscale habitat"
sbtt <- paste0("Variance within the dataset covers *",orders_of_magnitude_covered," orders of magnitude*.")

mtext(side=3, line = 1, at =-0.07, adj=0, cex = 1, ttl, font=2)
mtext(side=3, line = 0.25, at =-0.07, adj=0, cex = 0.7, sbtt)

dev.off()

#####################

# m2 <- manyglm(mv_wtmpord ~ wtmp$Year, family = "negative.binomial");summary(m2)
# saveRDS(m2, file = paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".rdat"))
m2 <- readRDS(paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".rdat"))
# m2out <- mvabund::anova.manyglm(m2,p.uni = "adjusted")
# saveRDS(m2out, file = paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".pw.rdat"))
m2out <- readRDS(paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".pw.rdat"))
m2tmp1 <- t(as.data.frame(m2out$uni.p))[,2]
names(m2tmp1) <- names(wtmpord)
m2tmp1[m2tmp1<0.056]; range(m2tmp1[m2tmp1<0.051])
print(paste0(length(names(m2tmp1[m2tmp1<0.056]))," 'significant' taxa"))

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
  labs(y="Taxon abundance",
       caption=paste0(unique(wtmp$BSH)[1]," BSH"))+
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
  labs(x="log(Taxon abundance (n+1))",
       caption=paste0(unique(wtmp$BSH)[1]," BSH"))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
    facet_wrap(.~Year)+
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
        axis.text.y = element_text(size=12,face="italic"),
        # axis.text.x = element_blank(),
        axis.title.x = element_text(face="bold"),
        strip.text.x = element_text(face="bold",size=12),
        plot.caption = element_text(face="bold",size=12)) -> pl3)

ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_relabund.pdf"),
       width = 14, height = 6, units="in",plot=pl2)
ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_relabund_ver2.pdf"),
       width = 14, height = 6, units="in",plot=pl3)
rm(wtmp,wtmpord,mv_wtmpord,m2,m2out,m2tx,m2txl,pl2,pl3)
#######
## A5.2 ####
wtmp <- dfwA52

## remove 'empty' columns
wtmp <- wtmp %>%
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

# create data for ordination
wtmp %>% 
  dplyr::select(-c(1:19)) -> wtmpord

wtmpord <- wtmpord %>% # convert -999 to 1
  mutate_all(~ ifelse(. < 0, 1, .))

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
mds_scores$Easting <- wtmp$Easting
mds_scores$Northing <- wtmp$Northing
mds_scores$Station.Code <- wtmp$Station.Code
mds_scores$GEAR <- wtmp$GEAR
mds_scores$GRAVELPCT <- wtmp$GRAVELPCT
mds_scores$SANDPCT <- wtmp$SANDPCT
mds_scores$MUDPCT <- wtmp$MUDPCT
mds_scores$FOLK <- wtmp$FOLK
mds_scores$BSH <- wtmp$BSH
mds_scores$BSH_CODE <- wtmp$BSH_CODE

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
             aes(x=NMDS1,y=NMDS2,colour=as.factor(Year),shape=as.factor(Year),
                 xend=BSH_centroid1,
                 yend=BSH_centroid2))+  
  geom_hline(yintercept = 0,lty=2, col="grey")+
  geom_vline(xintercept = 0,lty=2, col="grey")+
  geom_text(data=spp_scores, aes(x=NMDS1, y= NMDS2,label=species_sh),
            col="grey",
            alpha=0.5,
            inherit.aes = FALSE)+
  geom_segment(show.legend = FALSE,alpha=0.6)+## add 'spider legs' joining samples-centroids
  geom_point(show.legend = FALSE,colour=1,aes(fill=as.factor(Year)),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_Site,aes(x=NMDS1,y=NMDS2,label=Year,fill=as.factor(Year)),
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
  labs(caption=paste0(unique(wtmp$BSH)[1]," infaunal data"))
pl

ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_mds.pdf"),
       width = 14, height = 14, units="in",plot=pl)
rm(pl,mds_scores,centroid_Site,spp_scores)

## analyse for sig diffs between years ####
mv_wtmpord <- mvabund::as.mvabund(wtmpord)
mvabund::meanvar.plot(mv_wtmpord)
# m2 <- manyglm(mv_wtmpord ~ wtmp$Year, family = "negative.binomial");summary(m2)
# saveRDS(m2, file = paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".rdat"))
m2 <- readRDS(paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".rdat"))
# m2out <- mvabund::anova.manyglm(m2,p.uni = "adjusted")
# saveRDS(m2out, file = paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".pw.rdat"))
m2out <- readRDS(paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".pw.rdat"))
m2tmp1 <- t(as.data.frame(m2out$uni.p))[,2]
names(m2tmp1) <- names(wtmpord)
m2tmp1[m2tmp1<0.056]; range(m2tmp1[m2tmp1<0.051])
print(paste0(length(names(m2tmp1[m2tmp1<0.056]))," 'significant' taxa"))

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
  labs(y="Taxon abundance",
       caption=paste0(unique(wtmp$BSH)[1]," BSH"))+
  theme(axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) ->pl2

m2txl$Year <- as.factor(m2txl$Year)
ggplot(m2txl,aes(x=log(value+1),y=name, fill=Year, shape = Year))+
  geom_jitter(height = 0.15,size=3, alpha = 0.4) +
  scale_shape_manual(values = c(21:24))+
  labs(x="log(Taxon abundance (n+1))",
       caption=paste0(unique(wtmp$BSH)[1]," BSH"))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  theme(axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl3

ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_relabund.pdf"),
       width = 14, height = 6, units="in",plot=pl2)
ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_relabund_ver2.pdf"),
       width = 14, height = 6, units="in",plot=pl3)
rm(wtmp,wtmpord,mv_wtmpord,m2,m2out,m2tx,m2txl,pl2,pl3)
#######

## A5.3 ####
wtmp <- dfwA53

## remove 'empty' columns
wtmp <- wtmp %>%
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

# create data for ordination
wtmp %>% 
  dplyr::select(-c(1:19)) -> wtmpord

wtmpord <- wtmpord %>% # convert -999 to 1
  mutate_all(~ ifelse(. < 0, 1, .))

## how many unique taxa?
print(paste0(ncol(wtmpord)," unique taxa"))

## by year(?)
vegan::specnumber(wtmpord, groups=wtmp$Year)

## remove "Protodorvillea kefersteini" (skews ordination)
wtmpord <- wtmpord %>% 
  dplyr::select(-c("Protodorvillea kefersteini"))

kp <- vegan::specnumber(wtmpord) >0
wtmpord <- wtmpord[kp,]

## run ordination ##
set.seed(80);tmpord <- metaMDS(wtmpord, try = 30, trymax = 200)
plot(tmpord)

wtmp <- wtmp[kp,]

### ordination plot using ggplot ####
### extract site scores and groups
mds_scores <- as_tibble(as.data.frame(scores(tmpord,"site")))
mds_scores$Year <- wtmp$Year
mds_scores$Easting <- wtmp$Easting
mds_scores$Northing <- wtmp$Northing
mds_scores$Station.Code <- wtmp$Station.Code
mds_scores$GEAR <- wtmp$GEAR
mds_scores$GRAVELPCT <- wtmp$GRAVELPCT
mds_scores$SANDPCT <- wtmp$SANDPCT
mds_scores$MUDPCT <- wtmp$MUDPCT
mds_scores$FOLK <- wtmp$FOLK
mds_scores$BSH <- wtmp$BSH
mds_scores$BSH_CODE <- wtmp$BSH_CODE

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
             aes(x=NMDS1,y=NMDS2,colour=as.factor(Year),shape=as.factor(Year),
                 xend=BSH_centroid1,
                 yend=BSH_centroid2))+  
  geom_hline(yintercept = 0,lty=2, col="grey")+
  geom_vline(xintercept = 0,lty=2, col="grey")+
  geom_text(data=spp_scores, aes(x=NMDS1, y= NMDS2,label=species_sh),
            col="grey",
            alpha=0.5,
            inherit.aes = FALSE)+
  geom_segment(show.legend = FALSE,alpha=0.6)+## add 'spider legs' joining samples-centroids
  geom_point(show.legend = FALSE,colour=1,aes(fill=as.factor(Year)),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_Site,aes(x=NMDS1,y=NMDS2,label=Year,fill=as.factor(Year)),
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
  labs(caption=paste0(unique(wtmp$BSH)[1]," infaunal data\nOne site containing only a single taxon from 2011 removed"))
pl

ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_mds.pdf"),
       width = 14, height = 14, units="in",plot=pl)
rm(pl,mds_scores,centroid_Site,spp_scores)

## analyse for sig diffs between years ####
mv_wtmpord <- mvabund::as.mvabund(wtmpord)
mvabund::meanvar.plot(mv_wtmpord)
# m2 <- manyglm(mv_wtmpord ~ wtmp$Year, family = "negative.binomial");summary(m2)
# saveRDS(m2, file = paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".rdat"))
m2 <- readRDS(paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".rdat"))
# m2out <- mvabund::anova.manyglm(m2,p.uni = "adjusted")
# saveRDS(m2out, file = paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".pw.rdat"))
m2out <- readRDS(paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".pw.rdat"))
m2tmp1 <- t(as.data.frame(m2out$uni.p))[,2]
names(m2tmp1) <- names(wtmpord)
m2tmp1[m2tmp1<0.056]; range(m2tmp1[m2tmp1<0.051])
print(paste0(length(names(m2tmp1[m2tmp1<0.056]))," 'significant' taxa"))

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
  labs(y="Taxon abundance",
       caption=paste0(unique(wtmp$BSH)[1]," BSH"))+
  theme(axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) ->pl2

m2txl$Year <- as.factor(m2txl$Year)
ggplot(m2txl,aes(x=log(value+1),y=name, fill=Year, shape = Year))+
  geom_jitter(height = 0.15,size=3, alpha = 0.4) +
  scale_shape_manual(values = c(21:24))+
  labs(x="log(Taxon abundance (n+1))",
       caption=paste0(unique(wtmp$BSH)[1]," BSH"))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  theme(axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl3

ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_relabund.pdf"),
       width = 14, height = 6, units="in",plot=pl2)
ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_relabund_ver2.pdf"),
       width = 14, height = 6, units="in",plot=pl3)
rm(wtmp,wtmpord,mv_wtmpord,m2,m2out,m2tx,m2txl,pl2,pl3)
#######
## A5.4 ####
wtmp <- dfwA54

## remove 'empty' columns
wtmp <- wtmp %>%
  dplyr::select_if(where( ~ !is.numeric(.) || sum(.) !=0)) ## remove numeric variables that sum to 0

# create data for ordination
wtmp %>% 
  dplyr::select(-c(1:19)) -> wtmpord

wtmpord <- wtmpord %>% # convert -999 to 1
  mutate_all(~ ifelse(. < 0, 1, .))

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
mds_scores$Easting <- wtmp$Easting
mds_scores$Northing <- wtmp$Northing
mds_scores$Station.Code <- wtmp$Station.Code
mds_scores$GEAR <- wtmp$GEAR
mds_scores$GRAVELPCT <- wtmp$GRAVELPCT
mds_scores$SANDPCT <- wtmp$SANDPCT
mds_scores$MUDPCT <- wtmp$MUDPCT
mds_scores$FOLK <- wtmp$FOLK
mds_scores$BSH <- wtmp$BSH
mds_scores$BSH_CODE <- wtmp$BSH_CODE

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
             aes(x=NMDS1,y=NMDS2,colour=as.factor(Year),shape=as.factor(Year),
                 xend=BSH_centroid1,
                 yend=BSH_centroid2))+  
  geom_hline(yintercept = 0,lty=2, col="grey")+
  geom_vline(xintercept = 0,lty=2, col="grey")+
  geom_text(data=spp_scores, aes(x=NMDS1, y= NMDS2,label=species_sh),
            col="grey",
            alpha=0.5,
            inherit.aes = FALSE)+
  geom_segment(show.legend = FALSE,alpha=0.6)+## add 'spider legs' joining samples-centroids
  geom_point(show.legend = FALSE,colour=1,aes(fill=as.factor(Year)),
             size = 2.5)+
  geom_textbox(#size=6,
    data=centroid_Site,aes(x=NMDS1,y=NMDS2,label=Year,fill=as.factor(Year)),
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
  labs(caption=paste0(unique(wtmp$BSH)[1]," infaunal data"))
pl

ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_mds.pdf"),
       width = 14, height = 14, units="in",plot=pl)
rm(pl,mds_scores,centroid_Site,spp_scores)

## analyse for sig diffs between years ####
mv_wtmpord <- mvabund::as.mvabund(wtmpord)
mvabund::meanvar.plot(mv_wtmpord)
# m2 <- manyglm(mv_wtmpord ~ wtmp$Year, family = "negative.binomial");summary(m2)
# saveRDS(m2, file = paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".rdat"))
m2 <- readRDS(paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".rdat"))
# m2out <- mvabund::anova.manyglm(m2,p.uni = "adjusted")
# saveRDS(m2out, file = paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".pw.rdat"))
m2out <- readRDS(paste0("outputs/mvabund.inf.",unique(wtmp$BSH_CODE)[1],".pw.rdat"))
m2tmp1 <- t(as.data.frame(m2out$uni.p))[,2]
names(m2tmp1) <- names(wtmpord)
m2tmp1[m2tmp1<0.056]; range(m2tmp1[m2tmp1<0.051])
print(paste0(length(names(m2tmp1[m2tmp1<0.056]))," 'significant' taxa"))

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
  labs(y="Taxon abundance",
       caption=paste0(unique(wtmp$BSH)[1]," BSH"))+
  theme(axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) ->pl2

m2txl$Year <- as.factor(m2txl$Year)
ggplot(m2txl,aes(x=log(value+1),y=name, fill=Year, shape = Year))+
  geom_jitter(height = 0.15,size=3, alpha = 0.4) +
  scale_shape_manual(values = c(21:24))+
  labs(x="log(Taxon abundance (n+1))",
       caption=paste0(unique(wtmp$BSH)[1]," BSH"))+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  theme(axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        strip.text = element_text(face="bold")) -> pl3

ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_relabund.pdf"),
       width = 14, height = 14, units="in",plot=pl2)
ggsave(filename = paste0("figs/infauna_",unique(wtmp$BSH)[1],"_relabund_ver2.pdf"),
       width = 14, height = 14, units="in",plot=pl3)
rm(wtmp,wtmpord,mv_wtmpord,m2,m2out,m2tx,m2txl,pl2,pl3)
#######

# tidy up ####
detach("package:vegan",unload = TRUE)
detach("package:lmerTest",unload=TRUE)
detach("package:ggthemes",unload=TRUE)
detach("package:MASS",unload=TRUE)
detach("package:report",unload=TRUE)
detach("package:ggpubr",unload=TRUE)
detach("package:rstatix",unload=TRUE)
detach("package:ggtext",unload=TRUE)
detach("package:ggpmisc",unload=TRUE)
detach("package:mvabund",unload=TRUE)
detach("package:tidyverse",unload=TRUE)

rm(list = ls(pattern = "^df"))
rm(list = ls(pattern = "^m"))
rm(list = ls(pattern = "^tmp"))
rm(centroid_BSH,ord,sp_tx,cbPalette,datfol,kp,kptx,libfolder,nit,ppi,datfol_Rproj)
