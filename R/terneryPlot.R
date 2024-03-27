### terneryPlot.R ####
# produce ternery plot for sediment data ####

# Set up ####
# source("R/datfol.R")
source("R/metadata.R")

## load required packages ####
ld_pkgs <- c("tidyverse", "Ternary") # what packages do we need to load?
vapply(ld_pkgs, library, logical(1L), # load them and display TRUE/FALSE if loaded
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

### load data
df0 <- openxlsx::read.xlsx(paste0(datfol_Rproj,"data_in/2021_PSA_tern_plot.xlsx"),sheet="2021") #sample data

# dflns <- openxlsx::read.xlsx("data_in/2021_PSA_tern_plot.xlsx",sheet="lines") #gridlines

#toy example####
# x <- c(100,0,0)
# y <- c(0,100,0)
# z <- c(0,0,100)
# xpos <- 0.5 * (2*y + z) / (x+y+z)
# ypos <- sqrt(3) / 2 * z / (x+y+z)
# lb <- c("MUD","SAND","GRAVEL")
# df_dummy <- data.frame(x,y,z,xpos,ypos,lb)
# 
# ggplot(df_dummy,aes(xpos, ypos)) +
#   geom_point(aes(colour=lb,shape=lb),size=3) +
#   scale_colour_manual(values=cbPalette)+
#   annotate("path", x = c(0, 0.5, 1, 0), y = c(0,sqrt(3)/2,0,0)) +
#   coord_equal()+
#   xlim(-.10,1.1)+ylim(0,NA)+
#   # geom_segment(aes(x=.55, y=0.85, xend=1.05,yend=0),
#   #              arrow=arrow(length = unit(0.5, "cm")))+
#   theme(axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         panel.border = element_blank(),
#         legend.title = element_blank())


# manual version in ggplot2 ####
#### plot sediment data ####
### assign mud, sand and gravel %s to x, y, and z labels
x <- df0$mud
y <- df0$sand
z <- df0$gravel

## manually calculate locations in ternery space & append to data frame
xpos <- 0.5 * (2*y + z) / (x+y+z)
ypos <- sqrt(3) / 2 * z / (x+y+z)

### generate gridlines for plot
# xl <- dflns$Mud
# yl <- dflns$Snd
# zl <- dflns$Grv
# 
# xlnpos <- 0.5 * (2*yl + zl) / (xl+yl+zl)
# ylnpos <- sqrt(3) / 2 * zl / (xl+yl+zl)
# dflns <- as_tibble(data.frame(xlnpos,xlnpos))
df <- as_tibble(cbind(df0,xpos,ypos))

(pl <- ggplot(df,aes(xpos, ypos)) +
  # geom_point(aes(colour=BSH_CODE,shape=BSH_CODE),size=3) +
    geom_point(aes(colour=BSH,shape=BSH),size=3) +
    scale_colour_manual(values=cbPalette)+
    annotate("path", x = c(0, 0.5, 1, 0), y = c(0,sqrt(3)/2,0,0)) +coord_equal() +
    xlim(-.10,1.1) +
    geom_segment(aes(x=.55, y=0.85, xend=1.05,yend=0),#SAND
               arrow=arrow(length = unit(0.5, "cm")), linewidth = 1.5)+
    geom_segment(aes(x=-.05, y=0.005, xend=.45,yend=0.85),#GRAVEL
                 arrow=arrow(length = unit(0.5, "cm")), linewidth = 1.5) +
    geom_segment(aes(x=.95, y=-0.05, xend=.05,yend=-0.05),#MUD
                 arrow=arrow(length = unit(0.5, "cm")), linewidth = 1.5) +
    annotate(geom="text", x = 0.15, y=0.45, label = "Gravel", angle = 60, size=5)+
    annotate(geom="text", x = 0.85, y=0.45, label = "Sand", angle = 302, size=5)+
    annotate(geom="text", x = 0.5, y=-0.08, label = "Mud", angle = 0, size=5)+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.85,0.8))
  )
ggsave(plot=pl,filename = "figs/SedimentTernery.png")

# Using pkg::Ternary ####
# Ternary::TernaryApp()

par(mar = rep(0.3, 4))
# TernaryPlot("% Gravel", "% Sand", "% Mud")
# TernaryPoints(df0[,7:9], col = df0$BSH_CODE)
# TernaryPoints(df0[df0$BSH_CODE=="A5.1",][,7:9], col = 1, bg = 2, pch=21, cex=2)
# TernaryPoints(df0[df0$BSH_CODE=="A5.2",][,7:9], col = 1, bg = 3, pch=22, cex=2)
# TernaryPoints(df0[df0$BSH_CODE=="A5.3",][,7:9], col = 1, bg = 4, pch=23, cex=2)
# TernaryPoints(df0[df0$BSH_CODE=="A5.4",][,7:9], col = 1, bg = 5, pch=24, cex=2)
# legend("right", 
#        legend = c("A5.1", "A5.2", "A5.3", "A5.4"),
#        cex = 0.8, bty = "n", pch = 21:24, pt.cex = 1.8,
#        pt.bg = c(2:5))

#### ver 2 ####
pdf("figs/SedimentTernery.pdf", width=14,height = 14)
par(mar = rep(0.3, 4))
TernaryPlot(
  atip = "",
  btip = "",
  ctip = "",
  alab = "% Gravel",
  blab = "% Sand",
  clab = "% Silt/clay",
  #lab.offset = 0.16,
  lab.offset = 0.1,
  lab.col = "#000000",
  point = "up",
  clockwise = TRUE,
  xlim = NULL,
  ylim = NULL,
  # lab.cex = 1,
  # lab.cex = 1.5,
  lab.cex = 2.5, #axis label size
  lab.font = 2, #1 = normal, 2= bold
  tip.cex = 1,
  tip.font = 1,
  tip.col = "#000000",
  isometric = TRUE,
  atip.rotate = NULL,
  btip.rotate = NULL,
  ctip.rotate = NULL,
  atip.pos = NULL,
  btip.pos = NULL,
  ctip.pos = NULL,
  padding = 0.08,
  col = "#FFFFFF",
  grid.lines = 10,
  grid.col = "#A9A9A9",
  grid.lty = "solid",
  grid.lwd = 1,
  grid.minor.lines = 0,
  grid.minor.col = "#D3D3D3",
  grid.minor.lty = "solid",
  grid.minor.lwd = 1,
  axis.lty = "solid",
  axis.labels = TRUE,
  # axis.cex = 0.8,
  # axis.cex = 1.1,
  axis.cex = 2.5, #axis tick label size
  axis.font = 1,
  axis.rotate = TRUE,
  axis.tick = TRUE,
  axis.lwd = 1,
  ticks.lwd = 1,
  ticks.length = 0.025,
  axis.col = "#000000",
  ticks.col = "#A9A9A9"
)

TernaryPoints(df0[df0$BSH_CODE=="A5.1",][,7:9], col = 1, bg = cbPalette[2], pch=21, cex=3)
TernaryPoints(df0[df0$BSH_CODE=="A5.2",][,7:9], col = 1, bg = cbPalette[3], pch=22, cex=3)
TernaryPoints(df0[df0$BSH_CODE=="A5.3",][,7:9], col = 1, bg = cbPalette[4], pch=23, cex=3)
TernaryPoints(df0[df0$BSH_CODE=="A5.4",][,7:9], col = 1, bg = cbPalette[5], pch=24, cex=3)

legend(
  # x=.2,y=.905,
  x=.175,y=.905,
  # "right", 
       legend = c("A5.1 Subtidal coarse sediment",
                  "A5.2 Subtidal sand",
                  "A5.3 Subtidal mud",
                  "A5.4 Subtidal mixed sediments"),
       cex = 1.75, bty = "n", pch = 21:24, pt.cex = 2,
       pt.bg = c(cbPalette[2],
                 cbPalette[3],
                 cbPalette[4],
                 cbPalette[5]))

dev.off()

# tidy up ####
rm(df,df0,pl, cbPalette,ppi,x,xpos,y,ypos,z, libfolder)

detach("package:tidyverse",unload=TRUE)
detach("package:Ternary",unload=TRUE)
