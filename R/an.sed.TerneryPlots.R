#######################################################################
###########           MCZ Monitoring Reporting              ###########
###########              Ternary plots                      ###########
#######################################################################


#### Load libraries ####
require(xlsx)
require(openxlsx)
require(ggplot2)
require(ggtern)


#### Load the Baseplot.RData workspace saved together with this script ####
load(file.choose())


#### Read data ####

# Data should be in an Excel sheet that includes the columns Gravel, Sand and Mud

# Choose file path
stfile = file.choose()
filename = strsplit(stfile,'\\',fixed=TRUE)[[1]][length(strsplit(stfile,'\\',fixed=TRUE)[[1]])]
dir = substr(stfile, 1, nchar(stfile)-nchar(filename))

# Read and inspect data sheet names
stsheet <- getSheetNames(stfile)
stsheet

# Load data from the selected sheet
df <- read.xlsx(stfile, sheet = stsheet[1], # sheet = 'NameOfSheet' or sheet = stsheet[number of sheet]
                detectDates = TRUE,
                colNames = TRUE)
# Check column names include Gravel, Sand and Mud
names(df)


#### Produce ternary plots in the same directory as input data ####

jpeg(filename = paste(dir,'TernaryPlot.jpg',sep=""), width = 14, height = 12,units = 'cm',res = 600)
Baseplot +
  geom_point(data = df,size=1.5,shape=21) 
dev.off()

jpeg(filename = paste(dir,'TernaryPlot2.jpg',sep=""), width = 14, height = 12,units = 'cm',res = 600)
Baseplot2 +
  geom_point(data = df,size=1.5,shape=21) 
dev.off()


