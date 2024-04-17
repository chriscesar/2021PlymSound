# 2021PlymSound
Analysis of data for Plymouth Sound

## TO DO list ##
* ~~Ternary plot: increase font size for axes and legend~~
* Increase text size on nMDS plots
* Relative abundance plots: Increase all text sizes and contrast plot colours
  * use updated data in openxlsx::readxlsx(paste0(fol,"Data/2021 data/Infauna/Updated/"),
sheet = "Matched sites"); These data have been matched up by location (using 2021 site names). with Tamar & Yealm sites removed
* "Site" column indicates which samples have replicates. e.g., "PLYMS11a,b and c"
  * Before analysis, will need to calculate mean abundances (including 0s!) across reps
* reproduce MDS plots using values in "2021 matched site" to display sites instead of points
  * need to generate a new site name that includes the sampling year
* For MDS plots, instead of linking samples to a centroid (e.g., by year), link samples by sample site between the different years
* Perform ANOSIM, using site ID as a predictor (might need to chat with customer)
  * current version analyses with mvabund.  Consider retaining that analysis and discuss the difference (if any)?
* Need to see how sediments are changing at specific sites
  * maybe reproduce Ternary plots? Linking samples from same station, but different years?
  * could also run RDA plot with PSA data as a predictor?
  * could run gllvm with PSA data as a predictor? (or possibly incorporate taxon data and HMSC?)
