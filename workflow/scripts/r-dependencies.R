## This script installs any packages that we could not install with conda. Note that because the
## conda environment is currently active, these will be installed into the corresponding R library.

if (!"exactextractr" %in% installed.packages()) {
  install.packages("exactextractr", repos='https://cloud.r-project.org')
}

if (!"qmap" %in% installed.packages()) { 
  install.packages("qmap", repos='https://cloud.r-project.org')
}
## if (!"terra" %in% installed.packages()) {
##   install.packages("terra", repos='https://cloud.r-project.org')
## }

## if (!"sf" %in% installed.packages()) {
##   install.packages("sf", repos='https://cloud.r-project.org')
## }

## if (!"raster" %in% installed.packages()) {
##   install.packages("raster", repos='https://cloud.r-project.org')
## }
