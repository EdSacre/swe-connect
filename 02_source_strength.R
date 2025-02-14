# Import packages
library(terra)
library(rio)
library(Matrix)
options(scipen = 999)

# This script can safely be looped as it is much quicker to run than the previous scripts
for(i in 1:16){
  
  # Define extent
  ex <- c(350000, 1000000, 6091807, 7424396) # full extent
  
  # Specify output directory
  outdir <- "outputs/connectivity_source_strength/"
  if(!dir.exists(outdir)){dir.create(outdir)}
  
  # Import species spreadsheet
  meta <- rio::import("spec_list.xlsx", sheet = "species")
  specname <- meta$code[i]
  grid <- rast("grid.tif")
  
  # Import species habitats
  spec <- rast(paste0("./inputs/habs_standard/", specname, ".tif"))
  spec[is.na(grid)] <- NA
  
  # Import connectivity matrix without pressures
  p1 <- paste0(".*con_", specname, ".*\\.mtx$")
  l1 <- list.files("outputs/connectivity_matrices", pattern = p1, full.names = TRUE, recursive = T)
  l1 <- l1[grepl("presFALSE", l1)]
  if(identical(l1, character(0))){next}
  mat <- readMM(l1)
  
  # Calculate source strength over X generations
  mat <- mat %*% mat %*% mat # modify the number of generations here
  source_str <- rowSums(mat)
  source_str[source_str < 0] <- 0 # Remove small rounding errors
  source_str <- source_str / max(source_str) # Scale from 0 to 1
  
  # Add to spatial raster
  spec_source_str <- spec
  spec_source_str[!is.na(spec_source_str)] <- source_str
  spec_source_str[is.na(spec_source_str)] <- 0
  spec_source_str[is.na(grid)] <- NA
  
  # Export files
  outfile <- paste0(outdir, "/", specname, "_source_strength.tif")
  writeRaster(spec_source_str, outfile, overwrite = TRUE)
  cat("\nIteration", i, "complete - Species:", meta$english[i], " - Species code: ", meta$code[i])
}