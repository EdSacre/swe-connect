# NOTE: These models are extremely time consuming, so it is only recommended on
# a powerful computer or supercomputer cluster. If you use the smaller extent for testing
# it will override the model outputs currently in the outputs folder, so be sure not to
# override these if you want to test the subsequent scripts on the full extent.
# I would recommend running all subsequent scripts (02 - 10) before testing these models
# to avoid breaking the code. Models take 1 or two days to run for a single species
# on a single core on a normal computer, depending on the species.

# Packages
library(gdistance)
library(raster)
library(rio)
library(Matrix)
options(scipen = 999)

# Define extent
ex <- c(350000, 1000000, 6091807, 7424396) # full extent
#ex <- c(660000, 720000, 6560000, 6620000) # small extent for testing

# Specify output directory
outdir <- "outputs/connectivity_matrices/"
if(!dir.exists(outdir)){dir.create(outdir)}

# Import species spreadsheet
meta <- rio::import("spec_list.xlsx", sheet = "species")

# Set run parameters
#n.cores <- 4 # Number of cores to use for parallel processing
usepres <- TRUE # Set pressures as TRUE or FALSE, TRUE = include pressures in cost layer
habsonly <- TRUE # Only calculate connectivity to other habitat cells

# Select species (iterate in background for pseudo-parallel, i.e. 1 through 16)
# You can do this using "Source as background job" in RStudio and changing k 
# to 1,2,3...16 for each background job
# Here I have set it up so you can place it within a for loop if you like,
# i.e. for(k in 1:16){...}
# Or run 16 background processes in parallel for faster results
k <- 1

# Specify input and output files
specname <- meta$code[k]
max.dist <- meta$dispersal_distance_km[k]*1000
dir.create(paste0(outdir, specname))
namerun <- paste0(specname,"_250m")
outfile <- paste0(outdir, specname, "/con_", namerun,
                  "_", max.dist/1000, "km_pres", usepres,
                  "_habsonly", habsonly, ".mtx")
input_habs <- paste0("./inputs/habs_standard/", specname, ".tif")

# Import cost/destination/goal raster
base.cost <- raster::raster("grid.tif")+1
base.cost[base.cost == 0] <- NA
base.cost <- raster::crop(base.cost, ex)
base.cost <- raster::extend(base.cost, ex)
plot(base.cost)

# Import species habitats
species <- raster::raster(input_habs)
species <- raster::crop(species, ex)
species <- raster::extend(species, ex)
species[species == 0] <- NA
species[is.na(base.cost) == TRUE] <- NA

# Import and create physical disturbance raster
pres <- raster::raster("inputs/pres_standard/metria_connectivity_disturbance_250m.tif")
pres[is.na(pres)] <- 1
pres[is.na(base.cost)] <- NA
pres <- raster::crop(pres, ex)
pres <- raster::extend(pres, ex)
if(usepres){base.cost <- pres}

# Convert origin and goal rasters to points
species.o <- raster::rasterToPoints(species, spatial = TRUE) # Set the origin points
species.g <- raster::rasterToPoints(base.cost, spatial = TRUE) # Set the destination points
if(habsonly){species.g <- species.o} # If only calculating connectivity to other habitat cells
names(species.o) <- "habvalue"
species.o$habvalue <- species.o$habvalue / max(species.o$habvalue) # Scale habitat values to a max of 1
length(species.o[]) # Number of habitat cells
length(species.g[]) # Number of total cells


# Create transition layer for cost distance
species.tr <- gdistance::transition(base.cost, transitionFunction = min, directions = 8) # Transition file
species.tr <- gdistance::geoCorrection(species.tr, type = "c") # Geo-correct transition file

time1 <- Sys.time() # Record start processing time

# Specify parameters for the model
alpha <- 0.3 # This defines how steep the dispersal kernel is - a value of 1 is almost linear
d <- max.dist * alpha
a <- 1/d # this specifies that the dispersal kernel should taper towards 0 around the max dist

ncells <- length(species.g[])
nhabs <- length(species.o[])
speco <- species.o
specg <- species.g
spectr <- species.tr
source_only <- habsonly

if(exists("cd_mat")){rm(cd_mat)}
for(i in 1:nhabs){
  species.cd <- gdistance::costDistance(spectr, speco[i,], specg)
  cd.temp <- species.cd[1,]
  cd.temp[cd.temp > max.dist] <- NA
  cd.temp[is.na(cd.temp) == FALSE] <- exp(-a*(cd.temp[is.na(cd.temp) == FALSE]))
  #cd.temp[cd.temp > 0] <- 1 # Add this line to set all connectivity values equal (no dispersal kernel)
  cd.temp[is.na(cd.temp) == TRUE] <- 0
  cd.temp <- cd.temp * speco[i,]$habvalue # Add this line to give higher connectivity values to higher habitat quality/coverage
  cd.temp <- Matrix(cd.temp, sparse = TRUE, nrow = 1)
  
  if(exists("cd_mat")){cd_mat <- rbind(cd_mat, cd.temp)}
  if(!exists("cd_mat")){cd_mat <- cd.temp}
  
  cat("\nIteration ", i, "of ", nhabs, " done - Species:", specname)
}

writeMM(cd_mat, file = outfile)
