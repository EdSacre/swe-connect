# You will need to install the Gurobi optimizer to run this code
# Intructions to do so are here:
# https://prioritizr.net/articles/gurobi_installation_guide.html

# Load packages
library(terra)
library(prioritizr)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rio)
library(ggthemes)

# Create necessary directories
if(!dir.exists("outputs/prioritization_expansion_raster")){dir.create("outputs/prioritization_expansion_raster", recursive = T)}
if(!dir.exists("outputs/prioritization_expansion_table")){dir.create("outputs/prioritization_expansion_table", recursive = T)}

# Create base planning units
pu <- rast("grid.tif")
pu <- pu + 1

# Prioritization features
meta <- rio::import("spec_list.xlsx", sheet = "species")
l1 <- list.files("outputs/connectivity_source_strength", pattern = "\\.tif$", full.names = TRUE)
l2 <- list.files("outputs/connectivity_vulnerability", pattern = "\\.tif$", full.names = TRUE)
l3 <- list.files("inputs/habs_standard", pattern = "\\.tif$", full.names = TRUE)
l <- c(l1,l2, l3)

n1 <- list.files("outputs/connectivity_source_strength", pattern = "\\.tif$")
n2 <- list.files("outputs/connectivity_vulnerability", pattern = "\\.tif$")
n3 <- list.files("inputs/habs_standard", pattern = "\\.tif$")
n <- c(n1,n2,n3)
n <- tools::file_path_sans_ext(n)

feat <- rast(l)
names(feat) <- n
feat <- feat * pu

# Import MPA network
mpa <- rast("outputs/mpas/mpas.tif")
mpa <- app(mpa, fun = "sum")
mpa <- mpa*pu

# Determine the budget (in terms of area)
# Incrementally increasing the budget in increments of 1% until all
# objectives are achieved (100% protection of all features)
perc <- seq(0.11, 1.0, by = 0.01)
exp_budget <- sum(values(pu, na.rm = T)) * perc

# Out of MPA planning units
temp <- mpa
temp[temp == 1] <- NA
temp[temp == 0] <- 1
out_mpa <- temp 

# Run expansion prioritization scenario
for(i in 1:length(exp_budget)){
  p1 <- problem(x = pu, features = feat) %>%
    add_min_shortfall_objective(budget = exp_budget[i]) %>%
    add_relative_targets(1.0) %>%
    add_locked_in_constraints(mpa) %>% # Lock in existing protected areas
    add_gurobi_solver(gap = 0.01, verbose = FALSE)
  
  p1_sol <- solve(p1, force = T, run_checks = F)
  df <- eval_feature_representation_summary(p1, p1_sol)
  colnames(df)[1] <- "species"
  df$species <- rep(meta$english, 3)
  df$feature <- NA
  df$feature[1:16] <- "Source connectivity"
  df$feature[17:32] <- "Vulnerability"
  df$feature[33:48] <- "Habitat"
  write.csv(df, paste0("outputs/prioritization_expansion_table/feature_representation_", perc[i]*100, "percent.csv"), row.names = FALSE)
  
  p1_sol <- p1_sol * out_mpa
  p1_sol[p1_sol == 0] <- NA
  writeRaster(p1_sol, paste0("outputs/prioritization_expansion_raster/prioritization_expansion_", perc[i]*100, "percent.tif"), overwrite = TRUE)
  
  if(i == 1){
    p_inc <- p1_sol
  } else {
    p_inc <- c(p_inc, p1_sol)
  }
  
  if(sum(df$relative_held) == 48){break}
}

# Priority map based on expansion
p_inc[is.na(p_inc)] <- 0
prio <- app(p_inc, fun = "sum")
prio[is.na(pu)] <- NA
writeRaster(prio, "outputs/prioritization_expansion_raster/prioritization_expansion_priority.tif", overwrite = TRUE)
