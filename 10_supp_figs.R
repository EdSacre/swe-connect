# Import packages
library(terra)
library(rio)
library(Matrix)
library(ggplot2)
library(gridExtra)
options(scipen = 999)

# This script can safely be looped as it is much quicker to run than the previous scripts

plots <- list()

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
  if(nrow(mat) <= 10000){
    mat <- mat %*% mat %*% mat # modify the number of generations here
    source_str <- rowSums(mat)
    
    mat <- readMM(l1)
    mat[mat > 0] <- 1
    mat <- mat %*% mat %*% mat # modify the number of generations here
    source_deg <- rowSums(mat)
  }
  
  if(nrow(mat) > 10000){
    ind <- sample(1:nrow(mat), 10000)
    mats <- mat[ind, ind]
    
    mats <- mats %*% mats %*% mats # modify the number of generations here
    source_str <- rowSums(mats)
    
    mats <- mat[ind, ind]
    mats[mats > 0] <- 1
    mats <- mats %*% mats %*% mats # modify the number of generations here
    source_deg <- rowSums(mats)
  }
  
  #source_deg[source_deg < 0] <- 0 # Remove small rounding errors
  #source_str[source_str < 0] <- 0 # Remove small rounding errors
  df <- data.frame(source_deg = source_deg, source_str = source_str)
  
  plots[[i]] <- ggplot(df, aes(x = source_deg/1000000, y = source_str/1000000)) +
    geom_point(alpha = 0.1) +
    ggtitle(paste(meta$english[i])) +
    xlab("No. of connections") +
    ylab("Connectivity exports") +
    theme_minimal()
  
  cat("Iteration", i, "done\n")
}

jpeg(filename = "figures/connection_number.jpg", width=7000, height=7000, res=600)
grid.arrange(grobs = plots, nrow = 4, ncol = 4)
dev.off()
