
#------------------------------------------------------------------------------------------------------
#   MCC-HEATPTB FIRST STAGE ANALYSIS – FAKE DATA DEMO
#------------------------------------------------------------------------------------------------------
# This script performs the first-stage modeling for the MCC-HEATPTB project using
# example data (dlistFake). The core function is 'fit_fs'—see details in 'functionsTS.R'.
# It checks for errors during model fitting and saves the results.
# Core function: fit_fs
# Input data: dlistFake (list of Time series); mcitiesFake (metadata)
# Output data: fs_res
#-------------------------------------------------------------------------------------------------------



# 0 # REMOVE OLD OBJECTS AND LOAD PACKAGES -------------------------------------------------------------
rm(list = ls())  # Clean the workspace

library(pacman)
pacman::p_load(tidyverse, MASS, splines, mgcv,gnm, dlnm, mixmeta)



# 1 # LOAD FUNCTIONS AND DATA --------------------------------------------------------------------------

source("functionsTS.R")   # Load time-series functions (including fit_fs)
load("dataFake.RData")    # Load example dataset

dlist    <- dlistFake    # Assign example list of data
mcities0 <- mcitiesFake  # Assign example metadata



# 2 # FIRST STAGE MODEL FITTING AND ERROR CHECK --------------------------------------------------------

# Fit first-stage models using 'fit_fs', wrapped with 'safely' to detect errors.
fs_resProv <- map(dlist, safely(fit_fs))

## a) Check for errors
errors <- map(fs_resProv, function(x) {
  df <- data.frame(success = is.null(x$error))
  return(df)}) %>% bind_rows(.id = "id")

errors %>% filter(success == FALSE)  # Display any failed fits

## b) If there are no errors (<0 rows>), extract the results and remove provisional objects
fs_res <- map_depth(fs_resProv, 1, function(x) return(x$result))



# 4 # SAVE RESULTS ------------------------------------------------------------------------------------
# Extract 'mcities' information and save results
mcities_fs <- map(fs_res, function(x) { return(x$mcities) }) %>%
  bind_rows()

rm(list = setdiff(ls(), c("dlistFake", "fs_res", "mcities_fs")))
save.image("FSFake.RData")



#-------------------------------------------------------------------------------------------------------
# END OF SCRIPT