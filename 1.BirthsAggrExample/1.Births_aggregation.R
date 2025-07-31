#------------------------------------------------------------------------------------------------------
#   MCC-HEATPTB SCRIPT FOR AGGREGATION – COUNTRY: Co1
#------------------------------------------------------------------------------------------------------
# This script prepares time series of preterm births (PTB) for the basic (standard PTB) and additional 
# analyses (other endpoints of gestational duration and stratified analyses by fetal and maternal characteristics).
# Core function: faggr2
# Input data: dBCo1
# Output data: dTSCo1
#------------------------------------------------------------------------------------------------------



# 0 # REMOVE OLD OBJECTS AND LOAD PACKAGES -------------------------------------------------------------
rm(list = ls())  # Clean the workspace

library(pacman)
pacman::p_load(tidyverse, readr, fs, tsModel, lubridate, tsibble, foreign)



# 1 # LOAD DATA AND FUNCTIONS -------------------------------------------------------------------------

load("dBCo1.RData")    # Load main dataset
source("faggr.R")      # Load aggregation functions

ls()  # Expected: "dBCo1" "faggr2" 



# 2 # DATA PREPARATION -------------------------------------------------------------------------------
# We start with a list called 'dlist', where each element is a data.frame corresponding to one city/location.
# This section harmonizes variable names and creates subgroups.

Conam <- "Co1"
dlist <- split(dBCo1, factor(dBCo1$Ci))

# ---- FETAL GENDER (type: character) ----
# Original variable: Fetal.gender; values: "2.female" and "1.male"
dlistSEXp1 <- map(dlist, function(x) {  
  x$SEX <- x$Fetal.gender
  lab <- "1.male"
  x <- x[complete.cases(x$SEX), ]
  x$SEX <- factor(ifelse(x$SEX == lab, "2.Male", "1.Female"))
  x <- x %>% select(BDATE, GWEEK, SEX)
  l <- split(x, x$SEX)
  l[[3]] <- x
  names(l)[3] <- "allSEX"   
  return(l)
})

dlistSEXp2 <- map_depth(dlistSEXp1, 1, bind_rows, .id = "category") %>%
  bind_rows(.id = "city") %>%
  mutate(city = factor(city),
         category = factor(category))

dlistSEX <- split(dlistSEXp2, dlistSEXp2$category) %>%
  map(function(x) split(x, x$city))

# ---- MATERNAL AGE (type: numeric) ----
# Original variable: Mother.age
dlistMAGEAp1 <- map(dlist, function(x) {
  x$MAGE <- x$Mother.age
  x <- x[complete.cases(x$MAGE), ]
  breaks <- unique(c(min(min(x$MAGE), 24), 24, 34, max(max(x$MAGE), 34))) # <25, 25-34, ≥35
  x$MAGE_3CA <- cut(x$MAGE, breaks, include.lowest = TRUE)
  x <- x %>% select(BDATE, GWEEK, MAGE_3CA)
  l <- split(x, x$MAGE_3CA)
  l[[4]] <- x  
  return(l)
})

dlistMAGEAp1 <- lapply(dlistMAGEAp1, function(x) {
  names(x) <- c("1.24-", "2.[25,34]", "3.35+", "allMAGEA")
  return(x)
})

dlistMAGEAp2 <- map_depth(dlistMAGEAp1, 1, bind_rows, .id = "category") %>%
  bind_rows(.id = "city") %>%
  mutate(city = factor(city),
         category = factor(category))

dlistMAGEA <- split(dlistMAGEAp2, dlistMAGEAp2$category) %>%
  map(function(x) split(x, x$city))



# 3 # AGGREGATION ---------------------------------------------------------------------------------------
# We aggregate the data into time series per city for basic and extended analyses.

dlist1 <- map(dlist, function(x) x %>% select(BDATE, GWEEK))

# ---- BASIC ANALYSIS: Preterm [22,37) weeks
dpret37 <- map(dlist1, faggr2)
names(dpret37) <- names(dlist1)

# ---- OTHER GESTATIONAL ENDPOINTS
dpret28    <- map(dlist1, faggr2, wfdef = 27, widef = 22)   # [22,28) weeks
dpret28_32 <- map(dlist1, faggr2, wfdef = 31, widef = 28)   # [28,32) weeks
dpret32_37 <- map(dlist1, faggr2, wfdef = 36, widef = 32)   # [32,37) weeks
dpret37_39 <- map(dlist1, faggr2, wfdef = 38, widef = 37)   # [37,39) weeks
dfullt     <- map(dlist1, faggr2, wfdef = 42, widef = 39)   # ≥39 weeks

names(dpret28) <- names(dpret28_32) <-names(dpret32_37) <-names(dpret37_39) <-names(dfullt) <-names(dlist1)


# ---- STRATIFIED BY MATERNAL AND FETAL CHARACTERISTICS
dpret37_bySEX   <- map_depth(dlistSEX, 2, faggr2)    # By fetal gender
dpret37_byMAGEA <- map_depth(dlistMAGEA, 2, faggr2)  # By maternal age

# 4 # ORGANIZE AND SAVE -------------------------------------------------------------------------------------
# Create a nested list for extended analyses
dpretEP <- list(
  w28    = dpret28,
  w28_32 = dpret28_32,
  w32_37 = dpret32_37,
  w37    = dpret37,
  w37_39 = dpret37_39,
  wfullt = dfullt
)

dpretExtended <- list(
  EP   = dpretEP,  
  SEX  = dpret37_bySEX, 
  AGEA = dpret37_byMAGEA
)


# ---- Flatten nested list into a single-level dataframe
dlong <- map_depth(dpretExtended, 2, function(x) bind_rows(x, .id = "Ci")) %>%
  map_depth(1, function(x) bind_rows(x, .id = "Sub")) %>%
  bind_rows(.id = "Cat")

dTSCo1 <- dlong %>%
  mutate(Co = Conam,
         Cat.Sub.Co.Ci = paste(Cat, Sub, Co, Ci, sep = ".")) %>%
  select(Cat.Sub.Co.Ci, Cat, Sub, Co, Ci, everything())

save.image("dTSCo1.RData")



#-----------------------------------------------------------------------------------------------------------------------------------
#  NEXT STEPS:
#
#  1) Merge 'dTSCo1' with temperature data, define the 5 hottest months
#     (Mysumm_ini, Mysumm_end), and create a summer grouping variable (vgroup),
#     even when summer spans two calendar years:
#        aux <- ifelse(x$mm %in% (Mysumm_ini-1):12, x$yy + 1, x$yy)
#        vgroup <- if (Mysumm_end > 5) x$yy else aux
#
#  2) 'dTSCo1' now contains the time series required for first-stage
#     analysis in all available cities of Co1.
#
#  3) Standardize the 'mcitiesCo1' metadata dataframe.
#
#  4) Combine all 'dTSCoi' and split as list 'dlist'.
#     Bind all 'mcitiesCi' objects and save them as 'data.RData'.
#-----------------------------------------------------------------------------------------------------------------------------------

# END OF SCRIPT