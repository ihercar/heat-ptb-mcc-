
#------------------------------------------------------------------------------------------------------
# MCC-HEATPTB SECOND STAGE META-ANALYSIS (FAKE DATA)
#------------------------------------------------------------------------------------------------------
# Meta-analysis for standard preterm birth (PTB): overall curve, also summarized by Country & Climate 
# zone. Provides curves and RR for moderate and extreme heat.
# Core functions: fpredmv, plot_curves, plot_RR --see details in functionsTS.R.
# Core data: - Results from the first stage (1. FirstStage -> FSFake.RData). Concretely those for standard PTB
#            - CB parameters (parameters.R)
#------------------------------------------------------------------------------------------------------

# 0 # CLEAN ENVIRONMENT AND LOAD PACKAGES -------------------------------------------------------------
rm(list = ls())
library(pacman)
pacman::p_load(tidyverse, splines, mgcv, dlnm, gnm, mixmeta, ggplot2, ggthemes, ggtext, patchwork, scales)



# 1 # LOAD FUNCTIONS AND PARAMETERS ------------------------------------------------------------------
source("functionsTS.R")
source("parameters.R")



# 2 # LOAD FIRST STAGE RESULTS -----------------------------------------------------------------------
load("FSFake.RData")



# 3 # PREPARE META-ANALYSIS DATA ----------------------------------------------------------------------
x <- mcities_fs %>% filter(Cat == "EP", Sub == "w37")
coefall <- as.matrix(x %>% dplyr::select(starts_with("coef")))
vcovall <- as.matrix(x %>% dplyr::select(starts_with("vcov")))
dPTB <- x %>% dplyr::select(-starts_with("coef"), -starts_with("vcov"))



# 4 # META-ANALYTICAL MODEL ----------------------------------------------------------------------------
mvall <- mixmeta(coefall ~ rgtmean_s, vcovall, dPTB, random = ~1 | Co / Ci, 
                 control = list(igls.inititer = 20))
fwald(mvall, "rgtmean_s")  # p-value for slope (fake data)
summary(mvall)



# 5 # EXTRACT MODEL OBJECTS ----------------------------------------------------------------------------
gof_mvallPTB <- tibble(         #goodness of fit summary
  I2 = summary(mvall)$i2stat[[".all"]],
  Q = summary(mvall)$qstat$Q[[".all"]],
  pQ = summary(mvall)$qstat$pvalue[[".all"]]
)

coef_mvallPTB <- list(mvcoef = coef(mvall), mvvcov = vcov(mvall))  #coefs and vcovs

blup_mvallPTB <- blup(mvall,vcov=T)  #blup
names(blup_mvallPTB)<-dPTB$Co.Ci



# 6 # PREDICTION FROM MIXMETA (FIXED EFFECTS)-----------------------------------------------------------

#____ 6.1 # all 
mall <- dPTB%>%summarize(avgtmean_s=median(avgtmean_s), rgtmean_s=median(rgtmean_s))
EEi <- predict(mvall, newdata=mall, vcov=T)
EEall <- fpredmv(EEi); rownames(EEall)<-NULL
EEall <- EEall%>%mutate( Cat="All", Sub="All",)%>%dplyr::select( Cat, Sub, x,y, lowy, highy, everything())


#____ 6.2 # by country 
metaCo <- dPTB%>%group_by(Co)%>% summarize(rgtmean_s=median(rgtmean_s))
nCo <- unique(dPTB$Co)  
lisCo <- split(nCo, factor(nCo))
EECo <- map( lisCo, function(x){metavi<-filter(metaCo,metaCo$Co==x)
          EEi <- predict(mvall, newdata=metavi, vcov=T)
          fpredmv(EEi)})%>%bind_rows(.id="Sub")
EECo <- EECo%>%mutate(Cat="countryname")%>%dplyr::select(Cat, Sub, x,y, lowy, highy, everything())
rownames(EECo)<-NULL

#____ 6.3 # by kgclzone1
metaKz <- dPTB%>%group_by(kgclzone1)%>% summarize(rgtmean_s=median(rgtmean_s))
nKz <- unique(dPTB$kgclzone1)  
lisKz <- split(nKz, factor(nKz))
EEKz <- map( lisKz, function(x){metavi<-filter(metaKz,metaKz$kgclzone1==x)
          EEi <- predict(mvall, newdata=metavi, vcov=T)
          fpredmv(EEi)})%>%bind_rows(.id="Sub")
EEKz <- EEKz%>%mutate(Cat="kgclzone1") %>%dplyr::select( Cat, Sub, x,y, lowy, highy, everything())
rownames(EEKz)<-NULL


# 7 # PLOTTING ----------------------------------------------------------------------------------------
theme_set(theme_cc2())

#____ 7.1 # Overall curves
p<-plot_curve(EEall)+ ggtitle("Heat-PTB association (Fake data).")  #All
p
ggsave("figures/1.PTB/CurvePTB.png", height = 7, width = 9, bg = "white", type = "cairo")

p<-plot_curve(EECo)+facet_wrap(Sub ~ ., nrow = 2, scales = "free_y") #By country
p<-p+ggtitle("Heat-PTB association by country (Fake data).")
p
ggsave("figures/1.PTB/CurvePTBbyCo.png", height = 7, width = 9, bg = "white", type = "cairo")

p<-plot_curve(EEKz)+facet_wrap(Sub ~ ., nrow = 2, scales = "free_y") #By climate zone
p<-p+ggtitle("Heat-PTB association by climate zone (Fake data).")
p
ggsave("figures/1.PTB/CurvePTBbyKz.png", height = 7, width = 9, bg = "white", type = "cairo")


#____ 7.2 # Relative Risks 

heat_points <- tibble(x = c(0.75, 0.95), labheat = c("Moderate (%75)", "Extreme (%95)"))
EE <- bind_rows(EEall, EEKz, EECo)
RR <- filter(EE, x %in% c(0.75, 0.95))
RRall_co_cli <- left_join(RR, heat_points, by = "x") %>%
  mutate(
      Sub = case_when(Sub == "All" ~ "All",
      Sub == "A" ~ "Tropical (A)",
      Sub == "B" ~ "Arid (B)",
      Sub == "C" ~ "Temperate (C)",
      Sub == "D" ~ "Continental (D)",
      TRUE ~ Sub),
      Cat = factor(Cat, levels = c("All", "countryname", "kgclzone1"), 
                   labels = c("", "By country", "By climate zone")),
      labheat = fct_rev(labheat))

p<-plot_RR(RRall_co_cli)
p<-p+facet_grid(Cat ~ ., scales = "free_y", space = "free_y") + ggtitle("Heat-PTB association (Fake data).")
p
ggsave("figures/1.PTB/RRAll_Co_Kz.png", height = 7, width = 9, bg = "white", type = "cairo")


# CLEANUP AND SAVE ------------------------------------------------------------------------------------
mvallPTB<-mvall
rm(list = setdiff(ls(), 
                  c("EE", "RR", "dPTB", "gof_mvallPTB", "mvallPTB", "coef_mvallPTB","blup_mvallPTB")))
save.image("ResPTBFake.RData")



#------------------------------------------------------------------------------------------------------------
# END OF SCRIPT
#------------------------------------------------------------------------------------------------------------
