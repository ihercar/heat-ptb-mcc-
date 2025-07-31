
#------------------------------------------------------------------------------------------------------
# MCC-HEATPTB SECOND STAGE META-ANALYSIS (FAKE DATA)
#------------------------------------------------------------------------------------------------------
# Joint Meta-analysis for stratified analyses in the categories of each maternal/fetal characteristic. 
# Core functions: fpredmv, plot_RR 
# Core data: Results from 1. FirstStage -> FSFake.RData. Concretely those for SEX and MAGEA
#------------------------------------------------------------------------------------------------------



# 0 # CLEAN ENVIRONMENT AND LOAD PACKAGES -------------------------------------------------------------
rm(list = ls())
library(pacman)
pacman::p_load(tidyverse, splines, mgcv, dlnm, gnm, mixmeta, ggplot2, ggthemes, ggtext, patchwork, scales)


# 1 # LOAD FUNCTIONS AND PARAMETERS -------------------------------------------------------------------
source("functionsTS.R")
source("parameters.R")



# 2 # LOAD FIRST STAGE RESULTS ------------------------------------------------------------------------
load("FSFake.RData")



# 3 # ORGANIZE EACH CHARACTERISTIC WITHIN A LIST   ----------------------------------------------------
dCHAR <- mcities_fs%>%filter(Cat!="EP" & !grepl("all", Sub))
rownames(dCHAR) <- dCHAR$Ca.Su.Co.Ci
lCHAR <- split(dCHAR,f=factor(dCHAR$Cat))%>%map(function(x){ x<-x%>%mutate(CoSub=interaction(Co, Sub)) 
                                                             return(x) })



# 4 # # EXTRACT COEFS AND FIT META-ANALYSIS (JOINT CATEGORICAL) FOR EACH CHARACTERISTIC ---------------
mvallCHARg <- map(lCHAR, function(x){
                coefall <- as.matrix(x[,grep("coef", names(x))])
                vcovall <- as.matrix(x[,grep("vcov", names(x))])
                mvallCHAR <- mixmeta(coefall~rgtmean_s*Sub,vcovall,x,
                                     random=~1|CoSub/Ci,control=list(showiter=F, igls.inititer=15))
                return(mvallCHAR)})



# 5 # P-VALUES AND GOODNESS OF FIT SUMMARY ------------------------------------------------------------
metaCHAR.gof <- map(lCHAR, function(x){
                  coefall <- as.matrix(x[,grep("coef", names(x))])
                  vcovall <- as.matrix(x[,grep("vcov", names(x))])
                  mvallCHAR <- mixmeta(coefall~rgtmean_s*Sub,vcovall,x,
                                      random=~1|CoSub/Ci,control=list(showiter=F, igls.inititer=15))
                  sg <- summary(mvallCHAR)
                  sgml <- summary(update(mvallCHAR, method="ml"))
                  df <- data.frame(Subcat="All", type="jointbycat1", I2=sg$i2stat[[".all"]], 
                  pQ=sg$qstat$pvalue[[".all"]], prgtmean=fwald(mvallCHAR,"rgtmean"), AIC=sgml$AIC,
                  pcomplete=fwald(mvallCHAR,"Sub"), pmain=fwald2(mvallCHAR,"Sub",discount="rgtmean_s"))
                  return(df)})%>%bind_rows(.id="Cat") %>% mutate( fe="rgtmean_s*Sub",re="1|CoSub/city")



# 6 # PREDICTION FROM MIXMETA (JOINT CATEGORICAL) -----------------------------------------------------
EECHARg <- imap(mvallCHARg, function(x,name){
                   mSub<-dCHAR%>%filter(Cat == name)%>%group_by(Sub)%>% summarize(rgtmean_s=median(rgtmean_s), 
                                                                                  Sub=unique(Sub),
                                                                                  Cat=unique(Cat))
                   EEi<-predict(x, newdata=mSub, vcov=T)
                   names(EEi)<-mSub$Sub
                   EE <- map(EEi, function(x){ fpredmv(x)})%>% bind_rows(.id = "Sub")})%>%bind_rows(.id="Cat")
EECHARg <- EECHARg%>%mutate(Type="jointbycat") 



# 7 # PLOTTING -----------------------------------------------------------------------------------------
theme_set(theme_cc2())

RRCHAR0 <- filter(EECHARg, x %in% c(0.75, 0.95))%>%dplyr::select(Cat:highy)

fake <- expand.grid(Cat=unique(RRCHAR0$Cat), x=c(0.75,0.95))%>%mutate(Sub=Cat, y=NA, lowy=1, highy=1)
fake <- fake[,names(RRCHAR0)]
RRCHAR1 <- bind_rows(RRCHAR0,fake)%>%arrange(Cat)

heat_points <- tibble(x = c(0.75, 0.95), labheat = c("Moderate (%75)", "Extreme (%95)"))
RRCHAR2 <- left_join(RRCHAR1, heat_points, by = "x")

lev <- c("1.24-","2.[25,34]","3.35+","AGEA","1.Female","2.Male","SEX")
lab <- c("24-","[25,34]","35+","***AGE***","Female","Male","***GENDER***") 
RRCHAR3 <- RRCHAR2%>%mutate(Sub=factor(Sub, levels=lev, labels=lab),dummy = str_detect(Sub, "\\*") )

RRCHAR4 <- arrange(RRCHAR3, -as.numeric(Sub)) %>% mutate(dummy = ifelse(dummy == TRUE, as.character(Sub), NA) 
           %>% zoo::na.locf() %>% str_remove_all("\\*"),labheat = fct_rev(labheat))

RRCHAR <- RRCHAR4 %>% filter(str_detect(Sub, "\\*", negate = TRUE))


p<-plot_RR(RRCHAR)
p<-p+facet_grid(dummy ~ ., scales = "free_y", space = "free_y") + 
      ggtitle("Heat-related RR by mother age and fetal Gender (Fake data).")
p
ggsave("figures/3.CHAR/RRCHAR.png", height = 7, width = 9, bg = "white", type = "cairo")
 
 

# 8 # CLEANUP AND SAVE --------------------------------------------------------------------------------
rm(list = setdiff(ls(), c("dCHAR","mvallCHARg","metaCHAR.gof","EECHARg", "RRCHAR")))
save.image("ResCHARfake.RData")



#------------------------------------------------------------------------------------------------------
# END OF SCRIPT
#------------------------------------------------------------------------------------------------------

