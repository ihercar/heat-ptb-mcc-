
#------------------------------------------------------------------------------------------------------
# MCC-HEATPTB SECOND STAGE META-ANALYSIS (FAKE DATA)
#------------------------------------------------------------------------------------------------------
# Longitudinal Joint Meta-analysis for all of Endpoints of PTB and ATB. Aimed of: 
# 1) exploring critical windows of exposure 
# 2) calculating RR for each Endpoint at moderate and extreme heat.
# Core functions:fpredmv 
# Core data: - Results from the first stage (1. FirstStage -> FSFake.RData). Concretely those for each EP
#            - CB parameters (parameters.R) 
#------------------------------------------------------------------------------------------------------


# 0 # CLEAN ENVIRONMENT AND LOAD PACKAGES -------------------------------------------------------------
rm(list = ls())
library(pacman)
pacman::p_load(tidyverse, splines, mgcv, dlnm, gnm, mixmeta, ggplot2, ggthemes, ggtext, patchwork, scales)



# 1 # LOAD FUNCTIONS AND PARAMETERS -------------------------------------------------------------------
source("functionsTS.R")
source("parameters.R")



# 2 # LOAD FIRST STAGE RESULTS -------------------------------------------------------------------------
load("FSFake.RData")



# 3 # EXTRACT COEFS AND METAPREDICTORS -----------------------------------------------------------------
x <- mcities_fs %>% filter(Cat == "EP", Sub!= "w37") 
coefall <- as.matrix(x %>% dplyr::select(starts_with("coef")))
vcovall <- as.matrix(x %>% dplyr::select(starts_with("vcov")))
dEP <- x %>% dplyr::select(-starts_with("coef"), -starts_with("vcov"))
rownames(dEP) <- rownames(coefall) <- rownames(vcovall) <- dEP$Ca.Su.Co.Ci



# 4 # ASSIGN THE VALUE OF GESTATIONAL AGE FOR EACH ENDPOINT AND THE INTERACTION COUNTRY-ENDPOINT--------
sini <- c(22,28,32,37,39)
sfin <- c(27,31,36,38,42)
meanl <- (sini+sfin)/2

tkey <- data.frame(Sub=c("w28", "w28_32", "w32_37", "w37_39", "wfullt"), meanl=meanl)
dEP <- left_join(dEP, tkey, by="Sub")
dEP <- dEP%>%mutate(CoSub=interaction(Co, Sub),Sub=factor(Sub)) %>% mutate(Sub= relevel(Sub, ref = "wfullt"))
rownames(dEP) <- dEP$Ca.Su.Co.Ci  


# 5 # META-ANALYSIS (JOINT LONGITUDINAL)---------------------------------------------------------------
mvallEPl <- mixmeta(coefall~rgtmean_s+ns(meanl, knots=32, Boundary.knots =c(22,42)),vcovall,dEP,
                    random=~1|CoSub/Ci,control=list(showiter=F, igls.inititer=15))


# 6 # P-VALUES AND GOODNESS OF FIT SUMMARY -------------------------------------------------------------

#____ 6.1  calculating some p-values:

## For Gestational length as continuous
(pEPl <- fwald(mvallEPl, "meanl"))    # p-value for ns(meanl)

## For EP as categorical variable (from an alternative model)
mvallEPg <- mixmeta(coefall~rgtmean_s*Sub,vcovall,dEP,random=~1|CoSub/Ci,
                    control=list(showiter=F, igls.inititer=15))
(pEPg_c <- fwald(mvallEPg, "Sub"))    #p-value for the complete effect (Sub+Sub:rgtmean)
(pEPg_m <- fwald2(mvallEPg, "Sub",discount="rgtmean_s"))  #p-value for the main effect (Sub)

## Comparing each Endpoint with wfullt (alternative model)
(pEPw28   <- fwald(mvallEPg, "Subw28"))    
(pEPw28_32 <- fwald(mvallEPg, "Subw28_32")) 
(pEPw32_37 <- fwald(mvallEPg, "Subw32_37")) 
(pEPw37_39 <- fwald(mvallEPg, "Subw37_39")) 


#____ 6.2  Summary of our meta-model:

sl <- summary(mvallEPl); slml <- summary(update(mvallEPl, method="ml"))

metaEP.gof <-data.frame(Type="jointlong", fe="rgtmean_s+ns(meanl,32)",re="1|CoSub/city",
                        I2=sl$i2stat[[".all"]], pQ=sl$qstat$pvalue[[".all"]],
                        AIC=slml$AIC,
                        prgtmean=fwald(mvallEPl,"rgtmean"),pcomplete=pEPl)



# 7 # PREDICTION FROM MIXMETA (JOINT LONGITUDINAL) -----------------------------------------------------

#____ 7.1 # For each week 
mWeek <- expand.grid(rgtmean_s=median(dEP$rgtmean_s), meanl=22:42)
EEi <- predict(mvallEPl, newdata=mWeek, vcov=T)
names(EEi) <- mWeek$meanl 
EEWeek <- map(EEi, function(x){ fpredmv(x)})%>% bind_rows(.id = "Week")

#____ 7.2 # For each endpoint (Subcat)
mSub <- dEP%>%group_by(Sub)%>% summarize(rgtmean_s=median(rgtmean_s),meanl=unique(meanl))
EEi <- predict(mvallEPl, newdata=mSub, vcov=T)
names(EEi) <- mSub$Sub 
EEEPl <- map(EEi, function(x){ fpredmv(x)})%>% bind_rows(.id = "Sub")



# 8 # PLOTTING -----------------------------------------------------------------------------------------
theme_set(theme_cc2())
heat_points <- tibble(x = c(0.75, 0.95), labheat = c("Moderate (%75)", "Extreme (%95)"))

#____ 7.1 # RR of delivery by week at moderate and extreme heat
RRWeek <- filter(EEWeek, x %in% c(0.75, 0.95))
RRWeek <- left_join(RRWeek, heat_points, by = "x")%>%mutate(labheat = fct_rev(labheat)) 

RRWeek%>%mutate(Week=as.numeric(Week))%>%
  ggplot(aes(Week, y, ymin =lowy , ymax = highy, fill=labheat)) + 
  geom_ribbon(aes(ymin = lowy, ymax = highy), alpha=0.10, show.legend = F) +
  geom_line(aes(color=labheat), linewidth=0.5) +
  geom_hline(yintercept = 1, linewidth = 0.5, color="gray70", linetype="dashed") +
  scale_x_continuous(breaks=seq(22,42,2))+
  scale_color_manual( name  ="Heat intensity:",values=c("#feb24c","#e31a1c"))+
  scale_fill_manual(  values=c("#feb24c","#e31a1c"))+
  labs(x = "Week of gestation", y = "RR ", title="Heat-related birhts by gestational week (Fake data)")
ggsave("figures/2.EP/RRWeek.png", height = 7, width = 9, bg="white",type = "cairo")

#____ 7.2 # RR of delivery by week at moderate and extreme heat (statistically significance)
RRWeek %>% mutate(sig=factor(ifelse(lowy>1,1,0)))%>%
  ggplot(aes(Week, y, ymin =lowy , ymax = highy, color=sig)) + 
  geom_hline(yintercept = 1, linewidth = 0.5, linetype="dashed") +
  geom_linerange(aes(x = Week, y = y, ymin = lowy, ymax = highy, color=sig), linewidth = .5, show.legend = F) +
  geom_point(aes(color = sig), shape = 21, fill = "white", size = 2,show.legend = F) +
  scale_color_brewer(palette="Set1",direction = -1) +
  labs(x = "Week of gestation", y = "RR ",title="Heat-related Births by gestational week (Fake data).")+
  facet_wrap(labheat~., nrow = 1)
ggsave("figures/2.EP/RRWeek_sig.png", height = 7, width = 9, bg="white",type = "cairo")

WinM<-RRWeek %>% mutate(sig=factor(ifelse(lowy>1,1,0)))%>%arrange(Week)%>%filter(sig==1& x==0.75)
WinE<-RRWeek %>% mutate(sig=factor(ifelse(lowy>1,1,0)))%>%arrange(Week)%>%filter(sig==1& x==0.95)


#____ 7.3 # RR of be born at each Endpoint at moderate and extreme heat 
RREP0 <- filter(EEEPl, x %in% c(0.75, 0.95))
RREP1 <- left_join(RREP0, heat_points, by = "x")%>%mutate(labheat = fct_rev(labheat))

lev <- c("w28","w28_32","w32_37","w37_39","wfullt")
lab <- c("Extreme PTB","Very PTB","Late PTB", "Early ATB","Full ATB") 
RREP <- RREP1%>%mutate(Sub=factor(Sub,levels=lev, labels=lab),PT =as.numeric(str_detect(Sub, "PTB"))) %>%
         mutate(PT=factor(PT, levels=c(0,1)))
RREP%>% 
  ggplot(aes(x=Sub, y, ymin =lowy , ymax = highy, shape = labheat, group = labheat))+ 
  geom_hline(yintercept = 1, size = 0.5, linetype="dashed") +
  geom_linerange(aes(x = Sub, y = y, ymin = lowy, ymax = highy, colour = PT), show.legend = F,
                 position = position_dodge2(width = 0.5, padding=0.3)) +
  geom_point(aes(fill = PT, shape = labheat),colour = "white", size = 2,
             position = position_dodge2(width = 0.5, padding = 0.3)) +
  scale_color_manual(values = c("#008B8B", "#67001F"), guide = "none") +
  scale_fill_manual(values = c("#008B8B", "#67001F"),  guide = "none") +
  scale_shape_manual(values = c(21, 24),labels = c("Moderate (75%)", "Extreme (95%)"),
                     guide = guide_legend(override.aes = list(size = 3, fill = "black")))+
  labs(x = "", y = "RR ", shape = "Heat intensity:",  title="Heat-related births by endpoint (Fake data).") +
  coord_flip(ylim=c(0.9,1.1))         
ggsave("figures/2.EP/RREndpoints.png", height = 7, width = 9, bg="white",type = "cairo")


# 9 # CLEANUP AND SAVE ---------------------------------------------------------------------------------

rm(list = setdiff(ls(), c("dEP","mvallEPl","EEWeek","RRWeek", "EEEPl", "RREP", "metaEP.gof", "pEPl", 
                          "mvallEPg","pEPg_c","pEPg_m","pEPw28","pEPw28_32","pEPw32_37","pEPw37_39")))
save.image("ResEPFake.RData")



#--------------------------------------------------------------------------------------------------------
# END OF SCRIPT
#--------------------------------------------------------------------------------------------------------

