
#------------------------------------------------------------------------------------------------------
# MCC-HEATPTB AF and TAF CALCULATION (FAKE DATA)
#------------------------------------------------------------------------------------------------------
# To calculate  fractions of births attributable to heat
# Core functions: funANsim --see details in functionsTS.R.
# Core data: 1) TS Data, e.g. from 1. FirstStage -> FSFake.RData. 
#            2) Results from 2.SecondStage->ResPTBFake.RData
#------------------------------------------------------------------------------------------------------


# 0 # CLEAN ENVIRONMENT AND LOAD PACKAGES -------------------------------------------------------------
rm(list = ls())
library(pacman)
pacman::p_load(tidyverse, splines, mgcv, dlnm, gnm, mixmeta, ggplot2, ggthemes, ggtext, patchwork, scales,
               RColorBrewer,flextable, officer)


# 1 # CHARGE FUNCTIONS AND PARAMETERS ------------------------------------------------------------------
source("functionsTS.R")
source("parameters.R")



# 2 # LOAD DATA AND SECOND STAGE RESULTS ---------------------------------------------------------------
load("FSFake.RData")
load("ResPTBFake.RData")



# 3 # GET BLUPs AND ORGANIZE BY COUNTRY------------------------------------------------------------------
blup <- blup_mvallPTB

nCo <- unique(dPTB$Co)  
lisCo <- split(nCo, factor(nCo))
lisblupCo <- map(lisCo, function(z){ CoCii<-dPTB$Co.Ci[dPTB$Co==z]
              blupCoi<-blup[CoCii]
              l1<-map_depth(blupCoi,1, function(x){
                     coef<-x$blup
                     vcov<-vechMat(x$vcov, diag=T)
                     names(vcov)<-paste("vcov", 1:length(vcov),sep="")
                     matcoef<-data.frame(matrix(c(coef, vcov), nrow=1))
                     names(matcoef)<-c(names(coef), names(vcov))
                     return(matcoef)}) 
              return(l1)})

# 4 # GET PTB and ptmean time series BY COUNTRY
nCo <- unique(dPTB$Co)  
lisCo <- split(nCo, factor(nCo))
lisdPTBCo <- map(lisCo, function(z){ dPTBCoi<-dPTB[dPTB$Co==z,]
              l2<-split(dPTBCoi$Co.Ci, factor(dPTBCoi$Co.Ci))
              l3<-map(l2, function(x){ namx<-paste("EP.w37.",x, sep="")
                     fs_resij<-fs_res[[namx]]
                     dij<-fs_resij$mod$data
                     return(dij)})
              return(l3)})


# 5 # CALCULATE AN (eCI) IN EACH CITY BY COUNTRY, FUNCTION funANsim 
attrsubsep <- map(lisCo, function(x) { 
                       listBlupsubco<- lisblupCo[[x]]
                       listDBsubco  <- lisdPTBCo[[x]]
                       listDBsubco<-listDBsubco[names(listBlupsubco)] 
                       map2(listBlupsubco, listDBsubco, funANsim)})    



# 6 # CALCULATE AF and TAF (eCI) FROM AN, by City, Country and Overall 

#____ 6.1 # By City 
ressub_city <- map_depth(attrsubsep, 2, ~ .x$AFres)%>% imap(~ bind_rows(.x, .id = "city")) %>%
                 bind_rows(.id = "country")
ressub_city <- ressub_city%>%mutate(TAF=(AF*totPT*10^4)/totB,  
                                    lowTAF=(lowAF*totPT*10^4)/totB,  
                                    hiTAF= (hiAF*totPT*10^4)/totB)  #AF in %, thus TAF is rate by 10^6
rownames(ressub_city)<-NULL


#____ 6.2 # By Country 
ressub_country <- map(attrsubsep, function(x){
                    AN_sum <- reduce(map_depth(x, 1, ~ .x$ANsim),`+`)
                    Matres_ci<-map_depth(x, 1, ~ .x$AFres)%>%bind_rows(.id="city")
                    Matres_co<- summarize(Matres_ci, 
                        N=sum(N),
                        heat=mean(heat),
                        totB=sum( totB),
                        totPT=sum(totPT),
                        lowAF=100*quantile(AN_sum,0.025)/totPT,
                        AF=100*quantile(AN_sum,0.5)/totPT,
                        hiAF=100*quantile(AN_sum,0.975)/totPT,
                        TAF=(AF*totPT*10^4)/totB,  
                        lowTAF=(lowAF*totPT*10^4)/totB,  
                        hiTAF= (hiAF*totPT*10^4)/totB)
                     return(Matres_co)})%>%bind_rows(.id="country")


#____ 6.3 # Overall
AN_sum_co <- map_depth(attrsubsep, 1, function(x){AN_sum<-reduce(map_depth(x, 1, ~ .x$ANsim),`+`)
                                                return(AN_sum)})
AN_sum <- reduce(map_depth(AN_sum_co, 1, ~ .x),`+`)
ressub_all <- summarize(ressub_city,  
                        N=sum(N),
                        heat=mean(heat),
                        totB=sum( totB),
                        totPT=sum(totPT),
                        lowAF=100*quantile(AN_sum,0.025)/totPT,
                        AF=100*quantile(AN_sum,0.5)/totPT,
                        hiAF=100*quantile(AN_sum,0.975)/totPT,
                        TAF=(AF*totPT*10^4)/totB,  
                        lowTAF=(lowAF*totPT*10^4)/totB,  
                        hiTAF= (hiAF*totPT*10^4)/totB )
  
TAFw37sep <- list(AF_city=ressub_city, AF_country=ressub_country, AF_all=ressub_all)



# 7 # TABULATING ----------------------------------------------------------------------------------------


#____ 7.1 # function get_color
color_palette <- brewer.pal(7,"YlOrRd")
get_color <- function(values) {
  cut(values, 
      breaks = quantile(values, probs = seq(0, 1, length.out = 8), na.rm = TRUE), 
      labels = color_palette, 
      include.lowest = TRUE)
}

#____ 7.2 # Prepare the Table
TAFall <- TAFw37sep$AF_all %>% mutate(
             eAF  = paste0(round(AF, 2), " (95%CI: ", round(lowAF, 2), ", ", round(hiAF, 2), ")"),
             eTAF = paste0(round(TAF),   " (95%CI: ", round(lowTAF),   ", ", round(hiTAF),   ")"))

tTAF <- TAFw37sep$AF_country %>% arrange(desc(TAF)) %>% bind_rows(
             TAFw37sep$AF_all %>% mutate(country = "All") %>% select(names(.)))

TAF1 <- tTAF%>%mutate(Country=country,
                    AFIC=paste( round(AF,2)," (", round(lowAF,2),", ",round(hiAF,2),") ", sep=""),
                    Color_AF="",
                    TAFIC=paste( round(TAF)," (", round(lowTAF),", ",round(hiTAF),") ", sep=""),
                    Color_RAF="")%>% dplyr::select(Country:Color_RAF)

#____ 7.3 # Flextable
Color_AF <- as.character(get_color(tTAF$AF))
Color_RAF <- as.character(get_color(tTAF$TAF))

ft <- qflextable(TAF1) %>%
  align(align = "center", part = "all") %>%                      
  bg(j = 3, bg = Color_AF,  part = "body") %>%                   
  bg(j = 5, bg = Color_RAF, part = "body") %>%                   
  width(j = 3, width = 0.07) %>%                                 
  width(j = 5, width = 0.07) %>%                                 
  hline(i = 5, border = fp_border(color = "#666666"), part = "body") %>%  
  add_footer_lines(
    "AF: PTB fraction attributable to heat(%), RAF: PTB Rate (per million births) attributable to heat. 
    Color scales classify estimates into six increasing-magnitude categories, from light to dark red"
  ) %>%                                                          
  color(part = "footer", color = "#666666") %>%                  
  set_header_labels(
    Country = "Country", AFIC = "AF (CI95%)", Color_AF = "", 
    TAFIC = "RAF (CI95%)", Color_RAF = ""
  )
ft



# 8 # CLEANUP AND SAVE ----------------------------------------------------------------------------------


rm(list=setdiff(ls(), c("lisblupCo","lisdPTBCo","TAFw37sep","ft")))
save.image("ResTAFFake.RData")



#--------------------------------------------------------------------------------------------------------
# END OF SCRIPT
#--------------------------------------------------------------------------------------------------------

