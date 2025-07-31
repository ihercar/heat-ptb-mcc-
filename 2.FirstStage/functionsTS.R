#---------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------
# Description:Function to fit a conditional Poisson DLNM model for each outcome in each city. Main details:
#   - CB parameters are fixed by default but can be updated for sensitivity analysis
#   - Analysis restricted to summer (from month summ_ini to summ_end)
#   - At percentile scale
#   - Using gnm with strata: year, month
#   - Using the Pregnancy-at-risk approach: control: log(W), offset: Z  
# Requires:
#   - x: data.frame with the time series
#   - mcities0: metadata
# Returns:
#   - A list with several useful results. The main element is 'mcities', which includes
#     the original metadata for the city, descriptive summaries, and the set of 
#     coefficients and vcov for the cross-basis in the fitted model.
#-----------------------------------------------------------

fit_fs <- function(x, nlag = 4, breaksl=1, varper=c(50,90), varfun="ns", mcities=mcities0, prob = c(0.75, 0.95)){
  
  # basic id components
  
  dfid <- x%>%dplyr::select(Cat.Sub.Co.Ci:Ci)%>%summarize(across(everything(), ~ unique(.)))
  
  
  # summer definition
  
  mcitiesi <- filter(mcities, Co.Ci==unique(x$Co.Ci))
  mcitiesi <- data.frame(mcitiesi)
  summ_ini <- mcitiesi$mmmin_s
  summ_end <- mcitiesi$mmmax_s
  
  # select summer months including several days before to don't loose n due to dlnm
  
  if(summ_ini<summ_end) xs <- x[x$mm %in% summ_ini:summ_end |x$mm==summ_ini-1 & x$dd>=max(x$dd[x$mm==summ_ini-1])-nlag,]
  
  if(summ_ini>summ_end) xs <- x[x$mm %in% summ_ini:12 | x$mm %in% 1:summ_end |x$mm==summ_ini-1 & x$dd>=max(x$dd[x$mm==summ_ini-1])-nlag,]
  
  
  # moving absolute temperature to percentile order
  
  Fd <- ecdf(xs$tmean)
  xs$ptmean <- Fd(xs$tmean)
  
  # knots 
  
  kper <- varper/100
  
  # prediction space
  
  ppred <- unique(c((1:99)/100, prob))
  pprede <- quantile(xs$ptmean, probs=ppred, na.rm=T)
  
  
  # general basis
  
  xlag <- 0:nlag
  argvarm <- list(fun=varfun, knots=kper)
  arglagm <- list(fun="strata", breaks=breaksl)
  cb <- crossbasis(xs$ptmean, lag=nlag, argvar=argvarm, arglag=arglagm, group=xs$vgroup)
  bvar <- do.call("onebasis",c(list(x=ppred),attr(cb,"argvar")))
  blag <- do.call("onebasis",c(list(x=xlag),attr(cb,"arglag")))
  argvar <- list(fun=varfun,knots=kper)
  
  # strata generation
  
  xs <- xs[order(xs$date),]
  xs$month <- as.factor(months(xs$date))
  xs$year <- as.factor(xs$vgroup)
  xs$stratum <- as.factor(xs$year:xs$month)
  
  # day of week
  xs$dow <- as.factor(weekdays(xs$date))
  
  
  # dummy for non-zero strata 
  
  ind <- tapply(xs$PT,xs$stratum,sum)[xs$stratum]
  
  # seasonal control and offset 
  
  xs$lZ <- log(xs$Z)
  xs$lW <- log(xs$W)
  
  modelPTi <- gnm(PT~ cb+lW+dow, offset=lZ, family=quasipoisson,na.action="na.exclude", data=xs,subset=ind>0, eliminate=factor(stratum))
  
  
  # overall effects, extracting  coefs and vcovs as data.frame
  
  redi <- crossreduce(cb, modelPTi, at = ppred)
  
  coefi <- t(redi$coef)
  coefi <- as.data.frame(coefi)
  names(coefi) <- paste("coef", 1:dim(coefi)[2],sep="")
  
  vcovi <- t(vechMat(redi$vcov, diag=T))
  vcovi <- as.data.frame(vcovi)
  names(vcovi) <- paste("vcov", 1:dim(vcovi)[2],sep="")
  
  dfcc <- bind_cols(coefi, vcovi)
  
  
  # calculating mmt, restricted to [0.1, 0.5], and its CI
  
  dfcc$mmti <- findmin(cb,model=modelPTi,from=0.01,to=0.5)
  red_mmti <- crossreduce(cb,  modelPTi, at = ppred, cen = dfcc$mmti)
  
  # calculating effects by lag for extreme and moderate heat
  
  redvarHeatm  <- crossreduce(cb, modelPTi, type = "var", value = prob[1], cen = dfcc$mmti)
  redvarHeatE  <- crossreduce(cb, modelPTi, type = "var",  value = prob[2], cen = dfcc$mmti)
  
  # calculating some gof measures
  
  phii <- summary(modelPTi)$dispersion
  LogLi <- sum(dpois(modelPTi$y, modelPTi$fitted.values, log=TRUE))
  dfi <- summary(modelPTi)$df[3]
  QAICi <- -2*LogLi/phii + 2*dfi
  
  gofi <- data.frame(phi=phii,LogL=LogLi, df=dfi, QAIC=QAICi)
  
  # extending mcitiesi with  additional descriptives, coefs, vcov, and mmt
  
  dfdesc <- x%>%summarize(mintmean=min(tmean, na.rm=T), maxtmean=max(tmean, na.rm=T),
                        avgtmean=mean(tmean, na.rm=T),rgtmean=maxtmean-mintmean, 
                        nB=sum(B), nPTB=sum(PT))
  
  dfdesc_s <- xs%>%summarize(mintmean_s=min(tmean, na.rm=T), maxtmean_s=max(tmean, na.rm=T),
                           avgtmean_s=mean(tmean, na.rm=T),rgtmean_s=maxtmean_s-mintmean_s, 
                           nB_s=sum(B), nPTB_s=sum(PT))
  
  
  mcitiesi <- bind_cols(dfid, mcitiesi, dfdesc,dfdesc_s,dfcc)
  
  
  #list to be return
  
  output <- list(mod = modelPTi, red_mmt = red_mmti, redlagHeatE = redvarHeatE, redlagHeatM = redvarHeatE,
                 gof=gofi, mcities=mcitiesi)
  
  return(output)
  
}
#---------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------
# Function: fpredmv
# Description:
#   Performs detailed predictions from a mixmeta model's coefficients and variance-covariance matrix.
#   Calculates relative risks (RR) with confidence intervals centered at a chosen percentile.
#   Supports centering at the minimum risk percentile ("m") or a fixed percentile ("f").
#   Relies on global variables bvar (basis), ppred (prediction points), cen (centering method), 
#   and cenf (fixed centile).
# Requires:
#   - EEi: A list containing 'fit' (coefficients) and 'vcov' (variance-covariance matrix) from the meta-model.
# Returns:
#   - A dataframe with prediction points, relative risks, confidence intervals, centering information, and other metrics.
#-----------------------------------------------------------

fpredmv<-function(EEi){
  predi <- crosspred(bvar,coef= EEi$fit,vcov= EEi$vcov,at=ppred,model.link="log",cen=0.01)
  mmti<-findmin(bvar, coef=EEi$fit, vcov=EEi$vcov, by=0.01,from=0.01, to=0.5)
  if(cen=="m"){ceni<-mmti
  predicen <- crosspred(bvar,coef=EEi$fit,vcov=EEi$vcov,at=ppred,model.link="log",by=0.01, cen=ceni)}
  if(cen=="f"){ceni<-cenf}
  predicen <- crosspred(bvar,coef=EEi$fit,vcov=EEi$vcov,at=ppred,model.link="log",by=0.01, cen=ceni)  #centrado para el dibujo
  ch =(exp(predicen$allfit)-1)*100
  lowch   =(exp((predicen$allfit-1.96*predicen$allse))-1)*100
  highch   =(exp((predicen$allfit+1.96*predicen$allse))-1)*100
  predi<-data.frame(x=predicen$predvar, 
                    y= predicen$allRRfit,lowy=predicen$allRRlow,highy=predicen$allRRhigh,
                    ch=(exp(predicen$allfit)-1)*100, lowch=lowch, highch=highch, 
                    mmt=mmti, cen=ceni)
  return(predi)
}
#---------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------
# function to start the plot of curves for Heat-PTB association using ggplot2, with consistent styling
#-----------------------------------------------------------

plot_curve <- function(data) {
  p <- ggplot(data, aes(x, y)) +
    geom_hline(yintercept = 1, linewidth = 0.8, linetype = "dashed") +
    geom_vline(xintercept = 0.75, color = "#feb24c", alpha = 0.7, linewidth = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.95, color = "#e31a1c", alpha = 0.7, linewidth = 0.5, linetype = "dashed") +
    annotate("point", x = 0.75, y = 0.98, shape = 21, size = 2, color = "#feb24c", fill = "#feb24c") +
    annotate("point", x = 0.95, y = 0.98, shape = 24, size = 2, color = "#e31a1c", fill = "#e31a1c") +
    geom_ribbon(aes(ymin = lowy, ymax = highy), fill = "grey85", alpha = 0.3) +
    geom_line(colour = "#008B8B", linewidth = 1) +
    scale_x_continuous(breaks = seq(0, 1, 0.25),labels = ~number(., scale = 100))+
    labs(x = "Percentile of daily mean temperature", y = "RR")+
    theme_get()
  return(p)}

#-----------------------------------------------------------------------
#-----------------------------------------------------------
# function to start the RR plot for Moderate and Extreme heat using ggplot2 with consistent styling
#-----------------------------------------------------------

plot_RR <- function(data) {
  p <- ggplot(data,aes(y = fct_reorder(Sub, y), x = y, xmin = lowy, xmax = highy, shape = labheat)) +
    geom_vline(xintercept = 1, linewidth = 0.8) +
    geom_linerange(aes(colour = labheat), linewidth = 0.7, position = position_dodge2(width = 0.7, padding = 0.3)) +
    geom_point(aes(fill = labheat), colour = "white", size = 2.2, position = position_dodge2(width = 0.7, padding = 0.3)) +
    guides(color = "none") +
    scale_color_manual(values = c("#feb24c", "#e31a1c")) +
    scale_fill_manual(values = c("#feb24c", "#e31a1c")) +
    scale_shape_manual(values = c(21, 24)) +
    scale_x_continuous(breaks = pretty_breaks(n = 7)) +
    labs(y = NULL, x = "RR",shape = "Heat intensity:", colour = "Heat intensity:", fill = "Heat intensity:") 
  theme_get()
  return(p)}

#---------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------
# FUNCTION TO ESTIMATE MINIMUM OF A EXPOSURE-RESPONSE FUNCTION FROM A FITTED MODEL
# Tobias, Armstrong, Gasparrini. Epidemiol, 2017. 
#-----------------------------------------------------------

findmin <- function(basis,model=NULL,coef=NULL,vcov=NULL,at=NULL,from=NULL,to=NULL,by=NULL,sim=FALSE,nsim=5000) {
  
  #-----------------------------------------------------------
  # ARGUMENTS:
  # - basis: A SPLINE OR OTHER BASIS FOR AN EXPOSURE x CREATED BY DLNM FUNCTION CROSSBASIS OR ONEBASIS
  # - model: THE FITTED MODEL
  # - coef AND vcov: COEF AND VCOV FOR basis IF model IS NOT PROVIDED 
  # - at: A NUMERIC VECTOR OF x VALUES OVER WHICH THE MINIMUM IS SOUGHT 
  # OR
  # - from, to: RANGE OF x VALUES OVER WHICH THE MINIMUM IS SOUGHT.
  # - by: INCREMENT OF THE SEQUENCES x VALUES OVER WHICH THE MINIMUM IS SOUGHT 
  # - sim: IF BOOTSTRAP SIMULATION SAMPLES SHOULD BE RETURNED
  # - nsim: NUMBER OF SIMULATION SAMPLES 
  #-----------------------------------------------------------

   # CREATE THE BASIS AND EXTRACT COEF-VCOV #
    # CHECK AND DEFINE BASIS 
    if(!any(class(basis)%in%c("crossbasis","onebasis")))
      stop("the first argument must be an object of class 'crossbasis' or 'onebasis'") 
  #
    # INFO
    one <- any(class(basis)%in%c("onebasis"))
    attr <- attributes(basis)
    range <- attr(basis,"range")
    if(is.null(by)) by <- 0.1
    lag <- if(one) c(0,0) else cb=attr(basis,"lag") 
    if(is.null(model)&&(is.null(coef)||is.null(vcov)))
      stop("At least 'model' or 'coef'-'vcov' must be provided")
    name <- deparse(substitute(basis))
    cond <- if(one) paste(name,"[[:print:]]*b[0-9]{1,2}",sep="") else
      paste(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}",sep="") 
  #
    # SET COEF, VCOV CLASS AND LINK
    if(!is.null(model)) {
      model.class <- class(model)
      coef <- dlnm:::getcoef(model,model.class)
      ind <- grep(cond,names(coef))
      coef <- coef[ind]
      vcov <- dlnm:::getvcov(model,model.class)[ind,ind,drop=FALSE] 
      model.link <- dlnm:::getlink(model,model.class)
  } else model.class <- NA
  #
    # CHECK
    if(length(coef)!=ncol(basis) || length(coef)!=dim(vcov)[1] || any(is.na(coef)) || any(is.na(vcov)))
      stop("model or coef/vcov not consistent with basis")
  #
    # DEFINE at
  at <- dlnm:::mkat(at,from,to,by,range,lag,bylag=1) 
  predvar <- if(is.matrix(at)) rownames(at) else at 
  predlag <- dlnm:::seqlag(lag,by=1)
  #
    # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON TYPE) 
    type <- if(one) "one" else "cb"
    Xpred <- dlnm:::mkXpred(type,basis,at,predvar,predlag,cen=NULL)
    Xpredall <- 0
    for(i in seq(length(predlag))) {
      ind <- seq(length(predvar))+length(predvar)*(i-1)
      Xpredall <- Xpredall + Xpred[ind,,drop=FALSE] 
    }
  #
   # FIND THE MINIMUM
  #
    pred <- drop(Xpredall%*%coef) 
    ind <- which.min(pred)
    min <- predvar[ind]
  # 
   # SIMULATIONS
    #
    if(sim) {
      # SIMULATE COEFFICIENTS
      k <- length(coef)
      eigen <- eigen(vcov)
      X <- matrix(rnorm(length(coef)*nsim),nsim)
      coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X) 
      # COMPUTE MINIMUM
      minsim <- apply(coefsim,2,function(coefi) { 
        pred <- drop(Xpredall%*%coefi)
        ind <- which.min(pred) 
        return(predvar[ind])
      })
    }
  # 
   res <- if(sim) minsim else min
  #
    return(res)
  }

#--------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------
# Functions to calculate p-Wald for covs in mveta, by default test the complete effect when an 
# interaction is in the model                                                       
#-----------------------------------------------------------

# WALD TEST by default test the complete effect when an interaction is in the model
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)))
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}


#--------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------
# Functions to calculate p-Wald for covs in mveta, excluding the interaction in case is present                                                       
#-----------------------------------------------------------

fwald2 <- function(model,var, discount) {
  ind   <- grep(var,names(coef(model)))
  indno <- grep(discount,names(coef(model)))
  indsi<-setdiff(ind, indno)
  #names(coef(model))[indsi]
  coef <- coef(model)[indsi]
  vcov <- vcov(model)[indsi,indsi]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}

#--------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------
# function to compute AN to heat. Data list is .y and blup list is .x. 
# Requires the function attrdlr
#-----------------------------------------------------------

funANsim <- function(.x,.y, nsim = 2000){
  
  coefallx <- as.numeric(.x[grep("coef", names(.x))])
  mat<-matrix(0, nrow=length(coefallx), ncol=length(coefallx))
  mat[lower.tri(mat, diag = TRUE)] <- as.numeric(.x[grep("vcov", names(.x))])
  vcovallx<-mat+t(mat)-diag(diag(mat))
  colnames(vcovallx)<-rownames(vcovallx)<-colnames(coefallx)
  
  
  if(cen=="m") {mmti<-findmin(bvar, coef=coefallx, vcov=vcovallx, by=0.01,from=0.01, to=0.5)
  }
  if(cen=="f") {mmti[j]<-0.01} # get descriptives
  N<-length(.y$ptmean)
  heat<-quantile(.y$tmean, probs=mmti, na.rm=T)
  Nheat<- length(.y$ptmean > mmti)
  totPT<-sum(.y$PT, na.rm=T)
  totB<-sum(.y$B, na.rm=T)
  
  cb <- crossbasis(.y$ptmean, lag=nlag, argvar=argvarm, arglag=arglagm, group=.y$vgroup)
  expi<-.y$ptmean
  bvari <- do.call("onebasis",c(list(x=expi),attr(cb,"argvar")))
  ran<-c(mmti,max(expi,na.rm=T))
  respi <- .y$PT
  BMAyi <-  zoo::rollmean(respi, nlag, fill = NA)
  attri<-attrdlr(expi, cb, respi, coef=coefallx, vcov=vcovallx,
                 type="an", dir="for", tot=TRUE, 
                 cen=mmti, range=ran,sim=T, nsim=1000)
  AF<-100*quantile(attri, c(0.025,0.5,0.975))/totPT
  res<-data.frame(N=N, heat=heat,Nheat=Nheat,totB=totB, totPT=totPT, mmti=mmti, lowAF=AF[1], AF=AF[2], hiAF=AF[3] )
  return(list(AFres=res, ANsim=attri))
}

#---------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------
# FUNCTION FOR COMPUTING ATTRIBUTABLE MEASURES FROM DLNM
# ADAPTED FROM FUNCTION DEVELOPED BY Antonio Gasparrini 2014
# https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata/blob/master/attrdl.R
#-----------------------------------------------------------

attrdlr <- function(x,basis,cases,model=NULL,coef=NULL,vcov=NULL,type="af",
                    dir="back",tot=TRUE,cen,range=NULL,sim=FALSE,nsim=5000) {
  
  #
  # CHECK VERSION OF THE DLNM PACKAGE
  if(packageVersion("dlnm")<"2.2.0") 
    stop("update dlnm package to version >= 2.2.0")
  #
  # EXTRACT NAME AND CHECK type AND dir
  name <- deparse(substitute(basis))
  type <- match.arg(type,c("an","af"))
  dir <- match.arg(dir,c("back","forw"))
  #
  # DEFINE CENTERING
  if(missing(cen) && is.null(cen <- attr(basis,"argvar")$cen))
    stop("'cen' must be provided")
  if(!is.numeric(cen) && length(cen)>1L) stop("'cen' must be a numeric scalar")
  attributes(basis)$argvar$cen <- NULL
  #  
  # SELECT RANGE (FORCE TO CENTERING VALUE OTHERWISE, MEANING NULL RISK)
  if(!is.null(range)) x[x<range[1]|x>range[2]] <- cen
  #
  # COMPUTE THE MATRIX OF
  #   - LAGGED EXPOSURES IF dir="back"
  #   - CONSTANT EXPOSURES ALONG LAGS IF dir="forw"
  lag <- attr(basis,"lag")
  if(NCOL(x)==1L) {
    at <- if(dir=="back") tsModel:::Lag(x,seq(lag[1],lag[2])) else 
      matrix(rep(x,diff(lag)+1),length(x))
  } else {
    if(dir=="forw") stop("'x' must be a vector when dir='forw'")
    if(ncol(at <- x)!=diff(lag)+1) 
      stop("dimension of 'x' not compatible with 'basis'")
  }
  #
  # NUMBER USED FOR THE CONTRIBUTION AT EACH TIME IN FORWARD TYPE
  #   - IF cases PROVIDED AS A MATRIX, TAKE THE ROW AVERAGE
  #   - IF PROVIDED AS A TIME SERIES, COMPUTE THE FORWARD MOVING AVERAGE
  #   - THIS EXCLUDES MISSING ACCORDINGLY
  # ALSO COMPUTE THE DENOMINATOR TO BE USED BELOW
  if(NROW(cases)!=NROW(at)) stop("'x' and 'cases' not consistent")
  if(NCOL(cases)>1L) {
    if(dir=="back") stop("'cases' must be a vector if dir='back'")
    if(ncol(cases)!=diff(lag)+1) stop("dimension of 'cases' not compatible")
    den <- sum(rowMeans(cases,na.rm=TRUE),na.rm=TRUE)
    cases <- rowMeans(cases)
  } else {
    den <- sum(cases,na.rm=TRUE) 
    if(dir=="forw") 
      cases <- rowMeans(as.matrix(tsModel:::Lag(cases,-seq(lag[1],lag[2]))))
  }
  #
  
  #
  # EXTRACT COEF AND VCOV IF MODEL IS PROVIDED
  if(!is.null(model)) {
    cond <- paste0(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}")
    if(ncol(basis)==1L) cond <- name
    model.class <- class(model)
    coef <- dlnm:::getcoef(model,model.class)
    ind <- grep(cond,names(coef))
    coef <- coef[ind]
    vcov <- dlnm:::getvcov(model,model.class)[ind,ind,drop=FALSE]
    model.link <- dlnm:::getlink(model,model.class)
    if(model.link!="log") stop("'model' must have a log link function")
  }
  #
  # IF REDUCED ESTIMATES ARE PROVIDED
  typebasis <- ifelse(length(coef)!=ncol(basis),"one","cb")
  #
  
  # PREPARE THE ARGUMENTS FOR TH BASIS TRANSFORMATION
  predvar <- if(typebasis=="one") x else seq(NROW(at))
  predlag <- if(typebasis=="one") 0 else dlnm:::seqlag(lag)
  #  
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON typebasis)
  if(typebasis=="cb") {
    Xpred <- dlnm:::mkXpred(typebasis,basis,at,predvar,predlag,cen)
    Xpredall <- 0
    for (i in seq(length(predlag))) {
      ind <- seq(length(predvar))+length(predvar)*(i-1)
      Xpredall <- Xpredall + Xpred[ind,,drop=FALSE]
    }
  } else {
    basis <- do.call(onebasis,c(list(x=x),attr(basis,"argvar")))
    Xpredall <- dlnm:::mkXpred(typebasis,basis,x,predvar,predlag,cen)
  }
  #  
  # CHECK DIMENSIONS  
  if(length(coef)!=ncol(Xpredall))
    stop("arguments 'basis' do not match 'model' or 'coef'-'vcov'")
  if(any(dim(vcov)!=c(length(coef),length(coef)))) 
    stop("arguments 'coef' and 'vcov' do no match")
  if(typebasis=="one" && dir=="back")
    stop("only dir='forw' allowed for reduced estimates")
  #
  
  # COMPUTE AF AND AN 
  af <- 1-exp(-drop(as.matrix(Xpredall%*%coef)))
  an <- af*cases
  #
  # TOTAL
  #   - SELECT NON-MISSING OBS CONTRIBUTING TO COMPUTATION
  #   - DERIVE TOTAL AF
  #   - COMPUTE TOTAL AN WITH ADJUSTED DENOMINATOR (OBSERVED TOTAL NUMBER)
  if(tot) {
    isna <- is.na(an)
    af <- sum(an[!isna])/sum(cases[!isna])
    an <- af*sum(cases[!isna])  
  }
  #
  
  # EMPIRICAL CONFIDENCE INTERVALS
  if(!tot && sim) {
    sim <- FALSE
    warning("simulation samples only returned for tot=T")
  }
  if(sim) {
    # SAMPLE COEF
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    # RUN THE LOOP
    # pre_afsim <- (1 - exp(- Xpredall %*% coefsim)) * cases # a matrix
    # afsim <- colSums(pre_afsim,na.rm=TRUE) / sum(cases[!isna],na.rm=TRUE)
    afsim <- apply(coefsim,2, function(coefi) {
      ani <- (1-exp(-drop(Xpredall%*%coefi)))*cases
      sum(ani[!is.na(ani)])/sum(cases[!is.na(ani)])
    })
    ansim <- afsim*den
  }
  #
  
  #
  res <- if(sim) {
    if(type=="an") ansim else afsim
  } else {
    if(type=="an") an else af    
  }
  #
  return(res)
}


#--------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------
# gglot themes

theme_cc<-function(){ 
  theme_minimal() %+replace%
    theme(
      plot.margin = margin(10, 10, 10, 10),
      panel.border = element_rect(colour= "black", linewidth= 1,fill = NA,linetype=1),
      panel.grid=element_line(color="white"),
      plot.title = element_text(hjust = 0.5, size=12, face="bold"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8, vjust = 0.5),
      axis.title=element_text(hjust = 0.5, size=9),
      panel.background = element_blank(), 
      panel.spacing = unit(1, "lines"),
      legend.position="bottom"
    )}


theme_cc2<-function(){ 
  theme_bw(base_size = 12) +
    theme(
      panel.grid=element_line(color="white"),
      strip.text = element_text(hjust = 0),
      legend.position = "bottom",
      legend.justification = "left",
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size=13, face="bold")
    )}    
#-----------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
