

#summary(dlist1[[1]])
#wfdef<-31
#widef<-28
#faggr2(dlist1[[1]], wfdef=31,widef=28)

faggr2<-function(dataIN,wfdef=36,widef=22){
  
  
  if(widef>wfdef) STOP("the last week should be greater than the first one")
  
	      wmin<-22                  #first viable week (according who)
      	wmax<-42       						#last week to deliver
      
	dataIN<-arrange(dataIN, BDATE)%>%
              select(BDATE,GWEEK)%>% mutate_at("GWEEK", round)%>%
              filter(GWEEK >= wmin, GWEEK <= wmax)



	NDATE<-data.frame(BDATE=seq(min(dataIN$BDATE),max(dataIN$BDATE),by="day"))
   

      # DAILY COUNT OF PT BABIES (<=wfdef,>=vmef) 
  
	dataOUT <- group_by(dataIN, BDATE) %>% summarise(PT=length(BDATE[GWEEK<=wfdef& GWEEK>=widef]),
									  B=length(BDATE))
   

   
      dataOUT<-left_join(NDATE,dataOUT, by="BDATE")
      dataOUT[is.na(dataOUT)==T]<-0

	
      #  M AND W^j BY WEEK
      
	ndays<-nrow(dataOUT)	    				#n days in ts
	nw<-floor(ndays/7)	    				  #n complete weeks in the study
	nwt<-ifelse(ndays-nw*7==0,nw,nw+1)  		#total n of weeks
     
      cat("\n","########", "\n", "\t", "number of days:", ndays,
				 "\n", "\t", "number of births:", sum(dataOUT$B),
				 "\n", "\t", "number of births born at [",widef,",", wfdef+1,") g weeks:", sum(dataOUT$PT),
	"\n","########","\n"
      )

      
	dataOUT$week<-factor(sort(c(rep( 1:nw,7),rep(nwt, ndays-nw*7))),levels=1:nwt)

	dataIN<-left_join(dataIN,dataOUT,by="BDATE")%>% 
                  select(-c(PT,B))%>%
          		mutate(fSEM=factor(GWEEK,levels=widef:wmax, labels=widef:wmax))

	listM<-dataIN %>%
  			group_by(week,.drop=F) %>%
  			group_map(~ {table(.x$fSEM)})

	M  <- as.matrix(as.data.frame(do.call(rbind, listM)))
      
	nsj<-apply(M, 2, sum)
      n<-sum(nsj)   #n? fetuses
	W <-as.matrix(nsj/(n-cumsum(nsj)+nsj),ncol=1)
      
	Wj<-W[rownames(W) %in% widef:wfdef,]


	# Z_i, W_i BY WEEK
	
      k<-wmax-widef+1
      d<-wfdef-widef+1

	if(k>nwt) STOP("the study period is too much short")
      
      fZij<-function(i){
			Mi<-M[i:(i+k-1),]
			Zj<-rep(NA,length=d)
			for(j in 0:(d-1)){
 				kj<-k-j
		      	Zj[j+1]<-sum(diag(Mi[1:kj,(j+1):k]))}
			return(Zj)
			}
      
	listZij<- map(1:(nwt-k),fZij)

      Zij<- t(as.matrix(do.call(cbind, listZij)))
      
      Z<-apply(Zij,1,sum) 
	    W<-Zij%*%Wj/Z       
      W[Z==0]<-NA
      ZW<-data.frame(week=as.character(1:nrow(Zij)), Z=Z,W=W)
  
      dataOUT<-left_join(dataOUT,ZW,by="week")%>%select(-week)


	#RETURN AGGREGATED DATA
	return(dataOUT)}