#' Population reconstruction code
#' 
#' Computes popn statistics and optionally plots smoothed net change, pup counts, total population size, non-pup population size over time 
#' Some of these functions are not used in the current OSP analysis.
#' 
#' @param pup_counts vector of pup counts
#' @param model survival model object from marked
#' @param Phiformula formula for survival model
#' @param Phidata dataframe with values for survival predictions
#' @param Phi_coef fixed effect coefficients for Phi
#' @param random if TRUE, add random component for survival calculation in dfu$u
#' @param df portion of design data for survival for which calculatoins are to be performed with Phidata values
#' @param dfu matrix with normal random variate for each prediction
#' @param Phi_sigma vector of estimates for log(sigma)
#' @param begin year to initiate population projection
#' @param plot if TRUE, project_pop produces plots; otherwise, only returns stats
#' @param r assumed rate of increase at beginning of time series to create initial stable age distribution
#' @param nrep number of replicate simulations
#' @return vector containing year of maximum smoothed net increase, smoothed population size at MNPL and current smoothed N, all for non-pups.
#' @author Jeff Laake
#' @aliases popn_reconstruction popn_reconstruction_error population_reconstruction make_predictions make_pup_counts population_reconstruction_error
#' @export popn_reconstruction popn_reconstruction_error population_reconstruction make_predictions make_pup_counts population_reconstruction_error
popn_reconstruction=function(pup_counts,model,Phiformula,Phidata,Phi_coef=NULL,random=FALSE,marked=FALSE,df=NULL,dfu=NULL,
		begin=1975,plot=TRUE,r=0.41,Phi_sigma=NULL)
{
	if(marked)
	{
		Phidm=model.matrix(Phiformula,Phidata)
		if(is.null(Phi_coef))Phi_coef=model$results$coeflist$phi_beta
#   create estimates
		if(random)
			survivals=plogis(Phidm%*%Phi_coef+dfu%*%exp(Phi_sigma))
		else
			survivals=plogis(Phidm%*%Phi_coef)
	} else {
		if(is.null(df))stop("\nMust specify df for RMark Burnham model")
		Phidata=Phidata[order(Phidata$sex,Phidata$Time,Phidata$Age),]
		if(is.null(Phi_coef))
			survivals=predict_real(model,df,"S",data=Phidata,beta=Phi_coef)$estimate^(1/Phidata$time.intervals)
		else
			survivals=predict_real(model,df,"S",data=Phidata,beta=Phi_coef,se=FALSE,vcv=FALSE)$estimate^(1/Phidata$time.intervals)	
	}
	Sf=t(tapply(survivals[Phidata$sex=="F"],list(Phidata$Time[Phidata$sex=="F"],Phidata$Age[Phidata$sex=="F"]),sum))
	Sm=t(tapply(survivals[Phidata$sex=="M"],list(Phidata$Time[Phidata$sex=="M"],Phidata$Age[Phidata$sex=="M"]),sum))
# take 11/12ths of the pup survival to project to age 1.
	Sm[1,]=Sm[1,]^(11/12)
	Sf[1,]=Sf[1,]^(11/12)	
	Sm=Sm[-nrow(Sm),]
	Sf=Sf[-nrow(Sf),]
# project population
	results=project_pop(pup_counts,Sf, Sm, r=r,begin=1975,plot=plot)
	return(results)
}

popn_reconstruction_error=function(nrep,pup_counts,model,Phiformula,Phidata,marked=FALSE,df=NULL,random=FALSE,dfu=NULL,begin=1975,r=0.041,Phi_sigma=NULL)
{
	if(random)stop("Cannot use random effect models at present")
	if(marked)
	{
#       extract the variance-covariance matrix for the betas
		vcv = model$results$beta.vcv
		if(is.null(vcv)|all(is.na(vcv)))stop("Missing vcv for beta")
		vcv_names = rownames(model$results$beta.vcv)
#       extract the portion of the v-c matrix for the Phi parameters
		Phi_vcv = vcv[grep("Phi", vcv_names),grep("Phi", vcv_names)]
	}
# compute nrep multivariate normal variables with means Phi_coef
# and variance-covariance matrix Phi_vcv
	if(marked)
	{
		Phi_coef=model$results$coeflist$phi_beta
	}else
	{
		which.beta=grep("S:",rownames(model$results$beta))
		Phi_coef=model$results$beta$estimate[which.beta]
		Phi_vcv= model$results$beta.vcv[which.beta,which.beta]
	}
	Phi_rep = rmvnorm(nrep, Phi_coef, Phi_vcv)
# loop over each and project pop and save summary statistics
	mnpl_t=NULL
	for (i in 1:nrep)
	{
		if(!is.list(pup_counts))
		{
			stats=popn_reconstruction(pup_counts,model,Phiformula,Phidata,Phi_coef=Phi_rep[i,],df=df,random=FALSE,marked=marked,begin=begin,plot=FALSE,r=r,Phi_sigma=Phi_sigma)$stats
		} else {
			stats=popn_reconstruction(pup_counts[[i]],model,Phiformula,Phidata,Phi_coef=Phi_rep[i,],df=df,random=FALSE,marked=marked,begin=begin,plot=FALSE,r=r,Phi_sigma=Phi_sigma)$stats
		}
		mnpl_t=rbind(mnpl_t,stats)
	}
	rownames(mnpl_t)=NULL
	return(mnpl_t)
}	
population_reconstruction=function(predictions,reml.formula,key.estimates,pup_counts,r=0.041)
{
	z=predictions$real
	z$py=z$pup+z$yearling
	z=droplevels(z)
	dm=model.matrix(reml.formula,z)
	mod.reml=var.components.reml(z$real,design=dm,predictions$vcv, rdesign=model.matrix(~-1+py:time,z))
	
# get estimates of reals for entire 1975-2014 time series
	data(Phidata)
	Phidata$Time=Phidata$Time+2
	Phidata=Phidata[order(Phidata$sex,Phidata$Time,Phidata$Age),]
	
	Phidata$time=factor(cut(Phidata$Time+1987,1986:2015,include.lowest=TRUE))
	Phidata$pup=ifelse(Phidata$Age==0,1,0)
	Phidata$yearling=ifelse(Phidata$Age==1,1,0)
	Phidata$twothree=ifelse(Phidata$Age%in%2:3, 1,0)
	Phidata$threeplus=ifelse(Phidata$Age>=3,1,0)
	Phidata$fourplus=ifelse(Phidata$Age>=4,1,0)
	
# indicator variables for <1993 (pre93) and >=1993 (p93plus); >=1989 post (excludes
# 1989 if non-pup).
	Phidata$p87plus=ifelse(Phidata$Time>=0&Phidata$time!="1999",1,0)
	Phidata$p88plus=ifelse(Phidata$Time>=1&Phidata$time!="1999",1,0)
	Phidata$p89plus=ifelse(Phidata$Time>=2&Phidata$time!="1999",1,0)
	Phidata$p90plus=ifelse(Phidata$Time>=3&Phidata$time!="1999",1,0)
	Phidata$p91plus=ifelse(Phidata$Time>=4&Phidata$time!="1999",1,0)
	Phidata$p94plus=ifelse(Phidata$Time>=7,1,0)
	Phidata$timebin=cut(Phidata$Time,c(0,4,6:28),include.lowest=TRUE,labels=c("1988-1991","1992-1993",1994:2015))
	Phidata$plus4time=cut(Phidata$Time,c(0,7:28),include.lowest=TRUE,labels=c("1988-1994",1995:2015))
	Phidata$AgeS=cut(Phidata$Age,c(0,4,7,10,13,17,30),right=FALSE)
	Phidata$age=factor(Phidata$Age)
	Phidata$key=paste(Phidata$sex,Phidata$Time,Phidata$Age)
	
	which.reals=Phidata$key%in%key.estimates
	
	dm=model.matrix(reml.formula,Phidata)
	Phidata=cbind(Phidata,S=dm%*%mod.reml$beta$Estimate)
# maybe add process error to pup/yearling estimates and iid to all
# fix any to 1 that exceeded 1
	Phidata$S[Phidata$S>1]=1
	Phidata$S[which.reals]=predictions$real$real
# convert to Sf, Sm matrices
	Sf=matrix(Phidata$S[Phidata$sex=="F"],ncol=40,nrow=25)
	Sm=matrix(Phidata$S[Phidata$sex=="M"],ncol=40,nrow=25)
	Sf[1,]=Sf[1,]^(11/12)
	Sm[1,]=Sm[1,]^(11/12)
	Sm=Sm[-25,]
	Sf=Sf[-25,]
	results=project_pop(pup_counts[-41],Sf[,-40],Sm[,-40],r=r)
	return(results)
}
make_predictions=function(model,zc.ddl,beta=NULL,vcv=TRUE)
{
	data(Phidata)
	df=zc.ddl$S[zc.ddl$S$wt=="(-13,-4.5]",]
	Phidata$Time=Phidata$Time+2
	Phidata=Phidata[Phidata$Time>=0,]
	Phidata=Phidata[order(Phidata$sex,Phidata$Time,Phidata$Age),]
	
	Phidata$key=paste(Phidata$sex,Phidata$Time,Phidata$Age)
	df$key=paste(df$sex,df$Time,df$Age)
	Phidata=Phidata[Phidata$key %in% df$key,]
	df=df[df$key %in% Phidata$key,]
	df=df[order(df$sex,df$Time,df$Age),]
	
	Phidata$time=factor(Phidata$Time+1987)
#	Phidata$time=factor(cut(Phidata$Time+1987,1986:2015,include.lowest=TRUE))
	Phidata$pup=ifelse(Phidata$Age==0,1,0)
	Phidata$yearling=ifelse(Phidata$Age==1,1,0)
	Phidata$twothree=ifelse(Phidata$Age%in%2:3, 1,0)
	Phidata$threeplus=ifelse(Phidata$Age>=3,1,0)
	Phidata$fourplus=ifelse(Phidata$Age>=4,1,0)
	
# indicator variables for <1993 (pre93) and >=1993 (p93plus); >=1989 post (excludes
# 1989 if non-pup).
	Phidata$p87plus=ifelse(Phidata$Time>=0&Phidata$time!="1999",1,0)
	Phidata$p88plus=ifelse(Phidata$Time>=1&Phidata$time!="1999",1,0)
	Phidata$p89plus=ifelse(Phidata$Time>=2&Phidata$time!="1999",1,0)
	Phidata$p90plus=ifelse(Phidata$Time>=3&Phidata$time!="1999",1,0)
	Phidata$p91plus=ifelse(Phidata$Time>=4&Phidata$time!="1999",1,0)
	Phidata$p94plus=ifelse(Phidata$Time>=7,1,0)
	Phidata$timebin=cut(Phidata$Time,c(0,4,6:28),include.lowest=TRUE,labels=c("1988-1991","1992-1993",1994:2015))
    #Phidata$timebin=cut(Phidata$Time,c(0,4,6:28),include.lowest=TRUE,labels=c("1988-1991","1992-1993",1994:2015))
	Phidata$plus4time=cut(Phidata$Time,c(0,7:28),include.lowest=TRUE,labels=c("1988-1994",1995:2015))
	Phidata$AgeS=cut(Phidata$Age,c(0,1,2,4,7,10,13,17,30),right=FALSE)
	Phidata$time.intervals=1
	Phidata$time.intervals[Phidata$Age==0]=1/c(0.779,0.773,0.689,0.661,0.694,0.629,0.697,0.742,0.753,0.758,0.766,0.684,0.757,0.753,0.735,0.756,0.758,0.718,0.703,0.772,0.778,0.775,0.683,0.764,0.773,0.762,.767,.767)
	key.estimates=Phidata$key

    if(vcv)
	   pp=predict_real(model,df,"S",data=Phidata,beta=beta,vcv=vcv)
   else
	   pp=predict_real(model,df,"S",data=Phidata,beta=beta,vcv=FALSE,se=FALSE)
	return(list(pp=pp,key.estimates=key.estimates))
}

make_pup_counts=function(nreps)
{
	data(PupCounts)
	save_PupCounts=PupCounts
	incyears=c(1975:1977,1981:1984,1990:1997)
	snicount=PupCounts[PupCounts$Year%in%incyears,c("SNI_PupCnt")]
	nonsnicount=rowSums(PupCounts[PupCounts$Year%in%incyears,c("SMI_PupCnt","SCI_PupCnt","SBI_PupCnt")])
	mod=lm(log(snicount)~log(nonsnicount))
	sni_pred=predict(mod,newdata=data.frame(nonsnicount=rowSums(PupCounts[PupCounts$Year%in%1985:1989,c("SMI_PupCnt","SCI_PupCnt","SBI_PupCnt")])),se.fit=TRUE)
	
	PupCounts$SNI_PupCnt[PupCounts$Year%in%1985:1989]=sni_pred$fit
	PupCounts$US_PupCnt[PupCounts$Year%in%1985:1989]=rowSums(PupCounts[PupCounts$Year%in%1985:1989,c("SMI_PupCnt","SNI_PupCnt","SCI_PupCnt","SBI_PupCnt","Other_PupCnt")])
# now predict number with linear regression between 1978 to 1980
	incyears=c(1975:1977,1981:1982,1990:1997,1999:2001)
	mod=lm(log(US_PupCnt)~Year,data=PupCounts[PupCounts$Year%in%incyears,])
	total_pred=predict(mod,newdata=data.frame(Year=1978:1980),se.fit=TRUE)
	PupCounts$US_PupCnt[PupCounts$Year%in%1978:1980]=exp(total_pred$fit+total_pred$residual.scale^2/2)
	pup_counts=PupCounts$US_PupCnt
# finally compute vaues for 2009, 2010 and 2015 from ground counts at SMI and SNI(partial) 
	smi_counts=c(788,873,508)+c(12018,14258,15998)
	sni_counts=c(6328,4906,3064)-c(74,62,0)
# compute proportion in area 32-54 to expand to total from ground count
	data(SNIGroundProp)
	ground_prop=SNIGroundProp$Area32.54/SNIGroundProp$Total
	Year=I(SNIGroundProp$Year-1990)
	mod=glm(ground_prop~poly(Year,3),family=gaussian(link="log"))
	df= summary(mod)$df[2]
	snitotal_pred=predict(mod,newdata=data.frame(Year=c(2009,2010,2015)-1990),type="response",se.fit=TRUE)
# from SMI count and SNI estimated total compute total US population from SNISMI-Total model
	SMISNI_PupCnt=PupCounts$SMI_PupCnt+PupCounts$SNI_PupCnt
	mod_SNISMI=lm(US_PupCnt~SMISNI_PupCnt,data=PupCounts)
	fromSNISMI_pred=predict(mod_SNISMI,newdata=data.frame(SMISNI_PupCnt=sni_counts/snitotal_pred$fit+smi_counts),se.fit=TRUE)
	pup_counts[PupCounts$Year%in%c(2009,2010,2015)]=fromSNISMI_pred$fit

# for error in counts need std errors from sni_pred, total_pred, snitotal_pred and fromSNISMI_pred
# for prediction of 5 missing SNI totals
	slope_error=rt(nreps,df=sni_pred$df)
	errors_snipred=t(t(matrix(slope_error,nrow=nreps,ncol=5))*sni_pred$se.fit)+matrix(rnorm(5*nreps,0,sni_pred$residual.scale),ncol=5)
# for prediction of 1978-1980
	slope_error=rt(nreps,df=total_pred$df)
	errors_totalpred=t(t(matrix(slope_error,nrow=nreps,ncol=3))*total_pred$se.fit)+matrix(rnorm(3*nreps,0,total_pred$residual.scale),ncol=3)
# both for prediction of 2009,2010,2015 from ground counts
	slope_error=rt(nreps,df=df)
	errors_snitotalpred=t(t(matrix(slope_error,nrow=nreps,ncol=3))*snitotal_pred$se.fit)+matrix(rnorm(3*nreps,0,snitotal_pred$residual.scale),ncol=3)
	errors_fromSNISMI=matrix(0,nrow=nreps,ncol=3)
	for(i in 1:nreps)
	{
		xx=predict(mod_SNISMI,newdata=data.frame(SMISNI_PupCnt=sni_counts/(snitotal_pred$fit+errors_snitotalpred[i,])+smi_counts),se.fit=TRUE)
		slope_error=rt(1,df=fromSNISMI_pred$df)
		errors_fromSNISMI[i,]=slope_error*xx$se.fit+rnorm(3,0,xx$residual.scale)
	}
	Nlist=vector("list",length=nreps)
	for(i in 1:nreps)
	{
		PC=save_PupCounts
		PC$SNI_PupCnt[PC$Year%in%1985:1989]=exp(sni_pred$fit+errors_snipred[i,])
		PC$US_PupCnt[PC$Year%in%1985:1989]=rowSums(PC[PC$Year%in%1985:1989,c("SMI_PupCnt","SNI_PupCnt","SCI_PupCnt","SBI_PupCnt","Other_PupCnt")])
		PC$US_PupCnt[PC$Year%in%1978:1980]=exp(total_pred$fit+errors_totalpred[i,]+total_pred$residual.scale^2/2)
		PC$US_PupCnt[PC$Year%in%c(2009,2010,2015)]=fromSNISMI_pred$fit+errors_fromSNISMI[i,]
		Nlist[[i]]=PC$US_PupCnt
	}
	return(list(pup_counts=pup_counts,Nlist=Nlist))
}

population_reconstruction_error=function(nreps,pup_counts,model,zc.ddl,begin=1975,r=0.041,reml.formula)
{
	    model=RMark:::load.model(model)
		which.beta=grep("S:",rownames(model$results$beta))
		Phi_coef=model$results$beta$estimate[which.beta]
		Phi_vcv= model$results$beta.vcv[which.beta,which.beta]
	    Phi_rep = rmvnorm(nreps, Phi_coef, Phi_vcv)
		Phi_vcv=make_predictions(model,zc.ddl,vcv=TRUE)$pp$vcv
		mnpl_t=NULL
		for (i in 1:nreps)
		{
			pred.list=make_predictions(model,zc.ddl,beta=Phi_rep[i,],vcv=FALSE)
			pred.list$pp=list(real=pred.list$pp,vcv=pred.list$vcv)
			pred.list$pp$vcv=Phi_vcv
			stats=population_reconstruction(pred.list$pp,reml.formula,pred.list$key.estimates,pup_counts[[i]],r=r)$stats
			mnpl_t=rbind(mnpl_t,stats)
		}
		rownames(mnpl_t)=NULL
		return(mnpl_t)
}	
