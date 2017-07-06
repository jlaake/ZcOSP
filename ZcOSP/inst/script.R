library(ZcOSP)
library(RMark)
library(splines)
library(lme4)
library(mgcv)
library(optimx)
library(mvtnorm)
##############################################################################################################
# Define Functions
##############################################################################################################

get_Scoef=function(nreps,model)
{
	which.beta=grep("S:",rownames(model$results$beta))
	# remove 2015 which are not estimable &remove 1996 2+ male and then replace with 1995 value
	rep_index=c(which(rownames(model$results$beta)[which.beta]=="S:timebin1996:p89plus:male:twoplus"),grep("2015",rownames(model$results$beta)[which.beta]))

	Phi_coef=model$results$beta$estimate[which.beta][-rep_index]
	Phi_vcv= model$results$beta.vcv[which.beta,which.beta][-rep_index,-rep_index]
	Phi_rep=matrix(NA,nrow=nreps,ncol=length(which.beta))
	Phi_rep[,which.beta[-rep_index]] = rmvnorm(nreps, Phi_coef, Phi_vcv)
	Phi_rep[,rep_index]=Phi_rep[,rep_index-1]
	return(Phi_rep)
}

get_estimates_pop=function(z,zcpop,fixed,init_model,pup_counts,stochastic=FALSE,which=NULL,reps=10)
{
	K<-250000  
	z$N=NULL
	Ndf=data.frame(time=factor(1987:2014),N=apply(zcpop,3,sum)[-(1:12)])
	z$seq=1:nrow(z)
	z=merge(z,Ndf,by="time",all.x=TRUE)
	z=z[order(z$seq),]
	z$seq=NULL
	z$NK=z$N/K
	models=fitmixed(fixed.f=fixed,random.f=list(init_model$best.r),data=z[z$time!=2014,])
	tt=1:7
	r=coef(lm(log(pup_counts[1:7])~tt))[2]
	
	# loop over model that may have N/K in it
	for(i in 1:reps)
	{
		tt=1:7
		r=coef(lm(log(pup_counts[1:7])~tt))[2]
		z$N=NULL
		Ndf=data.frame(time=factor(1987:2014),N=apply(zcpop,3,sum)[-(1:12)])
		z$seq=1:nrow(z)
		z=merge(z,Ndf,by="time",all.x=TRUE)
		z=z[order(z$seq),]
		z$seq=NULL
		z$NK=z$N/K
		# fit selected lme model	
		mod=lmer(paste(models$best.f,"+(", models$best.r,")",sep=""), data = z[z$time!=2014,], 
				REML=TRUE)
		if(stochastic)
		{
			if(i==1)
			{
				sd=sqrt(diag(VarCorr(popn$mod)$time))
				puperr=rnorm(39,0,sd[1])
				yrerr=rnorm(39,0,sd[2])
				twopluserr=rnorm(39,0,sd[3])
				errdf=data.frame(time=as.character(1975:2013),Age=rep(0,39),err=puperr)
				errdf=rbind(errdf,data.frame(time=as.character(1975:2013),Age=rep(1,39),err=yrerr))
				for(i in 2:24)
					errdf=rbind(errdf,data.frame(time=as.character(1975:2013),Age=rep(i,39),err=twopluserr))
				errdf$key=paste(as.character(errdf$time),errdf$Age,sep="")
			}
		} else
			errdf=NULL
		
		####################################################################################################	
#	 get predictions for all survivals - use level 0 (fixed effect estimates) for missing survivals
		####################################################################################################
		Phidata=SavePhiData[SavePhiData$Time!=25,]
		Phidata$Time=Phidata$Time+2
		Phidata=Phidata[order(Phidata$sex,Phidata$Time,Phidata$Age),]
		Phidata$time=factor(Phidata$Time+1987)
		Phidata$pup=ifelse(Phidata$Age==0,1,0)
		Phidata$yearling=ifelse(Phidata$Age==1,1,0)
		Phidata$twothree=ifelse(Phidata$Age%in%2:3, 1,0)
		Phidata$threeplus=ifelse(Phidata$Age>=3,1,0)
		Phidata$fourplus=ifelse(Phidata$Age>=4,1,0)
		Phidata$timebin=cut(Phidata$Time,c(0,4,6:28),include.lowest=TRUE,labels=c("1988-1991","1992-1993",1994:2015))
		Phidata$plus4time=cut(Phidata$Time,c(0,7:28),include.lowest=TRUE,labels=c("1988-1994",1995:2015))
		Phidata$AgeS=cut(Phidata$Age,c(0,4,7,10,13,17,30),right=FALSE)
		Phidata$age=factor(Phidata$Age)
		Phidata$key=paste(Phidata$sex,Phidata$Time,Phidata$Age)
		Phidata$py=Phidata$pup+Phidata$yearling
		Ndf=data.frame(time=factor(1975:2014),N=apply(zcpop,3,sum))
		Phidata$seq=1:nrow(Phidata)
		Phidata=merge(Phidata,Ndf,by="time",all.x=TRUE)
		Phidata=Phidata[order(Phidata$seq),]
		Phidata$key=paste(as.character(Phidata$time),Phidata$Age,sep="")
		if(!is.null(errdf))
		{
			Phidata=merge(Phidata,errdf[,c("key","err")],by="key",all.x=TRUE)
			Phidata=Phidata[order(Phidata$seq),]
		} else
			Phidata$err=0
		Phidata$seq=NULL
		Phidata$NK=Phidata$N/K
		level1=plogis(predict(mod))
		level0=plogis(predict(mod,newdata=Phidata,re.form=NA)+Phidata$err)
		Phidata$key=paste(Phidata$sex,Phidata$Time,Phidata$Age)
		which.reals=Phidata$key%in%key.estimates
		Phidata$S=NA
		Phidata$S[which.reals]=level1
		Phidata$S[!which.reals]=level0[!which.reals]
		
		######################################################################
#   Reconstruct population
		######################################################################		
		Sf=matrix(Phidata$S[Phidata$sex=="F"],ncol=39,nrow=25)
		Sm=matrix(Phidata$S[Phidata$sex=="M"],ncol=39,nrow=25)
		Sf[1,]=Sf[1,]^(11/12)
		Sm[1,]=Sm[1,]^(11/12)
		Sm=Sm[-25,]
		Sf=Sf[-25,]
		initial=age_dist(Sf,Sm,pups=pup_counts[1],r=r)	
		birth_param=get_estimates_br(zcpop,DDT=get_DDT(),which=which,predict=TRUE)
#		if(debug)cat("\n",birth_param$br_coef)
#		f4plus=apply(zcpop[1,5:24,],2,sum)
#		pup_counts[4:6]=birth_param$br[4:6]*f4plus[4:6]
		res=project_pop(pup_counts[-41],Sf,Sm,stochastic=FALSE,initial=initial,r=r)
		zcpop=res$zcpop
		K =res$stats["maxsmN"]
		if(debug)cat("likelihood = ",logLik(mod)," K = ", K," r =",r)
	}
	if(stochastic)
		res=project_pop(pup_counts[-41],Sf,Sm,stochastic=TRUE,initial=initial,r=r)
	return(list(res=res,mod=mod,birth_vc=birth_param$vc))
}

get_DDT=function(stochastic=FALSE)
{
	Time=c(1970,1972, 1991, 2002)
	DDT=c((981*1303/(1303+6061)+121*6061/(1303+6061)),(752*1769/(1769+6133)+129*6133/(1769+6133)), (40.4*552/(552+15696)+9.46*(15696/(552+15696))),10.8)
	mod=lm(log(DDT)~I(Time-1971))
	predlist=predict(mod,newdata=data.frame(Time=1975:2014),se.fit=TRUE)
	if(!stochastic)
		DDT = exp(predlist$fit+summary(mod)$sigma^2/2)
	else
	{
		coef=as.vector(rmvt(1,vcov(mod),df=summary(mod)$df[1]))+coef(mod)
		DDT = exp(model.matrix(~Time,data.frame(Time=1975:2014-1971))%*%coef+summary(mod)$sigma^2/2)
	}   
	return(DDT)
}		

get_estimates_br=function(zcpop,DDT,stochastic=FALSE,which=NULL,predict=FALSE,model.average=FALSE)
{	
	# compute predicted DDT values over time from values measured for females in 1971,1991 and 2002
	
	data(Phidata)
	SST=Phidata$JulytoJuneSSTAnomalies[Phidata$sex=="F"&Phidata$Age==0]
	natality=data.frame(time=1975:2014,Pups=floor(apply(zcpop[,1,],2,sum)),Nf=floor(apply(zcpop[1,5:25,],2,sum)),NK=floor(apply(zcpop,3,sum))/K,SST=c(-.8,SST[-40]),DDT=DDT)
#	natality$Time=natality$time-1974
#   natality$seq=1:nrow(natality)
#	fishx=fish
#	fishx$Time=fishx$Time+1
#	fishx$time=factor(fishx$Time)
#	natality=merge(natality,fishx,by="time",all.x=TRUE)
#	natality=natality[order(natality$seq),]
#	natality$seq=NULL
#	rm(fishx)
	
	#SST
	par=c(8,1,0)
	nat_formula=~SST
	
	mod_SST=suppressWarnings(optimx(par=par,br_fit,method=c("nlminb"),nat_formula=nat_formula,natality=natality))
	mod_SST$AICc=AICvalue
	mod_SST$npar=npar
	if(debug&(is.null(which)||which==1))cat("birth neglnl=",br_fit(c(mod_SST$p1,mod_SST$p2,mod_SST$p3),nat_formula,natality))
	#SST +DDT
	par=c(mod_SST$p1,mod_SST$p2,mod_SST$p3,0 )
	nat_formula=~SST+DDT
	mod_SSTDDT=suppressWarnings(optimx(par=par,br_fit,method=c("nlminb"),nat_formula=nat_formula,natality=natality))
	mod_SSTDDT$AICc=AICvalue
	mod_SSTDDT$npar=npar
	if(debug&(is.null(which)||which==2))cat("birth neglnl=",br_fit(c(mod_SSTDDT$p1,mod_SSTDDT$p2,mod_SSTDDT$p3,mod_SSTDDT$p4),nat_formula,natality))
	
	#SST + NK
	par=c(mod_SST$p1,mod_SST$p2,mod_SST$p3,0 )
	nat_formula=~SST+NK
	mod_SSTNK=suppressWarnings(optimx(par,br_fit,method=c("nlminb"),nat_formula=nat_formula,natality=natality))
	mod_SSTNK$AICc=AICvalue
	mod_SSTNK$npar=npar
	if(debug&(is.null(which)||which==3))cat("birth neglnl=",br_fit(c(mod_SSTNK$p1,mod_SSTNK$p2,mod_SSTNK$p3,mod_SSTNK$p4),nat_formula,natality))
	
	#  SST +DDT +NK
	par=c(mod_SSTDDT$p1,mod_SSTDDT$p2,mod_SSTDDT$p3,mod_SSTDDT$p4,0 )
	nat_formula=~SST+DDT+NK
	mod_SSTDDTNK=suppressWarnings(optimx(par,br_fit,method=c("nlminb"),nat_formula=nat_formula,natality=natality))
	mod_SSTDDTNK$AICc=AICvalue
	mod_SSTDDTNK$npar=npar
	if(debug&(is.null(which)||which==4))cat("birth neglnl=",br_fit(c(mod_SSTDDTNK$p1,mod_SSTDDTNK$p2,mod_SSTDDTNK$p3,mod_SSTDDTNK$p4,mod_SSTDDTNK$p5),nat_formula,natality))
	AICcvalues=c(mod_SST$AICc,mod_SSTDDT$AICc,mod_SSTNK$AICc,mod_SSTDDTNK$AICc)
	whichbest=which.min(AICcvalues)
	wts=AICcvalues-min(AICcvalues)
	wts=exp(-.5*wts)/sum(exp(-.5*wts))
	if(!is.null(which))whichbest=which
	if(whichbest==1)
	{
		br_coef=c(mod_SST$p1,mod_SST$p2,mod_SST$p3,0,0)
		vc=solve(attr(mod_SST,"details")[,"nhatend"][[1]])
	}
	if(whichbest==2)
	{
		br_coef=c(mod_SSTDDT$p1,mod_SSTDDT$p2,mod_SSTDDT$p3,mod_SSTDDT$p4,0)
		vc=solve(attr(mod_SSTDDT,"details")[,"nhatend"][[1]])
	}
	if(whichbest==3)
	{
		br_coef=c(mod_SSTNK$p1,mod_SSTNK$p2,mod_SSTNK$p3,0,mod_SSTNK$p4)
		vc=solve(attr(mod_SSTNK,"details")[,"nhatend"][[1]])
	}
	if(whichbest==4)
	{
		br_coef=c(mod_SSTDDTNK$p1,mod_SSTDDTNK$p2,mod_SSTDDTNK$p3,mod_SSTDDTNK$p4,mod_SSTDDTNK$p5)
		vc=solve(attr(mod_SSTDDTNK,"details")[,"nhatend"][[1]])
	}
	mavg=c(mod_SST$p1,mod_SST$p2,mod_SST$p3,0,0)*wts[1]+c(mod_SSTDDT$p1,mod_SSTDDT$p2,mod_SSTDDT$p3,mod_SSTDDT$p4,0)*wts[2]+
		c(mod_SSTNK$p1,mod_SSTNK$p2,mod_SSTNK$p3,0,mod_SSTNK$p4)*wts[3]+c(mod_SSTDDTNK$p1,mod_SSTDDTNK$p2,mod_SSTDDTNK$p3,mod_SSTDDTNK$p4,mod_SSTDDTNK$p5)*wts[4]
	if(!predict)
		if(model.average)
		    return(br_coef=mavg)
	    else
			return(br_coef=br_coef)
	else
	{
	    if(model.average) br_coef=mavg
		  gamma=br_coef[1]
	    beta_nat=br_coef[2:length(br_coef)]
	    dm_nat=model.matrix(nat_formula,natality)
	    br=plogis(dm_nat%*%beta_nat)
		return(list(br_coef=br_coef,br=br,vc=vc))
	}
}


br_fit=function(par,nat_formula,natality)
{
	gamma=par[1]
	beta_nat=par[2:length(par)]
	dm_nat=model.matrix(nat_formula,natality)
	br<-plogis(dm_nat%*%beta_nat)
	sigmasq=natality$Nf*br*(1-br)*(exp(gamma)+1)
	neglnl=sum((natality$Pups-natality$Nf*br)^2/(2*sigmasq) +log(sqrt(sigmasq)))
	npar<<-length(par)
	AICvalue<<-2*as.numeric(neglnl)+2*npar + 2*npar*(npar+1)/(nrow(dm_nat)-npar+1)
	return(as.numeric(neglnl))
}


seRatio=function(y,x)
{
	R=sum(y)/sum(x)
	return(sqrt((var(y)+R^2*var(x)-2*R*cov(x,y))/length(x)/mean(x)^2))
}

get_removals=function(bycatch,stochastic=FALSE)
{
	bycatch_rate_prereg=sum(bycatch$ObservedBycatch[bycatch$Year<=1994])/sum(bycatch$ObservedEffort[bycatch$Year<=1994])
	bycatch_rate_postreg=sum(bycatch$ObservedBycatch[bycatch$Year>1994])/sum(bycatch$ObservedEffort[bycatch$Year>1994])
	effort_mod=with(bycatch,gam(log(FishingEffort)~s(Year)))
	Effort=c(exp(predict(effort_mod,newdata=data.frame(Year=1975:1980))),bycatch$FishingEffort[-c(33:34)],exp(predict(effort_mod,newdata=data.frame(Year=2013:2014))))
	removals=Effort*bycatch_rate_prereg
	removals[21:40]=Effort[21:40]*bycatch_rate_postreg
	if(stochastic)
	{
		n_prereg=length(bycatch$ObservedBycatch[bycatch$ObservedEffort>0&bycatch$Year<=1994])
		se_preReg=seRatio(bycatch$ObservedBycatch[bycatch$ObservedEffort>0&bycatch$Year<=1994],bycatch$ObservedEffort[bycatch$ObservedEffort>0&bycatch$Year<=1994])
		n_postreg=length(bycatch$ObservedBycatch[bycatch$ObservedEffort>0&bycatch$Year>1994])
		se_postReg=seRatio(bycatch$ObservedBycatch[bycatch$ObservedEffort>0&bycatch$Year>1994],bycatch$ObservedEffort[bycatch$ObservedEffort>0&bycatch$Year>1994])
		removals=Effort*(bycatch_rate_prereg+rt(length(Effort),n_prereg-1)*se_preReg)
		removals[21:40]=Effort[21:40]*(bycatch_rate_postreg+rt(length(removals[21:40]),n_postreg-1)*se_postReg)
		removals[removals<0]=0
	}
	return(removals)
}

get_logistic_fit=function(zcpop,start=NULL,removals=NULL)
{
	data(Phidata)
	SST=Phidata$JulytoJuneSSTAnomalies[Phidata$sex=="F"&Phidata$Age==0]
	N=apply(zcpop,3,sum)
	if(is.null(removals))
		removals=rep(0,40)
	panel=list(y.obs=N,time.index=1:40,removals=removals)
	if(is.null(start))
	{
		panel$z=2
		panel$Rm=0.05
		panel$n0=90000
		panel$K=267000
		panel$slope=0
	}
	else
		panel=c(panel,start)
	panel$SST=SST
	panel$proportional.errors=TRUE
	panel$lower=list(z=1,Rm=-.5,n0=20000,K=200000,slope=-5)
	panel$upper=list(z=20,Rm=.2,n0=120000,K=600000,slope=5)
	Nmod=NLS.PT(panel,upper=panel$upper)
	logLik(Nmod)
	par_N=coef(Nmod)
	pred=project.PT(par_N[1],par_N[2],par_N[3],par_N[4],par_N[5],1:40,SST=panel$SST,removals=removals)
	plot(1975:2014,pred,type="l",ylim=c(min(c(pred,panel$y.obs)),max(c(pred,panel$y.obs))),ylab="Total population size",xlab="Year")
	points(1975:2014,panel$y.obs)
	return(par_N)
}


plot_logistic_fit=function(zcpop,par_N,boot,removals=NULL)
{
	data(Phidata)
	SST=Phidata$JulytoJuneSSTAnomalies[Phidata$sex=="F"&Phidata$Age==0]
	N=apply(zcpop,3,sum)
	if(is.null(removals))
		removals=rep(0,40)
	panel=list(y.obs=N,time.index=1:40,removals=removals)
	panel$SST=SST
	pred=project.PT(par_N[1],par_N[2],par_N[3],par_N[4],par_N[5],1:40,SST=panel$SST,removals=removals)
	plot(1975:2014,pred,type="l",ylab="Total population size",xlab="Year",xaxp=c(1975,2014,39),ylim=c(60000,320000))
	points(1975:2014,panel$y.obs)
	bounds=t(sapply(logistic.boot, function(x) project.PT(x[1],x[2],x[3],x[4],x[5],1:40,SST=panel$SST,removals=removals)))
	lowerbounds=apply(bounds,2, function(x) sort(x)[.025*nreps])
	upperbounds=apply(bounds,2, function(x) sort(x)[.975*nreps])
	lines(1975:2014,lowerbounds,lty=2)
	lines(1975:2014,upperbounds,lty=2)
	invisible()
}

plot_error_bars=function(popn.boot)
{
	limits=apply(do.call(rbind,lapply(popn.boot,function(x)apply(x$res$zcpop,3,sum))),2,function(x)c(sort(x)[.025*nreps],sort(x)[.975*nreps]))
	for(i in 1:ncol(limits))
		lines(c(1974+i,1974+i),limits[,i])
}

fitmixed=function(fixed.f,random.f,data,save.model=TRUE,...)
{
	if(class(fixed.f)=="formula") fixed.f=list(fixed.f)
	if(length(fixed.f)>1)
	{
		if(length(random.f)>1) stop("cannot specify multiple sets of fixed and random formula")
		REML=FALSE
	} else
	{
		if(length(random.f)>1)
			REML=TRUE
	}
	results=vector("list",length(fixed.f)*length(random.f))
	model.table=data.frame(fixed=rep(NA,length(fixed.f)*length(random.f)),random=rep(NA,length(fixed.f)*length(random.f)))
	i=0
	for(f in fixed.f)
		for(r in random.f)
		{
			i=i+1
			results[[i]]=lmer(formula=as.formula(paste(f,"+(",r,")",sep="")),data=data,REML=REML,...)
			if(!save.model)
			{
				results[[i]]$data=NULL
				results[[i]]$fitted=NULL
				results[[i]]$residuals=NULL
			}
			cf=as.character(f)
			cf[1]=cf[2]
			cf[2]="~"
			model.table[i,]=cbind(fixed=f,random=r)
		}
	model.table$AIC=sapply(results,AIC)
	results$model.table=model.table[order(model.table$AIC),]
	results$best=as.numeric(row.names(results$model.table)[1]) 
	if(length(fixed.f)==1)
		results$best.f=fixed.f[[1]]
	else
		results$best.f=fixed.f[[results$best]]
	if(length(random.f)==1)
		results$best.r=random.f[[1]]
	else
		results$best.r=random.f[[results$best]]
	return(results)
}


project.PT <-function (z,Rm,n0,K,slope,times,SST,removals=NULL)
{
	n=vector(mode="numeric",length=max(times))
	n[1]  = n0
	for(i in 2:max(times))
		if(all(removals==0))
			n[i] = n[i-1]*(1+(Rm+slope*SST[i-1])*(1-(n[i-1]/(K))^z))
		else  
			n[i] = n[i-1]*(1+(Rm+slope*SST[i-1])*(1-(n[i-1]/(K))^z))-removals[times[i-1]]
	n
}
nls.PT.mod <-function (z,Rm,n0,K,slope,y.obs,time.index,SST,removals=NULL)
{
	if(debug) cat("z= ",z," Rm= ",Rm," n0= ",n0," K= ",K,"slope=",slope,"\n")
	
	if(!is.null(removals))
		n=project.PT (z,Rm,n0,K,slope,time.index,SST,removals[!is.na(removals)])
	else
		n=project.PT (z,Rm,n0,K,slope,time.index,SST)
	resid=n[time.index]-y.obs
	resid[y.obs==0]=0
	resid
}
nls.PT.mod.pe <-function (z,Rm,n0,K,slope,y.obs,time.index,SST,removals)
{
	if(debug) cat("z= ",z," Rm= ",Rm," n0= ",n0," K= ",K,"slope=",slope,"\n")
	if(!is.null(removals))
		n=project.PT (z,Rm,n0,K,slope,time.index,SST,removals[!is.na(removals)])
	else
		n=project.PT (z,Rm,n0,K,slope,time.index,SST)
	resid=(n[time.index]-y.obs)/n[time.index]
	resid[y.obs==0]=0
	resid
}

NLS.PT<-function (panel,MaxIter=5000,upper=NULL,...)
{
	y.obs=panel$y.obs
	time.index=panel$time.index
	removals=panel$removals
	SST=panel$SST
	dfm = data.frame(y.obs=y.obs,time.index=time.index,removals=removals,SST=SST)
	dfm$y.obs[is.na(dfm$y.obs)]=0
	z = panel$z ; Rm = panel$Rm
	n0 = panel$n0 ; K = panel$K
	slope=panel$slope
	if(is.null(upper))
		if(panel$proportional.errors)
			nls(~nls.PT.mod.pe(z,Rm,n0,K,slope,y.obs,time.index,SST,removals), data=dfm, start=list(z=z,Rm=Rm,n0=n0,K=K,slope=slope),
					nls.control(maxit=MaxIter,...))
		else
			nls(~nls.PT.mod(z,Rm,n0,K,slope,y.obs,time.index,SST,removals), data=dfm, start=list(z=z,Rm=Rm,n0=n0,K=K,slope=slope),
					nls.control(maxit=MaxIter,...))
	else
	if(panel$proportional.errors)
		nls(~nls.PT.mod.pe(z,Rm,n0,K,slope,y.obs,time.index,SST,removals), data=dfm,start=list(z=z,Rm=Rm,n0=n0,K=K,slope=slope),algorithm="port",
				upper=upper,lower=panel$lower,  nls.control(maxit=MaxIter,...))
	else
		nls(~nls.PT.mod(z,Rm,n0,K,slope,y.obs,time.index,SST,removals), data=dfm,start=list(z=z,Rm=Rm,n0=n0,K=K0,slope=slope),algorithm="port",
				upper=upper,lower=panel$lower,  nls.control(maxit=MaxIter,...))
}



########################End of Functions #######################################################################

nreps=1000
debug=FALSE
#attach bycatch data
data(bycatch)
# get PupCounts 
pup_count_list=make_pup_counts(nreps)
pup_counts=pup_count_list$pup_counts
# This is nreps bootstraps of the imputed counts with the actual counts
Nlist=pup_count_list$Nlist
# Get estimated S values and all estimates
data(best) # loads into model
# modify estimates to handle boundary parameters
which.beta=grep("S:",rownames(model$results$beta))
save.which=which.beta
# remove 2015 which are not estimable and not used; replace with 2014
which.beta[which.beta%in%grep("2015",rownames(model$results$beta)[which.beta])]=grep("2015",rownames(model$results$beta)[which.beta])-1
# remove 1996 2+ male and then replace with 1995 value in Phi_rep
which.beta[which.beta==104]=103
model$results$beta[save.which,]=model$results$beta[which.beta,]
model$results$beta.vcv[save.which,save.which]=model$results$beta.vcv[which.beta,which.beta]
Phi_coef=model$results$beta$estimate[which.beta]
Phi_vcv= model$results$beta.vcv[which.beta,which.beta]

data(SST)
data(zcddl)
data(zcpop)
predlist=make_predictions(model,zc.ddl)
pp=predlist$pp
key.estimates=predlist$key.estimates
z=pp$real
z$py=z$pup+z$yearling
z$oneplus=ifelse(z$age==0,0,1)
z$link=RMark:::compute.links.from.reals(pp$real$real,model,vcv.real=pp$vcv.real)$estimates
z$N=NULL
Ndf=data.frame(time=factor(1987:2014),N=apply(zcpop,3,sum)[-(1:12)])
z$seq=1:nrow(z)
z=merge(z,Ndf,by="time",all.x=TRUE)
z=z[order(z$seq),]
z$seq=NULL
K=250000
z$NK=z$N/K
z=z[z$time!=2014,]
z=droplevels(z)
data(Phidata)
SavePhiData=Phidata
randomlist=list("-1+pup+yearling+twoplus|time","-1+pup+oneplus|time","1|time")
fixed=list(	
		"link~sex * bs(Age) + OcttoJuneSSTAnomalies:pup + yearling:JulytoJuneSSTAnomalies + twothree:JulytoJuneSSTAnomalies + fourplus:JulytoJuneSSTAnomalies + pup:NK +   yearling:NK + twothree:NK + fourplus:NK + pup:weight + yearling:weight",
		"link~sex * bs(Age) + OcttoJuneSSTAnomalies:pup + yearling:JulytoJuneSSTAnomalies + twothree:JulytoJuneSSTAnomalies+  fourplus:JulytoJuneSSTAnomalies + pup:NK + yearling:NK + twothree:NK + pup:weight + yearling:weight",
		"link~sex * bs(Age) + OcttoJuneSSTAnomalies:pup + yearling:JulytoJuneSSTAnomalies + twothree:JulytoJuneSSTAnomalies + pup:NK + yearling:NK + twothree:NK + pup:weight + yearling:weight",
		"link~sex * bs(Age) + OcttoJuneSSTAnomalies:pup + yearling:JulytoJuneSSTAnomalies + twothree:JulytoJuneSSTAnomalies + pup:NK + yearling:NK + pup:weight + yearling:weight",
		"link~sex * bs(Age) + OcttoJuneSSTAnomalies:pup + yearling:JulytoJuneSSTAnomalies +  pup:NK + yearling:NK + pup:weight + yearling:weight")
#fixed=list(	
#		"link~ JulytoJuneSSTAnomalies + NK","link~ JulytoJuneSSTAnomalies + NK")
#randomlist=list("1|time")

# load pup counts into population structure
zcpop[1,1,]=pup_counts[1:40]/2
zcpop[2,1,]=pup_counts[1:40]/2

init_model=fitmixed(fixed=fixed[1],random=randomlist,z)
popn=get_estimates_pop(z,zcpop,fixed,init_model,pup_counts,which=4)
birth=get_estimates_br(popn$res$zcpop,DDT=get_DDT(),model.average=TRUE)
birth_predictions=get_estimates_br(popn$res$zcpop,DDT=get_DDT(),model.average=TRUE,predict=T)
logistic=get_logistic_fit(popn$res$zcpop,removals=get_removals(bycatch))

popn.boot=vector("list",length=nreps)
Phi_rep=get_Scoef(nreps,model)
for(i in 1:nreps)
{
	cat("\n",i)
	predlist=make_predictions(model,zc.ddl,beta=Phi_rep[i,],vcv=FALSE)
	z=predlist$pp
	key.estimates=predlist$key.estimates
	z$py=z$pup+z$yearling
	z=droplevels(z)
	z$link=RMark:::compute.links.from.reals(z$real,model,vcv.real=pp$vcv.real)$estimates
	z$N=NULL
	Ndf=data.frame(time=factor(1987:2014),N=apply(zcpop,3,sum)[-(1:12)])
	z$seq=1:nrow(z)
	z=merge(z,Ndf,by="time",all.x=TRUE)
	z=z[order(z$seq),]
	z$seq=NULL
	z$NK=z$N/K
	popn.boot[[i]]=get_estimates_pop(z,zcpop,fixed,init_model,pup_counts=Nlist[[i]],stochastic=TRUE,reps=7)
}


br.boot=vector("list",length=nreps)
for(i in 1:nreps)
{
	cat("\n",i)
	br.boot[[i]]=get_estimates_br(popn.boot[[i]]$res$zcpop,DDT=get_DDT(stochastic=TRUE),stochastic=TRUE)
	cat("\n",br.boot[[i]])
}

br.boot.mavg=vector("list",length=nreps)
for(i in 1:nreps)
{
	cat("\n",i)
	br.boot.mavg[[i]]=get_estimates_br(popn.boot[[i]]$res$zcpop,DDT=get_DDT(stochastic=TRUE),stochastic=TRUE,model.average=TRUE)
	cat("\n",br.boot.mavg[[i]])
}

br.pred.mavg=vector("list",length=nreps)
for(i in 1:nreps)
{
  cat("\n",i)
  br.pred.mavg[[i]]=get_estimates_br(popn.boot[[i]]$res$zcpop,DDT=get_DDT(stochastic=TRUE),stochastic=TRUE,model.average=TRUE,predict=TRUE)$br
}

logistic.boot=vector("list",length=nreps)
for(i in 1:nreps)
{	   
	cat("\n",i)
	removals=get_removals(bycatch,stochastic=TRUE)
	logistic.boot[[i]]=try(get_logistic_fit(popn.boot[[i]]$res$zcpop,removals=removals))
	if(class(logistic.boot[[i]])=="try-error") 	    
		logistic.boot[[i]]=try(get_logistic_fit(popn.boot[[i]]$res$zcpop,start=list(z=3,Rm=0.09,K=300000,n0=60000,slope=0),removals=removals))
}

#plot of sex/age trends
win.metafile("plots.wmf",width=8,height=4.5)
par(mfrow=c(2,2),cex.lab=1.25,mar=c(2,4,2,2))
plot(1975:2014,apply(popn$res$zcpop[1,1:2,],2,sum),pch="F",type="b",ylim=c(0,60000),xlab="Year",ylab="Pups & Yearlings")
points(1975:2014,apply(popn$res$zcpop[2,1:2,],2,sum),pch="M",typ="b")
plot(1975:2014,apply(popn$res$zcpop[1,3:4,],2,sum),pch="F",type="b",ylim=c(0,60000),xlab="Year",ylab="Ages 2-3")
points(1975:2014,apply(popn$res$zcpop[2,3:4,],2,sum),pch="M",typ="b")
plot(1975:2014,apply(popn$res$zcpop[1,5:8,],2,sum),pch="F",type="b",ylim=c(0,60000),xlab="Year",ylab="Ages 4-7")
points(1975:2014,apply(popn$res$zcpop[2,5:8,],2,sum),pch="M",typ="b")
plot(1975:2014,apply(popn$res$zcpop[1,9:24,],2,sum),pch="F",type="b",ylim=c(0,60000),xlab="Year",ylab="Ages 8 and older")
points(1975:2014,apply(popn$res$zcpop[2,9:24,],2,sum),pch="M",typ="b")
dev.off()


# bootstrap values for MNPL
boots=do.call("rbind",logistic.boot)
MNPL_gl=apply(boots,1,function(x) x[4]*(x[1]+1)^-(1/x[1]))
MNPL_glcl=c(sort(MNPL_gl)[floor(nreps*.025)],sort(MNPL_gl)[floor(nreps*.975)])
N_MNPL=sapply(popn.boot,function(x) sum(x$res$zcpop[,,40]))
currentN_MNPL=sum(N_MNPL>MNPL_gl)/nreps
MNPL_gl=logistic[4]*(logistic[1]+1)^-(1/logistic[1])

