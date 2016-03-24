#' Population projection plots
#' 
#' Computes popn statistics and optionally plots smoothed net change, pup counts, total population size, non-pup population size over time 
#'  
#' @param zcpop  population array of numbers at each sex, age and time
#' @param begin year to initiate population projection
#' @param plot if TRUE, produces plots; otherwie, only returns stats
#' @return vector containing year of maximum smoothed net increase, smoothed population size at MNPL and current smoothed N, all for non-pups.
#' @author Jeff Laake
#' @import mgcv
#' @export 

plot_pop=function(zcpop,begin=1975,plot=TRUE,span=1/3)
{
	year=begin+1:(dim(zcpop)[3]-1)
	pop=apply(zcpop,3,sum)
	popr=(pop[3:dim(zcpop)[3]]-pop[1:(dim(zcpop)[3]-2)])
	tpop=colSums(zcpop[1,,])+colSums(zcpop[2,,])
	sm=gam(popr~s(year[-length(year)]))
	pred.sm=predict(sm,se=TRUE)
	smnp=gam(pop~s(c(begin,year)))
	pred.pop=predict(smnp)
	smtp=gam(tpop~s(c(begin,year)))
	pred.tpop=predict(smtp,se.fit=TRUE)
	if(plot)
	{
		par(mfrow=c(3,2),cex.lab=1.25,cex.axis=1.25,bg="transparent")
		years=year[-length(year)]
		plot(years,popr,ylab="Net change",xlab="Year")
		lines(years,pred.sm$fit)
		lines(years,pred.sm$fit+2*pred.sm$se.fit,lty=2)
		lines(years,pred.sm$fit-2*pred.sm$se.fit,lty=2)
		plot(c(begin,year),tpop,ylab="Total population size",xlab="Year")
		lines(c(begin,year),pred.tpop$fit)
		lines(c(begin,year),pred.tpop$fit+2*pred.tpop$se.fit,lty=2)
		lines(c(begin,year),pred.tpop$fit-2*pred.tpop$se.fit,lty=2)
		plot(year-1,diff(pred.tpop$fit),ylab="Smoothed net change",xlab="Year",type="b")
#		plot(smnp$x,pop,ylab="Non-pup population size",xlab="Year")
#		lines(smnp$x,smnp$y)
		plot(c(1975,year),zcpop[1,1,]+zcpop[2,1,],ylab="Pup count",xlab="Year")
		plot(c(1975,year),(zcpop[1,1,]+zcpop[2,1,])/colSums(zcpop[1,5:25,]),xlab="Year",ylab="Pups/4+ females")
		abline(h=1)
		plot(c(1975,year),tpop/(zcpop[1,1,]+zcpop[2,1,]),xlab="Year",ylab="Implied correction factor")
	}
	maxN= max(apply(zcpop,3,sum))
	names(maxN)=names(which.max(apply(zcpop,3,sum)))
	return(c(mnpl_year=year[which.max(diff(pred.tpop$fit))]-1,mnpl_N=sum(zcpop[,,which.max(diff(pred.tpop$fit))]),current_N=sum(zcpop[,,length(smnp$y)]),maxN=maxN,maxsmN=max(pred.tpop$fit)))
#	return(c(mnpl_year=year[which.max(sm$y)]-1,mnpl_N=sum(zcpop[,,which.max(sm$y)]),current_N=sum(zcpop[,,length(smnp$y)]),maxN=maxN))
}



