#' Population projection
#' 
#' Projects a population with a 25 year age structure from a series of 
#' annual pup counts and a sex- and age-specific survival rate matrix.
#' Animals can reach age 24 and then they disappear (S from age 24 to 25 is assumed to be 0).
#' 
#' @param pup_counts vector of annual pup counts
#' @param Sf a matrix of age by year survival rates for females
#' @param Sm a matrix of age by year survival rates for males
#' @param r assumed rate of increase at beginning of time series to create initial stable age distribution
#' @param begin year to initiate population projection
#' @param plot if TRUE, produces plots; otherwie, only returns stats
#' @param stochastic if TRUE, includes demographic stochasticity
#' @return  list with population array of numbers at each sex, age and time and vector of statistics for MNPL etc
#' @author Jeff Laake
#' @export 
project_pop=function(pup_counts,Sf, Sm, r=0.041,begin=1975,plot=TRUE,span=1/3,stochastic=TRUE,initial=age_dist(Sf,Sm,r,pup_counts[1]))
{
	# compute initial stable age distribution and population size from pup count in first year
	#xx=age_dist(Sf,Sm,r,pup_counts[1])
	nyears=length(pup_counts)
	# initialize population matrix with 2 sexes, 25 age classes and nyears
	zcpop=array(0,dim=c(2,25,nyears))
	# initialize first year with state age distribution
	zcpop[1,1:25,1]=initial$N_f[1:25]
	zcpop[2,1:25,1]=initial$N_m[1:25]
	# put in pup counts for each year assuming a 50:50 sex ratio
	zcpop[1,1,]=pup_counts/2
	zcpop[2,1,]=pup_counts/2
	# loop over years projecting each pup cohort over time with a binomial distribution for the 
	# number of survivors from size of previous years cohort and the survival rate for each age-time
	for(year in 1:(nyears-1))
	{
		for(y in year:(nyears-1))
		{
			if(any(is.na(Sf[,year])))browser()
			if(any(is.na(Sm[,year])))browser()
			if(any(is.na(zcpop[1,1:24,y])) | any(zcpop[1,1:24,y]<0) )browser()
			if(any(is.na(zcpop[2,1:24,y]))| any(zcpop[2,1:24,y]<0))browser()	
			if(!stochastic)
			{
				zcpop[1,2:25,y+1]=zcpop[1,1:24,y]*Sf[,year]
				zcpop[2,2:25,y+1]=zcpop[2,1:24,y]*Sm[,year]
			}else
			{
				zcpop[1,2:25,y+1]=rbinom(24,floor(zcpop[1,1:24,y]),Sf[,year])
				zcpop[2,2:25,y+1]=rbinom(24,floor(zcpop[2,1:24,y]),Sm[,year])
			}
		}
	}
	# return the population array of numbers
	dimnames(zcpop)=list(sex=c("Female","Male"),age=0:24,years=begin:(begin+nyears-1))
	stats=plot_pop(zcpop,begin=begin,plot=plot,span=span)
	return(list(zcpop=zcpop,stats=stats))
}

