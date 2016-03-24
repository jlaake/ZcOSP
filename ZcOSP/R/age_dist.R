#' Stable age distribution
#' 
#' Computes sex-specific stable age distribution and abundance from pup count
#' from survival for females and males (Sf and Sm) and an assumed rate of growth.
#' The age proportions is the turned into abundance by computing N=pups/c(0) where c(0) 
#' is the proportion that are pups in the age distribution. Actually it is
#' pups/2/cf(0) + pups/2/cm(0) where cf(0) and cm(0) are the proportions that are pups for
#' females and males respectively. Then N is split into the sex and age categories based on
#' the stable age distribution.
#' 
#' @param Sf a matrix of age by year survival rates for females
#' @param Sm a matrix of age by year survival rates for males
#' @param r assumed rate of increase at beginning of time series to create initial stable age distribution
#' @param pups number of pups born in the initial year
#' @param proportions if TRUE, returns cx rather than N
#' @return  list with abundance by age for males (N_m), females (N_f) and total N.
#' @author Jeff Laake
#' @export 
#' 
age_dist=function(Sf,Sm,r=0.041,pups=NULL,proportions=FALSE)
{
	cxf=cumprod(c(1,Sf[,1],0))*exp(-r*(0:(length(Sf[,1])+1)))/sum(cumprod(c(1,Sf[,1],0))*exp(-r*(0:(length(Sf[,1])+1))))
	cxm=cumprod(c(1,Sm[,1],0))*exp(-r*(0:(length(Sm[,1])+1)))/sum(cumprod(c(1,Sm[,1],0))*exp(-r*(0:(length(Sm[,1])+1))))
	if(!proportions)
	{
		Nm=pups/2/cxm[1]
		Nf=pups/2/cxf[1]
	    N_m=Nm*cxm
	    N_f=Nf*cxf
	    return(list(N_m=N_m,N_f=N_f,N=Nm+Nf))
	}else
	{
		return(list(cxf=cxf,cxm=cxm))
	}
}

