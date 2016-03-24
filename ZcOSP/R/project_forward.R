#' Forward Population projection
#' 
#' Projects an expected population with a 25 year age structure forward from the endpoint of
#' a population array with an assumed birth rate of reproductive females. It uses survival rate 
#' from final year in projecting forward.
#'  
#' @param population array of numbers at each sex, age and time
#' @param birth assumed rate of birth (pups/ repro female) 
#' @param nyears number of years forward to project
#' @return  population array of numbers at each sex, age and time extended foward
#' @author Jeff Laake
#' @export 
project_forward=function(zcpop,birth,nyears)
{
	newpop=array(0,dim=c(2,25,dim(zcpop)[3]+nyears))
	newpop[,,1:dim(zcpop)[3]]=zcpop
	for(i in 1:nyears)
	{
		newpop[1,2:25,dim(zcpop)[3]+i]=newpop[1,1:24,dim(zcpop)[3]+i-1]*Sf[,ncol(Sf)]
		newpop[2,2:25,dim(zcpop)[3]+i]=newpop[2,1:24,dim(zcpop)[3]+i-1]*Sm[,ncol(Sm)]
		repro_females=sum(newpop[1,5:25,dim(zcpop)[3]+i])
		pups=repro_females*birth
		newpop[,1,dim(zcpop)[3]+i]=pups/2
	}
	return(newpop)
}


