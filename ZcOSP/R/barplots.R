#' Population projection barplots
#' 
#' Plots population age structure over time 
#'  
#' @param zcpop  population array of numbers at each sex, age and time
#' @param begin year to initiate population projection
#' @param range a range of years to be plotted
#' @return list with cummulative female and male age structure by year
#' @author Jeff Laake
#' @export
barplot_pop=function(zcpop,begin=1975,range=NULL)
{
par(mfrow=c(2,1))
female_agestructure=apply(zcpop[1,,],2,function(x) c(x[1:2],sum(x[3:4]),sum(x[5:11]),sum(x[12:25])))
male_agestructure=apply(zcpop[2,,],2,function(x) c(x[1:2],sum(x[3:4]),sum(x[5:11]),sum(x[12:25])))
rownames(female_agestructure)=c("Pups","Yrling","Two-Three","Four-Ten","Eleven+")
rownames(male_agestructure)=c("Pups","Yrling","Two-Three","Four-Ten","Eleven+")
colnames(female_agestructure)=begin:(begin+dim(zcpop)[3])
colnames(male_agestructure)=begin:(begin+dim(zcpop)[3])
if(!is.null(range))
{
	female_agestructure=female_agestructure[,colnames(female_agestructure)%in%range,drop=FALSE]
	male_agestructure=male_agestructure[,colnames(male_agestructure)%in%range,drop=FALSE]
}
maxy=max(apply(female_agestructure,2,sum))
if(is.null(range))
{
	barplot(male_agestructure,legend=T,args.legend=list(x=13,y=1.1*maxy),col=c("red","orange","green","yellow","blue"),main="Male",ylim=c(0,maxy))
	barplot(female_agestructure,legend=T,args.legend=list(x=13,y=1.1*maxy),col=c("red","orange","green","yellow","blue"),main="Female",ylim=c(0,maxy))
}else
{
	barplot(male_agestructure,col=c("red","orange","green","yellow","blue"),main="Male",ylim=c(0,maxy))
	barplot(female_agestructure,col=c("red","orange","green","yellow","blue"),main="Female",ylim=c(0,maxy))
}
return(list(female=apply(t(t(female_agestructure)/colSums(female_agestructure)),2,cumsum),
male=apply(t(t(male_agestructure)/colSums(male_agestructure)),2,cumsum)))
}

