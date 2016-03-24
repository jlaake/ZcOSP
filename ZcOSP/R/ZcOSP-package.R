#' Survival prediction data
#' 
#' @name Phidata
#' @docType data
#' @format A dataframe containing each age-sex-year combination and covariate data used to predict 
#' age-sex-time specific survival estimates from a mark-recapture model to be used in the 
#' ZcOSP model of the population dynamics.
#' @keywords datasets
#' @examples
#' #Create estimates from model for each age for 1975-2014
#' Phidata=expand.grid(Time=1975:2014,Age=0:24,sex=c("F","M"),weight=0)
#' # These are the observed average weights expressed as an anomaly from long term sex-specific mean
#' fwt=c(3.11615384615385,-1.68294117647059,0.0641176470588221,-2.10719298245614,-1.60045454545455,2.07952380952381,2.16555555555556,1.11939759036144,-4.49133333333333,-4.02772727272727,-0.510851063829787,1.20013698630137,-0.177272727272726,0.356315789473683,4.32555555555556,1.40959183673469,0.843333333333334,-2.72875,-0.422307692307694,-0.823169398907105,0.587791411042943,0.936134185303516,-3.07438040345821,-3.40334963325183,0.355761589403972,-0.556234567901235,-1.76404255319149,-0.753952095808383,0.165928753180662,3.58868852459016,3.50409836065574,0.449636363636362,0.435909090909089,-0.599743589743589,-1.22697986577181,-1.39880341880342,-3.7555230125523,-4.67135021097046,-2.68744680851064,-4.03443349753695)
#' mwt=c(3.85,-1.55387755102041,-0.166122448979593,-2.15069767441861,-1.64214285714286,2.58827586206896,1.30956521739131,1.35212121212121,-5.03164179104478,-3.31818181818182,1.59396825396825,1.37875,0.917777777777779,0.275294117647061,5.31777777777778,1.85496062992126,1.45750741839763,-2.65941860465116,-1.09893442622951,-0.64358208955224,1.13770114942529,1.19434782608696,-3.70278350515464,-4.52962457337884,0.295999999999999,-1.18459016393443,-2.55223300970874,-0.980555555555554,0.585454545454546,4.63950248756219,4.6175,1.15082251082251,1.04,-0.412173913043478,-2.96694444444444,-1.46421052631579,-3.64518518518518,-5.86983606557377,-2.84129032258064,-5.29898305084746)
#' Phidata$weight[Phidata$sex=="F"&Phidata$Age==0]=fwt
#' Phidata$weight[Phidata$sex=="M"&Phidata$Age==0]=mwt
#' Phidata$weight[Phidata$sex=="F"&Phidata$Age==1]=c(0,fwt[1:(length(fwt)-1)])
#' Phidata$weight[Phidata$sex=="M"&Phidata$Age==1]=c(0,mwt[1:(length(mwt)-1)])
#' Phidata$pup=ifelse(Phidata$Age==0,1,0)
#' Phidata$yearling=ifelse(Phidata$Age==1,1,0)
#' Phidata$twoplus=ifelse(Phidata$Age>=2,1,0)
#'# merge Phidata with SST data
#' library(CIPinnipedAnalysis)
#' lastyear=2014
#' source(file.path(system.file(package="CIPinnipedAnalysis"),"CreateAnomalies.r"))
#' env.data=data.frame(cohort=1975:lastyear,ApriltoSeptSSTAnomalies=AprtoSeptAnomalies[4:numyears],OcttoJuneSSTAnomalies=OcttoJuneAnomalies[4:numyears],JulytoJuneSSTAnomalies=JulytoJuneAnomalies[4:numyears],
#' 		MEI=LaggedMEIAprtoSept[-1],MEI1=LaggedMEIOcttoFeb[-1],MEI2=LaggedMEIJunetoFeb[-1],UWI33=UWImeansAprtoSept[1,-(1:6)],UWI36=UWImeansAprtoSept[2,-(1:6)],
#' 		UWI331=UWImeansOcttoFeb[1,-(1:6)],UWI361=UWImeansOcttoFeb[2,-(1:6)],UWI332=UWImeansJunetoFeb[1,-(1:6)],UWI362=UWImeansJunetoFeb[2,-(1:6)],
#' 		SLH=AprtoSeptSLH,SLH1=OcttoFebSLH)
#' Phidata=merge(Phidata,env.data,by.x="Time",by.y="cohort")
#' Phidata=Phidata[order(Phidata$Time,Phidata$Age),]
#' Phidata$Time=Phidata$Time-1989
#' save(Phidata,file="Phidata.rda")

NULL
