#plot method for "summary.mutPPAs.list" objects
plot.summary.mutPPAs.list<-function(x,prior=NULL,entropy=c("max","mean"), ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="summary.mutPPAs.list") stop("'x' is not a 'summary.mutPPAs.list' object")
	for(i in 1:length(x)) plot(x[[i]],prior=prior,entropy=entropy)
}

