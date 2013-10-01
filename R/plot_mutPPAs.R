#plot method for "mutPPAs" objects
plot.mutPPAs<-function(x,thresh=0.5,digits=2,prior=NULL,entropy=c("max","mean"), ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="mutPPAs") stop("'x' is not a 'mutPPAs' object")
	plot(summary.mutPPAs(x,thresh=thresh,digits=digits),prior=prior,entropy=entropy)
}

