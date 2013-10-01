#print method for "mutPPAs" objects
print.mutPPAs<-function(x,thresh=0.5,digits=2, ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="mutPPAs") stop("'x' is not a 'mutPPAs' object")
	print(summary.mutPPAs(x,thresh=thresh,digits=digits, ...))
}

