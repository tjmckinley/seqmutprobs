#print method for "mutPPAs.list" objects
print.mutPPAs.list<-function(x,thresh=0.5,digits=2, ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="mutPPAs.list") stop("'x' is not a 'mutPPAs.list' object")
	sapply(x,function(x) if(class(x)!="mutPPAs") stop("Elements of 'x' are not 'mutPPAs' objects"))
	print(summary(x,thresh=thresh,digits=digits, ...))
}

