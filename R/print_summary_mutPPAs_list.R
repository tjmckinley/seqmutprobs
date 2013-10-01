#print method for "summary.mutPPAs.list" objects
print.summary.mutPPAs.list<-function(x, ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="summary.mutPPAs.list") stop("'x' is not a 'summary.mutPPAs.list' object")
	for(i in 1:length(x)) print(x[[i]])
}

