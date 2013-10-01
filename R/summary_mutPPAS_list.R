#summary method for "mutPPAs.list" objects
summary.mutPPAs.list<-function(object,thresh=0.5,digits=2, ...)
{
	if(missing(object)) stop("'object' argument missing")
	if(class(object)!="mutPPAs.list") stop("'object' is not a 'mutPPAs.list' object")
	sapply(object,function(x) if(class(x)!="mutPPAs") stop("Elements of 'x' are not 'mutPPAs' objects"))
	object.out<-list(NULL)
	for(i in 1:length(object)) object.out[[i]]<-summary(object[[i]],thresh=thresh,digits=digits, ...)
	class(object.out)<-"summary.mutPPAs.list"
	names(object.out)<-names(object)
	object.out
}

