#summary method for "mutPPAs" objects
summary.mutPPAs<-function(object,thresh=0.5,digits=2, ...)
{
	#check input
	x<-object
	if(missing(x)) stop("'object' argument missing")
	if(class(x)!="mutPPAs") stop("'object' is not a 'mutPPAs' object")
	if(!is.numeric(thresh)|length(thresh)>1|thresh>1|thresh<0) stop("'thresh' parameter not in correct format")
	if(!is.numeric(digits)|length(digits)>1|digits<0) stop("'digits' parameter not in correct format")
	
	hyps<-x$hyp_PPAs
	#extract criteria by prior
	hyps<-lapply(as.list(1:2),function(i,hyps) lapply(as.list(1:length(hyps)),function(j,hyps,i) hyps[[j]][[i]],hyps=hyps,i=i),hyps=hyps)
	hyps<-lapply(hyps,function(hyps) do.call("cbind",hyps))
	hyps<-lapply(hyps,function(hyps) hyps<-cbind(inds=as.numeric(rownames(hyps)),hyps))
	hyp_output<-NA
	hyp_names_orig<-c("Less stringent","Stringent")
	if(is.na(x$PPAs$nprior[2]))
	{
		hyp_output<-"LESS STRINGENT criterion was not evaluated"
		hyps<-hyps[2]
	}
	if(is.na(x$PPAs$nprior[3]))
	{
		hyp_output<-"STRINGENT criterion was not evaluated"
		hyps<-hyps[1]
	}

	#extract subset of sites by threshold
	sitesofinterest<-lapply(hyps,function(hyps,thresh) hyps[!is.nan(hyps[,2])&!is.na(hyps[,2])&hyps[,2]>=thresh,],thresh=thresh)
	sitesofinterest[sapply(sitesofinterest,length)==0]<-NA

	#set hypothesis ID for printing
	hyp_id<-(1:length(sitesofinterest))
	hyp_names<-NA
	#remove any criteria for which there are no sites of interest
	if(length(sitesofinterest[is.na(sitesofinterest)])>0)
	{
		hyp_id<-hyp_id[!is.na(sitesofinterest)]
		sitesofinterest<-sitesofinterest[!is.na(sitesofinterest)]
	}
	if(length(sitesofinterest)>0)
	{
		nhyp<-length(hyp_id)
		#expand to include all sites
		sitesofinterest<-lapply(sitesofinterest,function(sites,nucind)
		{
			if(is.null(ncol(sites))) sites<-matrix(sites,ncol=length(sites))
			for(i in 1:length(sites[,1]))
			{
				temp<-which(nucind==sites[i,1])
				if(i==1) output<-cbind(temp,matrix(rep(sites[i,2:ncol(sites)],length(temp)),ncol=(ncol(sites)-1),byrow=T))
				else output<-rbind(output,cbind(temp,matrix(rep(sites[i,2:ncol(sites)],length(temp)),ncol=(ncol(sites)-1),byrow=T)))
			}
			output
		},nucind=x$nucind)
		#sort into ascending order by site
		sitesofinterest<-lapply(sitesofinterest,function(x){x<-x[sort.list(x[,1]),];if(!is.matrix(x)) x<-matrix(x,1);x})
		#merge output from different hypotheses together in order to print
		allsites<-matrix(unique(do.call("c",lapply(sitesofinterest,function(x) x[,1]))),nrow=1)
		allsites<-sort(allsites)
		sitesofinterest<-lapply(sitesofinterest,function(sites,unisites)
		{
			output<-matrix(NA,length(unisites),ncol(sites)-1)
			output[unisites %in% sites[,1],]<-sites[,2:ncol(sites)]
			output
		},unisites=allsites)
		sitesofinterest<-do.call("cbind",sitesofinterest)
		colnames(sitesofinterest)<-rep(paste("(",x$priorPA,")",sep=""),nhyp)
		sitesofinterest<-cbind(sites=allsites,round(sitesofinterest,digits=digits))
		nprior<-length(x$priorPA)
		nrowsites<-nrow(sitesofinterest)
		ncolsites<-ncol(sitesofinterest)
	
		#sort and output top models
		sitesofinterest<-sitesofinterest[sort.list(sitesofinterest[,2],decreasing=T),]
		if(!is.matrix(sitesofinterest))
		{
			temp_names<-names(sitesofinterest)
			sitesofinterest<-matrix(sitesofinterest,1)
			colnames(sitesofinterest)<-temp_names
		}
		sitesofinterest[,2:ncol(sitesofinterest)]<-round(sitesofinterest[,2:ncol(sitesofinterest)],digits=digits)
		hyp_names<-character(nprior*length(hyp_id))
		j<-1
		for(i in 1:length(hyp_names))
		{
			if((i-1)%%nprior==0)
			{
				hyp_names[i]<-hyp_names_orig[hyp_id[j]]
				j<-j+1
			}
			else hyp_names[i]<-""
		}
		#calculate normalised entropy measure
		dists<-x$basedist[,x$nucind[sitesofinterest[,1]]]
		if(is.null(ncol(dists))) dists<-matrix(dists,ncol=1)
		entropy<-apply(dists,2,function(x)
		{
			dists<-matrix(x,nrow=4)
			ans<-numeric(ncol(dists)-1)
			#calculate entropy
			for(j in 2:ncol(dists))
			{
				#produce normalised entropy
				s2<-sum(dists[,j])
				s2<-s2*diag(4)
				ans1<-apply(s2,2,function(x,dists) entropy.fn(dists[,1],x),dists=dists)
				ans[j-1]<-entropy.fn(dists[,1],dists[,j])/max(ans1)
			}
			ans
		})
		if(nrow(dists)>8) entropy<-t(entropy)
		else entropy<-matrix(entropy,ncol=1)
		entropy<-cbind(entropy,apply(entropy,1,mean),apply(entropy,1,max))
		entropy<-round(entropy,digits=digits)
		colnames(entropy)<-c(paste(2:(ncol(entropy)-1),":1",sep=""),"Mean","Max")
	}
	else
	{
		sitesofinterest<-NA
		hyp_id<-NA
		entropy<-NA
	}

	#remove elements of x that are no longer necessary
	x$sitesofinterest<-sitesofinterest
	x$entropy<-entropy
	x$hyp_output<-hyp_output
	x$hyp_names<-hyp_names
	x$hyp_id<-hyp_id
	x$thresh<-thresh
	class(x)<-"summary.mutPPAs"
	x
}
