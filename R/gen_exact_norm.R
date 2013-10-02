#function to generate EXACT normalising constants for use in efficient EXACT search routine
gen_exact_norm<-function(seqdata,pstar,lpriorPA,mc.cores=1, ...)
{	
	#calculate models
	nsamp<-nrow(seqdata)/4
	if(nsamp==1) stop("Need more than one sample")
	
	#generate total number of columns required to store intermediate calcs
	ntotcol<-sum(choose(nsamp,1:nsamp))
	#calculate normalising constants for each sample
	seqdata1<-apply(seqdata,2,function(x) list(x))
	seqdata1<-lapply(seqdata1,function(x) x[[1]])
	lnorms<-mclapply(seqdata1,function(seqs,nsamp,pstar,ntotcol,lpriors)
	{
		#generate (10 x ncol)-matrix of intermediate values for calculating PPAs
		lPDM_int_mat<-.C("calc_lPDM_int_fn",as.integer(nsamp),as.integer(seqs),as.double(pstar),lPDM_int_mat=as.double(numeric(10*ntotcol)), PACKAGE = "seqmutprobs")$lPDM_int_mat
		#now generate normalising constants
		lnorms<-apply(lpriors,1,function(x,nsamp,ntotcol,lPDM_int_mat) lnorm<-.Call("genmodels_exact",nsamp,ntotcol,x,lPDM_int_mat, PACKAGE = "seqmutprobs"),nsamp=nsamp,ntotcol=ntotcol,lPDM_int_mat=lPDM_int_mat)
		lnorms
	},nsamp=nsamp,pstar=pstar,ntotcol=ntotcol,lpriors=lpriorPA,mc.cores=mc.cores)
	lnorms<-do.call("cbind",lapply(lnorms,as.numeric))
	#output object
	lnorms
}

