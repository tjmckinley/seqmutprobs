#function for calculating log[P'(D|M)]s for full model set
mutlPDMs_full<-function(seqdata,mc.cores=1, ...)
{	
	#extract elements from 'seqdata'
	seqind<-seqdata$seqind
	pstar<-seqdata$pstar
	rem_sites<-seqdata$rem_sites
	full_cons<-seqdata$full_cons
	nnuc<-seqdata$nnuc
	nseq<-seqdata$nseq
	nsamp<-seqdata$nsamp
	sites<-seqdata$sites
	samp_names<-seqdata$samp_names
	genes<-seqdata$genes
	seqdata<-seqdata$seqdata
	
	#generate models and log-likelihoods

	#calculate total number of feasible models
	models_num<-.Call("genmodels",nsamp, PACKAGE = "seqmutprobs")
	models_num<-matrix(models_num,ncol=2*nsamp)
	totmods<-nrow(models_num)
	
	#generate total number of columns required to store intermediate calcs
	ntotcol<-sum(choose(nsamp,1:nsamp))
	
	#calculate log[P'(D|M)] values for each sample
	seqdata1<-apply(seqdata,2,function(x) list(x))
	seqdata1<-lapply(seqdata1,function(x) x[[1]])
	lPDM_mat<-mclapply(seqdata1,function(seqs,nsamp,totmods,pstar,ntotcol,models_num)
	{
		#generate (10 x ncol)-matrix of intermediate values for calculating PPAs
		lPDM_int_mat<-.C("calc_lPDM_int_fn",as.integer(nsamp),as.integer(seqs),as.double(pstar),lPDM_int_mat=as.double(numeric(10*ntotcol)), PACKAGE = "seqmutprobs")$lPDM_int_mat
		#calculate log[P'(D|M)] values for samples
		lPDM_mat<-.C("calc_lPDM_fn",as.integer(nsamp),as.integer(totmods),as.integer(models_num),as.double(lPDM_int_mat),lPDM_mat=as.double(numeric(totmods)), PACKAGE = "seqmutprobs")$lPDM_mat
		lPDM_mat
	},nsamp=nsamp,totmods=totmods,pstar=pstar,ntotcol=ntotcol,models_num=models_num,mc.cores=mc.cores)
	lPDM_mat<-do.call("cbind",lPDM_mat)
	#output list
	output<-list(nseq=nseq,nnuc=nnuc,pstar=pstar,samp_names=samp_names,basedist=seqdata,nucind=seqind,cons=full_cons,lPDM=lPDM_mat,rem_sites=rem_sites,test_sites=sites,models_num=models_num,genes=genes)
	#return list
	output
}



