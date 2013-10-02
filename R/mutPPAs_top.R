#function for selecting top models
mutPPAs_top<-function(seqdata,c=20,priorPA=c(0.001,0.05,0.01),criteria=c("both","stringent","less"),ppas=TRUE,justsetup=FALSE,mc.cores=mc.cores, ...)
{	
	#sort priorPA (important for summary)
	priorPA<-sort(priorPA)
		
	#extract elements from 'seqdata'
	seqind<-seqdata$seqind
	pstar<-seqdata$pstar
	rem_sites<-seqdata$rem_sites
	full_cons<-seqdata$full_cons
	nnuc<-seqdata$nnuc
	nseq<-seqdata$nseq
	nsamp<-seqdata$nsamp
	samp_names<-seqdata$samp_names
	sites<-seqdata$sites
	genes<-seqdata$genes
	seqdata<-seqdata$seqdata
		
	#produce prior specifications by generating (but not recording) each potential model sequentially
	structure<-.Call("genmodels_priors",nsamp, PACKAGE = "seqmutprobs")
	nstructure<-length(structure)
	totmods<-structure[nstructure-2]
	nalt_ls<-structure[nstructure-1]
	nalt_s<-structure[nstructure]
	nnull_ls<-totmods-nalt_ls
	nnull_s<-totmods-nalt_s
	structure<-structure[1:(nstructure-3)]
	#generate priorPAs for different hypotheses
	lpriornull_ls<-(1-priorPA)/nnull_ls
	lprioralt_ls<-priorPA/nalt_ls
	lpriornull_s<-(1-priorPA)/nnull_s
	lprioralt_s<-priorPA/nalt_s
	lpriors<-cbind(lpriornull_ls,lprioralt_ls,lpriornull_s,lprioralt_s)
	lpriors<-log(lpriors)
	
	if(criteria[1]=="less") nalt_s<-NA
	if(criteria[1]=="stringent") nalt_ls<-NA
	nprior<-c(nmod=totmods,nalt_less=nalt_ls,nalt_string=nalt_s)
	#calculate log-threshold
	logc<-log(c)
	
	#generate total number of columns required to store intermediate calcs
	ntotcol<-sum(choose(nsamp,1:nsamp))
		
	#calculate PPAs for each sample
	if(justsetup==FALSE)
	{
		seqdata1<-apply(seqdata,2,function(x) list(x))
		seqdata1<-lapply(seqdata1,function(x) x[[1]])
		PPA_mat<-mclapply(seqdata1,function(seqs,nsamp,pstar,ntotcol,logc,lpriors,structure,criteria,priorPA,ppas)
		{
			#generate (10 x ncol)-matrix of intermediate values for calculating PPAs
			lPDM_int_mat<-.C("calc_lPDM_int_fn",as.integer(nsamp),as.integer(seqs),as.double(pstar),lPDM_int_mat=as.double(numeric(10*ntotcol)), PACKAGE = "seqmutprobs")$lPDM_int_mat
			lPPAs<-apply(lpriors,1,function(x,seqs,nsamp,pstar,ntotcol,logc,structure,criteria,lPDM_int_mat,ppas)
			{
				#sort out and remove duplicates in 'lPDM_int_mat'
				lPDM_int_mat1<-lPDM_int_mat
				lPDM_int_mat<-matrix(lPDM_int_mat,10)
				uni<-matrix(0,10,nsamp*10)
				uni_ind<-matrix(0,10,nsamp)
				for(i in 1:nsamp)
				{
					z<-unique(lPDM_int_mat[duplicated(lPDM_int_mat[,i]),i])
					if(length(z)>0)
					{
						for(j in z)
						{
							y<-which(lPDM_int_mat[,i]==j)
							y1<-y[1]
							y<-y[2:length(y)]
							lPDM_int_mat[y,i]<-(-1e10)
							uni[y1,(i-1)*10+(1:length(y))]<-y-1
							uni_ind[y1,i]<-length(y)
						}
					}
				}
				lPDM_int_mat<-as.numeric(lPDM_int_mat)
				uni<-as.numeric(uni)
				uni_ind<-as.numeric(uni_ind)
				#now calculate model sets
				if(criteria=="both"|criteria=="less")
				{
					#now calculate PPAs according to approximation routine for LESS-STRINGENT criteria
					if(ppas==TRUE) lPPA_mat_ls<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[1],x[2],lPDM_int_mat1,0,structure,length(structure)/nsamp,1,uni,uni_ind, PACKAGE = "seqmutprobs")
					else lPPA_mat_ls<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[1],x[2],lPDM_int_mat,0,structure,length(structure)/nsamp,0,uni,uni_ind, PACKAGE = "seqmutprobs")
					totmods<-lPPA_mat_ls[length(lPPA_mat_ls)]
					lPPA_mat_ls<-lPPA_mat_ls[1:(length(lPPA_mat_ls)-1)]
					models_num<-lPPA_mat_ls[1:(2*nsamp*totmods)]
					models_num<-matrix(models_num,nrow=totmods,byrow=T)
					hyp<-lPPA_mat_ls[(2*nsamp*totmods+1):length(lPPA_mat_ls)]
					lPPA_mat_ls<-hyp[(totmods+1):length(hyp)]
					multfact<-lPPA_mat_ls[(totmods+1):length(lPPA_mat_ls)]
					hyp<-hyp[1:totmods]
					lPPA_mat_ls<-lPPA_mat_ls[1:totmods]
					if(ppas==TRUE) lPPA_mat_ls<-list(models_num=models_num,hyps=hyp,lPPA=lPPA_mat_ls)
					else
					{
						#calculate final lPPAs and normalising constant
						norm<-multfact*exp(lPPA_mat_ls)
						#if necessary then use "Rmpfr" package to calculate result to multiple precision
						if(length(which(!is.finite(log(norm))))>0)
						{
							prec<-60
							while(length(which(!is.finite(log(norm))))>0 & prec<=240)
							{
								prec<-prec*2
								norm<-exp(mpfr(lPPA_mat_ls,prec))
								if(length(which(!is.finite(log(norm))))==0)
								{
									multfact<-mpfr(multfact,prec)
									norm<-sum(multfact*norm)
									lPPA_final<-lPPA_mat_ls[hyp==1]
									if(length(lPPA_final)>0)
									{
										lPPA_final<-mpfr(lPPA_final,prec)
										multfact<-multfact[hyp==1]
										lPPA_final<-lPPA_final-log(norm)
										lPPA_final<-as.numeric(sum(multfact*exp(lPPA_final)))
									}
									else lPPA_final<-as.numeric(0)
								}
							}
							if(length(which(!is.finite(log(norm))))>0) stop("Precision issue with normalising constant")
							norm<-NA
						}
						else
						{
							norm<-sum(norm)	
							#calculate final lPPAs
							lPPA_final<-lPPA_mat_ls[hyp==1]
							multfact<-multfact[hyp==1]
							if(length(lPPA_final)>0)
							{
								lPPA_final<-lPPA_final-log(norm)
								lPPA_final<-sum(multfact*exp(lPPA_final))
							}
							else lPPA_final<-0
						}

						#save out results
						lPPA_mat_ls<-list(lPPA_final=lPPA_final,norm=norm)
					}
				}
				if(criteria=="both"|criteria=="stringent")
				{
					#now calculate PPAs according to approximation routine for LESS-STRINGENT criteria
					if(ppas==TRUE) lPPA_mat_s<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[3],x[4],lPDM_int_mat1,1,structure,length(structure)/nsamp,1,uni,uni_ind, PACKAGE = "seqmutprobs")
					else lPPA_mat_s<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[3],x[4],lPDM_int_mat,1,structure,length(structure)/nsamp,0,uni,uni_ind, PACKAGE = "seqmutprobs")
					totmods<-lPPA_mat_s[length(lPPA_mat_s)]
					lPPA_mat_s<-lPPA_mat_s[1:(length(lPPA_mat_s)-1)]
					models_num<-lPPA_mat_s[1:(2*nsamp*totmods)]
					models_num<-matrix(models_num,nrow=totmods,byrow=T)
					hyp<-lPPA_mat_s[(2*nsamp*totmods+1):length(lPPA_mat_s)]
					lPPA_mat_s<-hyp[(totmods+1):length(hyp)]
					multfact<-lPPA_mat_s[(totmods+1):length(lPPA_mat_s)]
					hyp<-hyp[1:totmods]
					lPPA_mat_s<-lPPA_mat_s[1:totmods]
					if(ppas==TRUE) lPPA_mat_s<-list(models_num=models_num,hyps=hyp,lPPA=lPPA_mat_s)
					else
					{
						#calculate final lPPAs and normalising constant
						norm<-multfact*exp(lPPA_mat_s)
						#if necessary then use "Rmpfr" package to calculate result to multiple precision
						if(length(which(!is.finite(log(norm))))>0)
						{
							prec<-60
							while(length(which(!is.finite(log(norm))))>0 & prec<=240)
							{
								prec<-prec*2
								norm<-exp(mpfr(lPPA_mat_s,prec))
								if(length(which(!is.finite(log(norm))))==0)
								{
									multfact<-mpfr(multfact,prec)
									norm<-sum(multfact*norm)
									lPPA_final<-lPPA_mat_s[hyp==1]
									if(length(lPPA_final)>0)
									{
										lPPA_final<-mpfr(lPPA_final,prec)
										multfact<-multfact[hyp==1]
										lPPA_final<-lPPA_final-log(norm)
										lPPA_final<-as.numeric(sum(multfact*exp(lPPA_final)))
									}
									else lPPA_final<-as.numeric(0)
								}
							}
							if(length(which(!is.finite(log(norm))))>0) stop("Precision issue with normalising constant")
							norm<-NA
						}
						else
						{
							norm<-sum(norm)	
							#calculate final lPPAs
							lPPA_final<-lPPA_mat_s[hyp==1]
							multfact<-multfact[hyp==1]
							if(length(lPPA_final)>0)
							{
								lPPA_final<-lPPA_final-log(norm)
								lPPA_final<-sum(multfact*exp(lPPA_final))
							}
							else lPPA_final<-0
						}

						#save out results
						lPPA_mat_s<-list(lPPA_final=lPPA_final,norm=norm)
					}
				}
				if(criteria=="less") lPPA_mat_s<-NA
				if(criteria=="stringent") lPPA_mat_ls<-NA			
				#output lists
				list(less=lPPA_mat_ls,stringent=lPPA_mat_s)
			},seqs=seqs,nsamp=nsamp,pstar=pstar,ntotcol=ntotcol,logc=logc,structure=structure,criteria=criteria,lPDM_int_mat=lPDM_int_mat,ppas=ppas)
			names(lPPAs)<-paste("priorPA_",priorPA,sep="")
			lPPAs
		},nsamp=nsamp,pstar=pstar,ntotcol=ntotcol,logc=logc,lpriors=lpriors,structure=structure,criteria=criteria[1],priorPA=priorPA,ppas=ppas,mc.cores=mc.cores)
	}
	else PPA_mat<-NA
	
	#output object
	output<-list(nseq=nseq,nnuc=nnuc,pstar=pstar,samp_names=samp_names,basedist=seqdata,nucind=seqind,cons=full_cons,rem_sites=rem_sites,test_sites=sites,priorPA=priorPA,PPAs=PPA_mat,lpriors=lpriors,seqdata=seqdata,nprior=nprior,genes=genes)
	output
}
