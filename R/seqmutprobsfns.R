#R script containing functions to run screening process
#for searching for mutations in sequence data using the
#method of McKinley et al., PLoS Comp Biol, (2011) and
#McKinley et al. (2012), in preparation.

#generate consensus sequence from a set of samples
calc_cons<-function(seq_matrix,reference)
{
	cons<-sapply(as.list(1:ncol(seq_matrix)),function(i,x,ref)
	{
		x<-x[,i]
		y<-table(x)
		y<-y[!is.na(y)]
		if(length(y)!=0)
		{
			x<-names(y)[y==max(y)]
			if(length(x)>1)
			{
				if(is.na(match(ref[i],x)))
				{
					cat(paste("Can't generate a consensus for site,",i,"\nand can't match to reference in a sensible manner,\nso this site has been removed:\nREF = ",ref[i],"\n"))
					y<-matrix(y[1:4],ncol=1)
					rownames(y)<-c("a","c","g","t")
					cat(print(y))
					cat("\n")
					x<-NA
				}
				else x<-ref[i]
			}
		}
		else x<-NA
		x
	},x=seq_matrix,ref=reference)
	cons
}

#calculate distribution of bases and sort so that consensus base is in position 4
calc_basedist<-function(seqs,cons)
{
	basedist<-apply(seqs,2,function(x){x<-x[!is.na(x)];c(length(x[x==1]),length(x[x==2]),length(x[x==3]),length(x[x==4]))})
	basedist<-rbind(basedist,cons)
	rows<-1:4
	basedist<-apply(basedist,2,function(x,rows)
	{
		cons<-x[length(x)]
		x<-x[1:(length(x)-1)]
		x<-x[c(rows[-cons],cons)]
		x
	},rows=rows)
	names(basedist)<-1:ncol(basedist)
	basedist
}

#function to turn numerical representation of a single model
#to character form for printing
numtochar<-function(model,nsamp)
{
	#extract models and model indicators
	inds<-model[(nsamp+1):(2*nsamp)]
	#convert any multiple M0s to correct indicator for printing
	mods<-model[1:nsamp]
#	inds[mods==0]<-inds[mods==0][1]
	#convert to character
	mods_char<-as.character(mods)
	#generate counts from indictators
	counts<-table(inds)
	#produce character representation of models
	if(length(counts)!=nsamp)
	{
		counts<-as.numeric(names(counts[counts>1]))
		suffix<-letters[1:length(counts)]
		if(length(suffix[is.na(suffix)])>0) stop("Too many combinations to produce character suffixes for models.")
		for(i in 1:length(counts)) mods_char[inds==counts[i]]<-paste(mods_char[inds==counts[i]],suffix[i],sep="")
		mods_char<-paste(mods_char,collapse=" ")
	}
	else mods_char<-paste(mods_char,collapse=" ")
	#return character representation
	mods_char
}

#function to convert alignment objects into pile-up tables
convert_align_pile<-function(seqdata, pstar=NULL, sites=NA, samp_names, genes, reference)
{
	#'seqdata' must be a list of 'alignment' objects from the 'seqinr' package
	if(!is.list(seqdata)) stop("'seqdata' not a list of alignment objects")
	if(sum(sapply(lapply(seqdata,class),function(x) ifelse(x=="alignment",0,1)))>0) stop("'seqdata' not a list of alignment objects")
	
	nsamp<-length(seqdata)
	if(nsamp==1) stop("Need more than one sample")
		
	#extract number of sequences for summary
	nseq<-sapply(seqdata,function(x) x$nb)
	
	#now convert all alignments into character matrix format
	seqdata<-lapply(seqdata,function(x) 
	{
		x<-as.matrix(x)
		x<-tolower(x)
	})
	#extract number of nucleotide sites for summary
	nnuc<-sapply(seqdata,function(x) ncol(x))
	
	#check that the number of nucleotide sites is the same across all samples
	if(length(unique(nnuc))!=1) stop("Different samples do not contain the same number of nucleotide sites")
	nnuc<-nnuc[1]
	
	#calculate consensus from initial sample
	cons<-calc_cons(seqdata[[1]],reference)
	#check for dodgy sites
	cons_tab<-names(table(cons))
	if(length(cons_tab[cons_tab!="a"&cons_tab!="c"&cons_tab!="g"&cons_tab!="t"&cons_tab!="-"])>0) stop("A nucleotide site exists in the initial sample that does not have a valid consensus (e.g. none of 'a', 'c', 'g', 't' or '-')")
		
	#check for and remove insertions if necessary
	ins_loc<-grep("-",cons)
	full_cons<-cons
	if(length(ins_loc)>0) 
	{
#		cat(paste("Positions:",paste(ins_loc,collapse=" "),"removed as possible insertion sites (i.e. initial sample consensus of '-')\n"))
		seqdata<-lapply(seqdata,function(x,ins_loc) x[,-ins_loc],ins_loc=ins_loc)
		cons<-cons[-ins_loc]
	}
	else ins_loc<-NA
	
	#now convert to numeric format (s.t. {acgt}->{1234}) for quicker processing
	seqdata<-lapply(seqdata,function(x) t(apply(x,1,s2n))+1)
	cons<-s2n(cons)+1

	#now calculate distribution of bases at each site and bind together
	seqdata<-do.call("rbind",lapply(seqdata,function(x,cons) calc_basedist(x,cons),cons=cons))
	
	#if 'pstar' not specified then calculate
	if(is.null(pstar)) pstar<-sum(seqdata[(1:nrow(seqdata))[(1:nrow(seqdata))%%4!=0],])/sum(seqdata)
	
	#now remove sites that aren't selected for testing
	if(!is.na(sites[1]))
	{
		rem_sites<-(1:length(full_cons))[-sites]
		rems<-(as.numeric(colnames(seqdata)) %in% rem_sites)
		rems<-which(rems==T)
		seqdata<-seqdata[,-rems]
		if(is.null(ncol(seqdata))) seqdata<-matrix(seqdata,ncol=1)
		cons<-cons[-rems]
		if(!is.na(ins_loc[1]))
		{
			rem_sites<-c(rem_sites,ins_loc)
			rem_sites<-unique(rem_sites)
		}
	}
	else
	{
		if(!is.na(ins_loc[1])) rem_sites<-ins_loc
		else rem_sites<-NA
	}
	
	#remove duplicate columns and create indicator to reduce memory requirements
	seqdata_temp<-seqdata
	seqdata<-seqdata[,!duplicated(seqdata,MARGIN=2)]
	if(is.null(ncol(seqdata))) seqdata<-matrix(seqdata,ncol=1)
	colnames(seqdata)<-1:ncol(seqdata)
	seqind_temp<-apply(seqdata_temp,2,function(x,y) (1:ncol(y))[apply(y,2,function(y,x) all(x==y),x=x)],y=seqdata)
	#add back in insertion information so that 'seqind' is the same length as the original sequence
	seqind<-numeric(nnuc)
	if(!is.na(rem_sites[1]))
	{
		seqind[rem_sites]<-NA
		seqind[-rem_sites]<-seqind_temp
	}
	else seqind<-seqind_temp
	rm(seqdata_temp,seqind_temp)
	
	#output data
	seqs<-list(seqdata=seqdata,seqind=seqind,pstar=pstar,rem_sites=rem_sites,full_cons=full_cons,nnuc=nnuc,nseq=nseq,nsamp=nsamp,samp_names=samp_names,sites=sites,genes=genes)
}

#function to create pileup tables from data files
create.pileup<-function(filenames,format=c("mase","clustal","phylip","fasta","msf","bam","pileup"),pstar=NULL,genes=NA,sites=NA,samp_names=NULL,cov_thresh=5,reference=NULL,ref_file,mc.cores=1,...)
{
	#'filenames' is a list of characters each containing paths to files we wish to import
	#'format' describes the formats of the files we want to import (all files must be the same format)
	#'pstar' is overall mutation probability 
	#'genes' specifies which segments to extract from any BAM files (only relevant for BAM files)
	#'sites' specify the sites to extract from the alignments
	#'samp_names' is vector of sample names
	#'cov_thresh' is the coverage threshold for extracting sites from NGS data
	#'reference' is a either a single character vector or a list containing the reference sequences
	#	for different Gene segments if necessary
	#'ref_file' is path to reference file
	#'mc.cores' is the number of cores to use for parallel processing
		
	#read in sequences
	if(format[1]!="bam" & format[1]!="pileup")
	{
		if(!is.null(pstar[1])&length(pstar)>1) stop("'pstar' specified incorrectly (must be of length 1 or NULL")
		if(!is.character(reference)) stop("'reference' in wrong format")
		#read in a list of alignment objects
		seqs<-lapply(filenames,function(files,format) read.alignment(files,format=format),format=format[1])
		#generate sample names
		if(is.null(samp_names[1]))
		{
			names(seqs)<-unlist(filenames)
			samp_names<-names(seqs)
		}
		else
		{
			if(length(samp_names)!=length(seqs)) stop("'samp_names' do not match 'filenames'")
			else names(seqs)<-samp_names
		}
		seqnames<-lapply(seqs,names)
		names(seqnames)<-samp_names
		#generate Gene segment names
		if(is.na(genes[1])) genes<-"Segment 1"
		else if(length(genes)!=1) stop("'genes' argument not of length 1")		
		#convert list of 'alignment' objects to pile-up tables
		seqs<-convert_align_pile(seqs,pstar,sites,samp_names,genes,reference)
	}
	else
	{
		if(is.na(cov_thresh[1])) stop("No coverage threshold given")
		#generate pile-up tables using "samtools"
		if(format[1]=="bam")
		{
			if(!file.exists(paste(ref_file,".fai",sep=""))) system(paste("samtools faidx",ref_file))
			pathtoperl<-paste(system.file(package = "seqmutprobs"),"/Perl/pileup2csv.pl",sep="")
			seqs<-mclapply(filenames,function(files,ref) read.table(pipe(paste("samtools mpileup -BQ0 -d10000000 -f ",ref," ",files," | perl ",pathtoperl,sep="")),header=TRUE,quote = "\"", sep="\t"),ref=ref_file,mc.cores=mc.cores)
		}
		else
		{
			pathtoperl<-paste(system.file(package = "seqmutprobs"),"/Perl/pileup2csv.pl",sep="")
			seqs<-lapply(filenames,function(files) read.table(pipe(paste("perl",pathtoperl,files,sep=" ")),header=TRUE,quote = "\"", sep="\t"))
		}
		
		if(is.null(samp_names[1]))
		{
			names(seqs)<-unlist(filenames)
			samp_names<-names(seqs)
		}
		else
		{
			if(length(samp_names)!=length(seqs)) stop("'samp_names' do not match 'filenames'")
			else names(seqs)<-samp_names
		}
		seqnames<-lapply(seqs,function(x) as.character(unique(x$NAME)))
		names(seqnames)<-samp_names		
		
		#extract relevant genes from BAM files
		if(is.na(genes[1])) genes<-seqnames[[1]]
		#check 'pstar' is either NULL, length 1 (in which case applied to all segments) or of the same length as 'genes'
		if(!is.null(pstar[1])&length(pstar)>1&length(pstar)!=length(genes)) stop("'pstar' is not the correct length (i.e. NULL, length 1 or the same length as 'genes'")
		
		matches<-lapply(seqnames,function(nam,genes) match(genes,nam),genes=genes)
		#remove any Gene segments not present in initial sample
		matchesini<-genes[!is.na(matches[[1]])]
		if(!is.null(pstar[1])&length(pstar)>1) pstar<-pstar[!is.na(matches[[1]])]
		mismatches<-genes[is.na(matches[[1]])]
		#exit if no matches found in initial sample
		if(length(matchesini)<1) stop("No matching segments in initial sample")
		#otherwise carry on
		if(length(mismatches)>0)
		{
			if(length(mismatches)>2) cat(paste(paste(paste(mismatches[1:(length(mismatches)-1)],collapse=", "),mismatches[length(mismatches)],collapse=" and "),"not present in initial sample and so are removed\n"))
			else cat(paste(paste(mismatches,collapse="and"),"are only present in initial sample and so are removed\n"))
		}
		matches<-lapply(matches,function(x,matchesini) x[!is.na(matchesini)],matchesini=matches[[1]])
		seqnames<-lapply(as.list(1:length(seqnames)),function(i,nam,matches)
		{
			nam<-nam[[i]]
			matches<-matches[[i]]
			nam[matches[!is.na(matches)]]
		},nam=seqnames,matches=matches)
		#check at least two samples for each segment
		matches<-do.call("c",seqnames)
		matches<-table(matches)
		mismatches<-matches[matches<2]
		if(!is.null(pstar[1])&length(pstar)>1) pstar<-pstar[matches[matches<2]]
		matches<-matches[matches>=2]
		if(length(matches)<1) stop("No segments have more than one sample")
		if(length(mismatches)>0)
		{
			if(length(mismatches)>2) cat(paste(paste(paste(names(mismatches[1:(length(mismatches)-1)]),collapse=", "),names(mismatches[length(mismatches)]),collapse=" and "),"are present in less than two samples and so are removed\n"))
			else cat(paste(paste(names(mismatches),collapse="and"),"are present in less than two samples and so are removed\n"))
		}
		seqnames<-lapply(seqnames,function(nam,mismatch) 
		{
			for(i in mismatch) nam<-nam[nam!=i]
			nam
		},mismatch=names(mismatches))
		names(seqnames)<-names(seqs)
		seqdata<-seqs
		rm(seqs)
		
		#remove sites with a low coverage threshold
		seqdata<-lapply(seqdata,function(x,cov_thresh) x[apply(x[,5:8],1,sum)>cov_thresh,],cov_thresh=cov_thresh)
				
		#merge pile-up tables across individuals for each gene segment
		genes<-seqnames[[1]]
		seqs<-list(length(genes))
		ref_match<-match(genes,names(reference))
		if(length(ref_match[is.na(ref_match)])>0) stop("Reference sequence not available for some gene segments")
		#extract and sort reference sequences into correct order
		reference<-reference[ref_match]
		#now extract information by Gene segment
		for(i in 1:length(genes))
		{
			#extract correct segment for each sample
			seqs[[i]]<-lapply(seqdata,function(seqs,nam) seqs[seqs$NAME==nam,],nam=genes[i])
			#remove any missing samples
			seqs[[i]]<-seqs[[i]][sapply(seqs[[i]],length)>0]
			#figure out consensus and merge pileup tables together
			seqs[[i]]<-lapply(seqs[[i]],function(seqs,nref)
			{
				output<-matrix(NA,4,nref)
				output[,seqs$POS]<-t(seqs[,5:8])
				output
			},nref=length(reference[[i]]))
			cons<-seqs[[i]][[1]]
			cons<-sapply(as.list(1:ncol(cons)),function(i,x)
			{
				ref<-which(c("a","c","g","t")==x[5,i])
				x<-as.numeric(x[1:4,i])
				if(!is.na(sum(x)))
				{
					y<-which(x==max(x))
					if(length(y)>1)
					{
						if(length(!is.na(match(y,ref)))==1) cons<-ref
						else
						{
							cat(paste("Can't generate a consensus for site ",i,",\nand can't match to reference in a sensible manner,\nso this site has been removed:\nREF = ",c("A","C","G","T")[ref],"\n",sep=""))
							cons<-NA
						}
					}
					else cons<-y
				}
				else cons<-NA
				cons		
			},x=rbind(cons,reference[[i]]))
			
			#calculate length of nucleotide sequence and number of samples
			nnuc<-length(reference[[i]])
			nsamp<-length(seqs[[i]])
			
			#merge samples together
			seqs[[i]]<-do.call("rbind",seqs[[i]])
			seq_nam<-1:ncol(seqs[[i]])
			
			#check for and remove invalid sites
			ins_loc<-which(is.na(apply(seqs[[i]],2,sum)) | is.na(cons))
			full_cons<-cons
			full_cons<-n2s(full_cons-1)
			if(length(ins_loc)>0) 
			{
				seqs[[i]]<-seqs[[i]][,-ins_loc]
				seq_nam<-seq_nam[-ins_loc]
				cons<-cons[-ins_loc]
			}
			else ins_loc<-NA
			
			if(length(cons)==0) seqs[[i]]<-NA
			else
			{
				#record range of number of sequences
				nseq<-list(nsamp)
				for(j in 1:nsamp) nseq[[j]]<-range(apply(seqs[[i]][1:4+(j-1)*4,],2,sum))
			
				#reorder rows so that consensus base is in row 4
				seqs[[i]]<-sapply(1:ncol(seqs[[i]]),function(i,x)
				{
					x<-x[,i]
					cons<-x[length(x)]
					x<-x[-length(x)]
					x<-matrix(x,4)
					x<-rbind(x[-cons,],x[cons,])
					x<-as.numeric(x)
					x
				},x=rbind(seqs[[i]],cons))
						
				#if 'pstar' not specified then calculate
				if(is.null(pstar[1])) pstar1<-sum(seqs[[i]][(1:nrow(seqs[[i]]))[(1:nrow(seqs[[i]]))%%4!=0],])/sum(seqs[[i]])
				else
				{
					if(length(pstar)==1) pstar1<-pstar
					else pstar1<-pstar[i]
				}
				#now remove sites that aren't selected for testing
				if(!is.na(sites[1]))
				{
					rem_sites<-(1:length(full_cons))[-sites]
					rems<-(seq_nam %in% rem_sites)
					rems<-which(rems==T)
					seq_nam<-seq_nam[-rems]
					seqs[[i]]<-seqs[[i]][,-rems]
					if(is.null(ncol(seqs[[i]]))) seqs[[i]]<-matrix(seqs[[i]],ncol=1)
					cons<-cons[-rems]
					if(!is.na(ins_loc[1]))
					{
						rem_sites<-c(rem_sites,ins_loc)
						rem_sites<-unique(rem_sites)
					}
					colnames(seqs[[i]])<-seq_nam
				}
				else
				{
					if(!is.na(ins_loc[1])) rem_sites<-ins_loc
					else rem_sites<-NA
				}
	
				#remove duplicate columns and create indicator to reduce memory requirements
				seqdata_temp<-seqs[[i]]
				seqs[[i]]<-seqs[[i]][,!duplicated(seqs[[i]],MARGIN=2)]
				if(is.null(ncol(seqs[[i]]))) seqs[[i]]<-matrix(seqs[[i]],ncol=1)
				colnames(seqs[[i]])<-1:ncol(seqs[[i]])
				seqind_temp<-apply(seqdata_temp,2,function(x,y) (1:ncol(y))[apply(y,2,function(y,x) all(x==y),x=x)],y=seqs[[i]])
				#add back in insertion information so that 'seqind' is the same length as the original sequence
				seqind<-numeric(nnuc)
				if(!is.na(rem_sites[1]))
				{
					seqind[rem_sites]<-NA
					seqind[-rem_sites]<-seqind_temp
				}
				else seqind<-seqind_temp
				rm(seqdata_temp,seqind_temp)
			
				seqs[[i]]<-list(seqdata=seqs[[i]],seqind=seqind,pstar=pstar1,rem_sites=rem_sites,full_cons=full_cons,nnuc=nnuc,nseq=nseq,nsamp=nsamp,sites=sites,samp_names=samp_names,genes=genes[i])
			}
		}
		names(seqs)<-genes
	}
	seqs
}

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
	models_num<-.Call("genmodels",nsamp)
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
		lPDM_int_mat<-.C("calc_lPDM_int_fn",as.integer(nsamp),as.integer(seqs),as.double(pstar),lPDM_int_mat=as.double(numeric(10*ntotcol)))$lPDM_int_mat
		#calculate log[P'(D|M)] values for samples
		lPDM_mat<-.C("calc_lPDM_fn",as.integer(nsamp),as.integer(totmods),as.integer(models_num),as.double(lPDM_int_mat),lPDM_mat=as.double(numeric(totmods)))$lPDM_mat
		lPDM_mat
	},nsamp=nsamp,totmods=totmods,pstar=pstar,ntotcol=ntotcol,models_num=models_num,mc.cores=mc.cores)
	lPDM_mat<-do.call("cbind",lPDM_mat)
	#output list
	output<-list(nseq=nseq,nnuc=nnuc,pstar=pstar,samp_names=samp_names,basedist=seqdata,nucind=seqind,cons=full_cons,lPDM=lPDM_mat,rem_sites=rem_sites,test_sites=sites,models_num=models_num,genes=genes)
	#return list
	output
}

#function for calculating PPAs from output of 'mutlPDMs_full'
mutPPAs_full<-function(mutlPDM,priorPA=c(0.001,0.01,0.05),criteria=c("both","stringent","less"),mc.cores=1, ...)
{	
	#sort priorPA (important for summary)
	priorPA<-sort(priorPA)
	
	#generate null and alternative hypotheses based on model structure and criteria
	if(criteria[1]=="both")
	{
		hyp_lessstring<-.C("calc_lessstring_fn",as.integer(length(mutlPDM$nseq)),as.integer(nrow(mutlPDM$models_num)),as.integer(mutlPDM$models_num),hyp=as.integer(numeric(nrow(mutlPDM$models_num))))$hyp
		hyp_string<-.C("calc_string_fn",as.integer(length(mutlPDM$nseq)),as.integer(nrow(mutlPDM$models_num)),as.integer(mutlPDM$models_num),hyp=as.integer(numeric(nrow(mutlPDM$models_num))))$hyp
	}
	else
	{
		if(criteria[1]=="stringent")
		{
			hyp_lessstring<-NA
			hyp_string<-.C("calc_string_fn",as.integer(length(mutlPDM$nseq)),as.integer(nrow(mutlPDM$models_num)),as.integer(mutlPDM$models_num),hyp=as.integer(numeric(nrow(mutlPDM$models_num))))$hyp
		}
		else
		{
			hyp_lessstring<-.C("calc_lessstring_fn",as.integer(length(mutlPDM$nseq)),as.integer(nrow(mutlPDM$models_num)),as.integer(mutlPDM$models_num),hyp=as.integer(numeric(nrow(mutlPDM$models_num))))$hyp
			hyp_string<-NA
		}
	}
	#produce priorPAs for different criteria
	if(criteria[1]=="both"|criteria[1]=="less")
	{
		#calculate how priorPAs are split across models for less stringent hypothesis
		nalt<-sum(hyp_lessstring)
		nnull<-length(hyp_lessstring)-nalt
		p_lessstring<-apply(matrix(priorPA,1),2,function(p,nalt,nnull,hyp)
		{
			palt<-p/nalt
			pnull<-(1-p)/nnull
			probs<-hyp
			probs[hyp==1]<-palt
			probs[hyp==0]<-pnull
			probs
		},nalt=nalt,nnull=nnull,hyp=hyp_lessstring)
		p_lessstring<-log(p_lessstring)
		p_lessstring<-matrix(p_lessstring,ncol=length(priorPA))
		#calculate PPAs for less stringent hypothesis
		PPAs_lessstring<-mclapply(as.list(1:ncol(p_lessstring)),function(i,pPA,lPDM)
		{
			pPA<-pPA[,i]			
			PPA<-apply(lPDM,2,function(l,p)
			{
				l<-l+p
				ppa<-exp(l)
				if(length(which(!is.finite(log(ppa))))>0)
				{
#					print("Hullo")
					require(Rmpfr)
					prec<-60
					while(length(which(!is.finite(log(ppa))))>0 & prec<=240)
					{
						prec<-prec*2
						ppa<-mpfr(l,prec)
						ppa<-exp(ppa)
						if(length(which(!is.finite(log(ppa))))==0)
						{
							norm<-sum(ppa)
							ppa<-as.numeric(exp(log(ppa)-log(norm)))
							return(list(ppa=ppa,norm=NA))
						}
					}
					stop("Precision issue with normalising constant")
				}
				else
				{
					norm<-sum(ppa)
					ppa<-exp(l-log(norm))
					return(list(ppa=ppa,norm=norm))
				}
			},p=pPA)
			norm<-do.call("c",lapply(PPA,function(ppa) ppa$norm))
			PPA<-do.call("cbind",lapply(PPA,function(ppa) ppa$ppa))
			PPA<-list(norm=norm,PPAs=PPA)
			PPA
		},pPA=p_lessstring,lPDM=mutlPDM$lPDM,mc.cores=mc.cores)
		#retain intermediate prior calculations to include in 'mutPPAs' class
		nalt_less<-nalt
	}
	if(criteria[1]=="both"|criteria[1]=="stringent")
	{
		#calculate how priorPAs are split across models for stringent hypothesis
		nalt<-sum(hyp_string)
		nnull<-length(hyp_string)-nalt
		p_string<-apply(matrix(priorPA,1),2,function(p,nalt,nnull,hyp)
		{
			palt<-p/nalt
			pnull<-(1-p)/nnull
			probs<-hyp
			probs[hyp==1]<-palt
			probs[hyp==0]<-pnull
			probs
		},nalt=nalt,nnull=nnull,hyp=hyp_string)
		p_string<-log(p_string)
		p_string<-matrix(p_string,ncol=length(priorPA))
		#calculate PPAs for stringent hypothesis
		PPAs_string<-mclapply(as.list(1:ncol(p_string)),function(i,pPA,lPDM)
		{
			pPA<-pPA[,i]			
			PPA<-apply(lPDM,2,function(l,p)
			{
				l<-l+p
				ppa<-exp(l)
				if(length(which(!is.finite(log(ppa))))>0)
				{
#					print("Hullo")
					require(Rmpfr)
					prec<-60
					while(length(which(!is.finite(log(ppa))))>0 & prec<=240)
					{
						prec<-prec*2
						ppa<-mpfr(l,prec)
						ppa<-exp(ppa)
						if(length(which(!is.finite(log(ppa))))==0)
						{
							norm<-sum(ppa)
							ppa<-as.numeric(exp(log(ppa)-log(norm)))
							return(list(ppa=ppa,norm=NA))
						}
					}
					stop("Precision issue with normalising constant")
				}
				else
				{
					norm<-sum(ppa)
					ppa<-exp(l-log(norm))
					return(list(ppa=ppa,norm=norm))
				}
			},p=pPA)
			norm<-do.call("c",lapply(PPA,function(ppa) ppa$norm))
			PPA<-do.call("cbind",lapply(PPA,function(ppa) ppa$ppa))
			PPA<-list(norm=norm,PPAs=PPA)
			PPA
		},pPA=p_string,lPDM=mutlPDM$lPDM,mc.cores=mc.cores)
		#retain intermediate prior calculations to include in 'mutPPAs' class
		nalt_string<-nalt
	}
	if(criteria[1]=="less") nalt_string<-NA
	if(criteria[1]=="stringent") nalt_less<-NA
	#condense into correct format
	if(criteria[1]=="both") PPAs<-lapply(as.list(1:length(PPAs_lessstring)),function(i,ppas_ls,ppas_s) list(less=ppas_ls[[i]],stringent=ppas_s[[i]]),ppas_ls=PPAs_lessstring,ppas_s=PPAs_string)
	else
	{
		if(criteria[1]=="less")	PPAs<-lapply(as.list(1:length(PPAs_lessstring)),function(i,ppas_ls) list(less=ppas_ls[[i]],stringent=list(norm=NA,PPAs=NA)),ppas_ls=PPAs_lessstring)
		else PPAs<-lapply(as.list(1:length(PPAs_string)),function(i,ppas_s) list(less=list(norm=NA,PPAs=NA),stringent=ppas_s[[i]]),ppas_s=PPAs_string)
	}
	#output information
	norms<-lapply(PPAs,function(ppas)
	{
		norms<-lapply(ppas,function(ppas) ppas$norm)
		norms
	})
	names(norms)<-paste("priorPA_",priorPA,sep="")
	PPAs<-lapply(PPAs,function(ppas)
	{
		PPAs<-lapply(ppas,function(ppas) ppas$PPAs)
		PPAs
	})
	names(PPAs)<-paste("priorPA_",priorPA,sep="")
	#re-arrange to match output from 'mutPPAs_top'
	PPAs<-lapply(as.list(1:ncol(mutlPDM$basedist)),function(i,ppas)
	{
		lapply(ppas,function(ppas,i)
		{
			lapply(ppas,function(ppas,i)
			{
				if(!is.na(ppas[1])) ppas<-ppas[,i]
				else ppas<-NA
				ppas
			},i=i)
		},i=i)
	},ppas=PPAs)
	names(PPAs)<-colnames(mutlPDM$basedist)
	#re-arrange to match output from 'mutPPAs_top'
	norms<-lapply(as.list(1:ncol(mutlPDM$basedist)),function(i,norms)
	{
		lapply(norms,function(norms,i)
		{
			lapply(norms,function(norms,i)
			{
				if(!is.na(norms[1])) norms<-as.numeric(norms[i])
				else norms<-NA
				norms
			},i=i)
		},i=i)
	},norms=norms)
	names(norms)<-colnames(mutlPDM$basedist)
	hyps<-list(hyp_lessstring,hyp_string)
	nmods<-sapply(hyps,length)
	if(criteria[1]=="less") nmods[2]<-NA
	if(criteria[1]=="stringent") nmods[1]<-NA
	nprior<-c(nmod=nmods[!is.na(nmods)][1],nalt_less=nalt_less,nalt_string=nalt_string)
	PPAs<-list(models_num=mutlPDM$models_num,hyps=hyps,nprior=nprior,PPAs=PPAs,norms=norms)
	names(PPAs$hyps)<-c("less","stringent")
	
	#reset names of sequences
	names(PPAs$PPAs)<-as.character(1:length(PPAs$PPAs))
	names(PPAs$norms)<-as.character(1:length(PPAs$norms))
		
	#output object of class "mutPPAs"
	output<-mutlPDM
	output["lPDM"]<-NULL
	output["models_num"]<-NULL
	output$PPAs<-PPAs
	output$priorPA<-priorPA
	output$estimate<-"full"
	output$supp_output<-FALSE
	output<-output[c(1:8,10,9,11:length(output))]

	#return "mutPPAs" object
	output
}

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
	structure<-.Call("genmodels_priors",nsamp)
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
			lPDM_int_mat<-.C("calc_lPDM_int_fn",as.integer(nsamp),as.integer(seqs),as.double(pstar),lPDM_int_mat=as.double(numeric(10*ntotcol)))$lPDM_int_mat
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
					if(ppas==TRUE) lPPA_mat_ls<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[1],x[2],lPDM_int_mat1,0,structure,length(structure)/nsamp,1,uni,uni_ind)
					else lPPA_mat_ls<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[1],x[2],lPDM_int_mat,0,structure,length(structure)/nsamp,0,uni,uni_ind)
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
								require(Rmpfr)
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
					if(ppas==TRUE) lPPA_mat_s<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[3],x[4],lPDM_int_mat1,1,structure,length(structure)/nsamp,1,uni,uni_ind)
					else lPPA_mat_s<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[3],x[4],lPDM_int_mat,1,structure,length(structure)/nsamp,0,uni,uni_ind)
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
								require(Rmpfr)
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
		lPDM_int_mat<-.C("calc_lPDM_int_fn",as.integer(nsamp),as.integer(seqs),as.double(pstar),lPDM_int_mat=as.double(numeric(10*ntotcol)))$lPDM_int_mat
		#now generate normalising constants
		lnorms<-apply(lpriors,1,function(x,nsamp,ntotcol,lPDM_int_mat) lnorm<-.Call("genmodels_exact",nsamp,ntotcol,x,lPDM_int_mat),nsamp=nsamp,ntotcol=ntotcol,lPDM_int_mat=lPDM_int_mat)
		lnorms
	},nsamp=nsamp,pstar=pstar,ntotcol=ntotcol,lpriors=lpriorPA,mc.cores=mc.cores)
	lnorms<-do.call("cbind",lapply(lnorms,as.numeric))
	#output object
	lnorms
}

#create a generic function and default method for producing PPAs
seqtoPPAs<-function(filenames,ref_file,format=c("fasta","clustal","phylip","mase","msf","bam","pileup"),ref_format=c("fasta","clustal","phylip","mase","msf"),estimate=c("top","full"),criteria=c("both","stringent","less"),supp_output=TRUE,priorPA=c(0.001,0.01,0.05),c=20,samp_names=NULL,pstar=NULL,sites=NA,genes=NA,cov_thresh=5,nswitch_to_supp_output=5, mc.cores=1, ...) UseMethod("seqtoPPAs")

#function that will take sequences directly and output "mutPPAs" object
seqtoPPAs.default<-function(filenames,ref_file,format=c("fasta","clustal","phylip","mase","msf","bam","pileup"),ref_format=c("fasta","clustal","phylip","mase","msf"),estimate=c("top","full"),criteria=c("both","stringent","less"),supp_output=TRUE,priorPA=c(0.001,0.01,0.05),c=20,samp_names=NULL,pstar=NULL,sites=NA,genes=NA,cov_thresh=5,nswitch_to_supp_output=5, mc.cores=1, ...)
{
	#check that inputs are in correct format
	if(missing(filenames)) stop("'filenames' argument missing")
	if(!is.list(filenames)) stop("'filenames' not in list format")
	#check input format
	if(format[1]!="mase"&format[1]!="clustal"&format[1]!="phylip"&format[1]!="fasta"&format[1]!="msf"&format[1]!="bam"&format[1]!="pileup") stop("File for initial sample not one of either 'fasta', 'mase', 'clustal', 'phylip', 'msf', 'bam' or 'pileup' formats")
	format<-format[1]
	sapply(filenames,function(x,format1)
	{
		if(length(x)!=1) stop("Elements of 'filenames' are too long")
		if(!is.character(x)) stop("Elements of 'filenames' are not characters")
		if(!file.exists(x)) stop(paste(x,"does not exist!"))
		#create index files for .bam files if necessary
		if(format1=="bam") if(!file.exists(paste(x,".bai",sep=""))) system(paste("samtools index",x))
	},format1=format)
	if(length(filenames)<=1) stop("Not enough samples to conduct analysis")
	if(missing(ref_file)) stop("'ref_file' argument missing")
	if(!is.character(ref_file)) stop("'ref_file' argument not a character")
	if(length(ref_file)>1) stop("'ref_file' argument must be of length 1")
	if(is.null(samp_names))
	{
		samp_names<-names(filenames)
		if(is.null(samp_names))
		{
			samp_names<-character(length(filenames))
			for(i in 1:length(filenames)) samp_names[i]<-paste("Sample",i)
		}
	}
	if(!is.character(samp_names)) stop("'samp_names' is not a character vector")
	if(length(samp_names)!=length(filenames)) stop("Length of 'samp_names' and 'filenames' don't match")
	
	#check format for reference file
	if(ref_format[1]!="mase"&ref_format[1]!="clustal"&ref_format[1]!="phylip"&ref_format[1]!="fasta"&ref_format[1]!="msf") stop("File for reference sequence not one of either 'fasta', 'mase', 'clustal', 'phylip' or 'msf' formats")
	if((format=="bam" | format=="pileup") & ref_format[1]!="fasta") stop("Reference file needs to be in 'fasta' format when reading in data from BAM files or pileup tables")
	if(!is.na(genes[1])&!is.character(genes)) stop("'genes' in incorrect format")
	if(!is.na(cov_thresh[1])&(!is.numeric(cov_thresh[1])|length(cov_thresh)>1)) stop("'cov_thresh' is not a scalar quantity")
	if(!is.null(pstar[1])&!is.numeric(pstar[1])) stop("'pstar' is not a numeric value")
	if(!is.numeric(priorPA)) stop("'priorPA' not a numeric vector")
	if(!is.numeric(c)|length(c)>1) stop("Wrong input for 'c'")
	if(estimate[1]!="full"&estimate[1]!="top") stop("Wrong value for 'estimate'")
	if(criteria[1]!="both"&criteria[1]!="stringent"&criteria[1]!="less") stop("Wrong value for 'criteria'")
	if(!is.logical(supp_output)) stop("Wrong value for 'supp_output'")
	if(sum(apply(matrix(priorPA,nrow=1),2,function(x) ifelse(x<=0|x>=1,1,0)))>0) stop("'priorPA' incorrectly specified")
	if((!is.na(sites[1])&!is.numeric(sites))|(!is.na(sites[1])&length(sites[sites<0])>0)) stop("'sites' argument incorrect")
	if(!is.numeric(nswitch_to_supp_output)|length(nswitch_to_supp_output)>1) stop("Wrong input for 'nswitch_to_supp_output'")
	if((!is.numeric(mc.cores)&!is.na(mc.cores))|length(mc.cores)>1) stop("Wrong input for 'mc.cores'")
	if(is.na(mc.cores)) mc.cores<-multicore:::detectCores()
	if(!is.finite(mc.cores)) mc.cores<-1
	cat(paste("No. of cores set to",mc.cores,"\n"))
	
	#check whether the number of samples is greater than or equal to 'nswitch_to_supp_output'
	if(length(filenames)>=nswitch_to_supp_output&supp_output==FALSE)
	{
		supp_output<-TRUE
		cat("Number of samples is greater than or equal to 'nswitch_to_supp_output',\nand so individual model outputs are suppressed\n")
	}
	
	#check that the reference file and initial sample matches up
	if(format!="bam" & format!="pileup")
	{
		ini<-ncol(as.matrix(read.alignment(filenames[[1]],format=format)))
		ini_ref<-read.alignment(ref_file,format=ref_format[1])
		if(length(ini_ref$seq)!=1) stop("Reference file contains more than one sequence")
		ini_ref<-strsplit(ini_ref$seq[[1]],"")[[1]]
		ini_ref<-ini_ref[ini_ref!="\r"]
		if(ini!=length(ini_ref)) stop("The length of the sequences in the initial sample and the reference file do not match")
		rm(ini)
	}
	else
	{
		if(format=="bam") ini<-scanBamHeader(filenames[[1]])[[1]]$targets
		else
		{
			pathtoperl<-paste(system.file(package = "seqmutprobs"),"/Perl/pileup2csv.pl",sep="")
			seqs<-mclapply(filenames,function(files) read.table(pipe(paste("perl",pathtoperl,files,sep=" ")),header=TRUE,quote = "\"", sep="\t"),mc.cores=mc.cores)
			ini<-unique(as.character(seqs[[1]]$NAME))
			names(ini)<-ini
		}	
		ini_ref<-read.alignment(ref_file,format=ref_format[1])
		ini_ref_nam<-ini_ref$nam
		#remove excess "\r"s if necessary
		ini_ref_nam<-sapply(as.list(ini_ref_nam),function(x) paste(strsplit(x,"\r")[[1]],collapse=""))
		if(length(ini_ref$seq)!=length(ini)) stop("Reference file does not contain the correct number of Gene segments compared to the initial sample")
		ini_ref<-lapply(ini_ref$seq,function(x)
		{
			x<-strsplit(x,"")[[1]]
			x<-x[x!="\r"]
			x
		})
		#check Gene segments are in initial sample
		names(ini_ref)<-ini_ref_nam
		if(length(which(is.na(match(ini_ref_nam,names(ini)))))!=0) stop("Some Gene segments in initial sample are not in reference file")
#		for(i in 1:length(ini)) if(ini_ref$nam[i]!=names(ini)[i]) stop("Names and/or order of Gene segments in reference file do not match those in the initial sample")
		if(format=="bam") for(i in length(ini)) if(length(ini_ref[[i]])!=ini[i]) stop("The length of the sequences in the initial sample and the reference file do not match")
		rm(ini,ini_ref_nam)
	}
	
	#create pile-up tables dependent on the format of the data
	seqdata<-create.pileup(filenames,format=format,pstar=pstar,genes=genes,sites=sites,samp_names=samp_names,cov_thresh=cov_thresh,reference=ini_ref,ref_file=ref_file,mc.cores=mc.cores)
	if(format!="bam" & format!="pileup") seqdata1<-list(seqdata)
	else
	{
		#check for invalid Gene segments
		seqdata1<-names(seqdata[sapply(seqdata,function(x) ifelse(is.null(names(x)[1]),1,0))==1])
		seqdata<-seqdata[sapply(seqdata,function(x) ifelse(is.null(names(x)[1]),1,0))==0]
		if(length(seqdata)==0) stop("None of the Gene segments contain any valid sites")
		else
		{
			if(length(seqdata1)>0)
			{
				if(length(seqdata1)>2)
				{
					cat(paste("Gene segments ",paste(paste(seqdata1[1:(length(seqdata1)-1)],collapse=", "),seqdata1[length(seqdata1)],collapse=" and ")," do not contain any valid sites\nand have been removed.\n",sep=""))
				}
				else
				{
					if(length(seqdata1)==2) cat(paste("Gene segments ",paste(seqdata1,collapse=" and ")," do not contain any valid sites\nand have been removed.\n",sep=""))
					else cat(paste("Gene segment ",seqdata1," does not contain any valid sites\nand has been removed.\n",sep=""))
				}
			}
		}
		seqdata1<-seqdata
	}
	
	PPAs_output<-list(length(seqdata))
	for(q in 1:length(seqdata1))
	{
		#extract correct element of seqdata
		seqdata<-seqdata1[[q]]
		cat(paste("Analysing segment",names(seqdata1)[q],"\n"))
		#check estimate type to decide subsequent analyses
		if(estimate[1]=="full")
		{
			#check output type to decide subsequent analyses
			if(supp_output==FALSE)
			{
				#calculate log[P'(D|M)]s for all models
				lPDMs<-mutlPDMs_full(seqdata,mc.cores=mc.cores)
				#calculate PPAs based on lPDMs
				PPAs<-mutPPAs_full(lPDMs,priorPA,criteria,mc.cores=mc.cores)		
				rm(lPDMs)
				#calculate final PPAs based on criteria
				PPAs$hyp_PPAs<-lapply(PPAs$PPAs$PPAs,function(ppas,hyps)
				{
					#cycle through priors
					hyp_ppas<-lapply(ppas,function(ppas,hyps)
					{
						#cycle through criteria
						hyp_ppas<-lapply(as.list(1:length(ppas)),function(i,ppas,hyps)
						{
							hyps<-hyps[[i]]
							ppas<-ppas[[i]]
							hyp_ppas<-sum(ppas[hyps==1])
							hyp_ppas
						},ppas=ppas,hyps=hyps)
						names(hyp_ppas)<-c("less","stringent")
						hyp_ppas
					},hyps=hyps)
					hyp_ppas
				},hyps=PPAs$PPAs$hyps)
				#convert to correct format
				#cycle across priors
				PPAs$hyp_PPAs<-lapply(as.list(1:length(priorPA)),function(i,hyps)
				{
					#cycle across criteria
					hyps<-lapply(as.list(1:2),function(j,hyps,i)
					{
						hyps<-sapply(hyps,function(hyps,i,j) hyps[[i]][[j]],i=i,j=j)
						hyps
					},hyps=hyps,i=i)
					names(hyps)<-c("less","stringent")
					hyps
				},hyps=PPAs$hyp_PPAs)
				names(PPAs$hyp_PPAs)<-paste("priorPA_",priorPA,sep="")
				PPAs$hyp_PPAs<-lapply(PPAs$hyp_PPAs,function(hyps) lapply(hyps,function(hyps){names(hyps)<-as.character(1:length(hyps));hyps}))
			}
			else
			{
				#manipulate the data and generate prior information
				PPAs<-mutPPAs_top(seqdata,c,priorPA,criteria,justsetup=TRUE,mc.cores=mc.cores)
				#produce normalising constant by generating (but not recording) each potential model sequentially
				norms<-gen_exact_norm(PPAs$seqdata,PPAs$pstar,PPAs$lpriors,mc.cores=mc.cores)
				if(is.null(ncol(norms))) norms<-matrix(norms,ncol=1)
				
				#extract precision indicators
				precs<-which((1:nrow(norms))%%5==0)
				precs<-norms[precs,]
				if(is.null(ncol(precs))) precs<-matrix(precs,ncol=1)
				precs<-apply(precs,1,function(x) list(which(x==1)))
				precs<-lapply(precs,function(x) x[[1]])
				#remove sites that require multiple precision
#				precs1<-lapply(as.list(1:length(precs)),function(i,x,prior,ppa)
#				{
#					x<-x[[i]]
#					prior<-prior[i]
#					if(length(x)>0)
#					{
#						x<-match(ppa$nucind,x)
#						x<-which(!is.na(x))
##						cat(paste("WARNING: for prior PA = ",prior,", sites ",paste(x,collapse=", "),"\nmight have a precision issue and have NOT been evaluated\n",sep=""))
#					}
#					else x<-NA
#					x<-list(x,prior)
#				},x=precs,prior=priorPA,ppa=PPAs)
				
				#now extract sites that need rerunning and rerun using full model
				#before discarding irrelevant output (INEFFICIENT BUT AVOIDS USING GMP/MPFR IN C)
				if(sum(sapply(precs,length))!=0)
				{
					#extract unique sites that require rerunning
					precs<-do.call("c",precs)
					precs<-unique(precs)
					precs<-sort(precs)
					#create a list containing these sites in order to run full model
					seqdata_temp<-seqdata
					seqdata_temp$seqdata<-seqdata_temp$seqdata[,precs]
					if(length(precs)==1) seqdata_temp$seqdata<-matrix(seqdata_temp$seqdata,ncol=1)
					colnames(seqdata_temp$seqdata)<-1:ncol(seqdata_temp$seqdata)
					seqdata_temp$rem_sites<-(1:seqdata_temp$nnuc)[is.na(match(seqdata_temp$seqind,precs))]
					seqdata_temp$rem_sites<-sort(seqdata_temp$rem_sites)
					seqdata_temp$seqind<-match(seqdata_temp$seqind,precs)
					
					#now generate full models for subset of sites				
					lPDMs<-mutlPDMs_full(seqdata_temp,mc.cores=mc.cores)
					#calculate PPAs based on lPDMs
					PPAs_temp<-mutPPAs_full(lPDMs,priorPA,criteria,mc.cores=mc.cores)		
					rm(lPDMs)
#					#calculate final PPAs based on criteria
					PPAs_temp$hyp_PPAs<-lapply(PPAs_temp$PPAs$PPAs,function(ppas,hyps)
					{
						#cycle through priors
						hyp_ppas<-lapply(ppas,function(ppas,hyps)
						{
							#cycle through criteria
							hyp_ppas<-lapply(as.list(1:length(ppas)),function(i,ppas,hyps)
							{
								hyps<-hyps[[i]]
								ppas<-ppas[[i]]
								hyp_ppas<-sum(ppas[hyps==1])
								hyp_ppas
							},ppas=ppas,hyps=hyps)
							names(hyp_ppas)<-c("less","stringent")
							hyp_ppas
						},hyps=hyps)
						hyp_ppas
					},hyps=PPAs_temp$PPAs$hyps)
					#convert to correct format
					#cycle across priors
					PPAs_temp$hyp_PPAs<-lapply(as.list(1:length(priorPA)),function(i,hyps)
					{
						#cycle across criteria
						hyps<-lapply(as.list(1:2),function(j,hyps,i)
						{
							hyps<-sapply(hyps,function(hyps,i,j) hyps[[i]][[j]],i=i,j=j)
							hyps
						},hyps=hyps,i=i)
						names(hyps)<-c("less","stringent")
						hyps
					},hyps=PPAs_temp$hyp_PPAs)
#					names(PPAs_temp$hyp_PPAs)<-paste("priorPA_",priorPA,sep="")
#					PPAs_temp$hyp_PPAs<-lapply(PPAs$hyp_PPAs,function(hyps) lapply(hyps,function(hyps){names(hyps)<-as.character(1:length(hyps));hyps}))
					PPAs_temp$hyp_PPAs<-do.call("rbind",lapply(PPAs_temp$hyp_PPAs,function(x) do.call("rbind",x)))
					PPAs_temp<-PPAs_temp$hyp_PPAs
					PPAs_temp<-apply(PPAs_temp,2,function(x)
					{
						x<-matrix(x,ncol=length(x)/2)
						x<-rbind(matrix(NA,nrow(x),ncol(x)),x)
						x<-as.numeric(x)
						x
					})
					#now add in the missing PPAs
					rowinds<-(1:nrow(norms))[-which((1:nrow(norms))%%5==0)]
					norms[rowinds,precs]<-PPAs_temp
					rm(PPAs_temp)
				}
				
				#now calculate PPAs for those sites not requiring multiple precision
				norms<-norms[-which((1:nrow(norms))%%5==0),]
				if(is.null(ncol(norms))) norms<-matrix(norms,ncol=1)
				#extract PPAs based on criteria
				hyp_PPAs<-apply(norms,2,function(x,nprior)
				{
					x<-matrix(x,ncol=nprior)
					n<-nrow(x)
					x<-x[-(1:(n/2)),]
					x
				},nprior=length(priorPA))
				#extract norms
				norms<-apply(norms,2,function(x,nprior)
				{
					x<-matrix(x,ncol=nprior)
					n<-nrow(x)
					x<-x[1:(n/2),]
					x
				},nprior=length(priorPA))
				if(criteria[1]=="less") norms[which((1:nrow(norms))%%2==1)+1,]<-NA
				if(criteria[1]=="stringent") norms[which((1:nrow(norms))%%2==1),]<-NA
				if(criteria[1]=="less") hyp_PPAs[which((1:nrow(hyp_PPAs))%%2==1)+1,]<-NA
				if(criteria[1]=="stringent") hyp_PPAs[which((1:nrow(hyp_PPAs))%%2==1),]<-NA
				
#				#adjust to remove those sites with a precision issue
#				for(i in 1:length(precs)) hyp_PPAs[(i-1)*2+(1:2),precs[[i]]]<-NA
#				for(i in 1:length(precs)) norms[(i-1)*2+(1:2),precs[[i]]]<-NA
				
				#convert to correct format for output
				#cycle across priors
				hyp_PPAs<-lapply(as.list(1:length(priorPA)),function(i,hyps)
				{
					#cycle across criteria
					hyps<-lapply(as.list(1:2),function(j,hyps,i)
					{
						hyps<-apply(hyps,2,function(hyps,i,j) hyps[(i-1)*2+j],i=i,j=j)
						hyps
					},hyps=hyps,i=i)
					names(hyps)<-c("less","strict")
					hyps
				},hyps=hyp_PPAs)
				names(hyp_PPAs)<-paste("priorPA_",priorPA,sep="")
				PPAs$PPAs<-list(nprior=PPAs$nprior,norms=norms)
				
				hyp_PPAs<-lapply(hyp_PPAs,function(hyps) lapply(hyps,function(hyps){names(hyps)<-as.character(1:length(hyps));hyps}))

				#remove extraneous outputs and re-arrange to correct format
				PPAs["lpriors"]<-NULL
				PPAs["seqdata"]<-NULL
				PPAs["nprior"]<-NULL
				PPAs$estimate<-estimate[1]
				PPAs$supp_output<-supp_output
				PPAs$hyp_PPAs<-hyp_PPAs
#				PPAs$warning_sites<-precs1
			}
		}
		else
		{
			if(supp_output==FALSE)
			{
				#generate top models
				PPAs<-mutPPAs_top(seqdata,c,priorPA,criteria,mc.cores=mc.cores)
				#produce normalising constants and PPA_SIs
				norms<-mclapply(PPAs$PPAs,function(ppas)
				{
					#cycle through priors
					norms<-lapply(ppas,function(ppas)
					{
						#cycle through criteria
						norms<-lapply(ppas,function(ppas)
						{
							if(!is.na(ppas[1]))
							{
								ppa<-exp(ppas$lPPA)
								if(length(which(!is.finite(log(ppa))))>0)
								{
									norm<-ppa
									prec<-60
									while(length(which(!is.finite(log(norm))))>0 & prec<=240)
									{
										require(Rmpfr)
										prec<-prec*2
										norm<-exp(mpfr(ppas$lPPA,prec))
										if(length(which(!is.finite(log(norm))))==0)
										{
											norm<-sum(norm)
											ppa<-mpfr(ppas$lPPA,prec)
											ppa<-ppa-log(norm)
											ppa<-exp(ppa)
											hyp<-which(ppas$hyps==1)
											if(length(hyp)>0)
											{
												ppa<-sum(ppa[hyp])
												ppa<-as.numeric(ppa)
											}
											else ppa<-as.numeric(0)
										}
									}
									if(length(which(!is.finite(log(norm))))>0) stop("Precision issue with normalising constant")
									norm<-NA
									return(list(norm=norm,hyp_PPA=ppa))
								}
								else
								{
									norm<-sum(ppa)
									hyp_PPA<-exp(ppas$lPPA-log(norm))
									hyp_PPA<-sum(hyp_PPA[ppas$hyps==1])
									return(list(norm=norm,hyp_PPA=hyp_PPA))
								}
							}
							else return(NA)
						})
						norms
					})
					norms
				},mc.cores=mc.cores)
				#extract PPAs in correct order
				hyp_PPAs<-lapply(as.list(1:length(norms[[1]])),function(i,norms)
				{
					#cycle through criteria
					hyp_PPAs<-lapply(as.list(1:length(norms[[1]][[1]])),function(j,norms,i)
					{
						if(!is.na(norms[[1]][[i]][[j]][1]))
						{
							hyp_PPAs<-sapply(as.list(1:length(norms)),function(k,norms,i,j) norms[[k]][[i]][[j]]$hyp_PPA,norms=norms,i=i,j=j)
							names(hyp_PPAs)<-names(norms)
						}
						else hyp_PPAs<-NA
						hyp_PPAs
					},norms=norms,i=i)
					hyp_PPAs
				},norms=norms)
				names(hyp_PPAs)<-names(norms[[1]])
				names(hyp_PPAs[[1]])<-names(norms[[1]][[1]])
				#extract norms in correct order
				norms<-lapply(norms,function(norms)
				{
					#cycle through priors
					norms<-lapply(norms,function(norms)
					{
						#cycle through criteria
						lapply(norms,function(norms)
						{
							if(!is.na(norms[1])) return(norms$norm)
							else return(NA)
						})
					})
					norms
				})
				#extract outputs in correct order
				models_num<-lapply(PPAs$PPAs,function(ppas)
				{
					#cycle through priors
					models<-lapply(ppas,function(ppas)
					{
						#cycle through criteria
						lapply(ppas,function(ppas)
						{
							if(!is.na(ppas[1])) return(ppas$models_num)
							else return(NA)
						})
					})
					models
				})
				hyps<-lapply(PPAs$PPAs,function(ppas)
				{
					#cycle through priors
					hyps<-lapply(ppas,function(ppas)
					{
						#cycle through criteria
						lapply(ppas,function(ppas)
						{
							if(!is.na(ppas[1])) return(ppas$hyps)
							else return(NA)
						})
					})
					hyps
				})
				ppas<-PPAs$PPAs
				ppas<-mclapply(as.list(1:length(ppas)),function(i,ppas,norms)
				{
					#cycle through priors
					ppas<-lapply(as.list(1:length(ppas[[1]])),function(j,ppas,norms,i)
					{
						#cycle through criteria
						lapply(as.list(1:length(ppas[[1]][[1]])),function(k,ppas,norms,i,j)
						{
							if(!is.na(ppas[[i]][[j]][[k]][1]))
							{
								if(is.na(norms[[i]][[j]][[k]]))
								{
									norm<-exp(ppas[[i]][[j]][[k]]$lPPA)
									prec<-60
									while(length(which(!is.finite(log(norm))))>0 & prec<=240)
									{
										require(Rmpfr)
										prec<-prec*2
										norm<-exp(mpfr(ppas[[i]][[j]][[k]]$lPPA,prec))
										if(length(which(!is.finite(log(norm))))==0)
										{
											norm<-sum(norm)
											ppa<-mpfr(ppas[[i]][[j]][[k]]$lPPA,prec)
											ppa<-ppa-log(norm)
											ppa<-as.numeric(exp(ppa))
										}
									}
									if(length(which(!is.finite(log(norm))))>0) stop("Precision issue with normalising constant")
									return(ppa)
								}
								else return(exp(ppas[[i]][[j]][[k]]$lPPA-log(norms[[i]][[j]][[k]])))
							}
							else return(NA)
						},ppas=ppas,norms=norms,i=i,j=j)
					},ppas=ppas,norms=norms,i=i)
					ppas
				},ppas=ppas,norms=norms,mc.cores=mc.cores)
				ppas<-list(models_num=models_num,hyps=hyps,nprior=PPAs$nprior,PPAs=ppas,norms=norms)
		
				#remove extraneous outputs and re-arrange to correct format
				PPAs["lpriors"]<-NULL
				PPAs["seqdata"]<-NULL
				PPAs["nprior"]<-NULL
				PPAs$PPAs<-ppas
				PPAs$estimate<-estimate[1]
				PPAs$supp_output<-supp_output[1]
				PPAs$hyp_PPAs<-hyp_PPAs
			}
			else
			{
				PPAs<-mutPPAs_top(seqdata,c,priorPA,criteria,ppas=FALSE,mc.cores=mc.cores)
				#extract PPAs in correct order - cycle through priors
				hyp_PPAs<-lapply(as.list(1:length(PPAs$PPAs[[1]])),function(i,ppas)
				{
					#cycle through criteria
					ppas<-lapply(as.list(1:length(ppas[[1]][[i]])),function(j,ppas)
					{
						#extract correct element
						ppas<-sapply(ppas,function(ppas,i,j) ppas[[i]][[j]]$lPPA_final,i=i,j=j)
						ppas
					},ppas=ppas)
					ppas
				},ppas=PPAs$PPAs)
				names(hyp_PPAs)<-names(PPAs$PPAs[[1]])
				names(hyp_PPAs[[1]])<-names(PPAs$PPAs[[1]][[1]])
				
				#extract normalising constants in correct order - cycle through priors
				norms<-lapply(as.list(1:length(PPAs$PPAs[[1]])),function(i,ppas)
				{
					#cycle through criteria
					ppas<-lapply(as.list(1:length(ppas[[1]][[i]])),function(j,ppas)
					{
						#extract correct element
						ppas<-sapply(ppas,function(ppas,i,j) ppas[[i]][[j]]$norm,i=i,j=j)
						ppas
					},ppas=ppas)
					ppas
				},ppas=PPAs$PPAs)
				norms<-do.call("rbind",lapply(norms,function(norms) do.call("rbind",norms)))
				colnames(norms)<-names(PPAs$PPAs)				
				
				PPAs$PPAs<-list(nprior=PPAs$nprior,norms=norms)
				
				#remove extraneous outputs and re-arrange to correct format
				PPAs["lpriors"]<-NULL
				PPAs["seqdata"]<-NULL
				PPAs["nprior"]<-NULL
				PPAs$estimate<-estimate[1]
				PPAs$supp_output<-supp_output
				PPAs$hyp_PPAs<-hyp_PPAs
			}
		}
		#output class "mutPPAs" object
		class(PPAs)<-"mutPPAs"
		PPAs_output[[q]]<-PPAs
	}
	if(length(PPAs_output)==1) PPAs_output<-PPAs_output[[1]]
	else
	{
		class(PPAs_output)<-"mutPPAs.list"
		names(PPAs_output)<-names(seqdata1)
	}
	PPAs_output
}

#print method for "mutPPAs.list" objects
print.mutPPAs.list<-function(x,thresh=0.5,digits=2, ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="mutPPAs.list") stop("'x' is not a 'mutPPAs.list' object")
	sapply(x,function(x) if(class(x)!="mutPPAs") stop("Elements of 'x' are not 'mutPPAs' objects"))
	print(summary(x,thresh=thresh,digits=digits, ...))
}

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

#plot method for "mutPPAs.list" objects
plot.mutPPAs.list<-function(x,thresh=0.5,digits=2,prior=NULL,entropy=c("max","mean"),...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="mutPPAs.list") stop("'x' is not a 'mutPPAs.list' object")
	sapply(x,function(x) if(class(x)!="mutPPAs") stop("Elements of 'x' are not 'mutPPAs' objects"))
	plot(summary(x,thresh=thresh,digits=digits),prior=prior,entropy=entropy)
}

#print method for "summary.mutPPAs.list" objects
print.summary.mutPPAs.list<-function(x, ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="summary.mutPPAs.list") stop("'x' is not a 'summary.mutPPAs.list' object")
	for(i in 1:length(x)) print(x[[i]])
}

#plot method for "summary.mutPPAs.list" objects
plot.summary.mutPPAs.list<-function(x,prior=NULL,entropy=c("max","mean"), ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="summary.mutPPAs.list") stop("'x' is not a 'summary.mutPPAs.list' object")
	for(i in 1:length(x)) plot(x[[i]],prior=prior,entropy=entropy)
}

#print method for "mutPPAs" objects
print.mutPPAs<-function(x,thresh=0.5,digits=2, ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="mutPPAs") stop("'x' is not a 'mutPPAs' object")
	print(summary.mutPPAs(x,thresh=thresh,digits=digits, ...))
}

#plot method for "mutPPAs" objects
plot.mutPPAs<-function(x,thresh=0.5,digits=2,prior=NULL,entropy=c("max","mean"), ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="mutPPAs") stop("'x' is not a 'mutPPAs' object")
	plot(summary.mutPPAs(x,thresh=thresh,digits=digits),prior=prior,entropy=entropy)
}

#subset method for "mutPPAs" objects
"[.mutPPAs"<-function(x,i)
{	
	if(missing(i)) return(x)
	if(is.null(i)) return(x)
	subset<-i
	if(all(abs(subset)>0&abs(subset)<length(x$nucind)))
	{
		if(max(table(abs(subset)))>1) stop("Some identical negative and positive elements for subsetting")
		neg<-abs(subset[subset<0])
		pos<-subset[subset>0]
		#remove negative subsets
		if(length(neg)>0)
		{
			x$nucind[neg]<-NA
			x$rem_sites<-c(x$rem_sites,neg)
			x$rem_sites<-sort(unique(x$rem_sites))
		}
		if(length(pos)>0)
		{
			#retain positive sites
			temp<-which(!is.na(match(pos,x$rem_sites)))
			if(length(temp)>0)
			{
				if(length(temp)==1) temp.err<-paste("Site:",pos[temp],"already removed from analysis\n")
				else temp.err<-paste("Sites: ",paste(pos[temp[-length(temp)]],collapse=", "),pos[temp[length(temp)]]," already removed from analysis\n",sep="")
				pos<-pos[-temp]
			}
			if(length(pos)>0)
			{
				x$nucind[-pos]<-NA
				x$rem_sites<-(1:length(x$nucind))[-pos]
				if(exists("temp.err")) cat(temp.err)
			}
			else stop(temp.err)
		}
	}
	else stop("Some elements of 'subset' were not in correct range")
	return(x)
}

#functions for calculating normalised entropy
entropy.fn<-function(z1,z2,...)
{
	s1<-sum(z1)
	s2<-sum(z2)
	log(s1-1)+log(s1-z1[4]+2)+log(s1-z1[4]+1)-log(s2-1)-log(s2-z2[4]+2)-log(s2-z2[4]+1)+lgamma(s1+1)-sum(lgamma(z1+1))-lgamma(s2+1)+sum(lgamma(z2+1))+sum((z1-z2)*digamma(z1+1))-(s1-z1[4]-s2+z2[4])*(digamma(s1-z1[4]+3)-digamma(s1-z1[4]+1))-(s1-s2)*digamma(s1+2)
}

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
	
#print method for "summary.mutPPAs" objects
print.summary.mutPPAs<-function(x, ...)	
{
	#print summary to screen
	cat("\n################################################################################\n")
	if(x$estimate=="full") cat(paste("Results for ",x$genes," based on evaluating ALL models\n",sep=""))
	else cat(paste("Results for ",x$genes," based on evaluating TOP models\n",sep=""))
	cat("################################################################################\n\n")
	cat(paste(length(x$nseq),"sets of samples, containing sequences of length",x$nnuc,"nucleotides:\n"))
	cat("\n")
	if(is.list(x$nseq)) cat(paste("Initial sample (",x$samp_names[1],") ranges between ",x$nseq[[1]][1]," and ",x$nseq[[1]][2] ," sequences",sep=""))
	else cat(paste("Initial sample (",x$samp_names[1],"): ",x$nseq[1]," sequences",sep=""))
	cat("\n")
	for(i in 2:length(x$nseq))
	{
		if(is.list(x$nseq)) cat(paste("Sample ",i-1," (",x$samp_names[i],") ranges between ",x$nseq[[i]][1]," and ",x$nseq[[i]][2] ," sequences",sep=""))
		else cat(paste("Sample ",i-1," (",x$samp_names[i],"): ",x$nseq[i]," sequences",sep=""))
		cat("\n")
	}
	if(!is.list(x$nseq)) cat("\nWarning: some sites might contain less bases due to removal of incomplete data.\nCheck individual sites for details.\n")
	if(!is.na(x$rem_sites[1]))
	{
		test<-x$rem_sites
		if(length(test)>1)
		{
			test1<-test[2:length(test)]-test[1:(length(test)-1)]
			test1<-c(2,test1)
			test2<-NULL
			test3<-numeric(length(test))
			j<-0
			for(i in 1:length(test))
			{
				j<-ifelse(test1[i]==1,j+1,0)
				test3[i]<-j
			}
			test3<-which(test3==0)
			if(length(test3)==1) test2<-paste(min(test),"-",max(test),sep="")
			else
			{
				for(i in 2:length(test3))
				{
					if((test3[i]-test3[i-1])==1) test2<-c(test2,paste(test[test3[i-1]],sep=""))
					else test2<-c(test2,paste(test[test3[i-1]],"-",test[test3[i]-1],sep=""))
					if(i==length(test3))
					{
						if(test3[i]<length(test)) test2<-c(test2,paste(test[test3[i]],"-",test[length(test)],sep=""))
						if(test3[i]==length(test)) test2<-c(test2,paste(test[test3[i-1]],sep=""))
					}
				}
			}
			if(length(test2)>1)
			{
				test2<-as.list(test2)
				if(length(test2)>2) test2[1:(length(test2)-2)]<-lapply(test2[1:(length(test2)-2)],function(x) c(x,", "))
				test2[length(test2)-1][[1]]<-c(test2[length(test2)-1][[1]]," ")
				test2[[length(test2)]]<-c("and ",test2[[length(test2)]][1])
				if(length(test2)>6) test2[(1:length(test2))%%6==0]<-lapply(test2[(1:length(test2))%%6==0],function(x) c(x,"\n"))
				test2<-do.call("c",test2)
				test2<-paste(test2,collapse="")
			}
		}
		else test2<-as.character(test)
		cat(paste("\nPosition(s): ",test2," \nremoved either as possible insertion sites (i.e. initial sample consensus '-')\nor because they have too few samples and/or it is not possible to generate \na consensus base\n",sep=""))
	}
	if(is.na(x$test_sites[1])) cat("\n####### ALL valid sites tested #######\n")
	else
	{
		cat(paste("\n####### SUBSET of valid sites tested #######\n",sep=""))
		cat(paste(paste(x$test_sites,collapse=" "),"\n",sep=""))
	}
	cat(paste("\n",ncol(x$basedist),"/",x$nnuc," unique distributions of bases across the combined set of samples\n",sep=""))
	cat(sprintf("\np* = %.3g\n",x$pstar))
	
	#remove any criteria for which there are no sites of interest
	if(!is.na(x$hyp_id[1]))
	{
		hyp_names_temp<-c("less stringent","stringent")
		hyp_names_temp<-hyp_names_temp[-(x$hyp_id)]
		for(i in hyp_names_temp) cat(paste("\nNo nucleotide sites found that have PPAs > ",x$thresh," for ",i," criterion\n",sep=""))
		if(!is.na(x$sitesofinterest[1]))
		{
			cat(paste("\nAll nucleotide sites that have PPAs>",x$thresh," (for priorPA = ",x$priorPA[1],"):\n",sep=""))
			if(!is.na(x$hyp_output))
			{
				cat("\n############# ")
				cat(x$hyp_output)
				cat(" #############\n\n")
			}
			#print results to screen
			x$sitesofinterest<-rbind(colnames(x$sitesofinterest),x$sitesofinterest)
			x$sitesofinterest<-rbind(c("",x$hyp_names),x$sitesofinterest)
			x$sitesofinterest[is.na(x$sitesofinterest)]<-""
			write.table(format(x$sitesofinterest),quote=F,na="",row.names=F,col.names=F)
			#now print entropy results to the screen
			cat("\n")
			x$entropy<-rbind(colnames(x$entropy),x$entropy)
			x$entropy<-rbind(c("Entropy",rep("",ncol(x$entropy)-1)),x$entropy)
			x$entropy[is.na(x$entropy)]<-""
			x$entropy<-cbind(x$sitesofinterest[,1],x$entropy)
			write.table(format(x$entropy),quote=F,na="",row.names=F,col.names=F)
		}
	}
	else cat(paste("\nNo nucleotide sites found that have PPAs > ",x$thresh," for EITHER criterion\n",sep=""))
}

#plot method for "summary.mutPPAs" objects
plot.summary.mutPPAs<-function(x, prior=NULL, entropy=c("max","mean"), ...)
{
	if(class(x)!="summary.mutPPAs") stop("'x' is not a 'summary.mutPPAs' object")
	if(is.numeric(prior)==FALSE & is.null(prior)==FALSE) stop("'prior' argument is not a numeric scalar or NULL")
	if(length(prior)>1)
	{
		print("'prior' argument has length > 1, and so only first element is used")
		prior<-prior[1]
	}
	if(entropy[1]!="max" & entropy[1]!="mean") stop("Wrong input for 'entropy' argument")
	
	
	sites1<-x$sitesofinterest
	nsites<-nrow(sites1)
	if(nsites==0) stop("No sites-of-interest")
	if(is.null(prior))
	{
		#search for location of smallest prior
		prior<-unique(colnames(sites1)[2:ncol(sites1)])
		prior<-sapply(strsplit(prior,"\\("),function(x) x[[2]])
		prior<-sapply(strsplit(prior,"\\)"),function(x) x[[1]])
		prior1<-min(as.numeric(prior))
		prior<-paste("(",prior1,")",sep="")
		sites1<-sites1[,c(1,which(colnames(sites1)==prior))]
		if(!is.matrix(sites1)) sites1<-matrix(sites1,nrow=1)
	}
	else
	{
		prior1<-prior
		prior<-paste("(",prior,")",sep="")
		sites1<-sites1[,c(1,which(colnames(sites1)==prior))]
		if(!is.matrix(sites1)) sites1<-matrix(sites1,nrow=1)
		if(ncol(sites1)==1) stop("'prior' argument doesn't exist in 'x$sitesofinterest'")
	}
	colnames(sites1)[2:ncol(sites1)]<-x$hyp_names[x$hyp_names!=""]
	#set up colours and point symbols
	cols<-rainbow(n=nsites)
	pchs<-rep(21:25,ceiling(nsites/5))
	pchs<-pchs[1:nsites]
	ymin<-unlist(sites1[,2:ncol(sites1)])
	ymin<-ymin[!is.na(ymin)]
	ymax<-max(ymin)
	ymin<-min(ymin)
	if(entropy[1]=="max") entcol<-ncol(x$entropy)
	else entcol<-ncol(x$entropy)-1
	xmin<-unlist(x$entropy[,entcol])
	xmin<-xmin[!is.na(xmin)]
	xmax<-max(xmin)
	xmin<-min(xmin)
	#plot PPA versus entropy
	
	#calculate outer margin necessary for plot
	ncol.legend<-ceiling(nsites/25)
	
	#open new graphics device and set parameters
	if(dev.cur()==1) dev<-0
	else
	{
		dev<-names(dev.cur())
		if(length(grep("X11",dev))>0 | length(grep("quartz",dev))>0 | length(grep("windows",dev))>0) dev<-0
		else dev<-1
	}
	if(dev==0) dev.new(width=7*(ncol(sites1)-1)+0.5*ncol.legend,height=7)
	par(mfrow=c(1,ncol(sites1)-1),oma=c(0,0,2,4*ncol.legend))
	for(j in 2:ncol(sites1))
	{
		cols1<-cols[!is.na(sites1[,j])]
		pchs1<-pchs[!is.na(sites1[,j])]
		entropy1<-x$entropy[!is.na(sites1[,j]),entcol]
		sites2<-sites1[!is.na(sites1[,j]),]
		if(!is.matrix(sites2)) sites2<-matrix(sites2,nrow=1)
		plot(entropy1,sites2[,j],bg=cols1,pch=pchs1,ylab=paste("PPA (prior PA=",prior1,")",sep=""),xlab=paste(ifelse(entropy[1]=="max","Max.","Mean"),"entropy"),main=paste(colnames(sites2)[j]),ylim=c(ymin,ymax),xlim=c(xmin,xmax))
	}
	legend(par("usr")[2]*1.05,mean(par("usr")[3:4]),legend=sites1[,1],pch=pchs,pt.bg=cols,xpd=NA,yjust=0.5,ncol=ncol.legend)
	mtext(text=paste(x$genes),side=3,font=2,cex=1.5,outer=T)
}

#function to extract and print models relating to specific nucleotides
#to be run on "mutPPAs.list" objects
extract_site_info<-function(object,site,criteria=c("less","stringent"),num_mod=5,digits=2,c=20,draw=FALSE, ...)
{
	#check that inputs are in correct format
	if(missing(object)) stop("'object' argument missing")
	if(missing(site)) stop("'site' argument missing")
	if(class(object)!="mutPPAs") stop("'object' is not a 'mutPPAs' object")
	
	hyp_name<-criteria
	if(!is.character(hyp_name[1])) stop("'hyp_name' is not a character")
	if(hyp_name[1]!="less"&hyp_name[1]!="stringent") stop("'hyp_name is not either 'less' or 'stringent'")
	if(length(site)>1) stop("'site' is not a scalar")
	if(!is.numeric(site)|site<1) stop("'site' argument is incorrect")
	if(!is.numeric(num_mod)|length(num_mod)>1) stop("'num_mod' argument is wrong")
	if(!is.numeric(digits)|length(digits)>1|digits<0) stop("'digits' parameter not in correct format")
	if(!is.logical(draw)) stop("'draw' not a logical")
	
	hyp_name_orig<-c("Less stringent","Stringent")
	
	#for each specified nucleotide site extract and print relevant information
	for(i in site)
	{
		ind<-object$nucind[i]
		if(is.na(ind)) cat(paste("Site",i,"was removed from analysis and/or not evaluated\n"))
		else
		{
			if(object$supp_output==FALSE)
			{
				if(object$estimate=="full")
				{
					if(hyp_name[1]=="less")
					{
						#check whether criteria is missing or not
						if(!is.na(object$PPAs$nprior[2]))
						{
							n_samp_mod<-ncol(object$PPAs$models_num)
							output<-cbind(object$PPAs$models_num,sapply(object$PPAs$PPAs[[ind]],function(ppas,crit) ppas[[crit]],crit=1),object$PPAs$hyps[[1]])
							models<-output[,1:n_samp_mod]
							output<-output[,(n_samp_mod+1):ncol(output)]
							#now generate character representation of models
							models<-apply(models,1,function(mod,nsamp) numtochar(mod,nsamp),nsamp=ncol(models)/2)
							#append models to output data frame
							output<-data.frame(models,output)
							colnames(output)<-c("models",paste("(",object$priorPA,")",sep=""),"Null/Alt")
						}
						else stop(paste("Criteria '",hyp_name[1],"' was not evaluated in this case",sep=""))
					}
					else
					{
						#check whether criteria is missing or not
						if(!is.na(object$PPAs$nprior[3]))
						{
							n_samp_mod<-ncol(object$PPAs$models_num)
							output<-cbind(object$PPAs$models_num,sapply(object$PPAs$PPAs[[ind]],function(ppas,crit) ppas[[crit]],crit=2),object$PPAs$hyps[[2]])
							models<-output[,1:n_samp_mod]
							output<-output[,(n_samp_mod+1):ncol(output)]
							#now generate character representation of models
							models<-apply(models,1,function(mod,nsamp) numtochar(mod,nsamp),nsamp=ncol(models)/2)
							#append models to output data frame
							output<-data.frame(models,output)
							colnames(output)<-c("models",paste("(",object$priorPA,")",sep=""),"Null/Alt")
						}
						else stop(paste("Criteria '",hyp_name[1],"' was not evaluated in this case",sep=""))
					}
				}
				else
				{
					if(hyp_name[1]=="less")
					{
						#check whether criteria is missing or not
						if(!is.na(object$PPAs$nprior[2]))
						{
							#now extract the relevant models from the 'mutPPAs' object
							crit<-1
							n_samp_mod<-ncol(object$PPAs$models_num[[ind]][[1]][[crit]])
							models_num<-lapply(object$PPAs$models_num[[ind]],function(x,crit) x[[crit]],crit=crit)
							hyps<-lapply(object$PPAs$hyps[[ind]],function(x,crit) x[[crit]],crit=crit)
							ppas<-lapply(object$PPAs$PPAs[[ind]],function(x,crit) x[[crit]],crit=crit)
							#sort into the correct order and amalgamate for printing
							models_num<-lapply(models_num,function(models)
							{
								#generate character representation of models
								models<-apply(models,1,function(mod,nsamp) numtochar(mod,nsamp),nsamp=ncol(models)/2)
								models
							})
							if(length(models_num)>1)
							{
								#match up models
								models<-do.call("c",models_num)
								hyps<-do.call("c",hyps)
								hyps<-hyps[!duplicated(models)]
								models<-models[!duplicated(models)]
								models_num<-lapply(models_num,function(mods,compmods) match(compmods,mods),compmods=models)
								#now expand PPAs for each prior separately and fill in the blanks
								ppas<-lapply(as.list(1:length(ppas)),function(i,ppas,orders)
								{
									ppas<-ppas[[i]]
									orders<-orders[[i]]
									ppas<-ppas[orders[!is.na(orders)]]
									orders<-match(orders,orders[!is.na(orders)])
									ppas1<-rep(NA,length(orders))
									ppas1[orders[!is.na(orders)]]<-ppas
									ppas1
								},ppas=ppas,orders=models_num)
								models_num<-models
							}
							#append models to output data frame
							output<-data.frame(models_num,ppas,hyps)
							colnames(output)<-c("models",paste("(",object$priorPA,")",sep=""),"Null/Alt")
						}
						else stop(paste("Criteria '",hyp_name[1],"' was not evaluated in this case",sep=""))
					}
					else
					{
						#check whether criteria is missing or not
						if(!is.na(object$PPAs$nprior[3]))
						{
							#now extract the relevant models from the 'mutPPAs' object
							crit<-2
							n_samp_mod<-ncol(object$PPAs$models_num[[ind]][[1]][[crit]])
							models_num<-lapply(object$PPAs$models_num[[ind]],function(x,crit) x[[crit]],crit=crit)
							hyps<-lapply(object$PPAs$hyps[[ind]],function(x,crit) x[[crit]],crit=crit)
							ppas<-lapply(object$PPAs$PPAs[[ind]],function(x,crit) x[[crit]],crit=crit)
							#sort into the correct order and amalgamate for printing
							models_num<-lapply(models_num,function(models)
							{
								#generate character representation of models
								models<-apply(models,1,function(mod,nsamp) numtochar(mod,nsamp),nsamp=ncol(models)/2)
								models
							})
							if(length(models_num)>1)
							{
								#match up models
								models<-do.call("c",models_num)
								hyps<-do.call("c",hyps)
								hyps<-hyps[!duplicated(models)]
								models<-models[!duplicated(models)]
								models_num<-lapply(models_num,function(mods,compmods) match(compmods,mods),compmods=models)
								#now expand PPAs for each prior separately and fill in the blanks
								ppas<-lapply(as.list(1:length(ppas)),function(i,ppas,orders)
								{
									ppas<-ppas[[i]]
									orders<-orders[[i]]
									ppas<-ppas[orders[!is.na(orders)]]
									orders<-match(orders,orders[!is.na(orders)])
									ppas1<-rep(NA,length(orders))
									ppas1[orders[!is.na(orders)]]<-ppas
									ppas1
								},ppas=ppas,orders=models_num)
								models_num<-models
							}
							#append models to output data frame
							output<-data.frame(models_num,ppas,hyps)
							colnames(output)<-c("models",paste("(",object$priorPA,")",sep=""),"Null/Alt")
						}
						else stop(paste("Criteria '",hyp_name[1],"' was not evaluated in this case",sep=""))
					}
				}
			}
			else
			{
				#produce set of 'top' models for printing
				if(hyp_name[1]=="less")
				{
					#check whether criteria is missing or not
					if(!is.na(object$PPAs$nprior[2]))
					{
						#now search for "top" models
						seqdata<-matrix(object$basedist[,ind],ncol=1)
												
						#calculate models
						nsamp<-nrow(seqdata)/4
						if(nsamp==1) stop("Need more than one sample")
					
						priorPA<-object$priorPA

						#produce prior specifications by generating (but not recording) each potential model sequentially
						structure<-.Call("genmodels_priors",nsamp)
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

						if(hyp_name[1]=="less") nalt_s<-NA
						if(hyp_name[1]=="stringent") nalt_ls<-NA
						nprior<-c(nmod=totmods,nalt_less=nalt_ls,nalt_string=nalt_s)
						#calculate log-threshold
						logc<-log(c)

						#generate total number of columns required to store intermediate calcs
						ntotcol<-sum(choose(nsamp,1:nsamp))
						
						PPA_mat<-apply(seqdata,2,function(seqs,nsamp,pstar,ntotcol,logc,lpriors,structure,criteria)
						{
							#generate (10 x ncol)-matrix of intermediate values for calculating PPAs
							lPDM_int_mat<-.C("calc_lPDM_int_fn",as.integer(nsamp),as.integer(seqs),as.double(pstar),lPDM_int_mat=as.double(numeric(10*ntotcol)))$lPDM_int_mat
							lPPAs<-apply(lpriors,1,function(x,seqs,nsamp,pstar,ntotcol,logc,structure,criteria,lPDM_int_mat)
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
								if(criteria=="both"|criteria=="less")
								{
									#now calculate PPAs according to approximation routine for LESS-STRINGENT criteria
									lPPA_mat_ls<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[1],x[2],lPDM_int_mat1,0,structure,length(structure)/nsamp,1,uni,uni_ind)
									totmods<-lPPA_mat_ls[length(lPPA_mat_ls)]
									lPPA_mat_ls<-lPPA_mat_ls[1:(length(lPPA_mat_ls)-1)]
									models_num<-lPPA_mat_ls[1:(2*nsamp*totmods)]
									models_num<-matrix(models_num,nrow=totmods,byrow=T)
									hyp<-lPPA_mat_ls[(2*nsamp*totmods+1):length(lPPA_mat_ls)]
									lPPA_mat_ls<-hyp[(totmods+1):length(hyp)]
									multfact<-lPPA_mat_ls[(totmods+1):length(lPPA_mat_ls)]
									hyp<-hyp[1:totmods]
									lPPA_mat_ls<-lPPA_mat_ls[1:totmods]
									lPPA_mat_ls<-list(models_num=models_num,hyps=hyp,lPPA=lPPA_mat_ls)
								}
								if(criteria=="both"|criteria=="stringent")
								{
									#now calculate PPAs according to approximation routine for LESS-STRINGENT criteria
									lPPA_mat_s<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[3],x[4],lPDM_int_mat1,1,structure,length(structure)/nsamp,1,uni,uni_ind)
									totmods<-lPPA_mat_s[length(lPPA_mat_s)]
									lPPA_mat_s<-lPPA_mat_s[1:(length(lPPA_mat_s)-1)]
									models_num<-lPPA_mat_s[1:(2*nsamp*totmods)]
									models_num<-matrix(models_num,nrow=totmods,byrow=T)
									hyp<-lPPA_mat_s[(2*nsamp*totmods+1):length(lPPA_mat_s)]
									lPPA_mat_s<-hyp[(totmods+1):length(hyp)]
									multfact<-lPPA_mat_s[(totmods+1):length(lPPA_mat_s)]
									hyp<-hyp[1:totmods]
									lPPA_mat_s<-lPPA_mat_s[1:totmods]
									lPPA_mat_s<-list(models_num=models_num,hyps=hyp,lPPA=lPPA_mat_s)
								}
								if(criteria=="less") lPPA_mat_s<-NA
								if(criteria=="stringent") lPPA_mat_ls<-NA			
								#output lists
								list(less=lPPA_mat_ls,stringent=lPPA_mat_s)
							},seqs=seqs,nsamp=nsamp,pstar=pstar,ntotcol=ntotcol,logc=logc,structure=structure,criteria=criteria,lPDM_int_mat=lPDM_int_mat)
							names(lPPAs)<-paste("priorPA_",priorPA,sep="")
							lPPAs
						},nsamp=nsamp,pstar=object$pstar,ntotcol=ntotcol,logc=logc,lpriors=lpriors,structure=structure,criteria=hyp_name[1])
						#extract correct values from 'norms' if required
						norms<-matrix(object$PPAs$norms[,ind],ncol=1)
						#calculate subset of top models with exact PPAs
						PPAs.adjust<-lapply(as.list(1:ncol(norms)),function(i,ppas,norms,names.priors,estimate,crit)
						{
							#extract correct sequence
							ppas<-ppas[[i]]
							norms<-matrix(norms[,i],nrow=2)
							#cycle through priors
							ppas<-lapply(as.list(1:ncol(norms)),function(j,ppas,norms,estimate,crit)
							{
								ppas<-ppas[[j]]
								norms<-norms[crit,j]
								#cycle through criteria
								if(is.na(norms))
								{
									if(estimate=="top")
									{
										ppa<-exp(ppas[[crit]]$lPPA)
										if(length(which(!is.finite(log(ppa))))>0)
										{
											norm<-ppa
											prec<-60
											while(length(which(!is.finite(log(norm))))>0 & prec<=240)
											{
												require(Rmpfr)
												prec<-prec*2
												norm<-exp(mpfr(ppas[[crit]]$lPPA,prec))
												if(length(which(!is.finite(log(norm))))==0)
												{
													norm<-sum(norm)
													ppa<-mpfr(ppas[[crit]]$lPPA,prec)
													ppa<-ppa-log(norm)
													ppa<-as.numeric(exp(ppa))
												}
											}
											if(length(which(!is.finite(log(norm))))>0) stop("Precision issue with normalising constant")
											ppas[[crit]]$PPAs<-ppa
											ppas[[crit]]$norm<-NA
										}
									}
									else stop("Precision issue with normalising constant")
								}
								else
								{
									ppas[[crit]]$PPAs<-exp(ppas[[crit]]$lPPA-log(norms))
									ppas[[crit]]$norm<-norms
								}
								ppas["lPPA"]<-NULL
								ppas
							},ppas=ppas,norms=norms,estimate=estimate,crit=crit)
							names(ppas)<-names.priors
							ppas				
						},ppas=PPA_mat,norms=norms,names.priors=names(PPA_mat[[1]]),estimate=object$estimate,crit=1)
						PPAs.adjust<-PPAs.adjust[[1]]
						#output results
						crit<-1
						n_samp_mod<-ncol(PPAs.adjust[[1]][[crit]]$models_num)
						models_num<-lapply(PPAs.adjust,function(x,crit) x[[crit]]$models_num,crit=crit)
						hyps<-lapply(PPAs.adjust,function(x,crit) x[[crit]]$hyps,crit=crit)
						ppas<-lapply(PPAs.adjust,function(x,crit) x[[crit]]$PPAs,crit=crit)
						#sort into the correct order and amalgamate for printing
						models_num<-lapply(models_num,function(models)
						{
							#generate character representation of models
							models<-apply(models,1,function(mod,nsamp) numtochar(mod,nsamp),nsamp=ncol(models)/2)
							models
						})
						if(length(models_num)>1)
						{
							#match up models
							models<-do.call("c",models_num)
							hyps<-do.call("c",hyps)
							hyps<-hyps[!duplicated(models)]
							models<-models[!duplicated(models)]
							models_num<-lapply(models_num,function(mods,compmods) match(compmods,mods),compmods=models)
							#now expand PPAs for each prior separately and fill in the blanks
							ppas<-lapply(as.list(1:length(ppas)),function(i,ppas,orders)
							{
								ppas<-ppas[[i]]
								orders<-orders[[i]]
								ppas<-ppas[orders[!is.na(orders)]]
								orders<-match(orders,orders[!is.na(orders)])
								ppas1<-rep(NA,length(orders))
								ppas1[orders[!is.na(orders)]]<-ppas
								ppas1
							},ppas=ppas,orders=models_num)
							models_num<-models
						}
						#append models to output data frame
						output<-data.frame(models_num,ppas,hyps)
						colnames(output)<-c("models",paste("(",object$priorPA,")",sep=""),"Null/Alt")
					}
					else stop(paste("Criteria '",hyp_name[1],"' was not evaluated in this case",sep=""))
				}
				else
				{
					#check whether criteria is missing or not
					if(!is.na(object$PPAs$nprior[3]))
					{
						#now search for "top" models
						seqdata<-matrix(object$basedist[,ind],ncol=1)
												
						#calculate models
						nsamp<-nrow(seqdata)/4
						if(nsamp==1) stop("Need more than one sample")
					
						priorPA<-object$priorPA

						#produce prior specifications by generating (but not recording) each potential model sequentially
						structure<-.Call("genmodels_priors",nsamp)
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

						if(hyp_name[1]=="less") nalt_s<-NA
						if(hyp_name[1]=="stringent") nalt_ls<-NA
						nprior<-c(nmod=totmods,nalt_less=nalt_ls,nalt_string=nalt_s)
						#calculate log-threshold
						logc<-log(c)

						#generate total number of columns required to store intermediate calcs
						ntotcol<-sum(choose(nsamp,1:nsamp))
						
						PPA_mat<-apply(seqdata,2,function(seqs,nsamp,pstar,ntotcol,logc,lpriors,structure,criteria)
						{
							#generate (10 x ncol)-matrix of intermediate values for calculating PPAs
							lPDM_int_mat<-.C("calc_lPDM_int_fn",as.integer(nsamp),as.integer(seqs),as.double(pstar),lPDM_int_mat=as.double(numeric(10*ntotcol)))$lPDM_int_mat
							lPPAs<-apply(lpriors,1,function(x,seqs,nsamp,pstar,ntotcol,logc,structure,criteria,lPDM_int_mat)
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
								if(criteria=="both"|criteria=="less")
								{
									#now calculate PPAs according to approximation routine for LESS-STRINGENT criteria
									lPPA_mat_ls<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[1],x[2],lPDM_int_mat1,0,structure,length(structure)/nsamp,1,uni,uni_ind)
									totmods<-lPPA_mat_ls[length(lPPA_mat_ls)]
									lPPA_mat_ls<-lPPA_mat_ls[1:(length(lPPA_mat_ls)-1)]
									models_num<-lPPA_mat_ls[1:(2*nsamp*totmods)]
									models_num<-matrix(models_num,nrow=totmods,byrow=T)
									hyp<-lPPA_mat_ls[(2*nsamp*totmods+1):length(lPPA_mat_ls)]
									lPPA_mat_ls<-hyp[(totmods+1):length(hyp)]
									multfact<-lPPA_mat_ls[(totmods+1):length(lPPA_mat_ls)]
									hyp<-hyp[1:totmods]
									lPPA_mat_ls<-lPPA_mat_ls[1:totmods]
									lPPA_mat_ls<-list(models_num=models_num,hyps=hyp,lPPA=lPPA_mat_ls)
								}
								if(criteria=="both"|criteria=="stringent")
								{
									#now calculate PPAs according to approximation routine for LESS-STRINGENT criteria
									lPPA_mat_s<-.Call("calc_PPAs_approx_fn",nsamp,ntotcol,logc,x[3],x[4],lPDM_int_mat1,1,structure,length(structure)/nsamp,1,uni,uni_ind)
									totmods<-lPPA_mat_s[length(lPPA_mat_s)]
									lPPA_mat_s<-lPPA_mat_s[1:(length(lPPA_mat_s)-1)]
									models_num<-lPPA_mat_s[1:(2*nsamp*totmods)]
									models_num<-matrix(models_num,nrow=totmods,byrow=T)
									hyp<-lPPA_mat_s[(2*nsamp*totmods+1):length(lPPA_mat_s)]
									lPPA_mat_s<-hyp[(totmods+1):length(hyp)]
									multfact<-lPPA_mat_s[(totmods+1):length(lPPA_mat_s)]
									hyp<-hyp[1:totmods]
									lPPA_mat_s<-lPPA_mat_s[1:totmods]
									lPPA_mat_s<-list(models_num=models_num,hyps=hyp,lPPA=lPPA_mat_s)
								}
								if(criteria=="less") lPPA_mat_s<-NA
								if(criteria=="stringent") lPPA_mat_ls<-NA			
								#output lists
								list(less=lPPA_mat_ls,stringent=lPPA_mat_s)
							},seqs=seqs,nsamp=nsamp,pstar=pstar,ntotcol=ntotcol,logc=logc,structure=structure,criteria=criteria,lPDM_int_mat=lPDM_int_mat)
							names(lPPAs)<-paste("priorPA_",priorPA,sep="")
							lPPAs
						},nsamp=nsamp,pstar=object$pstar,ntotcol=ntotcol,logc=logc,lpriors=lpriors,structure=structure,criteria=hyp_name[1])
						#extract correct values from 'norms' if required
						norms<-matrix(object$PPAs$norms[,ind],ncol=1)
						#calculate subset of top models with exact PPAs
						PPAs.adjust<-lapply(as.list(1:ncol(norms)),function(i,ppas,norms,names.priors,estimate,crit)
						{
							#extract correct sequence
							ppas<-ppas[[i]]
							norms<-matrix(norms[,i],nrow=2)
							#cycle through priors
							ppas<-lapply(as.list(1:ncol(norms)),function(j,ppas,norms,estimate,crit)
							{
								ppas<-ppas[[j]]
								norms<-norms[crit,j]
								#cycle through criteria
								if(is.na(norms))
								{
									if(estimate=="top")
									{
										ppa<-exp(ppas[[crit]]$lPPA)
										if(length(which(!is.finite(log(ppa))))>0)
										{
											norm<-ppa
											prec<-60
											while(length(which(!is.finite(log(norm))))>0 & prec<=240)
											{
												require(Rmpfr)
												prec<-prec*2
												norm<-exp(mpfr(ppas[[crit]]$lPPA,prec))
												if(length(which(!is.finite(log(norm))))==0)
												{
													norm<-sum(norm)
													ppa<-mpfr(ppas[[crit]]$lPPA,prec)
													ppa<-ppa-log(norm)
													ppa<-as.numeric(exp(ppa))
												}
											}
											if(length(which(!is.finite(log(norm))))>0) stop("Precision issue with normalising constant")
											ppas[[crit]]$PPAs<-ppa
											ppas[[crit]]$norm<-NA
										}
									}
									else stop("Precision issue with normalising constant")
								}
								else
								{
									ppas[[crit]]$PPAs<-exp(ppas[[crit]]$lPPA-log(norms))
									ppas[[crit]]$norm<-norms
								}
								ppas["lPPA"]<-NULL
								ppas
							},ppas=ppas,norms=norms,estimate=estimate,crit=crit)
							names(ppas)<-names.priors
							ppas				
						},ppas=PPA_mat,norms=norms,names.priors=names(PPA_mat[[1]]),estimate=object$estimate,crit=2)
						PPAs.adjust<-PPAs.adjust[[1]]
						#output results
						crit<-2
						n_samp_mod<-ncol(PPAs.adjust[[1]][[crit]]$models_num)
						models_num<-lapply(PPAs.adjust,function(x,crit) x[[crit]]$models_num,crit=crit)
						hyps<-lapply(PPAs.adjust,function(x,crit) x[[crit]]$hyps,crit=crit)
						ppas<-lapply(PPAs.adjust,function(x,crit) x[[crit]]$PPAs,crit=crit)
						#sort into the correct order and amalgamate for printing
						models_num<-lapply(models_num,function(models)
						{
							#generate character representation of models
							models<-apply(models,1,function(mod,nsamp) numtochar(mod,nsamp),nsamp=ncol(models)/2)
							models
						})
						if(length(models_num)>1)
						{
							#match up models
							models<-do.call("c",models_num)
							hyps<-do.call("c",hyps)
							hyps<-hyps[!duplicated(models)]
							models<-models[!duplicated(models)]
							models_num<-lapply(models_num,function(mods,compmods) match(compmods,mods),compmods=models)
							#now expand PPAs for each prior separately and fill in the blanks
							ppas<-lapply(as.list(1:length(ppas)),function(i,ppas,orders)
							{
								ppas<-ppas[[i]]
								orders<-orders[[i]]
								ppas<-ppas[orders[!is.na(orders)]]
								orders<-match(orders,orders[!is.na(orders)])
								ppas1<-rep(NA,length(orders))
								ppas1[orders[!is.na(orders)]]<-ppas
								ppas1
							},ppas=ppas,orders=models_num)
							models_num<-models
						}
						#append models to output data frame
						output<-data.frame(models_num,ppas,hyps)
						colnames(output)<-c("models",paste("(",object$priorPA,")",sep=""),"Null/Alt")
					}
					else stop(paste("Criteria '",hyp_name[1],"' was not evaluated in this case",sep=""))
				}
			}
			#check whether the correct number of models can be output (in case using the 'top' routine)
			if(num_mod>nrow(output)) num_mod<-nrow(output)
			#sort into decreasing order according to PPA and extract top models
			output<-output[order(-output[,2],output[,1]),]
			output<-output[1:num_mod,]
			#print base distribution
			cat(paste("\nDistribution of bases at site ",i,":\n",sep=""))
			basedist<-matrix(object$basedist[,ind],4)
			cons<-object$cons[i]
			bases<-c("a","c","g","t")
			bases<-c(bases[bases!=cons],bases[bases==cons])
			rownames(basedist)<-toupper(bases)
			colnames(basedist)<-object$samp_names
			basedistprop<-round(apply(basedist,2,function(x) x/sum(x)),digits=2)
			basedistprop.draw<-basedistprop
			basedistprop<-rbind(colnames(basedist),basedistprop)
			basedist<-rbind(colnames(basedist),basedist)
			write.table(format(basedist,justify="centre"),quote=F,na="",col.names=F)
			#print base distribution as proportions (easier to visualise for NGS data)
			cat(paste("\nDistribution of bases at site",i,"as proportions:\n"))
			write.table(format(basedistprop,justify="centre"),quote=F,na="",col.names=F)
			#print top num_mod models
			if(hyp_name[1]=="less") hyp_name<-"LESS STRINGENT"
			else hyp_name<-"STRINGENT"
			cat(paste("\nTop ",num_mod," models at site ",i," based on '",hyp_name,"' criterion:\n",sep=""))
			output[,2:ncol(output)]<-round(output[,2:ncol(output)],digits=digits)
#			output1<-output
			output<-apply(output,2,as.character)
			output<-rbind(colnames(output),output)
			output<-format(output,justify="centre")
			write.table(output,quote=F,na="",row.names=F,col.names=F)
			cat(paste("\nOverall PPA_SI for ",hyp_name," criterion:\n",sep=""))
			output<-sapply(object$hyp_PPAs,function(ppas,ind,crit)
			{
				if(crit=="LESS STRINGENT") x<-ppas[[1]][[ind]]
				else x<-ppas[[2]][[ind]]
				x
			},ind=ind,crit=hyp_name)
			output<-round(output,digits=digits)
			output<-as.character(output)
			output<-rbind(sapply(strsplit(names(object$hyp_PPAs),"_"),function(x) paste(x[1]," (",x[2],")",sep="")),output)
			output<-format(output,justify="centre")
			write.table(output,quote=F,na="",row.names=F,col.names=F)
			#add entropy calculations
			#calculate normalised entropy measure
			dists<-object$basedist[,ind]
			dists<-matrix(dists,nrow=4)
			entropy<-numeric(ncol(dists)-1)
			#calculate entropy
			for(j in 2:ncol(dists))
			{
				#produce normalised entropy
				s2<-sum(dists[,j])
				s2<-s2*diag(4)
				ans1<-apply(s2,2,function(x,dists) entropy.fn(dists[,1],x),dists=dists)
				entropy[j-1]<-entropy.fn(dists[,1],dists[,j])/max(ans1)
			}
			entropy<-round(entropy,digits=digits)
			entropy<-c(entropy,mean(entropy),max(entropy))
			entropy<-as.character(entropy)
			entropy<-matrix(entropy,nrow=1)
			entropy<-rbind(c(paste(2:(ncol(entropy)-1),":1",sep=""),"Mean","Max"),entropy)
			cat("\nEntropy measures\n")
			entropy<-format(entropy,justify="centre")
			write.table(entropy,quote=F,na="",row.names=F,col.names=F)
			#draw if required
			if(draw==TRUE)
			{
				basecols<-c("red","blue","green","yellow")
				basedistprop.draw<-apply(basedistprop.draw,2,rev)
				basecols<-basecols[match(rownames(basedistprop.draw),c("A","C","G","T"))]
				par(mar=c(3,4,3,3)+0.1,oma=c(0,0,0,2))
				barplot(basedistprop.draw,col=basecols,ylab="Proportion of bases",main=paste("Site",site,"in",object$genes))
				legend(par("usr")[2]*1.12,mean(par("usr")[3:4]),legend=c("A","C","G","T"),fill=c("red","blue","green","yellow"),xjust=1,xpd=T,yjust=0.5)
			}
		}
	}
#	if(!exists("output1")) output1<-NA
#	return(output1)
}
