#create a generic function and default method for producing PPAs
seqtoPPAs<-function(filenames,ref_file,format=c("fasta","clustal","phylip","mase","msf","bam","pileup"),ref_format=c("fasta","clustal","phylip","mase","msf"),estimate=c("top","full"),criteria=c("both","stringent","less"),supp_output=TRUE,priorPA=c(0.001,0.01,0.05),c=20,samp_names=NULL,pstar=NULL,sites=NA,genes=NA,cov_thresh=5,nswitch_to_supp_output=5, mc.cores=1, ...) UseMethod("seqtoPPAs")

