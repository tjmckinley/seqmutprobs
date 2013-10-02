#include "functions.h"
#include <R.h>

/*function to run recursive loop to calculate final lPPAs for Occam's window approach without recording all models*/
void final_lPPA_recur(int nsamp, int nmod, int nmodcol, int incmods, int *totmods, int currsamp, int *r1, int m1, int *currmods, double curr_lPPA, int curr_hyp, int **models_num, int **hyp, double **lPPA_mat, int **multfact, int *uni, int *uni_ind, int string, double lpriornull, double lprioralt)
{
	/*'nsamp' is the number of samples
	'nmod' is the number of models
	'nmodcol' is the number of samples plus model indicators (the number of columns of 'models_num')
	'incmods' is INITIAL arbitrary number of models corresponding to the amount of memory 
		made available for models_num 
	'totmods' is CURRENT arbitrary number of models corresponding to the amount of memory 
		made available for models_num (may need to be extended during run or
		cut down at the end)
	'currsamp' denotes the sample (i.e. column in 'lPDM_int_mat') that was changed
		at the previous level of the recursion
	'r1' denotes current number of models in adjusted model set
	'm1' is current model in model set that is being compared to
	'currmods' is model structure at the previous level of the recursion
	'curr_lPPA' is the log[P'(D|M)] for the current model at the previous
		level of the hierarchy
	'curr_hyp' is the null/alt status of the model at the previous
		level of the hierarchy
	'models_num' is matrix for storing model choices
	'hyp' is a binary vector of length 'incmods', with 0=null and 1=alt
	'lPPA_mat' is an output vector of log[P'(M|D)]s corresponding to the set of output models
	'multfact' is a vector of multiplication factors
	'uni' and 'uni_ind' are intermediate matrices for use in calculations
	'string' is binary recording criteria type
	'lprior*' represent the prior information for null and alternative models*/
		
	int i,j,k,l,q,d,temp_hyp;
	double temp_lPPA;
	int *currmods1 = (int *) Calloc(nmodcol,int);
	for(i=0;i<nmodcol;i++) currmods1[i]=currmods[i];
	
	//start loop
	for(i=currsamp;i<nsamp;i++)
	{
		j=currmods[i];
		if(uni_ind[index2(j,i,nmod)]>0)
		{
			for(k=0;k<uni_ind[index2(j,i,nmod)];k++)
			{
				q=0;
				for(l=0;l<nsamp;l++) if(currmods[l+nsamp]==currmods[i+nsamp]) q++;
				if(q==1)
				{
					for(l=0;l<nmodcol;l++) currmods1[l]=currmods[l];
					currmods1[i]=uni[index3(j,k,i,nmod,nmod)];
					temp_lPPA=curr_lPPA;
					if(string==0) temp_hyp=calc_lessstring_single(nsamp,currmods1);
					else temp_hyp=calc_string_single(nsamp,currmods1);
					if(temp_hyp!=curr_hyp)
					{
						if(temp_hyp==1) temp_lPPA+=(lprioralt-lpriornull);
						else temp_lPPA+=(lpriornull-lprioralt);
						//update model set
						for(l=0;l<nmodcol;l++) (*models_num)[index2_col(*r1,l,nmodcol)]=currmods1[l];
						(*lPPA_mat)[*r1]=temp_lPPA;
						(*hyp)[*r1]=temp_hyp;
						(*multfact)[*r1]=1;
						*r1=*r1+1;
						//increase size of output vectors if necessary
						if(*r1>(*totmods-1)) realloc_approx_new(nmodcol,totmods,incmods,models_num,hyp,lPPA_mat,multfact);
						//now enter further recursion if necessary
						if((i+1)<nsamp)
						{
							//just check next sample is independent
							q=i+1;
							d=2;
							while(q<nsamp&&d!=1)
							{
								d=0;
								for(l=0;l<nsamp;l++) if(currmods1[l+nsamp]==currmods1[q+nsamp]) d++;
								if(d==1) final_lPPA_recur(nsamp,nmod,nmodcol,incmods,totmods,q,r1,*r1,currmods1,temp_lPPA,temp_hyp,models_num,hyp,lPPA_mat,multfact,uni,uni_ind,string,lpriornull,lprioralt);
								q++;
							}
						}
					}
					else
					{
						(*multfact)[m1]++;
						//now enter further recursion if necessary
						if((i+1)<nsamp)
						{
							//just check next sample is independent
							q=i+1;
							d=2;
							while(q<nsamp&&d!=1)
							{
								d=0;
								for(l=0;l<nsamp;l++) if(currmods1[l+nsamp]==currmods1[q+nsamp]) d++;
								if(d==1) final_lPPA_recur(nsamp,nmod,nmodcol,incmods,totmods,q,r1,m1,currmods1,temp_lPPA,temp_hyp,models_num,hyp,lPPA_mat,multfact,uni,uni_ind,string,lpriornull,lprioralt);
								q++;
							}
						}
					}
				}
			}
		}
	}
	//free memory from the heap
	Free(currmods1);
	return;
}

