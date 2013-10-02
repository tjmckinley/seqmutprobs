#include "functions.h"
#include <R.h>

//function for calculating log[P'(D|M)]s for all unique models for a given set of base distributions
void calc_lPDM_fn(int *nsamp1, int *totmods1, int *models_num, double *lPDM_int_mat, double *lPDM_mat)
{
	int i,nsamp,totmods,nmod=10;
	nsamp=*nsamp1;
	totmods=*totmods1;
	//generate vector containing number of combinations in each subset
	//(for generating column indices for use in model simplification)
	int *ncomb_sub = (int *) R_alloc(nsamp,sizeof(int));
	ncomb_sub[0]=0;
	for(i=1;i<nsamp;i++) ncomb_sub[i]=choose_c(nsamp,i);
	for(i=1;i<nsamp;i++) ncomb_sub[i]=ncomb_sub[i]+ncomb_sub[i-1];
	//now declare intermediate vectors for model simplification
	int *sum = (int *) R_alloc(nsamp,sizeof(int));
	int *temp_index = (int *) R_alloc(nsamp,sizeof(int));
	int *mods = (int *) R_alloc(nsamp,sizeof(int));
	int *inds = (int *) R_alloc(nsamp,sizeof(int));
	
	//calculate log[P'(D|M)]s
	for(i=0;i<totmods;i++) lPDM_mat[i]=calc_lPDMs_single(i,nsamp,totmods,nmod,models_num,lPDM_int_mat,ncomb_sub,sum,temp_index,mods,inds);
	return;
}

