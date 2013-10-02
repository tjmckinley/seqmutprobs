#include "functions.h"
#include <R.h>

//recursive function
void genfullmodels_recur_exact(int nhier, int currhier1, int *temp, int *tempind, int q, int nrow, int ncol, int *structure, int *r, int *models_num, double *norm_ls, double *norm_s, double *palt_ls, double *palt_s, double lpriornull_ls, double lprioralt_ls, double lpriornull_s, double lprioralt_s, double *lPDM_int_mat, int *ncomb_sub, int *sum, int *temp_index, int *mods, int *inds, int *prec_ind)
{
	/*'nhier' is number of hierarchies
	'currhier1' is current hierarchy
	'temp' and 'tempind' are vectors recording structural information
	'q' is specific row of 'structure' we wish to extract
	'nrow' and 'ncol' are using for indexing
	'structure' is a matrix containing intermediate model structures
	'r' denotes current model
	'models_num' is a vector recording model structure
	'norm_ls' and 'norm_s' record normalising constants
	'palt_ls' and 'palt_s' record PPAs based on criteria
	'lprior...' are the log-priors for a single model for the two criteria
	'lPDM_int_mat' is intermediate vector for calculating normalising constants
	'sum', 'temp_index', 'mods' and 'inds' are all intermediate vectors for use during calculation*/
		
	//declare necessary variables
	int j,k,currhier,nsamp,nmod=10,priorind;
	double temp_norm,temp_norm1;
	nsamp=ncol;
	currhier=currhier1+1;
	//continue working through hierarchies in order and expand model set
	if(temp[currhier]==1)
	{
		for(j=0;j<nmod;j++)
		{
			models_num[tempind[currhier]]=j;
			if(currhier<(nhier-1)) genfullmodels_recur_exact(nhier,currhier,temp,tempind,q,nrow,ncol,structure,r,models_num,norm_ls,norm_s,palt_ls,palt_s,lpriornull_ls,lprioralt_ls,lpriornull_s,lprioralt_s,lPDM_int_mat,ncomb_sub,sum,temp_index,mods,inds,prec_ind);
			else
			{
				*r=*r+1;
				//now calculate prior specification based on model
				temp_norm=calc_lPDMs_single_vec(nsamp,nmod,models_num,lPDM_int_mat,ncomb_sub,sum,temp_index,mods,inds);
				priorind=calc_lessstring_single(nsamp,models_num);
				temp_norm1=exp(temp_norm+(priorind==1 ? lprioralt_ls:lpriornull_ls));
				if(temp_norm1==0.0) *prec_ind=1;
				*norm_ls=*norm_ls+temp_norm1;
				*palt_ls=*palt_ls+(priorind==1 ? temp_norm1:0.0);
				priorind=calc_string_single(nsamp,models_num);
				temp_norm1=exp(temp_norm+(priorind==1 ? lprioralt_s:lpriornull_s));
				if(temp_norm1==0.0) *prec_ind=1;
				*norm_s=*norm_s+temp_norm1;
				*palt_s=*palt_s+(priorind==1 ? temp_norm1:0.0);
			}
		}
	}
	else
	{
		for(j=0;j<nmod;j++)
		{
			//now fill in gaps and pass to recursive function if necessary
			for(k=0;k<nsamp;k++)
			{
				if(structure[index2(q,k,nrow)]==tempind[currhier])
				{
					models_num[k]=j;
					models_num[k+nsamp]=tempind[currhier];
				}
			}
			if(currhier<(nhier-1)) genfullmodels_recur_exact(nhier,currhier,temp,tempind,q,nrow,ncol,structure,r,models_num,norm_ls,norm_s,palt_ls,palt_s,lpriornull_ls,lprioralt_ls,lpriornull_s,lprioralt_s,lPDM_int_mat,ncomb_sub,sum,temp_index,mods,inds,prec_ind);
			else
			{
				*r=*r+1;
				//now calculate prior specification based on model
				temp_norm=calc_lPDMs_single_vec(nsamp,nmod,models_num,lPDM_int_mat,ncomb_sub,sum,temp_index,mods,inds);
				priorind=calc_lessstring_single(nsamp,models_num);
				temp_norm1=exp(temp_norm+(priorind==1 ? lprioralt_ls:lpriornull_ls));
				if(temp_norm1==0.0) *prec_ind=1;
				*norm_ls=*norm_ls+temp_norm1;
				*palt_ls=*palt_ls+(priorind==1 ? temp_norm1:0.0);
				priorind=calc_string_single(nsamp,models_num);
				temp_norm1=exp(temp_norm+(priorind==1 ? lprioralt_s:lpriornull_s));
				if(temp_norm1==0.0) *prec_ind=1;
				*norm_s=*norm_s+temp_norm1;
				*palt_s=*palt_s+(priorind==1 ? temp_norm1:0.0);
			}
		}
	}
	return;
}

