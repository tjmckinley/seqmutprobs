#include "functions.h"
#include <R.h>

/*function to run recursive loop to calculate all models with fully
independent structures that lie within threshold of maximum model -
to be called from "calc_PPAs_approx_fn"*/
void calc_approx_indep_recur(int nsamp, int nmod, int nmodcol, int incmods, int *totmods, int currsamp, int *r, int *currmods, int *currlocs, int *maxmods, int *maxlocs, double *max_lPPA, double curr_lPPA, double thresh, int **models_num, int **hyp, double *lPDM_int_mat, int *lPDM_int_ind, double **lPPA_mat)
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
	'r' denotes current number of models in model set
	'currmods' is model structure at the previous level of the recursion
	'currlocs' are row locations of 'lPDM_int_mat' corresponding to current model
	'maxmods' is model structure of the maximum model
	'maxlocs' are row locations of 'lPDM_int_mat' corresponding to maximum model
	'max_lPPA' is the log[P'(D|M)] for the current maximum model
	'curr_lPPA' is the log[P'(D|M)] for the current model at the previous
		level of the hierarchy
	'thresh' is the log-threshold value for acceptance relative to the maximum model
	'models_num' is matrix for storing model choices
	'hyp' is a binary vector of length 'incmods', with 0=null and 1=alt
	'lPDM_int_mat' is a matrix of log[P'(D|M)] values for use in intermediary calculations
		for calculating independent model structures
	'lPDM_int_ind' is a binary matrix of indicators denoting which model/sample combinations
		we can explore in the recursion
	'lPPA_mat' is an output vector of log[P'(M|D)]s corresponding to the set of output models*/
	
	/*DECLARE VARIABLES*/
	int i,j,valid;
	double temp_lPPA;
	int * currmods1 = (int *) Calloc(nsamp,int);
	int * currlocs1 = (int *) Calloc(nsamp,int);
	
	//run recursion through to calculate maximum fully independent model
	i=0;	
	valid=1;
	while(i<10&&valid==1)
	{
		//update model
		for(j=0;j<nsamp;j++) (*models_num)[index2_col(*r,j,nmodcol)]=currmods[j];
		(*models_num)[index2_col(*r,currsamp,nmodcol)]=lPDM_int_ind[index2(i,currsamp,nmod)];
		if((*models_num)[index2_col(*r,currsamp,nmodcol)]==currmods[currsamp])
		{
			//enter into recursion
			for(j=0;j<nsamp;j++)
			{
				currmods1[j]=currmods[j];
				currlocs1[j]=currlocs[j];
			}
			currlocs1[currsamp]=i;
			if((currsamp+1)<nsamp) calc_approx_indep_recur(nsamp,nmod,nmodcol,incmods,totmods,currsamp+1,r,currmods1,currlocs1,maxmods,maxlocs,max_lPPA,curr_lPPA,thresh,models_num,hyp,lPDM_int_mat,lPDM_int_ind,lPPA_mat);
			i++;
		}
		else
		{
			//calculate new log-PDM
			temp_lPPA=curr_lPPA-lPDM_int_mat[index2(currlocs[currsamp],currsamp,nmod)]+lPDM_int_mat[index2(i,currsamp,nmod)];
			//add to lPPA_mat output vector
			(*lPPA_mat)[*r]=temp_lPPA;
			//if new model is maximum, then update maximum model
			if(temp_lPPA<=(*max_lPPA))
			{
				//if new model is accepted then update model set
				if((*max_lPPA-temp_lPPA)<thresh)
				{
					//increase model counter
					*r=*r+1;
					//increase size of output vectors if necessary
					if(*r>(*totmods-1)) realloc_approx(nmodcol,totmods,incmods,models_num,hyp,lPPA_mat);
					//enter into recursion
					for(j=0;j<nsamp;j++)
					{
						currmods1[j]=(*models_num)[index2_col(*r-1,j,nmodcol)];
						currlocs1[j]=currlocs[j];
					}
					currlocs1[currsamp]=i;
					if((currsamp+1)<nsamp) calc_approx_indep_recur(nsamp,nmod,nmodcol,incmods,totmods,currsamp+1,r,currmods1,currlocs1,maxmods,maxlocs,max_lPPA,temp_lPPA,thresh,models_num,hyp,lPDM_int_mat,lPDM_int_ind,lPPA_mat);
				}
				else valid=0;
			}
			else
			{
				//increase model counter
				*r=*r+1;
				//increase size of output vectors if necessary
				if(*r>(*totmods-1)) realloc_approx(nmodcol,totmods,incmods,models_num,hyp,lPPA_mat);
				//update maximum model
				*max_lPPA=temp_lPPA;
				//set up vectors to record current maximum model plus locations
				maxmods[currsamp]=(*models_num)[index2_col(*r-1,currsamp,nmodcol)];
				maxlocs[currsamp]=i;
				//enter into recursion
				for(j=0;j<nsamp;j++)
				{
					currmods[j]=currmods[j];
					currlocs1[j]=currlocs[j];
				}
				currlocs1[currsamp]=i;
				if((currsamp+1)<nsamp) calc_approx_indep_recur(nsamp,nmod,nmodcol,incmods,totmods,currsamp+1,r,currmods1,currlocs1,maxmods,maxlocs,max_lPPA,temp_lPPA,thresh,models_num,hyp,lPDM_int_mat,lPDM_int_ind,lPPA_mat);
			}
			i++;
		}
	}
	//free temporary memory from the heap
	Free(currmods1);Free(currlocs1);
	return;
}

