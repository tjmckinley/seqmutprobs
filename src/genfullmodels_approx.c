#include "functions.h"
#include <R.h>

//function to generate top subset of fully and patially dependent models based on an intermediate set of structures
void genfullmodels_approx(int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int nmodcol, double *lPDM_int_mat, int *lPDM_int_ind,  int ntotcol, int **models_num, int **hyp, double **lPPA_mat, double logthresh, int *ncomb_sub)
{
	/*'q' is specific row of 'structure' we wish to extract
	'nrow' and 'ncol' are using for indexing (based on 'structure')
	'structure' is a matrix containing intermediate model structures
	'r' denotes current model in 'models_num'
	'totmods' is maximum number of models in 'models_num'
	'incmods' is arbitrary number of models to increase if reallocation necessary
	'nmodcol' is number of columns of 'models_num'
	'lPDM_int_mat' and 'lPDM_int_ind' are intermediate matrices (in vector form)
		containing log[P'(D|M)] information
	'ntotcol' is number of columns of 'lPDM_int_ind' and 'lPDM_int_mat'
	'models_num' is a matrix (in vector form) recording final models
	'hyp' is a binary vector of length 'incmods', with 0=null and 1=alt
	'lPPA_mat' is an output vector of log[P'(M|D)]s corresponding to the set of output models
	'logthresh' is the log-threshold relative to the maximum model
	'ncomb_sub' is intermediate vector used for indexing*/
	
	//declare necessary variables
	int i,j,k,nhier,nsamp,nmod=10,valid;
	nsamp=ncol;
	int * temp = (int *) Calloc(nsamp,int);
	int * tempind = (int *) Calloc(nsamp,int);
	int * temp_index = (int *) Calloc(nsamp,int);
	int * currmods1 = (int *) Calloc(nmodcol,int);
	//count up how many hierarchies
	for(i=0;i<nsamp;i++)
	{
		temp[i]=0;
		for(j=0;j<nsamp;j++) temp[i]+=(structure[index2(q,j,nrow)]==i ? 1:0);
	}
	nhier=0;
	for(j=0;j<nsamp;j++) nhier+=(temp[j]>0 ? 1:0);
	//order hierarchies
	bubble_sort_dec_int(temp,tempind,nsamp);
	
	//work out variables relating to maximum independent model
	double temp_lPPA,max_lPPA,max_indep_lPPA;
	max_lPPA=(*lPPA_mat)[0];
	max_indep_lPPA=(*lPPA_mat)[0];
	int corrcol;
	
	//now work through hierarchies in order and expand model set
	if(nhier==1)
	{
		corrcol=ntotcol-1;
		i=0;
		valid=0;
/*		Rprintf("r=%d\n",*r);		*/
		while(i<nmod&&valid==0)
		{
/*			Rprintf("i=%d\t%f\t%d\n",i,lPDM_int_mat[index2(i,corrcol,nmod)],lPDM_int_ind[index2(i,corrcol,nmod)]);*/
/*			if(lPDM_int_ind[index2(i,corrcol,nmod)]>0)*/
/*			{*/
/*				Rprintf("i=%d\n",i);*/
				for(j=0;j<nsamp;j++)
				{
					(*models_num)[index2_col(*r,j,nmodcol)]=lPDM_int_ind[index2(i,corrcol,nmod)];
					(*models_num)[index2_col(*r,j+nsamp,nmodcol)]=0;
				}
				//change model
				temp_lPPA=lPDM_int_mat[index2(i,corrcol,nmod)];
				//check if new model is selected
				if((max_lPPA-temp_lPPA)<logthresh)
				{
					(*lPPA_mat)[*r]=temp_lPPA;
					*r=*r+1;
					//increase size of output vectors if necessary
					if(*r>(*totmods-1)) realloc_approx(nmodcol,totmods,incmods,models_num,hyp,lPPA_mat);
					/*if new model is maximum model then reset maximum model*/
					if(temp_lPPA>max_lPPA) max_lPPA=temp_lPPA;
				}
				else valid=1;
/*			}*/
			i++;
		}
/*		Rprintf("r=%d\n",*r);*/
/*		getchar();*/
	}
	else
	{
		//calculate correct column of lPDM_int_mat based on structure of current hierarchy
		corrcol=calc_corr_col(q,nrow,structure,nsamp,ncomb_sub,tempind[0],temp[0],temp_index);
		//calculate PPAs
		for(i=0;i<nmod;i++)
		{
/*			if(lPDM_int_ind[index2(i,corrcol,nmod)]>0)*/
/*			{*/
				//initialise structure to best independent model
				for(k=0;k<nmodcol;k++) (*models_num)[index2_col(*r,k,nmodcol)]=(*models_num)[index2_col(0,k,nmodcol)];
				//now fill in gaps and initialise calculation of changes to lPDM
				temp_lPPA=max_indep_lPPA;
				for(k=0;k<nsamp;k++)
				{
					if(structure[index2(q,k,nrow)]==tempind[0])
					{
						(*models_num)[index2_col(*r,k,nmodcol)]=lPDM_int_ind[index2(i,corrcol,nmod)];
						(*models_num)[index2_col(*r,k+nsamp,nmodcol)]=tempind[0];
						//adjust lPDM to correct structure
						temp_lPPA-=lPDM_int_mat[index2(0,k,nmod)];
					}
				}
				//adjust lPDM to correct structure
				temp_lPPA+=lPDM_int_mat[index2(i,corrcol,nmod)];
				//enter into recursion
				for(j=0;j<nmodcol;j++) currmods1[j]=(*models_num)[index2_col(*r,j,nmodcol)];
				//enter further hierarchy
				genfullmodels_approx_recur(nhier,0,temp,tempind,currmods1,q,nrow,ncol,structure,r,totmods,incmods,nmodcol,lPDM_int_mat,lPDM_int_ind,ntotcol,models_num,hyp,lPPA_mat,logthresh,ncomb_sub,temp_lPPA,&max_lPPA);
/*			}*/
		}
	}
	//free memory from the heap
	Free(temp);Free(tempind);Free(temp_index);Free(currmods1);
	return;
}

