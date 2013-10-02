#include "functions.h"
#include <R.h>

//function to calculate intermediate values for log[P'(D|M)]
void calc_lPDM_int_fn(int *nsamp1, int *z, double *pstar1, double *lPDM_int_mat)
{	
	/*'nsamp1' is the number of samples
	'z' is a vector of length 4 giving the distribution of bases
		(such that the consensus base is at the last position)
	'pstar1' is the overall mutation rate
	'lPDM_int_mat' is a matrix (in vector form) recording the values of log[P'(D|M)]
		for different models and sample combinations*/
	int i,j,k,l,nsamp,nmod=10;
	double pstar;
	nsamp=*nsamp1;
	pstar=*pstar1;
	//set up vectors to produce lexicographic ordering of columns
	int *ncomb_sub = (int *) R_alloc(nsamp,sizeof(int));
	ncomb_sub[0]=0;
	for(i=1;i<nsamp;i++) ncomb_sub[i]=choose_c(nsamp,i);
	for(i=1;i<nsamp;i++) ncomb_sub[i]=ncomb_sub[i]+ncomb_sub[i-1];
/*	for(i=0;i<nsamp;i++) Rprintf("%d\t",ncomb_sub[i]);*/
/*	Rprintf("\n\n");*/
		
	//generate matrix containing log[P'(D|M)] values for each model/column combination
	int *comb = (int *) Calloc(nsamp,int);
	int *z1 = (int *) Calloc(5,int);
	for(i=0;i<nsamp;i++)
	{
		for(j=0;j<4;j++) z1[j]=z[index2(j,i,4)];
		z1[4]=0;
		for(j=0;j<3;j++) z1[4]+=z[index2(j,i,4)];
		for(j=0;j<nmod;j++) lPDM_int_mat[index2(j,i,nmod)]=lPDM_mod_fn(z1,j,pstar);
	}
	
	//set initial column indictator
	l=nsamp;
	for(i=2;i<=nsamp;i++)
	{
		//initialise combinations
		for(j=0;j<i;j++) comb[j]=j;
		//generate corresponding summed columns
		for(k=0;k<4;k++)
		{
			z1[k]=0;
			for(j=0;j<i;j++) z1[k]+=z[index2(k,comb[j],4)];
		}
		z1[4]=0;
		for(k=0;k<3;k++) z1[4]+=z1[k];
		//calculate log[P'(D|M)] values
		for(j=0;j<nmod;j++) lPDM_int_mat[index2(j,l,nmod)]=lPDM_mod_fn(z1,j,pstar);
		l++;
		//repeat until all combinations at this level are done
		while(next_comb(comb,i,nsamp)==1)
		{
			//generate corresponding summed columns
			for(k=0;k<4;k++)
			{
				z1[k]=0;
				for(j=0;j<i;j++) z1[k]+=z[index2(k,comb[j],4)];
			}
			z1[4]=0;
			for(k=0;k<3;k++) z1[4]+=z1[k];
			//calculate log[P'(D|M)] values
			for(j=0;j<nmod;j++) lPDM_int_mat[index2(j,l,nmod)]=lPDM_mod_fn(z1,j,pstar);
			l++;
		}
	}
	//now sort columns into ascending lexicographic ordering
	int *temp_order,*temp_order1,*temp_ind,*temp_ind1,*temp_val,*temp_samp;
	int nstruct=ncomb_sub[nsamp-1]+1;
	int *overall_order = (int *) Calloc (nstruct,int);
	for(i=0;i<nstruct;i++) overall_order[i]=i;
	int norder;//,temp_norder;
	for(i=2;i<nsamp;i++)
	{
		norder=(ncomb_sub[i]-ncomb_sub[i-1]);
		temp_order = (int *) Calloc(norder*i,int);
		temp_order1 = (int *) Calloc(norder*i,int);
		temp_ind = (int *) Calloc(norder,int);
		temp_ind1 = (int *) Calloc(i,int);
		temp_val = (int *) Calloc(norder,int);
		temp_samp = (int *) Calloc(nsamp,int);
		for(j=0;j<norder;j++) temp_ind[j]=j;
		k=0;
		//initialise combinations
		for(j=0;j<i;j++) comb[j]=j;
		for(j=0;j<i;j++) temp_order[index2_col(k,j,i)]=comb[j];
		while(next_comb(comb,i,nsamp)==1)
		{
			k++;
			for(j=0;j<i;j++) temp_order[index2_col(k,j,i)]=comb[j];
		}	
		//now order each row
		for(j=0;j<norder;j++) bubble_sort_dec_int(&temp_order[index2_col(j,0,i)],temp_ind1,i);
		//store in different format in order to bubble sort each column
		for(j=0;j<norder;j++) for(l=0;l<i;l++) temp_order1[index2(j,l,norder)]=temp_order[index2_col(j,l,i)];
		//now sort according to the first column
		bubble_sort_inc_int(&temp_order1[index2(0,0,norder)],temp_ind,norder);
		for(l=1;l<i;l++)
		{
			for(k=0;k<norder;k++) temp_val[k]=temp_order1[index2(k,l,norder)];
			for(k=0;k<norder;k++) temp_order1[index2(k,l,norder)]=temp_val[temp_ind[k]];
		}
		//adjust overall ordering
		for(k=0;k<norder;k++) temp_val[k]=overall_order[ncomb_sub[i-1]+k];
		for(k=0;k<norder;k++) overall_order[ncomb_sub[i-1]+k]=temp_val[temp_ind[k]];
		//now generate conditions for sequential sorting at next level
		for(k=0;k<nsamp;k++)
		{
			temp_samp[k]=0;
			for(j=0;j<norder;j++) temp_samp[k]+=(temp_order1[index2(j,0,norder)]==k ? 1:0);
		}
		//convert to cumulative counts
		for(k=1;k<nsamp;k++) temp_samp[k]+=temp_samp[k-1];
		lexico_recur(nsamp,1,i,temp_samp,temp_order1,norder,0,overall_order,ncomb_sub);
/*		Rprintf("\n");*/
/*		//print to screen as check		*/
/*		for(l=0;l<norder;l++)*/
/*		{*/
/*			for(j=0;j<i;j++) Rprintf("%d\t",temp_order1[index2(l,j,norder)]);*/
/*			Rprintf("\n");*/
/*		}*/
/*		Rprintf("\n");*/
		Free(temp_order);Free(temp_order1);Free(temp_ind);Free(temp_ind1);Free(temp_val);Free(temp_samp);
		temp_order=NULL;temp_order1=NULL;temp_ind=NULL;temp_ind1=NULL;temp_val=NULL;temp_samp=NULL;
	}
/*	Rprintf("\n");*/
/*	for(l=0;l<nstruct;l++) Rprintf("%d\n",overall_order[l]);*/
/*	getchar();*/
	//now order lPDM_int_mat into correct column ordering
	double * temp_lPDM = (double *) Calloc(nmod*nstruct,double);
	for(i=0;i<nmod;i++) for(j=0;j<nstruct;j++) temp_lPDM[index2(i,j,nmod)]=lPDM_int_mat[index2(i,j,nmod)];
	for(i=0;i<nmod;i++) for(j=0;j<nstruct;j++) lPDM_int_mat[index2(i,j,nmod)]=temp_lPDM[index2(i,overall_order[j],nmod)];
	//free memory from the heap (automatically sets pointer to NULL)
	Free(comb);Free(z1);Free(overall_order);Free(temp_lPDM);
	return;
}

