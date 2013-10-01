/*C code for generating models for sequence data, intended
to be called by associated R code*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

//function for calculating integer powers
int pow_int(int base, int exp)
{
	int i,ans;
	if(exp==0) ans=1;
	else
	{
		ans=base;
		for(i=0;i<(exp-1);i++) ans=ans*base;
	}
	return ans;
}

//indexing function for matrices
inline int index2(int i, int j, int nrow)
{
	return (nrow*j)+i;
}

//indexing function for matrices
inline int index2_col(int i, int j, int ncol)
{
	return (ncol*i)+j;
}

//indexing function for matrices
inline int index3(int i, int j, int k, int nrow, int ncol)
{
	return (k*nrow*ncol)+(nrow*j)+i;
}

//log-factorial function
double lfactorial(int n)
{
	int i;
	double ans=0.0;
	for(i=n;i>0;i--) ans+=log(i);
	return ans;
}

//factorial function
int factorial(int n)
{
	int i;
	int ans=1;
	for(i=n;i>0;i--) ans=ans*i;
	return ans;
}

//choose function
int choose_c(int n, int k)
{
	int ans;
	if(n<k) ans=0;
	else ans=factorial(n)/(factorial(k)*factorial(n-k));
	return ans;
}

/*This code to calculate combinations was adapted from a post by scvalex at
http://compprog.wordpress.com/2007/10/17/generating-combinations-1/.
It's very elegant and saved me a headache---many thanks!*/
int next_comb(int *comb, int k, int n)
{
	/*'comb' is vector of length 'k', initialised
		outside of this function
	'k' is number of combinations to choose
	'n' is number of elements to choose 'k' items from*/
	
	int i=k-1;
	comb[i]++;
	while((i>=0)&&(comb[i]>=n-k+1+i))
	{
		i--;
		if(i>=0) comb[i]++;
	}
	//end if no more combinations can be generated
	if(comb[0]>n-k) return 0;
	//else adjust comb vector
	for(i=i+1;i<k;i++) comb[i]=comb[i-1]+1;
	return 1;
}

/*bubble sort algorithm (sorts into decreasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_dec(double *values, int *ind, int n)
{
	/*'values' corresponds to unsorted vector
	'ind' corresponds to order indicators (used to sort other
		vectors in the same order as 'values')
	'n' is vector length*/
	int i,j,tempind;
	double temp;
	for(i=0;i<n;i++) ind[i]=i;
	i=n;
	while(i>0)
	{
		for(j=0;j<(i-1);j++)
		{
			temp=values[j];
			tempind=ind[j];
			if(temp<values[j+1])
			{
				values[j]=values[j+1];
				values[j+1]=temp;
				ind[j]=ind[j+1];
				ind[j+1]=tempind;
			}
		}
		i--;
	}
	return;
}

/*bubble sort algorithm for INTEGERS (sorts into decreasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_dec_int(int *values, int *ind, int n)
{
	/*'values' corresponds to unsorted vector
	'ind' corresponds to order indicators (used to sort other
		vectors in the same order as 'values')
	'n' is vector length*/
	int i,j,tempind,temp;
	for(i=0;i<n;i++) ind[i]=i;
	i=n;
	while(i>0)
	{
		for(j=0;j<(i-1);j++)
		{
			temp=values[j];
			tempind=ind[j];
			if(temp<values[j+1])
			{
				values[j]=values[j+1];
				values[j+1]=temp;
				ind[j]=ind[j+1];
				ind[j+1]=tempind;
			}
		}
		i--;
	}
	return;
}

/*bubble sort algorithm for INTEGERS (sorts into increasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_inc_int(int *values, int *ind, int n)
{
	/*'values' corresponds to unsorted vector
	'ind' corresponds to order indicators (used to sort other
		vectors in the same order as 'values')
	'n' is vector length*/
	int i,j,tempind,temp;
	for(i=0;i<n;i++) ind[i]=i;
	i=n;
	while(i>0)
	{
		for(j=0;j<(i-1);j++)
		{
			temp=values[j];
			tempind=ind[j];
			if(temp>values[j+1])
			{
				values[j]=values[j+1];
				values[j+1]=temp;
				ind[j]=ind[j+1];
				ind[j+1]=tempind;
			}
		}
		i--;
	}
	return;
}


//internal function to check for duplicated combinations
void check_duplicates(int start, int ncomb, int ncol, int *comb, int *ind)
{
	/*'start' corresponds to which element of 'comb' is the current
		comparison combination
	'ncomb' is the length of 'comb' (e.g. the number of rows)
	'ncol' is the number of columns of 'comb'
	'comb' is a matrix where each row corresponds to a unique
		combination of length 'ncol'
	'ind' is binary vector of length 'ncomb' recording unique combinations*/
		
	int i,j,k;
	
	for(i=(start+1);i<ncomb;i++)
	{
		ind[i]=1;
		for(j=0;j<ncol;j++) for(k=0;k<ncol;k++) if(comb[index2(i,j,ncomb)]==comb[index2(start,k,ncomb)]) ind[i]=0;
	}
	return;
}

//function to perform memory reallocation of intermediate structures if required
void realloc_maxmod(int nsamp, int *maxmodels, int incmodels, int **structure)
{
	/*'nsamp' is number of samples
	'maxmodels' records maximum number of models
	'incmodels' is value with which to increment 'maxmodels'
	'structure' is matrix requiring reallocation*/
	
	int i,j;
	int * temp_ptr = NULL;
	*maxmodels = *maxmodels+incmodels;
	temp_ptr = (int *) Calloc(*maxmodels*nsamp,int);
	for(i=0;i<(*maxmodels-incmodels);i++) for(j=0;j<nsamp;j++) temp_ptr[index2(i,j,*maxmodels)]=(*structure)[index2(i,j,*maxmodels-incmodels)];
	Free(*structure);
	*structure = temp_ptr;
	temp_ptr = NULL;
	return;
}

//function to perform memory reallocation of full model set if required
void realloc_maxmod_full(int nsamp, int *maxmodels, int incmodels, int **models_num)
{
	/*'nsamp' is number of samples
	'maxmodels' records maximum number of models
	'incmodels' is value with which to increment 'maxmodels'
	'models_num' is matrix requiring reallocation*/
	
	int i,j;
	int * temp_ptr = NULL;
	*maxmodels = *maxmodels+incmodels;
	temp_ptr = (int *) Calloc(*maxmodels*2*nsamp,int);
	for(i=0;i<(*maxmodels-incmodels);i++) for(j=0;j<(2*nsamp);j++) temp_ptr[index2(i,j,*maxmodels)]=(*models_num)[index2(i,j,*maxmodels-incmodels)];
	Free(*models_num);
	*models_num = temp_ptr;
	temp_ptr = NULL;
	return;
}

//recursive function to be called from 'genfullmodels'
void genfullmodels_recur(int nhier, int currhier1, int *temp, int *tempind, int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int **models_num)
{
	/*'nhier' is number of hierarchies
	'currhier1' is current hierarchy
	'temp' and 'tempind' are vectors recording structural information
	'q' is specific row of 'structure' we wish to extract
	'nrow' and 'ncol' are using for indexing
	'structure' is a matrix containing intermediate model structures
	'r' denotes current model in 'models_num'
	'totmods' is maximum number of models in 'models_num'
	'incmods' is arbitrary number of models to increase if reallocation necessary
	'models_num' is a matrix recording final models*/
	
	//declare necessary variables
	int i,j,k,currhier,nsamp,nmod=10;
	nsamp=ncol;
	currhier=currhier1+1;
	//continue working through hierarchies in order and expand model set
	if(temp[currhier]==1)
	{
		for(j=0;j<nmod;j++)
		{
			(*models_num)[index2(*r,tempind[currhier],*totmods)]=j;
			if(currhier<(nhier-1)) genfullmodels_recur(nhier,currhier,temp,tempind,q,nrow,ncol,structure,r,totmods,incmods,models_num);
			else
			{
				*r=*r+1;
				//perform memory reallocation if required
				if(*r>=*totmods) realloc_maxmod_full(nsamp,totmods,incmods,models_num);
				//initialise new model structure
				for(k=0;k<(2*nsamp);k++) (*models_num)[index2(*r,k,*totmods)]=(*models_num)[index2(*r-1,k,*totmods)];
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
					(*models_num)[index2(*r,k,*totmods)]=j;
					(*models_num)[index2(*r,k+nsamp,*totmods)]=tempind[currhier];
				}
			}
			if(currhier<(nhier-1)) genfullmodels_recur(nhier,currhier,temp,tempind,q,nrow,ncol,structure,r,totmods,incmods,models_num);
			else
			{
				*r=*r+1;
				//perform memory reallocation if required
				if(*r>=*totmods) realloc_maxmod_full(nsamp,totmods,incmods,models_num);
				//initialise new model structure
				for(k=0;k<(2*nsamp);k++) (*models_num)[index2(*r,k,*totmods)]=(*models_num)[index2(*r-1,k,*totmods)];
			}
		}
	}
	return;
}

//function to generate all possible models based on an intermediate set of structures
void genfullmodels(int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int **models_num)
{
	/*'q' is specific row of 'structure' we wish to extract
	'nrow' and 'ncol' are using for indexing
	'structure' is a matrix containing intermediate model structures
	'r' denotes current model in 'models_num'
	'totmods' is maximum number of models in 'models_num'
	'incmods' is arbitrary number of models to increase if reallocation necessary
	'models_num' is a matrix recording final models*/
	
	//declare necessary variables
	int i,j,k,l,nhier,nsamp,nmod=10;
	nsamp=ncol;
	int * temp = (int *) Calloc(nsamp,int);
	int * tempind = (int *) Calloc(nsamp,int);
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
	//now work through hierarchies in order and expand model set
	if(nhier==1)
	{
		for(i=0;i<nmod;i++)
		{
			for(j=0;j<nsamp;j++)
			{
				(*models_num)[index2(*r,j,*totmods)]=i;
				(*models_num)[index2(*r,j+nsamp,*totmods)]=0;
			}
			*r=*r+1;
			//perform memory reallocation if required
			if(*r>=*totmods) realloc_maxmod_full(nsamp,totmods,incmods,models_num);
		}
	}
	else
	{
		//temp[0] must be >1 in this loop
		for(j=0;j<nmod;j++)
		{
			//initialise structures
			for(k=0;k<nsamp;k++) (*models_num)[index2(*r,k+nsamp,*totmods)]=k;
			//now fill in gaps and pass to recursive function
			for(k=0;k<nsamp;k++)
			{
				if(structure[index2(q,k,nrow)]==tempind[0])
				{
					(*models_num)[index2(*r,k,*totmods)]=j;
					(*models_num)[index2(*r,k+nsamp,*totmods)]=tempind[0];
				}
			}
			genfullmodels_recur(nhier,0,temp,tempind,q,nrow,ncol,structure,r,totmods,incmods,models_num);
		}
	}
	//free memory from the heap (automatically sets pointer to NULL)
	Free(temp);Free(tempind);
	return;
}

//recursive function to be called from 'genmodels'
void genmodels_recur(int nsamp, int nstart, int *indexl, int *indexl_ind, int repeatind, int *comb_next, int ncomb_next, int ntotsamp, int *currstruct, int *q, int **structure, int *maxmodels, int incmodels)
{
	/*'nsamp' is number of remaining samples to choose from
	'nstart' is number of samples to start selecting
	'indexl' is vector recording which samples are remaining
	'indexl_ind' is a binary indictator denoting which samples
		can be selected or not
	'repeatind' is binary indicator 
		= 1 if repeated denominator, 0 otherwise
	'comb_next' is matrix recording available samples for special
		case of repeated denominators
	'ncomb_next' is number of rows of 'comb_next' (used for indexing)
	'ntotsamp' is total number of samples
	'currstruct' is vector recording initial structure at previous hierarchy
	'q' is current model structure indicator
	'structure' is matrix recording models structures
	'maxmodels' and 'incmodels' record the maximum number of models
		and the number to increment if reallocation necessary*/
	
	int i,j,k,l,m,nrem,ncomb,ncomb1,nstart1,valid,stop;
	//declare intermediate vectors for recursion
	int *indexl1 = (int *) Calloc(ntotsamp,int);
	int *indexl_ind1 = (int *) Calloc(ntotsamp,int);
	int *currstruct1 = (int *) Calloc(ntotsamp,int);
	int *comb,*comb_set,*comb_ind,*comb_next1;
	
	if(repeatind==1) stop=nstart-1;
	else stop=0;
	for(k=nstart;k>stop;k--)
	{
		if(k==1)
		{
			//increment models
			*q=*q+1;
			//initialise new model structure
			for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
			//perform memory reallocation if required
			if(*q>=*maxmodels) realloc_maxmod(ntotsamp,maxmodels,incmodels,structure);
		}
		else
		{
			/*calculate how many samples remaining, and if necessary
			set up for further loops*/
			nrem=nsamp-k;
			nstart1=(nrem>k ? k:nrem);
			if(nstart1>0)
			{
				if(repeatind==0)
				{
					//create set of combinations to use during looping
					ncomb = choose_c(nsamp,k);
					comb = (int *) Calloc(k,int);
					comb_set = (int *) Calloc(k*ncomb,int);
					comb_ind = (int *) Calloc(ncomb,int);
					//produce combinations
					for(j=0;j<k;j++) comb[j]=j;
					for(j=0;j<k;j++) comb_set[index2(0,j,ncomb)]=indexl[comb[j]];
					i=1;
					while(next_comb(comb,k,nsamp)==1)
					{
						for(j=0;j<k;j++) comb_set[index2(i,j,ncomb)]=indexl[comb[j]];
						i++;
					}
				}
				else
				{
					//create set of combinations to use during looping
					ncomb = ncomb_next;
					comb = (int *) Calloc(k,int);
					comb_set = (int *) Calloc(k*ncomb,int);
					comb_ind = (int *) Calloc(ncomb,int);
					//produce combinations
					for(j=0;j<ncomb;j++) for(m=0;m<k;m++) comb_set[index2(j,m,ncomb)]=comb_next[index2(j,m,ncomb)];
				}
				if(nstart1==k)
				{
					//now begin looping over combinations
					i=0;valid=1;
					while(i<(ncomb-1)&&valid==1)
					{
						for(j=0;j<ncomb;j++) comb_ind[j]=1;
						//extract subset that can be looped over further
						check_duplicates(i,ncomb,k,comb_set,comb_ind);
						ncomb1=0;
						for(j=(i+1);j<ncomb;j++) ncomb1+=comb_ind[j];
						if(ncomb1==0) valid=0;
						else
						{
							//reset indicators
							for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
							for(j=0;j<ntotsamp;j++)
							{
								indexl1[j]=indexl[j];
								indexl_ind1[j]=indexl_ind[j];
							}
							for(j=0;j<k;j++) (*structure)[index2(*q,comb_set[index2(i,j,ncomb)],*maxmodels)]=comb_set[index2(i,0,ncomb)];
							if(nrem>0)
							{
								if(nrem==1)
								{
									//increment models
									*q=*q+1;
									//initialise new model structure
									for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
									//perform memory reallocation if required
									if(*q>=*maxmodels) realloc_maxmod(ntotsamp,maxmodels,incmodels,structure);
								}
								else
								{
									//extract viable choices for next level
									comb_next1 = (int *) Calloc(ncomb1*k,int);
									l=0;
									for(j=(i+1);j<ncomb;j++)
									{
										if(comb_ind[j]==1)
										{
											for(m=0;m<k;m++) comb_next1[index2(l,m,ncomb1)]=comb_set[index2(j,m,ncomb)];
											l++;
										}
									}
									//extract remaining samples for loop
									for(j=0;j<k;j++) indexl_ind1[comb_set[index2(i,j,ncomb)]]=0;
									l=0;
									for(j=0;j<ntotsamp;j++)
									{
										if(indexl_ind1[j]==1)
										{
											indexl1[l]=j;
											l++;
										}
									}
									for(j=0;j<ntotsamp;j++) currstruct1[j]=(*structure)[index2(*q,j,*maxmodels)];
									genmodels_recur(nrem,nstart1,indexl1,indexl_ind1,1,comb_next1,ncomb1,ntotsamp,currstruct1,q,structure,maxmodels,incmodels);
									//free memory from heap
									Free(comb_next1);comb_next1=NULL;
								}
							}
							else
							{
								//increment models
								*q=*q+1;
								//initialise new model structure
								for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
								//perform memory reallocation if required
								if(*q>=*maxmodels) realloc_maxmod(ntotsamp,maxmodels,incmodels,structure);
							}
						}
						i++;
					}
					//if necessary, then loop over further hierarchies
					nstart1--;
					if(nstart1>1)
					{
						while(nstart1>1)
						{
							//loop over combinations
							for(i=0;i<ncomb;i++)
							{
								//reset indicators
								for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
								for(j=0;j<ntotsamp;j++)
								{
									indexl1[j]=indexl[j];
									indexl_ind1[j]=indexl_ind[j];
								}
								for(j=0;j<k;j++) (*structure)[index2(*q,comb_set[index2(i,j,ncomb)],*maxmodels)]=comb_set[index2(i,0,ncomb)];
								if(nrem>0)
								{
									//extract remaining samples for loop
									for(j=0;j<k;j++) indexl_ind1[comb_set[index2(i,j,ncomb)]]=0;
									l=0;
									for(j=0;j<ntotsamp;j++)
									{
										if(indexl_ind1[j]==1)
										{
											indexl1[l]=j;
											l++;
										}
									}
									if(nrem==1)
									{
										//increment models
										*q=*q+1;
										//initialise new model structure
										for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
										//perform memory reallocation if required
										if(*q>=*maxmodels) realloc_maxmod(ntotsamp,maxmodels,incmodels,structure);
									}
									else
									{
										for(j=0;j<ntotsamp;j++) currstruct1[j]=(*structure)[index2(*q,j,*maxmodels)];
										genmodels_recur(nrem,nstart1,indexl1,indexl_ind1,0,comb_next1,0,ntotsamp,currstruct1,q,structure,maxmodels,incmodels);
									}
								}
								else
								{
									//increment models
									*q=*q+1;
									//initialise new model structure
									for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
									//perform memory reallocation if required
									if(*q>=*maxmodels) realloc_maxmod(ntotsamp,maxmodels,incmodels,structure);
								}
							}
							nstart1--;
						}
					}
					else
					{
						//now repeat for non-repeated hierarchies - looping over combinations
						for(i=0;i<ncomb;i++)
						{
							for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
							for(j=0;j<k;j++) (*structure)[index2(*q,comb_set[index2(i,j,ncomb)],*maxmodels)]=comb_set[index2(i,0,ncomb)];
							//increment models
							*q=*q+1;
							//initialise new model structure
							for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
							//perform memory reallocation if required
							if(*q>=*maxmodels) realloc_maxmod(ntotsamp,maxmodels,incmodels,structure);
						}
					}
				}
				else
				{
					//loop over combinations
					for(i=0;i<ncomb;i++)
					{
						//reset indicators
						for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
						for(j=0;j<ntotsamp;j++)
						{
							indexl1[j]=indexl[j];
							indexl_ind1[j]=indexl_ind[j];
						}
						for(j=0;j<k;j++) (*structure)[index2(*q,comb_set[index2(i,j,ncomb)],*maxmodels)]=comb_set[index2(i,0,ncomb)];
						if(nrem>0)
						{
							//extract remaining samples for loop
							for(j=0;j<k;j++) indexl_ind1[comb_set[index2(i,j,ncomb)]]=0;
							l=0;
							for(j=0;j<ntotsamp;j++)
							{
								if(indexl_ind1[j]==1)
								{
									indexl1[l]=j;
									l++;
								}
							}
							if(nrem==1)
							{
								//increment models
								*q=*q+1;
								//initialise new model structure
								for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
								//perform memory reallocation if required
								if(*q>=*maxmodels) realloc_maxmod(ntotsamp,maxmodels,incmodels,structure);
							}
							else
							{
								for(j=0;j<ntotsamp;j++) currstruct1[j]=(*structure)[index2(*q,j,*maxmodels)];
								genmodels_recur(nrem,nstart1,indexl1,indexl_ind1,0,comb_next1,0,ntotsamp,currstruct1,q,structure,maxmodels,incmodels);
							}
						}
						else
						{
							//increment models
							*q=*q+1;
							//initialise new model structure
							for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
							//perform memory reallocation if required
							if(*q>=*maxmodels) realloc_maxmod(ntotsamp,maxmodels,incmodels,structure);
						}
					}
				}
				//free memory from the heap (automatically sets pointer to NULL)
				Free(comb);Free(comb_set);Free(comb_ind);
			}
			else
			{
				for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
				for(j=0;j<k;j++) (*structure)[index2(*q,indexl[j],*maxmodels)]=indexl[0];
				//increment models
				*q=*q+1;
				//initialise new model structure
				for(j=0;j<ntotsamp;j++) (*structure)[index2(*q,j,*maxmodels)]=currstruct[j];
				//perform memory reallocation if required
				if(*q>=*maxmodels) realloc_maxmod(ntotsamp,maxmodels,incmodels,structure);
			}
		}
	}
	//free memory from heap
	Free(indexl1);Free(indexl_ind1);Free(currstruct1);
	return;
}

/*function to generate all model structures in an efficient manner
(returns an R vector of a given size, though it is not necessary to
specify the size in advance)*/
SEXP genmodels(SEXP nsamp1)
{
	/*'nsamp1' is the number of samples*/
	
	//declare variables
	int i,j,k,l,m,q,r,nsamp,nrem,nstart,ncomb,ncomb1,valid;
	//extract integer parts of 'nsamp1'
	nsamp1 = coerceVector(nsamp1, INTSXP);
	nsamp=INTEGER(nsamp1)[0];
	//pointer for memory reallocation
	int * temp_ptr = NULL;
	
	/*this is arbitrary number of initial models
	(can be altered if necessary during runtime)*/
	int maxmodels,incmodels=10000;
	maxmodels=incmodels;
	/*set up pointer to array of pointers, each of which will point at an
	array holding relevant model structures (to be filled in during runtime)*/
	int *structure = (int *) Calloc(maxmodels*nsamp,int);
	//declare intermediate vectors for recursion
	int *indexl = (int *) Calloc(nsamp,int);
	int *indexl_ind = (int *) Calloc(nsamp,int);
	int *currstruct = (int *) Calloc(nsamp,int);
	int *comb,*comb_set,*comb_ind,*comb_next;
	
	//for fully DEPENDENT model structures
	q=0;
	for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=0;
	q++;
	
	//now complete for PARTIALLY DEPENDENT model structures
	for(k=(nsamp-1);k>1;k--)
	{
		//create set of combinations to use during looping
		ncomb = choose_c(nsamp,k);
		comb = (int *) Calloc(k,int);
		comb_set = (int *) Calloc(k*ncomb,int);
		comb_ind = (int *) Calloc(ncomb,int);
		//produce combinations
		for(j=0;j<k;j++) comb[j]=j;
		for(j=0;j<k;j++) comb_set[index2(0,j,ncomb)]=comb[j];
		i=1;
		while(next_comb(comb,k,nsamp)==1)
		{
			for(j=0;j<k;j++) comb_set[index2(i,j,ncomb)]=comb[j];
			i++;
		}
		/*calculate how many samples remaining, and if necessary
		set up for further loops*/
		nrem=nsamp-k;
		nstart=(nrem>k ? k:nrem);
		if(nstart>1)
		{
			if(nstart==k)
			{
				//now begin looping over combinations
				i=0;valid=1;
				while(i<(ncomb-1)&&valid==1)
				{
					//extract subset that can be looped over further
					check_duplicates(i,ncomb,k,comb_set,comb_ind);
					ncomb1=0;
					for(j=(i+1);j<ncomb;j++) ncomb1+=comb_ind[j];
					if(ncomb1==0) valid=0;
					else
					{
						//reset indicators
						for(j=0;j<nsamp;j++) indexl_ind[j]=1;
						//initialise structures
						for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
						for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
						//extract viable choices for next level
						comb_next = (int *) Calloc(ncomb1*k,int);
						l=0;
						for(j=(i+1);j<ncomb;j++)
						{
							if(comb_ind[j]==1)
							{
								for(m=0;m<k;m++) comb_next[index2(l,m,ncomb1)]=comb_set[index2(j,m,ncomb)];
								l++;
							}
						}
						//extract remaining samples for loop
						for(j=0;j<k;j++) indexl_ind[comb_set[index2(i,j,ncomb)]]=0;
						l=0;
						for(j=0;j<nsamp;j++)
						{
							if(indexl_ind[j]==1)
							{
								indexl[l]=j;
								l++;
							}
						}
						for(j=0;j<nsamp;j++) currstruct[j]=structure[index2(q,j,maxmodels)];
						genmodels_recur(nrem,nstart,indexl,indexl_ind,1,comb_next,ncomb1,nsamp,currstruct,&q,&structure,&maxmodels,incmodels);
						//free memory from heap
						Free(comb_next);comb_next=NULL;
					}
					i++;
				}
				//now enter further recursive loop if required
				nstart--;
				if(nstart>1)
				{
					while(nstart>1)
					{
						//loop over combinations
						for(i=0;i<ncomb;i++)
						{
							//reset indicators
							for(j=0;j<nsamp;j++) indexl_ind[j]=1;
							//initialise structures
							for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
							for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
							//extract remaining samples for loop
							for(j=0;j<k;j++) indexl_ind[comb_set[index2(i,j,ncomb)]]=0;
							l=0;
							for(j=0;j<nsamp;j++)
							{
								if(indexl_ind[j]==1)
								{
									indexl[l]=j;
									l++;
								}
							}
							for(j=0;j<nsamp;j++) currstruct[j]=structure[index2(q,j,maxmodels)];
							genmodels_recur(nrem,nstart,indexl,indexl_ind,0,comb_next,0,nsamp,currstruct,&q,&structure,&maxmodels,incmodels);
						}
						nstart--;	
					}
				}
				else
				{
					//now repeat for non-repeated hierarchies - looping over combinations
					for(l=0;l<ncomb;l++)
					{
						for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
						for(j=0;j<k;j++) structure[index2(q,comb_set[index2(l,j,ncomb)],maxmodels)]=comb_set[index2(l,0,ncomb)];
						q++;
						//perform memory reallocation if required
						if(q>=maxmodels) realloc_maxmod(nsamp,&maxmodels,incmodels,&structure);
					}
				}	
			}
			else
			{
				//loop over combinations
				for(i=0;i<ncomb;i++)
				{
					//reset indicators
					for(j=0;j<nsamp;j++) indexl_ind[j]=1;
					//initialise structures
					for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
					for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
					//extract remaining samples for loop
					for(j=0;j<k;j++) indexl_ind[comb_set[index2(i,j,ncomb)]]=0;
					l=0;
					for(j=0;j<nsamp;j++)
					{
						if(indexl_ind[j]==1)
						{
							indexl[l]=j;
							l++;
						}
					}
					for(j=0;j<nsamp;j++) currstruct[j]=structure[index2(q,j,maxmodels)];
					genmodels_recur(nrem,nstart,indexl,indexl_ind,0,comb_next,0,nsamp,currstruct,&q,&structure,&maxmodels,incmodels);
				}	
			}
		}
		else
		{
			//complete model structures
			for(i=0;i<ncomb;i++)
			{
				for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
				for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
				q++;
				//perform memory reallocation if required
				if(q>=maxmodels) realloc_maxmod(nsamp,&maxmodels,incmodels,&structure);
			}
		}
		//free memory from the heap (automatically sets pointer to NULL)
		Free(comb);Free(comb_ind);Free(comb_set);
	}
	//now reallocate memory to save space
	temp_ptr = (int *) Calloc(q*nsamp,int);
	for(i=0;i<nsamp;i++) for(j=0;j<q;j++) temp_ptr[index2(j,i,q)]=structure[index2(j,i,maxmodels)];
	Free(structure);
	structure = temp_ptr;
	temp_ptr = NULL;
	//now generate all possible model structures
	int totmods,incmods=100000;
	totmods=incmods;
	int * models_num = (int *) Calloc(totmods*2*nsamp,int);
	r=0;
	for(i=0;i<q;i++) genfullmodels(i,q,nsamp,structure,&r,&totmods,incmods,&models_num);
	
	//now append FULLY INDEPENDENT models
	int nmod=10;
	l=r+pow_int(nmod,nsamp);
	//now reallocate memory to save space
	temp_ptr = (int *) Calloc(l*2*nsamp,int);
	for(i=0;i<(2*nsamp);i++) for(j=0;j<r;j++) temp_ptr[index2(j,i,l)]=models_num[index2(j,i,totmods)];
	Free(models_num);
	models_num = temp_ptr;
	temp_ptr = NULL;
	totmods=l;
	//initialise indicator part of models_num
	for(i=r;i<totmods;i++) for(j=0;j<nsamp;j++) models_num[index2(i,j+nsamp,totmods)]=j;
	//expand grid of independent models into models_num
	int each,rep;
	for(i=0;i<nsamp;i++)
	{
		l=r;
		each=pow_int(nmod,i);
		rep=pow_int(nmod,nsamp-i-1);
		while(rep>0)
		{
			for(k=0;k<nmod;k++)
			{
				for(j=0;j<each;j++)
				{
					models_num[index2(l,i,totmods)]=k;
					l++;
				}
			}
			rep--;
		}
	}
	//now export into R and return
	SEXP models_num_R;
	PROTECT(models_num_R = allocVector(INTSXP,totmods*2*nsamp));
	for(i=0;i<totmods;i++) for(j=0;j<(2*nsamp);j++) INTEGER(models_num_R)[index2(i,j,totmods)] = models_num[index2(i,j,totmods)];
	UNPROTECT(1);
	//free memory from the heap (automatically sets pointer to NULL)
	Free(structure);Free(models_num);Free(indexl);Free(indexl_ind);Free(currstruct);
	return models_num_R;
}

//function to calculate l[P'(D|M)] for a given distribution of bases
double lPDM_mod_fn(int *z, int ind, double pstar)
{	
	/*'z' is a pointer to a vector of length 5, where the
		first 4 elements correspond to each base (with the consensus
		as the fourth element z[3]). The final element is 
		S-z[3]=sum(z[0:2])
	'ind' denotes which model (from 0:9) is to be calculated
	'pstar' is the overall mutation rate*/
	
	double lPDM;
/*	switch(ind)*/
/*	{*/
/*		//Null p1=p2=p3=p3*/
/*		case 0 :lPDM=z[4]*log(pstar/3.0)+z[3]*log(1.0-pstar);*/
/*			break;*/
/*		//Alt that one free pi is different: e.g. p1!=p2=p3 etc. but mutation rate constrained to p**/
/*		case 1 :lPDM=-(z[1]+z[2])*log(2.0)+lfactorial(z[0])+lfactorial(z[1]+z[2])-lfactorial(z[4]+1)+z[4]*log(pstar)+z[3]*log(1.0-pstar);*/
/*			break;*/
/*		case 2 :lPDM=-(z[0]+z[2])*log(2.0)+lfactorial(z[1])+lfactorial(z[0]+z[2])-lfactorial(z[4]+1)+z[4]*log(pstar)+z[3]*log(1.0-pstar);*/
/*			break;*/
/*		case 3 :lPDM=-(z[0]+z[1])*log(2.0)+lfactorial(z[2])+lfactorial(z[0]+z[1])-lfactorial(z[4]+1)+z[4]*log(pstar)+z[3]*log(1.0-pstar);*/
/*			break;*/
/*		//Alt that all pis are different but mutation rate constrained to p**/
/*		case 4 :lPDM=log(2.0)+lfactorial(z[0])+lfactorial(z[1])+lfactorial(z[2])-lfactorial(z[4]+2)+z[4]*log(pstar)+z[3]*log(1.0-pstar);*/
/*			break;*/
/*		//Alt that pis are uniform but not constrained to sum to p**/
/*		case 5 :lPDM=-z[4]*log(3.0)+lfactorial(z[4])+lfactorial(z[3])-lfactorial(z[3]+z[4]+1);*/
/*			break;*/
/*		//Alt that one free pi is different: e.g. p1!=p2=p3 etc.*/
/*		case 6 :lPDM=-(z[1]+z[2])*log(2.0)+lfactorial(z[1]+z[2])+lfactorial(z[0])+lfactorial(z[3])-log(z[4]+1)-lfactorial(z[4]+z[3]+1);*/
/*			break;*/
/*		case 7 :lPDM=-(z[0]+z[2])*log(2.0)+lfactorial(z[0]+z[2])+lfactorial(z[1])+lfactorial(z[3])-log(z[4]+1)-lfactorial(z[4]+z[3]+1);*/
/*			break;*/
/*		case 8 :lPDM=-(z[0]+z[1])*log(2.0)+lfactorial(z[0]+z[1])+lfactorial(z[2])+lfactorial(z[3])-log(z[4]+1)-lfactorial(z[4]+z[3]+1);*/
/*			break;*/
/*		//Alt that all pis are different*/
/*		case 9 :lPDM=log(2.0)+lfactorial(z[0])+lfactorial(z[1])+lfactorial(z[2])+lfactorial(z[3])-log(z[4]+2)-log(z[4]+1)-lfactorial(z[4]+z[3]+1);*/
/*			break;*/
/*	}*/
	int i;
	double z1[5];
	for(i=0;i<5;i++) z1[i]=(double) z[i];
	switch(ind)
	{
		//Null p1=p2=p3=p/3 where p<=p*
		case 0 :lPDM=-z1[4]*log(3.0)-log(pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,1,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		//Alt that one free pi is different: e.g. p1!=p2=p3 etc. but mutation rate constrained to be <= p*
		case 1 :lPDM=-(z[1]+z[2])*log(2.0)+lfactorial(z[0])+lfactorial(z[1]+z[2])-lfactorial(z[4]+1)-log(pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,1,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		case 2 :lPDM=-(z[0]+z[2])*log(2.0)+lfactorial(z[1])+lfactorial(z[0]+z[2])-lfactorial(z[4]+1)-log(pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,1,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		case 3 :lPDM=-(z[0]+z[1])*log(2.0)+lfactorial(z[2])+lfactorial(z[0]+z[1])-lfactorial(z[4]+1)-log(pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,1,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		//Alt that all pis are different but mutation rate constrained to be <= p*
		case 4 :lPDM=log(2.0)+lfactorial(z[0])+lfactorial(z[1])+lfactorial(z[2])-lfactorial(z[4]+2)-log(pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,1,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		//Alt p1=p2=p3=p/3 where p>p*
		case 5 :lPDM=-z1[4]*log(3.0)-log(1.0-pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,0,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		//Alt that one free pi is different: e.g. p1!=p2=p3 etc. but mutation rate constrained to be > p*
		case 6 :lPDM=-(z[1]+z[2])*log(2.0)+lfactorial(z[0])+lfactorial(z[1]+z[2])-lfactorial(z[4]+1)-log(1.0-pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,0,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		case 7 :lPDM=-(z[0]+z[2])*log(2.0)+lfactorial(z[1])+lfactorial(z[0]+z[2])-lfactorial(z[4]+1)-log(1.0-pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,0,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		case 8 :lPDM=-(z[0]+z[1])*log(2.0)+lfactorial(z[2])+lfactorial(z[0]+z[1])-lfactorial(z[4]+1)-log(1.0-pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,0,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		//Alt that all pis are different but mutation rate constrained to be > p*
		case 9 :lPDM=log(2.0)+lfactorial(z[0])+lfactorial(z[1])+lfactorial(z[2])-lfactorial(z[4]+2)-log(1.0-pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,0,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
	}
	return lPDM;
}

//recursive function to produce lexicographic orderings
void lexico_recur(int nsamp, int col, int totcol, int *temp_samp, int *temp_order1, int norder, int start, int *overall_order, int *ncomb_sub)
{
	/*'nsamp' is the number of samples
	'col' is current column
	'totcol' is total number of columns
	'temp_samp' is a vector recording where to start ordering in absolute terms
	'temp_order1' is a matrix (in vector form) containing current orderings
	'norder' is number of rows of 'temp_order1'
	'start' is the absolute starting point in 'temp_order1' used for indexing in recursion
	'overall_order' is vector containing absolute orderings across all possible dependencies
	'ncomb_sub' is intermediate vector used for indexing 'overall_order'*/
	
	int q,l,k,temp_norder;
	int * temp_ind = (int *) Calloc(norder,int);
	int * temp_val = (int *) Calloc(norder,int);
	int * temp_samp1 = (int *) Calloc(nsamp,int);
	
	for(q=1;q<nsamp;q++)
	{
		temp_norder=temp_samp[q]-temp_samp[q-1];
		if(temp_norder>1)
		{
			bubble_sort_inc_int(&temp_order1[index2(start+temp_samp[q-1],col,norder)],temp_ind,temp_norder);
/*			Rprintf("totcol=%d\tcol=%d\tq=%d\n",totcol,col,q);*/
/*			for(l=0;l<temp_norder;l++) Rprintf("%d\n",temp_ind[l]);*/
/*			Rprintf("\n");*/
			for(l=(col+1);l<totcol;l++)
			{
				for(k=0;k<temp_norder;k++) temp_val[k]=temp_order1[index2(start+temp_samp[q-1]+k,l,norder)];
				for(k=0;k<temp_norder;k++) temp_order1[index2(start+temp_samp[q-1]+k,l,norder)]=temp_val[temp_ind[k]];
			}
			//adjust overall ordering
			for(k=0;k<temp_norder;k++) temp_val[k]=overall_order[ncomb_sub[totcol-1]+start+temp_samp[q-1]+k];
			for(k=0;k<temp_norder;k++) overall_order[ncomb_sub[totcol-1]+start+temp_samp[q-1]+k]=temp_val[temp_ind[k]];
			//now generate conditions for sequential sorting at next level
			for(k=0;k<nsamp;k++)
			{
				temp_samp1[k]=0;
				for(l=(start+temp_samp[q-1]);l<(start+temp_samp[q-1]+temp_norder);l++) temp_samp1[k]+=(temp_order1[index2(l,col,norder)]==k ? 1:0);
			}
			//convert to cumulative counts
			for(l=1;l<nsamp;l++) temp_samp1[l]+=temp_samp1[l-1];
			//enter into further recursion
			if((col+1)<totcol) lexico_recur(nsamp,col+1,totcol,temp_samp1,temp_order1,norder,start+temp_samp[q-1],overall_order,ncomb_sub);
		}
	}
	Free(temp_ind);Free(temp_val);Free(temp_samp1);
	return;
}

//function to calculate intermediate values for log[P'(D|M)]
void calc_lPDM_int_fn(int *nsamp1, int *z, double *pstar1, double *lPDM_int_mat)
{	
	/*'nsamp1' is the number of samples
	'z' is a vector of length 4 giving the distribution of bases
		(such that the consensus base is at the last position)
	'pstar1' is the overall mutation rate
	'lPDM_int_mat' is a matrix (in vector form) recording the values of log[P'(D|M)]
		for different models and sample combinations*/
	int i,j,k,l,m,q,nsamp,nmod=10;
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
	int norder,temp_norder;
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

//function to simplify a single unique model and return the log[P'(D|M)]
double calc_lPDMs_single(int r, int nsamp, int totmods, int nmod, int *models_num, double *lPDM_int_mat, int *ncomb_sub, int *sum, int *temp_index, int *mods, int *inds)
{
	/*'r' denotes the current model (the r^th row of 'models_num')
	'nsamp' is the number of samples
	'totmods' is the number of rows of 'models_num'
	'nmod' is the number of models
	'models_num' is a matrix (in vector form) recording each unique model structure
	'lPDM_int_mat' is a matrix (in vector form) recording the values of log[P'(D|M)]
		for different models and sample combinations
	'ncomb_sub' is a vector recording numbers of combinations (to be used for selecting
		correct columns of 'lPDM_int_mat'
	'sum', 'temp_index', 'mods' and 'inds' are all intermediate vectors for use during calculation*/
	
	int j,k,l,m;
	double lPDM;
	
	//first part of code simplifies the model structures and indicators
	for(j=0;j<nsamp;j++)
	{
		//count number of times indicator j appears in model
		sum[j]=0;
		for(k=0;k<nsamp;k++) if(models_num[index2(r,k+nsamp,totmods)]==j) sum[j]++;
	}
	//now generate correct column indicators
	k=0;
	for(j=0;j<nsamp;j++)
	{
		if(sum[j]>0)
		{
			l=0;
			m=0;
			while(l<nsamp&&m<sum[j])
			{
				if(models_num[index2(r,l+nsamp,totmods)]==j)
				{
					if(m==0) mods[k]=models_num[index2(r,l,totmods)];
					temp_index[m]=l;
					m++;
				}
				l++;
			}
			//calculate correct column
			m=0;
			for(l=0;l<sum[j];l++) m+=choose_c(temp_index[l],l+1);
			m+=ncomb_sub[sum[j]-1];
			inds[k]=m;
			k++;
		}
	}
	//now generate lPDM
	lPDM=0.0;
	for(j=0;j<k;j++) lPDM+=lPDM_int_mat[index2(mods[j],inds[j],nmod)];
	return lPDM;
}

//function to simplify a single unique model and return the log[P'(D|M)]
//(EXACTLY THE SAME AS THE ABOVE FUNCTION ONLY FOR WHEN MODELS_NUM IS A SINGLE VECTOR)
double calc_lPDMs_single_vec(int nsamp, int nmod, int *models_num, double *lPDM_int_mat, int *ncomb_sub, int *sum, int *temp_index, int *mods, int *inds)
{
	/*'r' denotes the current model (the r^th row of 'models_num')
	'nsamp' is the number of samples
	'totmods' is the number of rows of 'models_num'
	'nmod' is the number of models
	'models_num' is a vector recording the unique model structure
	'lPDM_int_mat' is a matrix (in vector form) recording the values of log[P'(D|M)]
		for different models and sample combinations
	'ncomb_sub' is a vector recording numbers of combinations (to be used for selecting
		correct columns of 'lPDM_int_mat'
	'sum', 'temp_index', 'mods' and 'inds' are all intermediate vectors for use during calculation*/
	
	int j,k,l,m;
	double lPDM;
	
	//first part of code simplifies the model structures and indicators
	for(j=0;j<nsamp;j++)
	{
		//count number of times indicator j appears in model
		sum[j]=0;
		for(k=0;k<nsamp;k++) if(models_num[k+nsamp]==j) sum[j]++;
	}
	//now generate correct column indicators
	k=0;
	for(j=0;j<nsamp;j++)
	{
		if(sum[j]>0)
		{
			l=0;
			m=0;
			while(l<nsamp&&m<sum[j])
			{
				if(models_num[l+nsamp]==j)
				{
					if(m==0) mods[k]=models_num[l];
					temp_index[m]=l;
					m++;
				}
				l++;
			}
			//calculate correct column
			m=0;
			for(l=0;l<sum[j];l++) m+=choose_c(temp_index[l],l+1);
			m+=ncomb_sub[sum[j]-1];
			inds[k]=m;
			k++;
		}
	}
	//now generate lPDM
	lPDM=0.0;
	for(j=0;j<k;j++) lPDM+=lPDM_int_mat[index2(mods[j],inds[j],nmod)];
	return lPDM;
}		

//function for calculating log[P'(D|M)]s for all unique models for a given set of base distributions
void calc_lPDM_fn(int *nsamp1, int *totmods1, int *models_num, double *lPDM_int_mat, double *lPDM_mat)
{
	int i,j,nsamp,totmods,nmod=10;
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

/*function to discriminate null and alternative models
based on less-stringent hypothesis*/
void calc_lessstring_fn(int *nsamp1, int *totmods1, int *models_num, int *hyp)
{
	/*'nsamp1' is the number of samples
	'totmods1' is the total number of rows in 'models_num'
	'models_num' is a matrix (in vector form) containing all unique models
	'hyp' is binary vector of length 'totmods1' recording whether
		a model is a null (0) or alternative (1) model*/
	
	int i,j,nsamp,totmods,inimod,iniind;
	nsamp=*nsamp1;
	totmods=*totmods1;
	
	//calculate PPAs
	for(i=0;i<totmods;i++)
	{
		inimod=models_num[index2(i,0,totmods)];
		iniind=models_num[index2(i,nsamp,totmods)];
		hyp[i]=0;
		j=1;
		while(j<nsamp&&hyp[i]==0)
		{
			if(models_num[index2(i,j,totmods)]>=5) if(models_num[index2(i,j,totmods)]!=inimod||models_num[index2(i,j+nsamp,totmods)]!=iniind) hyp[i]=1;
			j++;
		}
	}
	return;
}

/*function to discriminate null and alternative models
based on stringent hypothesis*/
void calc_string_fn(int *nsamp1, int *totmods1, int *models_num, int *hyp)
{
	/*'nsamp1' is the number of samples
	'totmods1' is the total number of rows in 'models_num'
	'models_num' is a matrix (in vector form) containing all unique models
	'hyp' is binary vector of length 'totmods1' recording whether
		a model is a null (0) or alternative (1) model*/
	
	int i,j,nsamp,totmods,inimod,iniind;
	nsamp=*nsamp1;
	totmods=*totmods1;
	
	//calculate PPAs
	for(i=0;i<totmods;i++)
	{
		inimod=models_num[index2(i,0,totmods)];
		iniind=models_num[index2(i,nsamp,totmods)];
		hyp[i]=1;
		j=1;
		while(j<nsamp&&hyp[i]==1)
		{
			if(models_num[index2(i,j,totmods)]<5) hyp[i]=0;
			else{if(models_num[index2(i,j,totmods)]==inimod&&models_num[index2(i,j+nsamp,totmods)]==iniind) hyp[i]=0;}
			j++;
		}
	}
	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*THE FOLLOWING ARE FOR THE EFFICIENT SEARCH/APPROXIMATION ROUTINES*/

/*functions to generate priors (log[P(M)]) for null and alternative models
without having to record all model structures (for use in approximation routine)*/

/*function to discriminate null and alternative models based on 
less-stringent hypothesis for a single model structure*/
int calc_lessstring_single(int nsamp, int *models_num)
{
	/*'nsamp' is the number of samples
	'models_num' is a vector containing model structure*/
	
	int j,inimod,iniind,hyp;
	
	//calculate PPAs
	inimod=models_num[0];
	iniind=models_num[nsamp];
	hyp=0;
	j=1;
	while(j<nsamp&&hyp==0)
	{
		if(models_num[j]>=5) if(models_num[j]!=inimod||models_num[j+nsamp]!=iniind) hyp=1;
		j++;
	}
	//return binary indicator (1 if alternative, 0 if null)
	return hyp;
}

/*function to discriminate null and alternative models based on 
stringent hypothesis for a single model structure*/
int calc_string_single(int nsamp, int *models_num)
{
	/*'nsamp' is the number of samples
	'models_num' is a vector containing model structure*/
	
	int j,inimod,iniind,hyp;
	
	//calculate PPAs
	inimod=models_num[0];
	iniind=models_num[nsamp];
	hyp=1;
	j=1;
	while(j<nsamp&&hyp==1)
	{
		if(models_num[j]<5) hyp=0;
		else{if(models_num[j]==inimod&&models_num[j+nsamp]==iniind) hyp=0;}
		j++;
	}
	//return binary indicator (1 if alternative, 0 if null)
	return hyp;
}

//recursive function to be called from 'genfullmodels_priors'
void genfullmodels_recur_priors(int nhier, int currhier1, int *temp, int *tempind, int q, int nrow, int ncol, int *structure, int *r, int *models_num, int *nalt_ls, int *nalt_s)
{
	/*'nhier' is number of hierarchies
	'currhier1' is current hierarchy
	'temp' and 'tempind' are vectors recording structural information
	'q' is specific row of 'structure' we wish to extract
	'nrow' and 'ncol' are using for indexing
	'structure' is a matrix containing intermediate model structures
	'r' denotes current model
	'models_num' is a vector recording model structure
	'nalt_ls' records the number of alternative models 
		according to "less stringent" criteria
	'nalt_s' records the number of alternative models 
		according to "stringent" criteria*/
		
	//declare necessary variables
	int i,j,k,currhier,nsamp,nmod=10;
	nsamp=ncol;
	currhier=currhier1+1;
	//continue working through hierarchies in order and expand model set
	if(temp[currhier]==1)
	{
		for(j=0;j<nmod;j++)
		{
			models_num[tempind[currhier]]=j;
			if(currhier<(nhier-1)) genfullmodels_recur_priors(nhier,currhier,temp,tempind,q,nrow,ncol,structure,r,models_num,nalt_ls,nalt_s);
			else
			{
				*r=*r+1;
				//now calculate prior specification based on model
				*nalt_ls=*nalt_ls+calc_lessstring_single(nsamp,models_num);
				*nalt_s=*nalt_s+calc_string_single(nsamp,models_num);
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
			if(currhier<(nhier-1)) genfullmodels_recur_priors(nhier,currhier,temp,tempind,q,nrow,ncol,structure,r,models_num,nalt_ls,nalt_s);
			else
			{
				*r=*r+1;
				//now calculate prior specification based on model
				*nalt_ls=*nalt_ls+calc_lessstring_single(nsamp,models_num);
				*nalt_s=*nalt_s+calc_string_single(nsamp,models_num);
			}
		}
	}
	return;
}

/*function to generate all possible models based on an intermediate set of structures
(returns just prior information)*/
void genfullmodels_priors(int q, int nrow, int ncol, int *structure, int *r, int *models_num, int *nalt_ls, int *nalt_s)
{
	/*'q' is specific row of 'structure' we wish to extract
	'nrow' and 'ncol' are using for indexing
	'structure' is a matrix containing intermediate model structures
	'r' denotes current model
	'models_num' is a vector recording model structures
	'nalt_ls' records the number of alternative models 
		according to "less stringent" criteria
	'nalt_s' records the number of alternative models 
		according to "stringent" criteria*/
	
	//declare necessary variables
	int i,j,k,l,nhier,nsamp,nmod=10;
	nsamp=ncol;
	int * temp = (int *) Calloc(nsamp,int);
	int * tempind = (int *) Calloc(nsamp,int);
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
	//now work through hierarchies in order and expand model set
	if(nhier==1)
	{
		for(i=0;i<nmod;i++)
		{
			for(j=0;j<nsamp;j++)
			{
				models_num[j]=i;
				models_num[j+nsamp]=0;
			}
			//now calculate prior specification based on model
			*nalt_ls=*nalt_ls+calc_lessstring_single(nsamp,models_num);
			*nalt_s=*nalt_s+calc_string_single(nsamp,models_num);
			*r=*r+1;
		}
	}
	else
	{
		//temp[0] must be >1 in this loop
		for(j=0;j<nmod;j++)
		{
			//initialise structures
			for(k=0;k<nsamp;k++) models_num[k+nsamp]=k;
			//now fill in gaps and pass to recursive function
			for(k=0;k<nsamp;k++)
			{
				if(structure[index2(q,k,nrow)]==tempind[0])
				{
					models_num[k]=j;
					models_num[k+nsamp]=tempind[0];
				}
			}
			genfullmodels_recur_priors(nhier,0,temp,tempind,q,nrow,ncol,structure,r,models_num,nalt_ls,nalt_s);
		}
	}
	//free memory from the heap (automatically sets pointer to NULL)
	Free(temp);Free(tempind);
	return;
}

/*function to generate independent model structure recursively
(for use in approximation routine)*/
void genfullmodels_indep_priors(int currcol, int nsamp, int nmod, int *models_num, int *nalt_ls, int *nalt_s)
{
	/*'currcol' is the current sample in 'models_num'
	'nsamp' is the total number of samples of 'models_num'
	'nmod' is the number of model choices
	'models_num' is a vector recording model structures
	'nalt_ls' records the number of alternative models 
		according to "less stringent" criteria
	'nalt_s' records the number of alternative models 
		according to "stringent" criteria*/
		
	int i;
	
	for(i=0;i<nmod;i++)
	{
		models_num[currcol]=i;
		if(currcol==(nsamp-1))
		{
			//now calculate prior specification based on model
			*nalt_ls=*nalt_ls+calc_lessstring_single(nsamp,models_num);
			*nalt_s=*nalt_s+calc_string_single(nsamp,models_num);
		}
		else genfullmodels_indep_priors(currcol+1,nsamp,nmod,models_num,nalt_ls,nalt_s);
	}
	return;
}

/*function to generate all model structures in an efficient manner
and calculate prior information for use in efficient search routine*/
SEXP genmodels_priors(SEXP nsamp1)
{
	/*'nsamp1' is the number of samples*/
	
	//declare variables
	int i,j,k,l,m,q,r,nsamp,nrem,nstart,ncomb,ncomb1,valid,totmods;
	//extract integer parts of 'nsamp1'
	nsamp1 = coerceVector(nsamp1, INTSXP);
	nsamp=INTEGER(nsamp1)[0];
	
	/*this is arbitrary number of initial models
	(can be altered if necessary during runtime)*/
	int maxmodels,incmodels=10000;
	maxmodels=incmodels;
	/*set up pointer to array of pointers, each of which will point at an
	array holding relevant model structures (to be filled in during runtime)*/
	int *structure = (int *) Calloc(maxmodels*nsamp,int);
	//declare intermediate vectors for recursion
	int *indexl = (int *) Calloc(nsamp,int);
	int *indexl_ind = (int *) Calloc(nsamp,int);
	int *currstruct = (int *) Calloc(nsamp,int);
	int *comb,*comb_set,*comb_ind,*comb_next;
	
	//for fully DEPENDENT model structures
	q=0;
	for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=0;
	q++;
	
	//now complete for PARTIALLY DEPENDENT model structures
	for(k=(nsamp-1);k>1;k--)
	{
		//create set of combinations to use during looping
		ncomb = choose_c(nsamp,k);
		comb = (int *) Calloc(k,int);
		comb_set = (int *) Calloc(k*ncomb,int);
		comb_ind = (int *) Calloc(ncomb,int);
		//produce combinations
		for(j=0;j<k;j++) comb[j]=j;
		for(j=0;j<k;j++) comb_set[index2(0,j,ncomb)]=comb[j];
		i=1;
		while(next_comb(comb,k,nsamp)==1)
		{
			for(j=0;j<k;j++) comb_set[index2(i,j,ncomb)]=comb[j];
			i++;
		}
		/*calculate how many samples remaining, and if necessary
		set up for further loops*/
		nrem=nsamp-k;
		nstart=(nrem>k ? k:nrem);
		if(nstart>1)
		{
			if(nstart==k)
			{
				//now begin looping over combinations
				i=0;valid=1;
				while(i<(ncomb-1)&&valid==1)
				{
					//extract subset that can be looped over further
					check_duplicates(i,ncomb,k,comb_set,comb_ind);
					ncomb1=0;
					for(j=(i+1);j<ncomb;j++) ncomb1+=comb_ind[j];
					if(ncomb1==0) valid=0;
					else
					{
						//reset indicators
						for(j=0;j<nsamp;j++) indexl_ind[j]=1;
						//initialise structures
						for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
						for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
						//extract viable choices for next level
						comb_next = (int *) Calloc(ncomb1*k,int);
						l=0;
						for(j=(i+1);j<ncomb;j++)
						{
							if(comb_ind[j]==1)
							{
								for(m=0;m<k;m++) comb_next[index2(l,m,ncomb1)]=comb_set[index2(j,m,ncomb)];
								l++;
							}
						}
						//extract remaining samples for loop
						for(j=0;j<k;j++) indexl_ind[comb_set[index2(i,j,ncomb)]]=0;
						l=0;
						for(j=0;j<nsamp;j++)
						{
							if(indexl_ind[j]==1)
							{
								indexl[l]=j;
								l++;
							}
						}
						for(j=0;j<nsamp;j++) currstruct[j]=structure[index2(q,j,maxmodels)];
						genmodels_recur(nrem,nstart,indexl,indexl_ind,1,comb_next,ncomb1,nsamp,currstruct,&q,&structure,&maxmodels,incmodels);
						//free memory from heap
						Free(comb_next);
					}
					i++;
				}
				//now enter further recursive loop if required
				nstart--;
				if(nstart>1)
				{
					while(nstart>1)
					{
						//loop over combinations
						for(i=0;i<ncomb;i++)
						{
							//reset indicators
							for(j=0;j<nsamp;j++) indexl_ind[j]=1;
							//initialise structures
							for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
							for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
							//extract remaining samples for loop
							for(j=0;j<k;j++) indexl_ind[comb_set[index2(i,j,ncomb)]]=0;
							l=0;
							for(j=0;j<nsamp;j++)
							{
								if(indexl_ind[j]==1)
								{
									indexl[l]=j;
									l++;
								}
							}
							for(j=0;j<nsamp;j++) currstruct[j]=structure[index2(q,j,maxmodels)];
							genmodels_recur(nrem,nstart,indexl,indexl_ind,0,comb_next,0,nsamp,currstruct,&q,&structure,&maxmodels,incmodels);
						}
						nstart--;	
					}
				}
				else
				{
					//now repeat for non-repeated hierarchies - looping over combinations
					for(l=0;l<ncomb;l++)
					{
						for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
						for(j=0;j<k;j++) structure[index2(q,comb_set[index2(l,j,ncomb)],maxmodels)]=comb_set[index2(l,0,ncomb)];
						q++;
						//perform memory reallocation if required
						if(q>=maxmodels) realloc_maxmod(nsamp,&maxmodels,incmodels,&structure);
					}
				}	
			}
			else
			{
				//loop over combinations
				for(i=0;i<ncomb;i++)
				{
					//reset indicators
					for(j=0;j<nsamp;j++) indexl_ind[j]=1;
					//initialise structures
					for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
					for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
					//extract remaining samples for loop
					for(j=0;j<k;j++) indexl_ind[comb_set[index2(i,j,ncomb)]]=0;
					l=0;
					for(j=0;j<nsamp;j++)
					{
						if(indexl_ind[j]==1)
						{
							indexl[l]=j;
							l++;
						}
					}
					for(j=0;j<nsamp;j++) currstruct[j]=structure[index2(q,j,maxmodels)];
					genmodels_recur(nrem,nstart,indexl,indexl_ind,0,comb_next,0,nsamp,currstruct,&q,&structure,&maxmodels,incmodels);
				}	
			}
		}
		else
		{
			//complete model structures
			for(i=0;i<ncomb;i++)
			{
				for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
				for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
				q++;
				//perform memory reallocation if required
				if(q>=maxmodels) realloc_maxmod(nsamp,&maxmodels,incmodels,&structure);
			}
		}
		//free memory from the heap (automatically sets pointer to NULL)
		Free(comb);Free(comb_ind);Free(comb_set);
	}
	//now reallocate memory to save space
	int * temp_ptr = (int *) Calloc(q*nsamp,int);
	for(i=0;i<nsamp;i++) for(j=0;j<q;j++) temp_ptr[index2(j,i,q)]=structure[index2(j,i,maxmodels)];
	Free(structure);
	structure = temp_ptr;
	temp_ptr = NULL;
	/*now generate all possible model structures in turn, and count
	the total number of models and the number of alternative models for
	each hypothesis in turn*/
	int nalt_ls=0,nalt_s=0;
	int * models_num = (int *) Calloc(2*nsamp,int);
	r=0;
	for(i=0;i<q;i++) genfullmodels_priors(i,q,nsamp,structure,&r,models_num,&nalt_ls,&nalt_s);
	
	//now append FULLY INDEPENDENT models
	int nmod=10;
	l=r+pow_int(nmod,nsamp);
	totmods=l;
	//initialise indicator part of models_num
	for(j=0;j<nsamp;j++) models_num[j+nsamp]=j;
	//use recursive loop
	for(i=0;i<nmod;i++)
	{
		models_num[0]=i;
		genfullmodels_indep_priors(1,nsamp,nmod,models_num,&nalt_ls,&nalt_s);
	}
	//now export into R and return
	SEXP output_R;
	PROTECT(output_R = allocVector(INTSXP,3+nsamp*q));
	k=0;
	for(i=0;i<q;i++)
	{
		for(j=0;j<nsamp;j++)
		{
			INTEGER(output_R)[index2(i,j,q)] = structure[index2(i,j,q)];
			k++;
		}
	}
	INTEGER(output_R)[k]=totmods;
	INTEGER(output_R)[k+1]=nalt_ls;
	INTEGER(output_R)[k+2]=nalt_s;
	UNPROTECT(1);
	//free memory from the heap (automatically sets pointers to NULL)
	Free(structure);Free(models_num);Free(indexl);Free(indexl_ind);Free(currstruct);
	return output_R;
}

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
	int i,j,k,currhier,nsamp,nmod=10,priorind;
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

/*function to generate all possible models based on an intermediate set of structures
(returns just normalising constants for efficient EXACT routine)*/
void genfullmodels_exact(int q, int nrow, int ncol, int *structure, int *r, int *models_num, double *norm_ls, double *norm_s, double *palt_ls, double *palt_s, double lpriornull_ls, double lprioralt_ls, double lpriornull_s, double lprioralt_s, double *lPDM_int_mat, int *ncomb_sub, int *sum, int *temp_index, int *mods, int *inds, int *prec_ind)
{
	/*'q' is specific row of 'structure' we wish to extract
	'nrow' and 'ncol' are using for indexing
	'structure' is a matrix containing intermediate model structures
	'r' denotes current model
	'models_num' is a vector recording model structures
	'norm_ls' and 'norm_s' record normalising constants
	'palt_ls' and 'palt_s' record PPAs based on criteria
	'lprior...' are the log-priors for a single model for the two criteria
	'lPDM_int_mat' is intermediate vector for calculating normalising constants
	'sum', 'temp_index', 'mods' and 'inds' are all intermediate vectors for use during calculation
	'prec_ind' is an indicator to record whether there is a potential precision issue*/
	
	//declare necessary variables
	int i,j,k,l,nhier,nsamp,nmod=10,priorind;
	double temp_norm,temp_norm1;
	nsamp=ncol;
	int * temp = (int *) Calloc(nsamp,int);
	int * tempind = (int *) Calloc(nsamp,int);
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
	//now work through hierarchies in order and expand model set
	if(nhier==1)
	{
		for(i=0;i<nmod;i++)
		{
			for(j=0;j<nsamp;j++)
			{
				models_num[j]=i;
				models_num[j+nsamp]=0;
			}
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
			*r=*r+1;
		}
	}
	else
	{
		//temp[0] must be >1 in this loop
		for(j=0;j<nmod;j++)
		{
			//initialise structures
			for(k=0;k<nsamp;k++) models_num[k+nsamp]=k;
			//now fill in gaps and pass to recursive function
			for(k=0;k<nsamp;k++)
			{
				if(structure[index2(q,k,nrow)]==tempind[0])
				{
					models_num[k]=j;
					models_num[k+nsamp]=tempind[0];
				}
			}
			genfullmodels_recur_exact(nhier,0,temp,tempind,q,nrow,ncol,structure,r,models_num,norm_ls,norm_s,palt_ls,palt_s,lpriornull_ls,lprioralt_ls,lpriornull_s,lprioralt_s,lPDM_int_mat,ncomb_sub,sum,temp_index,mods,inds,prec_ind);
		}
	}
	//free memory from the heap (automatically sets pointer to NULL)
	Free(temp);Free(tempind);
	return;
}

/*function to generate independent model structure recursively
(for use in efficient EXACT routine)*/
void genfullmodels_indep_exact(int currcol, int nsamp, int nmod, int *models_num, double *norm_ls, double *norm_s, double *palt_ls, double *palt_s, double lpriornull_ls, double lprioralt_ls, double lpriornull_s, double lprioralt_s, double *lPDM_int_mat, int *ncomb_sub, int *sum, int *temp_index, int *mods, int *inds, int *prec_ind)
{
	/*'currcol' is the current sample in 'models_num'
	'nsamp' is the total number of samples of 'models_num'
	'nmod' is the number of model choices
	'models_num' is a vector recording model structures
	'norm_ls' and 'norm_s' record normalising constants
	'palt_ls' and 'palt_s' record PPAs based on criteria
	'lprior...' are the log-priors for a single model for the two criteria
	'lPDM_int_mat' is intermediate vector for calculating normalising constants
	'sum', 'temp_index', 'mods' and 'inds' are all intermediate vectors for use during calculation
	'prec_ind' is an indicator to record whether there is a potential precision issue*/
		
	int i,priorind;
	double temp_norm,temp_norm1;
	
	for(i=0;i<nmod;i++)
	{
		models_num[currcol]=i;
		if(currcol==(nsamp-1))
		{
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
		else genfullmodels_indep_exact(currcol+1,nsamp,nmod,models_num,norm_ls,norm_s,palt_ls,palt_s,lpriornull_ls,lprioralt_ls,lpriornull_s,lprioralt_s,lPDM_int_mat,ncomb_sub,sum,temp_index,mods,inds,prec_ind);
	}
	return;
}

/*function to generate all model structures in an efficient manner
and calculate normalising for use in EXACT efficient search routine*/
SEXP genmodels_exact(SEXP nsamp1, SEXP ntotcol1, SEXP lpriors, SEXP lPDM_int_mat1)
{
	/*'nsamp1' is the number of samples
	'ntotcol1' is number of columns of 'lPDM_int_mat1'
	'lpriors' is a vector of log-priors in order c(null_ls,alt_ls,null_s,alt_s)
	'lPDM_int_mat1' is a matrix for intermediate values for use in
		calculating normalising constant*/
	
	//declare variables
	int i,j,k,l,m,q,r,nsamp,nrem,nstart,ncomb,ncomb1,valid,totmods,nmod=10,ntotcol;
	double lpriornull_ls,lpriornull_s,lprioralt_ls,lprioralt_s;
	//extract integer parts of 'nsamp1'
	nsamp1 = coerceVector(nsamp1, INTSXP);
	nsamp=INTEGER(nsamp1)[0];
	//extract integer parts of 'ntotcol1'
	ntotcol1 = coerceVector(ntotcol1, INTSXP);
	ntotcol=INTEGER(ntotcol1)[0];
	//extract log-prior information
	lpriors = coerceVector(lpriors, REALSXP);
	lpriornull_ls=REAL(lpriors)[0];
	lprioralt_ls=REAL(lpriors)[1];
	lpriornull_s=REAL(lpriors)[2];
	lprioralt_s=REAL(lpriors)[3];
	//extract intermediate vector
	double *lPDM_int_mat = (double *) R_alloc(nmod*ntotcol,sizeof(double));
	lPDM_int_mat1 = coerceVector(lPDM_int_mat1, REALSXP);
	for(i=0;i<(nmod*ntotcol);i++) lPDM_int_mat[i]=REAL(lPDM_int_mat1)[i];
	
	//generate intermediate vectors for calculations
	int *ncomb_sub = (int *) R_alloc(nsamp,sizeof(int));
	ncomb_sub[0]=0;
	for(i=1;i<nsamp;i++) ncomb_sub[i]=choose_c(nsamp,i);
	for(i=1;i<nsamp;i++) ncomb_sub[i]=ncomb_sub[i]+ncomb_sub[i-1];
	//now declare intermediate vectors for model simplification
	int *sum = (int *) R_alloc(nsamp,sizeof(int));
	int *temp_index = (int *) R_alloc(nsamp,sizeof(int));
	int *mods = (int *) R_alloc(nsamp,sizeof(int));
	int *inds = (int *) R_alloc(nsamp,sizeof(int));
	
	/*this is arbitrary number of initial models
	(can be altered if necessary during runtime)*/
	int maxmodels,incmodels=10000;
	maxmodels=incmodels;
	/*set up pointer to array of pointers, each of which will point at an
	array holding relevant model structures (to be filled in during runtime)*/
	int *structure = (int *) Calloc(maxmodels*nsamp,int);
	//declare intermediate vectors for recursion
	int *indexl = (int *) Calloc(nsamp,int);
	int *indexl_ind = (int *) Calloc(nsamp,int);
	int *currstruct = (int *) Calloc(nsamp,int);
	int *comb,*comb_set,*comb_ind,*comb_next;
	
	//for fully DEPENDENT model structures
	q=0;
	for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=0;
	q++;
	
	//now complete for PARTIALLY DEPENDENT model structures
	for(k=(nsamp-1);k>1;k--)
	{
		//create set of combinations to use during looping
		ncomb = choose_c(nsamp,k);
		comb = (int *) Calloc(k,int);
		comb_set = (int *) Calloc(k*ncomb,int);
		comb_ind = (int *) Calloc(ncomb,int);
		//produce combinations
		for(j=0;j<k;j++) comb[j]=j;
		for(j=0;j<k;j++) comb_set[index2(0,j,ncomb)]=comb[j];
		i=1;
		while(next_comb(comb,k,nsamp)==1)
		{
			for(j=0;j<k;j++) comb_set[index2(i,j,ncomb)]=comb[j];
			i++;
		}
		/*calculate how many samples remaining, and if necessary
		set up for further loops*/
		nrem=nsamp-k;
		nstart=(nrem>k ? k:nrem);
		if(nstart>1)
		{
			if(nstart==k)
			{
				//now begin looping over combinations
				i=0;valid=1;
				while(i<(ncomb-1)&&valid==1)
				{
					//extract subset that can be looped over further
					check_duplicates(i,ncomb,k,comb_set,comb_ind);
					ncomb1=0;
					for(j=(i+1);j<ncomb;j++) ncomb1+=comb_ind[j];
					if(ncomb1==0) valid=0;
					else
					{
						//reset indicators
						for(j=0;j<nsamp;j++) indexl_ind[j]=1;
						//initialise structures
						for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
						for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
						//extract viable choices for next level
						comb_next = (int *) Calloc(ncomb1*k,int);
						l=0;
						for(j=(i+1);j<ncomb;j++)
						{
							if(comb_ind[j]==1)
							{
								for(m=0;m<k;m++) comb_next[index2(l,m,ncomb1)]=comb_set[index2(j,m,ncomb)];
								l++;
							}
						}
						//extract remaining samples for loop
						for(j=0;j<k;j++) indexl_ind[comb_set[index2(i,j,ncomb)]]=0;
						l=0;
						for(j=0;j<nsamp;j++)
						{
							if(indexl_ind[j]==1)
							{
								indexl[l]=j;
								l++;
							}
						}
						for(j=0;j<nsamp;j++) currstruct[j]=structure[index2(q,j,maxmodels)];
						genmodels_recur(nrem,nstart,indexl,indexl_ind,1,comb_next,ncomb1,nsamp,currstruct,&q,&structure,&maxmodels,incmodels);
						//free memory from heap
						Free(comb_next);
					}
					i++;
				}
				//now enter further recursive loop if required
				nstart--;
				if(nstart>1)
				{
					while(nstart>1)
					{
						//loop over combinations
						for(i=0;i<ncomb;i++)
						{
							//reset indicators
							for(j=0;j<nsamp;j++) indexl_ind[j]=1;
							//initialise structures
							for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
							for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
							//extract remaining samples for loop
							for(j=0;j<k;j++) indexl_ind[comb_set[index2(i,j,ncomb)]]=0;
							l=0;
							for(j=0;j<nsamp;j++)
							{
								if(indexl_ind[j]==1)
								{
									indexl[l]=j;
									l++;
								}
							}
							for(j=0;j<nsamp;j++) currstruct[j]=structure[index2(q,j,maxmodels)];
							genmodels_recur(nrem,nstart,indexl,indexl_ind,0,comb_next,0,nsamp,currstruct,&q,&structure,&maxmodels,incmodels);
						}
						nstart--;	
					}
				}
				else
				{
					//now repeat for non-repeated hierarchies - looping over combinations
					for(l=0;l<ncomb;l++)
					{
						for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
						for(j=0;j<k;j++) structure[index2(q,comb_set[index2(l,j,ncomb)],maxmodels)]=comb_set[index2(l,0,ncomb)];
						q++;
						//perform memory reallocation if required
						if(q>=maxmodels) realloc_maxmod(nsamp,&maxmodels,incmodels,&structure);
					}
				}	
			}
			else
			{
				//loop over combinations
				for(i=0;i<ncomb;i++)
				{
					//reset indicators
					for(j=0;j<nsamp;j++) indexl_ind[j]=1;
					//initialise structures
					for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
					for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
					//extract remaining samples for loop
					for(j=0;j<k;j++) indexl_ind[comb_set[index2(i,j,ncomb)]]=0;
					l=0;
					for(j=0;j<nsamp;j++)
					{
						if(indexl_ind[j]==1)
						{
							indexl[l]=j;
							l++;
						}
					}
					for(j=0;j<nsamp;j++) currstruct[j]=structure[index2(q,j,maxmodels)];
					genmodels_recur(nrem,nstart,indexl,indexl_ind,0,comb_next,0,nsamp,currstruct,&q,&structure,&maxmodels,incmodels);
				}	
			}
		}
		else
		{
			//complete model structures
			for(i=0;i<ncomb;i++)
			{
				for(j=0;j<nsamp;j++) structure[index2(q,j,maxmodels)]=j;
				for(j=0;j<k;j++) structure[index2(q,comb_set[index2(i,j,ncomb)],maxmodels)]=comb_set[index2(i,0,ncomb)];
				q++;
				//perform memory reallocation if required
				if(q>=maxmodels) realloc_maxmod(nsamp,&maxmodels,incmodels,&structure);
			}
		}
		//free memory from the heap (automatically sets pointer to NULL)
		Free(comb);Free(comb_ind);Free(comb_set);
	}
	//now reallocate memory to save space
	int * temp_ptr = (int *) Calloc(q*nsamp,int);
	for(i=0;i<nsamp;i++) for(j=0;j<q;j++) temp_ptr[index2(j,i,q)]=structure[index2(j,i,maxmodels)];
	Free(structure);
	structure = temp_ptr;
	temp_ptr = NULL;
	/*now generate all possible model structures in turn, and count
	the total number of models and the number of alternative models for
	each hypothesis in turn and figure out partial normalising constant*/
	double norm_ls=0.0,norm_s=0.0, palt_ls=0.0, palt_s=0.0;
	int * models_num = (int *) Calloc(2*nsamp,int);
	int prec_ind=0;
	r=0;
	for(i=0;i<q;i++) genfullmodels_exact(i,q,nsamp,structure,&r,models_num,&norm_ls,&norm_s,&palt_ls,&palt_s,lpriornull_ls,lprioralt_ls,lpriornull_s,lprioralt_s,lPDM_int_mat,ncomb_sub,sum,temp_index,mods,inds,&prec_ind);
	
	//now append FULLY INDEPENDENT models
	l=r+pow_int(nmod,nsamp);
	totmods=l;
	//initialise indicator part of models_num
	for(j=0;j<nsamp;j++) models_num[j+nsamp]=j;
	//use recursive loop
	for(i=0;i<nmod;i++)
	{
		models_num[0]=i;
		genfullmodels_indep_exact(1,nsamp,nmod,models_num,&norm_ls,&norm_s,&palt_ls,&palt_s,lpriornull_ls,lprioralt_ls,lpriornull_s,lprioralt_s,lPDM_int_mat,ncomb_sub,sum,temp_index,mods,inds,&prec_ind);
	}
	//now export into R and return
	SEXP output_R;
	PROTECT(output_R = allocVector(REALSXP,5));
	REAL(output_R)[0]=norm_ls;
	REAL(output_R)[1]=norm_s;
	REAL(output_R)[2]=palt_ls/norm_ls;
	REAL(output_R)[3]=palt_s/norm_s;
	REAL(output_R)[4]=prec_ind;
	UNPROTECT(1);
	//free memory from the heap (automatically sets pointers to NULL)
	Free(structure);Free(models_num);Free(indexl);Free(indexl_ind);Free(currstruct);
	return output_R;
}

//function to perform memory reallocation of full model set if required
void realloc_approx(int nmodcol, int *totmods, int incmods, int **models_num, int **hyp, double **lPPA_mat)
{
	/*'nmodcol' is number of columns of 'models_num'
	'totmods' records maximum number of models
	'incmods' is value with which to increment 'maxmodels'
	'models_num', 'hyp' and 'lPPA_mat' are matrices
		(and vectors) requiring reallocation*/
		
	*totmods=(*totmods)+incmods;
	*models_num = (int *) Realloc(*models_num,(*totmods)*nmodcol,int);
	*hyp = (int *) Realloc(*hyp,*totmods,int);
	*lPPA_mat = (double *) Realloc(*lPPA_mat,*totmods,double);
	return;
}

//function to perform memory reallocation of full model set if required
void realloc_approx_new(int nmodcol, int *totmods, int incmods, int **models_num, int **hyp, double **lPPA_mat, int **multfact)
{
	/*'nmodcol' is number of columns of 'models_num'
	'totmods' records maximum number of models
	'incmods' is value with which to increment 'maxmodels'
	'models_num', 'hyp', 'lPPA_mat' and 'multfact' are matrices
		(and vectors) requiring reallocation*/
		
	*totmods=(*totmods)+incmods;
	*models_num = (int *) Realloc(*models_num,(*totmods)*nmodcol,int);
	*hyp = (int *) Realloc(*hyp,*totmods,int);
	*lPPA_mat = (double *) Realloc(*lPPA_mat,*totmods,double);
	*multfact = (int *) Realloc(*multfact,*totmods,int);
	return;
}

//function to calculate correct column of lPDM_int_mat based on a single set of changes
int calc_corr_col(int q, int nrow, int *structure, int nsamp, int *ncomb_sub, int samp,  int sum, int *temp_index)
{
	/*'q' is relevant row of 'structure'
	'nrow' is number of rows of 'structure'
	'structure' is a vector recording dependency indicators
	'nsamp' is number of samples
	'ncomb_sub' is a vector recording numbers of combinations (to be used for selecting
		correct columns of 'lPDM_int_mat'
	'samp' corresponds to which indicators to search for
	'sum' corresponds to how many indicators of type 'samp' there are 
	'temp_index' is an intermediate vector for use during calculation*/
	
	int l,m;

	//now generate correct column indicators
	l=0;
	m=0;
	while(l<nsamp&&m<sum)
	{
		if(structure[index2(q,l,nrow)]==samp)
		{
			temp_index[m]=l;
			m++;
		}
		l++;
	}
	//calculate correct column
	m=0;
	for(l=0;l<sum;l++) m+=choose_c(temp_index[l],l+1);
	m+=ncomb_sub[sum-1];
	return m;
}

//function to generate top subset of fully and patially dependent models based on an intermediate set of structures
void genfullmodels_approx_recur(int nhier, int currhier1, int *temp, int *tempind, int *currmods, int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int nmodcol, double *lPDM_int_mat, int *lPDM_int_ind, int ntotcol, int **models_num, int **hyp, double **lPPA_mat, double logthresh, int *ncomb_sub, double curr_lPPA, double *max_lPPA)
{
	/*'nhier' is number of hierarchies
	'currhier1' is current hierarchy
	'temp' and 'tempind' are vectors recording structural information
	'currmods' is model structure for current model
	'q' is specific row of 'structure' we wish to extract
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
	'ncomb_sub' is intermediate vector used for indexing
	'curr_lPPA' is current log-PDM at previous hierarchy
	'max_lPPA' is current maximum log-PDM*/
	
	//declare necessary variables
	int i,j,k,l,m,nsamp,nmod=10,valid;
	nsamp=ncol;
		
	//work out variables relating to maximum independent model
	double temp_lPPA;
	int corrcol;
	int * currmods1 = (int *) Calloc(nmodcol,int);
	int * temp_index = (int *) Calloc(nsamp,int);
	
	int currhier=currhier1+1;
	//if only independent structures left to evaluate
	if(temp[currhier]==1)
	{
		corrcol=tempind[currhier];
		i=0;valid=0;
		while(i<nmod&&valid==0)
		{
			//reset model structure
			for(j=0;j<nmodcol;j++) (*models_num)[index2_col(*r,j,nmodcol)]=currmods[j];
			(*models_num)[index2_col(*r,corrcol,nmodcol)]=lPDM_int_ind[index2(i,corrcol,nmod)];
			//change model
			temp_lPPA=curr_lPPA-lPDM_int_mat[index2(0,corrcol,nmod)]+lPDM_int_mat[index2(i,corrcol,nmod)];
			//if final choice then record model and update system
			if(currhier==(nhier-1))
			{
				//check if new model is selected
				if(((*max_lPPA)-temp_lPPA)<logthresh)
				{
					(*lPPA_mat)[*r]=temp_lPPA;
					*r=*r+1;
					//increase size of output vectors if necessary
					if(*r>(*totmods-1)) realloc_approx(nmodcol,totmods,incmods,models_num,hyp,lPPA_mat);
					/*if new model is maximum model then recalculate log-acceptance threshold then reset maximum model*/
					if(temp_lPPA>(*max_lPPA)) *max_lPPA=temp_lPPA;
				}
				else valid=1;
			}
			else
			{
				for(j=0;j<nmodcol;j++) currmods1[j]=(*models_num)[index2_col(*r,j,nmodcol)];
				//enter further hierarchy
				genfullmodels_approx_recur(nhier,currhier,temp,tempind,currmods1,q,nrow,ncol,structure,r,totmods,incmods,nmodcol,lPDM_int_mat,lPDM_int_ind,ntotcol,models_num,hyp,lPPA_mat,logthresh,ncomb_sub,temp_lPPA,max_lPPA);
			}
			i++;
		}
	}
	else
	{
		//calculate correct column of lPDM_int_mat based on structure of current hierarchy
		corrcol=calc_corr_col(q,nrow,structure,nsamp,ncomb_sub,tempind[currhier],temp[currhier],temp_index);
		//calculate PPAs
		for(i=0;i<nmod;i++)
		{
/*			if(lPDM_int_ind[index2(i,corrcol,nmod)]>0)*/
/*			{*/
				//reset model structure
				for(j=0;j<nmodcol;j++) (*models_num)[index2_col(*r,j,nmodcol)]=currmods[j];
				//now fill in gaps and initialise calculation of changes to lPDM
				temp_lPPA=curr_lPPA;
				for(k=0;k<nsamp;k++)
				{
					if(structure[index2(q,k,nrow)]==tempind[currhier])
					{
						(*models_num)[index2_col(*r,k,nmodcol)]=lPDM_int_ind[index2(i,corrcol,nmod)];
						(*models_num)[index2_col(*r,k+nsamp,nmodcol)]=tempind[currhier];
						temp_lPPA-=lPDM_int_mat[index2(0,k,nmod)];
					}
				}
				//adjust lPDM to correct structure
				temp_lPPA+=lPDM_int_mat[index2(i,corrcol,nmod)];
				//if final choice then record model and update system
				if(currhier==(nhier-1))
				{
					//check if new model is selected
					if(((*max_lPPA)-temp_lPPA)<logthresh)
					{
						(*lPPA_mat)[*r]=temp_lPPA;
						*r=*r+1;
						//increase size of output vectors if necessary
						if(*r>(*totmods-1)) realloc_approx(nmodcol,totmods,incmods,models_num,hyp,lPPA_mat);
						/*if new model is maximum model then recalculate log-acceptance threshold then reset maximum model*/
						if(temp_lPPA>(*max_lPPA)) *max_lPPA=temp_lPPA;
					}
				}
				else
				{
					for(j=0;j<nmodcol;j++) currmods1[j]=(*models_num)[index2_col(*r,j,nmodcol)];
					//enter further hierarchy
					genfullmodels_approx_recur(nhier,currhier,temp,tempind,currmods1,q,nrow,ncol,structure,r,totmods,incmods,nmodcol,lPDM_int_mat,lPDM_int_ind,ntotcol,models_num,hyp,lPPA_mat,logthresh,ncomb_sub,temp_lPPA,max_lPPA);
				}
/*			}*/
		}
	}
	//free memory from the heap (automatically sets pointers to NULL)
	Free(currmods1);Free(temp_index);
	return;
}

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
	int i,j,k,l,m,nhier,nsamp,nmod=10,valid;
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
		
	int i,j,k,l,m,q,d,temp_hyp,temp_hyp1;
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

//function to calculate "top" set of models for use in approximation
SEXP calc_PPAs_approx_fn(SEXP nsamp1, SEXP ntotcol1, SEXP logc1, SEXP lpriornull1, SEXP lprioralt1, SEXP lPDM_int_mat1, SEXP string1, SEXP structure1, SEXP nstruct1, SEXP return_final1, SEXP uni1, SEXP uni_ind1)
{
	/*'nsamp1' is the number of samples
	'ntotcol1' is total number of columns in 'lPDM_int_mat'
	'logc1' is the log-threshold value for acceptance relative to the maximum model
	'lpriornull1' is the log[P(M)] for a single 'null' model
	'lprioralt1' is the log[P(M)] for a single 'alternative' model
	'lPDM_int_mat1' is a matrix of log[P'(D|M)] values for use in intermediary calculations
		for calculating independent model structures
	'string1' is a binary indicator denoting prior criterion, such that
		0=less stringent and 1=stringent
	'structure1' is a matrix (in vector form) containing structures generated from 
		'genmodels_priors'
	'nstruct1' is number of structures in 'structure'
	'return_final' is a binary value taking the value 1 if model sets are to be returned and 0 otherwise
	'uni1' and 'uni_ind1' are matrices (in vector form) for use in producing PPAs without recording the full model sets*/
	
	/*DECLARE VARIABLES*/
	int i,j,k,l,m,q,d,r,nmod=10,nsamp,incmods,totmods,ntotcol,nmodcol,string,nstruct,currsamp,valid,return_final;
	double lpriornull,lprioralt,lpriordiff,temp_lPPA,max_lPPA,logc,thresh;
	//extract relevant parts of inputs for use in C code
	nsamp1 = coerceVector(nsamp1, INTSXP);
	nsamp=INTEGER(nsamp1)[0];
	ntotcol1 = coerceVector(ntotcol1, INTSXP);
	ntotcol=INTEGER(ntotcol1)[0];
	string1 = coerceVector(string1, INTSXP);
	string=INTEGER(string1)[0];
	nstruct1 = coerceVector(nstruct1, INTSXP);
	nstruct=INTEGER(nstruct1)[0];
	logc1 = coerceVector(logc1, REALSXP);
	logc=REAL(logc1)[0];
	lpriornull1 = coerceVector(lpriornull1, REALSXP);
	lpriornull=REAL(lpriornull1)[0];
	lprioralt1 = coerceVector(lprioralt1, REALSXP);
	lprioralt=REAL(lprioralt1)[0];
	return_final1 = coerceVector(return_final1, INTSXP);
	return_final=INTEGER(return_final1)[0];
	nmodcol=nsamp*2;
	incmods=10000;
	totmods=incmods;
	//DECLARE AND INITIALISE FIXED-SIZE ARRAYS
	int *structure = (int *) R_alloc(nsamp*nstruct,sizeof(int));
	structure1 = coerceVector(structure1, INTSXP);
	for(i=0;i<nstruct;i++) for(j=0;j<nsamp;j++) structure[index2(i,j,nstruct)]=INTEGER(structure1)[index2(i,j,nstruct)];
	double *lPDM_int_mat = (double *) R_alloc(nmod*ntotcol,sizeof(double));
	lPDM_int_mat1 = coerceVector(lPDM_int_mat1, REALSXP);
	for(i=0;i<(nmod*ntotcol);i++) lPDM_int_mat[i]=REAL(lPDM_int_mat1)[i];
	int *uni = (int *) R_alloc(nsamp*nmod*nmod,sizeof(int));
	uni1 = coerceVector(uni1, INTSXP);
	for(i=0;i<(nsamp*nmod*nmod);i++) uni[i]=INTEGER(uni1)[i];
	int *uni_ind = (int *) R_alloc(nsamp*nmod,sizeof(int));
	uni_ind1 = coerceVector(uni_ind1, INTSXP);
	for(i=0;i<(nsamp*nmod);i++) uni_ind[i]=INTEGER(uni_ind1)[i];
	//print to screen as a check
/*	Rprintf("uni\n");*/
/*	for(i=0;i<nmod;i++)*/
/*	{*/
/*		for(j=0;j<nsamp;j++) for(k=0;k<nmod;k++) Rprintf("%d\t",uni[index3(i,k,j,nmod,nmod)]);*/
/*		Rprintf("\n");*/
/*	}*/
/*	Rprintf("uni_ind\n");*/
/*	for(i=0;i<nmod;i++)*/
/*	{*/
/*		for(j=0;j<nsamp;j++) Rprintf("%d\t",uni_ind[index2(i,j,nmod)]);*/
/*		Rprintf("\n");*/
/*	}*/
	
	//declare vectors for intermediate calculations
	int *lPDM_int_ind = (int *) R_alloc(nmod*ntotcol,sizeof(int));
	for(i=0;i<(nmod*ntotcol);i++) lPDM_int_ind[i]=0;
	int *currmods = (int *) R_alloc(nsamp,sizeof(int));
	int *currmods1 = (int *) R_alloc(nsamp,sizeof(int));
	int *currlocs = (int *) R_alloc(nsamp,sizeof(int));
	int *currlocs1 = (int *) R_alloc(nsamp,sizeof(int));
	int *maxmods = (int *) R_alloc(nsamp,sizeof(int));
	int *maxlocs = (int *) R_alloc(nsamp,sizeof(int));
	int *comb = (int *) R_alloc(nsamp,sizeof(int));
	//now generate intermediate vector for model simplification
	int *temp_index = (int *) R_alloc(nsamp,sizeof(int));
	/*generate vector containing number of combinations in each subset
	(for generating column indices for use in model simplification)*/
	int *ncomb_sub = (int *) R_alloc(nsamp+1,sizeof(int));
	ncomb_sub[0]=0;
	for(i=1;i<=nsamp;i++) ncomb_sub[i]=choose_c(nsamp,i);
	for(i=1;i<=nsamp;i++) ncomb_sub[i]=ncomb_sub[i]+ncomb_sub[i-1];
	
	//DECLARE AND INITIALISE VARIABLE-SIZED ARRAYS
	//generate model output vector
	double *lPPA_mat = (double *) Calloc(totmods,double);
	int *hyp = (int *) Calloc(totmods,int);
	int *models_num = (int *) Calloc(nmodcol*totmods,int);
	
	//order lPDM matrix and store correctly ordered indices
	for(i=0;i<nmod;i++) for(j=0;j<ntotcol;j++) lPDM_int_ind[index2(i,j,nmod)]=i;
	
	for(i=0;i<ntotcol;i++) bubble_sort_dec(&lPDM_int_mat[index2(0,i,nmod)],&lPDM_int_ind[index2(0,i,nmod)],nmod);
	//record current maximum model
	r=0;
	max_lPPA=0.0;
	for(j=0;j<nsamp;j++)
	{
		max_lPPA+=lPDM_int_mat[index2(0,j,nmod)];
		models_num[index2_col(r,j,nmodcol)]=lPDM_int_ind[index2(0,j,nmod)];
		models_num[index2_col(r,j+nsamp,nmodcol)]=j;
	}
	lPPA_mat[r]=max_lPPA;
	//set up vectors to record current maximum model plus locations
	for(j=0;j<nsamp;j++)
	{
		maxmods[j]=models_num[index2_col(r,j,nmodcol)];
		maxlocs[j]=0;
	}
	//update model counter
	r++;
	
	//now calculate maximum effect of changing prior
	lpriordiff=lpriornull-lprioralt;
	if(lpriordiff<0) lpriordiff=-lpriordiff;
	//update threshold
	thresh=logc+lpriordiff;
		
	//run recursion through to calculate maximum fully independent model
	i=0;
	valid=1;
	currsamp=0;
	while(i<10&&valid==1)
	{
		//update model
		for(j=0;j<nsamp;j++) models_num[index2_col(r,j,nmodcol)]=maxmods[j];
		models_num[index2_col(r,currsamp,nmodcol)]=lPDM_int_ind[index2(i,currsamp,nmod)];
		if(models_num[index2_col(r,currsamp,nmodcol)]==maxmods[currsamp])
		{
			//enter into recursion
			for(j=0;j<nsamp;j++)
			{
				currmods[j]=maxmods[j];
				currlocs[j]=maxlocs[j];
			}
			currlocs[currsamp]=i;
			if((currsamp+1)<nsamp) calc_approx_indep_recur(nsamp,nmod,nmodcol,incmods,&totmods,currsamp+1,&r,currmods,currlocs,maxmods,maxlocs,&max_lPPA,max_lPPA,thresh,&models_num,&hyp,lPDM_int_mat,lPDM_int_ind,&lPPA_mat);
			i++;
		}
		else
		{
			//calculate new log-PDM
			temp_lPPA=max_lPPA-lPDM_int_mat[index2(maxlocs[currsamp],currsamp,nmod)]+lPDM_int_mat[index2(i,currsamp,nmod)];
			//add to lPPA_mat output vector
			lPPA_mat[r]=temp_lPPA;
			//if new model is maximum, then update maximum model
			if(temp_lPPA<=max_lPPA)
			{
				//if new model is accepted then update model set
				if((max_lPPA-temp_lPPA)<thresh)
				{
					//increase model counter
					r++;
					//increase size of output vectors if necessary
					if(r>(totmods-1)) realloc_approx(nmodcol,&totmods,incmods,&models_num,&hyp,&lPPA_mat);
					//enter into recursion
					for(j=0;j<nsamp;j++)
					{
						currmods[j]=models_num[index2_col(r-1,j,nmodcol)];
						currlocs[j]=maxlocs[j];
					}
					currlocs[currsamp]=i;
					if((currsamp+1)<nsamp) calc_approx_indep_recur(nsamp,nmod,nmodcol,incmods,&totmods,currsamp+1,&r,currmods,currlocs,maxmods,maxlocs,&max_lPPA,temp_lPPA,thresh,&models_num,&hyp,lPDM_int_mat,lPDM_int_ind,&lPPA_mat);
				
				}
				else valid=0;
			}
			else
			{
				//increase model counter
				r++;
				//increase size of output vectors if necessary
				if(r>(totmods-1)) realloc_approx(nmodcol,&totmods,incmods,&models_num,&hyp,&lPPA_mat);
				//update maximum model
				max_lPPA=temp_lPPA;
				//set up vectors to record current maximum model plus locations
				maxmods[currsamp]=models_num[index2_col(r-1,currsamp,nmodcol)];
				maxlocs[currsamp]=i;
				//enter into recursion
				for(j=0;j<nsamp;j++)
				{
					currmods[j]=models_num[index2_col(r-1,j,nmodcol)];;
					currlocs[j]=maxlocs[j];
				}
				if((currsamp+1)<nsamp) calc_approx_indep_recur(nsamp,nmod,nmodcol,incmods,&totmods,currsamp+1,&r,currmods,currlocs,maxmods,maxlocs,&max_lPPA,temp_lPPA,thresh,&models_num,&hyp,lPDM_int_mat,lPDM_int_ind,&lPPA_mat);
			}
			i++;
		}
	}
		
	//sort all independent models into rank order according to log(PPA')
	int *temp_ind = (int *) Calloc(r,int);
	bubble_sort_dec(lPPA_mat,temp_ind,r);
	//generate sorted vectors/matrices
	int * new_models_num = (int *) Calloc(nmodcol*totmods,int);
	for(i=0;i<r;i++) for(j=0;j<nmodcol;j++) new_models_num[index2_col(i,j,nmodcol)]=models_num[index2_col(temp_ind[i],j,nmodcol)];
	Free(models_num);
	models_num = new_models_num;
	new_models_num = NULL;
	Free(temp_ind);
	temp_ind=NULL;
	
	
/*	for(i=0;i<r;i++)*/
/*	{*/
/*		for(j=0;j<nmodcol;j++) Rprintf("%d\t",models_num[index2_col(i,j,nmodcol)]);*/
/*		Rprintf("%f\n",lPPA_mat[i]);*/
/*	}*/
/*	Rprintf("return_final=%d\n",return_final);*/
/*	getchar();*/
		
	//fill in dependency structures
	for(i=0;i<r;i++) for(j=0;j<nsamp;j++) models_num[index2_col(i,j+nsamp,nmodcol)]=j;
	
	//now generate all viable dependent structures relative to maximum INDEPENDENT model
	for(i=0;i<nstruct;i++) genfullmodels_approx(i,nstruct,nsamp,structure,&r,&totmods,incmods,nmodcol,lPDM_int_mat,lPDM_int_ind,ntotcol,&models_num,&hyp,&lPPA_mat,thresh,ncomb_sub);
	
/*	for(i=0;i<r;i++)*/
/*	{*/
/*		for(j=0;j<nmodcol;j++) Rprintf("%d\t",models_num[index2_col(i,j,nmodcol)]);*/
/*		Rprintf("%f\n",lPPA_mat[i]);*/
/*	}*/
/*	Rprintf("return_final=%d\n",return_final);*/
/*	getchar();*/
	
/*	for(i=0;i<r;i++) Rprintf("%f\n",lPPA_mat[i]);*/
/*	Rprintf("\n");*/
/*	Rprintf("lpriornull=%f\tlprioralt=%f\n",lpriornull,lprioralt);*/
	
	//add prior information
	if(string==0)
	{
		for(i=0;i<r;i++)
		{
			hyp[i]=calc_lessstring_single(nsamp,&models_num[index2_col(i,0,nmodcol)]);
			lPPA_mat[i]+=(hyp[i]==0 ? lpriornull:lprioralt);
		}
	}
	else
	{
		for(i=0;i<r;i++)
		{
			hyp[i]=calc_string_single(nsamp,&models_num[index2_col(i,0,nmodcol)]);
			lPPA_mat[i]+=(hyp[i]==0 ? lpriornull:lprioralt);
		}
	}
	
	//sort all models into rank order according to log(PPA')
	temp_ind = (int *) Calloc(r,int);
	bubble_sort_dec(lPPA_mat,temp_ind,r);
/*	for(i=0;i<r;i++) Rprintf("%f\n",lPPA_mat[i]);*/
/*	Rprintf("\n");*/
	//generate sorted vectors/matrices
	new_models_num = (int *) Calloc(nmodcol*totmods,int);
	for(i=0;i<r;i++) for(j=0;j<nmodcol;j++) new_models_num[index2_col(i,j,nmodcol)]=models_num[index2_col(temp_ind[i],j,nmodcol)];
	Free(models_num);
	models_num = new_models_num;
	new_models_num = NULL;
	
	int *new_hyp = (int *) Calloc(totmods,int);
	for(i=0;i<r;i++) new_hyp[i]=hyp[temp_ind[i]];
	Free(hyp);
	hyp = new_hyp;
	new_hyp = NULL;
	Free(temp_ind);
	
	//remove any below threshold
	if(return_final==1)
	{
		i=1;
		while(i<r&&(lPPA_mat[0]-lPPA_mat[i])<logc) i++;
		r=i;
	}
	else
	{	
		i=1;
		while(i<r&&(lPPA_mat[0]-lPPA_mat[i])<thresh) i++;
		r=i;
	}
/*	for(i=0;i<nmod;i++)*/
/*	{*/
/*		for(j=0;j<nsamp;j++) Rprintf("%f\t",lPDM_int_mat[index2(i,j,nmod)]);*/
/*		Rprintf("\n");*/
/*	}*/
/*	for(i=0;i<r;i++)*/
/*	{*/
/*		for(j=0;j<nmodcol;j++) Rprintf("%d\t",models_num[index2_col(i,j,nmodcol)]);*/
/*		Rprintf("%f\n",lPPA_mat[i]);*/
/*	}*/
/*	Rprintf("return_final=%d\n",return_final);*/
/*	getchar();*/
	
	//now make adjustments to model set if output is to be suppressed
	int *multfact = (int *) Calloc(totmods,int);
	int r1=r;
	if(return_final==0)
	{
		//resort lPDM_int_mat
		for(i=0;i<totmods;i++) multfact[i]=1;
		double temp_lPPA;
		int temp_hyp;
		int *currmods = (int *) R_alloc(nmodcol,sizeof(int));
		//start loop
		for(i=0;i<nsamp;i++)
		{
			for(j=0;j<nmod;j++)
			{
				if(uni_ind[index2(j,i,nmod)]>0)
				{
					for(k=0;k<uni_ind[index2(j,i,nmod)];k++)
					{
						for(m=0;m<r;m++)
						{
							if(models_num[index2_col(m,i,nmodcol)]==j)
							{
								q=0;
								for(l=0;l<nsamp;l++) if(models_num[index2_col(m,l+nsamp,nmodcol)]==models_num[index2_col(m,i+nsamp,nmodcol)]) q++;
								if(q==1)
								{
									for(l=0;l<nmodcol;l++) currmods[l]=models_num[index2_col(m,l,nmodcol)];
									currmods[i]=uni[index3(j,k,i,nmod,nmod)];
									temp_lPPA=lPPA_mat[m];
									if(string==0) temp_hyp=calc_lessstring_single(nsamp,currmods);
									else temp_hyp=calc_string_single(nsamp,currmods);
									if(temp_hyp!=hyp[m])
									{
										if(temp_hyp==1) temp_lPPA+=(lprioralt-lpriornull);
										else temp_lPPA+=(lpriornull-lprioralt);
										//update model set
										for(l=0;l<nmodcol;l++) models_num[index2_col(r1,l,nmodcol)]=currmods[l];
										lPPA_mat[r1]=temp_lPPA;
										hyp[r1]=temp_hyp;
										multfact[r1]=1;
										r1++;
										//increase size of output vectors if necessary
										if(r1>(totmods-1)) realloc_approx_new(nmodcol,&totmods,incmods,&models_num,&hyp,&lPPA_mat,&multfact);
										//now enter further recursion if necessary
										if((i+1)<nsamp)
										{
											//just check next sample is independent
											q=i+1;
											d=2;
											while(q<nsamp&&d!=1)
											{
												d=0;
												for(l=0;l<nsamp;l++) if(currmods[l+nsamp]==currmods[q+nsamp]) d++;
												if(d==1) final_lPPA_recur(nsamp,nmod,nmodcol,incmods,&totmods,q,&r1,r1,currmods,temp_lPPA,temp_hyp,&models_num,&hyp,&lPPA_mat,&multfact,uni,uni_ind,string,lpriornull,lprioralt);
												q++;
											}
										}
									}
									else
									{
										multfact[m]++;
										//now enter further recursion if necessary
										if((i+1)<nsamp)
										{
											//just check next sample is independent
											q=i+1;
											d=2;
											while(q<nsamp&&d!=1)
											{
												d=0;
												for(l=0;l<nsamp;l++) if(currmods[l+nsamp]==currmods[q+nsamp]) d++;
												if(d==1) final_lPPA_recur(nsamp,nmod,nmodcol,incmods,&totmods,q,&r1,m,currmods,temp_lPPA,temp_hyp,&models_num,&hyp,&lPPA_mat,&multfact,uni,uni_ind,string,lpriornull,lprioralt);
												q++;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
/*		Rprintf("\n");*/
/*		for(i=0;i<r1;i++)*/
/*		{*/
/*			for(j=0;j<nmodcol;j++) Rprintf("%d\t",models_num[index2_col(i,j,nmodcol)]);*/
/*			Rprintf("%d\t%f\n",multfact[i],lPPA_mat[i]);*/
/*		}*/
/*		getchar();*/
		r=r1;
		//sort all models into rank order according to log(PPA')
		temp_ind = (int *) Calloc(r,int);
		bubble_sort_dec(lPPA_mat,temp_ind,r);
		//generate sorted vectors/matrices
		new_models_num = (int *) Calloc(nmodcol*totmods,int);
		for(i=0;i<r;i++) for(j=0;j<nmodcol;j++) new_models_num[index2_col(i,j,nmodcol)]=models_num[index2_col(temp_ind[i],j,nmodcol)];
		Free(models_num);
		models_num = new_models_num;
		new_models_num = NULL;
	
		new_hyp = (int *) Calloc(totmods,int);
		for(i=0;i<r;i++) new_hyp[i]=hyp[temp_ind[i]];
		Free(hyp);
		hyp = new_hyp;
		new_hyp = NULL;
		
		int *new_multfact = (int *) Calloc(totmods,int);
		for(i=0;i<r;i++) new_multfact[i]=multfact[temp_ind[i]];
		Free(multfact);
		multfact = new_multfact;
		new_multfact = NULL;
		Free(temp_ind);
		
		//now cycle through and remove models outside threshold
		i=1;
		while(i<r&&(lPPA_mat[0]-lPPA_mat[i])<logc) i++;
		r=i;
		
		//reallocate memory to save space
		realloc_approx_new(nmodcol,&totmods,r-totmods,&models_num,&hyp,&lPPA_mat,&multfact);
	}
	else
	{
		//reallocate memory to save space
		realloc_approx(nmodcol,&totmods,r-totmods,&models_num,&hyp,&lPPA_mat);
		for(i=0;i<r;i++) multfact[i]=1;
	}	
		
	//now export into R and return
	SEXP models_num_R;
	PROTECT(models_num_R = allocVector(REALSXP,totmods*nmodcol+3*totmods+1));
	k=0;
	for(i=0;i<totmods;i++)
	{
		for(j=0;j<nmodcol;j++)
		{
			REAL(models_num_R)[index2_col(i,j,nmodcol)] = (double) models_num[index2_col(i,j,nmodcol)];
			k++;
		}
	}
	for(j=0;j<totmods;j++)
	{
		REAL(models_num_R)[k] = (double) hyp[j];
		k++;
	}
	for(j=0;j<totmods;j++)
	{
		REAL(models_num_R)[k] = lPPA_mat[j];
		k++;
	}
	for(j=0;j<totmods;j++)
	{
		REAL(models_num_R)[k] = (double) multfact[j];
		k++;
	}
	REAL(models_num_R)[k]=(double) totmods;
	UNPROTECT(1);
	//free memory from the heap (automatically sets pointer to NULL)
	Free(models_num);Free(hyp);Free(lPPA_mat);Free(multfact);
	return models_num_R;
}


