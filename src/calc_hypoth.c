#include "functions.h"
#include <R.h>

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
