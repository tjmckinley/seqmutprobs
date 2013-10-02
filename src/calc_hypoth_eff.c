#include "functions.h"
#include <R.h>
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

