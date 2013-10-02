#include "functions.h"
#include <R.h>

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

