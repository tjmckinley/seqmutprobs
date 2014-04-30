#include "functions.h"
#include <R.h>

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

    for(i = 0; i < nmod; i++)
    {
        models_num[currcol] = i;
        if(currcol == (nsamp - 1))
        {
            //now calculate prior specification based on model
            *nalt_ls = *nalt_ls + calc_lessstring_single(nsamp, models_num);
            *nalt_s = *nalt_s + calc_string_single(nsamp, models_num);
        }
        else genfullmodels_indep_priors(currcol + 1, nsamp, nmod, models_num, nalt_ls, nalt_s);
    }
    return;
}

