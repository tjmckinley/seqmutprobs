#include "functions.h"
#include <R.h>

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
    int j, k, currhier, nsamp, nmod = 10;
    nsamp = ncol;
    currhier = currhier1 + 1;
    //continue working through hierarchies in order and expand model set
    if(temp[currhier] == 1)
    {
        for(j = 0; j < nmod; j++)
        {
            models_num[tempind[currhier]] = j;
            if(currhier < (nhier - 1)) genfullmodels_recur_priors(nhier, currhier, temp, tempind, q, nrow, ncol, structure, r, models_num, nalt_ls, nalt_s);
            else
            {
                *r = *r + 1;
                //now calculate prior specification based on model
                *nalt_ls = *nalt_ls + calc_lessstring_single(nsamp, models_num);
                *nalt_s = *nalt_s + calc_string_single(nsamp, models_num);
            }
        }
    }
    else
    {
        for(j = 0; j < nmod; j++)
        {
            //now fill in gaps and pass to recursive function if necessary
            for(k = 0; k < nsamp; k++)
            {
                if(structure[index2(q, k, nrow)] == tempind[currhier])
                {
                    models_num[k] = j;
                    models_num[k + nsamp] = tempind[currhier];
                }
            }
            if(currhier < (nhier - 1)) genfullmodels_recur_priors(nhier, currhier, temp, tempind, q, nrow, ncol, structure, r, models_num, nalt_ls, nalt_s);
            else
            {
                *r = *r + 1;
                //now calculate prior specification based on model
                *nalt_ls = *nalt_ls + calc_lessstring_single(nsamp, models_num);
                *nalt_s = *nalt_s + calc_string_single(nsamp, models_num);
            }
        }
    }
    return;
}

