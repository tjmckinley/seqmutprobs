#include "functions.h"
#include <R.h>

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
    int i, j, k, nhier, nsamp, nmod = 10;
    nsamp = ncol;
    int * temp = (int *) Calloc(nsamp, int);
    int * tempind = (int *) Calloc(nsamp, int);
    //count up how many hierarchies
    for(i = 0; i < nsamp; i++)
    {
        temp[i] = 0;
        for(j = 0; j < nsamp; j++) temp[i] += (structure[index2(q, j, nrow)] == i ? 1 : 0);
    }
    nhier = 0;
    for(j = 0; j < nsamp; j++) nhier += (temp[j] > 0 ? 1 : 0);
    //order hierarchies
    bubble_sort_dec_int(temp, tempind, nsamp);
    //now work through hierarchies in order and expand model set
    if(nhier == 1)
    {
        for(i = 0; i < nmod; i++)
        {
            for(j = 0; j < nsamp; j++)
            {
                models_num[j] = i;
                models_num[j + nsamp] = 0;
            }
            //now calculate prior specification based on model
            *nalt_ls = *nalt_ls + calc_lessstring_single(nsamp, models_num);
            *nalt_s = *nalt_s + calc_string_single(nsamp, models_num);
            *r = *r + 1;
        }
    }
    else
    {
        //temp[0] must be >1 in this loop
        for(j = 0; j < nmod; j++)
        {
            //initialise structures
            for(k = 0; k < nsamp; k++) models_num[k + nsamp] = k;
            //now fill in gaps and pass to recursive function
            for(k = 0; k < nsamp; k++)
            {
                if(structure[index2(q, k, nrow)] == tempind[0])
                {
                    models_num[k] = j;
                    models_num[k + nsamp] = tempind[0];
                }
            }
            genfullmodels_recur_priors(nhier, 0, temp, tempind, q, nrow, ncol, structure, r, models_num, nalt_ls, nalt_s);
        }
    }
    //free memory from the heap (automatically sets pointer to NULL)
    Free(temp);
    Free(tempind);
    return;
}

