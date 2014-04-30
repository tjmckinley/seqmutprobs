#include "functions.h"
#include <R.h>

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
    int i, j, k, nhier, nsamp, nmod = 10, priorind;
    double temp_norm, temp_norm1;
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
            temp_norm = calc_lPDMs_single_vec(nsamp, nmod, models_num, lPDM_int_mat, ncomb_sub, sum, temp_index, mods, inds);
            priorind = calc_lessstring_single(nsamp, models_num);
            temp_norm1 = exp(temp_norm + (priorind == 1 ? lprioralt_ls : lpriornull_ls));
            if(temp_norm1 == 0.0) *prec_ind = 1;
            *norm_ls = *norm_ls + temp_norm1;
            *palt_ls = *palt_ls + (priorind == 1 ? temp_norm1 : 0.0);
            priorind = calc_string_single(nsamp, models_num);
            temp_norm1 = exp(temp_norm + (priorind == 1 ? lprioralt_s : lpriornull_s));
            if(temp_norm1 == 0.0) *prec_ind = 1;
            *norm_s = *norm_s + temp_norm1;
            *palt_s = *palt_s + (priorind == 1 ? temp_norm1 : 0.0);
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
            genfullmodels_recur_exact(nhier, 0, temp, tempind, q, nrow, ncol, structure, r, models_num, norm_ls, norm_s, palt_ls, palt_s, lpriornull_ls, lprioralt_ls, lpriornull_s, lprioralt_s, lPDM_int_mat, ncomb_sub, sum, temp_index, mods, inds, prec_ind);
        }
    }
    //free memory from the heap (automatically sets pointer to NULL)
    Free(temp);
    Free(tempind);
    return;
}

