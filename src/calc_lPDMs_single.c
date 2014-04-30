#include "functions.h"
#include <R.h>

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

    int j, k, l, m;
    double lPDM;

    //first part of code simplifies the model structures and indicators
    for(j = 0; j < nsamp; j++)
    {
        //count number of times indicator j appears in model
        sum[j] = 0;
        for(k = 0; k < nsamp; k++) if(models_num[index2(r, k + nsamp, totmods)] == j) sum[j]++;
    }
    //now generate correct column indicators
    k = 0;
    for(j = 0; j < nsamp; j++)
    {
        if(sum[j] > 0)
        {
            l = 0;
            m = 0;
            while(l < nsamp && m < sum[j])
            {
                if(models_num[index2(r, l + nsamp, totmods)] == j)
                {
                    if(m == 0) mods[k] = models_num[index2(r, l, totmods)];
                    temp_index[m] = l;
                    m++;
                }
                l++;
            }
            //calculate correct column
            m = 0;
            for(l = 0; l < sum[j]; l++) m += choose_c(temp_index[l], l + 1);
            m += ncomb_sub[sum[j] - 1];
            inds[k] = m;
            k++;
        }
    }
    //now generate lPDM
    lPDM = 0.0;
    for(j = 0; j < k; j++) lPDM += lPDM_int_mat[index2(mods[j], inds[j], nmod)];
    return lPDM;
}

