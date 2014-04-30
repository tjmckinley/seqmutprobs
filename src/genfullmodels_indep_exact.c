#include "functions.h"
#include <R.h>

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

    int i, priorind;
    double temp_norm, temp_norm1;

    for(i = 0; i < nmod; i++)
    {
        models_num[currcol] = i;
        if(currcol == (nsamp - 1))
        {
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
        }
        else genfullmodels_indep_exact(currcol + 1, nsamp, nmod, models_num, norm_ls, norm_s, palt_ls, palt_s, lpriornull_ls, lprioralt_ls, lpriornull_s, lprioralt_s, lPDM_int_mat, ncomb_sub, sum, temp_index, mods, inds, prec_ind);
    }
    return;
}
