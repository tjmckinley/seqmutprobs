#include "functions.h"
#include <R.h>

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

    int l, m;

    //now generate correct column indicators
    l = 0;
    m = 0;
    while(l < nsamp && m < sum)
    {
        if(structure[index2(q, l, nrow)] == samp)
        {
            temp_index[m] = l;
            m++;
        }
        l++;
    }
    //calculate correct column
    m = 0;
    for(l = 0; l < sum; l++) m += choose_c(temp_index[l], l + 1);
    m += ncomb_sub[sum - 1];
    return m;
}
