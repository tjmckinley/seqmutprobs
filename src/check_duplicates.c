#include "functions.h"
#include <R.h>

//internal function to check for duplicated combinations
void check_duplicates(int start, int ncomb, int ncol, int *comb, int *ind)
{
    /*'start' corresponds to which element of 'comb' is the current
    	comparison combination
    'ncomb' is the length of 'comb' (e.g. the number of rows)
    'ncol' is the number of columns of 'comb'
    'comb' is a matrix where each row corresponds to a unique
    	combination of length 'ncol'
    'ind' is binary vector of length 'ncomb' recording unique combinations*/

    int i, j, k;

    for(i = (start + 1); i < ncomb; i++)
    {
        ind[i] = 1;
        for(j = 0; j < ncol; j++) for(k = 0; k < ncol; k++) if(comb[index2(i, j, ncomb)] == comb[index2(start, k, ncomb)]) ind[i] = 0;
    }
    return;
}

