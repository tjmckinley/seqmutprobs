#include "functions.h"
#include <R.h>

//recursive function to produce lexicographic orderings
void lexico_recur(int nsamp, int col, int totcol, int *temp_samp, int *temp_order1, int norder, int start, int *overall_order, int *ncomb_sub)
{
    /*'nsamp' is the number of samples
    'col' is current column
    'totcol' is total number of columns
    'temp_samp' is a vector recording where to start ordering in absolute terms
    'temp_order1' is a matrix (in vector form) containing current orderings
    'norder' is number of rows of 'temp_order1'
    'start' is the absolute starting point in 'temp_order1' used for indexing in recursion
    'overall_order' is vector containing absolute orderings across all possible dependencies
    'ncomb_sub' is intermediate vector used for indexing 'overall_order'*/

    int q, l, k, temp_norder;
    int * temp_ind = (int *) Calloc(norder, int);
    int * temp_val = (int *) Calloc(norder, int);
    int * temp_samp1 = (int *) Calloc(nsamp, int);

    for(q = 1; q < nsamp; q++)
    {
        temp_norder = temp_samp[q] - temp_samp[q - 1];
        if(temp_norder > 1)
        {
            bubble_sort_inc_int(&temp_order1[index2(start + temp_samp[q - 1], col, norder)], temp_ind, temp_norder);
            /*			Rprintf("totcol=%d\tcol=%d\tq=%d\n",totcol,col,q);*/
            /*			for(l=0;l<temp_norder;l++) Rprintf("%d\n",temp_ind[l]);*/
            /*			Rprintf("\n");*/
            for(l = (col + 1); l < totcol; l++)
            {
                for(k = 0; k < temp_norder; k++) temp_val[k] = temp_order1[index2(start + temp_samp[q - 1] + k, l, norder)];
                for(k = 0; k < temp_norder; k++) temp_order1[index2(start + temp_samp[q - 1] + k, l, norder)] = temp_val[temp_ind[k]];
            }
            //adjust overall ordering
            for(k = 0; k < temp_norder; k++) temp_val[k] = overall_order[ncomb_sub[totcol - 1] + start + temp_samp[q - 1] + k];
            for(k = 0; k < temp_norder; k++) overall_order[ncomb_sub[totcol - 1] + start + temp_samp[q - 1] + k] = temp_val[temp_ind[k]];
            //now generate conditions for sequential sorting at next level
            for(k = 0; k < nsamp; k++)
            {
                temp_samp1[k] = 0;
                for(l = (start + temp_samp[q - 1]); l < (start + temp_samp[q - 1] + temp_norder); l++) temp_samp1[k] += (temp_order1[index2(l, col, norder)] == k ? 1 : 0);
            }
            //convert to cumulative counts
            for(l = 1; l < nsamp; l++) temp_samp1[l] += temp_samp1[l - 1];
            //enter into further recursion
            if((col + 1) < totcol) lexico_recur(nsamp, col + 1, totcol, temp_samp1, temp_order1, norder, start + temp_samp[q - 1], overall_order, ncomb_sub);
        }
    }
    Free(temp_ind);
    Free(temp_val);
    Free(temp_samp1);
    return;
}

