#include "functions.h"
#include <R.h>

//function to generate all possible models based on an intermediate set of structures
void genfullmodels(int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int **models_num)
{
    /*'q' is specific row of 'structure' we wish to extract
    'nrow' and 'ncol' are using for indexing
    'structure' is a matrix containing intermediate model structures
    'r' denotes current model in 'models_num'
    'totmods' is maximum number of models in 'models_num'
    'incmods' is arbitrary number of models to increase if reallocation necessary
    'models_num' is a matrix recording final models*/

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
                (*models_num)[index2(*r, j, *totmods)] = i;
                (*models_num)[index2(*r, j + nsamp, *totmods)] = 0;
            }
            *r = *r + 1;
            //perform memory reallocation if required
            if(*r >= *totmods) realloc_maxmod_full(nsamp, totmods, incmods, models_num);
        }
    }
    else
    {
        //temp[0] must be >1 in this loop
        for(j = 0; j < nmod; j++)
        {
            //initialise structures
            for(k = 0; k < nsamp; k++) (*models_num)[index2(*r, k + nsamp, *totmods)] = k;
            //now fill in gaps and pass to recursive function
            for(k = 0; k < nsamp; k++)
            {
                if(structure[index2(q, k, nrow)] == tempind[0])
                {
                    (*models_num)[index2(*r, k, *totmods)] = j;
                    (*models_num)[index2(*r, k + nsamp, *totmods)] = tempind[0];
                }
            }
            genfullmodels_recur(nhier, 0, temp, tempind, q, nrow, ncol, structure, r, totmods, incmods, models_num);
        }
    }
    //free memory from the heap (automatically sets pointer to NULL)
    Free(temp);
    Free(tempind);
    return;
}

