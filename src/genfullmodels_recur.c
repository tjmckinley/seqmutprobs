#include "functions.h"
#include <R.h>

//recursive function to be called from 'genfullmodels'
void genfullmodels_recur(int nhier, int currhier1, int *temp, int *tempind, int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int **models_num)
{
    /*'nhier' is number of hierarchies
    'currhier1' is current hierarchy
    'temp' and 'tempind' are vectors recording structural information
    'q' is specific row of 'structure' we wish to extract
    'nrow' and 'ncol' are using for indexing
    'structure' is a matrix containing intermediate model structures
    'r' denotes current model in 'models_num'
    'totmods' is maximum number of models in 'models_num'
    'incmods' is arbitrary number of models to increase if reallocation necessary
    'models_num' is a matrix recording final models*/

    //declare necessary variables
    int j, k, currhier, nsamp, nmod = 10;
    nsamp = ncol;
    currhier = currhier1 + 1;
    //continue working through hierarchies in order and expand model set
    if(temp[currhier] == 1)
    {
        for(j = 0; j < nmod; j++)
        {
            (*models_num)[index2(*r, tempind[currhier], *totmods)] = j;
            if(currhier < (nhier - 1)) genfullmodels_recur(nhier, currhier, temp, tempind, q, nrow, ncol, structure, r, totmods, incmods, models_num);
            else
            {
                *r = *r + 1;
                //perform memory reallocation if required
                if(*r >= *totmods) realloc_maxmod_full(nsamp, totmods, incmods, models_num);
                //initialise new model structure
                for(k = 0; k < (2 * nsamp); k++) (*models_num)[index2(*r, k, *totmods)] = (*models_num)[index2(*r - 1, k, *totmods)];
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
                    (*models_num)[index2(*r, k, *totmods)] = j;
                    (*models_num)[index2(*r, k + nsamp, *totmods)] = tempind[currhier];
                }
            }
            if(currhier < (nhier - 1)) genfullmodels_recur(nhier, currhier, temp, tempind, q, nrow, ncol, structure, r, totmods, incmods, models_num);
            else
            {
                *r = *r + 1;
                //perform memory reallocation if required
                if(*r >= *totmods) realloc_maxmod_full(nsamp, totmods, incmods, models_num);
                //initialise new model structure
                for(k = 0; k < (2 * nsamp); k++) (*models_num)[index2(*r, k, *totmods)] = (*models_num)[index2(*r - 1, k, *totmods)];
            }
        }
    }
    return;
}
