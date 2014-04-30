#include "functions.h"
#include <R.h>

//recursive function to be called from 'genmodels'
void genmodels_recur(int nsamp, int nstart, int *indexl, int *indexl_ind, int repeatind, int *comb_next, int ncomb_next, int ntotsamp, int *currstruct, int *q, int **structure, int *maxmodels, int incmodels)
{
    /*'nsamp' is number of remaining samples to choose from
    'nstart' is number of samples to start selecting
    'indexl' is vector recording which samples are remaining
    'indexl_ind' is a binary indictator denoting which samples
    	can be selected or not
    'repeatind' is binary indicator
    	= 1 if repeated denominator, 0 otherwise
    'comb_next' is matrix recording available samples for special
    	case of repeated denominators
    'ncomb_next' is number of rows of 'comb_next' (used for indexing)
    'ntotsamp' is total number of samples
    'currstruct' is vector recording initial structure at previous hierarchy
    'q' is current model structure indicator
    'structure' is matrix recording models structures
    'maxmodels' and 'incmodels' record the maximum number of models
    	and the number to increment if reallocation necessary*/

    int i, j, k, l, m, nrem, ncomb, ncomb1, nstart1, valid, stop;
    //declare intermediate vectors for recursion
    int *indexl1 = (int *) Calloc(ntotsamp, int);
    int *indexl_ind1 = (int *) Calloc(ntotsamp, int);
    int *currstruct1 = (int *) Calloc(ntotsamp, int);
    int *comb, *comb_set, *comb_ind, *comb_next1;

    if(repeatind == 1) stop = nstart - 1;
    else stop = 0;
    for(k = nstart; k > stop; k--)
    {
        if(k == 1)
        {
            //increment models
            *q = *q + 1;
            //initialise new model structure
            for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
            //perform memory reallocation if required
            if(*q >= *maxmodels) realloc_maxmod(ntotsamp, maxmodels, incmodels, structure);
        }
        else
        {
            /*calculate how many samples remaining, and if necessary
            set up for further loops*/
            nrem = nsamp - k;
            nstart1 = (nrem > k ? k : nrem);
            if(nstart1 > 0)
            {
                if(repeatind == 0)
                {
                    //create set of combinations to use during looping
                    ncomb = choose_c(nsamp, k);
                    comb = (int *) Calloc(k, int);
                    comb_set = (int *) Calloc(k * ncomb, int);
                    comb_ind = (int *) Calloc(ncomb, int);
                    //produce combinations
                    for(j = 0; j < k; j++) comb[j] = j;
                    for(j = 0; j < k; j++) comb_set[index2(0, j, ncomb)] = indexl[comb[j]];
                    i = 1;
                    while(next_comb(comb, k, nsamp) == 1)
                    {
                        for(j = 0; j < k; j++) comb_set[index2(i, j, ncomb)] = indexl[comb[j]];
                        i++;
                    }
                }
                else
                {
                    //create set of combinations to use during looping
                    ncomb = ncomb_next;
                    comb = (int *) Calloc(k, int);
                    comb_set = (int *) Calloc(k * ncomb, int);
                    comb_ind = (int *) Calloc(ncomb, int);
                    //produce combinations
                    for(j = 0; j < ncomb; j++) for(m = 0; m < k; m++) comb_set[index2(j, m, ncomb)] = comb_next[index2(j, m, ncomb)];
                }
                if(nstart1 == k)
                {
                    //now begin looping over combinations
                    i = 0;
                    valid = 1;
                    while(i < (ncomb - 1) && valid == 1)
                    {
                        for(j = 0; j < ncomb; j++) comb_ind[j] = 1;
                        //extract subset that can be looped over further
                        check_duplicates(i, ncomb, k, comb_set, comb_ind);
                        ncomb1 = 0;
                        for(j = (i + 1); j < ncomb; j++) ncomb1 += comb_ind[j];
                        if(ncomb1 == 0) valid = 0;
                        else
                        {
                            //reset indicators
                            for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                            for(j = 0; j < ntotsamp; j++)
                            {
                                indexl1[j] = indexl[j];
                                indexl_ind1[j] = indexl_ind[j];
                            }
                            for(j = 0; j < k; j++) (*structure)[index2(*q, comb_set[index2(i, j, ncomb)], *maxmodels)] = comb_set[index2(i, 0, ncomb)];
                            if(nrem > 0)
                            {
                                if(nrem == 1)
                                {
                                    //increment models
                                    *q = *q + 1;
                                    //initialise new model structure
                                    for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                                    //perform memory reallocation if required
                                    if(*q >= *maxmodels) realloc_maxmod(ntotsamp, maxmodels, incmodels, structure);
                                }
                                else
                                {
                                    //extract viable choices for next level
                                    comb_next1 = (int *) Calloc(ncomb1 * k, int);
                                    l = 0;
                                    for(j = (i + 1); j < ncomb; j++)
                                    {
                                        if(comb_ind[j] == 1)
                                        {
                                            for(m = 0; m < k; m++) comb_next1[index2(l, m, ncomb1)] = comb_set[index2(j, m, ncomb)];
                                            l++;
                                        }
                                    }
                                    //extract remaining samples for loop
                                    for(j = 0; j < k; j++) indexl_ind1[comb_set[index2(i, j, ncomb)]] = 0;
                                    l = 0;
                                    for(j = 0; j < ntotsamp; j++)
                                    {
                                        if(indexl_ind1[j] == 1)
                                        {
                                            indexl1[l] = j;
                                            l++;
                                        }
                                    }
                                    for(j = 0; j < ntotsamp; j++) currstruct1[j] = (*structure)[index2(*q, j, *maxmodels)];
                                    genmodels_recur(nrem, nstart1, indexl1, indexl_ind1, 1, comb_next1, ncomb1, ntotsamp, currstruct1, q, structure, maxmodels, incmodels);
                                    //free memory from heap
                                    Free(comb_next1);
                                    comb_next1 = NULL;
                                }
                            }
                            else
                            {
                                //increment models
                                *q = *q + 1;
                                //initialise new model structure
                                for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                                //perform memory reallocation if required
                                if(*q >= *maxmodels) realloc_maxmod(ntotsamp, maxmodels, incmodels, structure);
                            }
                        }
                        i++;
                    }
                    //if necessary, then loop over further hierarchies
                    nstart1--;
                    if(nstart1 > 1)
                    {
                        while(nstart1 > 1)
                        {
                            //loop over combinations
                            for(i = 0; i < ncomb; i++)
                            {
                                //reset indicators
                                for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                                for(j = 0; j < ntotsamp; j++)
                                {
                                    indexl1[j] = indexl[j];
                                    indexl_ind1[j] = indexl_ind[j];
                                }
                                for(j = 0; j < k; j++) (*structure)[index2(*q, comb_set[index2(i, j, ncomb)], *maxmodels)] = comb_set[index2(i, 0, ncomb)];
                                if(nrem > 0)
                                {
                                    //extract remaining samples for loop
                                    for(j = 0; j < k; j++) indexl_ind1[comb_set[index2(i, j, ncomb)]] = 0;
                                    l = 0;
                                    for(j = 0; j < ntotsamp; j++)
                                    {
                                        if(indexl_ind1[j] == 1)
                                        {
                                            indexl1[l] = j;
                                            l++;
                                        }
                                    }
                                    if(nrem == 1)
                                    {
                                        //increment models
                                        *q = *q + 1;
                                        //initialise new model structure
                                        for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                                        //perform memory reallocation if required
                                        if(*q >= *maxmodels) realloc_maxmod(ntotsamp, maxmodels, incmodels, structure);
                                    }
                                    else
                                    {
                                        for(j = 0; j < ntotsamp; j++) currstruct1[j] = (*structure)[index2(*q, j, *maxmodels)];
                                        genmodels_recur(nrem, nstart1, indexl1, indexl_ind1, 0, comb_next1, 0, ntotsamp, currstruct1, q, structure, maxmodels, incmodels);
                                    }
                                }
                                else
                                {
                                    //increment models
                                    *q = *q + 1;
                                    //initialise new model structure
                                    for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                                    //perform memory reallocation if required
                                    if(*q >= *maxmodels) realloc_maxmod(ntotsamp, maxmodels, incmodels, structure);
                                }
                            }
                            nstart1--;
                        }
                    }
                    else
                    {
                        //now repeat for non-repeated hierarchies - looping over combinations
                        for(i = 0; i < ncomb; i++)
                        {
                            for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                            for(j = 0; j < k; j++) (*structure)[index2(*q, comb_set[index2(i, j, ncomb)], *maxmodels)] = comb_set[index2(i, 0, ncomb)];
                            //increment models
                            *q = *q + 1;
                            //initialise new model structure
                            for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                            //perform memory reallocation if required
                            if(*q >= *maxmodels) realloc_maxmod(ntotsamp, maxmodels, incmodels, structure);
                        }
                    }
                }
                else
                {
                    //loop over combinations
                    for(i = 0; i < ncomb; i++)
                    {
                        //reset indicators
                        for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                        for(j = 0; j < ntotsamp; j++)
                        {
                            indexl1[j] = indexl[j];
                            indexl_ind1[j] = indexl_ind[j];
                        }
                        for(j = 0; j < k; j++) (*structure)[index2(*q, comb_set[index2(i, j, ncomb)], *maxmodels)] = comb_set[index2(i, 0, ncomb)];
                        if(nrem > 0)
                        {
                            //extract remaining samples for loop
                            for(j = 0; j < k; j++) indexl_ind1[comb_set[index2(i, j, ncomb)]] = 0;
                            l = 0;
                            for(j = 0; j < ntotsamp; j++)
                            {
                                if(indexl_ind1[j] == 1)
                                {
                                    indexl1[l] = j;
                                    l++;
                                }
                            }
                            if(nrem == 1)
                            {
                                //increment models
                                *q = *q + 1;
                                //initialise new model structure
                                for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                                //perform memory reallocation if required
                                if(*q >= *maxmodels) realloc_maxmod(ntotsamp, maxmodels, incmodels, structure);
                            }
                            else
                            {
                                for(j = 0; j < ntotsamp; j++) currstruct1[j] = (*structure)[index2(*q, j, *maxmodels)];
                                genmodels_recur(nrem, nstart1, indexl1, indexl_ind1, 0, comb_next1, 0, ntotsamp, currstruct1, q, structure, maxmodels, incmodels);
                            }
                        }
                        else
                        {
                            //increment models
                            *q = *q + 1;
                            //initialise new model structure
                            for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                            //perform memory reallocation if required
                            if(*q >= *maxmodels) realloc_maxmod(ntotsamp, maxmodels, incmodels, structure);
                        }
                    }
                }
                //free memory from the heap (automatically sets pointer to NULL)
                Free(comb);
                Free(comb_set);
                Free(comb_ind);
            }
            else
            {
                for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                for(j = 0; j < k; j++) (*structure)[index2(*q, indexl[j], *maxmodels)] = indexl[0];
                //increment models
                *q = *q + 1;
                //initialise new model structure
                for(j = 0; j < ntotsamp; j++) (*structure)[index2(*q, j, *maxmodels)] = currstruct[j];
                //perform memory reallocation if required
                if(*q >= *maxmodels) realloc_maxmod(ntotsamp, maxmodels, incmodels, structure);
            }
        }
    }
    //free memory from heap
    Free(indexl1);
    Free(indexl_ind1);
    Free(currstruct1);
    return;
}

