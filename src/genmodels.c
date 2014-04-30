#include "functions.h"
#include <R.h>
#include <Rinternals.h>

/*function to generate all model structures in an efficient manner
(returns an R vector of a given size, though it is not necessary to
specify the size in advance)*/
SEXP genmodels(SEXP nsamp1)
{
    /*'nsamp1' is the number of samples*/

    //declare variables
    int i, j, k, l, m, q, r, nsamp, nrem, nstart, ncomb, ncomb1, valid;
    //extract integer parts of 'nsamp1'
    nsamp1 = coerceVector(nsamp1, INTSXP);
    nsamp = INTEGER(nsamp1)[0];
    //pointer for memory reallocation
    int * temp_ptr = NULL;

    /*this is arbitrary number of initial models
    (can be altered if necessary during runtime)*/
    int maxmodels, incmodels = 10000;
    maxmodels = incmodels;
    /*set up pointer to array of pointers, each of which will point at an
    array holding relevant model structures (to be filled in during runtime)*/
    int *structure = (int *) Calloc(maxmodels * nsamp, int);
    //declare intermediate vectors for recursion
    int *indexl = (int *) Calloc(nsamp, int);
    int *indexl_ind = (int *) Calloc(nsamp, int);
    int *currstruct = (int *) Calloc(nsamp, int);
    int *comb, *comb_set, *comb_ind, *comb_next;

    //for fully DEPENDENT model structures
    q = 0;
    for(j = 0; j < nsamp; j++) structure[index2(q, j, maxmodels)] = 0;
    q++;

    //now complete for PARTIALLY DEPENDENT model structures
    for(k = (nsamp - 1); k > 1; k--)
    {
        //create set of combinations to use during looping
        ncomb = choose_c(nsamp, k);
        comb = (int *) Calloc(k, int);
        comb_set = (int *) Calloc(k * ncomb, int);
        comb_ind = (int *) Calloc(ncomb, int);
        //produce combinations
        for(j = 0; j < k; j++) comb[j] = j;
        for(j = 0; j < k; j++) comb_set[index2(0, j, ncomb)] = comb[j];
        i = 1;
        while(next_comb(comb, k, nsamp) == 1)
        {
            for(j = 0; j < k; j++) comb_set[index2(i, j, ncomb)] = comb[j];
            i++;
        }
        /*calculate how many samples remaining, and if necessary
        set up for further loops*/
        nrem = nsamp - k;
        nstart = (nrem > k ? k : nrem);
        if(nstart > 1)
        {
            if(nstart == k)
            {
                //now begin looping over combinations
                i = 0;
                valid = 1;
                while(i < (ncomb - 1) && valid == 1)
                {
                    //extract subset that can be looped over further
                    check_duplicates(i, ncomb, k, comb_set, comb_ind);
                    ncomb1 = 0;
                    for(j = (i + 1); j < ncomb; j++) ncomb1 += comb_ind[j];
                    if(ncomb1 == 0) valid = 0;
                    else
                    {
                        //reset indicators
                        for(j = 0; j < nsamp; j++) indexl_ind[j] = 1;
                        //initialise structures
                        for(j = 0; j < nsamp; j++) structure[index2(q, j, maxmodels)] = j;
                        for(j = 0; j < k; j++) structure[index2(q, comb_set[index2(i, j, ncomb)], maxmodels)] = comb_set[index2(i, 0, ncomb)];
                        //extract viable choices for next level
                        comb_next = (int *) Calloc(ncomb1 * k, int);
                        l = 0;
                        for(j = (i + 1); j < ncomb; j++)
                        {
                            if(comb_ind[j] == 1)
                            {
                                for(m = 0; m < k; m++) comb_next[index2(l, m, ncomb1)] = comb_set[index2(j, m, ncomb)];
                                l++;
                            }
                        }
                        //extract remaining samples for loop
                        for(j = 0; j < k; j++) indexl_ind[comb_set[index2(i, j, ncomb)]] = 0;
                        l = 0;
                        for(j = 0; j < nsamp; j++)
                        {
                            if(indexl_ind[j] == 1)
                            {
                                indexl[l] = j;
                                l++;
                            }
                        }
                        for(j = 0; j < nsamp; j++) currstruct[j] = structure[index2(q, j, maxmodels)];
                        genmodels_recur(nrem, nstart, indexl, indexl_ind, 1, comb_next, ncomb1, nsamp, currstruct, &q, &structure, &maxmodels, incmodels);
                        //free memory from heap
                        Free(comb_next);
                        comb_next = NULL;
                    }
                    i++;
                }
                //now enter further recursive loop if required
                nstart--;
                if(nstart > 1)
                {
                    while(nstart > 1)
                    {
                        //loop over combinations
                        for(i = 0; i < ncomb; i++)
                        {
                            //reset indicators
                            for(j = 0; j < nsamp; j++) indexl_ind[j] = 1;
                            //initialise structures
                            for(j = 0; j < nsamp; j++) structure[index2(q, j, maxmodels)] = j;
                            for(j = 0; j < k; j++) structure[index2(q, comb_set[index2(i, j, ncomb)], maxmodels)] = comb_set[index2(i, 0, ncomb)];
                            //extract remaining samples for loop
                            for(j = 0; j < k; j++) indexl_ind[comb_set[index2(i, j, ncomb)]] = 0;
                            l = 0;
                            for(j = 0; j < nsamp; j++)
                            {
                                if(indexl_ind[j] == 1)
                                {
                                    indexl[l] = j;
                                    l++;
                                }
                            }
                            for(j = 0; j < nsamp; j++) currstruct[j] = structure[index2(q, j, maxmodels)];
                            genmodels_recur(nrem, nstart, indexl, indexl_ind, 0, comb_next, 0, nsamp, currstruct, &q, &structure, &maxmodels, incmodels);
                        }
                        nstart--;
                    }
                }
                else
                {
                    //now repeat for non-repeated hierarchies - looping over combinations
                    for(l = 0; l < ncomb; l++)
                    {
                        for(j = 0; j < nsamp; j++) structure[index2(q, j, maxmodels)] = j;
                        for(j = 0; j < k; j++) structure[index2(q, comb_set[index2(l, j, ncomb)], maxmodels)] = comb_set[index2(l, 0, ncomb)];
                        q++;
                        //perform memory reallocation if required
                        if(q >= maxmodels) realloc_maxmod(nsamp, &maxmodels, incmodels, &structure);
                    }
                }
            }
            else
            {
                //loop over combinations
                for(i = 0; i < ncomb; i++)
                {
                    //reset indicators
                    for(j = 0; j < nsamp; j++) indexl_ind[j] = 1;
                    //initialise structures
                    for(j = 0; j < nsamp; j++) structure[index2(q, j, maxmodels)] = j;
                    for(j = 0; j < k; j++) structure[index2(q, comb_set[index2(i, j, ncomb)], maxmodels)] = comb_set[index2(i, 0, ncomb)];
                    //extract remaining samples for loop
                    for(j = 0; j < k; j++) indexl_ind[comb_set[index2(i, j, ncomb)]] = 0;
                    l = 0;
                    for(j = 0; j < nsamp; j++)
                    {
                        if(indexl_ind[j] == 1)
                        {
                            indexl[l] = j;
                            l++;
                        }
                    }
                    for(j = 0; j < nsamp; j++) currstruct[j] = structure[index2(q, j, maxmodels)];
                    genmodels_recur(nrem, nstart, indexl, indexl_ind, 0, comb_next, 0, nsamp, currstruct, &q, &structure, &maxmodels, incmodels);
                }
            }
        }
        else
        {
            //complete model structures
            for(i = 0; i < ncomb; i++)
            {
                for(j = 0; j < nsamp; j++) structure[index2(q, j, maxmodels)] = j;
                for(j = 0; j < k; j++) structure[index2(q, comb_set[index2(i, j, ncomb)], maxmodels)] = comb_set[index2(i, 0, ncomb)];
                q++;
                //perform memory reallocation if required
                if(q >= maxmodels) realloc_maxmod(nsamp, &maxmodels, incmodels, &structure);
            }
        }
        //free memory from the heap (automatically sets pointer to NULL)
        Free(comb);
        Free(comb_ind);
        Free(comb_set);
    }
    //now reallocate memory to save space
    temp_ptr = (int *) Calloc(q * nsamp, int);
    for(i = 0; i < nsamp; i++) for(j = 0; j < q; j++) temp_ptr[index2(j, i, q)] = structure[index2(j, i, maxmodels)];
    Free(structure);
    structure = temp_ptr;
    temp_ptr = NULL;
    //now generate all possible model structures
    int totmods, incmods = 100000;
    totmods = incmods;
    int * models_num = (int *) Calloc(totmods * 2 * nsamp, int);
    r = 0;
    for(i = 0; i < q; i++) genfullmodels(i, q, nsamp, structure, &r, &totmods, incmods, &models_num);

    //now append FULLY INDEPENDENT models
    int nmod = 10;
    l = r + pow_int(nmod, nsamp);
    //now reallocate memory to save space
    temp_ptr = (int *) Calloc(l * 2 * nsamp, int);
    for(i = 0; i < (2 * nsamp); i++) for(j = 0; j < r; j++) temp_ptr[index2(j, i, l)] = models_num[index2(j, i, totmods)];
    Free(models_num);
    models_num = temp_ptr;
    temp_ptr = NULL;
    totmods = l;
    //initialise indicator part of models_num
    for(i = r; i < totmods; i++) for(j = 0; j < nsamp; j++) models_num[index2(i, j + nsamp, totmods)] = j;
    //expand grid of independent models into models_num
    int each, rep;
    for(i = 0; i < nsamp; i++)
    {
        l = r;
        each = pow_int(nmod, i);
        rep = pow_int(nmod, nsamp - i - 1);
        while(rep > 0)
        {
            for(k = 0; k < nmod; k++)
            {
                for(j = 0; j < each; j++)
                {
                    models_num[index2(l, i, totmods)] = k;
                    l++;
                }
            }
            rep--;
        }
    }
    //now export into R and return
    SEXP models_num_R;
    PROTECT(models_num_R = allocVector(INTSXP, totmods * 2 * nsamp));
    for(i = 0; i < totmods; i++) for(j = 0; j < (2 * nsamp); j++) INTEGER(models_num_R)[index2(i, j, totmods)] = models_num[index2(i, j, totmods)];
    UNPROTECT(1);
    //free memory from the heap (automatically sets pointer to NULL)
    Free(structure);
    Free(models_num);
    Free(indexl);
    Free(indexl_ind);
    Free(currstruct);
    return models_num_R;
}

