#include "functions.h"
#include <R.h>
#include <Rinternals.h>

/*function to generate all model structures in an efficient manner
and calculate normalising for use in EXACT efficient search routine*/
SEXP genmodels_exact(SEXP nsamp1, SEXP ntotcol1, SEXP lpriors, SEXP lPDM_int_mat1)
{
    /*'nsamp1' is the number of samples
    'ntotcol1' is number of columns of 'lPDM_int_mat1'
    'lpriors' is a vector of log-priors in order c(null_ls,alt_ls,null_s,alt_s)
    'lPDM_int_mat1' is a matrix for intermediate values for use in
    	calculating normalising constant*/

    //declare variables
    int i, j, k, l, m, q, r, nsamp, nrem, nstart, ncomb, ncomb1, valid, nmod = 10, ntotcol;
    double lpriornull_ls, lpriornull_s, lprioralt_ls, lprioralt_s;
    //extract integer parts of 'nsamp1'
    nsamp1 = coerceVector(nsamp1, INTSXP);
    nsamp = INTEGER(nsamp1)[0];
    //extract integer parts of 'ntotcol1'
    ntotcol1 = coerceVector(ntotcol1, INTSXP);
    ntotcol = INTEGER(ntotcol1)[0];
    //extract log-prior information
    lpriors = coerceVector(lpriors, REALSXP);
    lpriornull_ls = REAL(lpriors)[0];
    lprioralt_ls = REAL(lpriors)[1];
    lpriornull_s = REAL(lpriors)[2];
    lprioralt_s = REAL(lpriors)[3];
    //extract intermediate vector
    double *lPDM_int_mat = (double *) R_alloc(nmod * ntotcol, sizeof(double));
    lPDM_int_mat1 = coerceVector(lPDM_int_mat1, REALSXP);
    for(i = 0; i < (nmod * ntotcol); i++) lPDM_int_mat[i] = REAL(lPDM_int_mat1)[i];

    //generate intermediate vectors for calculations
    int *ncomb_sub = (int *) R_alloc(nsamp, sizeof(int));
    ncomb_sub[0] = 0;
    for(i = 1; i < nsamp; i++) ncomb_sub[i] = choose_c(nsamp, i);
    for(i = 1; i < nsamp; i++) ncomb_sub[i] = ncomb_sub[i] + ncomb_sub[i - 1];
    //now declare intermediate vectors for model simplification
    int *sum = (int *) R_alloc(nsamp, sizeof(int));
    int *temp_index = (int *) R_alloc(nsamp, sizeof(int));
    int *mods = (int *) R_alloc(nsamp, sizeof(int));
    int *inds = (int *) R_alloc(nsamp, sizeof(int));

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
    int * temp_ptr = (int *) Calloc(q * nsamp, int);
    for(i = 0; i < nsamp; i++) for(j = 0; j < q; j++) temp_ptr[index2(j, i, q)] = structure[index2(j, i, maxmodels)];
    Free(structure);
    structure = temp_ptr;
    temp_ptr = NULL;
    /*now generate all possible model structures in turn, and count
    the total number of models and the number of alternative models for
    each hypothesis in turn and figure out partial normalising constant*/
    double norm_ls = 0.0, norm_s = 0.0, palt_ls = 0.0, palt_s = 0.0;
    int * models_num = (int *) Calloc(2 * nsamp, int);
    int prec_ind = 0;
    r = 0;
    for(i = 0; i < q; i++) genfullmodels_exact(i, q, nsamp, structure, &r, models_num, &norm_ls, &norm_s, &palt_ls, &palt_s, lpriornull_ls, lprioralt_ls, lpriornull_s, lprioralt_s, lPDM_int_mat, ncomb_sub, sum, temp_index, mods, inds, &prec_ind);

    //now append FULLY INDEPENDENT models
    l = r + pow_int(nmod, nsamp);
    //initialise indicator part of models_num
    for(j = 0; j < nsamp; j++) models_num[j + nsamp] = j;
    //use recursive loop
    for(i = 0; i < nmod; i++)
    {
        models_num[0] = i;
        genfullmodels_indep_exact(1, nsamp, nmod, models_num, &norm_ls, &norm_s, &palt_ls, &palt_s, lpriornull_ls, lprioralt_ls, lpriornull_s, lprioralt_s, lPDM_int_mat, ncomb_sub, sum, temp_index, mods, inds, &prec_ind);
    }
    //now export into R and return
    SEXP output_R;
    PROTECT(output_R = allocVector(REALSXP, 5));
    REAL(output_R)[0] = norm_ls;
    REAL(output_R)[1] = norm_s;
    REAL(output_R)[2] = palt_ls / norm_ls;
    REAL(output_R)[3] = palt_s / norm_s;
    REAL(output_R)[4] = prec_ind;
    UNPROTECT(1);
    //free memory from the heap (automatically sets pointers to NULL)
    Free(structure);
    Free(models_num);
    Free(indexl);
    Free(indexl_ind);
    Free(currstruct);
    return output_R;
}

