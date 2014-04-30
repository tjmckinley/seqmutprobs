/*C code for generating models for sequence data, intended
to be called by associated R code*/

#include <R.h>
#include <Rinternals.h>
#include "functions.h"

//function to calculate "top" set of models for use in approximation
SEXP calc_PPAs_approx_fn(SEXP nsamp1, SEXP ntotcol1, SEXP logc1, SEXP lpriornull1, SEXP lprioralt1, SEXP lPDM_int_mat1, SEXP string1, SEXP structure1, SEXP nstruct1, SEXP return_final1, SEXP uni1, SEXP uni_ind1)
{
    /*'nsamp1' is the number of samples
    'ntotcol1' is total number of columns in 'lPDM_int_mat'
    'logc1' is the log-threshold value for acceptance relative to the maximum model
    'lpriornull1' is the log[P(M)] for a single 'null' model
    'lprioralt1' is the log[P(M)] for a single 'alternative' model
    'lPDM_int_mat1' is a matrix of log[P'(D|M)] values for use in intermediary calculations
    	for calculating independent model structures
    'string1' is a binary indicator denoting prior criterion, such that
    	0=less stringent and 1=stringent
    'structure1' is a matrix (in vector form) containing structures generated from
    	'genmodels_priors'
    'nstruct1' is number of structures in 'structure'
    'return_final' is a binary value taking the value 1 if model sets are to be returned and 0 otherwise
    'uni1' and 'uni_ind1' are matrices (in vector form) for use in producing PPAs without recording the full model sets*/

    /*DECLARE VARIABLES*/
    int i, j, k, l, m, q, d, r, nmod = 10, nsamp, incmods, totmods, ntotcol, nmodcol, string, nstruct, currsamp, valid, return_final;
    double lpriornull, lprioralt, lpriordiff, temp_lPPA, max_lPPA, logc, thresh;
    //extract relevant parts of inputs for use in C code
    nsamp1 = coerceVector(nsamp1, INTSXP);
    nsamp = INTEGER(nsamp1)[0];
    ntotcol1 = coerceVector(ntotcol1, INTSXP);
    ntotcol = INTEGER(ntotcol1)[0];
    string1 = coerceVector(string1, INTSXP);
    string = INTEGER(string1)[0];
    nstruct1 = coerceVector(nstruct1, INTSXP);
    nstruct = INTEGER(nstruct1)[0];
    logc1 = coerceVector(logc1, REALSXP);
    logc = REAL(logc1)[0];
    lpriornull1 = coerceVector(lpriornull1, REALSXP);
    lpriornull = REAL(lpriornull1)[0];
    lprioralt1 = coerceVector(lprioralt1, REALSXP);
    lprioralt = REAL(lprioralt1)[0];
    return_final1 = coerceVector(return_final1, INTSXP);
    return_final = INTEGER(return_final1)[0];
    nmodcol = nsamp * 2;
    incmods = 10000;
    totmods = incmods;
    //DECLARE AND INITIALISE FIXED-SIZE ARRAYS
    int *structure = (int *) R_alloc(nsamp * nstruct, sizeof(int));
    structure1 = coerceVector(structure1, INTSXP);
    for(i = 0; i < nstruct; i++) for(j = 0; j < nsamp; j++) structure[index2(i, j, nstruct)] = INTEGER(structure1)[index2(i, j, nstruct)];
    double *lPDM_int_mat = (double *) R_alloc(nmod * ntotcol, sizeof(double));
    lPDM_int_mat1 = coerceVector(lPDM_int_mat1, REALSXP);
    for(i = 0; i < (nmod * ntotcol); i++) lPDM_int_mat[i] = REAL(lPDM_int_mat1)[i];
    int *uni = (int *) R_alloc(nsamp * nmod * nmod, sizeof(int));
    uni1 = coerceVector(uni1, INTSXP);
    for(i = 0; i < (nsamp * nmod * nmod); i++) uni[i] = INTEGER(uni1)[i];
    int *uni_ind = (int *) R_alloc(nsamp * nmod, sizeof(int));
    uni_ind1 = coerceVector(uni_ind1, INTSXP);
    for(i = 0; i < (nsamp * nmod); i++) uni_ind[i] = INTEGER(uni_ind1)[i];
    //print to screen as a check
    /*	Rprintf("uni\n");*/
    /*	for(i=0;i<nmod;i++)*/
    /*	{*/
    /*		for(j=0;j<nsamp;j++) for(k=0;k<nmod;k++) Rprintf("%d\t",uni[index3(i,k,j,nmod,nmod)]);*/
    /*		Rprintf("\n");*/
    /*	}*/
    /*	Rprintf("uni_ind\n");*/
    /*	for(i=0;i<nmod;i++)*/
    /*	{*/
    /*		for(j=0;j<nsamp;j++) Rprintf("%d\t",uni_ind[index2(i,j,nmod)]);*/
    /*		Rprintf("\n");*/
    /*	}*/

    //declare vectors for intermediate calculations
    int *lPDM_int_ind = (int *) R_alloc(nmod * ntotcol, sizeof(int));
    for(i = 0; i < (nmod * ntotcol); i++) lPDM_int_ind[i] = 0;
    int *currmods = (int *) R_alloc(nsamp, sizeof(int));
    /*	int *currmods1 = (int *) R_alloc(nsamp,sizeof(int));*/
    int *currlocs = (int *) R_alloc(nsamp, sizeof(int));
    /*	int *currlocs1 = (int *) R_alloc(nsamp,sizeof(int));*/
    int *maxmods = (int *) R_alloc(nsamp, sizeof(int));
    int *maxlocs = (int *) R_alloc(nsamp, sizeof(int));
    /*	int *comb = (int *) R_alloc(nsamp,sizeof(int));*/
    //now generate intermediate vector for model simplification
    /*	int *temp_index = (int *) R_alloc(nsamp,sizeof(int));*/
    /*generate vector containing number of combinations in each subset
    (for generating column indices for use in model simplification)*/
    int *ncomb_sub = (int *) R_alloc(nsamp + 1, sizeof(int));
    ncomb_sub[0] = 0;
    for(i = 1; i <= nsamp; i++) ncomb_sub[i] = choose_c(nsamp, i);
    for(i = 1; i <= nsamp; i++) ncomb_sub[i] = ncomb_sub[i] + ncomb_sub[i - 1];

    //DECLARE AND INITIALISE VARIABLE-SIZED ARRAYS
    //generate model output vector
    double *lPPA_mat = (double *) Calloc(totmods, double);
    int *hyp = (int *) Calloc(totmods, int);
    int *models_num = (int *) Calloc(nmodcol * totmods, int);

    //order lPDM matrix and store correctly ordered indices
    for(i = 0; i < nmod; i++) for(j = 0; j < ntotcol; j++) lPDM_int_ind[index2(i, j, nmod)] = i;

    for(i = 0; i < ntotcol; i++) bubble_sort_dec(&lPDM_int_mat[index2(0, i, nmod)], &lPDM_int_ind[index2(0, i, nmod)], nmod);
    //record current maximum model
    r = 0;
    max_lPPA = 0.0;
    for(j = 0; j < nsamp; j++)
    {
        max_lPPA += lPDM_int_mat[index2(0, j, nmod)];
        models_num[index2_col(r, j, nmodcol)] = lPDM_int_ind[index2(0, j, nmod)];
        models_num[index2_col(r, j + nsamp, nmodcol)] = j;
    }
    lPPA_mat[r] = max_lPPA;
    //set up vectors to record current maximum model plus locations
    for(j = 0; j < nsamp; j++)
    {
        maxmods[j] = models_num[index2_col(r, j, nmodcol)];
        maxlocs[j] = 0;
    }
    //update model counter
    r++;

    //now calculate maximum effect of changing prior
    lpriordiff = lpriornull - lprioralt;
    if(lpriordiff < 0) lpriordiff = -lpriordiff;
    //update threshold
    thresh = logc + lpriordiff;

    //run recursion through to calculate maximum fully independent model
    i = 0;
    valid = 1;
    currsamp = 0;
    while(i < 10 && valid == 1)
    {
        //update model
        for(j = 0; j < nsamp; j++) models_num[index2_col(r, j, nmodcol)] = maxmods[j];
        models_num[index2_col(r, currsamp, nmodcol)] = lPDM_int_ind[index2(i, currsamp, nmod)];
        if(models_num[index2_col(r, currsamp, nmodcol)] == maxmods[currsamp])
        {
            //enter into recursion
            for(j = 0; j < nsamp; j++)
            {
                currmods[j] = maxmods[j];
                currlocs[j] = maxlocs[j];
            }
            currlocs[currsamp] = i;
            if((currsamp + 1) < nsamp) calc_approx_indep_recur(nsamp, nmod, nmodcol, incmods, &totmods, currsamp + 1, &r, currmods, currlocs, maxmods, maxlocs, &max_lPPA, max_lPPA, thresh, &models_num, &hyp, lPDM_int_mat, lPDM_int_ind, &lPPA_mat);
            i++;
        }
        else
        {
            //calculate new log-PDM
            temp_lPPA = max_lPPA - lPDM_int_mat[index2(maxlocs[currsamp], currsamp, nmod)] + lPDM_int_mat[index2(i, currsamp, nmod)];
            //add to lPPA_mat output vector
            lPPA_mat[r] = temp_lPPA;
            //if new model is maximum, then update maximum model
            if(temp_lPPA <= max_lPPA)
            {
                //if new model is accepted then update model set
                if((max_lPPA - temp_lPPA) < thresh)
                {
                    //increase model counter
                    r++;
                    //increase size of output vectors if necessary
                    if(r > (totmods - 1)) realloc_approx(nmodcol, &totmods, incmods, &models_num, &hyp, &lPPA_mat);
                    //enter into recursion
                    for(j = 0; j < nsamp; j++)
                    {
                        currmods[j] = models_num[index2_col(r - 1, j, nmodcol)];
                        currlocs[j] = maxlocs[j];
                    }
                    currlocs[currsamp] = i;
                    if((currsamp + 1) < nsamp) calc_approx_indep_recur(nsamp, nmod, nmodcol, incmods, &totmods, currsamp + 1, &r, currmods, currlocs, maxmods, maxlocs, &max_lPPA, temp_lPPA, thresh, &models_num, &hyp, lPDM_int_mat, lPDM_int_ind, &lPPA_mat);

                }
                else valid = 0;
            }
            else
            {
                //increase model counter
                r++;
                //increase size of output vectors if necessary
                if(r > (totmods - 1)) realloc_approx(nmodcol, &totmods, incmods, &models_num, &hyp, &lPPA_mat);
                //update maximum model
                max_lPPA = temp_lPPA;
                //set up vectors to record current maximum model plus locations
                maxmods[currsamp] = models_num[index2_col(r - 1, currsamp, nmodcol)];
                maxlocs[currsamp] = i;
                //enter into recursion
                for(j = 0; j < nsamp; j++)
                {
                    currmods[j] = models_num[index2_col(r - 1, j, nmodcol)];;
                    currlocs[j] = maxlocs[j];
                }
                if((currsamp + 1) < nsamp) calc_approx_indep_recur(nsamp, nmod, nmodcol, incmods, &totmods, currsamp + 1, &r, currmods, currlocs, maxmods, maxlocs, &max_lPPA, temp_lPPA, thresh, &models_num, &hyp, lPDM_int_mat, lPDM_int_ind, &lPPA_mat);
            }
            i++;
        }
    }

    //sort all independent models into rank order according to log(PPA')
    int *temp_ind = (int *) Calloc(r, int);
    bubble_sort_dec(lPPA_mat, temp_ind, r);
    //generate sorted vectors/matrices
    int * new_models_num = (int *) Calloc(nmodcol * totmods, int);
    for(i = 0; i < r; i++) for(j = 0; j < nmodcol; j++) new_models_num[index2_col(i, j, nmodcol)] = models_num[index2_col(temp_ind[i], j, nmodcol)];
    Free(models_num);
    models_num = new_models_num;
    new_models_num = NULL;
    Free(temp_ind);
    temp_ind = NULL;


    /*	for(i=0;i<r;i++)*/
    /*	{*/
    /*		for(j=0;j<nmodcol;j++) Rprintf("%d\t",models_num[index2_col(i,j,nmodcol)]);*/
    /*		Rprintf("%f\n",lPPA_mat[i]);*/
    /*	}*/
    /*	Rprintf("return_final=%d\n",return_final);*/
    /*	getchar();*/

    //fill in dependency structures
    for(i = 0; i < r; i++) for(j = 0; j < nsamp; j++) models_num[index2_col(i, j + nsamp, nmodcol)] = j;

    //now generate all viable dependent structures relative to maximum INDEPENDENT model
    for(i = 0; i < nstruct; i++) genfullmodels_approx(i, nstruct, nsamp, structure, &r, &totmods, incmods, nmodcol, lPDM_int_mat, lPDM_int_ind, ntotcol, &models_num, &hyp, &lPPA_mat, thresh, ncomb_sub);

    /*	for(i=0;i<r;i++)*/
    /*	{*/
    /*		for(j=0;j<nmodcol;j++) Rprintf("%d\t",models_num[index2_col(i,j,nmodcol)]);*/
    /*		Rprintf("%f\n",lPPA_mat[i]);*/
    /*	}*/
    /*	Rprintf("return_final=%d\n",return_final);*/
    /*	getchar();*/

    /*	for(i=0;i<r;i++) Rprintf("%f\n",lPPA_mat[i]);*/
    /*	Rprintf("\n");*/
    /*	Rprintf("lpriornull=%f\tlprioralt=%f\n",lpriornull,lprioralt);*/

    //add prior information
    if(string == 0)
    {
        for(i = 0; i < r; i++)
        {
            hyp[i] = calc_lessstring_single(nsamp, &models_num[index2_col(i, 0, nmodcol)]);
            lPPA_mat[i] += (hyp[i] == 0 ? lpriornull : lprioralt);
        }
    }
    else
    {
        for(i = 0; i < r; i++)
        {
            hyp[i] = calc_string_single(nsamp, &models_num[index2_col(i, 0, nmodcol)]);
            lPPA_mat[i] += (hyp[i] == 0 ? lpriornull : lprioralt);
        }
    }

    //sort all models into rank order according to log(PPA')
    temp_ind = (int *) Calloc(r, int);
    bubble_sort_dec(lPPA_mat, temp_ind, r);
    /*	for(i=0;i<r;i++) Rprintf("%f\n",lPPA_mat[i]);*/
    /*	Rprintf("\n");*/
    //generate sorted vectors/matrices
    new_models_num = (int *) Calloc(nmodcol * totmods, int);
    for(i = 0; i < r; i++) for(j = 0; j < nmodcol; j++) new_models_num[index2_col(i, j, nmodcol)] = models_num[index2_col(temp_ind[i], j, nmodcol)];
    Free(models_num);
    models_num = new_models_num;
    new_models_num = NULL;

    int *new_hyp = (int *) Calloc(totmods, int);
    for(i = 0; i < r; i++) new_hyp[i] = hyp[temp_ind[i]];
    Free(hyp);
    hyp = new_hyp;
    new_hyp = NULL;
    Free(temp_ind);

    //remove any below threshold
    if(return_final == 1)
    {
        i = 1;
        while(i < r && (lPPA_mat[0] - lPPA_mat[i]) < logc) i++;
        r = i;
    }
    else
    {
        i = 1;
        while(i < r && (lPPA_mat[0] - lPPA_mat[i]) < thresh) i++;
        r = i;
    }
    /*	for(i=0;i<nmod;i++)*/
    /*	{*/
    /*		for(j=0;j<nsamp;j++) Rprintf("%f\t",lPDM_int_mat[index2(i,j,nmod)]);*/
    /*		Rprintf("\n");*/
    /*	}*/
    /*	for(i=0;i<r;i++)*/
    /*	{*/
    /*		for(j=0;j<nmodcol;j++) Rprintf("%d\t",models_num[index2_col(i,j,nmodcol)]);*/
    /*		Rprintf("%f\n",lPPA_mat[i]);*/
    /*	}*/
    /*	Rprintf("return_final=%d\n",return_final);*/
    /*	getchar();*/

    //now make adjustments to model set if output is to be suppressed
    int *multfact = (int *) Calloc(totmods, int);
    int r1 = r;
    if(return_final == 0)
    {
        //resort lPDM_int_mat
        for(i = 0; i < totmods; i++) multfact[i] = 1;
        double temp_lPPA;
        int temp_hyp;
        int *currmods = (int *) R_alloc(nmodcol, sizeof(int));
        //start loop
        for(i = 0; i < nsamp; i++)
        {
            for(j = 0; j < nmod; j++)
            {
                if(uni_ind[index2(j, i, nmod)] > 0)
                {
                    for(k = 0; k < uni_ind[index2(j, i, nmod)]; k++)
                    {
                        for(m = 0; m < r; m++)
                        {
                            if(models_num[index2_col(m, i, nmodcol)] == j)
                            {
                                q = 0;
                                for(l = 0; l < nsamp; l++) if(models_num[index2_col(m, l + nsamp, nmodcol)] == models_num[index2_col(m, i + nsamp, nmodcol)]) q++;
                                if(q == 1)
                                {
                                    for(l = 0; l < nmodcol; l++) currmods[l] = models_num[index2_col(m, l, nmodcol)];
                                    currmods[i] = uni[index3(j, k, i, nmod, nmod)];
                                    temp_lPPA = lPPA_mat[m];
                                    if(string == 0) temp_hyp = calc_lessstring_single(nsamp, currmods);
                                    else temp_hyp = calc_string_single(nsamp, currmods);
                                    if(temp_hyp != hyp[m])
                                    {
                                        if(temp_hyp == 1) temp_lPPA += (lprioralt - lpriornull);
                                        else temp_lPPA += (lpriornull - lprioralt);
                                        //update model set
                                        for(l = 0; l < nmodcol; l++) models_num[index2_col(r1, l, nmodcol)] = currmods[l];
                                        lPPA_mat[r1] = temp_lPPA;
                                        hyp[r1] = temp_hyp;
                                        multfact[r1] = 1;
                                        r1++;
                                        //increase size of output vectors if necessary
                                        if(r1 > (totmods - 1)) realloc_approx_new(nmodcol, &totmods, incmods, &models_num, &hyp, &lPPA_mat, &multfact);
                                        //now enter further recursion if necessary
                                        if((i + 1) < nsamp)
                                        {
                                            //just check next sample is independent
                                            q = i + 1;
                                            d = 2;
                                            while(q < nsamp && d != 1)
                                            {
                                                d = 0;
                                                for(l = 0; l < nsamp; l++) if(currmods[l + nsamp] == currmods[q + nsamp]) d++;
                                                if(d == 1) final_lPPA_recur(nsamp, nmod, nmodcol, incmods, &totmods, q, &r1, r1, currmods, temp_lPPA, temp_hyp, &models_num, &hyp, &lPPA_mat, &multfact, uni, uni_ind, string, lpriornull, lprioralt);
                                                q++;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        multfact[m]++;
                                        //now enter further recursion if necessary
                                        if((i + 1) < nsamp)
                                        {
                                            //just check next sample is independent
                                            q = i + 1;
                                            d = 2;
                                            while(q < nsamp && d != 1)
                                            {
                                                d = 0;
                                                for(l = 0; l < nsamp; l++) if(currmods[l + nsamp] == currmods[q + nsamp]) d++;
                                                if(d == 1) final_lPPA_recur(nsamp, nmod, nmodcol, incmods, &totmods, q, &r1, m, currmods, temp_lPPA, temp_hyp, &models_num, &hyp, &lPPA_mat, &multfact, uni, uni_ind, string, lpriornull, lprioralt);
                                                q++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        /*		Rprintf("\n");*/
        /*		for(i=0;i<r1;i++)*/
        /*		{*/
        /*			for(j=0;j<nmodcol;j++) Rprintf("%d\t",models_num[index2_col(i,j,nmodcol)]);*/
        /*			Rprintf("%d\t%f\n",multfact[i],lPPA_mat[i]);*/
        /*		}*/
        /*		getchar();*/
        r = r1;
        //sort all models into rank order according to log(PPA')
        temp_ind = (int *) Calloc(r, int);
        bubble_sort_dec(lPPA_mat, temp_ind, r);
        //generate sorted vectors/matrices
        new_models_num = (int *) Calloc(nmodcol * totmods, int);
        for(i = 0; i < r; i++) for(j = 0; j < nmodcol; j++) new_models_num[index2_col(i, j, nmodcol)] = models_num[index2_col(temp_ind[i], j, nmodcol)];
        Free(models_num);
        models_num = new_models_num;
        new_models_num = NULL;

        new_hyp = (int *) Calloc(totmods, int);
        for(i = 0; i < r; i++) new_hyp[i] = hyp[temp_ind[i]];
        Free(hyp);
        hyp = new_hyp;
        new_hyp = NULL;

        int *new_multfact = (int *) Calloc(totmods, int);
        for(i = 0; i < r; i++) new_multfact[i] = multfact[temp_ind[i]];
        Free(multfact);
        multfact = new_multfact;
        new_multfact = NULL;
        Free(temp_ind);

        //now cycle through and remove models outside threshold
        i = 1;
        while(i < r && (lPPA_mat[0] - lPPA_mat[i]) < logc) i++;
        r = i;

        //reallocate memory to save space
        realloc_approx_new(nmodcol, &totmods, r - totmods, &models_num, &hyp, &lPPA_mat, &multfact);
    }
    else
    {
        //reallocate memory to save space
        realloc_approx(nmodcol, &totmods, r - totmods, &models_num, &hyp, &lPPA_mat);
        for(i = 0; i < r; i++) multfact[i] = 1;
    }

    //now export into R and return
    SEXP models_num_R;
    PROTECT(models_num_R = allocVector(REALSXP, totmods * nmodcol + 3 * totmods + 1));
    k = 0;
    for(i = 0; i < totmods; i++)
    {
        for(j = 0; j < nmodcol; j++)
        {
            REAL(models_num_R)[index2_col(i, j, nmodcol)] = (double) models_num[index2_col(i, j, nmodcol)];
            k++;
        }
    }
    for(j = 0; j < totmods; j++)
    {
        REAL(models_num_R)[k] = (double) hyp[j];
        k++;
    }
    for(j = 0; j < totmods; j++)
    {
        REAL(models_num_R)[k] = lPPA_mat[j];
        k++;
    }
    for(j = 0; j < totmods; j++)
    {
        REAL(models_num_R)[k] = (double) multfact[j];
        k++;
    }
    REAL(models_num_R)[k] = (double) totmods;
    UNPROTECT(1);
    //free memory from the heap (automatically sets pointer to NULL)
    Free(models_num);
    Free(hyp);
    Free(lPPA_mat);
    Free(multfact);
    return models_num_R;
}


