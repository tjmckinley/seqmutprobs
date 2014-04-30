#include "functions.h"
#include <R.h>

//function to generate top subset of fully and patially dependent models based on an intermediate set of structures
void genfullmodels_approx_recur(int nhier, int currhier1, int *temp, int *tempind, int *currmods, int q, int nrow, int ncol, int *structure, int *r, int *totmods, int incmods, int nmodcol, double *lPDM_int_mat, int *lPDM_int_ind, int ntotcol, int **models_num, int **hyp, double **lPPA_mat, double logthresh, int *ncomb_sub, double curr_lPPA, double *max_lPPA)
{
    /*'nhier' is number of hierarchies
    'currhier1' is current hierarchy
    'temp' and 'tempind' are vectors recording structural information
    'currmods' is model structure for current model
    'q' is specific row of 'structure' we wish to extract
    'nrow' and 'ncol' are using for indexing (based on 'structure')
    'structure' is a matrix containing intermediate model structures
    'r' denotes current model in 'models_num'
    'totmods' is maximum number of models in 'models_num'
    'incmods' is arbitrary number of models to increase if reallocation necessary
    'nmodcol' is number of columns of 'models_num'
    'lPDM_int_mat' and 'lPDM_int_ind' are intermediate matrices (in vector form)
    	containing log[P'(D|M)] information
    'ntotcol' is number of columns of 'lPDM_int_ind' and 'lPDM_int_mat'
    'models_num' is a matrix (in vector form) recording final models
    'hyp' is a binary vector of length 'incmods', with 0=null and 1=alt
    'lPPA_mat' is an output vector of log[P'(M|D)]s corresponding to the set of output models
    'logthresh' is the log-threshold relative to the maximum model
    'ncomb_sub' is intermediate vector used for indexing
    'curr_lPPA' is current log-PDM at previous hierarchy
    'max_lPPA' is current maximum log-PDM*/

    //declare necessary variables
    int i, j, k, nsamp, nmod = 10, valid;
    nsamp = ncol;

    //work out variables relating to maximum independent model
    double temp_lPPA;
    int corrcol;
    int * currmods1 = (int *) Calloc(nmodcol, int);
    int * temp_index = (int *) Calloc(nsamp, int);

    int currhier = currhier1 + 1;
    //if only independent structures left to evaluate
    if(temp[currhier] == 1)
    {
        corrcol = tempind[currhier];
        i = 0;
        valid = 0;
        while(i < nmod && valid == 0)
        {
            //reset model structure
            for(j = 0; j < nmodcol; j++) (*models_num)[index2_col(*r, j, nmodcol)] = currmods[j];
            (*models_num)[index2_col(*r, corrcol, nmodcol)] = lPDM_int_ind[index2(i, corrcol, nmod)];
            //change model
            temp_lPPA = curr_lPPA - lPDM_int_mat[index2(0, corrcol, nmod)] + lPDM_int_mat[index2(i, corrcol, nmod)];
            //if final choice then record model and update system
            if(currhier == (nhier - 1))
            {
                //check if new model is selected
                if(((*max_lPPA) - temp_lPPA) < logthresh)
                {
                    (*lPPA_mat)[*r] = temp_lPPA;
                    *r = *r + 1;
                    //increase size of output vectors if necessary
                    if(*r > (*totmods - 1)) realloc_approx(nmodcol, totmods, incmods, models_num, hyp, lPPA_mat);
                    /*if new model is maximum model then recalculate log-acceptance threshold then reset maximum model*/
                    if(temp_lPPA > (*max_lPPA)) *max_lPPA = temp_lPPA;
                }
                else valid = 1;
            }
            else
            {
                for(j = 0; j < nmodcol; j++) currmods1[j] = (*models_num)[index2_col(*r, j, nmodcol)];
                //enter further hierarchy
                genfullmodels_approx_recur(nhier, currhier, temp, tempind, currmods1, q, nrow, ncol, structure, r, totmods, incmods, nmodcol, lPDM_int_mat, lPDM_int_ind, ntotcol, models_num, hyp, lPPA_mat, logthresh, ncomb_sub, temp_lPPA, max_lPPA);
            }
            i++;
        }
    }
    else
    {
        //calculate correct column of lPDM_int_mat based on structure of current hierarchy
        corrcol = calc_corr_col(q, nrow, structure, nsamp, ncomb_sub, tempind[currhier], temp[currhier], temp_index);
        //calculate PPAs
        for(i = 0; i < nmod; i++)
        {
            /*			if(lPDM_int_ind[index2(i,corrcol,nmod)]>0)*/
            /*			{*/
            //reset model structure
            for(j = 0; j < nmodcol; j++) (*models_num)[index2_col(*r, j, nmodcol)] = currmods[j];
            //now fill in gaps and initialise calculation of changes to lPDM
            temp_lPPA = curr_lPPA;
            for(k = 0; k < nsamp; k++)
            {
                if(structure[index2(q, k, nrow)] == tempind[currhier])
                {
                    (*models_num)[index2_col(*r, k, nmodcol)] = lPDM_int_ind[index2(i, corrcol, nmod)];
                    (*models_num)[index2_col(*r, k + nsamp, nmodcol)] = tempind[currhier];
                    temp_lPPA -= lPDM_int_mat[index2(0, k, nmod)];
                }
            }
            //adjust lPDM to correct structure
            temp_lPPA += lPDM_int_mat[index2(i, corrcol, nmod)];
            //if final choice then record model and update system
            if(currhier == (nhier - 1))
            {
                //check if new model is selected
                if(((*max_lPPA) - temp_lPPA) < logthresh)
                {
                    (*lPPA_mat)[*r] = temp_lPPA;
                    *r = *r + 1;
                    //increase size of output vectors if necessary
                    if(*r > (*totmods - 1)) realloc_approx(nmodcol, totmods, incmods, models_num, hyp, lPPA_mat);
                    /*if new model is maximum model then recalculate log-acceptance threshold then reset maximum model*/
                    if(temp_lPPA > (*max_lPPA)) *max_lPPA = temp_lPPA;
                }
            }
            else
            {
                for(j = 0; j < nmodcol; j++) currmods1[j] = (*models_num)[index2_col(*r, j, nmodcol)];
                //enter further hierarchy
                genfullmodels_approx_recur(nhier, currhier, temp, tempind, currmods1, q, nrow, ncol, structure, r, totmods, incmods, nmodcol, lPDM_int_mat, lPDM_int_ind, ntotcol, models_num, hyp, lPPA_mat, logthresh, ncomb_sub, temp_lPPA, max_lPPA);
            }
            /*			}*/
        }
    }
    //free memory from the heap (automatically sets pointers to NULL)
    Free(currmods1);
    Free(temp_index);
    return;
}

