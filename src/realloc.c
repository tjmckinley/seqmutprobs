#include "functions.h"
#include <R.h>

//function to perform memory reallocation of intermediate structures if required
void realloc_maxmod(int nsamp, int *maxmodels, int incmodels, int **structure)
{
    /*'nsamp' is number of samples
    'maxmodels' records maximum number of models
    'incmodels' is value with which to increment 'maxmodels'
    'structure' is matrix requiring reallocation*/

    int i, j;
    int * temp_ptr = NULL;
    *maxmodels = *maxmodels + incmodels;
    temp_ptr = (int *) Calloc(*maxmodels * nsamp, int);
    for(i = 0; i < (*maxmodels - incmodels); i++) for(j = 0; j < nsamp; j++) temp_ptr[index2(i, j, *maxmodels)] = (*structure)[index2(i, j, *maxmodels - incmodels)];
    Free(*structure);
    *structure = temp_ptr;
    temp_ptr = NULL;
    return;
}

//function to perform memory reallocation of full model set if required
void realloc_maxmod_full(int nsamp, int *maxmodels, int incmodels, int **models_num)
{
    /*'nsamp' is number of samples
    'maxmodels' records maximum number of models
    'incmodels' is value with which to increment 'maxmodels'
    'models_num' is matrix requiring reallocation*/

    int i, j;
    int * temp_ptr = NULL;
    *maxmodels = *maxmodels + incmodels;
    temp_ptr = (int *) Calloc(*maxmodels * 2 * nsamp, int);
    for(i = 0; i < (*maxmodels - incmodels); i++) for(j = 0; j < (2 * nsamp); j++) temp_ptr[index2(i, j, *maxmodels)] = (*models_num)[index2(i, j, *maxmodels - incmodels)];
    Free(*models_num);
    *models_num = temp_ptr;
    temp_ptr = NULL;
    return;
}
