#include "functions.h"
#include <R.h>

/*bubble sort algorithm (sorts into decreasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_dec(double *values, int *ind, int n)
{
    /*'values' corresponds to unsorted vector
    'ind' corresponds to order indicators (used to sort other
    	vectors in the same order as 'values')
    'n' is vector length*/
    int i, j, tempind;
    double temp;
    for(i = 0; i < n; i++) ind[i] = i;
    i = n;
    while(i > 0)
    {
        for(j = 0; j < (i - 1); j++)
        {
            temp = values[j];
            tempind = ind[j];
            if(temp < values[j + 1])
            {
                values[j] = values[j + 1];
                values[j + 1] = temp;
                ind[j] = ind[j + 1];
                ind[j + 1] = tempind;
            }
        }
        i--;
    }
    return;
}

/*bubble sort algorithm for INTEGERS (sorts into decreasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_dec_int(int *values, int *ind, int n)
{
    /*'values' corresponds to unsorted vector
    'ind' corresponds to order indicators (used to sort other
    	vectors in the same order as 'values')
    'n' is vector length*/
    int i, j, tempind, temp;
    for(i = 0; i < n; i++) ind[i] = i;
    i = n;
    while(i > 0)
    {
        for(j = 0; j < (i - 1); j++)
        {
            temp = values[j];
            tempind = ind[j];
            if(temp < values[j + 1])
            {
                values[j] = values[j + 1];
                values[j + 1] = temp;
                ind[j] = ind[j + 1];
                ind[j + 1] = tempind;
            }
        }
        i--;
    }
    return;
}

/*bubble sort algorithm for INTEGERS (sorts into increasing order
and provides indicators in order to sort multiple columns)*/
void bubble_sort_inc_int(int *values, int *ind, int n)
{
    /*'values' corresponds to unsorted vector
    'ind' corresponds to order indicators (used to sort other
    	vectors in the same order as 'values')
    'n' is vector length*/
    int i, j, tempind, temp;
    for(i = 0; i < n; i++) ind[i] = i;
    i = n;
    while(i > 0)
    {
        for(j = 0; j < (i - 1); j++)
        {
            temp = values[j];
            tempind = ind[j];
            if(temp > values[j + 1])
            {
                values[j] = values[j + 1];
                values[j + 1] = temp;
                ind[j] = ind[j + 1];
                ind[j + 1] = tempind;
            }
        }
        i--;
    }
    return;
}

