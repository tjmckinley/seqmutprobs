#include "functions.h"
#include <R.h>

/*This code to calculate combinations was adapted from a post by scvalex at
http://compprog.wordpress.com/2007/10/17/generating-combinations-1/.
It's very elegant and saved me a headache---many thanks!*/
int next_comb(int *comb, int k, int n)
{
    /*'comb' is vector of length 'k', initialised
    	outside of this function
    'k' is number of combinations to choose
    'n' is number of elements to choose 'k' items from*/

    int i = k - 1;
    comb[i]++;
    while((i >= 0) && (comb[i] >= n - k + 1 + i))
    {
        i--;
        if(i >= 0) comb[i]++;
    }
    //end if no more combinations can be generated
    if(comb[0] > n - k) return 0;
    //else adjust comb vector
    for(i = i + 1; i < k; i++) comb[i] = comb[i - 1] + 1;
    return 1;
}

