#include "functions.h"
#include <R.h>

//function for calculating integer powers
int pow_int(int base, int exp)
{
    int i, ans;
    if(exp == 0) ans = 1;
    else
    {
        ans = base;
        for(i = 0; i < (exp - 1); i++) ans = ans * base;
    }
    return ans;
}

