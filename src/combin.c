#include "functions.h"
#include <R.h>

//log-factorial function
double lfactorial(int n)
{
    int i;
    double ans = 0.0;
    for(i = n; i > 0; i--) ans += log(i);
    return ans;
}

//factorial function
int factorial(int n)
{
    int i;
    int ans = 1;
    for(i = n; i > 0; i--) ans = ans * i;
    return ans;
}

//choose function
int choose_c(int n, int k)
{
    int ans;
    if(n < k) ans = 0;
    else ans = factorial(n) / (factorial(k) * factorial(n - k));
    return ans;
}

