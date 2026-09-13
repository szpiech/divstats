#include "binom.h"

long double fact(int x)
{
    return floor(0.5 + exp(factln(x)));
}

//From Numerical Recipes in C
//returns n choose k
//
//Two hazards, both of which matter to the SFS subsampler:
// - it materializes the coefficient, so it overflows to inf once the value
//   leaves the range of long double -- n ~ 1030 where long double is 8 bytes
//   (arm64 macOS), n ~ 16000 where it is 80-bit (x86-64 Linux);
// - the rounding is not a no-op outside 0 <= k <= n. factln() returns 0 for
//   negative arguments, so nCk(a,b) with b > a evaluates a!/b!, and
//   floor(0.5 + a!/b!) rounds UP to 1 for (a,b) in {(0,1),(0,2),(1,2)} where
//   the true coefficient is 0.
//Prefer lnCk below for anything needing a RATIO of coefficients.
long double nCk(int n, int k)
{
    return floor(0.5 + exp(factln(n) - factln(k) - factln(n - k)));
}

//returns ln(n choose k), or -infinity where the coefficient is zero
//(k outside 0..n). exp(-infinity) is 0, so callers summing weights need no
//special case, and unlike nCk the zero is exact rather than rounded.
//Never overflows: ln C(n, n/2) grows like n*ln(2), so n = 10^6 lands near
//7*10^5.
long double lnCk(int n, int k)
{
    if (n < 0 || k < 0 || k > n) return -INFINITY;
    return factln(n) - factln(k) - factln(n - k);
}

//From Numerical Recipes in C
//returns ln(gamma(xx)) for xx >0
long double gammln(long double xx)
{
    long double x, y, tmp, ser;
    static long double cof[6] = {76.18009172947146, -86.50532032941677,
                                 24.01409824083091, -1.231739572450155,
                                 0.1208650973866179e-2, -0.5395239384953e-5
                                };
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

//From Numerical Recipes in C
//returns ln(n!)
long double factln(int n)
{
    static long double a[101];

    if (n <= 1) return 0.0;
    if (n <= 100) return a[n] ? a[n] : (a[n] = gammln(n + 1.0));
    else return gammln(n + 1.0);
}

