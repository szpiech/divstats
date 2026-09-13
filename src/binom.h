#ifndef __BINOM_H__
#define __BINOM_H__
#include <cmath>

//Was declared as gammaln, but binom.cpp defines gammln. Nothing called the
//declared name, so the mismatch never failed to link.
long double gammln(long double);
long double factln(int);
long double nCk(int, int);
long double lnCk(int, int);
long double fact(int x);

#endif
