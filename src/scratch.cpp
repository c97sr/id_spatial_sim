#include <iostream>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

int
main (void)
{
    double x = 5.0;
    double y = gsl_sf_bessel_J0 (x);
    cerr << ("J0(%g) = %.18e\n", x, y);
    return 0;
}

