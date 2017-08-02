#include <iostream>
#include "Faddeeva.hpp"
#include "Faddeeva.cpp"

using namespace std;

extern "C"
{
  void cwrap_faddeeva_( double *re, double *im, double *f_re, double *f_im );
}

void cwrap_faddeeva_( double *re, double *im, double *f_re, double *f_im )
{
  std::complex<double> z(*re, *im);
  std::complex<double> w_iz = Faddeeva::erfcx(z);

  *f_re = w_iz.real();
  *f_im = w_iz.imag();
}
