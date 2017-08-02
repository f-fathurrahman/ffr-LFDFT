#include <iostream>
#include <cstdio>
#include "Faddeeva.hpp"
#include "Faddeeva.cpp"

using namespace std;

int main()
{
  std::complex<double> z( 1.0, 1.0);
  std::complex<double> w_iz = Faddeeva::erfcx(z);
  
  printf("z    = %18.10f %18.10f\n", z.real(), z.imag());
  printf("w_iz = %18.10f %18.10f\n", w_iz.real(), w_iz.imag());
  
  return 0;
}
