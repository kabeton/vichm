#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <functional>

double mpi(double x0, double eps, double (*meth)(double x), double (*F)(double y)) {
  double xp = x0;
  double xr = meth(x0);
  std::cout << F(xr) << std::endl;
  int iter = 1;
  while(fabs(F(xr)) > eps) {
    xp = xr;
    xr = meth(xp);
    iter++;
    std::cout << "iter: " << iter << " " << "guess: " << xr << " " << "func: " << F(xr) << std::endl;
  }
  return xr;
}

int main() {
  auto eq1{
    [](double x) -> double {
      return sqrt(0.5) * exp(x*x - 0.5)/2;
    }
  };
  
  auto F{
    [](double x) -> double {
      return x * exp(-x*x) - sqrt(0.5)*exp(-0.5)/2;
    }
  };
  
  auto eq2{
    [](double x) -> double {
      return sqrt(log(2*x) - log(sqrt(0.5)*exp(-0.5)));
    }
  };
  
  double eps = 10e-5;
  std::cout << std::setprecision(-log10(eps) + 1);
  double x1 = 1 / sqrt(8);
  double x2 = sqrt(2);
  double xl = mpi(x1, eps/2, eq1, F);
  double xr = mpi(x2, eps/2, eq2, F);
  std::cout << xl << " " << xr << std::endl;
  std::cout << xr - xl << std::endl;
  return 0;
}
