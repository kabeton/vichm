#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <functional>
#include "prog.h"
#include "spline.h"

//problem in general:
//y'' + p(x)y' + q(x)y = f(x)
//al1*y(a) + be1*y'(a) = g1
//al2*y(b) + be2*y'(b) = g2

template <class T>
void draw(std::vector<double> x, std::vector<double> y, T newt) {
  std::ofstream ofi("spline.dat");
  //std::ofstream giv("poly.dat");
  /*
  for(int i = 0; i < x.size(); i++) {
    giv << x[i] << " " << y[i] << std::endl;
  }
  */
  for(double i = x[0]; i <= x[x.size() - 1]; i += (x[1]-x[0])/10) {
    ofi << i << " " << newt.value(i) << std::endl;
  }

  FILE *gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set grid x y\n");
  fprintf(gp, "show grid\n");
  fprintf(gp, "plot 'spline.dat' with lines\n");
  fprintf(gp, "pause mouse close\n");
  pclose(gp);
}

double p(double x) {
  return x*x - 3;
}

double q(double x) {
  return (x*x - 3)*cos(x);
}

double f(double x) {
  return 2 - 6*x + (x*x - 3)*exp(x)*sin(x)*(1 + cos(x)) + cos(x)*(exp(x) + (x*x - 1) + x*x*x*x - 3*x*x);
}

void draw_vect(std::vector<double> y, double a, double b) {
  std::ofstream ofi("data.dat");
  double h = (b - a)/y.size();

  for(int i = 0; i < y.size(); i++) {
    ofi << (a + i*h) << " " << y[i] << std::endl;
  }

  FILE *gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set grid x y\n");
  fprintf(gp, "show grid\n");
  fprintf(gp, "plot 'data.dat' with lines\n");
  fprintf(gp, "pause mouse close\n");
  pclose(gp);
}

int main() {
  double pi = 3.1415926;
  double a = 0, b = pi;
  double al1 = 1, al2 = 1, be1 = 0, be2 = 0;
  double g1 = 0, g2 = pi*pi;

  Solver solver(p, q, f, al1, be1, al2, be2, g1, g2, a, b);

  double eps = 10e-3;
  std::vector<double> err;
  std::vector<double> sol;
  sol = solver.solve_prec(eps, err);
  //draw_vect(sol, a, b);
  std::vector<double> x(sol.size());
  double h = (b - a)/sol.size();

  for(int i = 0; i < sol.size(); i++) {
    x[i] = (a + i*h);
  }

  Spline s(x, sol);
  //draw(x, sol, s);

  std::cout << s.value(0.5) << " " << s.value(1) << " " << s.value(1.5) << " " << s.value(2) << " " << s.value(2.5) << " " << s.value(3) << std::endl;
  for(int i = 0; i < err.size(); i++) {
    std::cout << err[i] << std::endl;
  }
  
  std::vector<double> sums;
  for(int i = err.size() - 2; i >= 0; i--) {
    double s = 0;
    for(int j = err.size() - 1; j >= i; j--) {
      s += err[j];
    }
    sums.push_back(s);
  }
  std::cout << "------------------" << std::endl;
  for(int i = 0; i < sums.size(); i++) {
    std::cout << sums[i] << std::endl;
  }
  return 0;
}