#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "poly.h"
#include "spline.h"

std::vector<std::vector<double>> newton_div_sum(std::vector<double> x, std::vector<double> y) {
  int n = x.size();
  std::vector<std::vector<double>> ds(n, std::vector<double>(n));

  for (int k = 0; k < n; k++) {
    ds[0][k] = y[k];
  }

  for (int i = 1; i < n; i++) {
    for(int k = 0; k < n - i; k++) {
      ds[i][k] = (ds[i-1][k+1] - ds[i-1][k]) / (x[k+i] - x[k]);
    }
  }
  return ds;
}

Poly newton_interpol_alg(std::vector<double> x, std::vector<double> y) {
  int n = x.size();
  std::vector<std::vector<double>> ds = newton_div_sum(x, y);
  
  Poly inter, temp{1};
  for(int i = 0; i < n; i++) {
    inter = inter + ds[i][0] * temp;
    temp = temp * Poly{{-1 * x[i], 1}};
  }
  
  return inter;
}

template <class T>
void draw(std::vector<double> x, std::vector<double> y, T newt) {
  std::ofstream ofi("data.dat");
  //std::ofstream giv("poly.dat");
  /*
  for(int i = 0; i < x.size(); i++) {
    giv << x[i] << " " << y[i] << std::endl;
  }
  */
  for(double i = x[0]; i <= x[x.size() - 1] + 10; i += (x[1]-x[0])/10) {
    ofi << i << " " << newt.value(i) << std::endl;
  }

  FILE *gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set grid x y\n");
  fprintf(gp, "show grid\n");
  fprintf(gp, "plot 'data.dat' with lines, 'real.dat'\n");
  fprintf(gp, "pause mouse close\n");
  pclose(gp);
}

int main() 
{
  std::vector<double> x = {1910., 1920., 1930., 1940., 1950., 1960., 1970., 1980., 1990., 2000.};
  std::vector<double> y = {92228496., 106021537., 123202624., 132164569., 151325798., 179323175., 203211926., 226545805., 248709873., 281421906.};

  std::vector<double> lastx = {1970, 1980., 1990., 2000.};
  std::vector<double> lasty = {203211926., 226545805., 248709873., 281421906.};

  Poly newt = newton_interpol_alg(lastx, lasty);
  std::cout << newt;

  Spline s(x, y);
  std::cout << s;

  std::cout << newt(2010) << std::endl;
  std::cout << s.value(2010) << std::endl;

  draw(lastx, lasty, newt);
  draw(x, y, s);

  std::vector<std::vector<double>> ds = newton_div_sum(x, y);
  int n = x.size();
  FILE *out = fopen("newline.dat", "w");

  //calculating P for the point r

  for (double r = x[0]; r <= x.back()+10; r += (x[1]-x[0])/10) {
    double p = y[0];
    for (int i = 1; i < n; i++) {
      double add = ds[i][0];
      for (int k = 0; k < i; k++) {
        add *= r - x[k];
      }
      p += add;
    }
    fprintf(out, "%g %g\n", r, p);
  }
  fclose(out);

  FILE *gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set grid x y\n");
  fprintf(gp, "show grid\n");
  fprintf(gp, "plot 'newline.dat' with lines, 'real.dat'\n");
  fprintf(gp, "pause mouse close\n");
  pclose(gp);

  return 0;
}
