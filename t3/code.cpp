#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <functional>

std::vector<double> gen_grid(double x0, double xn, int points) {
  std::vector<double> res(points);
  double h = (xn - x0) / points;
  for(double i = 0; i < xn; i += h) {
    res.push_back(i);
  }
  return res;
}

double step(double a, double b, std::function<double(double)> f) {
  return (b - a) * (f(a) + 3*f((2*a+b)/3) + 3*f((a+2*b)/3) + f(b))/8;
}

double integr(std::vector<double> grid, std::function<double(double)> f) {
  double res = 0;
  for(int i = 0; i < grid.size() - 1; i++) {
    res += step(grid[i], grid[i+1], f);
  }
  return res;
}

void solve(double &res, double &err, double x0, double xn, double eps, std::function<double(double)> f) {
  int n = 10;
  double r1 = 0, r2 = 0;
  std::vector<double> grid1 = gen_grid(x0, xn, n);
  std::vector<double> grid2 = gen_grid(x0, xn, n*2);
  r1 = integr(grid1, f);
  r2 = integr(grid2, f);
  err = fabs(r1 - r2) / 15;
  while (err > eps) {
    n *= 2;
    grid1 = grid2;
    r1 = r2;
    std::vector<double> grid2 = gen_grid(x0, xn, n*2);
    r2 = integr(grid2, f);
    err = fabs(r1 - r2) / 15;
  }
  std::cout << "final n: " << n << std::endl;
  res = r2;
}

double eq(double x) {
  return sin(100*x) * exp(-x*x) * cos(2*x);
}

int main() {
  std::function<double(double)> f = eq;
  double eps = 10e-15;
  std::cout << std::setprecision(-log10(eps) + 1);
  double res = 0, err = 0;
  double x0 = 0, xn = 3;

  solve(res, err, x0, xn, eps, f);
  std::cout << res << " +- " << err << std::endl;

  return 0;
}