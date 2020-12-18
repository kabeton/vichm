#ifndef SPLINE_H
#define SPLINE_H 1

#include "poly.h"
#include <vector>

class Spline {
  std::vector<std::vector<double>> a;
  std::vector<double> grid;

  std::vector<double> solve_for_c(std::vector<double> a, std::vector<double> h) {
    std::vector<double> temp(a.size(), 0);
    double b = 0, f = 0;
    std::vector<double> p(a.size(), 0), q(a.size(), 0);
    for(int i = 1; i < a.size() - 1; i++) {
      b = 2 * (h[i-1] + h[i]);
      f = 3 * ((a[i+1] - a[i]) / h[i] - (a[i] - a[i-1]) / h[i-1]);
      p[i+1] = -h[i] / (h[i-1] * p[i] + b);
      q[i+1] = (f - h[i-1] * q[i]) / (h[i-1] * p[i] + b);
    }

    for(int i = a.size() - 2; i > 0; i--) {
      temp[i] = p[i+1]*temp[i+1] + q[i+1];
    }
    return temp;
  }

  std::vector<double> solve_for_b(std::vector<double> a, std::vector<double> c, std::vector<double> h) {
    std::vector<double> temp(a.size(), 0);
    for(int i = 1; i < a.size(); i++) {
      temp[i] = (a[i] - a[i-1])/h[i-1] + (2 * c[i] + c[i-1])*h[i-1]/3;
    }
    return temp;
  }

  std::vector<double> solve_for_d(std::vector<double> c, std::vector<double> h) {
    std::vector<double> temp(c.size(), 0);
    for(int i = 1; i < c.size(); i++) {
      temp[i] = (c[i] - c[i-1])/(3 * h[i-1]);
    }
    return temp;
  }

public:
  Spline();

  Spline(std::vector<double> x, std::vector<double> y) {
    grid = x;
    int n = x.size();
    std::vector<double> h(n-1);
    for(int i = 0; i < n-1; i++) {
      h[i] = x[i+1] - x[i];
    }
    std::vector<std::vector<double>> local(4);
    local[0] = y;
    local[2] = solve_for_c(local[0], h);
    local[1] = solve_for_b(local[0], local[2], h);
    local[3] = solve_for_d(local[2], h);
    for(int i = 0; i < y.size(); i++) {
      a.push_back({local[0][i], local[1][i], local[2][i], local[3][i]});
    }
  }

  Spline(Poly p, std::vector<double> x) {
    grid = x;
    Poly pd = p.derivative();
    for(int i = 0; i <= p.getdeg() - 1; i++) {
      double h = x[i+1] - x[i];
      std::vector<double> local(4);
      local[0] = (-1*pd(x[i+1])*x[i]*x[i]*x[i+1]*h + p(x[i+1])*x[i]*x[i]*(3*x[i+1] - x[i]))/(h*h*h) +
                 (p(x[i])*x[i+1]*x[i+1]*(x[i+1] - 3*x[i]) - pd(x[i])*x[i]*x[i+1]*x[i+1]*h)/(h*h*h);
      local[1] = (pd(x[i+1])*x[i]*(2*x[i+1] + x[i])*h - 6*(p(x[i+1]) - p(x[i]))*x[i]*x[i+1])/(h*h*h) + 
                 (pd(x[i])*x[i+1]*(x[i+1] + 2*x[i])*h)/(h*h*h);
      local[2] = (-1*pd(x[i+1])*h*(x[i+1] + 2*x[i]) + 3*(p(x[i+1]) - p(x[i]))*(x[i+1] + x[i]))/(h*h*h) -
                 (pd(x[i])*h*(x[i] + 2*x[i+1]))/(h*h*h);
      local[3] = (pd(x[i+1])*h - 2*(p(x[i+1]) - p(x[i])) + pd(x[i])*h)/(h*h*h);
      a.push_back(local);
    }
  }
  
  std::vector<double> get_ith_coef(int i) const {
    if(i < a.size()) return a[i];
    return {0};
  }

  int size() const {
    return a.size();
  }

  double value(double x) const {
    int i = 1;
    double start = grid[i];
    while(start < x && i < grid.size()-1) {
      start = grid[++i];
    }
    Poly c = Poly{get_ith_coef(i)};
    return c.value(x - grid[i]);
  }

  ~Spline() {};
};

std::ostream& operator<<(std::ostream &os, const Spline x) {
  for(int i = 0; i < x.size(); i++) {
    std::vector<double> c = x.get_ith_coef(i);
    os << "i=" << i << " " << c[0] << " " << c[1] << " x " << c[2] << " x^2 " << c[3] << " x^3" << std::endl;
  }
  return os;
}

#endif
