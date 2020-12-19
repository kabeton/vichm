#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>

//problem in general:
//y'' + p(x)y' + q(x)y = f(x)
//al1*y(a) + be1*y'(a) = g1
//al2*y(b) + be2*y'(b) = g2

class Solver {
private:
  std::function<double(double)> p, q, f;
  double al1, be1, al2, be2, g1, g2, a, b;
public:
  Solver(std::function<double(double)> _p, std::function<double(double)> _q, std::function<double(double)> _f, 
         double _al1, double _be1, double _al2, double _be2, double _g1, double _g2, double _a, double _b) {
    p = _p; q = _q; f = _f; a = _a; b = _b;
    al1 = _al1; al2 = _al2; be1 = _be1; be2 = _be2; g1 = _g1; g2 = _g2;
  } 

  std::vector<double> gen_r(int n) {
    std::vector<double> temp(n);
    double h = (b - a) / (n-1);
    temp[0] = g1;
    temp[n-1] = g2;
    for(int i = 1; i < n-1; i++) {
      temp[i] = f(a + i*h);
    }
    return temp;
  }

  std::vector<std::vector<double>> gen_matrix(int n) {
    std::vector<std::vector<double>> temp(n, std::vector<double>(n));
    double h = (b - a)/(n-1);
    temp[0][0] = al1 - 1/h;
    temp[0][1] = be1/h;
    temp[n-1][n-2] = -be2/h;
    temp[n-1][n-1] = al2 + be2/h;
    for(int i = 1; i < n-1; i++) {
      temp[i][i-1] = 1/(h*h) - p(a+ i*h)/(2*h);
      temp[i][i] = q(a + i*h) - 2/(h*h);
      temp[i][i+1] = 1/(h*h) + p(a + i*h)/(2*h);
    }
    return temp;
  }

  std::vector<double> solve(std::vector<std::vector<double>> m, std::vector<double> r) {
    int n = r.size();
    std::vector<double> x(n);
    for(int i = 1; i < n-1; i++) {
      double w = m[i][i-1] / m[i-1][i-1];
      m[i][i] = m[i][i] - w*m[i-1][i];
      r[i] = r[i] - w*r[i-1];
    }
    x[n-1] = r[n-1]/m[n-1][n-1];
    for(int i = n-2; i >= 0; i--) {
      x[i] = (r[i] - m[i][i+1]*x[i+1])/m[i][i];
    }
    return x;
  }

  double dist(std::vector<double> x, std::vector<double> y) {
    int n = x.size();
    std::vector<double> d(n);
    for(int i = 0; i < n-1; i++) {
      d.push_back(fabs(x[i] - y[i*2]));
    }
    return *std::max_element(d.begin(), d.end());
  }

  std::vector<double> solve_prec(double eps, std::vector<double> &err) {
    int n = 10;
    std::vector<double> r1 = gen_r(n);
    std::vector<double> r2 = gen_r(2*n);
    std::vector<std::vector<double>> m1 = gen_matrix(n);
    std::vector<std::vector<double>> m2 = gen_matrix(2*n);
    std::vector<double> res1(n, 0), res2(n, 0);
    res1 = solve(m1, r1);
    res2 = solve(m2, r2);
    double errl = dist(res1, res2);
    while(errl > eps) {
      n *= 2;
      r1 = r2; r2 = gen_r(2*n);
      m1 = m2; m2 = gen_matrix(2*n);
      res1 = res2; res2 = solve(m2, r2);
      err.push_back(dist(res1, res2));
      errl = dist(res1, res2);
    }
    return res2;
  }
};