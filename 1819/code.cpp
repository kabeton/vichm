#include <iostream>
#include <cmath>

double sin_ch(double x) {
  return x;
}

double fact(int number) {
  return (number <= 1) ? 1 : number * fact(number - 1);
}

double power(double a, int n) {
  if (n == 1) {
    return a;
  } else if (n == 0) {
    return 1;
  }
  return (n % 2 == 0) ? power(a, n / 2) * power(a, n / 2) : a * power(a, n - 1);
}

double sin_member(int degree, double point) {
  return power(-1, degree) * power(point, 2 * degree + 1) / fact(2 * degree + 1);
}

double exp_member(int degree, double point) {
  return power(point, degree) / fact(degree);
}

int main() {
  double p = 0.5, q = 1;  
  std::cout << "------------------------------------" << std::endl;
  std::cout << "p = 0.5 and q = 1: " << std::endl;
  std::cout << sin(p) << " " << sin(q) << std::endl;
  std::cout << exp(p) << " " << exp(q) << std::endl;
  std::cout << "------------------------------------" << std::endl;

  int i = 0;
  double ansi = sin_member(0, p), ansii = ansi + sin_member(1, p);
  while (fabs(ansi - ansii) > 10e-4) {
    ansi += sin_member(++i, p);
    ansii += sin_member(i + 1, p);
  }
  std::cout << i << std::endl;
  std::cout << ansi << " " << ansii << " " << ansi - ansii << std::endl;
  std::cout << "------------------------------------" << std::endl;

  i = 0;
  ansi = sin_member(0, q); ansii = ansi + sin_member(1, q);
  while (fabs(ansi - ansii) > 10e-4) {
    ansi += sin_member(++i, q);
    ansii += sin_member(i + 1, q);
  }
  std::cout << i << std::endl;
  std::cout << ansi << " " << ansii << " " << ansi - ansii << std::endl;
  std::cout << "------------------------------------" << std::endl;

  i = 0;
  ansi = exp_member(0, p); ansii = ansi + exp_member(1, p);
  while (fabs(ansi - ansii) > 10e-4) {
    ansi += exp_member(++i, p);
    ansii += exp_member(i + 1, p);
  }
  std::cout << i << std::endl;
  std::cout << ansi << " " << ansii << " " << ansi - ansii << std::endl;
  std::cout << "------------------------------------" << std::endl;

  i = 0;
  ansi = exp_member(0, q); ansii = ansi + exp_member(1, q);
  while (fabs(ansi - ansii) > 10e-4) {
    ansi += exp_member(++i, q);
    ansii += exp_member(i + 1, q);
  }
  std::cout << i << std::endl;
  std::cout << ansi << " " << ansii << " " << ansi - ansii << std::endl;
  std::cout << "------------------------------------" << std::endl;


  p = 10; q = 11;  
  std::cout << "------------------------------------" << std::endl;
  std::cout << "p = 10 and q = 11: " << std::endl;
  std::cout << sin(p) << " " << sin(q) << std::endl;
  std::cout << exp(p) << " " << exp(q) << std::endl;
  std::cout << "------------------------------------" << std::endl;

  i = 0;
  ansi = sin_member(0, p); ansii = ansi + sin_member(1, p);
  while (fabs(ansi - ansii) > 10e-4) {
    ansi += sin_member(++i, p);
    ansii += sin_member(i + 1, p);
  }
  std::cout << i << std::endl;
  std::cout << ansi << " " << ansii << " " << ansi - ansii << std::endl;
  std::cout << "------------------------------------" << std::endl;

  i = 0;
  ansi = sin_member(0, q); ansii = ansi + sin_member(1, q);
  while (fabs(ansi - ansii) > 10e-4) {
    ansi += sin_member(++i, q);
    ansii += sin_member(i + 1, q);
  }
  std::cout << i << std::endl;
  std::cout << ansi << " " << ansii << " " << ansi - ansii << std::endl;
  std::cout << "------------------------------------" << std::endl;

  i = 0;
  ansi = exp_member(0, p); ansii = ansi + exp_member(1, p);
  while (fabs(ansi - ansii) > 10e-4) {
    ansi += exp_member(++i, p);
    ansii += exp_member(i + 1, p);
  }
  std::cout << i << std::endl;
  std::cout << ansi << " " << ansii << " " << ansi - ansii << std::endl;
  std::cout << "------------------------------------" << std::endl;

  i = 0;
  ansi = exp_member(0, q); ansii = ansi + exp_member(1, q);
  while (fabs(ansi - ansii) > 10e-4) {
    ansi += exp_member(++i, q);
    ansii += exp_member(i + 1, q);
  }
  std::cout << i << std::endl;
  std::cout << ansi << " " << ansii << " " << ansi - ansii << std::endl;
  std::cout << "------------------------------------" << std::endl;

  return 0;
}