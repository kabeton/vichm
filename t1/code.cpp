#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <functional>
#define sol ;

struct point{
  double x, y, c, len;
  void update() {
    c = (x + y) / 2;
    len = y - x;
  }
};

std::ostream& operator<<(std::ostream& os, const point &p) {
#ifdef otd
  os << p.x << " " << p.y << std::endl;
#endif
#ifdef sol
  os << p.c << " " << p.len << std::endl;
#endif
  return os;
}

void dichotomy(double eps, double (*eq)(double x), std::vector<point> &roots) 
{
  for(int i = 0; i < roots.size(); i++) {
  int j = 0;
  while (roots[i].len > eps) {
    if(eq(roots[i].x) * eq(roots[i].c) < 0) {
      roots[i].y = roots[i].c;
    } else {
      roots[i].x = roots[i].c;
    }
    roots[i].update();
    j++;
  }
  std::cout << roots[i];
  std::cout << j << std::endl;
  }
}


int main() {
  //plot 2 graphs and locate roots
  //calculating 
#ifdef otd
  std::vector<point> gr1(4000);
  std::vector<point> gr2(2000);
  std::ofstream cf("circ.dat");
  std::ofstream tf("tang.dat");

  for(int i = 0; i < 2000; i++) {
    double x = -1 + (double)i/1000;
    gr1[i].x = x;
    gr1[i].y = sqrt(1-x*x);
    cf << gr1[i];
  }
  for(int i = 2000; i < 4000; i++) {
    double x = 1 - ((double)i-2000)/1000;
    gr1[i].x = x;
    gr1[i].y = -sqrt(1-x*x);
    cf << gr1[i];
  }
  for(int i = 0; i < 2000; i++) {
    double x = -1 + (double)i/1000;
    gr2[i].x = x;
    gr2[i].y = tan(x);
    tf << gr2[i];
  }

  //plotting
  FILE *gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set grid x y \n");
  fprintf(gp, "show grid\n");
  fprintf(gp, "plot 'circ.dat' with lines, 'tang.dat' with lines\n");
  fprintf(gp, "pause mouse close\n");
  
  pclose(gp);
#endif
  //solving the system
#ifdef sol
  std::vector<point> roots(2);
  roots[0].x = 0.6; roots[0].y = 0.7;
  roots[1].x = -0.7; roots[1].y = -0.6;
  roots[0].update(); roots[1].update();
  double eps = 0.000001;
  std::cout << std::setprecision(-log10(eps) + 1);

  auto eq{
    [](double x) -> double {
      return x*x + tan(x)*tan(x) - 1;
    }
  };

  dichotomy(eps, eq, roots);

#endif
  return 0;
}
