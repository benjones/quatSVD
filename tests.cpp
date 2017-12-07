#include "quatSVD.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

int main(){

  double matrix[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  /*double ds = 0;
  for(double d : matrix){
	ds += d*d;
  }

  /ds = std::sqrt(ds);
  for(double& d : matrix){
	d /= ds;
	}*/
  std::cout << std::setprecision(12);
  svd(matrix);

  return 0;
}
