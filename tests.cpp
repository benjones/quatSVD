#include "quatSVD.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>


std::array<double, 9> randomMatrix(){
  std::array<double,9> ret;
  
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0,1);

  double s2 = 0;
  for(int i = 0; i < 9; ++i){
	ret[i] = dis(gen);
	s2 += ret[i]*ret[i];
  }
  
  s2 = std::sqrt(s2);
  
  for(int i = 0; i < 9; ++i){
	ret[i] /= s2;
  }

  return ret;
}


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


  //  auto errs = computeErrors(matrix);
  //  std::cout << errs.first << ' ' << errs.second << std::endl;


  int N = 1000000;
  std::vector<std::array<double,9>> matrices;
  matrices.reserve(N);
  for(int i = 0; i < N; ++i){
	matrices.push_back(randomMatrix());
  }
  
  std::vector<double> maxErrors, totalErrors;
  maxErrors.reserve(N);
  totalErrors.reserve(N);

  auto start = std::chrono::high_resolution_clock::now();
	
  for(int i = 0; i < 1000000; ++i){
	auto errors = computeErrors(matrices[i].data());
	//	std::cout << errors.first << ' ' << errors.second << std::endl;
	totalErrors.push_back(errors.first);
	maxErrors.push_back(errors.second);
  }

  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = stop-start;
  std::cout << "time to compute " << N << " SVDs " << diff.count() << " s\n";
  
  
  std::sort(maxErrors.begin(), maxErrors.end());
  std::sort(totalErrors.begin(), totalErrors.end());

  std::cout << "median total errors " << totalErrors[totalErrors.size()/2] << std::endl;
  std::cout << "median max errors " << totalErrors[maxErrors.size()/2] << std::endl;

  std::cout << "95th percential total errors "
			<< totalErrors[static_cast<int>(totalErrors.size()*.9)] << std::endl;
  std::cout << "95th percential max errors "
			<< totalErrors[static_cast<int>(maxErrors.size()*.9)] << std::endl;
  
  std::cout << "max total errors " << totalErrors.back() << std::endl;
  std::cout << "max max errors " << totalErrors.back() << std::endl;


  
  
  
  return 0;
}
