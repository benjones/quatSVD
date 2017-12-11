

#ifdef INCLUDE_QUATSVD_EIGEN_API
#include <Eigen/Eigen>
#endif

#include "quatSVD.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>

template<typename T>
std::array<T, 9> randomMatrix(){
  std::array<T,9> ret;
  
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<T> dis(0,1);

  T s2 = 0;
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

template<typename T>
void correctnessTests(int N){
  std::vector<std::array<T,9>> matrices;
  matrices.reserve(N);
  for(int i = 0; i < N; ++i){
	matrices.push_back(randomMatrix<T>());
  }
  
  std::vector<T> maxErrors, totalErrors;
  maxErrors.reserve(N);
  totalErrors.reserve(N);

  auto start = std::chrono::high_resolution_clock::now();
	
  for(int i = 0; i < N; ++i){
	
	auto errors = QuatSVD::computeErrors(matrices[i].data());
	
	totalErrors.push_back(errors.first);
	maxErrors.push_back(errors.second);
  }

  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = stop-start;
  std::cout << "time to compute and check " << N << " SVDs " << diff.count() << " s\n";
  
  
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



}

template<typename T>
void performanceTests(int N){
  std::vector<std::array<T,9>> matrices;
  matrices.reserve(N);
  for(int i = 0; i < N; ++i){
	matrices.push_back(randomMatrix<T>());
  }

  auto start = std::chrono::high_resolution_clock::now();

  double dummy = 0;
  for(int i = 0; i < N; ++i){
	auto SVD = QuatSVD::svd(matrices[i].data());
	dummy += SVD.S[0];
  }

  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = stop-start;
  std::cout << "time to compute " << N << " SVDs " << diff.count() << " s\n";
  std::cout << "dummy value to outsmart optimizer: " << dummy << std::endl;


}


#ifdef INCLUDE_QUATSVD_EIGEN_API
template<typename T>
Eigen::Matrix<T, 3, 3> randomEigenMatrix(){
  Eigen::Matrix<T, 3, 3> ret;
  
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<T> dis(0,1);

  T s2 = 0;
  for(int i = 0; i < 9; ++i){
	ret.data()[i] = dis(gen);
	s2 += ret.data()[i]*ret.data()[i];
  }
  
  s2 = std::sqrt(s2);
  
  for(int i = 0; i < 9; ++i){
	ret.data()[i] /= s2;
  }

  return ret;

}

template<typename T>
void eigenCorrectnessTests(int N){

  std::vector<Eigen::Matrix<T,3,3>> matrices;
  matrices.reserve(N);
  for(int i = 0; i < N; ++i){
	matrices.push_back(randomEigenMatrix<T>());
  }
  
  std::vector<T> maxErrors, totalErrors;
  maxErrors.reserve(N);
  totalErrors.reserve(N);

  auto start = std::chrono::high_resolution_clock::now();
	
  for(int i = 0; i < N; ++i){
	auto errors = QuatSVD::computeErrors(matrices[i]);
	totalErrors.push_back(errors.first);
	maxErrors.push_back(errors.second);
  }

  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = stop-start;
  std::cout << "time to compute and check " << N << " SVDs " << diff.count() << " s\n";
  
  
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

  
}

template<typename T>
void eigenPerformanceTests(int N){
  std::vector<Eigen::Matrix<T,3,3>> matrices;
  matrices.reserve(N);
  for(int i = 0; i < N; ++i){
	matrices.push_back(randomEigenMatrix<T>());
  }

  auto start = std::chrono::high_resolution_clock::now();

  double dummy = 0;
  for(int i = 0; i < N; ++i){
	auto SVD = QuatSVD::svd(matrices[i]);
	dummy += SVD.S[0];
  }

  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = stop-start;
  std::cout << "time to compute " << N << " SVDs " << diff.count() << " s\n";
  std::cout << "dummy value to outsmart optimizer: " << dummy << std::endl;


}


#endif

int main(){

  std::cout << std::setprecision(12);


  //  auto errs = computeErrors(matrix);
  //  std::cout << errs.first << ' ' << errs.second << std::endl;

  
  int N = 100000;
  std::cout << "doubles" << std::endl;
  correctnessTests<double>(N);
  performanceTests<double>(N);
  std::cout << "floats" << std::endl;
  correctnessTests<float>(N);
  performanceTests<float>(N);
  

#ifdef INCLUDE_QUATSVD_EIGEN_API
  std::cout << "tests with eigen" << std::endl;
  std::cout << "doubles " << std::endl;
  eigenCorrectnessTests<double>(N);
  eigenPerformanceTests<double>(N);
  
  std::cout << "floats " << std::endl;
  eigenCorrectnessTests<float>(N);
  eigenPerformanceTests<float>(N);

  #endif
  
  
  return 0;
}
