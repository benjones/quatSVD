#pragma once

/* compute the SVD of a 3x3 matrix and return 
   U, V as quaternions, and Sigma as a Vec3

   Based on the uwisc tech report from here: http://pages.cs.wisc.edu/~sifakis/project_pages/svd.html


   


*/

#include <cmath>
#include <x86intrin.h>
#include <array>


//IMPORT USAGE NOTE.  IF YOU WANT TO USE THIS WITH EIGEN MATRICES/QUATERNIONS
//UNCOMMENT THE NEXT LINE (You could add a -D flag for your build system... if you want
//#define INCLUDE_QUATSVD_EIGEN_API
//You should #include eigen headers yourself before including this file

namespace QuatSVD{
  /* this is the API (functions implemented below)
  template<typename T>
  inline SVD<T> svd(T* matrix, T eps = (.5*std::numeric_limits<T>::epsilon()));

  template<typename Derived, typename T = typename Derived::Scalar>
  inline EigenSVD<T> svd(const Eigen::DenseBase<Derived>& matrix,
  T eps = (.5*std::numeric_limits<T>::epsilon()));


  If you have a matrix that's already decomposed into its polar decomp A = RS (a quat times a symmetric matrix)
  We can speed up SVD computation.  It only requires a jacobi diagonalization of S and 
  a few quat multiplications/conjugations.

  This is currently only implemented for Eigen quats/matrices due to laziness

  template<typename Derived, typename T = typename Derived::Scalar>
  inline EigenSVD<T> svdFromPolar(
	  const Eigen::Quaternion<T>& R,  //R as a quaternion
	  const Eigen::DenseBase<Derived>& S,  //S as a 3x3 Eigen matrix
	  T eps = (.5*std::numeric_limits<T>::epsilon())){
  */


  
  template<typename T>
  struct SVD{
	std::array<T, 3> S;
	std::array<T, 4> U, V;
  };

#ifdef INCLUDE_QUATSVD_EIGEN_API

  template<typename T>
  struct EigenSVD{
	Eigen::Matrix<T, 3, 1> S;
	Eigen::Quaternion<T> U,V;
  };
  
#endif
}


//internal details
namespace{
  
  //from https://www.sebastiansylvan.com/post/scalarsseintrinsics/
  inline float rsqrt_fast(float x)
  {
	return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(x)));
  }
  inline float rsqrt(float x)
  {
	//the above + 1 step of newton's method
	float xrsqrt_est = rsqrt_fast(x);
	return xrsqrt_est*(1.5f - x*0.5f*xrsqrt_est*xrsqrt_est); // NR iteration
  }

  template<typename T>
	inline T fast_sqrt(T x){
	return x* rsqrt(x);
  }


  
  //stored as xx, yy, zz, xy, xz, yz
  template<typename T>
	struct Symm3x3{
	std::array<T, 6> arr;

	Symm3x3(T* a){  //assume A is row major for now
	  //form the normal equations A^T A
	  /*
		[a11 a21 a31; a12 a22 a32; a13 a23 a33] * 
		[a11 a12 a13; a21 a22 a23; a31 a32 a33]
		
		11 -> 0, 
		12 -> 1, 
		13 -> 2,
		21 -> 3,
		22 -> 4,
		23 -> 5,
		31 -> 6,
		32 -> 7,
		33 -> 8
	  */
	  auto& me = *this;
	  me(0,0) = a[0]*a[0] + a[3]*a[3] + a[6]*a[6];
	  me(1,1) = a[1]*a[1] + a[4]*a[4] + a[7]*a[7];
	  me(2,2) = a[2]*a[2] + a[5]*a[5] + a[8]*a[8];

	  me(0,1) = a[0]*a[1] + a[3]*a[4] + a[6]*a[7];
	  me(0,2) = a[0]*a[2] + a[3]*a[5] + a[6]*a[8];

	  me(1,2) = a[1]*a[2] + a[4]*a[5] + a[7]*a[8];
	  
	}


	/*
	  If we already have a symmetric matrix, in the correct format, just copy it in
	 */
	struct AlreadySymmetric{};

	Symm3x3(T* a, AlreadySymmetric){
	  std::copy(a, a+6, arr.data());
	}

#ifdef INCLUDE_QUATSVD_EIGEN_API
	template<typename Derived>
	Symm3x3(const Eigen::DenseBase<Derived>& a){
	  //TODO do these correctly
	  //	  static_assert(a.RowsAtCompileTime == 3, "must be a 3x3 eigen matrix type");
	  //	  static_assert(a.ColsAtCompileTime == 3, "must be a 3x3 eigen matrix type");

	  auto& me = *this;
	  me(0,0) = a(0,0)*a(0,0) + a(1,0)*a(1,0) + a(2,0)*a(2,0);
	  me(1,1) = a(0,1)*a(0,1) + a(1,1)*a(1,1) + a(2,1)*a(2,1);
	  me(2,2) = a(0,2)*a(0,2) + a(1,2)*a(1,2) + a(2,2)*a(2,2);

	  me(0,1) = a(0,0)*a(0,1) + a(1,0)*a(1,1) + a(2,0)*a(2,1);
	  me(0,2) = a(0,0)*a(0,2) + a(1,0)*a(1,2) + a(2,0)*a(2,2);
	  
	  me(1,2) = a(0,1)*a(0,2) + a(1,1)*a(1,2) + a(2,1)*a(2,2);

	  
	}

	//a is already symmetric, just copy it in
	template<typename Derived>
	Symm3x3(const Eigen::DenseBase<Derived>& a, AlreadySymmetric){
	  auto& me = *this;

	  me(0,0) = a(0,0);
	  me(1,1) = a(1,1);
	  me(2,2) = a(2,2);
	  me(0,1) = a(0,1);
	  me(0,2) = a(0,2);
	  me(1,2) = a(1,2);
	}

	
#endif
	

	//r must be less than c !!!
	constexpr inline T& operator()(int r, int c){
	  if(r == c){ return arr[r]; }
	  else if(r == 0){  return arr[2 + c]; }
	  else return arr[5];
	}
	//r must be less than c !!!
	constexpr inline T operator()(int r, int c) const {
	  if(r == c){ return arr[r]; }
	  else if(r == 0){  return arr[2 + c]; }
	  else return arr[5];
	}


	/*	inline void dump() const{
	  std::cout << (*this)(0,0) << '\t' << (*this)(0,1) << '\t' << (*this)(0,2) << '\n'
				<< (*this)(0,1) << '\t' << (*this)(1,1) << '\t' << (*this)(1,2) << '\n'
				<< (*this)(0,2) << '\t' << (*this)(1,2) << '\t' << (*this)(2,2) << std::endl;
	  
				}*/

	inline T frobeneius() const{
	  return arr[0]*arr[0] + arr[1]*arr[1] + arr[2]*arr[2] +
		2*(arr[3]*arr[3] + arr[4]*arr[4] + arr[5]*arr[5]);
	}

	inline T diagMag() const{
	  return arr[0]*arr[0] + arr[1]*arr[1] + arr[2]*arr[2];
	}
	
	//compute Q^T S Q
	//where Q is represented as the quaterion with c as the scalar and s in the slot that's not p or q
	inline void quatConjugate01(T c, T s);
	inline void quatConjugate02(T c, T s);
	inline void quatConjugate12(T c, T s);
	
	inline void quatConjugateFull(T* q){
	  //assume q is unit
	  /*
		R = [a b c; d e f; g h i]
	  */

	  auto a = 1 - 2*(q[2]*q[2] + q[3]*q[3]);
	  auto b = 2*(q[1]*q[2] - q[3]*q[0]);
	  auto c = 2*(q[1]*q[3] + q[2]*q[0]);

	  auto d = 2*(q[1]*q[2] + q[3]*q[0]);
	  auto e = 1 - 2*(q[1]*q[1] + q[3]*q[3]);
	  auto f = 2*(q[2]*q[3] - q[1]*q[0]);

	  auto g = 2*(q[1]*q[3] - q[2]*q[0]);
	  auto h = 2*(q[2]*q[3] + q[1]*q[0]);
	  auto i = 1 - 2*(q[1]*q[1] + q[2]*q[2]);


	  //B = Q^T S
	  const auto& me = *this;
	  auto B11 = a*me(0,0) + d*me(0,1) + g*me(0,2);
	  auto B12 = a*me(0,1) + d*me(1,1) + g*me(1,2);
	  auto B13 = a*me(0,2) + d*me(1,2) + g*me(2,2);
	  
	  auto B21 = b*me(0,0) + e*me(0,1) + h*me(0,2);
	  auto B22 = b*me(0,1) + e*me(1,1) + h*me(1,2);
	  auto B23 = b*me(0,2) + e*me(1,2) + h*me(2,2);
	  
	  auto B31 = c*me(0,0) + f*me(0,1) + i*me(0,2);
	  auto B32 = c*me(0,1) + f*me(1,1) + i*me(1,2);
	  auto B33 = c*me(0,2) + f*me(1,2) + i*me(2,2);

	  auto new11 = a*B11 + d*B12 + g*B13;
	  auto new22 = b*B21 + e*B22 + h*B23;
	  auto new33 = c*B31 + f*B32 + i*B33;
	  
	  auto new12 = a*B21 + d*B22 + g*B23;
	  auto new13 = a*B31 + d*B32 + g*B33;
		
	  auto new23 = b*B31 + e*B32 + h*B33;

	  (*this)(0,0) = new11;
	  (*this)(1,1) = new22;
	  (*this)(2,2) = new33;
	  (*this)(0,1) = new12;
	  (*this)(0,2) = new13;
	  (*this)(1,2) = new23;

	  
	}
	
	
  }; //end Symm3x3<T>
  
  
  template<typename T>
	inline void Symm3x3<T>::quatConjugate01(T c, T s){
	//rotate about the z axis
	/*
	  Q^T * *this * Q
	  [ c s 0; -s c 0; 0 0 1] * [ S ] * [c -s 0; s c 0; 0 0 1]
	  
	  = [ c s11 + s s12,   c s12 + s s22,   c s13 + s s23 ]
	  [ -s s11 + c s12,  -s s12 + c s22, -s s13 + c s23 ]
	  [ s13          ,   s23           ,    s33         ]   * [c -s 0; s c 0; 0 0 1]
	  
	  =  [ c c s11 + c s s12 + c s s12 + s s s22;  -c s s11 - s s s12 + c c s12 + c s s22;  c s13 + s s23]
	  [ -c s s11 + c c s12 -s s s12 + s c s22;  s s s11 - c s s12 - c s s12 + c c s22;  -s s13 + c s23]
	  [ c s13 + s s23                        ;  -s s13 + c s23                       ;  s33           ]
	*/
	auto realC = c*c - s*s;
	auto realS = 2*s*c;
	
	auto cs = realS*realC; //used a bunch of places
	auto cc = realC*realC;
	auto ss = realS*realS;

	const auto& me = *this;
	auto newS11 = cc * me(0,0) + 2*cs*me(0,1) + ss*me(1,1);
	auto newS22 = ss * me(0,0) - 2*cs*me(0,1) + cc*me(1,1);
	auto newS12 = me(0,1)*(cc - ss) + cs*( me(1,1) - me(0,0) );
	auto newS13 = realC*me(0,2) + realS*me(1,2);
	auto newS23 = realC*me(1,2) - realS*me(0,2);
	  
	(*this)(0,0) = newS11;
	(*this)(1,1) = newS22;
	(*this)(0,1) = newS12;
	(*this)(0,2) = newS13;
	(*this)(1,2) = newS23;
	  
  }

  template<typename T>
	inline void Symm3x3<T>::quatConjugate02(T c, T s){
	//rotate about the y axis
	//quat looks like (ch, 0, sh, 0)
	//R = [ 1 - 2*s*s, 0, 2*c*s;  0, 1, 0;  -2*c*s, 0, 1 - 2*s*s]
	//or:  [C, 0, S; 0 1 0; -S, 0, C]   where C = 1 - 2*s*s == c*c - s*s, and S = 2*c*s
	/*
	  Q^T * *this * Q
	  [ c 0 -s; - 1 -; s 0 c] * [ S ] * [c 0 s; 0 1 0; -s 0 c]]
	  
	  = [ c s11 - s s13,   c s12 - s s23,   c s13 - s s33 ]
	  [ s12,  s22, s23 ]
	  [ s*s11 + c*s13, s*s12 + c*s23, s*s13 + c*s33]   * [c 0 s; 0 1 0; -s 0 c]
	  
	  =  [ c c s11 - c s s13 - c s s13 + s s s33; c s12 - s s23;    c s s11 - s s s13 + c c s13 - c s s33]
	  [ c s12 - s s23;                            s22;              s s12 + c s23]
	  [ c s s11 + c c s13 - s s s13 - c s s33;      s s12 + c s23;  s s s11 + c s s13 + c s s13 + c c s33  ]
	*/
	auto realC = c*c - s*s;
	auto realS = 2*s*c;
	
	auto cs = realS*realC; //used a bunch of places
	auto cc = realC*realC;
	auto ss = realS*realS;
	const auto& me = *this;

	auto newS11 = cc*me(0,0) - 2*cs*me(0,2) + ss*me(2,2);
	auto newS33 = ss*me(0,0) + 2*cs*me(0,2) + cc*me(2,2);
	auto newS12 = realC*me(0,1) - realS*me(1,2);
	auto newS13 = cs*(me(0,0) - me(2,2))  + (cc - ss)*me(0,2);
	auto newS23 = realS*me(0,1) + realC*me(1,2);
	  
	(*this)(0,0) = newS11;
	(*this)(2,2) = newS33;
	(*this)(0,1) = newS12;
	(*this)(0,2) = newS13;
	(*this)(1,2) = newS23;
	
  }
  

  template<typename T>
	inline void Symm3x3<T>::quatConjugate12(T c, T s){
	//rotate about the x axis

	//quat looks like (ch, sh, 0, 0)
	//R = [ 1, 0, 0;  0, 1 - 2*s*s, -2cs ;0,   2cs,  1 - 2s*s   ]
	//or:  [1, 0, 0; 0 C -S; 0, S, C]   where C = 1 - 2*s*s == c*c - s*s, and S = 2*c*s


	
	/*
	  Q^T * *this * Q
	  [ 1 0 0;  0 c s; 0 -s c] * [ S ] * [1 0 0; 0 c -s; 0 s c]
	  
	  = [ s11, s12, s13]
	  [ c s12 + s s13, c s22 + s s23, c s23 + s s33  ]
	  [ -s s12 + c s13, -s s22 + c s23, -s s23 + c s33 ]   * [1 0 0; 0 c -s; 0 s c]
	  
	  =  [ s11,         c s12 + s s13,                           -s s12 + c s13  ]
	  [ c s12 + s s13,  c c s22 + c s s23 + c s s23 + s s s33,   -c s s22 - s s s23 + c c s23 + c s s33   ]
	  [ -s s12 + c s13,  -c s s22 + c c s23 -s s s23 + c s s33,  s s s22 - c s s23 - c s s23 + c c s33   ]
	  
	  
	*/

	auto realC = c*c - s*s;	
	auto realS = 2*s*c;
	
	
	auto cs = realS*realC; //used a bunch of places
	auto cc = realC*realC;
	auto ss = realS*realS;
	const auto& me = *this;
	auto newS22 = cc*me(1,1) + 2*cs*me(1,2) + ss*me(2,2);
	auto newS33 = ss*me(1,1) - 2*cs*me(1,2) + cc*me(2,2);
	auto newS12 = realC*me(0,1) + realS*me(0,2);
	auto newS13 = -realS*me(0,1) + realC*me(0,2);
	auto newS23 = (cc - ss)*me(1,2) + cs*(me(2,2) - me(1,1));
	
	(*this)(1,1) = newS22;
	(*this)(2,2) = newS33;
	(*this)(0,1) = newS12;
	(*this)(0,2) = newS13;
	(*this)(1,2) = newS23;
	
  }

	


  //variable templates are cool
  template<typename T>
	static constexpr T gamma =3 + 2*T{M_SQRT2};
  template<typename T>
	static const T sinBackup = T{.5}*std::sqrt(T{2} - T{M_SQRT2});
  template<typename T>
	static const T cosBackup = T{.5}*std::sqrt(T{2} + T{M_SQRT2});


  //B is A^T A , which we'll represent as a Vec 6
  //returns c, s, the givens angles
  template<typename T, int p, int q>
	inline std::pair<T,T> givensAngles(const Symm3x3<T>& B){
	
	T ch;// = B(p,p) - B(q,q);
	//sign is different depending on which element we're trying to eliminate

	//note if you don't have a C++ 17 compiler,
	//I suspect that if you remove "constexpr" that the optimzer will probably
	//give you the same code since these values are all compile time constants
	if constexpr(p == 0 && q == 1){ ch = (B(p,p) - B(q, q));}
	if constexpr(p == 0 && q == 2){ ch = (B(q,q) - B(p, p));}
	if constexpr(p == 1 && q == 2){ ch = (B(p,p) - B(q, q));}
	
	T sh = .5*B(p, q);                 

	T omega = rsqrt(ch*ch + sh*sh);
	ch *= omega;
	sh *= omega;
	
	bool approxValid = gamma<T>*sh*sh < ch*ch;
	
	ch = approxValid ? ch : cosBackup<T>;
	sh = approxValid ? sh : sinBackup<T>;
	
	return {ch, sh};
  }

  template<typename T>
	inline void quatTimesEquals(T* lhs, T* rhs){

	//store as s, x, y, z storage order
	/* 
	   s1 * v2 + s2* v1 + v1 x v2

	*/
	auto newS = lhs[0]*rhs[0] - lhs[1]*rhs[1] - lhs[2]*rhs[2] - lhs[3]*rhs[3];

	auto newX = lhs[0]*rhs[1] + rhs[0]*lhs[1] + lhs[2]*rhs[3] - rhs[2]*lhs[3];
	auto newY = lhs[0]*rhs[2] + rhs[0]*lhs[2] + lhs[3]*rhs[1] - rhs[3]*lhs[1];
	auto newZ = lhs[0]*rhs[3] + rhs[0]*lhs[3] + lhs[1]*rhs[2] - rhs[1]*lhs[2];

	lhs[0] = newS;
	lhs[1] = newX;
	lhs[2] = newY;
	lhs[3] = newZ;
	
  }

  template<typename T, int i>
	inline void quatTimesEqualCoordinateAxis(T* lhs, T c, T s){
	//the quat we're multiplying by is (c, ? s ?)  where s is in slot i of the vector part,
	//and the other entries are 0
	
	auto newS = lhs[0]*c - lhs[i+1]*s;
	
	T newVals[3];
	//the s2*v1 part
	newVals[0] = c*lhs[1];  
	newVals[1] = c*lhs[2];
	newVals[2] = c*lhs[3];
	//the s1*v2 part
	newVals[i] += lhs[0]*s;
	//the cross product part
	newVals[(i+1)%3] += s*lhs[1 + ((i+2)%3)];
	newVals[(i+2)%3] -= s*lhs[1 + ((i+1)%3)];
	
	lhs[0] = newS;
	lhs[1] = newVals[0];
	lhs[2] = newVals[1];
	lhs[3] = newVals[2];
	
  }

  /*  template<typename T>
	inline void quatDump(T* q){
	std::cout << "quat: " << q[0] << ' ' << q[1] << ' '<< q[2] << ' ' << q[3] << std::endl;
  }

  template<typename T>
	inline void matDump(T* m){
	std::cout << m[0] << ' ' << m[1] << ' ' << m[2] << '\n'
			  << m[3] << ' ' << m[4] << ' ' << m[5] << '\n'
			  << m[6] << ' ' << m[7] << ' ' << m[8] << std::endl;
			  }*/
  

  template<typename T>
	inline std::array<T, 4> jacobiDiagonalize(Symm3x3<T>& ATA){

  	std::array<T, 4> V {{1,0,0,0}};

	for(int i = 0; i < 4; ++i){
	  auto givens = givensAngles<T,0,1>(ATA);
	  ATA.quatConjugate01(givens.first, givens.second);
	  quatTimesEqualCoordinateAxis<T,2>(V.data(), givens.first, givens.second);

	  givens = givensAngles<T,1,2>(ATA);
	  ATA.quatConjugate12(givens.first, givens.second);
	  quatTimesEqualCoordinateAxis<T,0>(V.data(), givens.first, givens.second);

	  givens = givensAngles<T,0,2>(ATA);
	  ATA.quatConjugate02(givens.first, givens.second);
	  quatTimesEqualCoordinateAxis<T,1>(V.data(), givens.first, givens.second);

	}

	return V;
	
  }

  template<typename T>
	inline std::array<T, 9> computeAV(T* matrix, T* V){
	//compute quaternion matrix from V
	//V ->  [[a, b, c], [d, e, f], [g, h, i]]
	auto a = 1 - 2*(V[2]*V[2] + V[3]*V[3]);
	auto b = 2*(V[1]*V[2] - V[3]*V[0]);
	auto c = 2*(V[1]*V[3] + V[2]*V[0]);
	
	auto d = 2*(V[1]*V[2] + V[3]*V[0]);
	auto e = 1 - 2*(V[1]*V[1] + V[3]*V[3]);
	auto f = 2*(V[2]*V[3] - V[1]*V[0]);
	
	auto g = 2*(V[1]*V[3] - V[2]*V[0]);
	auto h = 2*(V[2]*V[3] + V[1]*V[0]);
	auto i = 1 - 2*(V[1]*V[1] + V[2]*V[2]);

	std::array<T, 9> ret;

	ret[0] = a*matrix[0]  + d*matrix[1]   + g*matrix[2] ;
	ret[1] = b*matrix[0]  + e*matrix[1]   + h*matrix[2] ;
	ret[2] = c*matrix[0]  + f*matrix[1]   + i*matrix[2] ;

	ret[3] = a*matrix[3]  + d*matrix[4]   + g*matrix[5] ;
	ret[4] = b*matrix[3]  + e*matrix[4]   + h*matrix[5] ;
	ret[5] = c*matrix[3]  + f*matrix[4]   + i*matrix[5] ;
	
	ret[6] = a*matrix[6]  + d*matrix[7]   + g*matrix[8] ;
	ret[7] = b*matrix[6]  + e*matrix[7]   + h*matrix[8] ;
	ret[8] = c*matrix[6]  + f*matrix[7]   + i*matrix[8] ;

	return ret;
	
  }

#ifdef INCLUDE_QUATSVD_EIGEN_API
  template<typename T, typename Derived>
	inline std::array<T, 9> computeAV(const Eigen::DenseBase<Derived>& matrix, T* V){

	//compute quaternion matrix from V
	//V ->  [[a, b, c], [d, e, f], [g, h, i]]
	auto a = 1 - 2*(V[2]*V[2] + V[3]*V[3]);
	auto b = 2*(V[1]*V[2] - V[3]*V[0]);
	auto c = 2*(V[1]*V[3] + V[2]*V[0]);
	
	auto d = 2*(V[1]*V[2] + V[3]*V[0]);
	auto e = 1 - 2*(V[1]*V[1] + V[3]*V[3]);
	auto f = 2*(V[2]*V[3] - V[1]*V[0]);
	
	auto g = 2*(V[1]*V[3] - V[2]*V[0]);
	auto h = 2*(V[2]*V[3] + V[1]*V[0]);
	auto i = 1 - 2*(V[1]*V[1] + V[2]*V[2]);

	std::array<T, 9> ret;

	ret[0] = a*matrix(0,0)  + d*matrix(0,1)   + g*matrix(0,2) ;
	ret[1] = b*matrix(0,0)  + e*matrix(0,1)   + h*matrix(0,2) ;
	ret[2] = c*matrix(0,0)  + f*matrix(0,1)   + i*matrix(0,2) ;

	ret[3] = a*matrix(1,0)  + d*matrix(1,1)   + g*matrix(1,2) ;
	ret[4] = b*matrix(1,0)  + e*matrix(1,1)   + h*matrix(1,2) ;
	ret[5] = c*matrix(1,0)  + f*matrix(1,1)   + i*matrix(1,2) ;
	
	ret[6] = a*matrix(2,0)  + d*matrix(2,1)   + g*matrix(2,2) ;
	ret[7] = b*matrix(2,0)  + e*matrix(2,1)   + h*matrix(2,2) ;
	ret[8] = c*matrix(2,0)  + f*matrix(2,1)   + i*matrix(2,2) ;

	
	return ret;
  }
#endif
  
  template <typename T, int i, int j>
	inline void swapColsNeg(T* B){

	auto tmp = -B[i];
	B[i] = B[j];
	B[j] = tmp;

	tmp = -B[i +3];
	B[i+3] = B[j +3];
	B[j+3] = tmp;

	tmp = -B[i+6];
	B[i+6] = B[j+6];
	B[j+6] = tmp;
  
  }

  template<typename T>
	inline void permuteColumns(T* B, T* V){

	T mags[3];
	mags[0] = B[0]*B[0] + B[3]*B[3] + B[6]*B[6];
	mags[1] = B[1]*B[1] + B[4]*B[4] + B[7]*B[7];
	mags[2] = B[2]*B[2] + B[5]*B[5] + B[8]*B[8];

	if(mags[0] < mags[1]){
	  swapColsNeg<T,0,1>(B);
	  //swaps cols 0 and 1 in the corresponding matrix form, negates the new col 1
	  quatTimesEqualCoordinateAxis<T,2>(V, T{M_SQRT1_2}, T{M_SQRT1_2});
	  std::swap(mags[0], mags[1]);
	}

	if(mags[0] < mags[2]){
	  swapColsNeg<T,0,2>(B);
	  quatTimesEqualCoordinateAxis<T,1>(V, T{M_SQRT1_2}, T{-M_SQRT1_2} );
	  std::swap(mags[0], mags[2]);
	}

	if(mags[1] < mags[2]){
	  swapColsNeg<T,1,2>(B);
	  quatTimesEqualCoordinateAxis<T,0>(V, T{M_SQRT1_2}, T{M_SQRT1_2});
	  //don't bother swapping the mags anymore... don't care
	}
	
	
  }

  //returns the 2 components of the quaternion
  //such that Q^T * B has a 0 in element p, q
  template<typename T, int r, int c>
	std::pair<T,T> computeGivensQR(T* B, T eps){

	auto app = B[4*c];
	auto apq = B[3*r + c];

	auto rho = std::sqrt(app*app + apq*apq);
	T sh = rho > eps ? apq : 0;
	T ch = std::abs(app) + std::max(rho, eps);

	if(app < 0){
	  std::swap(sh, ch);
	}

	auto omega = rsqrt(ch*ch + sh*sh);
	ch *= omega;
	sh *= omega;
	
	return {ch, sh};
  }


  //Q is the rot matrix defined by quaternion (ch, . . . sh .. . ) where sh is coord i
  template<typename T>
	inline void givensQTB2(T* B, T ch, T sh){
	//quat is (ch, 0, 0, sh), rotation around Z axis
	auto c = ch*ch - sh*sh;
	auto s = 2*sh*ch;
	//Q = [ c -s 0; s c 0; 0 0 1]
	
	auto newb00 = B[0]*c + B[3]*s;
	auto newb01 = B[1]*c + B[4]*s;
	auto newb02 = B[2]*c + B[5]*s;

	auto newb10 = 0;//B[3]*c - B[0]*s; //should be 0... maybe don't compute?
	auto newb11 = B[4]*c - B[1]*s;
	auto newb12 = B[5]*c - B[2]*s;

	B[0] = newb00;
	B[1] = newb01;
	B[2] = newb02;

	B[3] = newb10;
	B[4] = newb11;
	B[5] = newb12;

	
  }

  //This will be called after givensQTB<2>, so we know that
  //B10 is 0... which actually doesn't matter since that row won't change
  template<typename T>
	inline void givensQTB1(T* B, T ch, T sh){

	auto c = ch*ch - sh*sh;
	auto s = 2*sh*ch;
	//Q = [c 0 s; 0 1 0; -s 0 c];
	auto newb00 = B[0]*c - B[6]*s;
	auto newb01 = B[1]*c - B[7]*s;
	auto newb02 = B[2]*c - B[8]*s;

	auto newb20 = 0;// B[0]*s + B[6]*c; //should be 0... maybe don't compute?
	auto newb21 = B[1]*s + B[7]*c; 
	auto newb22 = B[2]*s + B[8]*c;

	B[0] = newb00;
	B[1] = newb01;
	B[2] = newb02;

	B[6] = newb20;
	B[7] = newb21;
	B[8] = newb22;
  }

  //B10 and B20 are 0, so don't bother filling in/computing them :)
  template<typename T>
	inline void givensQTB0(T* B, T ch, T sh){

	auto c = ch*ch - sh*sh;
	auto s = 2*ch*sh;


	/* we may not need to compute the off diags since B should be diagonal
	   after this step
	*/
	auto newb11 = B[4]*c + B[7]*s;
	//auto newb12 = B[5]*c + B[8]*s; 

	//auto newb21 = B[7]*c - B[4]*s;
	auto newb22 = B[8]*c - B[5]*s;

	B[4] = newb11;
	//B[5] = newb12;

	//B[7] = newb21;
	B[8] = newb22;
  }


  template<typename T>
	inline  std::pair<std::array<T, 3>, std::array<T, 4>>
	QRFactorize(T* AV, T* V, T eps){

	permuteColumns(AV, V);
	
	//perform QR decomposition, which is diagonalization here!
	std::array<T,4> U {{1, 0, 0, 0}};

	{
	  auto givens = computeGivensQR<T,1,0>(AV, eps);

	  givensQTB2(AV, givens.first, givens.second);
	  quatTimesEqualCoordinateAxis<T,2>(U.data(), givens.first, givens.second);
	}

	{
	  auto givens = computeGivensQR<T,2,0>(AV, eps);

	  //the sign of the sine should be swapped for this axis:
	  givensQTB1(AV, givens.first, -givens.second);
	  quatTimesEqualCoordinateAxis<T,1>(U.data(), givens.first, -givens.second);
	}
  
	{
	  auto givens = computeGivensQR<T,2,1>(AV, eps);
	
	  givensQTB0(AV, givens.first, givens.second);
	  quatTimesEqualCoordinateAxis<T,0>(U.data(), givens.first, givens.second);
	}

	std::array<T, 3> S{AV[0], AV[4], AV[8]};
	return  {S, U};
	
  }
}


namespace QuatSVD{

  template<typename T>
  inline SVD<T> svd(T* matrix, T eps = (.5*std::numeric_limits<T>::epsilon())){

	Symm3x3<T> ATA(matrix);

	auto V = jacobiDiagonalize(ATA);

	auto AV = computeAV(matrix, V.data());

	auto SU = QRFactorize(AV.data(), V.data(), eps);
	return SVD<T>{SU.first, SU.second, V};
  
  }

#ifdef INCLUDE_QUATSVD_EIGEN_API
  template<typename Derived, typename T = typename Derived::Scalar>
  inline EigenSVD<T> svd(const Eigen::DenseBase<Derived>& matrix,
	  T eps = (.5*std::numeric_limits<T>::epsilon())){
	
	Symm3x3<T> ATA(matrix);
	auto V = jacobiDiagonalize(ATA);
	auto AV = computeAV(matrix, V.data());

	auto SU = QRFactorize(AV.data(), V.data(), eps);

	EigenSVD<T> ret;
	ret.S.x() = SU.first[0];
	ret.S.y() = SU.first[1];
	ret.S.z() = SU.first[2];

	ret.U = Eigen::Quaternion<T>(SU.second[0], SU.second[1], SU.second[2], SU.second[3]);
	ret.V = Eigen::Quaternion<T>(V[0], V[1], V[2], V[3]);

	return ret;
  }

  /*
	If we have the polar decomposition of A, compute SVD from that more efficiently
	
   */
  template<typename Derived, typename T = typename Derived::Scalar>
  inline EigenSVD<T> svdFromPolar(
	  const Eigen::Quaternion<T>& R,
	  const Eigen::DenseBase<Derived>& S,
	  T eps = (.5*std::numeric_limits<T>::epsilon())){

	Symm3x3<T> SSym(S, typename Symm3x3<T>::AlreadySymmetric{});

	auto V = jacobiDiagonalize(SSym);  //S = V Sigma V^T, SSYM holds sigma

	//so, now RS = R * V * sigma * V^T
	EigenSVD<T> ret;
	ret.S(0,0) = SSym(0,0);
	ret.S(1,0) = SSym(1,1);
	ret.S(2,0) = SSym(2,2);

	ret.V = Eigen::Quaternion<T>(V[0], V[1], V[2], V[3]);
	ret.U = R*ret.V;

	return ret;
  }

  
#endif
  
  //U, V as quaternions in s, x, y, z order
  //S as a vec3
  template<typename T>
  inline std::array<T,9> reconstructMatrix(const SVD<T>& _svd){

	const T* S = _svd.S.data();
	const T* U = _svd.U.data();
	const T* V = _svd.V.data();
	
	//compute S * V^T 
  
	auto a = 1 - 2*(V[2]*V[2] + V[3]*V[3]);
	auto b = 2*(V[1]*V[2] - V[3]*V[0]);
	auto c = 2*(V[1]*V[3] + V[2]*V[0]);
  
	auto d = 2*(V[1]*V[2] + V[3]*V[0]);
	auto e = 1 - 2*(V[1]*V[1] + V[3]*V[3]);
	auto f = 2*(V[2]*V[3] - V[1]*V[0]);
  
	auto g = 2*(V[1]*V[3] - V[2]*V[0]);
	auto h = 2*(V[2]*V[3] + V[1]*V[0]);
	auto i = 1 - 2*(V[1]*V[1] + V[2]*V[2]);

	std::array<T,9> SVT;
	SVT[0] = S[0]*a;
	SVT[1] = S[0]*d;
	SVT[2] = S[0]*g;

	SVT[3] = S[1]*b;
	SVT[4] = S[1]*e;
	SVT[5] = S[1]*h;

	SVT[6] = S[2]*c;
	SVT[7] = S[2]*f;
	SVT[8] = S[2]*i;

	//compute U 
  
	a = 1 - 2*(U[2]*U[2] + U[3]*U[3]);
	b = 2*(U[1]*U[2] - U[3]*U[0]);
	c = 2*(U[1]*U[3] + U[2]*U[0]);
  
	d = 2*(U[1]*U[2] + U[3]*U[0]);
	e = 1 - 2*(U[1]*U[1] + U[3]*U[3]);
	f = 2*(U[2]*U[3] - U[1]*U[0]);
  
	g = 2*(U[1]*U[3] - U[2]*U[0]);
	h = 2*(U[2]*U[3] + U[1]*U[0]);
	i = 1 - 2*(U[1]*U[1] + U[2]*U[2]);
  
	std::array<T, 9> ret;

	ret[0] = SVT[0]*a + SVT[3]*b + SVT[6]*c;
	ret[1] = SVT[1]*a + SVT[4]*b + SVT[7]*c;
	ret[2] = SVT[2]*a + SVT[5]*b + SVT[8]*c;
  
	ret[3] = SVT[0]*d + SVT[3]*e + SVT[6]*f;
	ret[4] = SVT[1]*d + SVT[4]*e + SVT[7]*f;
	ret[5] = SVT[2]*d + SVT[5]*e + SVT[8]*f;
  
	ret[6] = SVT[0]*g + SVT[3]*h + SVT[6]*i;
	ret[7] = SVT[1]*g + SVT[4]*h + SVT[7]*i;
	ret[8] = SVT[2]*g + SVT[5]*h + SVT[8]*i;
  
	return ret;
  
  }

#ifdef INCLUDE_QUATSVD_EIGEN_API
  template<typename T>
  Eigen::Matrix<T, 3, 3> reconstructMatrix(const EigenSVD<T>& _svd){
	return _svd.U.toRotationMatrix()*_svd.S.asDiagonal()*
	  (_svd.V.toRotationMatrix().transpose());
  }
#endif

  
  template<typename T>
  inline std::pair<T,T> computeErrors(T* matrix){

	auto res = svd(matrix);

	auto reconstruction = reconstructMatrix(res);

	T errSquared = 0;
	T maxError = 0;
	for(int i = 0; i < 9; ++i){
	  auto e = matrix[i] - reconstruction[i];
	  e *= e;
	  errSquared += e;
	  if(e > maxError){ maxError = e;}
	}

	return {errSquared, maxError};
  
  }

#ifdef INCLUDE_QUATSVD_EIGEN_API
  template<typename Derived, typename T = typename Derived::Scalar>
  inline std::pair<T, T> computeErrors(const Eigen::DenseBase<Derived>& matrix){

	auto res = svd(matrix);
	auto reconstruction = reconstructMatrix(res);

	T errSquared = 0;
	T maxError = 0;
	for(int r = 0; r < 3; ++r){
	  for(int c = 0; c < 3; ++c){
		
		auto e = matrix(r,c) - reconstruction(r,c);
		e *= e;
		errSquared += e;
		if(e > maxError){ maxError = e;}
	  }
	}

	return {errSquared, maxError};
	
	
  }

  template<typename Derived, typename T = typename Derived::Scalar>
  inline std::pair<T, T> computeErrorsPolar(const Eigen::DenseBase<Derived>& matrix){

	auto res = svd(matrix);
	Eigen::Matrix<T,3,3> matCopy = matrix;
	Eigen::Quaternion<T> R = res.U*res.V.conjugate();
	Eigen::Matrix<T,3,3> S = res.V.toRotationMatrix()*res.S.asDiagonal()*res.V.conjugate().toRotationMatrix();

	Eigen::Matrix<T,3,3> reconstructedPolar = R.toRotationMatrix()*S;
	Eigen::Matrix<T,3,3> originalReconstructionError = reconstructedPolar - matCopy;
	
	auto pSVD = svdFromPolar(R, S);
	
	auto reconstruction = reconstructMatrix(pSVD);

	T errSquared = 0;
	T maxError = 0;
	for(int r = 0; r < 3; ++r){
	  for(int c = 0; c < 3; ++c){
		
		auto e = matrix(r,c) - reconstruction(r,c);
		e *= e;
		errSquared += e;
		if(e > maxError){ maxError = e;}
	  }
	}
	return {errSquared, maxError};
	
	
  }

  
#endif
  
}
