#pragma once

/* compute the SVD of a 3x3 matrix and return 
   U, V as quaternions, and Sigma as a Vec3

   Based on the uwisc tech report from here: http://pages.cs.wisc.edu/~sifakis/project_pages/svd.html

*/

#include <cmath>
#include <x86intrin.h>
#include <array>
#include <cassert>
#include <iostream>


//from https://www.sebastiansylvan.com/post/scalarsseintrinsics/
inline float rsqrt_fast(float x)
{
  return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(x)));
}
inline float rsqrt(float x)
{
  float xrsqrt_est = rsqrt_fast(x);
  return xrsqrt_est*(1.5f - x*0.5f*xrsqrt_est*xrsqrt_est); // NR iteration
}



//internal details
namespace{

  //stored as xx, yy, zz, xy, xz, yz
  struct Symm3x3{
	std::array<double, 6> arr;

	Symm3x3(double* a){  //assume A is row major for now
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


	
	constexpr inline double& operator()(int r, int c){
	  assert(r <= c);
	  if(r == c){ return arr[r]; }
	  else if(r == 0){  return arr[2 + c]; }
	  else return arr[5];
	}

	constexpr inline double operator()(int r, int c) const {
	  assert(r <= c);
	  if(r == c){ return arr[r]; }
	  else if(r == 0){  return arr[2 + c]; }
	  else return arr[5];
	}


	inline void dump() const{
	  std::cout << (*this)(0,0) << '\t' << (*this)(0,1) << '\t' << (*this)(0,2) << '\n'
				<< (*this)(0,1) << '\t' << (*this)(1,1) << '\t' << (*this)(1,2) << '\n'
				<< (*this)(0,2) << '\t' << (*this)(1,2) << '\t' << (*this)(2,2) << std::endl;

	}

	inline double frobeneius() const{
	  return arr[0]*arr[0] + arr[1]*arr[1] + arr[2]*arr[2] +
		2*(arr[3]*arr[3] + arr[4]*arr[4] + arr[5]*arr[5]);
	}

	inline double diagMag() const{
	  return arr[0]*arr[0] + arr[1]*arr[1] + arr[2]*arr[2];
	}
	
	//compute Q^T S Q
	//where Q is represented as the quaterion with c as the scalar and s in the slot that's not p or q
	template<int p, int q>
	void quatConjugate(double c, double s);

	inline void quatConjugateFull(double* q){
	  //assume q is unit
	  /*
		R = [a b c; d e f; g h i]
	  */

	  //	  std::cout << "qnorm in full conju: " <<
	  //		(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]) << std::endl;
		
	  
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
	
	
  };
  
  
  template<>
	inline void Symm3x3::quatConjugate<0,1>(double c, double s){
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

  template<>
	inline void Symm3x3::quatConjugate<0,2>(double c, double s){
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
  

  template<>
	inline void Symm3x3::quatConjugate<1,2>(double c, double s){
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

	


  
  static constexpr double gamma =3 + 2*M_SQRT2;
  static const double sinBackup = std::sin(M_PI/8);
  static const double cosBackup = std::cos(M_PI/8);


  //B is A^T A , which we'll represent as a Vec 6
  //returns c, s, the givens angles
  template<int p, int q>
	inline std::pair<double,double> givensAngles(const Symm3x3& B){
	std::cout << "computing angles for p: " << p << " q: " << q << std::endl;
	
	double ch;// = B(p,p) - B(q,q);
	//sign is different depending on which element we're trying to eliminate
	if constexpr(p == 0 && q == 1){ ch = (B(p,p) - B(q, q));}
	if constexpr(p == 0 && q == 2){ ch = (B(q,q) - B(p, p));}
	if constexpr(p == 1 && q == 2){ ch = (B(p,p) - B(q, q));}
	
	double sh = .5*B(p, q);                 

	std::cout << "guesses for ch " << ch << " sh " << sh
			  << " norm: " << ch*ch + sh*sh << std::endl;

	double omega = rsqrt(ch*ch + sh*sh);
	ch *= omega;
	sh *= omega;
	
	bool approxValid = gamma*sh*sh < ch*ch;
	
	ch = approxValid ? ch : cosBackup;
	sh = approxValid ? sh : sinBackup;
	
	std::cout << "approx valid? " << approxValid << std::endl;
	std::cout << "ch " << ch << " sh " << sh << " norm: " << ch*ch + sh*sh << std::endl;
	
	return {ch, sh};
  }


  inline void quatTimesEquals(double* lhs, double* rhs){

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

  template<int i>
	inline void quatTimesEqualCoordinateAxis(double* lhs, double c, double s){
	//the quat we're multiplying by is (c, ? s ?)  where s is in slot i of the vector part,
	//and the other entries are 0
	
	auto newS = lhs[0]*c - lhs[i+1]*s;
	
	float newVals[3];
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
  
  inline void quatNormalize(double* q){
	double qSum = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
	auto norm = std::sqrt(qSum);
	q[0] /= norm;
	q[1] /= norm;
	q[2] /= norm;
	q[3] /= norm;
  }

  inline void quatDump(double* q){
	std::cout << "quat: " << q[0] << ' ' << q[1] << ' '<< q[2] << ' ' << q[3] << std::endl;
  }


  inline std::array<double, 4> jacobiDiagonalize(double* matrix){
	Symm3x3 ATA(matrix);
	ATA.dump();
	std::cout << ATA.frobeneius() << std::endl;
	std::cout << ATA.diagMag() << std::endl;
  	std::array<double, 4> V {{1,0,0,0}};

	for(int i = 0; i < 4; ++i){


	  std::cout << std::endl << std::endl;
	  auto givens = givensAngles<0,1>(ATA);
	  ATA.quatConjugate<0,1>(givens.first, givens.second);
	  ATA.dump();
	  std::cout << ATA.frobeneius() << std::endl;
	  std::cout << ATA.diagMag() << std::endl;

	  //	std::cout << "quat before" << std::endl;
	  //	quatDump(V.data());
	  quatTimesEqualCoordinateAxis<2>(V.data(), givens.first, givens.second);
	  //	std::cout << "quat after" << std::endl;
	  quatDump(V.data());
	  /*	{
			std::array<double,4> qUndo = V;
			quatTimesEqualCoordinateAxis<2>(qUndo.data(), givens.first, -givens.second);
			std::cout << "times quatInverse: " << std::endl;
			quatDump(qUndo.data());
			}*/

	  /*	{
			std::cout << "unconjugating" << std::endl;
			Symm3x3 check = ATA;
			check.quatConjugate<0,1>(givens.first, -givens.second);
			check.dump();
			std::cout << std::endl;
			}*/

	
	  /*{
		Symm3x3 check(matrix);
		check.quatConjugateFull(V.data());
		std::cout << "check" << std::endl;
		check.dump();
		std::cout << "check frob: " << check.frobeneius() << std::endl;
		std::cout << "diag: " << check.diagMag() << std::endl;
		std::cout << std::endl;
		}*/

	  std::cout << std::endl << std::endl;
	  givens = givensAngles<1,2>(ATA);
	  ATA.quatConjugate<1,2>(givens.first, givens.second);
	  ATA.dump();
	  std::cout << ATA.frobeneius() << std::endl;
	  std::cout << ATA.diagMag() << std::endl;

	  //	std::cout << "quat before" << std::endl;
	  //	quatDump(V.data());

	  quatTimesEqualCoordinateAxis<0>(V.data(), givens.first, givens.second);
	  std::cout << "quat after " << std::endl;
	  quatDump(V.data());

	  /*{
		std::array<double,4> qUndo = V;
		quatTimesEqualCoordinateAxis<0>(qUndo.data(), givens.first, -givens.second);
		std::cout << "times quatInverse: " << std::endl;;
	
		quatDump(qUndo.data());
		}
		{
		Symm3x3 check(matrix);
		check.quatConjugateFull(V.data());
		std::cout << "check" << std::endl;
		check.dump();
		std::cout << "check frob: " << check.frobeneius() << std::endl;
		std::cout << "diag: " << check.diagMag() << std::endl;
	  
		std::cout << std::endl;
		}*/
	  std::cout << std::endl << std::endl;
	  givens = givensAngles<0,2>(ATA);
	  ATA.quatConjugate<0,2>(givens.first, givens.second);
	  ATA.dump();
	  std::cout << ATA.frobeneius() << std::endl;
	  std::cout << ATA.diagMag() << std::endl;


	  //	std::cout << "quat before" << std::endl;
	  //	quatDump(V.data());

	  quatTimesEqualCoordinateAxis<1>(V.data(), givens.first, givens.second);
	  std::cout << "quat after " << std::endl;
	  quatDump(V.data());
	  /*	{
			std::array<double,4> qUndo = V;
			quatTimesEqualCoordinateAxis<1>(qUndo.data(), givens.first, -givens.second);
			std::cout << "times quatInverse: " << std::endl;;
			quatDump(qUndo.data());
			}
			{
			Symm3x3 check(matrix);
			check.quatConjugateFull(V.data());
			std::cout << "check" << std::endl;
			check.dump();
			std::cout << "check frob: " << check.frobeneius() << std::endl;
			std::cout << "diag: " << check.diagMag() << std::endl;
	  
			std::cout << std::endl;
			}*/
	}

	Symm3x3 check(matrix);
	check.quatConjugateFull(V.data());
	std::cout << "final check" << std::endl;
	check.dump();

	return V;
	
  }

  inline std::array<double, 9> computeAV(double* matrix, double* V){
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

	std::array<double, 9> ret;

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


  template <int i, int j>
	inline void swapColsNeg(double* B){

	double tmp = -B[i];
	B[i] = B[j];
	B[j] = tmp;

	tmp = -B[i +3];
	B[i+3] = B[j +3];
	B[j+3] = tmp;

	tmp = -B[i+6];
	B[i+6] = B[j+6];
	B[j+6] = tmp;
  
  }
  
  inline void permuteColumns(double* B, double* V){

	double mags[3];
	mags[0] = B[0]*B[0] + B[3]*B[3] + B[6]*B[6];
	mags[1] = B[1]*B[1] + B[4]*B[4] + B[7]*B[7];
	mags[2] = B[2]*B[2] + B[5]*B[5] + B[8]*B[8];

	if(mags[0] < mags[1]){
	  swapColsNeg<0,1>(B);
	  //swaps cols 0 and 1 in the corresponding matrix form, negates the new col 1
	  quatTimesEqualCoordinateAxis<2>(V, M_SQRT1_2, M_SQRT1_2);
	  std::swap(mags[0], mags[1]);
	}

	if(mags[0] < mags[2]){
	  swapColsNeg<0,2>(B);
	  quatTimesEqualCoordinateAxis<1>(V, M_SQRT1_2, -M_SQRT1_2 );
	  std::swap(mags[0], mags[2]);
	}

	if(mags[1] < mags[2]){
	  swapColsNeg<1,2>(B);
	  quatTimesEqualCoordinateAxis<0>(V, M_SQRT1_2, M_SQRT1_2);
	  //don't bother swapping the mags anymore... don't care
	}
	
	
  }

  
}


inline void svd(double* matrix){

  auto VQuat = jacobiDiagonalize(matrix);

  auto AV = computeAV(matrix, VQuat.data());


  std::cout << AV[0] << ' ' << AV[1] << ' ' << AV[2] << '\n'
			<< AV[3] << ' ' << AV[4] << ' ' << AV[5] << '\n'
			<< AV[6] << ' ' << AV[7] << ' ' << AV[8] << std::endl;

  std::cout << AV[0]*AV[1] + AV[3]*AV[4] + AV[6]*AV[7] << std::endl;
  std::cout << AV[0]*AV[2] + AV[3]*AV[5] + AV[6]*AV[8] << std::endl;
  std::cout << AV[2]*AV[1] + AV[5]*AV[4] + AV[8]*AV[7] << std::endl;

  permuteColumns(AV.data(), VQuat.data());

  std::cout << "after permute columns " << std::endl;
  
  std::cout << AV[0] << ' ' << AV[1] << ' ' << AV[2] << '\n'
			<< AV[3] << ' ' << AV[4] << ' ' << AV[5] << '\n'
			<< AV[6] << ' ' << AV[7] << ' ' << AV[8] << std::endl;

  std::cout << AV[0]*AV[1] + AV[3]*AV[4] + AV[6]*AV[7] << std::endl;
  std::cout << AV[0]*AV[2] + AV[3]*AV[5] + AV[6]*AV[8] << std::endl;
  std::cout << AV[2]*AV[1] + AV[5]*AV[4] + AV[8]*AV[7] << std::endl;


  auto AVCheck = computeAV(matrix, VQuat.data());
  std::cout << "AVCheck: " << std::endl;

  std::cout << AVCheck[0] << ' ' << AVCheck[1] << ' ' << AVCheck[2] << '\n'
			<< AVCheck[3] << ' ' << AVCheck[4] << ' ' << AVCheck[5] << '\n'
			<< AVCheck[6] << ' ' << AVCheck[7] << ' ' << AVCheck[8] << std::endl;

  
}