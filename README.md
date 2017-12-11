# quatSVD
An implementation of the fast3x3SVD code from [Wisconsin](http://pages.cs.wisc.edu/~sifakis/project_pages/svd.html) , that's comprehensible

## Usage
This library exposes an QuatSVD::svd(matrix) function that returns a struct with the singular values, S as a vec3, and U and V as quaternions

Should be very easy to use.  Just include quatSVD.hpp There are 2 interfaces you can use:  
### Row major matrix, as pointer
The struct contains 3 std::array<T>'s for S, U, and V.  The Quaternions are in w, x, y, z order

### Passing an Eigen matrix
you need to \#include the appropriate Eigen headers before you include quatSVD.hpp \#define INCLUDE_QUATSVD_EIGEN_API
The returned struct returns an Eigen::Vec3<T> for the singular values, and U and V as Eigen::Quaternion<T>'s 

### Performance/Correctness
The included tests.cpp file runs this code on a bunch of random matrices (with frobeneius norm 1) and prints some error/timing info.
On my machine (2015 rMBP 13") this code seems to be about 2x as fast as Eigen's JacobiSVD for doubles, and slightly less than 2x faster for floats

I suspect that even if you add in the cost of turning the quaternions into rotation matrices (if you need that), this will be faster than Eigen.

I wrote this because I needed to do SVDs that give U and V as quaternions and I knew of the UW code from above.  
That code is really, really hard to grok, so I wrote this with the hope of being easier to understand (while being hopefully almost as fast, although I didn't measure against their code).


