/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.23  2005/03/15 11:40:56  cignoni
Added operator*=( std::vector<PointType> ...) to apply a matrix to a vector of vertexes (replacement of the old style mesh.Apply(tr).

Revision 1.22  2004/12/15 18:45:50  tommyfranken
*** empty log message ***

Revision 1.21  2004/10/22 14:41:30  ponchio
return in operator+ added.

Revision 1.20  2004/10/18 15:03:14  fiorin
Updated interface: all Matrix classes have now the same interface

Revision 1.19  2004/10/07 14:23:57  ganovelli
added function to take rows and comlumns. Added toMatrix and fromMatrix to comply
RotationTYpe prototype in Similarity.h

Revision 1.18  2004/05/28 13:01:50  ganovelli
changed scalar to ScalarType

Revision 1.17  2004/05/26 15:09:32  cignoni
better comments in set rotate

Revision 1.16  2004/05/07 10:05:50  cignoni
Corrected abuse of for index variable scope

Revision 1.15  2004/05/04 23:19:41  cignoni
Clarified initial comment, removed vector*matrix operator (confusing!)
Corrected translate and Rotate, removed gl stuff.

Revision 1.14  2004/05/04 02:34:03  ganovelli
wrong use of operator [] corrected

Revision 1.13  2004/04/07 10:45:54  cignoni
Added: [i][j] access, V() for the raw float values, constructor from T[16]

Revision 1.12  2004/03/25 14:57:49  ponchio

****************************************************************************/

#ifndef __VCGLIB_MATRIX44
#define __VCGLIB_MATRIX44

#include <memory.h>
#include <vcg/math/base.h>
#include <vcg/space/point3.h>
#include <vcg/space/point4.h>
#include <vector>


namespace vcg {

  /*
	Annotations:
Opengl stores matrix in  column-major order. That is, the matrix is stored as:

	a0  a4  a8  a12
	a1  a5  a9  a13
	a2  a6  a10 a14
	a3  a7  a11 a15

  Usually in opengl (see opengl specs) vectors are 'column' vectors 
  so usually matrix are PRE-multiplied for a vector.
  So the command glTranslate generate a matrix that 
  is ready to be premultipled for a vector:

	1 0 0 tx
	0 1 0 ty 
	0 0 1 tz
	0 0 0  1

Matrix44 stores matrix in row-major order i.e.

	a0  a1  a2  a3
	a4  a5  a6  a7
	a8  a9  a10 a11
	a12 a13 a14 a15

So for the use of that matrix in opengl with their supposed meaning you have to transpose them before feeding to glMultMatrix.
This mechanism is hidden by the templated function defined in wrap/gl/math.h; 
If your machine has the ARB_transpose_matrix extension it will use the appropriate;
The various gl-like command SetRotate, SetTranslate assume that you are making matrix
for 'column' vectors.

*/

  /** This class represent a 4x4 matrix. T is the kind of element in the matrix.
    */
template <class T> class Matrix44 {  
protected:
  T _a[16];

public:	
  typedef T ScalarType;

///@{

  /** $name Contrutors
    *  No automatic casting and default constructor is empty
    */
	Matrix44() {};	
	~Matrix44() {};
  Matrix44(const Matrix44 &m);
  Matrix44(const T v[]);

	///	Number of columns
	inline unsigned int ColumnsNumber() const
	{
		return 4;
	};

	/// Number of rows
	inline unsigned int RowsNumber() const
	{
		return 4;
	};

  T &ElementAt(const int row, const int col);
  T ElementAt(const int row, const int col) const;
  //T &operator[](const int i);
  //const T &operator[](const int i) const;
  T *V();
  const T *V() const ;
  
  T *operator[](const int i);
  const T *operator[](const int i) const;

	// return a copy of the i-th column
	Point4<T> GetColumn(const int& i)const{
		assert(i >=0);
		assert(i<4);
		int first = i<<2;
		return Point4<T>(_a[first],_a[first+1],_a[first+2],_a[first+3]);
	}

	// return the i-th row
	Point4<T> & GetColumn4(const int& i)const{
		assert(i >=0);
		assert(i<4);
		int first = i<<2;
		return Point4<T>(_a[first],_a[first+4],_a[first+8],_a[first+12]);
	}

	// return a copy of the i-th row
	Point4<T> Row4(const int& i)const{
		assert(i >=0);
		assert(i<4);
		return *((Point4<T>*)(&_a[i<<2]));
	}

	Point3<T> GetColumn3(const int& i)const{
		assert(i >=0);
		assert(i<4);
		int first  = i <<2;
		return Point3<T>(_a[first],_a[first+4],_a[first+8]);
	}	
	
	// return a copy of the i-th row
	Point3<T> Row3(const int& i)const{
		assert(i >=0);
		assert(i<4);
		return *((Point3<T>*)(&_a[i<<2]));
	}


  Matrix44 operator+(const Matrix44 &m) const;
  Matrix44 operator-(const Matrix44 &m) const;
  Matrix44 operator*(const Matrix44 &m) const;
  Point4<T> operator*(const Point4<T> &v) const;  

  bool operator==(const  Matrix44 &m) const;
  bool operator!= (const  Matrix44 &m) const;

  Matrix44 operator-() const;			
  Matrix44 operator*(const T k) const;
  void operator+=(const Matrix44 &m);
  void operator-=(const Matrix44 &m);	
  void operator*=( const Matrix44 & m );	
  void operator*=( const T k );
	
  template <class Matrix44Type>
	void ToMatrix(Matrix44Type & m) const {for(int i = 0; i < 16; i++) m.V()[i]=V()[i];}
	template <class Matrix44Type>
	void FromMatrix(const Matrix44Type & m){for(int i = 0; i < 16; i++) V()[i]=m.V()[i];}
	void SetZero();
  void SetIdentity();
  void SetDiagonal(const T k);
  Matrix44 &SetScale(const T sx, const T sy, const T sz);	
	Matrix44 &SetTranslate(const Point3<T> &t);
	Matrix44 &SetTranslate(const T sx, const T sy, const T sz);
  ///use radiants for angle.
  Matrix44 &SetRotate(T AngleRad, const Point3<T> & axis); 

  T Determinant() const;

  template <class Q> void Import(const Matrix44<Q> &m) {
    for(int i = 0; i < 16; i++) 
      _a[i] = (T)(m.V()[i]);
  }
};


/** Class for solving A * x = b. */
template <class T> class LinearSolve: public Matrix44<T> {
public:
  LinearSolve(const Matrix44<T> &m);
  Point4<T> Solve(const Point4<T> &b); // solve A · x = b 
  ///If you need to solve some equation you can use this function instead of Matrix44 one for speed.
  T Determinant() const;
protected:
  ///Holds row permutation.
  int index[4]; //hold permutation
  ///Hold sign of row permutation (used for determinant sign)
  T d;
  bool Decompose();
};

/*** Postmultiply */
//template <class T> Point3<T> operator*(const Point3<T> &p, const Matrix44<T> &m);

///Premultiply 
template <class T> Point3<T> operator*(const Matrix44<T> &m, const Point3<T> &p);

template <class T> Matrix44<T> &Transpose(Matrix44<T> &m);
//return NULL matrix if not invertible
template <class T> Matrix44<T> &Invert(Matrix44<T> &m);   
template <class T> Matrix44<T> Inverse(const Matrix44<T> &m);   

typedef Matrix44<short>  Matrix44s;
typedef Matrix44<int>    Matrix44i;
typedef Matrix44<float>  Matrix44f;
typedef Matrix44<double> Matrix44d;



template <class T> Matrix44<T>::Matrix44(const Matrix44<T> &m) {
  memcpy((T *)_a, (T *)m._a, 16 * sizeof(T));   
}

template <class T> Matrix44<T>::Matrix44(const T v[]) {
  memcpy((T *)_a, v, 16 * sizeof(T));   
}

template <class T> T &Matrix44<T>::ElementAt(const int row, const int col) {
  assert(row >= 0 && row < 4);
  assert(col >= 0 && col < 4);
  return _a[(row<<2) + col];
}

template <class T> T Matrix44<T>::ElementAt(const int row, const int col) const {
  assert(row >= 0 && row < 4);
  assert(col >= 0 && col < 4);
  return _a[(row<<2) + col];
}

//template <class T> T &Matrix44<T>::operator[](const int i) {
//  assert(i >= 0 && i < 16);
//  return ((T *)_a)[i];
//}
//
//template <class T> const T &Matrix44<T>::operator[](const int i) const {
//  assert(i >= 0 && i < 16);
//  return ((T *)_a)[i];
//}
template <class T> T *Matrix44<T>::operator[](const int i) {
  assert(i >= 0 && i < 16);
  return _a+i*4;
}

template <class T> const T *Matrix44<T>::operator[](const int i) const {
  assert(i >= 0 && i < 4);
  return _a+i*4;
}
template <class T>  T *Matrix44<T>::V()  { return _a;} 
template <class T> const T *Matrix44<T>::V() const { return _a;} 


template <class T> Matrix44<T> Matrix44<T>::operator+(const Matrix44 &m) const {
  Matrix44<T> ret;
  for(int i = 0; i < 16; i++) 
    ret.V()[i] = V()[i] + m.V()[i];  
  return ret;
}

template <class T> Matrix44<T> Matrix44<T>::operator-(const Matrix44 &m) const {
  Matrix44<T> ret;
  for(int i = 0; i < 16; i++) 
    ret.V()[i] = V()[i] - m.V()[i];  
  return ret;
}

template <class T> Matrix44<T> Matrix44<T>::operator*(const Matrix44 &m) const {
  Matrix44 ret;       
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++) {
      T t = 0.0;
      for(int k = 0; k < 4; k++)
        t += ElementAt(i, k) * m.ElementAt(k, j);
      ret.ElementAt(i, j) = t;
    }
  return ret;
}

template <class T> Point4<T> Matrix44<T>::operator*(const Point4<T> &v) const {
  Point4<T> ret;     
  for(int i = 0; i < 4; i++){
    T t = 0.0;
    for(int k = 0; k < 4; k++)
      t += ElementAt(i,k) * v[k];
    ret[i] = t;
   }
   return ret;
}


template <class T> bool Matrix44<T>::operator==(const  Matrix44 &m) const {
  for(int i = 0 ; i < 16; i++)
    if(operator[](i) != m[i])
      return false;
  return true;
}
template <class T> bool Matrix44<T>::operator!=(const  Matrix44 &m) const {
  for(int i = 0 ; i < 16; i++)
    if(operator[](i) != m[i])
      return true;
  return false;
}

template <class T> Matrix44<T> Matrix44<T>::operator-() const {
  Matrix44<T> res;
  for(int i = 0; i < 16; i++)
    res.V()[i] = -V()[i];
  return res;
}

template <class T> Matrix44<T> Matrix44<T>::operator*(const T k) const {
  Matrix44<T> res;
  for(int i = 0; i < 16; i++)
    res.V()[i] =V()[i] * k;
  return res;
}

template <class T> void Matrix44<T>::operator+=(const Matrix44 &m) {
  for(int i = 0; i < 16; i++)
    V()[i] += m.V()[i];  
}
template <class T> void Matrix44<T>::operator-=(const Matrix44 &m) {
  for(int i = 0; i < 16; i++)
    V()[i] -= m.V()[i];
}
template <class T> void Matrix44<T>::operator*=( const Matrix44 & m ) {
  *this = *this *m;
  
  /*for(int i = 0; i < 4; i++) { //sbagliato
    Point4<T> t(0, 0, 0, 0);
    for(int k = 0; k < 4; k++) {
      for(int j = 0; j < 4; j++) {
        t[k] += ElementAt(i, k) * m.ElementAt(k, j);
      }
    }
    for(int l = 0; l < 4; l++)
      ElementAt(i, l) = t[l];
  } */
}

template < class PointType , class T > void operator*=( std::vector<PointType> &vert, const Matrix44<T> & m ) {
  typename std::vector<PointType>::iterator ii;
  for(ii=vert.begin();ii!=vert.end();++ii)
    (*ii).P()=m * (*ii).P();
}

template <class T> void Matrix44<T>::operator*=( const T k ) {
  for(int i = 0; i < 4; i++)
    operator[](i) *= k;
}


template <class T> void Matrix44<T>::SetZero() {
  memset((T *)_a, 0, 16 * sizeof(T));
}

template <class T> void Matrix44<T>::SetIdentity() { 
  SetDiagonal(1);  
}

template <class T> void Matrix44<T>::SetDiagonal(const T k) {
  SetZero();
  ElementAt(0, 0) = k;
  ElementAt(1, 1) = k;
  ElementAt(2, 2) = k;
  ElementAt(3, 3) = 1;    
}

template <class T> Matrix44<T> &Matrix44<T>::SetScale(const T sx, const T sy, const T sz) {
  SetZero();
  ElementAt(0, 0) = sx;
  ElementAt(1, 1) = sy;
  ElementAt(2, 2) = sz;
  ElementAt(3, 3) = 1;
  return *this;
}

template <class T> Matrix44<T> &Matrix44<T>::SetTranslate(const Point3<T> &t) {
  SetTranslate(t[0], t[1], t[2]);
  return *this;
}
template <class T> Matrix44<T> &Matrix44<T>::SetTranslate(const T sx, const T sy, const T sz) {
  SetIdentity();
	ElementAt(0, 3) = sx;
  ElementAt(1, 3) = sy;
  ElementAt(2, 3) = sz;
  return *this;
}
template <class T> Matrix44<T> &Matrix44<T>::SetRotate(T AngleRad, const Point3<T> & axis) {  
  //angle = angle*(T)3.14159265358979323846/180; e' in radianti!
  T c = math::Cos(AngleRad);
  T s = math::Sin(AngleRad);
	T q = 1-c;  
	Point3<T> t = axis;
	t.Normalize();
	ElementAt(0,0) = t[0]*t[0]*q + c;
	ElementAt(0,1) = t[0]*t[1]*q - t[2]*s;
	ElementAt(0,2) = t[0]*t[2]*q + t[1]*s;
	ElementAt(0,3) = 0;
	ElementAt(1,0) = t[1]*t[0]*q + t[2]*s;
	ElementAt(1,1) = t[1]*t[1]*q + c;
	ElementAt(1,2) = t[1]*t[2]*q - t[0]*s;
	ElementAt(1,3) = 0;
	ElementAt(2,0) = t[2]*t[0]*q -t[1]*s;
	ElementAt(2,1) = t[2]*t[1]*q +t[0]*s;
	ElementAt(2,2) = t[2]*t[2]*q +c;
	ElementAt(2,3) = 0;
	ElementAt(3,0) = 0;
	ElementAt(3,1) = 0;									
	ElementAt(3,2) = 0;
	ElementAt(3,3) = 1;	
  return *this;
}


template <class T> T Matrix44<T>::Determinant() const {  
  LinearSolve<T> solve(*this);
  return solve.Determinant();
}


template <class T> Point3<T> operator*(const Matrix44<T> &m, const Point3<T> &p) {
  T w;
  Point3<T> s;
  s[0] = m.ElementAt(0, 0)*p[0] + m.ElementAt(0, 1)*p[1] + m.ElementAt(0, 2)*p[2] + m.ElementAt(0, 3);
  s[1] = m.ElementAt(1, 0)*p[0] + m.ElementAt(1, 1)*p[1] + m.ElementAt(1, 2)*p[2] + m.ElementAt(1, 3);
  s[2] = m.ElementAt(2, 0)*p[0] + m.ElementAt(2, 1)*p[1] + m.ElementAt(2, 2)*p[2] + m.ElementAt(2, 3);
     w = m.ElementAt(3, 0)*p[0] + m.ElementAt(3, 1)*p[1] + m.ElementAt(3, 2)*p[2] + m.ElementAt(3, 3);
	if(w!= 0) s /= w;
  return s;
}

//template <class T> Point3<T> operator*(const Point3<T> &p, const Matrix44<T> &m) {
//  T w;
//  Point3<T> s;
//  s[0] = m.ElementAt(0, 0)*p[0] + m.ElementAt(1, 0)*p[1] + m.ElementAt(2, 0)*p[2] + m.ElementAt(3, 0);
//  s[1] = m.ElementAt(0, 1)*p[0] + m.ElementAt(1, 1)*p[1] + m.ElementAt(2, 1)*p[2] + m.ElementAt(3, 1);
//  s[2] = m.ElementAt(0, 2)*p[0] + m.ElementAt(1, 2)*p[1] + m.ElementAt(2, 2)*p[2] + m.ElementAt(3, 2);
//  w    = m.ElementAt(0, 3)*p[0] + m.ElementAt(1, 3)*p[1] + m.ElementAt(2, 3)*p[2] + m.ElementAt(3, 3);
//	if(w != 0) s /= w;
//  return s;
//}

template <class T> Matrix44<T> &Transpose(Matrix44<T> &m) {
  for(int i = 1; i < 4; i++)
    for(int j = 0; j < i; j++) {
      T t = m.ElementAt(i, j); 
      m.ElementAt(i, j) = m.ElementAt(j, i);
      m.ElementAt(j, i) = t;
    }
  return m;
}




template <class T> Matrix44<T> &Invert(Matrix44<T> &m) {        
  LinearSolve<T> solve(m);

  for(int j = 0; j < 4; j++) { //Find inverse by columns.
    Point4<T> col(0, 0, 0, 0);
    col[j] = 1.0;
    col = solve.Solve(col);
    for(int i = 0; i < 4; i++) 
      m.ElementAt(i, j) = col[i];
  }  
  return m;
}

template <class T> Matrix44<T> Inverse(const Matrix44<T> &m) {        
  LinearSolve<T> solve(m);
  Matrix44<T> res;
  for(int j = 0; j < 4; j++) { //Find inverse by columns.
    Point4<T> col(0, 0, 0, 0);
    col[j] = 1.0;
    col = solve.Solve(col);
    for(int i = 0; i < 4; i++) 
      res.ElementAt(i, j) = col[i];
  }  
  return res;
}



/* LINEAR SOLVE IMPLEMENTATION */

template <class T> LinearSolve<T>::LinearSolve(const Matrix44<T> &m): Matrix44<T>(m) {
  if(!Decompose()) {
    for(int i = 0; i < 4; i++)
      index[i] = i;
    SetZero();
  }
}


template <class T> T LinearSolve<T>::Determinant() const {
  T det = d;
  for(int j = 0; j < 4; j++) 
    det *= ElementAt(j, j);   
  return det;
}


/*replaces a matrix by its LU decomposition of a rowwise permutation.
d is +or -1 depeneing of row permutation even or odd.*/
#define TINY 1e-100

template <class T> bool LinearSolve<T>::Decompose() {  
  
 /* Matrix44<T> A;
  for(int i = 0; i < 16; i++)
    A[i] = operator[](i);  
  SetIdentity();    
  Point4<T> scale;
  // Set scale factor, scale(i) = max( |a(i,j)| ), for each row
  for(int i = 0; i < 4; i++ ) {
    index[i] = i;			  // Initialize row index list
    T scalemax = (T)0.0;
    for(int j = 0; j < 4; j++) 
      scalemax = (scalemax > math::Abs(A.ElementAt(i,j))) ? scalemax : math::Abs(A.ElementAt(i,j));
    scale[i] = scalemax;
  }

  // Loop over rows k = 1, ..., (N-1)
  int signDet = 1;
  for(int k = 0; k < 3; k++ ) {
	  // Select pivot row from max( |a(j,k)/s(j)| )
    T ratiomax = (T)0.0;
	  int jPivot = k;
    for(int i = k; i < 4; i++ ) {
      T ratio = math::Abs(A.ElementAt(index[i], k))/scale[index[i]];
      if(ratio > ratiomax) {
        jPivot = i;
        ratiomax = ratio;
      }
    }
	  // Perform pivoting using row index list
	  int indexJ = index[k];
	  if( jPivot != k ) {	          // Pivot
      indexJ = index[jPivot];
      index[jPivot] = index[k];   // Swap index jPivot and k
      index[k] = indexJ;
	    signDet *= -1;			  // Flip sign of determinant
	  }
	  // Perform forward elimination
    for(int i=k+1; i < 4; i++ ) {
      T coeff = A.ElementAt(index[i],k)/A.ElementAt(indexJ,k);
      for(int j=k+1; j < 4; j++ )
        A.ElementAt(index[i],j) -= coeff*A.ElementAt(indexJ,j);
      A.ElementAt(index[i],k) = coeff;
      //for( j=1; j< 4; j++ ) 
      //  ElementAt(index[i],j) -= A.ElementAt(index[i], k)*ElementAt(indexJ, j);
    }
  }
  for(int i = 0; i < 16; i++)
    operator[](i) = A[i];

  d = signDet; 
  // this = A;
  return true;  */

  d = 1; //no permutation still
    
  T scaling[4];
  int i,j,k;
  //Saving the scvaling information per row
  for(i = 0; i < 4; i++) { 
    T largest = 0.0;
    for(j = 0; j < 4; j++) {
      T t = math::Abs(ElementAt(i, j));
      if (t > largest) largest = t;
    }

    if (largest == 0.0) { //oooppps there is a zero row!
      return false;
    }    
    scaling[i] = (T)1.0 / largest; 
  }

  int imax;
  for(j = 0; j < 4; j++) { 
    for(i = 0; i < j; i++) {
      T sum = ElementAt(i,j);
      for(int k = 0; k < i; k++) 
        sum -= ElementAt(i,k)*ElementAt(k,j);
      ElementAt(i,j) = sum;
    }
    T largest = 0.0; 
    for(i = j; i < 4; i++) { 
      T sum = ElementAt(i,j);
      for(k = 0; k < j; k++)
        sum -= ElementAt(i,k)*ElementAt(k,j);
      ElementAt(i,j) = sum;
      T t = scaling[i] * math::Abs(sum);
      if(t >= largest) { 
        largest = t;
        imax = i;
      }
    }
    if (j != imax) { 
      for (int k = 0; k < 4; k++) { 
        T dum = ElementAt(imax,k);
        ElementAt(imax,k) = ElementAt(j,k);
        ElementAt(j,k) = dum;
      }
      d = -(d);
      scaling[imax] = scaling[j]; 
    }
    index[j]=imax;
    if (ElementAt(j,j) == 0.0) ElementAt(j,j) = (T)TINY;
    if (j != 3) { 
      T dum = (T)1.0 / (ElementAt(j,j));
      for (i = j+1; i < 4; i++) 
        ElementAt(i,j) *= dum;
    }
  }
  return true; 
}


template <class T> Point4<T> LinearSolve<T>::Solve(const Point4<T> &b) { 
  Point4<T> x(b);
  int first = -1, ip;  
  for(int i = 0; i < 4; i++) { 
    ip = index[i];
    T sum = x[ip];
    x[ip] = x[i];
    if(first!= -1)
      for(int j = first; j <= i-1; j++) 
        sum -= ElementAt(i,j) * x[j];
    else 
      if(sum) first = i; 
    x[i] = sum;
  }
  for (int i = 3; i >= 0; i--) { 
    T sum = x[i];
    for (int j = i+1; j < 4; j++) 
      sum -= ElementAt(i, j) * x[j];
    x[i] = sum / ElementAt(i, i); 
  }
  return x;
}

} //namespace
#endif


