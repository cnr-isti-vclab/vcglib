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

$LOG$

****************************************************************************/

#ifndef __VCGLIB_MATRIX44
#define __VCGLIB_MATRIX44

#include <string.h>
#include <vcg/space/point3.h>
#include <vcg/space/point4.h>


namespace vcg {

  /*
	Annotations:
Opengl stores matrix in  column-major order. That is, the matrix is stored as:

	a0  a4  a8  a12
	a1  a5  a9  a13
	a2  a6  a10 a14
	a3  a7  a11 a15

e.g. glTranslate generate a matrix that is

	1 0 0 tx
	0 1 0 ty 
	0 0 1 tz
	0 0 0  1

Matrix44 stores matrix in row-major order i.e.

	a0  a1  a2  a3
	a4  a5  a6  a7
	a8  a9  a10 a11
	a12 a13 a14 a15

and the Translate Function generate:

	1  0  0  0
	0  1  0  0
	0  0  1  0
	tx ty tz 1

*/

  /** This class represent a 4x4 matrix. T is the kind of element in the matrix.
    */
template <class T> class Matrix44 {  
protected:
  T _a[16];

public:	
  typedef T scalar;

///@{

  /** $name Contrutors
    *  No automatic casting and default constructor is empty
    */
	Matrix44() {};	
	~Matrix44() {};
  Matrix44(const Matrix44 &m);
  Matrix44(const T v[]);

  T &element(const int row, const int col);
  T element(const int row, const int col) const;
  T &operator[](const int i);
  T operator[](const int i) const;

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

	void SetZero();
  void SetIdentity();
  void SetDiagonal(const T k);
  void SetScale(const T sx, const T sy, const T sz);	
	void SetTranslate(const Point3<T> &t);
	void SetTranslate(const T sx, const T sy, const T sz);
  ///use radiants for angle.
  void SetRotate(T angle, const Point3<T> & axis); 

  T Determinant() const;
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
  void Decompose();
};

///	Apply   POST moltiplica la matrice al vettore (e.g. la traslazione deve stare nell'ultima riga)
///	Project PRE  moltiplica la matrice al vettore (e.g. la traslazione deve stare nell'ultima colonna)

/*** Postmultiply (old Apply in the old interface). 
  * SetTranslate, SetScale, SetRotate prepare the matrix for this.
  */
template <class T> Point4<T> operator*(const Point4<T> &p, const Matrix44<T> &m);

///Premultiply  (old Project in the old interface)
template <class T> Point4<T> operator*(const Matrix44<T> &m, const Point4<T> &p);

template <class T> Matrix44<T> &Transpose(Matrix44<T> &m);
template <class T> Matrix44<T> &Invert(Matrix44<T> &m);

typedef Matrix44<short>  Matrix44s;
typedef Matrix44<int>	 Matrix44i;
typedef Matrix44<float>  Matrix44f;
typedef Matrix44<double> Matrix44d;



template <class T> Matrix44<T>::Matrix44(const Matrix44<T> &m) {
  memcpy((T *)_a, (T *)m._a, 16 * sizeof(T));   
}

template <class T> T &Matrix44<T>::element(const int row, const int col) {
  assert(row >= 0 && row < 4);
  assert(col >= 0 && col < 4);
  return _a[(row<<2) + col];
}

template <class T> T Matrix44<T>::element(const int row, const int col) const {
  assert(row >= 0 && row < 4);
  assert(col >= 0 && col < 4);
  return _a[(row<<2) + col];
}

template <class T> T &Matrix44<T>::operator[](const int i) {
  assert(i >= 0 && i < 16);
  return ((T *)_a)[i];
}

template <class T> T Matrix44<T>::operator[](const int i) const {
  assert(i >= 0 && i < 16);
  return ((T *)_a)[i];
}

template <class T> Matrix44<T> Matrix44<T>::operator+(const Matrix44 &m) const {
  Matrix44<T> ret;
  for(int i = 0; i < 16; i++) 
    ret[i] = operator[](i) + m[i];  
}

template <class T> Matrix44<T> Matrix44<T>::operator-(const Matrix44 &m) const {
  Matrix44<T> ret;
  for(int i = 0; i < 16; i++) 
    ret[i] = operator[](i) - m[i];  
}

template <class T> Matrix44<T> Matrix44<T>::operator*(const Matrix44 &m) const {
  Matrix44 ret;       
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++) {
      T t = 0.0;
      for(int k = 0; k < 4; k++)
        t += element[i][k] * m.element[k][j];
      ret.element[i][j] = t;
    }
  return ret;
}

template <class T> Point4<T> Matrix44<T>::operator*(const Point4<T> &v) const {
  Point4<T> ret;     
  for(int i = 0; i < 4; i++){
    T t = 0.0;
    for(int k = 0; k < 4; k++)
      t += element[i][k] * v[k];
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
    res[i] = -operator[](i);
  return res;
}

template <class T> Matrix44<T> Matrix44<T>::operator*(const T k) const {
  Matrix44<T> res;
  for(int i = 0; i < 16; i++)
    res[i] = operator[](i) * k;
  return res;
}

template <class T> void Matrix44<T>::operator+=(const Matrix44 &m) {
  for(int i = 0; i < 16; i++)
    operator[](i) += m[i];  
}
template <class T> void Matrix44<T>::operator-=(const Matrix44 &m) {
  for(int i = 0; i < 16; i++)
    operator[](i) -= m[i];
}
template <class T> void Matrix44<T>::operator*=( const Matrix44 & m ) {
  for(int i = 0; i < 4; i++) {
    Point4<T> t(0, 0, 0, 0);
    for(int k = 0; k < 4; k++) {
      for(int j = 0; j < 4; j++) {
        t[k] += element(i, k) * m.element(k, j);
      }
    }
    for(int l = 0; l < 4; l++)
      element(i, l) = t[l];
  }
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
  element(0, 0) = k;
  element(1, 1) = k;
  element(2, 2) = k;
  element(3, 3) = 1;    
}

template <class T> void Matrix44<T>::SetScale(const T sx, const T sy, const T sz) {
  SetZero();
  element(0, 0) = sx;
  element(1, 1) = sy;
  element(2, 2) = sz;
  element(3, 3) = 1;
}

template <class T> void Matrix44<T>::SetTranslate(const Point3<T> &t) {
  SetTranslate(t[0], t[1], t[2]);
}
template <class T> void Matrix44<T>::SetTranslate(const T sx, const T sy, const T sz) {
  SetIdentity();
	element(3, 0) = sx;
  element(3, 1) = sy;
  element(3, 2) = sz;
}
template <class T> void Matrix44<T>::SetRotate(T angle, const Point3<T> & axis) {  
  //angle = angle*(T)3.14159265358979323846/180; e' in radianti!
  T c = Math<T>::Cos(angle);
  T s = Math<T>::Sin(angle);
	T q = 1-c;  
	Point3<T> t = axis;
	t.Normalize();
	element(0,0) = t[0]*t[0]*q + c;
	element(1,0) = t[0]*t[1]*q - t[2]*s;
	element(2,0) = t[0]*t[2]*q + t[1]*s;
	element(3,0) = 0;
	element(0,1) = t[1]*t[0]*q + t[2]*s;
	element(1,1) = t[1]*t[1]*q + c;
	element(2,1) = t[1]*t[2]*q - t[0]*s;
	element(3,1) = 0;
	element(0,2) = t[2]*t[0]*q -t[1]*s;
	element(1,2) = t[2]*t[1]*q +t[0]*s;
	element(2,2) = t[2]*t[2]*q +c;
	element(3,2) = 0;
	element(0,3) = 0;
	element(1,3) = 0;									
	element(2,3) = 0;
	element(3,3) = 1;	
}


template <class T> T Matrix44<T>::Determinant() const {  
  LinearSolve<T> solve(*this);
  return solve.Determinant();
}


template <class T> Point4<T> operator*(const Matrix44<T> &m, const Point4<T> &p) {
  T w;
  Point4<T> s;
  s.x() = a[0][0]*p.x() + a[0][1]*p.y() + a[0][2]*p.z() + a[0][3];
  s.y() = a[1][0]*p.x() + a[1][1]*p.y() + a[1][2]*p.z() + a[1][3];
  s.z() = a[2][0]*p.x() + a[2][1]*p.y() + a[2][2]*p.z() + a[2][3];
  s.w() = a[3][0]*p.x() + a[3][1]*p.y() + a[3][2]*p.z() + a[3][3];
	if(s.w()!= 0) s /= s.w();
  return s;
}

template <class T> Point4<T> operator*(const Point4<T> &p, const Matrix44<T> &m) {
  T w;
  Point4<T> s;
  s.x() = a[0][0]*p.x() + a[0][1]*p.y() + a[0][2]*p.z() + a[0][3];
  s.y() = a[1][0]*p.x() + a[1][1]*p.y() + a[1][2]*p.z() + a[1][3];
  s.z() = a[2][0]*p.x() + a[2][1]*p.y() + a[2][2]*p.z() + a[2][3];
  s.w() = a[3][0]*p.x() + a[3][1]*p.y() + a[3][2]*p.z() + a[3][3];
	if(s.w() != 0) s /= s.w();
  return s;
}

template <class T> Matrix44<T> &Transpose(Matrix44<T> &m) {
  for(int i = 1; i < 4; i++)
    for(int j = 0; j < i; j++) {
      T t = element(i, j); 
      element(i, j) = element(j, i);
      element(j, i) = t;
    }
  return *this;
}




template <class T> Matrix44<T> &Invert(Matrix44<T> &m) {        
  LinearSolve<T> solve(m);

  for(int j = 0; j < 4; j++) { //Find inverse by columns.
    Point4<T> col(0, 0, 0, 0);
    col[j] = 1.0;
    col = solve.Solve(col);
    for(int i = 0; i < 4; i++) 
      m.element(i, j) = col[i];
  }  
}



/* LINEAR SOLVE IMPLEMENTATION */

template <class T> LinearSolve<T>::LinearSolve(const Matrix44<T> &m): Matrix44<T>(m) {
  Decompose();  
}


template <class T> T LinearSolve<T>::Determinant() const {
  T det = d;
  for(int j = 0; j < 4; j++) 
    det *= element(j, j);   
  return det;
}


/*replaces a matrix by its LU decomposition of a rowwise permutation.
d is +or -1 depeneing of row permutation even or odd.*/
#define TINY 1e-100

template <class T> void LinearSolve<T>::Decompose() {
  d = 1; //no permutation still
    
  T scaling[4];
  
  //Saving the scvaling information per row
  for(int i = 0; i < 4; i++) { 
    T largest = 0.0;
    for(int j = 0; j < 4; j++) {
      T t = fabs(element(i, j));
      if (t > largest) largest = t;
    }

    if (largest == 0.0) { //oooppps there is a zero row!
      return;
    }    
    scaling[i] = 1.0 / largest; 
  }

  int imax;
  for(int j = 0; j < 4; j++) { 
    for(int i = 0; i < j; i++) {
      T sum = element(i,j);
      for(int k = 0; k < i; k++) 
        sum -= element(i,k)*element(k,j);
      element(i,j) = sum;
    }
    T largest = 0.0; 
    for(int i = j; i < 4; i++) { 
      T sum = element(i,j);
      for(int k = 0; k < j; k++)
        sum -= element(i,k)*element(k,j);
      element(i,j) = sum;
      T t = scaling[i] * fabs(sum);
      if(t >= largest) { 
        largest = t;
        imax = i;
      }
    }
    if (j != imax) { 
      for (int k = 0; k < 4; k++) { 
        T dum = element(imax,k);
        element(imax,k) = element(j,k);
        element(j,k) = dum;
      }
      d = -(d);
      scaling[imax] = scaling[j]; 
    }
    index[j]=imax;
    if (element(j,j) == 0.0) element(j,j) = TINY;
    if (j != 3) { 
      T dum = 1.0 / (element(j,j));
      for (i = j+1; i < 4; i++) 
        element(i,j) *= dum;
    }
  }
}


template <class T> Point4<T> LinearSolve<T>::Solve(const Point4<T> &b) { 
  Point4<T> x(b);
  int first = 0, ip;  
  for(int i = 0; i < 4; i++) { 
    ip = index[i];
    T sum = x[ip];
    x[ip] = x[i];
    if(first)
      for(int j = first; j <= i-1; j++) 
        sum -= element(i,j) * x[j];
    else 
      if (sum) first = i; 
  }
  for (int i = 3; i >= 0; i--) { 
    T sum = x[i];
    for (int j = i+1; j < 4; j++) 
      sum -= element(i, j) * x[j];
    x[i] = sum / element(i, i); 
  }
  return x;
}

} //namespace
#endif




/*#ifdef __GL_H__
// Applicano la trasformazione intesa secondo la Apply
void glMultMatrix(Matrix44<double> const & M)   { 	glMultMatrixd(&(M.a[0][0]));}
void glMultMatrix(Matrix44<float>  const & M)   { 	glMultMatrixf(&(M.a[0][0]));}
void glLoadMatrix(Matrix44<double> const & M)   { 	glLoadMatrixd(&(M.a[0][0]));}
void glLoadMatrix(Matrix44<float>  const & M)   { 	glLoadMatrixf(&(M.a[0][0]));}
#endif*/


