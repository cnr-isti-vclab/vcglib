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
Revision 1.2  2004/07/13 06:48:26  cignoni
removed uppercase references in include

Revision 1.1  2004/05/28 13:09:05  ganovelli
created

Revision 1.1  2004/05/28 13:00:39  ganovelli
created


****************************************************************************/


#ifndef __VCGLIB_MATRIX33_H
#define __VCGLIB_MATRIX33_H

#include <stdio.h>
#include <vcg/space/point3.h>
#include <vector>

namespace vcg {

template<class S>
/** @name Matrix33
	Class Matrix33.
    This is the class for definition of a matrix 3x3.	
	@param S (Templete Parameter) Specifies the ScalarType field.
*/
class Matrix33
{
public:
	typedef S ScalarType;
	
	/// Default constructor
	inline Matrix33() {}

	/// Copy constructor
	Matrix33( const Matrix33 & m )
	{
		for(int i=0;i<9;++i)
			a[i] = m.a[i];
	}

	/// create from array
	Matrix33( const S * v )
	{
		for(int i=0;i<9;++i) a[i] = v[i];
	}

	///	Number of columns
	inline unsigned int ColumnsNumber() const
	{
		return 3;
	};

	/// Number of rows
	inline unsigned int RowsNumber() const
	{
		return 3;
	};

	/// Assignment operator
	Matrix33 & operator = ( const Matrix33 & m )
	{
		for(int i=0;i<9;++i)
			a[i] = m.a[i];
		return *this;
	}


	
	/// Operatore di indicizzazione
	inline S * operator [] ( const int i )
	{
		return a+i*3;
	}
	/// Operatore const di indicizzazione
	inline const S * operator [] ( const int i ) const
	{
		return a+i*3;
	}
	

	/// Modificatore somma per matrici 3x3
	Matrix33 & operator += ( const Matrix33  &m )
	{
		for(int i=0;i<9;++i)
			a[i] += m.a[i];
		return *this;
	}

	/// Modificatore sottrazione per matrici 3x3
	Matrix33 & operator -= ( const Matrix33 &m )
	{
		for(int i=0;i<9;++i)
			a[i] -= m.a[i];
		return *this;
	}

	/// Modificatore divisione per scalare
	Matrix33 & operator /= ( const S &s )
	{
		for(int i=0;i<9;++i)
			a[i] /= s;
		return *this;
	}


	/// Modificatore prodotto per matrice
	Matrix33 operator * ( const Matrix33< S> & t ) const
	{
		Matrix33<S> r;

		int i,j;
		for(i=0;i<3;++i)
			for(j=0;j<3;++j)
					r[i][j] = (*this)[i][0]*t[0][j] + (*this)[i][1]*t[1][j] + (*this)[i][2]*t[2][j]; 

		return r;
	}

	/// Modificatore prodotto per costante
	Matrix33 & operator *= ( const S t )
	{
		for(int i=0;i<9;++i)
			a[i] *= t;
		return *this;
	}

	/// Operatore prodotto per costante
	Matrix33 operator * ( const S t )
	{
		Matrix33<S> r;
		for(int i=0;i<9;++i)
			r.a[i] = a[i]* t;

		return r;
	}

	/// Operatore sottrazione per matrici 3x3
	Matrix33  operator - ( const Matrix33 &m )
	{
		Matrix33<S> r;
		for(int i=0;i<9;++i)
			r.a[i] = a[i] - m.a[i];

		return r;
	}

	/** Operatore per il prodotto matrice-vettore.
		@param v A point in $R^{3}$
		@return Il vettore risultante in $R^{3}$
	*/
	Point3<S> operator * ( const Point3<S> & v ) const
	{
		Point3<S> t;

		t[0] = a[0]*v[0] + a[1]*v[1] + a[2]*v[2];
		t[1] = a[3]*v[0] + a[4]*v[1] + a[5]*v[2];
		t[2] = a[6]*v[0] + a[7]*v[1] + a[8]*v[2];
		return t;
	}

	void OuterProduct(Point3<S> const &p0, Point3<S> const &p1) {
		Point3<S> row;
		row = p1*p0[0];
		a[0] = row[0];a[1] = row[1];a[2] = row[2];
		row = p1*p0[1];
		a[3] = row[0]; a[4] = row[1]; a[5] = row[2];
		row = p1*p0[2];
		a[6] = row[0];a[7] = row[1];a[8] = row[2];
	}

	void SetZero()	{
		for(int i=0;i<9;++i) a[i] =0;
	}
	void SetIdentity()	{
		for(int i=0;i<9;++i) a[i] =0;
		a[0]=a[4]=a[8]=1.0;
	}

	void Rotate(S angle, const Point3<S> & axis )
	{
		angle = angle*3.14159265358979323846/180;
		double c = cos(angle);
		double s = sin(angle);
		double q = 1-c;
		Point3<S> t = axis;
		t.Normalize();
		a[0] = t[0]*t[0]*q + c;
		a[1] = t[0]*t[1]*q - t[2]*s;
		a[2] = t[0]*t[2]*q + t[1]*s;
		a[3] = t[1]*t[0]*q + t[2]*s;
		a[4] = t[1]*t[1]*q + c;
		a[5] = t[1]*t[2]*q - t[0]*s;
		a[6] = t[2]*t[0]*q -t[1]*s;
		a[7] = t[2]*t[1]*q +t[0]*s;
		a[8] = t[2]*t[2]*q +c;
	}
	/// Funzione per eseguire la trasposta della matrice
	Matrix33 & Transpose()
	{
		swap(a[1],a[3]);
		swap(a[2],a[6]);
		swap(a[5],a[7]);
		return *this;
	}

	/// Funzione per costruire una matrice diagonale dati i tre elem.
	Matrix33 & SetDiagonal(S *v)
	{int i,j;
		for(i=0;i<3;i++)
      for(j=0;j<3;j++)
        if(i==j) (*this)[i][j] = v[i];
        else     (*this)[i][j] = 0;
    return *this;
	}


	/// Assegna l'n-simo vettore colonna
	void SetColumn(const int n, S* v){
		assert( (n>=0) && (n<3) );
		a[n]=v[0]; a[n+3]=v[1]; a[n+6]=v[2];
	};

	/// Assegna l'n-simo vettore riga
	void SetRow(const int n, S* v){
		assert( (n>=0) && (n<3) );
		int m=n*3;
		a[m]=v[0]; a[m+1]=v[1]; a[m+2]=v[2];
	};

	/// Assegna l'n-simo vettore colonna
	void SetColumn(const int n, const Point3<S> v){
		assert( (n>=0) && (n<3) );
		a[n]=v[0]; a[n+3]=v[1]; a[n+6]=v[2];
	};

	/// Assegna l'n-simo vettore riga
	void SetRow(const int n, const Point3<S> v){
		assert( (n>=0) && (n<3) );
		int m=n*3;
		a[m]=v[0]; a[m+1]=v[1]; a[m+2]=v[2];
	};

	/// Restituisce l'n-simo vettore colonna
	Point3<S> GetColumn(const int n) const {
		assert( (n>=0) && (n<3) );
		Point3<S> t;
		t[0]=a[n]; t[1]=a[n+3]; t[2]=a[n+6];
		return t;
	};

	/// Restituisce l'n-simo vettore riga
	Point3<S> GetRow(const int n) const {
		assert( (n>=0) && (n<3) );
		Point3<S> t;
		int m=n*3;
		t[0]=a[m]; t[1]=a[m+1]; t[2]=a[m+2];
		return t;
	};



	/// Funzione per il calcolo del determinante
	S Determinant() const
	{
		return a[0]*(a[4]*a[8]-a[5]*a[7]) -
			     a[1]*(a[3]*a[8]-a[5]*a[6]) +
					 a[2]*(a[3]*a[7]-a[4]*a[6]) ;
	}

	Matrix33 & Invert()
	{
			// Maple produsse:
		S t4  = a[0]*a[4];
		S t6  = a[0]*a[5];
		S t8  = a[1]*a[3];
		S t10 = a[2]*a[3];
		S t12 = a[1]*a[6];
		S t14 = a[2]*a[6];
		S t17 = 1/(t4*a[8]-t6*a[7]-t8*a[8]+t10*a[7]+t12*a[5]-t14*a[4]);
		S a0  = a[0];
		S a1  = a[1];
		S a3  = a[3];
		S a4  = a[4];
		a[0]  =  (a[4]*a[8]-a[5]*a[7])*t17;
		a[1]  = -(a[1]*a[8]-a[2]*a[7])*t17;
		a[2]  =  (a1  *a[5]-a[2]*a[4])*t17;
		a[3]  = -(a[3]*a[8]-a[5]*a[6])*t17;
		a[4]  =  (a0  *a[8]-t14      )*t17;
		a[5]  = -(t6 - t10)*t17;
		a[6]  =  (a3  *a[7]-a[4]*a[6])*t17;
		a[7]  = -(a[0]*a[7]-t12)*t17;
		a[8]  =  (t4-t8)*t17;

		return *this;
	}

	void show(FILE * fp)
	{
		for(int i=0;i<3;++i)
		    printf("| %g \t%g \t%g |\n",a[3*i+0],a[3*i+1],a[3*i+2]);
	}

// return the Trace of the matrix i.e. the sum of the diagonal elements
S Trace() const
{
	return a[0]+a[4]+a[8];
}

/* 
compute the matrix generated by the product of a * b^T
*/
void ExternalProduct(const Point3<S> &a, const Point3<S> &b) 
{
	for(int i=0;i<3;++i)
		for(int j=0;j<3;++j)
			 (*this)[i][j] = a[i]*b[j];
}

/* 
It compute the cross covariance matrix of two set of 3d points P and X;
it returns also the barycenters of P and X.
fonte:

Besl, McKay
A method for registration o f 3d Shapes 
IEEE TPAMI Vol 14, No 2 1992

*/ 
template <class STLPOINTCONTAINER >
void CrossCovariance(const STLPOINTCONTAINER &P, const STLPOINTCONTAINER &X, 
										 Point3<S> &bp, Point3<S> &bx) 
{
	Zero();
	assert(P.size()==X.size());
	bx.Zero();
	bp.Zero();
	Matrix33<S> tmp;
	typename std::vector <Point3<S> >::const_iterator pi,xi;
	for(pi=P.begin(),xi=X.begin();pi!=P.end();++pi,++xi){
		bp+=*pi;
		bx+=*xi;
		tmp.ExternalProduct(*pi,*xi);
		(*this)+=tmp;
	}
	bp/=P.size();
	bx/=X.size();
	(*this)/=P.size();
	tmp.ExternalProduct(bp,bx);
	(*this)-=tmp;
}

template <class STLPOINTCONTAINER, class STLREALCONTAINER>
void WeightedCrossCovariance(const STLREALCONTAINER &  weights,
							 const STLPOINTCONTAINER &P, 
							 const STLPOINTCONTAINER &X, 
							 Point3<S> &bp, 
							 Point3<S> &bx) 
{
	Zero();
	assert(P.size()==X.size());
	bx.Zero();
	bp.Zero();
	Matrix33<S> tmp;
	typename std::vector <Point3<S> >::const_iterator pi,xi;
	typename STLREALCONTAINER::const_iterator pw;

	for(pi=P.begin(),xi=X.begin();pi!=P.end();++pi,++xi){
		bp+=(*pi);
		bx+=(*xi);
	}
	bp/=P.size();
	bx/=X.size();

	for(pi=P.begin(),xi=X.begin(),pw = weights.begin();pi!=P.end();++pi,++xi,++pw){
		
		tmp.ExternalProduct(((*pi)-(bp)),((*xi)-(bp)));

		(*this)+=tmp*(*pw);
	}
}

private:
	S a[9];
};


/// 
typedef Matrix33<short>  Matrix33s;
typedef Matrix33<int>	 Matrix33i;
typedef Matrix33<float>  Matrix33f;
typedef Matrix33<double> Matrix33d;

} // end of namespace

#endif
