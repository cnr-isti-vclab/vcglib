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

#warning You are including deprecated math stuff
/*!
*	\deprecated use cols()
*	Number of columns
*/
EIGEN_DEPRECATED inline unsigned int ColumnsNumber() const
{
	return cols();
};


/*!
*	\deprecated use rows()
*	Number of rows
*/
EIGEN_DEPRECATED inline unsigned int RowsNumber() const
{
	return rows();
};

/*
*	\deprecated use *this.isApprox(m) or *this.cwise() == m
*	Equality operator.
*	\param m
*	\return true iff the matrices have same size and its elements have same values.
*/
// template<typename OtherDerived>
// EIGEN_DEPRECATED bool operator==(const MatrixBase<OtherDerived> &m) const
// {
// 	return (this->cwise() == m);
// }

/*
*	\deprecated use !*this.isApprox(m) or *this.cwise() != m
*	Inequality operator
*	\param m
*	\return true iff the matrices have different size or if their elements have different values.
*/
// template<typename OtherDerived>
// EIGEN_DEPRECATED bool operator!=(const MatrixBase<OtherDerived> &m) const
// {
// 	return (this->cwise() != m);
// };

/*!
*	\deprecated use *this(i,j) (or *this.coeff(i,j))
* Return the element stored in the <I>i</I>-th rows at the <I>j</I>-th column
*	\param i the row index
*	\param j the column index
*	\return the element
*/
EIGEN_DEPRECATED inline Scalar ElementAt(unsigned int i, unsigned int j)
{
	return (*this)(i,j);
};

/*!
*	\deprecated use *this.determinant() (or *this.lu().determinant() for large matrices)
*	Calculate and return the matrix determinant (Laplace)
*	\return	the matrix determinant
*/
EIGEN_DEPRECATED Scalar Determinant() const
{
	return determinant();
};

/*!
*	Return the cofactor <I>A<SUB>i,j</SUB></I> of the <I>a<SUB>i,j</SUB></I> element
*	\return	...
*/
EIGEN_DEPRECATED Scalar Cofactor(unsigned int i, unsigned int j) const
{
	assert(rows() == cols());
	assert(rows()>2);
	return (((i+j)%2==0) ? 1. : -1.) * minor(i,j).determinant();
};

/*!
*	\deprecated use *this.col(j)
*	Get the <I>j</I>-th column on the matrix.
*	\param j	the column index.
*	\return		the reference to the column elements. This pointer must be deallocated by the caller.
*/
EIGEN_DEPRECATED ColXpr GetColumn(const unsigned int j)
{
	return col(j);
};

/*!
*	\deprecated use *this.row(i)
*	Get the <I>i</I>-th row on the matrix.
*	\param i	the column index.
*	\return		the reference to the row elements. This pointer must be deallocated by the caller.
*/
EIGEN_DEPRECATED RowXpr GetRow(const unsigned int i)
{
	return row(i);
};

/*!
*	\deprecated use m1.col(i).swap(m1.col(j));
* Swaps the values of the elements between the <I>i</I>-th and the <I>j</I>-th column.
* \param i the index of the first column
* \param j the index of the second column
*/
EIGEN_DEPRECATED void SwapColumns(const unsigned int i, const unsigned int j)
{
	if (i==j)
		return;

	col(i).swap(col(j));
};

/*!
*	\deprecated use m1.col(i).swap(m1.col(j))
* Swaps the values of the elements between the <I>i</I>-th and the <I>j</I>-th row.
* \param i the index of the first row
* \param j the index of the second row
*/
EIGEN_DEPRECATED void SwapRows(const unsigned int i, const unsigned int j)
{
	if (i==j)
		return;

	row(i).swap(row(j));
};


/*!
*	\deprecated use *this.cwise() += k
*	(Modifier) Add to each element of this matrix the scalar constant <I>k</I>.
* \param k	the scalar constant
*	\return		the modified matrix
*/
EIGEN_DEPRECATED Derived& operator+=(const Scalar k)
{
	cwise() += k;
	return *this;
};

/*!
*	\deprecated use *this.cwise() -= k
*	(Modifier) Subtract from each element of this matrix the scalar constant <I>k</I>.
* \param k	the scalar constant
*	\return		the modified matrix
*/
EIGEN_DEPRECATED Derived& operator-=(const Scalar k)
{
	cwise() -= k;
	return *this;
};

/*!
*	\deprecated use *this.dot
*	Matrix multiplication: calculates the cross product.
*	\param	reference to the matrix to multiply by
*	\return the matrix product
*/
// template <int N,int M>
// EIGEN_DEPRECATED void DotProduct(Point<N,Scalar> &m,Point<M,Scalar> &result)
// {
// 	unsigned int i, j;
// 	for (i=0; i<M; i++)
// 	{ result[i]=0;
// 		for (j=0; j<N; j++)
// 			result[i]+=(*this)[i][j]*m[j];
// 	}
// };

/*!
*	\deprecated use (*this) * vec.asDiagonal() or (*this) * mat.mark<Diagonal>()
*	Matrix multiplication by a diagonal matrix
*/
// EIGEN_DEPRECATED Matrix<Scalar> operator*(const MatrixDiagBase &m) const
// {
// 	assert(_columns == _rows);
// 	assert(_columns == m.Dimension());
// 	int i,j;
// 	Matrix<Scalar> result(_rows, _columns);
// 
// 	for (i=0; i<result._rows; i++)
// 		for (j=0; j<result._columns; j++)
// 			result[i][j]*= m[j];
// 
// 	return result;
// };

/*!
*	\deprecated use *this = a * b.transpose()
*	Matrix from outer product.
*/
template <typename OtherDerived1, typename OtherDerived2>
EIGEN_DEPRECATED void OuterProduct(const MatrixBase<OtherDerived1>& a, const MatrixBase<OtherDerived2>& b)
{
	*this = a * b.transpose();
}

typedef CwiseUnaryOp<ei_scalar_add_op<Scalar>, Derived> ScalarAddReturnType;

/*!
*	\deprecated use *this.cwise() + k
*	Scalar sum.
*	\param k
*	\return		the resultant matrix
*/
EIGEN_DEPRECATED const ScalarAddReturnType operator+(const Scalar k) { return cwise() + k; }

/*!
*	\deprecated use *this.cwise() - k
*	Scalar difference.
*	\param k
*	\return		the resultant matrix
*/
EIGEN_DEPRECATED const ScalarAddReturnType operator-(const Scalar k) { return cwise() - k; }


/*!
*	\deprecated use *this.setZero() or *this = MatrixType::Zero(rows,cols), etc.
*	Set all the matrix elements to zero.
*/
EIGEN_DEPRECATED void SetZero()
{
	setZero();
};

/*!
*	\deprecated use *this.setIdentity() or *this = MatrixType::Identity(rows,cols), etc.
*	Set the matrix to identity.
*/
EIGEN_DEPRECATED void SetIdentity()
{
	setIdentity();
};

/*!
*	\deprecated use *this.col(j) = expression
*	Set the values of <I>j</I>-th column to v[j]
*	\param j	the column index
*	\param v	...
*/
EIGEN_DEPRECATED void SetColumn(const unsigned int j, Scalar* v)
{
	col(j) = Map<Matrix<Scalar,RowsAtCompileTime,1> >(v,cols(),1);
};

/*!
*	\deprecated use *this.row(i) = expression
*	Set the elements of the <I>i</I>-th row to v[j]
*	\param i	the row index
*	\param v	...
*/
EIGEN_DEPRECATED void SetRow(const unsigned int i, Scalar* v)
{
	row(i) = Map<Matrix<Scalar,1,ColsAtCompileTime> >(v,1,rows());
};

/*!
*	\deprecated use *this.diagonal() = expression
*	Set the diagonal elements <I>v<SUB>i,i</SUB></I> to v[i]
*	\param v
*/
EIGEN_DEPRECATED void SetDiagonal(Scalar *v)
{
	assert(rows() == cols());
	diagonal() = Map<Matrix<Scalar,RowsAtCompileTime,1> >(v,cols(),1);
}


/*!
*	\deprecated use *this = *this.transpose()
*/
// Transpose already exist
// EIGEN_DEPRECATED void Transpose()
// {
// 	assert(0 && "dangerous use of deprecated Transpose function, please use: m = m.transpose();");
// };


/*!
*	\deprecated use ostream << *this or ostream << *this.withFormat(...)
*	Print all matrix elements
*/
EIGEN_DEPRECATED void Dump()
{
	unsigned int i, j;
	for (i=0; i<rows(); ++i)
	{
		printf("[\t");
		for (j=0; j<cols(); j++)
			printf("%f\t", coeff(i,j));
		printf("]\n");
	}
	printf("\n");
}


