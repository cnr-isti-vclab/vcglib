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

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include <algorithm>


namespace vcg
{
	 /** \addtogroup math */
   /* @{ */

/*!
 * This class represent a generic <I>m</I>×<I>n</I> matrix. The class is templated over the scalar type field.
 * @param TYPE (Templete Parameter) Specifies the ScalarType field.
 */
	template<class TYPE> 
	class Matrix
		{
			
		public:
			typedef TYPE ScalarType;

			/*!
			*	Default constructor
			* All the elements are initialized to zero.
			*	\param m the number of matrix rows
			* \param n the number of matrix columns
			*/
			Matrix(unsigned int m, unsigned int n)
			{
				_rows = m;
				_columns = n;
				_data = new ScalarType[m*n];
				memset(_data, ScalarType(0.0), m*n*sizeof(ScalarType));
			};

			/*!
			*	Constructor
			* The matrix elements are initialized with the values of the elements in \i values.
			*	\param m the number of matrix rows
			* \param n the number of matrix columns
			*	\param values the values of the matrix elements
			*/
			Matrix(unsigned int m, unsigned int n, TYPE *values)
			{
				_rows = m;
				_columns = n;
				_data = new ScalarType[m*n];
				unsigned int i;
        for (i=0; i<_rows*_columns; i++)
					_data[i] = values[i];
			};

			/*!
			*	Copy constructor
			*	The matrix elements are initialized with the value of the corresponding element in \i m
			* \param m the matrix to be copied
			*/
			Matrix(const Matrix<TYPE> &m)
			{
				_rows = m._rows;
				_columns = m._columns;
				_data = new ScalarType[_rows*_columns];
				for (unsigned int i=0; i<_rows*_columns; i++)
					_data[i] = m._data[i];
			};

			/*!
			*	Default destructor
			*/
			~Matrix()
			{
				delete []_data;
			};

			/*!
			*	Number of columns
			*/
			inline unsigned int ColumnsNumber()
			{
				return _columns;
			};


			/*!
			*	Number of rows
			*/
			inline unsigned int RowsNumber()
			{
				return _rows;
			};

						/*!
			*	Equality operator.
			*	\param m
			*	\return true iff the matrices have same size and its elements have same values.
			*/
			bool operator==(const Matrix<TYPE> &m) const
			{
				if (_rows==m._rows && _columns==m._columns)
				{
					bool result = true;
					for (unsigned int i=0; i<_rows*_columns && result; i++)
						result = (_data[i]==m._data[i]);
					return result;
				}
				return false;
			};

			/*!
			*	Inequality operator
			*	\param m
			*	\return true iff the matrices have different size or if their elements have different values.
			*/
			bool operator!=(const Matrix<TYPE> &m) const
			{
				if (_rows==m._rows && _columns==m._columns)
				{
					bool result = false;
					for (unsigned int i=0; i<_rows*_columns && !result; i++)
						result = (_data[i]!=m._data[i]);
					return result;
				}
				return true;
			};

			/*!
			* Return the element stored in the <I>i</I>-th rows at the <I>j</I>-th column
			*	\param i the row index
			*	\param j the column index
			*	\return the element 
			*/
			inline TYPE ElementAt(unsigned int i, unsigned int j)
			{
				assert(i>=0 && i<_rows);
				assert(j>=0 && j<_columns);
				return _data[i*_columns+j];
			};

			/*!
			*	Calculate and return the matrix determinant (Laplace)
			*	\return	the matrix determinant
			*/
			TYPE Determinant() const
			{
				assert(_rows == _columns);
				switch (_rows)
				{
				case 2:
					{
						return _data[0]*_data[3]-_data[1]*_data[2];
						break;
					};
				case 3:
					{
						return	_data[0]*(_data[4]*_data[8]-_data[5]*_data[7]) - 
										_data[1]*(_data[3]*_data[8]-_data[5]*_data[6]) + 
										_data[2]*(_data[3]*_data[7]-_data[4]*_data[6]) ;
						break;
					};
				default:
					{
						// da migliorare: si puo' cercare la riga/colonna con maggior numero di zeri
						ScalarType det = 0;
						for (unsigned int j=0; j<_columns; j++)
							if (_data[j]!=0)
								det += _data[j]*this->Cofactor(0, j);

						return det;
					}
				};
			};

			/*!
			*	Return the cofactor <I>A<SUB>i,j</SUB></I> of the <I>a<SUB>i,j</SUB></I> element
			*	\return	...
			*/
			TYPE Cofactor(unsigned int i, unsigned int j) const
			{
				assert(_rows == _columns);
				assert(_rows>2);
				TYPE *values = new TYPE[(_rows-1)*(_columns-1)];
				unsigned int u, v, p, q, s, t;
				for (u=0, p=0, s=0, t=0; u<_rows; u++, t+=_rows)
				{
					if (i==u)
						continue;

					for (v=0, q=0; v<_columns; v++)
					{
						if (j==v)
							continue;
						values[s+q] = _data[t+v];
						q++;
					}
					p++;
					s+=(_rows-1);
				}
				Matrix<TYPE> temp(_rows-1, _columns-1, values);
				return (pow(-1, i+j)*temp.Determinant());
			};

			/*!
			*	Subscript operator: 
			* \param i	the index of the row
			*	\return a reference to the <I>i</I>-th matrix row
			*/
			inline TYPE* operator[](const unsigned int i)
			{
				assert(i>=0 && i<_columns);
				return _data + i*_columns;
			};

			/*!
			*	Const subscript operator
			* \param i	the index of the row
			*	\return a reference to the <I>i</I>-th matrix row
			*/
			inline const TYPE* operator[](const unsigned int i) const 
			{
				assert(i>=0 && i<_columns);
				return _data + i*_columns;
			};

						/*!
			*	Get the <I>j</I>-th column on the matrix.
			*	\param j	the column index.
			*	\return		the reference to the column elements.
			*/
			TYPE* GetColumn(const unsigned int j)
			{
				assert(j>=0 && j<_columns);
				ScalarType *v = new ScalarType[_columns];
				unsigned int i, p;
				for (i=0, p=j; i<_rows; i++, p+=_columns)
					v[i] = _data[p];
				return v;
			};

			/*!
			*	Get the <I>i</I>-th row on the matrix.
			*	\param i	the column index.
			*	\return		the reference to the row elements.
			*/
			TYPE* GetRow(const unsigned int i)
			{
				assert(i>=0 && i<_rows);
				ScalarType *v = new ScalarType[_rows];
				unsigned int j, p;
				for (j=0, p=i*_columns; j<_columns; j++, p++)
					v[j] = _data[p];
				return v;
			};

			/*!
			*	Assignment operator
			*	\param m ...
			*/
			Matrix<TYPE>& operator=(const Matrix<TYPE> &m)
			{
				if (this != &m)
				{
					assert(_rows == m._rows);
					assert(_columns == m._columns);
					for (unsigned int i=0; i<_rows*_columns; i++)
						_data[i] = m._data[i];
				}
				return *this;
			};

			/*!
			*	Adds a matrix <I>m</I> to this matrix.
			*	\param m  reference to matrix to add to this  
			*	\return		the matrix sum. 
			*/
			Matrix<TYPE>& operator+=(const Matrix<TYPE> &m)
			{
				assert(_rows == m._rows);
				assert(_columns == m._columns);
				for (unsigned int i=0; i<_rows*_columns; i++)
					_data[i] += m._data[i];
				return *this;
			};

			/*!
			*	Subtracts a matrix <I>m</I> to this matrix. 
			*	\param m  reference to matrix to subtract 
			*	\return		the matrix difference. 
			*/
			Matrix<TYPE>& operator-=(const Matrix<TYPE> &m)
			{
				assert(_rows == m._rows);
				assert(_columns == m._columns);
				for (unsigned int i=0; i<_rows*_columns; i++)
					_data[i] -= m._data[i];
				return *this;
			};

			/*!
			*	(Modifier) Add to each element of this matrix the scalar constant <I>k</I>. 
			* \param k	the scalar constant
			*	\return		the modified matrix
			*/
			Matrix<TYPE>& operator+=(const TYPE k)
			{
				for (unsigned int i=0; i<_rows*_columns; i++)
					_data[i] += k;
				return *this;
			};

			/*!
			*	(Modifier) Subtract from each element of this matrix the scalar constant <I>k</I>. 
			* \param k	the scalar constant
			*	\return		the modified matrix
			*/
			Matrix<TYPE>& operator-=(const TYPE k)
			{
				for (unsigned int i=0; i<_rows*_columns; i++)
					_data[i] -= k;
				return *this;
			};

			/*!
			*	(Modifier) Multiplies each element of this matrix by the scalar constant <I>k</I>. 
			* \param k	the scalar constant
			*	\return		the modified matrix
			*/
			Matrix<TYPE>& operator*=(const TYPE k)
			{
				for (unsigned int i=0; i<_rows*_columns; i++)
					_data[i] *= k;
				return *this;
			};

			/*!
			*	(Modifier) Divides each element of this matrix by the scalar constant <I>k</I>. 
			* \param k	the scalar constant
			*	\return		the modified matrix
			*/
			Matrix<TYPE>& operator/=(const TYPE k)
			{
				assert(k!=0);
				for (unsigned int i=0; i<_rows*_columns; i++)
					_data[i] /= k;
				return *this;
			};

			/*!
			*	Matrix multiplication: calculates the cross product.
			*	\param	reference to the matrix to multiply by 
			*	\result the matrix product
			*/
			Matrix<TYPE> operator*(const Matrix<TYPE> &m)
			{
				assert(_columns == m._rows);
				Matrix<TYPE> result(_rows, m._columns);
				unsigned int i, j, k, p, q, r;
				for (i=0, p=0, r=0; i<result._rows; i++, p+=_columns, r+=result._columns)
					for (j=0; j<result._columns; j++)
					{
						ScalarType temp = 0;
						for (k=0, q=0; k<_columns; k++, q+=m._columns)
							temp+=(_data[p+k]*m._data[q+j]);
						result._data[r+j] = temp;
					}
				
				return result;
			};

			/*!
			*	Scalar sum.  
			*	\param k	
			*	\return		the resultant matrix 
			*/
			Matrix<TYPE> operator+(const TYPE k)
			{
				Matrix<TYPE> result(_rows, _columns);
				for (unsigned int i=0; i<_rows*_columns; i++)
					result._data[i] =  _data[i]+k;
				return result;
			};

			/*!
			*	Scalar difference.  
			*	\param k	
			*	\return		the resultant matrix 
			*/
			Matrix<TYPE> operator-(const TYPE k)
			{
				Matrix<TYPE> result(_rows, _columns);
				for (unsigned int i=0; i<_rows*_columns; i++)
					results._data[i] =  _data[i]-k;
				return result;
			};

			/*!
			*	Negate all matrix elements
			*	\return ...
			*/
			Matrix<TYPE> operator-() const
			{
				Matrix<TYPE> result(_rows, _columns, _data);
				for (unsigned int i=0; i<_columns*_rows; i++)
					result._data[i] = -1*_data[i];
				return result;
			};

			/*!
			*	Scalar multiplication.  
			*	\param k	value to multiply every member by 
			*	\return		the resultant matrix 
			*/
			Matrix<TYPE> operator*(const TYPE k)
			{
				Matrix<TYPE> result(_rows, _columns);
				for (unsigned int i=0; i<_rows*_columns; i++)
					results._data[i] =  _data[i]*k;
				return result;
			};

			/*!
			*	Scalar division.  
			*	\param k	value to divide every member by 
			*	\return		the resultant matrix 
			*/
			Matrix<TYPE> operator/(const TYPE k)
			{
				Matrix<TYPE> result(_rows, _columns);
				for (unsigned int i=0; i<_rows*_columns; i++)
					results._data[i] =  _data[i]/k;
				return result;
			};


			/*!
			*	Set all the matrix elements to zero.
			*/
			void SetZero()
			{
				for (unsigned int i=0; i<_rows*_columns; i++)
					_data[i] = ScalarType(0.0);
			};

			/*!
			*	Set the values of <I>j</I>-th column to v[j]
			*	\param j	the column index
			*	\param v	...
			*/
			void SetColumn(const unsigned int j, TYPE* v)
			{
				assert(j>=0 && j<_columns);
				unsigned int i, p;
				for (i=0, p=0; i<_rows; i++, p+=_columns)
					_data[p] = v[i];
			};

			/*!
			*	Set the elements of the <I>i</I>-th row to v[j] 
			*	\param i	the row index
			*	\param v	...
			*/
			void SetRow(const unsigned int i, TYPE* v)
			{
				assert(i>=0 && i<_rows);
				unsigned int j, p;
				for (j=0, p=i*_rows; j<_columns; j++, p++)
					_data[p] = v[j];
			};

			/*!
			*	Resize the current matrix.
			*	\param m the number of matrix rows.
			* \param n the number of matrix columns.
			*/
			void Resize(const unsigned int m, const unsigned int n)
			{
				assert(m>=2);
				assert(n>=2);
				_rows = m;
				_columns = n;
				delete []_data;
				_data = new ScalarType[m*n];
				for (unsigned int i=0; i<m*n; i++)
					_data[i] = 0;
			};


			/*!
			*	Matrix transposition operation: set the current matrix to its transpose
			*/
			void Transpose()
			{
				ScalarType *temp = new ScalarType[_rows*_columns];
				unsigned int i, j, p, q;
				for (i=0, p=0; i<_rows; i++, p+=_columns)
					for (j=0, q=0; j<_columns; j++, q+=_rows)
						temp[q+i] = _data[p+j];
				
				std::swap(_columns, _rows);
				std::swap(_data, temp);
				delete []temp;
			};
		
			/*!
			*	Print all matrix elements
			*/
			void Dump()
			{
				unsigned int i, j, p;
				for (i=0, p=0; i<_rows; i++, p+=_columns)
				{
					printf("[\t");
					for (j=0; j<_columns; j++)
						printf("%g\t", _data[p+j]);
					
					printf("]\n");
				}
				printf("\n");
			};

		protected:
			///	the number of matrix rows
			unsigned int _rows;

			///	the number of matrix rows
			unsigned int _columns;

			/// the matrix elements 
			ScalarType *_data;
		};

  /*! @} */
}; // end of namespace