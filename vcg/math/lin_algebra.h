#include <vcg/math/matrix44.h>

namespace vcg
{
	/** \addtogroup math */
	/* @{ */

	/*!
	*	Computes all eigenvalues and eigenvectors of a real symmetric matrix .
	*	On output, elements of the input matrix above the diagonal are destroyed. 
	* \param d  returns the eigenvalues of a. 
	* \param v  is a matrix whose columns contain, the normalized eigenvectors 
	* \param nrot returns the number of Jacobi rotations that were required. 
	*/
	template <typename TYPE>
	static void Jacobi(Matrix44<TYPE> &w, Point4<TYPE> &d, Matrix44<TYPE> &v, int &nrot) 
	{ 
		int j,iq,ip,i; 
		//assert(w.IsSymmetric());
		TYPE tresh, theta, tau, t, sm, s, h, g, c; 
		Point4<TYPE> b, z; 

		v.SetIdentity();

		for (ip=0;ip<4;++ip)			//Initialize b and d to the diagonal of a. 
		{		
			b[ip]=d[ip]=w[ip][ip]; 
			z[ip]=0.0;							//This vector will accumulate terms of the form tapq as in equation (11.1.14). 
		}
		nrot=0; 
		for (i=0;i<50;i++) 
		{ 
			sm=0.0; 
			for (ip=0;ip<3;++ip)		// Sum off diagonal elements
			{
				for (iq=ip+1;iq<4;++iq) 
					sm += fabs(w[ip][iq]); 
			} 
			if (sm == 0.0)					//The normal return, which relies on quadratic convergence to machine underflow. 
			{				
				return; 
			} 
			if (i < 4) 	
				tresh=0.2*sm/(4*4); //...on the first three sweeps. 
			else 		
				tresh=0.0;				//...thereafter. 
			for (ip=0;ip<4-1;++ip) 
			{  
				for (iq=ip+1;iq<4;iq++) 
				{ 
					g=100.0*fabs(w[ip][iq]); 
					//After four sweeps, skip the rotation if the off-diagonal element is small. 
					if(i>4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip]) && (float)(fabs(d[iq])+g) == (float)fabs(d[iq])) 
						w[ip][iq]=0.0; 
					else if (fabs(w[ip][iq]) > tresh) 
					{ 
						h=d[iq]-d[ip]; 
						if ((float)(fabs(h)+g) == (float)fabs(h)) 
							t=(w[ip][iq])/h; //t =1/(2#) 
						else 
						{ 
							theta=0.5*h/(w[ip][iq]); //Equation (11.1.10). 
							t=1.0/(fabs(theta)+sqrt(1.0+theta*theta)); 
							if (theta < 0.0) t = -t; 
						} 
						c=1.0/sqrt(1+t*t); 
						s=t*c; 
						tau=s/(1.0+c); 
						h=t*w[ip][iq]; 
						z[ip] -= h; 
						z[iq] += h; 
						d[ip] -= h; 
						d[iq] += h; 
						w[ip][iq]=0.0; 
						for (j=0;j<=ip-1;j++) { //Case of rotations 1 <= j < p. 
							JacobiRotate<TYPE>(w,s,tau,j,ip,j,iq) ;
						} 
						for (j=ip+1;j<=iq-1;j++) { //Case of rotations p < j < q. 
							JacobiRotate<TYPE>(w,s,tau,ip,j,j,iq);
						} 
						for (j=iq+1;j<4;j++) { //Case of rotations q< j <= n. 
							JacobiRotate<TYPE>(w,s,tau,ip,j,iq,j);
						} 
						for (j=0;j<4;j++) { 
							JacobiRotate<TYPE>(v,s,tau,j,ip,j,iq);
						} 
						++nrot; 
					} 
				} 
			} 
			for (ip=0;ip<4;ip++) 
			{ 
				b[ip] += z[ip]; 
				d[ip]=b[ip]; //Update d with the sum of ta_pq , 
				z[ip]=0.0; //and reinitialize z. 
			} 
		} 
	};

	/*!
	*
	*/
	template< typename TYPE >
		void JacobiRotate(Matrix44<TYPE> &A, TYPE s, TYPE tau, int i,int j,int k,int l)
	{
		TYPE g=A[i][j];
		TYPE h=A[k][l];
		A[i][j]=g-s*(h+g*tau);
		A[k][l]=h+s*(g-h*tau); 
	};

	/*!
	*	Given a matrix <I>A<SUB>m×n</SUB></I>, this routine computes its singular value decomposition,
	*	i.e. <I>A=U·W·V<SUP>T</SUP></I>. The matrix <I>A</I> will be destroyed!
	*	\param A	...
	*	\param W	the diagonal matrix of singular values <I>W</I>, stored as a vector <I>W[1...N]</I>
	*	\param V	the matrix <I>V</I> (not the transpose <I>V<SUP>T</SUP></I>)
	*/
	template <typename MATRIX_TYPE>
		static void SingularValueDecomposition(MATRIX_TYPE &A, typename MATRIX_TYPE::ScalarType *W, MATRIX_TYPE &V)
	{
		typedef typename MATRIX_TYPE::ScalarType ScalarType;

	};  


	/*!
	*	Solves A·X = B for a vector X, where A is specified by the matrices <I>U<SUB>m×n</SUB></I>, 
	*	<I>W<SUB>n×1</SUB></I> and <I>V<SUB>n×n</SUB></I> as returned by <CODE>SingularValueDecomposition</CODE>.
	*	No input quantities are destroyed, so the routine may be called sequentially with different b’s.
	*	\param x	is the output solution vector (<I>x<SUB>n×1</SUB></I>)
	*	\param b	is the input right-hand side (<I>b<SUB>n×1</SUB></I>)
	*/
	template <typename MATRIX_TYPE>
		static void SingularValueBacksubstitution(const MATRIX_TYPE												&U,
		const typename MATRIX_TYPE::ScalarType	*W,
		const MATRIX_TYPE												&V,
		typename MATRIX_TYPE::ScalarType	*x,
		const typename MATRIX_TYPE::ScalarType *b)
	{
		unsigned int jj, j, i;
		ScalarType s;
		ScalarType tmp	=	new ScalarType[U._columns];
		for (j=0; j<U._columns; j++) //Calculate U^T * B.
		{			
			s = 0;
			if (W[j]!=0)							//Nonzero result only if wj is nonzero.
			{ 
				for (i=0; i<U._rows; i++) 
					s += U[i][j]*b[i];
				s /= w[j];							//This is the divide by wj .
			}
			tmp[j]=s;
		}
		for (j=0;j<U._columns;j++)	//Matrix multiply by V to get answer.
		{			
			s = 0;
			for (jj=0; jj<U._columns; jj++) 
				s += V[j][jj]*tmp[jj];
			x[j]=s;
		}
		delete []tmp;
	};

	// Computes (a^2 + b^2)^(1/2) without destructive underflow or overflow.
	template <typename TYPE>
		inline TYPE pythagora(TYPE a, TYPE b)
	{
		TYPE abs_a = fabs(a);
		TYPE abs_b = fabs(b);
		if (abs_a > abs_b) 
			return abs_a*sqrt(1.0+sqr(abs_b/abs_a));
		else 
			return (abs_b == 0.0 ? 0.0 : abs_b*sqrt(1.0+sqr(abs_a/abs_b)));
	};

	/*! @} */
}; // end of namespace