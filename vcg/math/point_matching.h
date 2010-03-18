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

$Log: point_matching.h,v $

****************************************************************************/
#ifndef _VCG_MATH_POINTMATCHING_H
#define _VCG_MATH_POINTMATCHING_H

#include <vcg/math/matrix33.h>
#include <vcg/math/quaternion.h>
#include <vcg/math/lin_algebra.h>
namespace vcg
{
template<class ScalarType> 
class PointMatching 
{
public:
  typedef Point3<ScalarType> Point3x;
  typedef Matrix33<ScalarType> Matrix33x;
  typedef Matrix44<ScalarType> Matrix44x;
  typedef Quaternion<ScalarType> Quaternionx;


/*
Compute a similarity matching (rigid + uniform scaling)
simply create a temporary point set with the correct scaling factor


*/ 
static bool ComputeSimilarityMatchMatrix(		Matrix44x &res,
                                    std::vector<Point3x> &Pfix,		// vertici corrispondenti su fix (rossi)
						std::vector<Point3x> &Pmov) 		// normali scelti su mov (verdi)
{
	Quaternionx qtmp;
	Point3x tr;
	
	std::vector<Point3x> Pnew(Pmov.size());
	
	ScalarType scalingFactor=0;
	
	for(size_t i=0;i<( Pmov.size()-1);++i)
	{
			scalingFactor += Distance(Pmov[i],Pmov[i+1])/ Distance(Pfix[i],Pfix[i+1]);
#ifdef _DEBUG
			printf("Scaling Factor is %f",scalingFactor/(i+1));
#endif
	}
	scalingFactor/=(Pmov.size()-1);

	for(size_t i=0;i<Pmov.size();++i)
		Pnew[i]=Pmov[i]/scalingFactor;		
		
	
	bool ret=ComputeRigidMatchMatrix(res,Pfix,Pnew,qtmp,tr);
	if(!ret) return false;
	Matrix44x scaleM; scaleM.SetDiagonal(1.0/scalingFactor);
	
	res = res * scaleM;
	return true;
}



static bool ComputeRigidMatchMatrix(		Matrix44x &res,
                                    std::vector<Point3x> &Pfix,		// vertici corrispondenti su fix (rossi)
						std::vector<Point3x> &Pmov) 		// normali scelti su mov (verdi)
{
	Quaternionx qtmp;
	Point3x tr;
	return ComputeRigidMatchMatrix(res,Pfix,Pmov,qtmp,tr);
}


/* 
Calcola la matrice di rototraslazione 
che porta i punti Pmov su Pfix

Basata sul paper 

Besl, McKay
A method for registration o f 3d Shapes 
IEEE TPAMI Vol 14, No 2 1992


	Esempio d'uso 
			const int np=1000;
			std::vector<Point3x> pfix(np),pmov(np);

			Matrix44x Rot,Trn,RotRes;
			Rot.Rotate(30,Point3x(1,0,1));
			Trn.Translate(0,0,100);
			Rot=Trn*Rot;
			
			for(int i=0;i<np;++i){
				pfix[i]=Point3x(-150+rand()%1000,-150+rand()%1000,0);
				pmov[i]=Rot.Apply(pfix[i]);
			}
			
			ComputeRigidMatchMatrix(RotRes,pfix,pmov);
      
			RotRes.Invert();
			assert( RotRes==Rot);  
			assert( RotRes.Apply(pmov[i]) == pfix[i] );
			
*/
static
bool ComputeWeightedRigidMatchMatrix(Matrix44x &res,
                  std::vector<Point3x> &Pfix,		
									std::vector<Point3x> &Pmov,
									std::vector<ScalarType> weights,
									Quaternionx &q,
									Point3x &tr
									) 	
{

  Matrix33x ccm; 
	Point3x bfix,bmov; // baricenter of src e trg
	ccm.WeightedCrossCovariance(weights,Pmov,Pfix,bmov,bfix);
	Matrix33x cyc; // the cyclic components of the cross covariance matrix.

	cyc=ccm - ccm.transpose();

	Matrix44x QQ;
	QQ.SetZero();
	Point3x D(cyc[1][2],cyc[2][0],cyc[0][1]);

  Matrix33x RM;
	RM.SetZero();
	RM[0][0]=-ccm.Trace();
  RM[1][1]=-ccm.Trace();
  RM[2][2]=-ccm.Trace();
  RM += ccm + ccm.transpose();

	QQ[0][0] = ccm.Trace();
	QQ[0][1] = D[0]; QQ[0][2] = D[1]; QQ[0][3] = D[2];
	QQ[1][0] = D[0]; QQ[2][0] = D[1];	QQ[3][0] = D[2];

	int i,j;
  for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			QQ[i+1][j+1]=RM[i][j];

//  printf(" Quaternion Matrix\n");
//	print(QQ);
	Point4d d;
  Matrix44x v;
	int nrot;
	Jacobi(QQ,d,v,nrot);
//	printf("Done %i iterations\n %f %f %f %f\n",nrot,d[0],d[1],d[2],d[3]);
//	print(v);
	// Now search the maximum eigenvalue
	double maxv=0;
	int maxind=-1;
  for(i=0;i<4;i++)
		if(maxv<fabs(d[i])) {
			q=Quaternionx(v[0][i],v[1][i],v[2][i],v[3][i]);
			maxind=i;
			maxv=d[i];
		}
  // The corresponding eigenvector define the searched rotation,
		Matrix44x Rot;
	q.ToMatrix(Rot);
  // the translation (last row) is simply the difference between the transformed src barycenter and the trg baricenter
	tr= (bfix - Rot *bmov);
	//res[3][0]=tr[0];res[3][1]=tr[1];res[3][2]=tr[2];
	Matrix44x Trn;
	Trn.SetTranslate(tr);
		
	res=Rot*Trn;
	return true;
}

static
bool ComputeRigidMatchMatrix(Matrix44x &res,
 						std::vector<Point3x> &Pfix,		
						std::vector<Point3x> &Pmov,
							Quaternionx &q,
							Point3x &tr) 	
{

  Matrix33x ccm; 
	Point3x bfix,bmov; // baricenter of src e trg
	ccm.CrossCovariance(Pmov,Pfix,bmov,bfix);
	Matrix33x cyc; // the cyclic components of the cross covariance matrix.

	cyc=ccm-ccm.transpose();

	Matrix44x QQ;
	QQ.SetZero();
	Point3x D(cyc[1][2],cyc[2][0],cyc[0][1]);

  Matrix33x RM;
	RM.SetZero();
	RM[0][0]=-ccm.Trace();
  RM[1][1]=-ccm.Trace();
  RM[2][2]=-ccm.Trace();
  RM += ccm + ccm.transpose();

	QQ[0][0] = ccm.Trace();
	QQ[0][1] = D[0]; QQ[0][2] = D[1]; QQ[0][3] = D[2];
	QQ[1][0] = D[0]; QQ[2][0] = D[1];	QQ[3][0] = D[2];

	int i,j;
  for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			QQ[i+1][j+1]=RM[i][j];

//  printf(" Quaternion Matrix\n");
//	print(QQ);
	Point4d d;
  Matrix44x v;
	int nrot;
	//QQ.Jacobi(d,v,nrot);
	Jacobi(QQ,d,v,nrot);
//	printf("Done %i iterations\n %f %f %f %f\n",nrot,d[0],d[1],d[2],d[3]);
//	print(v);
	// Now search the maximum eigenvalue
	double maxv=0;
	int maxind=-1;
  for(i=0;i<4;i++)
		if(maxv<fabs(d[i])) {
			q=Quaternionx(v[0][i],v[1][i],v[2][i],v[3][i]);
			maxind=i;
			maxv=d[i];
		}
  // The corresponding eigenvector define the searched rotation,
	Matrix44x Rot;
	q.ToMatrix(Rot);
  // the translation (last row) is simply the difference between the transformed src barycenter and the trg baricenter
	tr= (bfix - Rot*bmov);
	//res[3][0]=tr[0];res[3][1]=tr[1];res[3][2]=tr[2];
	Matrix44x Trn;
	Trn.SetTranslate(tr);
		
	res=Trn*Rot;
	return true;
}

// Dati due insiemi di punti e normali corrispondenti calcola la migliore trasformazione 
// che li fa corrispondere
static bool ComputeMatchMatrix(		Matrix44x &res,
 						std::vector<Point3x> &Ps,		// vertici corrispondenti su src (rossi)
						std::vector<Point3x> &Ns, 		// normali corrispondenti su src (rossi)
						std::vector<Point3x> &Pt)		// vertici scelti su trg (verdi) 
//						vector<Point3x> &Nt) 		// normali scelti su trg (verdi)
{
  assert(0);
  // Da qui in poi non compila che ha bisogno dei minimiquadrati
#if 0
  int sz=Ps.size();

	Matrix<double> A(sz,12);
	Vector<double> b(sz);
	Vector<double> x(12);

	//inizializzo il vettore per minimi quadrati
	// la matrice di trasf che calcolo con LeastSquares cerca avvicinare il piu' 
	// possibile le coppie di punti che trovo ho scelto  
	// Le coppie di punti sono gia' trasformate secondo la matrice <In> quindi come scelta iniziale 
	// per il metodo minimiquadrati scelgo l'identica (e.g. se ho allineato a mano perfettamente e 
	// le due mesh sono perfettamente uguali DEVE restituire l'identica)
	
	res.SetIdentity();
	int i,j,k;
	for(i=0; i<=2; ++i)
		for(j=0; j<=3; ++j)
			x[i*4+j] = res[i][j];


	//costruzione della matrice
	for(i=0;i<sz;++i)
	{
		for(j=0;j<3;++j)
			for(k=0;k<4;++k)
				if(k<3)
				{
					A[i][k+j*4] = Ns[i][j]*Pt[i][k];
				}
				else
				{
					A[i][k+j*4] = Ns[i][j];
				}
		b[i] = Ps[i]*Ns[i];
	}
	const int maxiter = 4096;
	int iter;
	LSquareGC(x,A,b,1e-16,maxiter,iter);
	
	TRACE("LSQ Solution");
	for(int ind=0; ind<12; ++ind) {
		if((ind%4)==0) TRACE("\n");
		TRACE("%8.5lf ", x[ind]); 
	} TRACE("\n");

	if(iter==maxiter)
	{
		TRACE("I minimi quadrati non convergono!!\n");
		return false;
	}
	else { TRACE("Convergenza in %d passi\n",iter); }

	//Devo riapplicare la matrice di trasformazione globale a 
	//trg inserendo il risultato nel vettore trgvert contenente 
	//copia dei suoi vertici
	Matrix44x tmp;
	for(i=0; i<=2; ++i)
		for(j=0; j<=3; ++j)
			res[j][i] = x[i*4+j];
	res[0][3] = 0.0;
	res[1][3] = 0.0;
	res[2][3] = 0.0;
	res[3][3] = 1.0;
	/*
	res.Transpose();
	Point3x scv,shv,rtv,trv;
	res.Decompose(scv,shv,rtv,trv);
	vcg::print(res);
	printf("Scale %f %f %f\n",scv[0],scv[1],scv[2]);
	printf("Shear %f %f %f\n",shv[0],shv[1],shv[2]);
	printf("Rotat %f %f %f\n",rtv[0],rtv[1],rtv[2]);
	printf("Trans %f %f %f\n",trv[0],trv[1],trv[2]);
	
	printf("----\n"); res.Decompose(scv,shv,rtv,trv);
	vcg::print(res);
	printf("Scale %f %f %f\n",scv[0],scv[1],scv[2]);
	printf("Shear %f %f %f\n",shv[0],shv[1],shv[2]);
	printf("Rotat %f %f %f\n",rtv[0],rtv[1],rtv[2]);
	printf("Trans %f %f %f\n",trv[0],trv[1],trv[2]);
	
	res.Transpose();
	*/
#endif
	return true;
}

/*
****** Questa parte per compilare ha bisogno di leastsquares e matrici generiche 
****** Da controllare meglio


static void CreatePairMatrix( Matrix<double> & A2, const Point3x & p, const Point3x & n, double d )
{	
	double t1 = p[0]*p[0];
	double t2 = n[0]*n[0];
	double t4 = t1*n[0];
	double t5 = t4*n[1];
	double t6 = t4*n[2];
	double t7 = p[0]*t2;
	double t8 = t7*p[1];
	double t9 = p[0]*n[0];
	double t10 = p[1]*n[1];
	double t11 = t9*t10;
	double t12 = p[1]*n[2];
	double t13 = t9*t12;
	double t14 = t7*p[2];
	double t15 = p[2]*n[1];
	double t16 = t9*t15;
	double t17 = p[2]*n[2];
	double t18 = t9*t17;
	double t19 = t9*n[1];
	double t20 = t9*n[2];
	double t21 = t9*d;
	double t22 = n[1]*n[1];
	double t25 = t1*n[1]*n[2];
	double t26 = p[0]*t22;
	double t27 = t26*p[1];
	double t28 = p[0]*n[1];
	double t29 = t28*t12;
	double t30 = t26*p[2];
	double t31 = t28*t17;
	double t32 = t28*n[2];
	double t33 = t28*d;
	double t34 = n[2]*n[2];

	double t36 = p[0]*t34;
	double t41 = p[1]*p[1]; double t43 = t41*n[0];
	double t46 = p[1]*t2;   double t48 = p[1]*n[0];
	double t49 = t48*t15;   double t50 = t48*t17;
	double t51 = t48*n[1];  double t52 = t48*n[2];
	double t57 = p[1]*t22;  double t59 = t10*t17;
	double t60 = t10*n[2];  double t63 = p[1]*t34;
	double t66 = p[2]*p[2]; double t68 = t66*n[0];
	double t72 = p[2]*n[0]; double t73 = t72*n[1];
	double t74 = t72*n[2];	double t80 = t15*n[2];
	
	A2[0][0] = t1*t2; A2[0][1] = t5;  A2[0][2] = t6;
	A2[0][3] = t8;    A2[0][4] = t11; A2[0][5] = t13;
	A2[0][6] = t14;   A2[0][7] = t16; A2[0][8] = t18;
	A2[0][9] = t7;   A2[0][10] = t19; A2[0][11] = t20;
	A2[0][12] = -t21;
	
	A2[1][1] = t1*t22; A2[1][2]  = t25; A2[1][3] = t11;
	A2[1][4] = t27;    A2[1][5]  = t29; A2[1][6] = t16;
	A2[1][7] = t30;    A2[1][8]  = t31; A2[1][9] = t19;
	A2[1][10] = t26;   A2[1][11] = t32; A2[1][12] = -t33;
	
	A2[2][2] = t1*t34; A2[2][3] = t13; A2[2][4] = t29;
	A2[2][5] = t36*p[1];    A2[2][6] = t18; A2[2][7] = t31;
	A2[2][8] = t36*p[2];    A2[2][9] = t20; A2[2][10] = t32;
	A2[2][11] = t36;   A2[2][12] = -p[0]*n[2]*d;
	
	A2[3][3] = t41*t2; A2[3][4] = t43*n[1]; A2[3][5] = t43*n[2];
	A2[3][6] = t46*p[2]; A2[3][7] = t49;  A2[3][8] = t50;
    A2[3][9] = t46; A2[3][10] = t51; A2[3][11] = t52;
    A2[3][12] = -t48*d;
	
	A2[4][4]  = t41*t22;  A2[4][5]  = t41*n[1]*n[2]; A2[4][6] = t49;
	A2[4][7]  = t57*p[2]; A2[4][8]  = t59; A2[4][9] = t51;
	A2[4][10] = t57;      A2[4][11] = t60; A2[4][12] = -t10*d;
	
	A2[5][5]  = t41*t34;  A2[5][6] = t50; A2[5][7] = t59;
	A2[5][8]  = t63*p[2]; A2[5][9] = t52; A2[5][10] = t60;
	A2[5][11] = t63;      A2[5][12] = -t12*d;
	
	A2[6][6]  = t66*t2;  A2[6][7] = t68*n[1]; A2[6][8] = t68*n[2];
	A2[6][9]  = p[2]*t2; A2[6][10] = t73;     A2[6][11] = t74;
	A2[6][12] = -t72*d;
	
	A2[7][7] = t66*t22;   A2[7][8] = t66*n[1]*n[2]; A2[7][9] = t73;
	A2[7][10] = p[2]*t22; A2[7][11] = t80;          A2[7][12] = -t15*d;
	
	A2[8][8] = t66*t34;   A2[8][9] = t74; A2[8][10] = t80;
	A2[8][11] = p[2]*t34; A2[8][12] = -t17*d;
	
	A2[9][9]   = t2;        A2[9][10]  = n[0]*n[1];
	A2[9][11]  = n[0]*n[2]; A2[9][12]  = -n[0]*d;
	
	A2[10][10] = t22; A2[10][11] = n[1]*n[2]; A2[10][12] = -n[1]*d;
	A2[11][11] = t34; A2[11][12] = -n[2]*d;
	A2[12][12] = d*d;
}

// Dati due insiemi di punti e normali corrispondenti calcola la migliore trasformazione 
// che li fa corrispondere
static bool ComputeMatchMatrix2(		Matrix44x &res,
 						std::vector<Point3x> &Ps,		// vertici corrispondenti su src (rossi)
						std::vector<Point3x> &Ns, 		// normali corrispondenti su src (rossi)
						std::vector<Point3x> &Pt)		// vertici scelti su trg (verdi) 
						//vector<Point3x> &Nt) 		// normali scelti su trg (verdi)
{
	const int N = 13;
	int i,j,k;

	Matrixd AT(N,N);
	Matrixd TT(N,N);
		// Azzeramento matrice totale (solo tri-superiore)
	for(i=0;i<N;++i)
		for(j=i;j<N;++j)
			AT[i][j] = 0;
		// Calcolo matrici locali e somma
	for(k=0;k<Ps.size();++k)
	{		
		CreatePairMatrix(TT,Pt[k],Ns[k],Ps[k]*Ns[k]);
		for(i=0;i<N;++i)
			for(j=i;j<N;++j)
				AT[i][j] += TT[i][j];
	}

	for(i=0;i<N;++i)
		for(j=0;j<i;++j)
				AT[i][j] = AT[j][i];

	std::vector<double> q;
	double error;
	affine_ls2(AT,q,error);
	//printf("error: %g \n",error);
	res[0][0] = q[0];
	res[0][1] = q[1];
	res[0][2] = q[2];
	res[0][3] = 0;
	res[1][0] = q[3];
	res[1][1] = q[4];
	res[1][2] = q[5];
	res[1][3] = 0;
	res[2][0] = q[6];
	res[2][1] = q[7];
	res[2][2] = q[8];
	res[2][3] = 0;
	res[3][0] = q[9];
	res[3][1] = q[10];
	res[3][2] = q[11];
	res[3][3] = q[12];

	return true;
}
*/
};
} // end namespace

#endif
