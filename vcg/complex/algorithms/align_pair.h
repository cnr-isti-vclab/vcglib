/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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
#ifndef VCG_ALIGN_PAIR_H
#define VCG_ALIGN_PAIR_H

#include <ctime>
#include <stdio.h>
#include <vcg/math/histogram.h>
#include <vcg/math/matrix44.h>
#include <vcg/math/random_generator.h>
#include <vcg/math/gen_normal.h>

#include <vcg/space/point_matching.h>
#include <vcg/space/index/grid_static_ptr.h>

#include <vcg/simplex/face/component_ep.h>

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/component_ep.h>
#include <vcg/complex/algorithms/update/position.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/point_matching_scale.h>


namespace vcg {
/**
Class for aligning pair of meshes using ICP (iterated closest point)
*/

class AlignPair {
public:

	AlignPair()
	{
		clear();
		myrnd.initialize(time(NULL));
	}

	enum ErrorCode {
		SUCCESS,
		NO_COMMON_BBOX,
		TOO_FEW_POINTS,
		LSQ_DIVERGE,
		TOO_MUCH_SHEAR,
		TOO_MUCH_SCALE,
		FORBIDDEN,
		INVALID,
		UNKNOWN_MODE };


	/*********************** Utility classes ****************************/

	class A2Vertex;

	class A2Face;

	class A2UsedTypes:
			public vcg::UsedTypes < vcg::Use<A2Vertex>::AsVertexType,
			vcg::Use<A2Face  >::AsFaceType >{};

	class A2Vertex : public vcg::Vertex<A2UsedTypes,vcg::vertex::Coord3d,vcg::vertex::Normal3d,vcg::vertex::BitFlags> {};
	class A2Face   : public vcg::Face< A2UsedTypes,vcg::face::VertexRef, vcg::face::Normal3d,vcg::face::Mark,vcg::face::BitFlags> {};

	class A2Mesh   : public vcg::tri::TriMesh< std::vector<A2Vertex>, std::vector<A2Face> >
	{
	public:
		//bool Import(const char *filename) { Matrix44d Tr; Tr.SetIdentity(); return Import(filename,Tr);}
		//bool Import(const char *filename, Matrix44d &Tr);

		inline bool initVert(const Matrix44d &Tr) {
			Matrix44d Idn; Idn.SetIdentity();
			if (Tr != Idn)
				tri::UpdatePosition<A2Mesh>::Matrix(*this, Tr);
			tri::UpdateNormal<A2Mesh>::NormalizePerVertex(*this);
			tri::UpdateBounding<A2Mesh>::Box(*this);
			return true;
		}
		inline bool init(const Matrix44d &Tr) {
			Matrix44d Idn; Idn.SetIdentity();
			tri::Clean<A2Mesh>::RemoveUnreferencedVertex(*this);
			if (Tr != Idn)
				tri::UpdatePosition<A2Mesh>::Matrix(*this, Tr);
			tri::UpdateNormal<A2Mesh>::PerVertexNormalizedPerFaceNormalized(*this);
			tri::UpdateFlags<A2Mesh>::FaceBorderFromNone(*this);
			tri::UpdateBounding<A2Mesh>::Box(*this);

			return true;
		}
	};

	typedef A2Mesh::FaceContainer FaceContainer;
	typedef A2Mesh::FaceType      FaceType;
	typedef GridStaticPtr<FaceType, double > A2Grid;
	typedef GridStaticPtr<A2Mesh::VertexType, double > A2GridVert;

	class Stat
	{
	public:

		class IterInfo
		{
		public:
			IterInfo()
			{
				memset ( (void *) this, 0, sizeof(IterInfo));
			}

			double MinDistAbs;
			int DistanceDiscarded;
			int AngleDiscarded;
			int BorderDiscarded;
			int SampleTested;  // how many points have been tested 
			int SampleUsed;    // how many points have been actually used to compute the transformation
			double pcl50;
			double pclhi;
			double AVG;
			double RMS;
			double StdDev;
			int Time;  // Ending time of this iteration

		};

		std::vector<IterInfo> I;

		double lastPcl50() const
		{
			return I.back().pcl50;
		}

		int lastSampleUsed() const {
			return I.back().SampleUsed;
		}

		int MovVertNum;
		int FixVertNum;
		int FixFaceNum;

		int totTime() {
			return I.back().Time-StartTime;
		}

		int iterTime(unsigned int i) const
		{
			const int clock_per_ms = std::max<int>(CLOCKS_PER_SEC / 1000,1);
			assert(i<I.size());
			if(i==0) return  (I[i].Time-StartTime )/clock_per_ms;
			else return (I[i].Time - I[i-1].Time)/clock_per_ms ;
		}
		int StartTime;

		inline void clear()
		{
			I.clear();
			StartTime = 0;
			MovVertNum = 0;
			FixVertNum = 0;
			FixFaceNum = 0;
		}

		inline void dump(FILE *fp)
		{
			if (I.size() == 0) {
				fprintf(fp, "Empty AlignPair::Stat\n");
				return;
			}
			fprintf(fp, "Final Err %8.5f In %i iterations Total Time %ims\n", lastPcl50(), (int)I.size(), totTime());
			fprintf(fp, "Mindist   Med   Hi    Avg  RMS   StdDev   Time Tested Used  Dist Bord Angl \n");
			for (unsigned int qi = 0; qi < I.size(); ++qi)
				fprintf(
							fp,
							"%5.2f (%6.3f:%6.3f) (%6.3f %6.3f %6.3f) %4ims %5i %5i %4i+%4i+%4i\n",
							I[qi].MinDistAbs,
							I[qi].pcl50, I[qi].pclhi,
							I[qi].AVG, I[qi].RMS, I[qi].StdDev,
							iterTime(qi),
							I[qi].SampleTested, I[qi].SampleUsed, I[qi].DistanceDiscarded, I[qi].BorderDiscarded, I[qi].AngleDiscarded);
		}

		//! Write a HTML table with the values
		inline void htmlDump(FILE *fp)
		{
			fprintf(fp, "Final Err %8.5f In %i iterations Total Time %ims\n", lastPcl50(), (int)I.size(), totTime());
			fprintf(fp, "<table border>\n");
			fprintf(fp, "<tr> <th>Mindist</th><th>    50ile </th><th>  Hi </th><th>   Avg  </th><th> RMS </th><th>  StdDev  </th><th> Time </th><th> Tested </th><th> Used </th><th> Dist </th><th> Bord </th><th> Angl \n");
			for (unsigned int qi = 0; qi < I.size(); ++qi)
				fprintf(
						fp, "<tr> <td> %8.5f </td><td align=\"right\"> %9.6f </td><td align=\"right\"> %8.5f </td><td align=\"right\"> %5.3f </td><td align=\"right\"> %8.5f </td><td align=\"right\"> %9.6f </td><td align=\"right\"> %4ims </td><td align=\"right\"> %5i </td><td align=\"right\"> %5i </td><td align=\"right\"> %4i </td><td align=\"right\"> %4i </td><td align=\"right\">%4i </td><td align=\"right\"></tr>\n",
							I[qi].MinDistAbs,
							I[qi].pcl50, I[qi].pclhi,
							I[qi].AVG, I[qi].RMS, I[qi].StdDev,
							iterTime(qi),
							I[qi].SampleTested, I[qi].SampleUsed, I[qi].DistanceDiscarded, I[qi].BorderDiscarded, I[qi].AngleDiscarded);
			fprintf(fp, "</table>\n");
		}

		// Restituisce true se nelle ultime <lastiter> iterazioni non e' diminuito
		// l'errore
		inline bool stable(int lastiter)
		{
			if (I.empty())
				return false;
			int parag = int(I.size()) - lastiter;

			if (parag < 0)
				parag = 0;
			if (I.back().pcl50 < I[parag].pcl50)
				return false; // se siamo diminuiti non e' stabile

			return true;
		}

	};


	class Param
	{
	public:
		enum MatchModeEnum  {MMSimilarity, MMRigid};
		enum SampleModeEnum {SMRandom, SMNormalEqualized};

		Param()
		{
			SampleNum         = 2000;
			MaxPointNum       = 100000;
			MinPointNum       = 30;
			MaxIterNum        = 75;
			TrgDistAbs        = 0.005f;
			MinDistAbs        = 10;
			MaxAngleRad       = math::ToRad(45.0);
			MaxShear          = 0.5;
			MaxScale          = 0.5; // means that the scale must be between 1-MaxScale and 1+MaxScale
			PassHiFilter      = 0.75;
			ReduceFactorPerc  = 0.80;
			MinMinDistPerc    = 0.01;
			EndStepNum        = 5;
			MatchMode         = MMRigid;
			SampleMode        = SMNormalEqualized;
			UGExpansionFactor = 10;
		}

		int SampleNum;        //! The intial number of samples that are chosen on the fix mesh.
		int MaxPointNum;      //! Maximum number off point used to compute the transforamtion matrix (rarely used)
		int MinPointNum;      //! Minimal number of point that we have to find to consider the resulting alignment valid
		double MinDistAbs;    //! Minimal distance for ICP points. Only points that are closer than this threshold are chosen. 
		                      //! This threshold is iteratively lowered by the ReduceFactor parameter 
		
		double MaxAngleRad;	  //! Max angle (in radiant) for ICP points. Only points whose normal differ less than this threshold are considered
		
		int MaxIterNum;       //! Maximum number of iteration to be done during aligning
		double TrgDistAbs;    //! Target distance (half of the samples should be below this distance to consider the icp reach convergence)
		int EndStepNum;       //! max number of ICP iterations
		  
		double PassHiFilter;     //! Percentile filtering threshold. Usually we discard the farthest quartile of the matched points
		double ReduceFactorPerc; //! At each step we discard the points farther than a given threshold. The threshold is iterativeley reduced;
		                         //! StartMinDist= min(StartMinDist, 5.0*H.Percentile(ap.ReduceFactorPerc))


		double MinMinDistPerc;	//! Ratio between initial starting distance (MinDistAbs) and what can reach by the application of the ReduceFactor.

		int UGExpansionFactor;  //! Size of the uniform grid as a ration of the fix mesh size

		bool UseVertexOnly;     //! if true all the Alignment pipeline ignores faces and works over point clouds.

		double MaxShear;
		double MaxScale;
		MatchModeEnum MatchMode;
		SampleModeEnum SampleMode;
		//void Dump(FILE *fp,double BoxDiag);

	};

  // Class to store the result of an alignment between two meshes
  // the points are intended in the initial reference system of the two meshes.
  //
  // if the meshes have a basic transformation in input,
  // this appears only during the A2Mesh::Import and then is forever forgotten.
  // These points are therefore in the reference systems built during the Import
  // the matrix Tr that which
	// Tr*Pmov[i]== Pfix


	class Result
	{
	public:
		int MovName;
		int FixName;

		Matrix44d Tr;
		std::vector<Point3d> Pfix;  // Corresponding Points on the Fix Mesh (red)
		std::vector<Point3d> Nfix;  // Corresponding Normals  on the Fix Mesh (red)
		std::vector<Point3d> Pmov;  // Chosen Points on the Mov Mesh (green) before the transformation
		std::vector<Point3d> Nmov;  // Chosen Normals on the Mov Mesh (green)
		Histogramf H;
		Stat as;
		Param ap;
		ErrorCode status;
		bool isValid()
		{
			return status==SUCCESS;
		}
		double err;
		float area; // the overlapping area, a percentage as computed in Occupancy Grid.

		bool operator <  (const Result & rr) const {return (err< rr.err);}
		bool operator <= (const Result & rr) const {return (err<=rr.err);}
		bool operator >  (const Result & rr) const {return (err> rr.err);}
		bool operator >= (const Result & rr) const {return (err>=rr.err);}
		bool operator == (const Result & rr) const {return (err==rr.err);}
		bool operator != (const Result & rr) const {return (err!=rr.err);}

		std::pair<double,double> computeAvgErr() const
		{
			double sum_before=0;
			double sum_after=0;
			for(unsigned int ii=0;ii<Pfix.size();++ii) {
				sum_before+=Distance(Pfix[ii],   Pmov[ii]);
				sum_after+=Distance(Pfix[ii], Tr*Pmov[ii]);
			}
			return std::make_pair(sum_before/double(Pfix.size()),sum_after/double(Pfix.size()) ) ;
		}

	};

	/******************* End utility classes ************************/

	static inline const char* errorMsg(ErrorCode code)
	{
		switch (code){
		case SUCCESS:
			return "Success";
		case NO_COMMON_BBOX:
			return "No Common BBox";
		case TOO_FEW_POINTS:
			return "Too few points";
		case LSQ_DIVERGE:
			return "LSQ not converge";
		case TOO_MUCH_SHEAR:
			return "Too much shear";
		case TOO_MUCH_SCALE:
			return "Too much scale";
		case UNKNOWN_MODE:
			return "Unknown mode ";
		default:
			assert(0);
			return "Catastrophic Error";
		}
		return 0;
	}

	void clear()
	{
		status=SUCCESS;
	}

	/******* Data Members *********/

	std::vector<A2Vertex> *mov;
	A2Mesh *fix;

	ErrorCode status;
	AlignPair::Param ap;

	math::SubtractiveRingRNG myrnd;

	/**** End Data Members *********/

	template < class MESH >
	void convertMesh(MESH &M1, A2Mesh &M2)
	{
		tri::Append<A2Mesh,MESH>::MeshCopy(M2,M1);
	}

	template < class VERTEX >
	void convertVertex(const std::vector<VERTEX> &vert1, std::vector<A2Vertex> &vert2, Box3d *Clip=0)
	{
		vert2.clear();
		typename std::vector<VERTEX>::const_iterator vi;
		A2Vertex tv;
		Box3<typename VERTEX::ScalarType> bb;
		if(Clip){
			bb.Import(*Clip);
			for(vi=vert1.begin();vi<vert1.end();++vi)
				if(!(*vi).IsD() && bb.IsIn((*vi).cP())){
					tv.P().Import((*vi).cP());
					tv.N().Import((*vi).cN());
					vert2.push_back(tv);
				}
		}
		else {
			for(vi=vert1.begin();vi<vert1.end();++vi) {
				if(!(*vi).IsD()){
					tv.P().Import((*vi).cP());
					tv.N().Import((*vi).cN());
					vert2.push_back(tv);
				}
			}
		}
	}

	inline bool sampleMovVert(
			std::vector<A2Vertex> &vert,
			int sampleNum,
			AlignPair::Param::SampleModeEnum sampleMode)
	{
		switch (sampleMode)
		{
		case AlignPair::Param::SMRandom:
			return SampleMovVertRandom(vert, sampleNum);
		case AlignPair::Param::SMNormalEqualized:
			return SampleMovVertNormalEqualized(vert, sampleNum);
		default:
			assert(0);
			return false;
		}
	}

	inline bool SampleMovVertRandom(std::vector<A2Vertex> &vert, int sampleNum)
	{
		if (int(vert.size()) <= sampleNum)
			return true;
		for (int i = 0; i < sampleNum; ++i) {
			int pos = myrnd.generate(vert.size());
			assert(pos >= 0 && pos < int(vert.size()));
			std::swap(vert[i], vert[pos]);
		}
		vert.resize(sampleNum);
		return true;
	}

	bool SampleMovVertNormalEqualized(std::vector<A2Vertex> &vert, int sampleNum)
	{
		std::vector<Point3d> NV;
		if (NV.size() == 0) {
			GenNormal<double>::Fibonacci(30, NV);
			printf("Generated %i normals\n", int(NV.size()));
		}
		// Bucket vector dove, per ogni normale metto gli indici
		// dei vertici ad essa corrispondenti
		std::vector<std::vector <int> > BKT(NV.size());
		for (size_t i = 0; i < vert.size(); ++i) {
			int ind = GenNormal<double>::BestMatchingNormal(vert[i].N(), NV);
			BKT[ind].push_back(int(i));
		}

		// vettore di contatori per sapere quanti punti ho gia' preso per ogni bucket
		std::vector <int> BKTpos(BKT.size(), 0);

		if (sampleNum >= int(vert.size()))
			sampleNum = vert.size() - 1;

		for (int i = 0; i < sampleNum;) {
			int ind = myrnd.generate(BKT.size()); // Scelgo un Bucket
			int &CURpos = BKTpos[ind];
			std::vector<int> &CUR = BKT[ind];

			if (CURpos<int(CUR.size())) {
				std::swap(CUR[CURpos], CUR[CURpos + myrnd.generate(BKT[ind].size() - CURpos)]);
				std::swap(vert[i], vert[CUR[CURpos]]);
				++BKTpos[ind];
				++i;
			}
		}
		vert.resize(sampleNum);

		return true;
	}

	/*
	This function is used to choose remove outliers after each ICP iteration.
	All the points with a distance over the given Percentile are discarded.
	It uses two parameters
	MaxPointNum an (unused) hard limit on the number of points that are chosen
	MinPointNum the minimum number of points that have to be chosen to be usable
	*/
	inline bool choosePoints(
		std::vector<Point3d> &ps, // vertici corrispondenti su fix (rossi)
		std::vector<Point3d> &ns, // normali corrispondenti su fix (rossi)
		std::vector<Point3d> &pt, // vertici scelti su mov (verdi)
		std::vector<Point3d> &opt,		// vertici scelti su mov (verdi)
		//vector<Point3d> &Nt, 		// normali scelti su mov (verdi)
		double passHi,
		Histogramf &h)
	{
		const int N = ap.MaxPointNum;
		double newmaxd = h.Percentile(float(passHi));
		int sz = int(ps.size());
		int fnd = 0;
		int lastgood = sz - 1;
		math::SubtractiveRingRNG myrnd;
		while (fnd < N && fnd < lastgood) {
			int index = fnd + myrnd.generate(lastgood - fnd);
			double dd = Distance(ps[index], pt[index]);
			if (dd <= newmaxd){
				std::swap(ps[index], ps[fnd]);
				std::swap(ns[index], ns[fnd]);
				std::swap(pt[index], pt[fnd]);
				std::swap(opt[index], opt[fnd]);
				++fnd;
			}
			else {
				std::swap(ps[index], ps[lastgood]);
				std::swap(ns[index], ns[lastgood]);
				std::swap(pt[index], pt[lastgood]);
				std::swap(opt[index], opt[lastgood]);
				lastgood--;
			}
		}
		ps.resize(fnd);
		ns.resize(fnd);
		pt.resize(fnd);
		opt.resize(fnd);

		if ((int)ps.size() < ap.MinPointNum){
			printf("Troppi pochi punti!\n");
			ps.clear();
			ns.clear();
			pt.clear();
			opt.clear();
			return false;
		}
		return true;
	}

/*
Minimal example of code for using the align.

AlignPair::A2Mesh Mov,Fix;                   // The two meshes to be alligned. Mov will moved. 
vector<AlignPair::A2Vertex> MovVert;         // Points chosen on Mov Mesh to compute the alignment
Matrix44d In;	In.SetIdentity();            // Initial transformation to be applied to Mov mesh to bring it over the Fix mesh.

AlignPair aa;                                // The main align class
AlignPair::Param ap;
UGrid< AlignPair::A2Mesh::face_container > UG;

Fix.LoadPly("FixMesh.ply");                             // Standard ply loading
Mov.LoadPly("MovMesh.ply");                        
Fix.Init( Ident, false);                                // Basic init (computation of normals)
Mov.Init( Ident, false);                                
                                                   
AlignPair::InitFix(&Fix, ap, UG);                       // Init UG for quick search
                                                   
aa.ConvertVertex(Mov.vert,MovVert);                     // Sample the Mov vertex in order to find good samples to be used. 
aa.SampleMovVert(MovVert, ap.SampleNum, ap.SampleMode);

aa.mov=&MovVert;                                        // basic init of the align class with the chosen vert and the fix mesh
aa.fix=&Fix;                                        
aa.ap = ap;                                         

aa.Align(In,UG,res);                                    // Main ICP algorithm
                                                        // the matrix containing the computed alignment is in res.Tr;

res.as.Dump(stdout);
*/

	bool align(const Matrix44d &in, A2Grid &UG, A2GridVert &UGV, Result &res)
	{
		res.ap=ap;

		bool ret=align(UG, UGV, in, res.Tr, res.Pfix, res.Nfix, res.Pmov, res.Nmov, res.H, res.as);

		res.err=res.as.lastPcl50();
		res.status=status;
		return ret;
	}

	double abs2Perc(double val, Box3d bb) const
	{
		return val/bb.Diag();
	}

	double perc2Abs(double val, Box3d bb) const
	{
		return val*bb.Diag();
	}

/************************************************************************************
Versione Vera della Align a basso livello.

Si assume che la mesh fix sia gia' stata messa nella ug u con le debite trasformazioni.
in

************************************************************************************/

	/**
	The Main ICP alignment Function:
	It assumes that:
	we have two meshes:
	- Fix the mesh that does not move and stays in the spatial indexing structure.
	- Mov the mesh we 'move' e.g. the one for which we search the transforamtion.

	requires normalize normals for vertices AND faces
	*/
	inline bool align(
			A2Grid &u,
			A2GridVert &uv,
			const Matrix44d &in,			// starting transformation that matches mov points to fix mesh
			Matrix44d &out,					// computed transformation
			std::vector<Point3d> &pfix,		// (red) corresponding vertices on src
			std::vector<Point3d> &nfix, 	// (red) corresponding normals on src
			std::vector<Point3d> &opmov,	// chosen vertices on trg (verdi) before the input transormation (Original Point Target)
			std::vector<Point3d> &onmov, 	// chosen normals on trg (verdi)
			Histogramf &h,
			Stat &as)
	{
		std::vector<char> beyondCntVec;  // flag vector to set the movverts that we should not use
		                                 // every time that a vertex is at a distance beyound max dist, its counter is incremented;
		                                 // movverts that has been discarded more than MaxCntDist times will not be considered anymore
		const int maxBeyondCnt = 3;
		std::vector< Point3d > movvert;
		std::vector< Point3d > movnorm;
		std::vector<Point3d> pmov;       // vertices chosen after the transformation
		status = SUCCESS;
		int tt0 = clock();

		out = in;

		int i;

		double cosAngleThr = cos(ap.MaxAngleRad);
		double startMinDist = ap.MinDistAbs;
		int tt1 = clock();
		int ttsearch = 0;
		int ttleast = 0;
		int nc = 0;
		as.clear();
		as.StartTime = clock();

		beyondCntVec.resize(mov->size(), 0);

		/**************** BEGIN ICP LOOP ****************/
		do {
			Stat::IterInfo ii;
			Box3d movbox;
			initMov(movvert, movnorm, movbox, out);
			h.SetRange(0.0f, float(startMinDist), 512, 2.5f);
			pfix.clear();
			nfix.clear();
			pmov.clear();
			opmov.clear();
			onmov.clear();
			int tts0 = clock();
			ii.MinDistAbs = startMinDist;
			int LocSampleNum = std::min(ap.SampleNum, int(movvert.size()));
			Box3d fixbox;
			if (u.Empty())
				fixbox = uv.bbox;
			else
				fixbox = u.bbox;
			for (i = 0; i < LocSampleNum; ++i) {
				if (beyondCntVec[i] < maxBeyondCnt) {
					if (fixbox.IsIn(movvert[i])) {
						double error = startMinDist;
						Point3d closestPoint, closestNormal;
						double maxd = startMinDist;
						ii.SampleTested++;
						if (u.Empty()) {  // using the point cloud grid 
							A2Mesh::VertexPointer vp = tri::GetClosestVertex(*fix, uv, movvert[i], maxd, error);
							if (error >= startMinDist) {
								ii.DistanceDiscarded++; ++beyondCntVec[i]; continue;
							}
							if (movnorm[i].dot(vp->N()) < cosAngleThr) {
								ii.AngleDiscarded++; continue;
							}
							closestPoint = vp->P();
							closestNormal = vp->N();
						}
						else {			// using the standard faces and grid
							A2Mesh::FacePointer f = vcg::tri::GetClosestFaceBase<vcg::AlignPair::A2Mesh, vcg::AlignPair::A2Grid >(*fix, u, movvert[i], maxd, error, closestPoint);
							if (error >= startMinDist) {
								ii.DistanceDiscarded++; ++beyondCntVec[i]; continue;
							}
							if (movnorm[i].dot(f->N()) < cosAngleThr) {
								ii.AngleDiscarded++; continue;
							}
							Point3d ip;
							InterpolationParameters<A2Face, double>(*f, f->N(), closestPoint, ip);
							const double IP_EPS = 0.00001;
							// If ip[i] == 0 it means that we are on the edge opposite to i
							if ((fabs(ip[0]) <= IP_EPS && f->IsB(1)) || (fabs(ip[1]) <= IP_EPS && f->IsB(2)) || (fabs(ip[2]) <= IP_EPS && f->IsB(0))){
								ii.BorderDiscarded++;  continue;
							}
							closestNormal = f->N();
						}
						// The sample was accepted. Store it.
						pmov.push_back(movvert[i]);
						opmov.push_back((*mov)[i].P());
						onmov.push_back((*mov)[i].N());
						nfix.push_back(closestNormal);
						pfix.push_back(closestPoint);
						h.Add(float(error));
						ii.SampleUsed++;
					}
					else {
						beyondCntVec[i] = maxBeyondCnt + 1;
					}
				}
			} // End for each pmov
			int tts1 = clock();
			printf("Found %d pairs\n",(int)pfix.size());
			if (!choosePoints(pfix, nfix, pmov, opmov, ap.PassHiFilter, h)) {
				if (int(pfix.size()) < ap.MinPointNum){
					status = TOO_FEW_POINTS;
					ii.Time = clock();
					as.I.push_back(ii);
					return false;
				}
			}
			Matrix44d newout;
			switch (ap.MatchMode){
			case AlignPair::Param::MMSimilarity:
				vcg::PointMatchingScale::computeRotoTranslationScalingMatchMatrix(newout, pfix, pmov);
				break;
			case AlignPair::Param::MMRigid:
				ComputeRigidMatchMatrix(pfix, pmov, newout);
				break;
			default:
				status = UNKNOWN_MODE;
				ii.Time = clock();
				as.I.push_back(ii);
				return false;
			}

			//double sum_before=0;
			//double sum_after=0;
			//for(unsigned int iii=0;iii<pfix.size();++iii){
			//	sum_before+=Distance(pfix[iii], out*OPmov[iii]);
			//	sum_after+=Distance(pfix[iii], newout*OPmov[iii]);
			//}
			//printf("Distance %f -> %f\n",sum_before/double(pfix.size()),sum_after/double(pfix.size()) ) ;

			// the following tuns will use as a initial transformation, the one that has been just found.
			// in the next loops the starting matrix will be this one.
			out = newout * out;

			assert(pfix.size() == pmov.size());
			int tts2 = clock();
			ttsearch += tts1 - tts0;
			ttleast += tts2 - tts1;
			ii.pcl50 = h.Percentile(.5);
			ii.pclhi = h.Percentile(float(ap.PassHiFilter));
			ii.AVG = h.Avg();
			ii.RMS = h.RMS();
			ii.StdDev = h.StandardDeviation();
			ii.Time = clock();
			as.I.push_back(ii);
			nc++;
			// The distance of the next points to be considered is lowered according to the <ReduceFactor> parameter.
			// We use 5 times the <ReduceFactor> percentile of the found points.
			if (ap.ReduceFactorPerc<1)
				startMinDist = std::max(ap.MinDistAbs*ap.MinMinDistPerc, std::min(startMinDist, 5.0*h.Percentile(float(ap.ReduceFactorPerc))));
			//as.dump(stderr);
		} while (
				nc <= ap.MaxIterNum &&
				h.Percentile(.5) > ap.TrgDistAbs &&
				(nc<ap.EndStepNum + 1 || !as.stable(ap.EndStepNum)) );

		/**************** END ICP LOOP ****************/
		int tt2 = clock();
		out[3][0] = 0; out[3][1] = 0; out[3][2] = 0; out[3][3] = 1;
		Matrix44d ResCopy = out;
		Point3d scv, shv, rtv, trv;
		Decompose(ResCopy, scv, shv, rtv, trv);
		if ((ap.MatchMode == vcg::AlignPair::Param::MMRigid) && (math::Abs(1 - scv[0])>ap.MaxScale || math::Abs(1 - scv[1]) > ap.MaxScale || math::Abs(1 - scv[2]) > ap.MaxScale)) {
			status = TOO_MUCH_SCALE;
			return false;
		}
		if (shv[0] > ap.MaxShear || shv[1] > ap.MaxShear || shv[2] > ap.MaxShear) {
			status = TOO_MUCH_SHEAR;
			return false;
		}
		printf("Grid %i %i %i - fn %i\n", u.siz[0], u.siz[1], u.siz[2], fix->fn);
		printf("Init %8.3f Loop %8.3f Search %8.3f least sqrt %8.3f\n",
		float(tt1 - tt0) / CLOCKS_PER_SEC, float(tt2 - tt1) / CLOCKS_PER_SEC,
		float(ttsearch) / CLOCKS_PER_SEC, float(ttleast) / CLOCKS_PER_SEC);

		return true;
	}

	/*
	 * Function called by Align at every cycle.
	 * It fills the <MovVert> and <MovNorm> vectors with the coordinates and normals
	 * taken from the the vertex vector "mov" of the mesh to move according to the
	 * matrix <In>.
	 * It computes also the new bounding box of the transformed vertices
	*/
	inline bool initMov(
			std::vector< Point3d > &movvert,
			std::vector< Point3d > &movnorm,
			Box3d &movbox,
			const Matrix44d &in	)
	{
		Point3d pp, nn;

		movvert.clear();
		movnorm.clear();
		movbox.SetNull();

		A2Mesh::VertexIterator vi;
		for (vi = mov->begin(); vi != mov->end(); vi++) {
			pp = in*(*vi).P();
			nn = in*Point3d((*vi).P() + (*vi).N()) - pp;
			nn.Normalize();
			movvert.push_back(pp);
			movnorm.push_back(nn);
			movbox.Add(pp);
		}
		return true;
	}

	static inline bool InitFixVert(
			A2Mesh *fm,
			AlignPair::Param &pp,
			A2GridVert &u,
			int preferredGridSize=0)
	{
		Box3d bb2 = fm->bbox;
		double minDist = pp.MinDistAbs*1.1;
		//the bbox of the grid should be enflated of the mindist used in the ICP search
		bb2.Offset(Point3d(minDist, minDist, minDist));

		u.SetBBox(bb2);

		//Inserisco la src nella griglia
		if (preferredGridSize == 0)
			preferredGridSize = int(fm->vert.size())*pp.UGExpansionFactor;
		u.Set(fm->vert.begin(), fm->vert.end());//, PreferredGridSize);
		printf("UG %i %i %i\n", u.siz[0], u.siz[1], u.siz[2]);
		return true;
	}

	static inline bool initFix(
			A2Mesh *fm,
			AlignPair::Param &pp,
			A2Grid &u,
			int preferredGridSize=0)
	{
		tri::InitFaceIMark(*fm);

		Box3d bb2 = fm->bbox;
		//	double MinDist= fm->bbox.Diag()*pp.MinDistPerc;
		double minDist = pp.MinDistAbs*1.1;
		//gonfio della distanza utente il BBox della seconda mesh
		bb2.Offset(Point3d(minDist, minDist, minDist));

		u.SetBBox(bb2);

		//Inserisco la src nella griglia
		if (preferredGridSize == 0)
			preferredGridSize = int(fm->face.size())*pp.UGExpansionFactor;
		u.Set(fm->face.begin(), fm->face.end(), preferredGridSize);
		printf("UG %i %i %i\n", u.siz[0], u.siz[1], u.siz[2]);
		return true;
	}

}; // end class

} // end namespace vcg

#endif
