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

#ifndef __PICK______H
#define __PICK______H

// Assumes OpenGL already included

#include <vector>
#include <algorithm>

#include <vcg/complex/complex.h>

namespace vcg{

template <class MESH_TYPE>
class GLPickTri
{
public:
	typedef typename MESH_TYPE::ScalarType     ScalarType;
	typedef typename MESH_TYPE::CoordType      CoordType;
	typedef typename MESH_TYPE::FaceIterator   FaceIterator;
	typedef typename MESH_TYPE::VertexIterator VertexIterator;
	typedef typename MESH_TYPE::FacePointer    FacePointer;
	typedef typename MESH_TYPE::VertexPointer  VertexPointer;
	typedef typename MESH_TYPE::VertexType     VertexType;

	static CoordType glProject(const Eigen::Matrix<ScalarType,4,4> &M, const ScalarType * viewport, const CoordType &p)
	{
		const ScalarType vx=viewport[0];
		const ScalarType vy=viewport[1];
		const ScalarType vw2=viewport[2]/ScalarType(2.0);
		const ScalarType vh2=viewport[3]/ScalarType(2.0);
		Eigen::Matrix<ScalarType,4,1> vp(p[0],p[1],p[2],ScalarType(1.0));
		Eigen::Matrix<ScalarType,4,1> vpp = M*vp;
		Eigen::Matrix<ScalarType,4,1> ndc = vpp/vpp[3];

		CoordType sc(vw2*ndc[0] + vx+vw2,
		             vh2*ndc[1] + vy+vh2,
		             ndc[2]);

		return sc;
	}

	static void FillProjectedVector(MESH_TYPE &m, std::vector<CoordType> &pVec, const Eigen::Matrix<ScalarType,4,4> &M, const ScalarType * viewportF)
	{
		pVec.resize(m.vert.size());
		for(size_t i=0;i<m.vert.size();++i) if(!m.vert[i].IsD())
		{
			pVec[i] = GLPickTri<MESH_TYPE>::glProject(M, viewportF,CoordType::Construct(m.vert[i].P()));
		}
	}

	static void glGetMatrixAndViewport(Eigen::Matrix<ScalarType,4,4> &M, ScalarType *viewportF)
	{
		Eigen::Matrix4d mp,mm;

		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT,viewport);
		for(int i=0;i<4;++i) viewportF[i]=viewport[i];

		glGetDoublev(GL_PROJECTION_MATRIX, mp.data());
		glGetDoublev(GL_MODELVIEW_MATRIX,  mm.data());

		M = (mp*mm).cast<ScalarType>();
	}

	// compute a bbox in Device Coordinate (with the z without the near far normalization and ranged in -1 1)
	static Box3<ScalarType> ComputeDCBox(int x, int y, int width, int height)
	{
		Box3<ScalarType> bb;
		bb.SetNull();
		bb.Add(CoordType(x-width/ScalarType(2.0),y-height/ScalarType(2.0),ScalarType(-1.0)));
		bb.Add(CoordType(x+width/ScalarType(2.0),y+height/ScalarType(2.0), ScalarType(1.0)));
		return bb;
	}

	static bool PickClosestFace(int x, int y, MESH_TYPE &m, FacePointer &fp,int width=4, int height=4)
	{
		Eigen::Matrix<ScalarType,4,4> M;
		ScalarType viewportF[4];
		glGetMatrixAndViewport(M,viewportF);
		Box3<ScalarType> reg=ComputeDCBox(x,y,width,height);

		ScalarType bzmin = std::numeric_limits<ScalarType>::max();
		fp=0;
		for(size_t i=0;i<m.face.size();++i) if(!m.face[i].IsD())
		{
			CoordType bary = vcg::Barycenter(m.face[i]);
			CoordType bz = glProject(M, viewportF,bary);

			if(bz[2]<bzmin && reg.IsIn(bz))
			{
				bzmin=bz[2];
				fp = &m.face[i];
			}
		}
		return fp!=0;
	}

	static bool PickClosestVert(int x, int y, MESH_TYPE &m, VertexPointer &vp,int width=4, int height=4)
	{
		Eigen::Matrix<ScalarType,4,4> M;
		ScalarType viewportF[4];
		glGetMatrixAndViewport(M,viewportF);
		ScalarType bzmin = std::numeric_limits<ScalarType>::max();
		vp=0;

		Box3<ScalarType> reg=ComputeDCBox(x,y,width,height);

		for(size_t i=0;i<m.vert.size();++i) if(!m.vert[i].IsD())
		{
			CoordType bz = glProject(M, viewportF,m.vert[i].P());
			if(bz[2]<bzmin && reg.IsIn(bz))
			{
				bzmin=bz[2];
				vp = &m.vert[i];
			}
		}
		return vp!=0;
	}

	static int PickVert(int x, int y, MESH_TYPE &m, std::vector<VertexPointer> &result, int width=4, int height=4)
	{
		result.clear();
		static Eigen::Matrix<ScalarType,4,4> lastM;
		static MESH_TYPE *lastm=0;
		static std::vector<CoordType> pVec;

		Eigen::Matrix<ScalarType,4,4> M;
		ScalarType viewportF[4];
		glGetMatrixAndViewport(M,viewportF);

		Box3<ScalarType> reg =ComputeDCBox(x,y,width,height);

		if ((M != lastM) || (&m != lastm) || (pVec.size() != m.VN()))
		{
			FillProjectedVector(m,pVec,M,viewportF);
			lastM = M;
			lastm = &m;
		}

		for(size_t i=0;i<m.vert.size();++i) if(!m.vert[i].IsD())
		{
			if(reg.IsIn(pVec[i]))
				result.push_back(&m.vert[i]);
		}
		return result.size();
	}

	static int PickFace(int x, int y, MESH_TYPE &m, std::vector<FacePointer> &result, int width = 4, int height = 4)
	{
		static Eigen::Matrix<ScalarType, 4, 4> lastM;
		static MESH_TYPE *lastm = 0;
		static std::vector<CoordType> pVec;

		ScalarType viewportF[4];
		Eigen::Matrix<ScalarType, 4, 4> M;
		glGetMatrixAndViewport(M, viewportF);
		result.clear();
		Box3<ScalarType> reg;
		reg.Add(CoordType(x - width / ScalarType(2.0), y - height / ScalarType(2.0), ScalarType(-1.0)));
		reg.Add(CoordType(x + width / ScalarType(2.0), y + height / ScalarType(2.0), ScalarType(1.0)));

		if ((M != lastM) || (&m != lastm) || (pVec.size() != m.VN()))
		{
			FillProjectedVector(m, pVec, M, viewportF);
			lastM = M;
			lastm = &m;
		}

		for (size_t i = 0; i < m.face.size(); ++i)
		{
			if (!m.face[i].IsD())
			{
				const CoordType &p0 = pVec[tri::Index(m, m.face[i].V(0))];
				const CoordType &p1 = pVec[tri::Index(m, m.face[i].V(1))];
				const CoordType &p2 = pVec[tri::Index(m, m.face[i].V(2))];

				if (!(abs(p0[2]) > 1 || abs(p1[2]) > 1 || abs(p2[2]) > 1) && IntersectionTriangleBox(reg, p0, p1, p2))
					result.push_back(&m.face[i]);
			}
		}
		return result.size();
	}

	// Same of above but it also assumes that you want only visible faces.
	// Visibility is computed according to the current depth buffer.
	static int PickVisibleFace(int x, int y, MESH_TYPE &m, std::vector<FacePointer> &resultZ, int width=4, int height=4)
	{
		ScalarType vp[4];
		Eigen::Matrix<ScalarType,4,4> M;
		glGetMatrixAndViewport(M,vp);

		int screenW = (int)(vp[2]-vp[0]);
		int screenH = (int)(vp[3]-vp[1]);

		GLfloat *buffer = new GLfloat[screenW*screenH];

		//glReadPixels does NOT accept GL_DOUBLE
		//GLenum err = glGetError();
		glReadPixels(vp[0],vp[1],vp[2],vp[3],GL_DEPTH_COMPONENT,GL_FLOAT,buffer);
		//err = glGetError();
		std::vector<FacePointer> result;
		PickFace(x,y,m,result,width,height);
		ScalarType LocalEpsilon(ScalarType(0.001));
		for(size_t i =0;i<result.size();++i)
		{
			CoordType p = glProject(M,vp,CoordType::Construct(Barycenter(*(result[i])))) ;
			if(p[0] >=0 && p[0]<screenW && p[1] >=0 && p[1]<screenH)
			{
				ScalarType bufZ(buffer[int(p[0])+int(p[1])*screenW]);
				if(bufZ + LocalEpsilon >= ScalarType(p[2]+1.0)/2.0)
					resultZ.push_back(result[i]);
			}
		}

		delete [] buffer;
		return resultZ.size();
	}

};

} // end namespace vcg

#endif
