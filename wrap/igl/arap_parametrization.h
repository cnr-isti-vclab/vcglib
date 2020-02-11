/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2017                                                \/)\/    *
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

/* Optimizes given UV-mapping with
 * [ARAP parametrization]
 * (minimizes area and angle distortions).
 *
 * Needs:
 * (-) per-vertex texture coords
 * (-) per-vertex flags to fix boundaries
 *     Fixed vertices are the flagged ones.
 *     By default: BORDER or SELECTED verts are fixed.
 *     (use fixedMask parameter to customize)
 *
 * Example of usage:
 *    MeshType m;
 *    ...
 *    vcg::tri::UpdateBounding<MeshType>::Box(m);
 *    vcg::tri::UpdateFlags<MeshType>::Clear(m);
 *    vcg::tri::UpdateFlags<MeshType>::VertexBorderFromNone(m);
 *    vcg::tri::OptimizeUV_ARAP(m);
 *
 */

#ifndef __VCG_IGL_ARAP_PARAMETRIZATION
#define __VCG_IGL_ARAP_PARAMETRIZATION

#include <cmath>

#include <igl/arap.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <wrap/igl/lscm_parametrization.h>

namespace vcg {
namespace tri {

template<class MeshType>
void OptimizeUV_ARAP(
        MeshType& m,
        unsigned int iterations = 100,
        unsigned int fixedMask = MeshType::VertexType::BORDER | MeshType::VertexType::SELECTED,
        bool generateInitialGuess = true)
{

	// check requirements
	vcg::tri::RequirePerVertexTexCoord(m);
	vcg::tri::RequirePerVertexFlags   (m);
	vcg::tri::RequireCompactness      (m);

	if (m.vert.size() <= 1 || m.face.size() == 0)
	{
		return;
	}

	// build fixed points data
	size_t nFixed = 0;
	if (fixedMask != 0)
	{
		for (size_t i=0; i<m.vert.size(); i++)
		{
			if (m.vert[i].Flags() & fixedMask) nFixed++;
		}
	}

	// all fixed, nothing to do? get out to avoid crashes
	if (nFixed == m.vert.size())
	{
		return;
	}

	if (generateInitialGuess)
	{
		// if not enough vertices are fixed, initialize manually fixed points
		// else initialize with the provided fixed values
		InitializeArapWithLSCM(m, (nFixed < 2) ? 0 : fixedMask);
	}

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::VectorXi b;
	Eigen::MatrixXd bc;
	Eigen::MatrixXd V_uv;
	vcg::tri::MeshToMatrix<MeshType>::GetTriMeshData(m, F, V);
	vcg::tri::MeshToMatrix<MeshType>::GetUVData(m, V_uv);

	b.resize(nFixed);
	bc.resize(nFixed,2);
	for (size_t i=0,k=0; i<m.vert.size(); i++)
	{
		if (m.vert[i].Flags() & fixedMask)
		{
			b(k) = i;
			bc(k,0) = m.vert[i].T().P()[0];
			bc(k,1) = m.vert[i].T().P()[1];
			k++;
		}
	}

	// Add dynamic regularization to avoid to specify boundary conditions
	::igl::ARAPData arap_data;
	arap_data.with_dynamics = true;
	arap_data.max_iter = iterations;

	// compute ARAP parametrization
	::igl::arap_precomputation(V, F, 2, b, arap_data);
	::igl::arap_solve(bc, arap_data, V_uv);

	// copy results back to mesh
	for (size_t i=0; i<m.vert.size(); i++)
	{
		m.vert[i].T().P()[0] = V_uv(i,0);
		m.vert[i].T().P()[1] = V_uv(i,1);
	}
}

template <class MeshType>
void InitializeArapWithLSCM(MeshType & m, unsigned int fixedMask = 0)
{
	typedef typename MeshType::ScalarType                          ScalarType;
	typedef typename MeshType::VertexType::TexCoordType::PointType TexPointType;
	typedef typename TexPointType::ScalarType                      TexScalarType;

	if (fixedMask == 0)
	{
		// automatically select 2 vertices to fix
		vcg::tri::UpdateFlags<MeshType>::Clear(m);

		int fixed0, fixed1 = -1;
		auto p = m.bbox.Center();
		ScalarType maxDist = -1;
		for (size_t i=0; i<m.vert.size(); i++)
		{
			// farthest point from the center
			const ScalarType dist =(m.vert[i].cP() - p).Norm();
			if (dist > maxDist)
			{
				fixed0 = i;
				maxDist = dist;
			}
		}
		maxDist = -1;
		p = m.vert[fixed0].cP();
		for (size_t i=0; i<m.vert.size(); i++)
		{
			// farthest point from the previous
			const ScalarType dist =(m.vert[i].cP() - p).Norm();
			if (dist > maxDist)
			{
				fixed1 = i;
				maxDist = dist;
			}
		}

		assert(fixed0 >= 0);
		assert(fixed1 >= 0);
		assert(fixed0 != fixed1);

		//then select them
		m.vert[fixed0].SetS();
		m.vert[fixed1].SetS();
		m.vert[fixed0].T().P() = TexPointType(0,0);
		m.vert[fixed1].T().P() = TexPointType(1,1);

		fixedMask = MeshType::VertexType::SELECTED;
	}

	vcg::tri::OptimizeUV_LSCM(m, fixedMask);

	// Rescale the parametrization to match the 3D area
	ScalarType meshArea2D = 0;
	ScalarType meshArea3D = 0;

	for (size_t i=0; i<m.face.size(); i++)
	{
		vcg::Triangle2<TexScalarType> t2(m.face[i].V(0)->T().P(),
		                                 m.face[i].V(1)->T().P(),
		                                 m.face[i].V(2)->T().P());
		meshArea2D += ScalarType(fabs(((t2.P(1) - t2.P(0)) ^ (t2.P(2) - t2.P(0)))/2));
		meshArea3D += vcg::DoubleArea(m.face[i])/2;
	}

	ScalarType scaleFact = std::sqrt(meshArea3D / meshArea2D);

	for (size_t i=0; i<m.vert.size(); i++)
	{
		TexPointType & UVCoord = m.vert[i].T().P();
		UVCoord *= scaleFact;
	}
}

}} // namespaces

#endif // __VCG_IGL_ARAP_PARAMETRIZATION
