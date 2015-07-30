/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
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

#ifndef __ML_SCENE_RENDERER_H
#define __ML_SCENE_RENDERER_H

#include <GL/glew.h>
#include "mesh.h"
#include <wrap/gl/gl_mesh_attributes_feeder.h>

#include <QObject>
#include <QMap>
#include <QReadWriteLock>



class MLThreadSafeMemoryInfo;


class MLThreadSafeGLMeshAttributesFeeder : public vcg::GLMeshAttributesFeeder<CMeshO>
{
public:
	struct MLThreadSafeTextureNamesContainer
	{
		MLThreadSafeTextureNamesContainer();
		~MLThreadSafeTextureNamesContainer();

		void push_back(GLuint textid);
		size_t size() const;
		bool empty() const;
		void clear();
		GLuint& operator[](size_t ii) {return _tmid[ii];};
		inline std::vector<GLuint>& textId() {return _tmid;};
	private:
		std::vector<GLuint> _tmid;
		mutable QReadWriteLock _lock;
	};


	MLThreadSafeGLMeshAttributesFeeder(CMeshO& mesh,MLThreadSafeMemoryInfo& gpumeminfo,size_t perbatchtriangles);
	~MLThreadSafeGLMeshAttributesFeeder() {};

	void setPerBatchTriangles(size_t perbatchtriangles);
	
	size_t perBatchTriangles() const;
	
	bool renderedWithBO() const;
	
	GLuint bufferObjectHandle() const;

	void meshAttributesUpdated(int mask);

	vcg::GLFeederInfo::ReqAtts setupRequestedAttributes(const vcg::GLFeederInfo::ReqAtts& rq,bool& allocated);

	void deAllocateBO();
	
	void drawWire(vcg::GLFeederInfo::ReqAtts& rq);

	void drawFlatWire(vcg::GLFeederInfo::ReqAtts& rq);

	void drawPoints(vcg::GLFeederInfo::ReqAtts& rq);

	void drawTriangles(vcg::GLFeederInfo::ReqAtts& rq);

	void drawBBox(vcg::GLFeederInfo::ReqAtts& rq);

	inline CMeshO& mesh() {return _mesh;}

	inline MLThreadSafeTextureNamesContainer& textureIDContainer() {return _textids;}

	void buffersDeAllocationRequested();
private:
	mutable QReadWriteLock _lock;
	MLThreadSafeTextureNamesContainer _textids;
};

#endif
