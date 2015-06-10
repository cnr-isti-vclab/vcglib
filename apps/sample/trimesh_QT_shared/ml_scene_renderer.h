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

#include <QObject>
#include <QMap>
#include <QReadWriteLock>

#include "mesh.h"

#include <GL/glew.h>
#include <wrap/gl/gl_mesh_attributes_feeder.h>



class MLThreadSafeMemoryInfo;


class MLThreadSafeGLMeshAttributesFeeder : public vcg::GLMeshAttributesFeeder<CMesh>
{
public:
	MLThreadSafeGLMeshAttributesFeeder(CMesh& mesh,MLThreadSafeMemoryInfo& gpumeminfo,size_t perbatchtriangles);
	~MLThreadSafeGLMeshAttributesFeeder() {};
	void setPerBatchTriangles(size_t perbatchtriangles);
	size_t perBatchTriangles() const;
	bool renderedWithBO() const;
	void update(int mask);
	void drawWire(GLuint& vaohandlespecificperopenglcontext,vcg::GLFeedEnum::NORMAL_MODALITY nm,vcg::GLFeedEnum::COLOR_MODALITY cm);

	void drawFlatWire(GLuint& vaohandlespecificperopenglcontext,vcg::GLFeedEnum::COLOR_MODALITY cm,vcg::GLFeedEnum::TEXTURE_MODALITY tm);

	void drawPoints(GLuint& vaohandlespecificperopenglcontext,vcg::GLFeedEnum::COLOR_MODALITY cm);

	void drawTriangles(GLuint& vaohandlespecificperopenglcontext,vcg::GLFeedEnum::NORMAL_MODALITY nm,vcg::GLFeedEnum::COLOR_MODALITY cm,vcg::GLFeedEnum::TEXTURE_MODALITY tm);
private:
	mutable QReadWriteLock _lock;
};

#endif
