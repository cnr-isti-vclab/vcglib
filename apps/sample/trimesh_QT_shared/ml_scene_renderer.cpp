#include "ml_scene_renderer.h"
#include "ml_thread_safe_memory_info.h"

MLThreadSafeGLMeshAttributesFeeder::MLThreadSafeGLMeshAttributesFeeder(CMesh& mesh,MLThreadSafeMemoryInfo& gpumeminfo,size_t perbatchtriangles)
	:GLMeshAttributesFeeder<CMesh>(mesh,gpumeminfo,perbatchtriangles),_lock(QReadWriteLock::Recursive)
{
}

void MLThreadSafeGLMeshAttributesFeeder::setPerBatchTriangles( size_t perbatchtriangles )
{
	QWriteLocker locker(&_lock);
	GLMeshAttributesFeeder<CMesh>::setPerBatchTriangles(perbatchtriangles);
}

size_t MLThreadSafeGLMeshAttributesFeeder::perBatchTriangles() const
{
	QReadLocker locker(&_lock);
	return GLMeshAttributesFeeder<CMesh>::perBatchTriangles();
	
}

bool MLThreadSafeGLMeshAttributesFeeder::renderedWithBO() const
{
	QReadLocker locker(&_lock);
	return GLMeshAttributesFeeder<CMesh>::renderedWithBO();
}

void MLThreadSafeGLMeshAttributesFeeder::update( int mask )
{
	QWriteLocker locker(&_lock);
	GLMeshAttributesFeeder<CMesh>::update(mask);
}

void MLThreadSafeGLMeshAttributesFeeder::drawWire(GLuint& vaohandlespecificperopenglcontext,vcg::GLFeedEnum::NORMAL_MODALITY nm,vcg::GLFeedEnum::COLOR_MODALITY cm )
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	drawTriangles(vaohandlespecificperopenglcontext,nm,cm,vcg::GLFeedEnum::TX_NONE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MLThreadSafeGLMeshAttributesFeeder::drawFlatWire(GLuint& vaohandlespecificperopenglcontext,vcg::GLFeedEnum::COLOR_MODALITY cm,vcg::GLFeedEnum::TEXTURE_MODALITY tm )
{
	glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT );
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1);
	drawTriangles(vaohandlespecificperopenglcontext,vcg::GLFeedEnum::NR_PERFACE,cm,tm);

	glDisable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
	glColor3f(.3f,.3f,.3f);

	drawWire(vaohandlespecificperopenglcontext,vcg::GLFeedEnum::NR_PERFACE,vcg::GLFeedEnum::CL_NONE);
	glPopAttrib();
}

void MLThreadSafeGLMeshAttributesFeeder::drawPoints(GLuint& vaohandlespecificperopenglcontext,vcg::GLFeedEnum::COLOR_MODALITY cm )
{
	QWriteLocker locker(&_lock);
	GLMeshAttributesFeeder<CMesh>::passPointsToOpenGL(vaohandlespecificperopenglcontext,vcg::GLFeedEnum::NR_PERVERT,cm);
}

void MLThreadSafeGLMeshAttributesFeeder::drawTriangles(GLuint& vaohandlespecificperopenglcontext,vcg::GLFeedEnum::NORMAL_MODALITY nm,vcg::GLFeedEnum::COLOR_MODALITY cm,vcg::GLFeedEnum::TEXTURE_MODALITY tm )
{
	QWriteLocker locker(&_lock);
	GLMeshAttributesFeeder<CMesh>::passTrianglesToOpenGL(vaohandlespecificperopenglcontext,nm,cm,tm);
}


