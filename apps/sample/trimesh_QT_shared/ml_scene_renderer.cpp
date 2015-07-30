#include "ml_scene_renderer.h"
#include "ml_thread_safe_memory_info.h"
#include "glarea.h"

#include <vector>

MLThreadSafeGLMeshAttributesFeeder::MLThreadSafeGLMeshAttributesFeeder(CMeshO& mesh,MLThreadSafeMemoryInfo& gpumeminfo,size_t perbatchtriangles)
	:GLMeshAttributesFeeder<CMeshO>(mesh,gpumeminfo,perbatchtriangles),_lock(QReadWriteLock::Recursive)
{
}

void MLThreadSafeGLMeshAttributesFeeder::setPerBatchTriangles( size_t perbatchtriangles )
{
	QWriteLocker locker(&_lock);
	GLMeshAttributesFeeder<CMeshO>::setPerBatchPrimitives(perbatchtriangles);
}

size_t MLThreadSafeGLMeshAttributesFeeder::perBatchTriangles() const
{
	QReadLocker locker(&_lock);
	return GLMeshAttributesFeeder<CMeshO>::perBatchPrimitives();

}

bool MLThreadSafeGLMeshAttributesFeeder::renderedWithBO() const
{
	QReadLocker locker(&_lock);
	return GLMeshAttributesFeeder<CMeshO>::isPossibleToUseBORendering();
}

void MLThreadSafeGLMeshAttributesFeeder::meshAttributesUpdated( int mask )
{
	QWriteLocker locker(&_lock);
	GLMeshAttributesFeeder<CMeshO>::meshAttributesUpdated(mask);
}

vcg::GLFeederInfo::ReqAtts MLThreadSafeGLMeshAttributesFeeder::setupRequestedAttributes(const vcg::GLFeederInfo::ReqAtts& rq,bool& allocated )
{
	QWriteLocker locker(&_lock);
	return GLMeshAttributesFeeder<CMeshO>::setupRequestedAttributes(rq,allocated);
}


void MLThreadSafeGLMeshAttributesFeeder::drawWire(vcg::GLFeederInfo::ReqAtts& rq)
{
	glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	drawTriangles(rq);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glPopAttrib();
}

void MLThreadSafeGLMeshAttributesFeeder::drawFlatWire(vcg::GLFeederInfo::ReqAtts& rq)
{
	glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT );
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1);
	drawTriangles(rq);
	glDisable(GL_POLYGON_OFFSET_FILL);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);

	ReqAtts tmp = rq;
	tmp[ATT_VERTCOLOR] = false;
	tmp[ATT_FACECOLOR] = false;
	tmp[ATT_MESHCOLOR] = false;

	glColor3f(.3f,.3f,.3f);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	QReadLocker locker(&_lock);
	GLMeshAttributesFeeder<CMeshO>::drawTriangles(tmp,_textids.textId());
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glPopAttrib();
}

void MLThreadSafeGLMeshAttributesFeeder::drawPoints(vcg::GLFeederInfo::ReqAtts& rq)
{
	QReadLocker locker(&_lock);
	GLMeshAttributesFeeder<CMeshO>::drawPoints(rq);
}

void MLThreadSafeGLMeshAttributesFeeder::drawTriangles(vcg::GLFeederInfo::ReqAtts& rq)
{
	QReadLocker locker(&_lock);
    GLMeshAttributesFeeder<CMeshO>::drawTriangles(rq,_textids.textId());
}

void MLThreadSafeGLMeshAttributesFeeder::drawBBox(vcg::GLFeederInfo::ReqAtts& rq)
{
	QReadLocker locker(&_lock);
	
	vcg::Box3f& b = _mesh.bbox;
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	GLuint bbhandle;
	glGenBuffers(1,&bbhandle);
	std::vector<vcg::Point3f> bbox(12 * 2);

	//0
	bbox[0] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.min[2]);
	bbox[1] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.min[2]);

	//1
	bbox[2] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.min[2]);
	bbox[3] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.min[2]);

	//2
	bbox[4] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.min[2]);
	bbox[5] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.min[2]);

	//3
	bbox[6] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.min[2]);
	bbox[7] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.min[2]);

	//4
	bbox[8] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.min[2]);
	bbox[9] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.max[2]);

	//5
	bbox[10] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.max[2]);
	bbox[11] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.max[2]);

	//6
	bbox[12] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.max[2]);
	bbox[13] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.min[2]);

	//7
	bbox[14] = vcg::Point3f((float)b.max[0],(float)b.min[1],(float)b.max[2]);
	bbox[15] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.max[2]);

	//8
	bbox[16] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.max[2]);
	bbox[17] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.min[2]);
	
	//9
	bbox[18] = vcg::Point3f((float)b.max[0],(float)b.max[1],(float)b.max[2]);
	bbox[19] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.max[2]);

	//10
	bbox[20] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.max[2]);
	bbox[21] = vcg::Point3f((float)b.min[0],(float)b.min[1],(float)b.max[2]);

	//11
	bbox[22] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.max[2]);
	bbox[23] = vcg::Point3f((float)b.min[0],(float)b.max[1],(float)b.min[2]);

	glColor3f(1.0f,1.0f,1.0f);
	glBindBuffer(GL_ARRAY_BUFFER,bbhandle);
	glBufferData(GL_ARRAY_BUFFER, 12 * 2 * sizeof(vcg::Point3f), &(bbox[0]), GL_STATIC_DRAW);
	glVertexPointer(3,GL_FLOAT,0,0);
	glBindBuffer(GL_ARRAY_BUFFER,0);
	glEnableClientState(GL_VERTEX_ARRAY);
	glDrawArrays(GL_LINES,0,24);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDeleteBuffers(1,&bbhandle);
	glPopAttrib();
}


void MLThreadSafeGLMeshAttributesFeeder::deAllocateBO()
{
	QWriteLocker locker(&_lock);
	GLMeshAttributesFeeder<CMeshO>::buffersDeAllocationRequested();
}

MLThreadSafeGLMeshAttributesFeeder::MLThreadSafeTextureNamesContainer::MLThreadSafeTextureNamesContainer()
	:_tmid(),_lock(QReadWriteLock::Recursive)
{

}

MLThreadSafeGLMeshAttributesFeeder::MLThreadSafeTextureNamesContainer::~MLThreadSafeTextureNamesContainer()
{
	clear();
}

void MLThreadSafeGLMeshAttributesFeeder::MLThreadSafeTextureNamesContainer::push_back( GLuint textid )
{
	QWriteLocker locker(&_lock);
	_tmid.push_back(textid);
}

size_t MLThreadSafeGLMeshAttributesFeeder::MLThreadSafeTextureNamesContainer::size() const
{
	QReadLocker locker(&_lock);
	return _tmid.size(); 
}

bool MLThreadSafeGLMeshAttributesFeeder::MLThreadSafeTextureNamesContainer::empty() const
{
	QReadLocker locker(&_lock);
	return _tmid.empty(); 
}

void MLThreadSafeGLMeshAttributesFeeder::MLThreadSafeTextureNamesContainer::clear()
{
	QWriteLocker locker(&_lock);
	_tmid.clear();
}