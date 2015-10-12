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

#ifndef __VCG_GL_MESH_ATTRIBUTES_FEEDER
#define __VCG_GL_MESH_ATTRIBUTES_FEEDER

#include <queue>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <limits>

//#include <GL/glew.h>
#include <wrap/gl/space.h>
#include <wrap/gl/math.h>
#include <vcg/space/color4.h>
#include<wrap/system/memory_info.h>


namespace vcg
{
    struct GLFeederInfo
    {
        struct GLFeederException : public std::exception
        {
            GLFeederException(const char* text)
                :std::exception(),_text(text) {}

            ~GLFeederException() throw() {}
            inline const char* what() const throw() {return _text.c_str();}
        private:
            std::string _text;
        };

        enum ATT_NAMES
        {
            ATT_ALL = -1,
            ATT_VERTPOSITION = 0,
            ATT_VERTNORMAL = 1,
            ATT_FACENORMAL = 2,
            ATT_VERTCOLOR = 3,
            ATT_FACECOLOR = 4,
            ATT_MESHCOLOR = 5,
            ATT_VERTTEXTURE = 6,
            ATT_WEDGETEXTURE = 7,
            ATT_VERTINDEX = 8,
            ATT_NAMES_ARITY = 9
        };

        enum PRIMITIVE_MODALITY
        {
            PR_NONE = 0,
            PR_POINTS = 1,
            PR_TRIANGLES = 2,
            PR_QUADS = 3,
            PR_NAMES_ARITY = 4
        };

        typedef unsigned int ATT_BIT_MASK;

        static ATT_BIT_MASK attBitMask(ATT_NAMES attname)
        {
            if ((attname == ATT_ALL) || (attname == ATT_VERTINDEX))
                return 0xffffffff;

            if ((attname >= ATT_VERTPOSITION) && (attname <= ATT_WEDGETEXTURE))
            {
                static const ATT_BIT_MASK res[] = {
                    0x00000000, /*attVertPosition()*/
                    0x00000001, /*attVertNormal()*/
                    0x00000002, /*attFaceNormal()*/
                    0x00000004, /*attVertColor()*/
                    0x00000008, /*attFaceColor()*/
                    0x00000010, /*attVertTexture()*/
                    0x00000020}; /*attWedgeTexture()*/
                    return res[attname];
            }
            else
                throw GLFeederException("Out of range value\n");
            return 0;
        }

        class ReqAtts
        {
        public:
            ReqAtts()
            {
                reset();
            }


            const bool& operator[](ATT_NAMES att) const
            {
                size_t ii = static_cast<size_t>(att);
                if (ii > ATT_VERTINDEX)
                    throw GLFeederException("Out of range value\n");
                return _atts[ii];
            }

            bool& operator[](ATT_NAMES att)
            {
                size_t ii = static_cast<size_t>(att);
                if (ii > ATT_VERTINDEX)
                    throw GLFeederException("Out of range value\n");
                return _atts[ii];
            }

            virtual void reset()
            {
                for(size_t ii = 0;ii < _attssize;++ii)
                    _atts[ii] = false;
                _pm = PR_NONE;
            }

            inline PRIMITIVE_MODALITY& primitiveModality()
            {
                return _pm;
            }

            inline PRIMITIVE_MODALITY primitiveModality() const
            {
                return _pm;
            }

            inline static size_t possibleAttributesNumber()
            {
                return _attssize;
            }

            static ReqAtts setUnion(const ReqAtts& a,const ReqAtts& b)
            {
                ReqAtts res;
                for(size_t ii = 0; ii < ReqAtts::possibleAttributesNumber();++ii)
                {
                    ATT_NAMES name = static_cast<ATT_NAMES>(ii);
                    res[name] = a[name] || b[name];
                }

                res.primitiveModality() = a.primitiveModality();
                if ((unsigned int) res.primitiveModality() <= (unsigned int) b.primitiveModality())
                    res.primitiveModality() = b.primitiveModality();
                return res;
            }

            static ReqAtts setComplement(const ReqAtts& a,const ReqAtts& b)
            {
                /*TRUTH TABLE*/
                //this[ATT_NAMES] | rq[ATT_NAMES] | res
                //    true        |     true      | false
                //    true        |     false     | true
                //    false       |     true      | false
                //    false       |     false     | false

                ReqAtts res = a;
                for(size_t ii = 0; ii < ReqAtts::possibleAttributesNumber();++ii)
                {
                    ATT_NAMES name = static_cast<ATT_NAMES>(ii);
                    if (res[name])
                        res[name] = !(b[name]);
                }

                res.primitiveModality() = b.primitiveModality();
                return res;
            }

            static bool isReplicatedPipeline(const ReqAtts& rqatt)
            {
                return (rqatt[ATT_FACENORMAL] || rqatt[ATT_FACECOLOR] || rqatt[ATT_WEDGETEXTURE]);
            }

            static bool isVertexIndexingRequired(const ReqAtts& rqatt)
            {
                PRIMITIVE_MODALITY pm = rqatt.primitiveModality();
                return (!isReplicatedPipeline(rqatt) && (pm != PR_POINTS) && (pm != PR_NONE));
            }

            template<typename MESHTYPE>
            static void computeARequestedAttributesSetCompatibleWithMesh(ReqAtts& rqatt,const MESHTYPE& mesh)
            {
                if (mesh.VN() == 0)
                {
                    rqatt.reset();
                    return;
                }

                rqatt[ATT_VERTPOSITION] = true;
                rqatt[ATT_VERTNORMAL] = rqatt[ATT_VERTNORMAL] && vcg::tri::HasPerVertexNormal(mesh);
                rqatt[ATT_FACENORMAL] = rqatt[ATT_FACENORMAL] && vcg::tri::HasPerFaceNormal(mesh);
                rqatt[ATT_VERTCOLOR] = rqatt[ATT_VERTCOLOR] && vcg::tri::HasPerVertexColor(mesh);
                rqatt[ATT_FACECOLOR] = rqatt[ATT_FACECOLOR] && vcg::tri::HasPerFaceColor(mesh);
                rqatt[ATT_MESHCOLOR] = rqatt[ATT_MESHCOLOR];
                rqatt[ATT_VERTTEXTURE] = rqatt[ATT_VERTTEXTURE] && vcg::tri::HasPerVertexTexCoord(mesh);
                rqatt[ATT_WEDGETEXTURE] = rqatt[ATT_WEDGETEXTURE] && vcg::tri::HasPerWedgeTexCoord(mesh);
                rqatt[ATT_VERTINDEX] = isVertexIndexingRequired(rqatt);
            }

        protected:
            static const size_t _attssize = ATT_NAMES_ARITY;
            bool _atts[_attssize];
            PRIMITIVE_MODALITY _pm;
        };
    };

    //WARNING! member functions of this class should be called by the host application using concurrency
    template <typename MESHTYPE>
    class GLMeshAttributesFeeder : public GLFeederInfo
    {
    public:
        GLMeshAttributesFeeder(/*const*/ MESHTYPE& mesh,MemoryInfo& meminfo, size_t perbatchprimitives)
            :_mesh(mesh),_gpumeminfo(meminfo),_bo(ATT_NAMES_ARITY,NULL),_currallocatedboatt(),_lastfeedingusedreplicatedpipeline(false),_perbatchprim(perbatchprimitives),_chunkmap(),_borendering(false),_rendermodinitialized(false)
        {
            _bo[GLFeederInfo::ATT_VERTPOSITION] = new GLBufferObject(3,GL_FLOAT,GL_VERTEX_ARRAY,GL_ARRAY_BUFFER);
            _bo[GLFeederInfo::ATT_VERTNORMAL] = new GLBufferObject(3,GL_FLOAT,GL_NORMAL_ARRAY,GL_ARRAY_BUFFER);
            _bo[GLFeederInfo::ATT_FACENORMAL] = new GLBufferObject(3,GL_FLOAT,GL_NORMAL_ARRAY,GL_ARRAY_BUFFER);
            _bo[GLFeederInfo::ATT_VERTCOLOR] = new GLBufferObject(4,GL_UNSIGNED_BYTE,GL_COLOR_ARRAY,GL_ARRAY_BUFFER);
            _bo[GLFeederInfo::ATT_FACECOLOR] = new GLBufferObject(4,GL_UNSIGNED_BYTE,GL_COLOR_ARRAY,GL_ARRAY_BUFFER);
            /*MESHCOLOR has not a buffer object associated with it. It's just a call to glColor3f. it's anyway added to the _bo arrays for sake of coherence*/
            _bo[GLFeederInfo::ATT_MESHCOLOR] = NULL;
            _bo[GLFeederInfo::ATT_VERTTEXTURE] = new GLBufferObject(2,GL_FLOAT,GL_TEXTURE_COORD_ARRAY,GL_ARRAY_BUFFER);
            _bo[GLFeederInfo::ATT_WEDGETEXTURE] = new GLBufferObject(2,GL_FLOAT,GL_TEXTURE_COORD_ARRAY,GL_ARRAY_BUFFER);
            _bo[GLFeederInfo::ATT_VERTINDEX] = new GLBufferObject(3,GL_UNSIGNED_INT,GL_ELEMENT_ARRAY_BUFFER);
        }

        ~GLMeshAttributesFeeder()
        {
            for(size_t ii = 0;ii < _bo.size();++ii)
                delete _bo[ii];
            _bo.clear();
        }

        //void meshAttributesUpdated(int mask)
        //{
        //    if ((_mesh.VN() != _mesh.vert.size()) || (_mesh.FN() != _mesh.face.size()))
        //    {
        //        throw GLFeederException("The current mesh contains vertices/faces marked as deleted!\nPlease, call the compact vectors function in order to properly remove them\n");
        //        return;
        //    }

        //    if (!_rendermodinitialized)
        //    {
        //        throw GLFeederException("Required attributes for rendering were not properly initialized\nInvokation of setAttributesToBeRendered function is strictly required!\n");
        //        return;
        //    }

        //    bool boupdatedrequired = false;
        //    if (((mask & attBitMask(GLFeederInfo::ATT_VERTPOSITION)) || (mask & attBitMask(ATT_ALL))) && (_currallocatedboatt[ATT_VERTPOSITION]))
        //    {
        //        boupdatedrequired = true;
        //        _bo[GLFeederInfo::ATT_VERTPOSITION]->_isvalid = false;
        //    }

        //    if (((mask & attBitMask(GLFeederInfo::ATT_VERTNORMAL)) || (mask & attBitMask(ATT_ALL))) &&  (_currallocatedboatt[ATT_VERTNORMAL] && vcg::tri::HasPerVertexNormal(_mesh)))
        //    {
        //        boupdatedrequired = true;
        //        _bo[GLFeederInfo::ATT_VERTNORMAL]->_isvalid = false;
        //    }

        //    if (((mask & attBitMask(GLFeederInfo::ATT_FACENORMAL)) || (mask & attBitMask(ATT_ALL))) && (_currallocatedboatt[ATT_FACENORMAL] && vcg::tri::HasPerFaceNormal(_mesh)))
        //    {
        //        boupdatedrequired = true;
        //        _bo[GLFeederInfo::ATT_FACENORMAL]->_isvalid = false;
        //    }

        //    if (((mask & attBitMask(GLFeederInfo::ATT_VERTCOLOR)) || (mask & attBitMask(ATT_ALL))) && (_currallocatedboatt[ATT_VERTCOLOR] && vcg::tri::HasPerVertexColor(_mesh)))
        //    {
        //        boupdatedrequired = true;
        //        _bo[GLFeederInfo::ATT_VERTCOLOR]->_isvalid = false;
        //    }

        //    if (((mask & attBitMask(GLFeederInfo::ATT_FACECOLOR)) || (mask & attBitMask(ATT_ALL))) && (_currallocatedboatt[ATT_FACECOLOR] && vcg::tri::HasPerFaceColor(_mesh)))
        //    {
        //        boupdatedrequired = true;
        //        _bo[GLFeederInfo::ATT_FACECOLOR]->_isvalid = false;
        //    }


        //    if (((mask & attBitMask(GLFeederInfo::ATT_VERTTEXTURE)) || (mask & attBitMask(ATT_ALL))) && (_currallocatedboatt[ATT_VERTTEXTURE] && vcg::tri::HasPerVertexTexCoord(_mesh)))
        //    {
        //        boupdatedrequired = true;
        //        _bo[GLFeederInfo::ATT_VERTTEXTURE]->_isvalid = false;
        //    }

        //    if (((mask & attBitMask(GLFeederInfo::ATT_WEDGETEXTURE)) || (mask & attBitMask(ATT_ALL))) && (_currallocatedboatt[ATT_WEDGETEXTURE] && vcg::tri::HasPerWedgeTexCoord(_mesh)))
        //    {
        //        boupdatedrequired = true;
        //        _bo[GLFeederInfo::ATT_WEDGETEXTURE]->_isvalid = false;
        //    }

        //    if (mask & attBitMask(ATT_ALL))
        //    {
        //        boupdatedrequired = true;
        //        _bo[GLFeederInfo::ATT_VERTINDEX]->_isvalid = false;
        //    }

        //    if (boupdatedrequired)
        //        tryToAllocateAndCopyAttributesInBO();
        //}

        ReqAtts setupRequestedAttributes(const ReqAtts& rq,bool& allocated)
        {

            if (!_rendermodinitialized)
                _rendermodinitialized = true;

            try
            {
                ReqAtts tmp = rq;
                ReqAtts::computeARequestedAttributesSetCompatibleWithMesh(tmp,_mesh);
                tmp = ReqAtts::setUnion(_currallocatedboatt,tmp);
               
                allocated = tryToAllocateAndCopyAttributesInBO(tmp);
                return tmp;
            }
            catch (GLFeederException& e)
            {
                return ReqAtts();
            }
            return ReqAtts();
        }

        ReqAtts removeRequestedAttributes(const ReqAtts& rq)
        {
            return _currallocatedboatt = ReqAtts::setComplement(rq,_currallocatedboatt);
        }

        void invalidateRequestedAttributes(const ReqAtts& rq)
        {
            size_t ii = 0;
            for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
            {
                ATT_NAMES boname = static_cast<ATT_NAMES>(ii);
                if ((*it != NULL) && (rq[boname]))
                    (*it)->_isvalid = false;
                ++ii;
            }
        }

        void buffersDeAllocationRequested()
        {
            size_t ii = 0;
            for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
            {
                ATT_NAMES boname = static_cast<ATT_NAMES>(ii);
                if (*it != NULL)
                    bufferDeAllocationRequested(boname);
                ++ii;
            }
        }

        void buffersDeAllocationRequested(const ReqAtts& rq)
        {
            size_t ii = 0;
            for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
            {
                ATT_NAMES boname = static_cast<ATT_NAMES>(ii);
                if ((*it != NULL) && (rq[boname]))
                    bufferDeAllocationRequested(boname);
                ++ii;
            }
        }

        void draw(const ReqAtts& rq,const std::vector<GLuint> textid = std::vector<GLuint>())
        {

            if (isPossibleToUseBORendering())
            {
                switch(rq.primitiveModality())
                {
                case(PR_TRIANGLES):
                    drawTriangles(rq,textid);
                    break;
                case(PR_POINTS):
                    drawPoints(rq);
                    break;
                case (PR_QUADS):
                    break;
                default:
                    break;
                }
            }
            else
                immediateModeRendering(rq,textid);
        }

        void setPerBatchPrimitives(size_t perbatchprimitives)
        {
            _perbatchprim = perbatchprimitives;
        }

        size_t perBatchPrimitives() const
        {
            return _perbatchprim;
        }

        bool isPossibleToUseBORendering() const
        {
            return _borendering;
        }

 

        //void invalidateRequestedAttributes(ReqAtts& rq)
        //{
        //    size_t ii = 0;
        //    for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
        //    {
        //        ATT_NAMES boname = static_cast<ATT_NAMES>(ii);
        //        if (((*it) != NULL) && (rq[boname]))
        //            (*it)->_isvalid = false;
        //        ++ii;
        //    }
        //    //_currallocatedboatt = vcg::GLFeederInfo::ReqAtts::setComplement(_currallocatedboatt,rq);
        //}

    protected:
        struct GLBufferObject
        {
            GLBufferObject(size_t components,GLenum gltype,GLenum clientstatetag,GLenum target)
                :_size(0),_components(components),_isvalid(false),_gltype(gltype),_clientstatetag(clientstatetag),_target(target),_bohandle(0)
            {
            }

            GLBufferObject(size_t components,GLenum gltype,GLenum target)
                :_size(0),_components(components),_isvalid(false),_gltype(gltype),_clientstatetag(),_target(target),_bohandle(0)
            {
            }

            size_t getSizeOfGLType() const
            {
                switch(_gltype)
                {
                case(GL_FLOAT):
                    return sizeof(GLfloat);
                case(GL_INT):
                    return sizeof(GLint);
                case(GL_UNSIGNED_INT):
                    return sizeof(GLuint);
                case(GL_UNSIGNED_BYTE):
                    return sizeof(GLubyte);
                }
                return 0;
            }

            size_t _size;
            const size_t _components;
            bool _isvalid;
            const GLenum _gltype;
            const GLenum _target;

            /*WARNING!!!!!!!!!!!!!!!!! In openGL INDEX BO doesn't require to be enabled/disabled so has NOT a valid tag associated.
            In this case the client state tag remains not initialized and it's not meaningful */
            const GLenum _clientstatetag;
            /**********************************************************************************/

            GLuint _bohandle;
        };

        const GLBufferObject& getBufferObjectInfo(ATT_NAMES boname) const
        {
            return _bo[boname];
        }

        void bufferDeAllocationRequested(const ATT_NAMES att)
        {
            size_t ind = static_cast<size_t>(att);
            if ((ind < 0) || (ind >= _bo.size()))
                return;
            GLBufferObject* bobj = _bo[att];
            if (bobj == NULL)
                return;
           
            if ((att != ATT_VERTINDEX ) && (ATT_MESHCOLOR))
                glDisableClientState(bobj->_clientstatetag);

            glDeleteBuffers(1,&(bobj->_bohandle));
            bobj->_bohandle = 0;
            if (bobj->_size > 0)
                //we don't use dim cause dim is the value that is going to be allocated, instead use (*it)->_size * (*it)->getSizeOfGLType() is the value already in the buffer
                _gpumeminfo.releasedMemory(bobj->_size * bobj->getSizeOfGLType());
            bobj->_isvalid = false;
            bobj->_size = 0;
            _currallocatedboatt[att] = false; 
        }

        bool buffersAllocationFunction(const ReqAtts& req,std::vector<bool>& attributestobeupdated)
        {
            bool replicated = ReqAtts::isReplicatedPipeline(req);
            attributestobeupdated.clear();
            attributestobeupdated.resize(_bo.size());
            std::ptrdiff_t bomemoryrequiredbymesh = bufferObjectsMemoryRequired(req);
            bool generateindex = ReqAtts::isVertexIndexingRequired(req);
            unsigned int ii = 0;
            for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
            {
                ATT_NAMES boname = static_cast<ATT_NAMES>(ii);
                size_t sz = boExpectedSize(boname,replicated,generateindex);
                size_t dim = boExpectedDimension(boname,replicated,generateindex);

                if (((*it) != NULL) && (
                    /*a mesh attributes has been updated and the number of already allocated bo cells to contain the attribute values doesn't suit anymore
                    (i.e. if i change just the vertex positions without changing the vertex numbers i have not to reallocate the previous vertposition bo, it's just sufficient to update the vertex coordinates)*/
                    ((!(*it)->_isvalid) && (sz != (*it)->_size)) ||
                    /*if _lastfeedingusedreplicatedpipeline == false means that maybe there are valid per vertex attribute buffer objects that MUST be reallocated anyway cause we have to switch to the replicated attributes pipeline*/
                    (replicated && !_lastfeedingusedreplicatedpipeline && (isPerVertexAttribute(boname) || (boname == GLFeederInfo::ATT_VERTINDEX))) ||
                    /*we switched back from the replicated pipeline to the indexed one. All the bos have to be regenerated*/
                    (!replicated && _lastfeedingusedreplicatedpipeline) ||
                    /*the buffer object is valid but for same reason the number of cells of the bo don't suit anymore the required size. we have to reallocate the buffer object*/
                    (((*it)->_isvalid) && (sz != (*it)->_size))
                    //||
                    //the buffer is valid, but the attribute is not required to be displayed
                    /*(((*it)->_isvalid) && !isAttributeRequiredToBeDisplayed(boname)))*/
                    ))
                {
                    bufferDeAllocationRequested(boname); 
                }
                /*there are already mesh attributes properly allocated in memory, we don't need to allocate them again. there could be invalid values attributes but with properly memory space already allocated in memory
                (i.e. i changed the pervertex colors but the vertex numbers remained constant)*/
                if ((*it != NULL) && (/*((*it)->_isvalid) &&*/ (sz == (*it)->_size)))
                    bomemoryrequiredbymesh -= dim;
                ++ii;
            }

            if (!_gpumeminfo.isAdditionalMemoryAvailable(bomemoryrequiredbymesh))
            {
                std::cout << "no additional memory available!!! memory required: " << bomemoryrequiredbymesh << std::endl;
                ii = 0;
                for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
                {
                    ATT_NAMES boname = static_cast<ATT_NAMES>(ii);
                    size_t sz(boExpectedSize(boname,replicated,generateindex));
                    //there are already valid mesh attributes properly allocated in memory but there is not enough gpu memory for the remaining mesh.
                    //we have to deallocate the previously allocated mesh attributes
                    if ((*it != NULL) && ((sz == (*it)->_size)))
                    {
                        ////std::ptrdiff_t dim(boExpectedDimension(boname,replicated,_currallocatedboatt[GLFeederInfo::ATT_VERTINDEX]));
                        //std::ptrdiff_t dim(boExpectedDimension(boname,replicated,generateindex));
                        ////disableClientState(boname,importattribute);
                        //if ((*it)->_size > 0)
                        //{


                        //    /*WARNING! THIS CODE MUST BE INCAPSULATED INTO A DEALLOCATE FUNCTION IN A PROPER MADE BUFFER OBJECT CLASS
                        //    I DON'T INSERT IT INTO THE GLBufferObjectInfo CLASS CAUSE I CONSIDER IT A TEMPORARY PRIVATE STRUCT*/
                        //    glDeleteBuffers(1,&(*it)->_bohandle);
                        //    (*it)->_bohandle = 0;
                        //    _gpumeminfo.releasedMemory(dim);
                        //}
                        //(*it)->_isvalid = false;
                        //(*it)->_size = 0;
                        ///*********************************************************************************************************/
                        //_currallocatedboatt[boname] = false; 
                        bufferDeAllocationRequested(boname); 
                    }
                    ++ii;
                }
                _borendering = false;
                _lastfeedingusedreplicatedpipeline = false;
                return false;
            }
            else
            {

                //I have to update the invalid buffers requested to be imported
                for(size_t kk = 0;kk < attributestobeupdated.size();++kk)
                    attributestobeupdated[kk] = req[static_cast<ATT_NAMES>(kk)];
                bool failedallocation = false;
                size_t ii = 0;
                typename std::vector<GLBufferObject*>::iterator it = _bo.begin();
                while((it != _bo.end()) && (!failedallocation))
                {
                    ATT_NAMES boname = static_cast<ATT_NAMES>(ii);
                    GLBufferObject* cbo = _bo.at(boname);
                    bool importatt = req[boname];
                    //glBindVertexArray(vaohandlespecificperopenglcontext);

                    /*if a bo is not valid but at this point has a valid handle means that attribute values have been updated but the arity of the vertices/faces didn't change. i can use the already allocated space*/
                    bool notvalidandtoberegenerated = (cbo != NULL) && (!cbo->_isvalid) && (cbo->_bohandle == 0) && (importatt);
                    if (notvalidandtoberegenerated)
                    {
                        cbo->_size = boExpectedSize(boname,replicated,generateindex);
                        std::ptrdiff_t dim = boExpectedDimension(boname,replicated,generateindex);
                        
                        glGenBuffers(1, &cbo->_bohandle);
                        glBindBuffer(cbo->_target, cbo->_bohandle);
                        //we call glGetError BEFORE the glBufferData function in order to clean the error flag
                        GLenum err = glGetError();
                        glBufferData(cbo->_target, dim, NULL, GL_STATIC_DRAW);
                        err = glGetError();
                        //even if there according the MemoryInfo subclass there is enough space we were not able to allocate an attribute buffer object. We have to deallocate all the bos related to this mesh
                        failedallocation = (err == GL_OUT_OF_MEMORY);
                        if (!failedallocation)
                        {
                            //setBufferPointerEnableClientState(boname);
                            setBufferPointer(boname);
                            _gpumeminfo.acquiredMemory(dim);
                        }
                        attributestobeupdated[boname] = !failedallocation;
                        cbo->_isvalid = !failedallocation;
                        _borendering = !failedallocation;
                        glBindBuffer(cbo->_target, 0);
                        _currallocatedboatt[boname] = !failedallocation;
                    }
                    else
                    {
                        attributestobeupdated[boname] = false;
                        //values in the bo have to be updated but the arity of the attribute doesn't change. i can use the already allocated space without reallocating it
                        if ((cbo != NULL) && (!cbo->_isvalid) && (cbo->_bohandle != 0) && (importatt))
                        {
                            attributestobeupdated[boname] = true;
                            cbo->_isvalid = true;
                            _currallocatedboatt[boname] = true;
                        }
                        if ((cbo != NULL) && (!importatt))
                        {
                            cbo->_isvalid = false;
                            _currallocatedboatt[boname] = false;
                        }
                    }

                    //if ((cbo == NULL) || (!cbo->_isvalid))
                    //{
                    //	//glBindVertexArray(_vaohandle);
                    //	disableClientState(boname,requestedattributes);
                    //	//glBindVertexArray(0);
                    //}
                    ++it;
                    ++ii;
                    //glBindVertexArray(0);
                }
                if (failedallocation)
                   buffersDeAllocationRequested();
                _borendering = !failedallocation;
                _lastfeedingusedreplicatedpipeline = replicated;
                return _borendering;
            }
        }

        bool tryToAllocateAndCopyAttributesInBO(const ReqAtts& req)
        {
            std::vector<bool> attributestobeupdated;
            bool immediatemode = !(buffersAllocationFunction(req,attributestobeupdated));
            bool replicated = ReqAtts::isReplicatedPipeline(req);

            if (immediatemode)
                return false;

            bool somethingtoupdate = false;
            for(size_t hh = 0;hh < attributestobeupdated.size();++hh)
                somethingtoupdate = somethingtoupdate || attributestobeupdated[hh];
            if (somethingtoupdate)
            {
                if (replicated)
                {
                    //WARNING!In case we have to update the wedgetexture bo maybe (not always!) we must update also the other buffer already in memory
                    //cause the wedgetexture pipeline force a change in the order of the triangles in GPU.
                    //they are now ordered by the texture seam and not more by the triangle index!
                    if (attributestobeupdated[ATT_WEDGETEXTURE])
                    {
                        for(size_t jj = 0;jj < attributestobeupdated.size();++jj)
                        {
                            ATT_NAMES att = static_cast<ATT_NAMES>(jj);
                            attributestobeupdated[jj] = _currallocatedboatt[att] || attributestobeupdated[jj];
                        }
                    }
                    updateBuffersReplicatedPipeline(attributestobeupdated);
                }
                else
                    updateBuffersIndexedPipeline(attributestobeupdated);
                glFinish();
            }
            return true;
        }

        //bool tryToAllocateAttributesInBO(NORMAL_MODALITY nm,COLOR_MODALITY cm)
        //{
        //	std::vector<bool> importattribute(_bo.size());
        //	std::vector<bool> attributestobeupdated;
        //	attributesToBeImportedInPointBasedPipeline(importattribute, nm, cm);
        //	bool immediatemode = !(buffersAllocationFunction(nm,cm,TX_NONE,importattribute,attributestobeupdated));
        //	if (immediatemode)
        //		return false;
        //
        //	bool somethingtoupdate = false;
        //	for(size_t hh = 0;hh < attributestobeupdated.size();++hh)
        //		somethingtoupdate = somethingtoupdate || attributestobeupdated[hh];
        //	if (somethingtoupdate)
        //		updateBuffersIndexedPipeline(attributestobeupdated);
        //	return true;
        //}

        bool updateBuffersIndexedPipeline(const std::vector<bool>& attributestobeupdated)
        {
            _chunkmap.clear();
            size_t vn = _mesh.VN();
            size_t tn = _mesh.FN();

            size_t facechunk = std::min(size_t(tn),_perbatchprim);
            size_t vertexchunk = std::min(size_t(vn),_perbatchprim);

            std::vector<vcg::Point3f> pv; //position vector
            if (attributestobeupdated[GLFeederInfo::ATT_VERTPOSITION])
                pv.resize(vertexchunk);

            std::vector<vcg::Point3f> nv; //per vertex normal vector
            if (attributestobeupdated[GLFeederInfo::ATT_VERTNORMAL])
                nv.resize(vertexchunk);

            std::vector<vcg::Color4b> cv; // Per vertex color vector
            if (attributestobeupdated[GLFeederInfo::ATT_VERTCOLOR])
                cv.resize(vertexchunk);

            std::vector<float> tv; // per vertex texture coord vector
            if (attributestobeupdated[GLFeederInfo::ATT_VERTTEXTURE])
                tv.resize(vertexchunk * 2);

            size_t chunkingpu = 0;

            for(size_t i=0;i<vn;++i)
            {
                size_t chunkindex = i % vertexchunk;
                if (attributestobeupdated[GLFeederInfo::ATT_VERTPOSITION])
                    pv[chunkindex].Import(_mesh.vert[i].cP());

                if (attributestobeupdated[GLFeederInfo::ATT_VERTNORMAL])
                {
                    nv[chunkindex].Import(_mesh.vert[i].cN());
                    nv[chunkindex].Normalize();
                }

                if (attributestobeupdated[GLFeederInfo::ATT_VERTCOLOR])
                    cv[chunkindex] = _mesh.vert[i].cC();
                if (attributestobeupdated[GLFeederInfo::ATT_VERTTEXTURE])
                {
                    tv[chunkindex*2+0] = _mesh.vert[i].cT().U();
                    tv[chunkindex*2+1] = _mesh.vert[i].cT().V();
                }

                if((i == vn - 1) || (chunkindex == vertexchunk - 1))
                {
                    size_t chunksize = vertexchunk;
                    if (i == vn - 1)
                        chunksize = chunkindex + 1;

                    if (attributestobeupdated[GLFeederInfo::ATT_VERTPOSITION])
                    {
                        GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_VERTPOSITION];
                        glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                        glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&pv[0]);
                        //std::vector<vcg::Point3f> tmppv; //position vector
                        //if (attributestobeupdated[GLFeederInfo::ATT_VERTPOSITION])
                        //    tmppv.resize(vertexchunk);
                        //glGetBufferSubData(GL_ARRAY_BUFFER,0,buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&tmppv[0]);
                        glBindBuffer(GL_ARRAY_BUFFER, 0);
                      
                    }
                    if (attributestobeupdated[GLFeederInfo::ATT_VERTNORMAL])
                    {
                        GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_VERTNORMAL];
                        glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                        glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&nv[0]);
                        glBindBuffer(GL_ARRAY_BUFFER, 0);
                    }
                    if (attributestobeupdated[GLFeederInfo::ATT_VERTCOLOR])
                    {
                        GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_VERTCOLOR];
                        glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                        glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&cv[0]);
                        glBindBuffer(GL_ARRAY_BUFFER, 0);
                    }
                    if (attributestobeupdated[GLFeederInfo::ATT_VERTTEXTURE])
                    {
                        GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_VERTTEXTURE];
                        glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                        glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&tv[0]);
                        glBindBuffer(GL_ARRAY_BUFFER, 0);
                    }
                    glFinish();
                    ++chunkingpu;
                }
            }

            pv.clear();
            nv.clear();
            cv.clear();
            tv.clear();

            chunkingpu = 0;
            std::vector<GLuint> ti(facechunk * 3);
            for(size_t i=0;i<tn;++i)
            {
                size_t chunkindex = i % facechunk;

                ti[chunkindex * 3 + 0] = GLuint(vcg::tri::Index(_mesh,_mesh.face[i].V(0)));
                ti[chunkindex * 3 + 1] = GLuint(vcg::tri::Index(_mesh,_mesh.face[i].V(1)));
                ti[chunkindex * 3 + 2] = GLuint(vcg::tri::Index(_mesh,_mesh.face[i].V(2)));

                if((i == tn - 1) || (chunkindex == facechunk - 1))
                {
                    size_t chunksize = facechunk;
                    if (i == tn - 1)
                        chunksize = chunkindex + 1;

                    if (attributestobeupdated[GLFeederInfo::ATT_VERTINDEX])
                    {
                        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bo[GLFeederInfo::ATT_VERTINDEX]->_bohandle);
                        glBufferSubData(GL_ELEMENT_ARRAY_BUFFER,chunkingpu * facechunk *  _bo[GLFeederInfo::ATT_VERTINDEX]->_components *  _bo[GLFeederInfo::ATT_VERTINDEX]->getSizeOfGLType(),_bo[GLFeederInfo::ATT_VERTINDEX]->_components *  _bo[GLFeederInfo::ATT_VERTINDEX]->getSizeOfGLType() * chunksize,&ti[0]);
                        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
                    }
                    ++chunkingpu;
                }
            }
            return true;
        }

        bool updateBuffersReplicatedPipeline(const std::vector<bool>& attributestobeupdated)
        {
            size_t tn = _mesh.fn;

            size_t facechunk = std::min(size_t(tn),_perbatchprim);

            std::vector<vcg::Point3f> rpv; //position vector
            if (attributestobeupdated[GLFeederInfo::ATT_VERTPOSITION])
                rpv.resize(facechunk * 3);

            std::vector<vcg::Point3f> rnv; //per vertex normal vector
            if (attributestobeupdated[GLFeederInfo::ATT_VERTNORMAL])
                rnv.resize(facechunk * 3);

            std::vector<vcg::Point3f> rfnv; //per face normal vector
            if (attributestobeupdated[GLFeederInfo::ATT_FACENORMAL])
                rfnv.resize(facechunk * 3);

            std::vector<vcg::Color4b> rcv; // Per vertex color vector
            if (attributestobeupdated[GLFeederInfo::ATT_VERTCOLOR])
                rcv.resize(facechunk * 3);

            std::vector<vcg::Color4b> rfcv; // Per vertex color vector
            if (attributestobeupdated[GLFeederInfo::ATT_FACECOLOR])
                rfcv.resize(facechunk * 3);

            std::vector<float> rtv; // per vertex texture coord vector
            if (attributestobeupdated[GLFeederInfo::ATT_VERTTEXTURE])
                rtv.resize(facechunk * 3 * 2);

            std::vector<float> rwtv; // per wedge texture coord vector
            if (attributestobeupdated[GLFeederInfo::ATT_WEDGETEXTURE])
                rwtv.resize(facechunk * 3 * 2);

            size_t chunkingpu = 0;

            //it's a map containing for each texture seams n a vector of all the triangle index ranges having n has texture seam
            //Suppose that in a mesh we have
            //TXS_0{t0,t1,t2,t3}, TXS_4{t4,t5},TXS_0{t6},TXS_-1{t7,t8,t9},TXS_4{t10,t11}
            //so chunkMap will contain
            // -1 -> [<t7,t9>]
            //  0 -> [<t0,t3>,<t6,t6>]
            //  4 -> [<t4,t5>,<t10,t11>]
            //
            //if the map has no-texture coords at all in order to unify the code we fill the ChunkMap with texture seam -1 and a single triangle range going from face_0 to face_n-1


            if (attributestobeupdated[GLFeederInfo::ATT_WEDGETEXTURE] || attributestobeupdated[GLFeederInfo::ATT_VERTTEXTURE])
            {
                _chunkmap.clear();
                if (attributestobeupdated[GLFeederInfo::ATT_WEDGETEXTURE])
                    fillchunkMap(_chunkmap);
                else
                    if(attributestobeupdated[GLFeederInfo::ATT_VERTTEXTURE])
                        _chunkmap[0].push_back(std::make_pair(0,tn-1));
            }

            //default case: no texture is required to be rendered but a non texture attribute has to be updated
            //we have to init the _chunkmap with just one entry (-1...that means no texture) referring all the triangles in the mesh
            if ((!_currallocatedboatt[ATT_VERTTEXTURE] && !_currallocatedboatt[ATT_WEDGETEXTURE]) &&
                (attributestobeupdated[GLFeederInfo::ATT_VERTPOSITION] ||
                attributestobeupdated[GLFeederInfo::ATT_VERTNORMAL] || attributestobeupdated[GLFeederInfo::ATT_FACENORMAL] ||
                attributestobeupdated[GLFeederInfo::ATT_VERTCOLOR] || attributestobeupdated[GLFeederInfo::ATT_FACECOLOR]))
            {
                _chunkmap.clear();
                _chunkmap[-1].push_back(std::make_pair(0,tn-1));
            }

            int t = 0;
            if (attributestobeupdated[GLFeederInfo::ATT_WEDGETEXTURE] || attributestobeupdated[GLFeederInfo::ATT_VERTTEXTURE])
            {
                _texindnumtriangles.clear();
                _texindnumtriangles.resize(_chunkmap.size());
            }
            int i = 0;
            size_t chunkindex = i;
            GLuint triangles = 0;
            for(ChunkMap::const_iterator mit = _chunkmap.begin();mit != _chunkmap.end();++mit)
            {
                for (ChunkVector::const_iterator cit = mit->second.begin();cit != mit->second.end();++cit)
                {
                    for(size_t indf = cit->first;indf<=cit->second;++indf)
                    {
                        chunkindex = i % facechunk;
                        if (attributestobeupdated[GLFeederInfo::ATT_VERTPOSITION])
                        {
                            rpv[chunkindex*3+0].Import(_mesh.face[indf].V(0)->P());
                            rpv[chunkindex*3+1].Import(_mesh.face[indf].V(1)->P());
                            rpv[chunkindex*3+2].Import(_mesh.face[indf].V(2)->P());
                        }
                        if (attributestobeupdated[GLFeederInfo::ATT_VERTNORMAL])
                        {
                            rnv[chunkindex*3+0].Import(_mesh.face[indf].V(0)->N().Normalize());
                            rnv[chunkindex*3+1].Import(_mesh.face[indf].V(1)->N().Normalize());
                            rnv[chunkindex*3+2].Import(_mesh.face[indf].V(2)->N().Normalize());
                        }

                        if (attributestobeupdated[GLFeederInfo::ATT_FACENORMAL])
                        {
                            rfnv[chunkindex*3+0].Import(_mesh.face[indf].N().Normalize());
                            rfnv[chunkindex*3+1].Import(_mesh.face[indf].N().Normalize());
                            rfnv[chunkindex*3+2].Import(_mesh.face[indf].N().Normalize());
                        }

                        if ((attributestobeupdated[GLFeederInfo::ATT_VERTCOLOR]))
                        {
                            rcv[chunkindex*3+0] = _mesh.face[indf].V(0)->C();
                            rcv[chunkindex*3+1] = _mesh.face[indf].V(1)->C();
                            rcv[chunkindex*3+2] = _mesh.face[indf].V(2)->C();
                        }

                        if ((attributestobeupdated[GLFeederInfo::ATT_FACECOLOR]))
                        {
                            rfcv[chunkindex*3+0] = _mesh.face[indf].C();
                            rfcv[chunkindex*3+1] = _mesh.face[indf].C();
                            rfcv[chunkindex*3+2] = _mesh.face[indf].C();
                        }

                        if (attributestobeupdated[GLFeederInfo::ATT_VERTTEXTURE])
                        {
                            rtv[chunkindex*6+0]=float(_mesh.face[indf].V(0)->T().U());
                            rtv[chunkindex*6+1]=float(_mesh.face[indf].V(0)->T().V());
                            rtv[chunkindex*6+2]=float(_mesh.face[indf].V(1)->T().U());
                            rtv[chunkindex*6+3]=float(_mesh.face[indf].V(1)->T().V());
                            rtv[chunkindex*6+4]=float(_mesh.face[indf].V(2)->T().U());
                            rtv[chunkindex*6+5]=float(_mesh.face[indf].V(2)->T().V());
                        }

                        if (attributestobeupdated[GLFeederInfo::ATT_WEDGETEXTURE])
                        {
                            rwtv[chunkindex*6+0]=float(_mesh.face[indf].WT(0).U());
                            rwtv[chunkindex*6+1]=float(_mesh.face[indf].WT(0).V());
                            rwtv[chunkindex*6+2]=float(_mesh.face[indf].WT(1).U());
                            rwtv[chunkindex*6+3]=float(_mesh.face[indf].WT(1).V());
                            rwtv[chunkindex*6+4]=float(_mesh.face[indf].WT(2).U());
                            rwtv[chunkindex*6+5]=float(_mesh.face[indf].WT(2).V());
                        }

                        if((i == tn - 1) || (chunkindex == facechunk - 1))
                        {
                            size_t chunksize = facechunk;
                            if (i == tn - 1)
                                chunksize = chunkindex + 1;

                            if (attributestobeupdated[GLFeederInfo::ATT_VERTPOSITION])
                            {
                                GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_VERTPOSITION];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 * buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rpv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[GLFeederInfo::ATT_VERTNORMAL])
                            {
                                GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_VERTNORMAL];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 * buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rnv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[GLFeederInfo::ATT_FACENORMAL])
                            {
                                GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_FACENORMAL];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 * buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rfnv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[GLFeederInfo::ATT_VERTCOLOR])
                            {
                                GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_VERTCOLOR];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 *buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rcv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[GLFeederInfo::ATT_FACECOLOR])
                            {
                                GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_FACECOLOR];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 *buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rfcv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[GLFeederInfo::ATT_VERTTEXTURE])
                            {
                                GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_VERTTEXTURE];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 *buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rtv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            if (attributestobeupdated[GLFeederInfo::ATT_WEDGETEXTURE])
                            {
                                GLBufferObject* buffobj = _bo[GLFeederInfo::ATT_WEDGETEXTURE];
                                glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
                                glBufferSubData(GL_ARRAY_BUFFER,chunkingpu * facechunk * 3 *buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rwtv[0]);
                                glBindBuffer(GL_ARRAY_BUFFER, 0);
                            }

                            ++chunkingpu;
                        }
                        ++i;
                    }
                    triangles += cit->second - cit->first + 1;
                }

                if (attributestobeupdated[GLFeederInfo::ATT_WEDGETEXTURE] || attributestobeupdated[GLFeederInfo::ATT_VERTTEXTURE])
                    _texindnumtriangles[t] = std::make_pair(mit->first,triangles);
                ++t;
            }

            //return (k != tn)
            //    throw MeshLabException("Mesh has not been properly partitioned");
            return true;
        }

        bool immediateModeRendering(const ReqAtts& req,const std::vector<GLuint>& textureindex = std::vector<GLuint>())
        {
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            if(_mesh.fn==0)
                return false;

            if(req[ATT_MESHCOLOR])
                glColor(_mesh.C());

            //typename MESHTYPE::FaceContainer::iterator fp;
            typename MESHTYPE::FaceIterator fi = _mesh.face.begin();

            short curtexname=-1;
            if(req[ATT_WEDGETEXTURE])
            {
                curtexname=(*fi).WT(0).n();
                if ((curtexname >= 0) && (curtexname < (int)textureindex.size()))
                {
                    glEnable(GL_TEXTURE_2D);
                    glBindTexture(GL_TEXTURE_2D,textureindex[curtexname]);
                }
                else
                {
                    glDisable(GL_TEXTURE_2D);
                }
            }

            if(req[ATT_VERTTEXTURE] && !textureindex.empty()) // in the case of per vertex tex coord we assume that we have a SINGLE texture.
            {
                curtexname = 0;
                glEnable(GL_TEXTURE_2D);
                glBindTexture(GL_TEXTURE_2D,textureindex[curtexname]);
            }

            GLenum primitive = GL_TRIANGLES;
            if (req.primitiveModality() == vcg::GLFeederInfo::PR_POINTS)
                primitive = GL_POINTS;
            glBegin(primitive);

            while(fi!=_mesh.face.end())
            {
                typename MESHTYPE::FaceType & f = *fi;
                if(!f.IsD())
                {
                    if(req[ATT_WEDGETEXTURE])
                        if(f.WT(0).n() != curtexname)
                        {
                            curtexname=(*fi).WT(0).n();
                            glEnd();

                            if (curtexname >= 0)
                            {
                                glEnable(GL_TEXTURE_2D);
                                if(!textureindex.empty())
                                    glBindTexture(GL_TEXTURE_2D,textureindex[curtexname]);
                            }
                            else
                            {
                                glDisable(GL_TEXTURE_2D);
                            }

                            glBegin(GL_TRIANGLES);
                        }

                        if(req[ATT_FACENORMAL])
                            glNormal(f.cN());
                        if(req[ATT_VERTNORMAL])
                            glNormal(f.V(0)->cN());

                        if(req[ATT_FACECOLOR])
                            glColor(f.C());
                        if(req[ATT_VERTCOLOR])
                            glColor(f.V(0)->C());
                        if(req[ATT_VERTTEXTURE])
                            glTexCoord(f.V(0)->T().P());
                        if(req[ATT_WEDGETEXTURE])
                            glTexCoord(f.WT(0).t(0));
                        glVertex(f.V(0)->P());

                        if(req[ATT_VERTNORMAL])
                            glNormal(f.V(1)->cN());
                        if(req[ATT_VERTCOLOR])
                            glColor(f.V(1)->C());
                        if(req[ATT_VERTTEXTURE])
                            glTexCoord(f.V(1)->T().P());
                        if(req[ATT_WEDGETEXTURE])
                            glTexCoord(f.WT(1).t(0));
                        glVertex(f.V(1)->P());

                        if(req[ATT_VERTNORMAL])
                            glNormal(f.V(2)->cN());
                        if(req[ATT_VERTCOLOR])
                            glColor(f.V(2)->C());
                        if(req[ATT_VERTTEXTURE])
                            glTexCoord(f.V(2)->T().P());
                        if(req[ATT_WEDGETEXTURE])
                            glTexCoord(f.WT(2).t(0));
                        glVertex(f.V(2)->P());
                }
                ++fi;
            }

            glEnd();
            glPopAttrib();
            return true;
        }

        void drawTriangles(const ReqAtts& req,const std::vector<GLuint>& textureindex = std::vector<GLuint>())
        {
            //isBORenderingPossible(
            if((!isPossibleToUseBORendering()) || (_mesh.VN() == 0))
                return;
            updateClientState(req);
            bool replicated = ReqAtts::isReplicatedPipeline(_currallocatedboatt);

            if (replicated)
            {
                //qDebug("Replicated drawing");
                int firsttriangleoffset = 0;
                if(!req[ATT_VERTTEXTURE] && !req[ATT_WEDGETEXTURE])
                {
                    glDisable(GL_TEXTURE_2D);
                    glDrawArrays(GL_TRIANGLES,0,_mesh.fn * 3);
                }
                else
                {
                    glEnable(GL_TEXTURE_2D);
                    for(std::vector< std::pair<short,GLuint> >::const_iterator it = _texindnumtriangles.begin();it != _texindnumtriangles.end();++it)
                    {
                        if ((it->first != -1) && (it->first < textureindex.size()))
                            glBindTexture(GL_TEXTURE_2D,textureindex[it->first]);
                        else
                            glBindTexture(GL_TEXTURE_2D,0);
                        glDrawArrays(GL_TRIANGLES,firsttriangleoffset,it->second * 3 - firsttriangleoffset);
                        firsttriangleoffset = it->second * 3;
                    }
                    glBindTexture(GL_TEXTURE_2D,0);
                    glDisable(GL_TEXTURE_2D);
                }

            }
            else
            {
                if(req[ATT_VERTTEXTURE])
                {
                    if (textureindex.size() > 0)
                    {
                        glEnable(GL_TEXTURE_2D);
                        glBindTexture(GL_TEXTURE_2D,textureindex[0]);

                    }
                }
                else
                    glDisable(GL_TEXTURE_2D);


                if  (_bo[ATT_VERTINDEX]->_isvalid)
                {
                    //qDebug("Indexed drawing");
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,_bo[ATT_VERTINDEX]->_bohandle);
                    glDrawElements( GL_TRIANGLES, _mesh.fn * _bo[ATT_VERTINDEX]->_components,GL_UNSIGNED_INT ,NULL);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
                }

                glBindTexture(GL_TEXTURE_2D,0);
                glDisable(GL_TEXTURE_2D);
            }
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
            glBindBuffer(GL_ARRAY_BUFFER,0);
            int ii = 0;
            for(typename std::vector<GLBufferObject*>::const_iterator it = _bo.begin();it != _bo.end();++it)
            {
                ATT_NAMES boname = static_cast<ATT_NAMES>(ii);
                if ((boname != GLFeederInfo::ATT_VERTINDEX) && (boname != GLFeederInfo::ATT_MESHCOLOR))
                    disableClientState(boname,req);
                ++ii;
            }
            /*disable all client state buffers*/
            ReqAtts tmp;
            updateClientState(tmp);
        }

        void drawPoints(const ReqAtts& req)
        {
            if ((!isPossibleToUseBORendering()) || (_mesh.VN() == 0))
                return;
            updateClientState(req);
            glDisable(GL_TEXTURE_2D);
            //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bo[GLFeederInfo::ATT_VERTINDEX]->_bohandle);
            size_t pointsnum = _mesh.vn;
            if (ReqAtts::isReplicatedPipeline(_currallocatedboatt))
                pointsnum = _mesh.fn * 3;

            glDrawArrays(GL_POINTS,0,pointsnum);
            //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

            /*disable all client state buffers*/
            ReqAtts tmp;
            updateClientState(tmp);
        }

        void updateClientState(const ReqAtts& req)
        {
            int ii = 0;
            for(typename std::vector<GLBufferObject*>::const_iterator it = _bo.begin();it != _bo.end();++it)
            {
                ATT_NAMES boname = static_cast<ATT_NAMES>(ii);
                if ((boname != GLFeederInfo::ATT_VERTINDEX) && (boname != GLFeederInfo::ATT_MESHCOLOR))
                {
                    if (req[boname] && _currallocatedboatt[boname] && (*it != NULL))
                    {
                        glBindBuffer((*it)->_target, (*it)->_bohandle);
                        setBufferPointer(boname);
                        glEnableClientState((*it)->_clientstatetag);
                        glBindBuffer((*it)->_target, 0);
                    }
                    else
                        disableClientState(boname,req);
                }
                ++ii;
            }
        }

        //bool isAttributeRequiredToBeDisplayed(ATT_NAMES att)
        //{
        //    bool res = false;
        //    for(std::map<unsigned int,ReqAtts>::const_iterator it = _allreqattsmap.begin();it != _allreqattsmap.end();++it)
        //        res |= it->second[att];
        //    return res;
        //}

        void setBufferPointer( ATT_NAMES boname) const
        {
            if ((boname < GLFeederInfo::ATT_VERTPOSITION) || (boname > GLFeederInfo::ATT_VERTINDEX))
                return;
            GLBufferObject* cbo = _bo[boname];
            if (cbo == NULL)
                return;

            switch(boname)
            {
            case(GLFeederInfo::ATT_VERTPOSITION):
                {
                    glVertexPointer(cbo->_components, cbo->_gltype, 0, 0);
                    break;
                }
            case(GLFeederInfo::ATT_VERTNORMAL):
            case(GLFeederInfo::ATT_FACENORMAL):
                {
                    glNormalPointer(cbo->_gltype, 0, 0);
                    break;
                }
            case(GLFeederInfo::ATT_VERTCOLOR):
            case(GLFeederInfo::ATT_FACECOLOR):
                {
                    glColorPointer(cbo->_components, cbo->_gltype, 0, 0);
                    break;
                }
            case(GLFeederInfo::ATT_VERTTEXTURE):
            case(GLFeederInfo::ATT_WEDGETEXTURE):
                {
                    glTexCoordPointer(cbo->_components, cbo->_gltype, 0, 0);
                    break;
                }
            case(GLFeederInfo::ATT_VERTINDEX):
                {
                    break;
                }
            default : break;
            }
        }

        void disableClientState( ATT_NAMES boname,const ReqAtts& req) const
        {
            if ((boname < GLFeederInfo::ATT_VERTPOSITION) || (boname > GLFeederInfo::ATT_VERTINDEX))
                return;

            switch(boname)
            {
            case(GLFeederInfo::ATT_VERTPOSITION):
                {
                    glDisableClientState(GL_VERTEX_ARRAY);
                    break;
                }
            case(GLFeederInfo::ATT_VERTNORMAL):
            case(GLFeederInfo::ATT_FACENORMAL):
                {
                    if (!req[GLFeederInfo::ATT_VERTNORMAL] && !req[GLFeederInfo::ATT_FACENORMAL])
                        glDisableClientState(GL_NORMAL_ARRAY);
                    break;
                }
            case(GLFeederInfo::ATT_VERTCOLOR):
            case(GLFeederInfo::ATT_FACECOLOR):
                {
                    if (!req[GLFeederInfo::ATT_VERTCOLOR] && !req[GLFeederInfo::ATT_FACECOLOR])
                        glDisableClientState(GL_COLOR_ARRAY);
                    break;
                }
            case(GLFeederInfo::ATT_VERTTEXTURE):
            case(GLFeederInfo::ATT_WEDGETEXTURE):
                {
                    if (!req[GLFeederInfo::ATT_VERTTEXTURE] && !req[GLFeederInfo::ATT_WEDGETEXTURE])
                        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
                    break;
                }
            default:
                {
                    break;
                }

            }
        }

        std::ptrdiff_t bufferObjectsMemoryRequired(const ReqAtts& rqatt) const
        {
            bool replicated = ReqAtts::isReplicatedPipeline(rqatt);
            std::ptrdiff_t result(0);
            bool generateindex = ReqAtts::isVertexIndexingRequired(rqatt);

            for(unsigned int ii = 0;ii < (unsigned int)ATT_NAMES_ARITY;++ii)
            {
                ATT_NAMES nm = static_cast<ATT_NAMES>(ii);
                if (rqatt[nm])
                    result += (std::ptrdiff_t) boExpectedDimension(nm,replicated,generateindex);
            }
            return result;
        }

        //expected number of cells should have the required bo
        //generateindex is true when i have a triangle based mesh
        //				is false when i have a point based mesh
        size_t boExpectedSize(ATT_NAMES name,bool replicatedpipeline,bool generateindex) const
        {
            try
            {
                GLBufferObject& cbo = *_bo.at(name);
                size_t vertnum(_mesh.VN());
                size_t facenum(_mesh.FN());

                switch(name)
                {
                case(GLFeederInfo::ATT_VERTPOSITION):
                case(GLFeederInfo::ATT_VERTNORMAL):
                case(GLFeederInfo::ATT_VERTCOLOR):
                case(GLFeederInfo::ATT_VERTTEXTURE):
                    {
                        if (replicatedpipeline)
                            return facenum * 3 * cbo._components;
                        else
                            return vertnum * cbo._components;
                    }

                case(GLFeederInfo::ATT_FACENORMAL):
                case(GLFeederInfo::ATT_FACECOLOR):
                case(GLFeederInfo::ATT_WEDGETEXTURE):
                    {
                        if (replicatedpipeline)
                            return facenum * 3 * cbo._components;
                        else
                            return 0;
                    }
                case(GLFeederInfo::ATT_VERTINDEX):
                    {
                        if (replicatedpipeline || !generateindex)
                            return 0;
                        else
                            return facenum * cbo._components;
                    }
                default : break;
                }
            }
            catch(std::out_of_range& /*exc*/)
            {
                return 0;
            }
            return 0;
        }

        //generateindex is true when i have a triangle based mesh
        //				is false when i have a point based mesh
        size_t boExpectedDimension(ATT_NAMES name,bool replicatedpipeline,bool generateindex) const
        {
            try
            {
                size_t sz = boExpectedSize(name,replicatedpipeline,generateindex);
                GLBufferObject* cbo = _bo.at(name);
                if (cbo == NULL)
                    return 0;
                else
                    return sz * cbo->getSizeOfGLType();
            }
            catch(std::out_of_range& /*exc*/)
            {
                return 0;
            }
            return 0;
        }

        static bool isPerVertexAttribute(ATT_NAMES name)
        {
            return ((name == GLFeederInfo::ATT_VERTPOSITION) ||(name == GLFeederInfo::ATT_VERTNORMAL) || (name == GLFeederInfo::ATT_VERTCOLOR) || (name == GLFeederInfo::ATT_VERTTEXTURE));
        }

        //it's a map containing for each texture seams n a vector of all the triangle index ranges having n has texture seam
        //Suppose that in a mesh we have
        //TXS_0{t0,t1,t2,t3}, TXS_4{t4,t5},TXS_0{t6},TXS_-1{t7,t8,t9},TXS_4{t10,t11}
        //so chunkMap will contain
        // -1 -> [<t7,t9>]
        //  0 -> [<t0,t3>,<t6,t6>]
        //  4 -> [<t4,t5>,<t10,t11>]

        typedef std::vector< std::pair< GLuint,GLuint > > ChunkVector;
        typedef std::map< short, ChunkVector > ChunkMap;

        void fillchunkMap(ChunkMap& cmap)
        {
            if (!vcg::tri::HasPerWedgeTexCoord(_mesh))
                return;
            cmap.clear();
            typename MESHTYPE::FaceIterator infrange = _mesh.face.begin();
            short texind = std::numeric_limits<short>::max();
            int hh = 0;
            for(typename MESHTYPE::FaceIterator fit = _mesh.face.begin();fit != _mesh.face.end();++fit)
            {
                if (fit->WT(0).N() != texind)
                {
                    if ((texind != std::numeric_limits<short>::max()) || (fit == _mesh.face.end() - 1))
                    {
                        GLuint lowind = std::distance(_mesh.face.begin(),infrange);
                        GLuint topind = std::distance(_mesh.face.begin(),fit) - 1;
                        cmap[texind].push_back(std::make_pair(lowind,topind));
                        infrange = fit;
                    }
                    texind = fit->WT(0).N();
                }
                ++hh;
            }
            cmap[texind].push_back(std::make_pair(std::distance(_mesh.face.begin(),infrange),std::distance(_mesh.face.begin(),_mesh.face.end() - 1)));
        }

        //ideally this should be const. I'm not yet sure if VCGLib will allow me to declare it as constant
        MESHTYPE& _mesh;
        MemoryInfo& _gpumeminfo;

        std::vector<GLBufferObject*> _bo;

        /*_currboatt contains the union of all the requested attributes by each single view on the scene. At the end it represents the BOs allocated in the GPU memory*/
        ReqAtts _currallocatedboatt;
        ///*_allreqattmap contains a map of the requested atts by each single view. it's maintained for the actual rendering step*/
        //std::map<unsigned int,ReqAtts> _allreqattsmap;

        std::vector< std::pair<short,GLuint> > _texindnumtriangles;

        bool _lastfeedingusedreplicatedpipeline;
        bool _borendering;
        size_t _perbatchprim;
        bool _rendermodinitialized;
        ChunkMap  _chunkmap;
    };
}

#endif
