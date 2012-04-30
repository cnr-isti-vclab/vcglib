#ifndef GLW_CONTEXT_H
#define GLW_CONTEXT_H

#include "./noncopyable.h"
#include "./buffer.h"
#include "./renderbuffer.h"
#include "./vertexshader.h"
#include "./geometryshader.h"
#include "./fragmentshader.h"
#include "./program.h"
#include "./texture2d.h"
#include "./framebuffer.h"

#include <string>
#include <set>
#include <map>

namespace glw
{

class Context : public detail::NonCopyable
{
	friend class detail::SharedObjectBase<Object>;

	public:

		typedef detail::NonCopyable BaseType;
		typedef Context             ThisType;

		Context(void)
			: m_acquired     (false)
			, m_textureUnits (0)
		{
			;
		}

		virtual ~Context(void)
		{
			this->release();
		}

		bool acquire(void)
		{
			this->release();
			this->initializeTargets();
			this->m_acquired = true;
			return this->m_acquired;
		}

		void release(void)
		{
			if (!this->isAcquired()) return;
			this->m_acquired = false;
			this->terminateTargets();
			this->destroyAllObjects();
		}

		bool isAcquired(void) const
		{
			return this->m_acquired;
		}

		bool isValid(void) const
		{
			return this->isAcquired();
		}

		BufferHandle createBuffer(const BufferArguments & args)
		{
			BufferHandle handle = this->createHandle<Buffer>();
			handle.object()->create(args);
			return handle;
		}

		BoundBuffer bindVertexBuffer(BufferHandle & handle)
		{
			return this->bind(BindingTarget(GL_ARRAY_BUFFER, 0), handle);
		}

		void unbindVertexBuffer(void)
		{
			BufferHandle nullHandle;
			this->bindVertexBuffer(nullHandle);
		}

		BoundBuffer bindIndexBuffer(BufferHandle & handle)
		{
			return this->bind(BindingTarget(GL_ELEMENT_ARRAY_BUFFER, 0), handle);
		}

		void unbindIndexBuffer(void)
		{
			BufferHandle nullHandle;
			this->bindIndexBuffer(nullHandle);
		}

		RenderbufferHandle createRenderbuffer(const RenderbufferArguments & args)
		{
			RenderbufferHandle handle = this->createHandle<Renderbuffer>();
			handle.object()->create(args);
			return handle;
		}

		BoundRenderbuffer bindRenderbuffer(RenderbufferHandle & handle)
		{
			return this->bind(BindingTarget(GL_RENDERBUFFER, 0), handle);
		}

		void unbindRenderbuffer(void)
		{
			RenderbufferHandle nullHandle;
			this->bindRenderbuffer(nullHandle);
		}

		VertexShaderHandle createVertexShader(const VertexShaderArguments & args)
		{
			VertexShaderHandle handle = this->createHandle<VertexShader>();
			handle.object()->create(args);
			return handle;
		}

		BoundVertexShader bindVertexShader(VertexShaderHandle & handle)
		{
			return this->bind(BindingTarget(GL_VERTEX_SHADER, 0), handle);
		}

		void unbindVertexShader(void)
		{
			VertexShaderHandle nullHandle;
			this->bindVertexShader(nullHandle);
		}

		GeometryShaderHandle createGeometryShader(const GeometryShaderArguments & args)
		{
			GeometryShaderHandle handle = this->createHandle<GeometryShader>();
			handle.object()->create(args);
			return handle;
		}

		BoundGeometryShader bindGeometryShader(GeometryShaderHandle & handle)
		{
			return this->bind(BindingTarget(GL_GEOMETRY_SHADER, 0), handle);
		}

		void unbindGeometryShader(void)
		{
			GeometryShaderHandle nullHandle;
			this->bindGeometryShader(nullHandle);
		}

		FragmentShaderHandle createFragmentShader(const FragmentShaderArguments & args)
		{
			FragmentShaderHandle handle = this->createHandle<FragmentShader>();
			handle.object()->create(args);
			return handle;
		}

		BoundFragmentShader bindFragmentShader(FragmentShaderHandle & handle)
		{
			return this->bind(BindingTarget(GL_FRAGMENT_SHADER, 0), handle);
		}

		void unbindFragmentShader(void)
		{
			FragmentShaderHandle nullHandle;
			this->bindFragmentShader(nullHandle);
		}

		ProgramHandle createProgram(const ProgramArguments & args)
		{
			ProgramHandle handle = this->createHandle<Program>();
			handle.object()->create(args);
			return handle;
		}

		BoundProgram bindProgram(ProgramHandle & handle)
		{
			return this->bind(BindingTarget(GL_CURRENT_PROGRAM, 0), handle);
		}

		void unbindProgram(void)
		{
			ProgramHandle nullHandle;
			this->bindProgram(nullHandle);
		}

		Texture2DHandle createTexture2D(const Texture2DArguments & args)
		{
			Texture2DHandle handle = this->createHandle<Texture2D>();
			handle.object()->create(args);
			return handle;
		}

		BoundTexture2D bindTexture2D(GLint unit, Texture2DHandle & handle)
		{
			glActiveTexture(GL_TEXTURE0 + unit);
			return this->bind(BindingTarget(GL_TEXTURE_2D, unit), handle);
		}

		void unbindTexture2D(GLint unit)
		{
			Texture2DHandle nullHandle;
			this->bindTexture2D(unit, nullHandle);
		}

		FramebufferHandle createFramebuffer(const FramebufferArguments & args)
		{
			FramebufferHandle handle = this->createHandle<Framebuffer>();
			handle.object()->create(args);
			return handle;
		}

		BoundFramebuffer bindReadFramebuffer(FramebufferHandle & handle)
		{
			FramebufferHandle nullHandle;
			this->bind(BindingTarget(GL_FRAMEBUFFER, 0), nullHandle);
			return this->bind(BindingTarget(GL_READ_FRAMEBUFFER, 0), handle);
		}

		void unbindReadFramebuffer(void)
		{
			FramebufferHandle nullHandle;
			this->bindReadFramebuffer(nullHandle);
		}

		BoundFramebuffer bindDrawFramebuffer(FramebufferHandle & handle)
		{
			FramebufferHandle nullHandle;
			this->bind(BindingTarget(GL_FRAMEBUFFER, 0), nullHandle);
			return this->bind(BindingTarget(GL_DRAW_FRAMEBUFFER, 0), handle);
		}

		void unbindDrawFramebuffer(void)
		{
			FramebufferHandle nullHandle;
			this->bindDrawFramebuffer(nullHandle);
		}

		BoundFramebuffer bindFramebuffer(FramebufferHandle & handle)
		{
			FramebufferHandle nullHandle;
			this->bind(BindingTarget(GL_READ_FRAMEBUFFER, 0), nullHandle);
			this->bind(BindingTarget(GL_DRAW_FRAMEBUFFER, 0), nullHandle);
			return this->bind(BindingTarget(GL_FRAMEBUFFER, 0), handle);
		}

		void unbindFramebuffer(void)
		{
			FramebufferHandle nullHandle;
			this->bindFramebuffer(nullHandle);
		}

	private:

		typedef detail::SharedObjectBase<Object>          SharedObjectType;

		typedef void *                                    GenericPtr;
		typedef std::pair<GLenum, GLint>                  BindingTarget;

		typedef std::set<SharedObjectType *>              ShaderObjecPtrSet;
		typedef ShaderObjecPtrSet::const_iterator         ShaderObjecPtrSetConstIterator;
		typedef ShaderObjecPtrSet::iterator               ShaderObjecPtrSetIterator;

		typedef std::map<BindingTarget, GenericPtr>       SharedObjectBindingPtrMap;
		typedef SharedObjectBindingPtrMap::const_iterator SharedObjectBindingPtrConstIterator;
		typedef SharedObjectBindingPtrMap::iterator       SharedObjectBindingPtrIterator;
		typedef SharedObjectBindingPtrMap::value_type     SharedObjectBindingPtrValue;

		bool                      m_acquired;
		int                       m_textureUnits;
		ShaderObjecPtrSet         m_shareds;
		SharedObjectBindingPtrMap m_bindings;

		template <typename TObject>
		void initializeTarget(BindingTarget bt)
		{
			typedef TObject                          ObjectType;
			typedef detail::UnsafeHandle<ObjectType> UnsafeObjectType;

			UnsafeObjectType * unsafeObject = 0;
			this->m_bindings.insert(SharedObjectBindingPtrValue(bt, unsafeObject));
		}

		template <typename TObject>
		void terminateTarget(BindingTarget bt)
		{
			typedef TObject                        ObjectType;
			typedef detail::SafeHandle<ObjectType> SafeObjectType;

			SafeObjectType nullHandle;
			this->bind<ObjectType>(bt, nullHandle);
		}

		void initializeTargets(void)
		{
			{
				GLint texUnits = 0;
				glGetIntegerv(GL_MAX_TEXTURE_UNITS, &texUnits);
				this->m_textureUnits = int(texUnits);
			}
			
			this->initializeTarget<Buffer        >(BindingTarget(GL_ARRAY_BUFFER,         0));
			this->initializeTarget<Buffer        >(BindingTarget(GL_ELEMENT_ARRAY_BUFFER, 0));
			this->initializeTarget<Renderbuffer  >(BindingTarget(GL_RENDERBUFFER,         0));
			this->initializeTarget<VertexShader  >(BindingTarget(GL_VERTEX_SHADER,        0));
			this->initializeTarget<GeometryShader>(BindingTarget(GL_GEOMETRY_SHADER,      0));
			this->initializeTarget<FragmentShader>(BindingTarget(GL_FRAGMENT_SHADER,      0));
			this->initializeTarget<Program       >(BindingTarget(GL_CURRENT_PROGRAM,      0));
			for (int i=0; i<this->m_textureUnits; ++i)
			{
				this->initializeTarget<Texture2D>(BindingTarget(GL_TEXTURE_2D, GLint(i)));
			}
			this->initializeTarget<Framebuffer   >(BindingTarget(GL_READ_FRAMEBUFFER,     0));
			this->initializeTarget<Framebuffer   >(BindingTarget(GL_DRAW_FRAMEBUFFER,     0));
			this->initializeTarget<Framebuffer   >(BindingTarget(GL_FRAMEBUFFER,          0));
		}

		void terminateTargets(void)
		{
			this->terminateTarget<Buffer        >(BindingTarget(GL_ARRAY_BUFFER,         0));
			this->terminateTarget<Buffer        >(BindingTarget(GL_ELEMENT_ARRAY_BUFFER, 0));
			this->terminateTarget<Renderbuffer  >(BindingTarget(GL_RENDERBUFFER,         0));
			this->terminateTarget<VertexShader  >(BindingTarget(GL_VERTEX_SHADER,        0));
			this->terminateTarget<GeometryShader>(BindingTarget(GL_GEOMETRY_SHADER,      0));
			this->terminateTarget<FragmentShader>(BindingTarget(GL_FRAGMENT_SHADER,      0));
			this->terminateTarget<Program       >(BindingTarget(GL_CURRENT_PROGRAM,      0));
			for (int i=0; i<this->m_textureUnits; ++i)
			{
				this->terminateTarget<Texture2D>(BindingTarget(GL_TEXTURE_2D, GLint(i)));
			}
			this->terminateTarget<Framebuffer   >(BindingTarget(GL_READ_FRAMEBUFFER,     0));
			this->terminateTarget<Framebuffer   >(BindingTarget(GL_DRAW_FRAMEBUFFER,     0));
			this->terminateTarget<Framebuffer   >(BindingTarget(GL_FRAMEBUFFER,          0));

			this->m_textureUnits = 0;
		}

		template <typename TObject>
		detail::UnsafeHandle<TObject> bind(BindingTarget bt, detail::SafeHandle<TObject> & h)
		{
			typedef TObject                                 ObjectType;
			typedef detail::SharedObjectBinding<ObjectType> SharedObjectBindingType;
			typedef detail::UnsafeHandle<ObjectType>        UnsafeObjectType;

			SharedObjectBindingPtrIterator it = this->m_bindings.find(bt);
			GLW_ASSERT(it != this->m_bindings.end());

			ObjectType * Object = h.object();

			SharedObjectBindingType * currentBinding = static_cast<SharedObjectBindingType *>(it->second);
			if (currentBinding != 0)
			{
				GLW_ASSERT(!currentBinding->isNull());
				// WARNING: as state could have been changed outside GLW, uncommenting the following line may prevent correct binding.
				//if (currentBinding->object() == Object) return UnsafeObjectType(currentBinding);
				if (h.isNull()) currentBinding->object()->unbind();
				currentBinding->setNull();
				currentBinding->unref();
				currentBinding = 0;
				it->second = 0;
			}

			if (h.isNull()) return UnsafeObjectType();

			SharedObjectBindingType * newBinding = new SharedObjectBindingType(h.shared(), bt.first, bt.second);
			newBinding->ref();
			it->second = newBinding;

			Object->setBinding(bt.first, bt.second);
			Object->bind();

			return UnsafeObjectType(newBinding);
		}

		template <typename TObject>
		TObject * createObject(void)
		{
			typedef TObject ObjectType;
			ObjectType * Object = new ObjectType(this);
			return Object;
		}

		void destroyObject(Object * Object)
		{
			GLW_ASSERT(Object != 0);
			Object->destroy();
			delete Object;
		}

		void removeShared(SharedObjectType * shared)
		{
			GLW_ASSERT(shared != 0);
			GLW_ASSERT(this->m_shareds.count(shared) > 0);
			this->m_shareds.erase(shared);
			this->destroyObject(shared->object());
		}

		void destroyAllObjects(void)
		{
			for (ShaderObjecPtrSetIterator it=this->m_shareds.begin(); it!=this->m_shareds.end(); ++it)
			{
				SharedObjectType * shared = *it;
				Object *           Object = shared->object();
				shared->setNull();
				this->destroyObject(Object);
			}
		}

		template <typename TObject>
		detail::SafeHandle<TObject> createHandle(void)
		{
			typedef TObject                          ObjectType;
			typedef detail::SharedObject<ObjectType> SharedObjectType;
			typedef detail::SafeHandle<ObjectType>   SafeType;

			ObjectType *       Object = new ObjectType(this);
			SharedObjectType * shared = new SharedObjectType(this, Object);

			this->m_shareds.insert(shared);

			return SafeType(shared);
		}
};

namespace detail
{

template <typename TObject>
inline void SharedObjectBase<TObject>::signalDestruction(void)
{
	this->m_context->removeShared(this);
}

};

};

#endif // GLW_CONTEXT_H
