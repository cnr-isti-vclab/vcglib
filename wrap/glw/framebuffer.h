#ifndef GLW_FRAMEBUFFER_H
#define GLW_FRAMEBUFFER_H

#include "./object.h"
#include <vector>

namespace glw
{

class RenderTarget
{
	public:

		typedef void         BaseType;
		typedef RenderTarget ThisType;

		RenderableHandle target;
		GLint            level;
		GLint            layer;
		GLenum           face;

		RenderTarget(void)
		{
			this->clear();
		}

		RenderTarget(RenderableHandle & rTarget, GLint rLevel, GLint rLayer, GLenum rFace)
			: target (rTarget)
			, level  (rLevel)
			, layer  (rLayer)
			, face   (rFace)
		{
			;
		}

		void clear(void)
		{
			this->target.setNull();
			this->level = 0;
			this->layer = -1;
			this->face  = GL_NONE;
		}
};

typedef std::vector<RenderTarget> RenderTargetVector;

inline RenderTarget texture2DTarget(Texture2DHandle & handle, GLint level = 0)
{
	return RenderTarget(handle, level, 0, GL_NONE);
}

inline RenderTarget renderbufferTarget(RenderbufferHandle & handle)
{
	return RenderTarget(handle, 0, 0, GL_NONE);
}

class RenderTargetMapping
{
	public:

		typedef void                BaseType;
		typedef RenderTargetMapping ThisType;

		typedef std::map<GLuint, RenderTarget> Map;
		typedef Map::const_iterator            ConstIterator;
		typedef Map::iterator                  Iterator;
		typedef Map::value_type                Value;

		Map bindings;

		RenderTargetMapping(void)
		{
			this->clear();
		}

		void clear(void)
		{
			this->bindings.clear();
		}

		const RenderTarget & operator [] (GLuint attachmentIndex) const
		{
			return this->bindings.find(attachmentIndex)->second;
		}

		RenderTarget & operator [] (GLuint attachmentIndex)
		{
			return this->bindings[attachmentIndex];
		}
};

class RenderTargetBinding
{
	public:

		typedef void                BaseType;
		typedef RenderTargetBinding ThisType;

		typedef std::map<GLuint, GLuint> Map;
		typedef Map::const_iterator      ConstIterator;
		typedef Map::iterator            Iterator;
		typedef Map::value_type          Value;

		Map bindings;

		RenderTargetBinding(void)
		{
			this->clear();
		}

		void clear(void)
		{
			this->bindings.clear();
		}

		GLuint operator [] (GLuint attachmentIndex) const
		{
			return this->bindings.find(attachmentIndex)->second;
		}

		GLuint & operator [] (GLuint attachmentIndex)
		{
			return this->bindings[attachmentIndex];
		}
};

class FramebufferArguments : public ObjectArguments
{
	public:

		typedef ObjectArguments      BaseType;
		typedef FramebufferArguments ThisType;

		RenderTargetMapping  colorTargets;
		RenderTarget         depthTarget;
		RenderTarget         stencilTarget;
		RenderTarget         depthStencilTarget;
		RenderTargetBinding  targetInputs;

		FramebufferArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
			this->colorTargets       .clear();
			this->depthTarget        .clear();
			this->stencilTarget      .clear();
			this->depthStencilTarget .clear();
			this->targetInputs       .clear();
		}
};

class SafeFramebuffer : public virtual SafeObject
{
	public:

		typedef SafeObject      BaseType;
		typedef SafeFramebuffer ThisType;

	protected:

		SafeFramebuffer(Context * ctx)
			: BaseType(ctx)
		{
			;
		}
};

class Framebuffer : public Object, public SafeFramebuffer
{
	friend class Context;
	friend class detail::SharedObjectBinding<Framebuffer>;

	public:

		typedef Object          BaseType;
		typedef SafeFramebuffer SafeType;
		typedef Framebuffer     ThisType;

		virtual Type type(void) const
		{
			return FramebufferType;
		}

	protected:

		Framebuffer(Context * ctx)
			: SafeObject  (ctx)
			, SafeType    (ctx)
			, BaseType    (ctx)
		{
			;
		}

		virtual ~Framebuffer(void)
		{
			this->destroy();
		}

		bool create(const FramebufferArguments & args)
		{
			this->destroy();

			this->m_config = args;

			GLint boundNameDraw = 0;
			glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &boundNameDraw);

			GLint boundNameRead = 0;
			glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &boundNameRead);

			glGenFramebuffers(1, &(this->m_name));
			this->setBinding(GL_FRAMEBUFFER, 0);
			this->bind();
			this->configure();

			glBindFramebuffer(GL_DRAW_FRAMEBUFFER, boundNameDraw);
			glBindFramebuffer(GL_READ_FRAMEBUFFER, boundNameRead);

			return true;
		}

		virtual void doDestroy(Context * ctx, GLuint name)
		{
			(void)ctx;
			if (name == 0) return;
			glDeleteFramebuffers(1, &name);
			this->m_config.clear();
		}

		virtual void doBind(void)
		{
			glBindFramebuffer(this->m_target, this->m_name);
		}

		virtual void doUnbind(void)
		{
			glBindFramebuffer(this->m_target, 0);
		}

	private:

		FramebufferArguments m_config;

		void configure(void)
		{
			for (RenderTargetMapping::Iterator it=this->m_config.colorTargets.bindings.begin(); it!=this->m_config.colorTargets.bindings.end(); ++it)
			{
				this->attachTarget(GL_COLOR_ATTACHMENT0 + it->first, it->second);
			}

			if (this->m_config.depthStencilTarget.target)
			{
				this->attachTarget(GL_DEPTH_ATTACHMENT,   this->m_config.depthStencilTarget);
				this->attachTarget(GL_STENCIL_ATTACHMENT, this->m_config.depthStencilTarget);
			}
			else
			{
				this->attachTarget(GL_DEPTH_ATTACHMENT,   this->m_config.depthTarget  );
				this->attachTarget(GL_STENCIL_ATTACHMENT, this->m_config.stencilTarget);
			}

			if (this->m_config.colorTargets.bindings.empty())
			{
				glDrawBuffer(GL_NONE);
				glReadBuffer(GL_NONE);
			}
			else
			{
				std::vector<GLenum> drawBuffers;
				drawBuffers.reserve(this->m_config.targetInputs.bindings.size());
				for (RenderTargetBinding::Iterator it=this->m_config.targetInputs.bindings.begin(); it!=this->m_config.targetInputs.bindings.end(); ++it)
				{
					const GLuint fragOutput      = it->second;
					const GLuint attachmentIndex = GL_COLOR_ATTACHMENT0 + it->first;
					if (drawBuffers.size() <= size_t(fragOutput))
					{
						drawBuffers.resize(size_t(fragOutput + 1), GL_NONE);
					}
					drawBuffers[fragOutput] = attachmentIndex;
				}
				glDrawBuffers(GLsizei(drawBuffers.size()), &(drawBuffers[0]));
				glReadBuffer(drawBuffers[0]);
			}
		}

		bool attachTarget(GLenum attachment, RenderTarget & target)
		{
			RenderableHandle & handle = target.target;

			if (!handle)
			{
				glFramebufferRenderbuffer(this->m_target, attachment, GL_RENDERBUFFER, 0);
				return false;
			}

			switch (handle->type())
			{
				case RenderbufferType : glFramebufferRenderbuffer (this->m_target, attachment, GL_RENDERBUFFER, handle->name()              ); break;
				case Texture2DType    : glFramebufferTexture2D    (this->m_target, attachment, GL_TEXTURE_2D,   handle->name(), target.level); break;
				default               : GLW_ASSERT(0);                                                                                         break;
			}

			return true;
		}
};

typedef detail::SafeHandle   <Framebuffer> FramebufferHandle;
typedef detail::UnsafeHandle <Framebuffer> BoundFramebuffer;

};

#endif // GLW_FRAMEBUFFER_H
