#ifndef GLW_RENDERBUFFER_H
#define GLW_RENDERBUFFER_H

#include "./renderable.h"

namespace glw
{

class RenderbufferArguments : public RenderableArguments
{
	public:

		typedef RenderableArguments   BaseType;
		typedef RenderbufferArguments ThisType;

		GLsizei width;
		GLsizei height;

		RenderbufferArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
			this->width  = 0;
			this->height = 0;
		}
};

class SafeRenderbuffer : public virtual SafeRenderable
{
	public:

		typedef SafeRenderable   BaseType;
		typedef SafeRenderbuffer ThisType;

		GLsizei width(void) const
		{
			return this->m_width;
		}

		GLsizei height(void) const
		{
			return this->m_height;
		}

	protected:

		GLsizei m_width;
		GLsizei m_height;

		SafeRenderbuffer(Context * ctx)
			: SafeObject (ctx)
			, BaseType   (ctx)
			, m_width    (0)
			, m_height   (0)
		{
			;
		}
};

class Renderbuffer : public Renderable, public SafeRenderbuffer
{
	friend class Context;
	friend class detail::SharedObjectBinding<Renderbuffer>;

	public:

		typedef Renderable       BaseType;
		typedef SafeRenderbuffer SafeType;
		typedef Renderbuffer     ThisType;

		virtual Type type(void) const
		{
			return RenderbufferType;
		}

		virtual int imageDimensions(void) const
		{
			return 2;
		}

		virtual bool isArray(void) const
		{
			return false;
		}

	protected:

		Renderbuffer(Context * ctx)
			: SafeObject     (ctx)
			, SafeRenderable (ctx)
			, BaseType       (ctx)
			, SafeType       (ctx)
		{
			;
		}

		virtual ~Renderbuffer(void)
		{
			this->destroy();
		}

		bool create(const RenderbufferArguments & args)
		{
			this->destroy();

			GLint boundName = 0;
			glGetIntegerv(GL_RENDERBUFFER_BINDING, &boundName);

			glGenRenderbuffers(1, &(this->m_name));
			this->setBinding(GL_RENDERBUFFER, 0);
			this->bind();
			this->allocate(args.format, args.width, args.height);

			glBindRenderbuffer(GL_RENDERBUFFER, boundName);

			return true;
		}

		virtual void doDestroy(Context * ctx, GLuint name)
		{
			(void)ctx;
			if (name == 0) return;
			this->m_format = GL_NONE;
			this->m_width  = 0;
			this->m_height = 0;
			glDeleteRenderbuffers(1, &name);
		}

		virtual void doBind(void)
		{
			glBindRenderbuffer(this->m_target, this->m_name);
		}

		virtual void doUnbind(void)
		{
			glBindRenderbuffer(this->m_target, 0);
		}

		void allocate(GLenum format, GLsizei width, GLsizei height)
		{
			GLW_ASSERT(this->isValid());
			glRenderbufferStorage(this->m_target, format, width, height);
			this->m_format = format;
			this->m_width  = width;
			this->m_height = height;
		}
};

typedef detail::SafeHandle   <Renderbuffer> RenderbufferHandle;
typedef detail::UnsafeHandle <Renderbuffer> BoundRenderbuffer;

} // end namespace glw

#endif // GLW_RENDERBUFFER_H
