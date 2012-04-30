#ifndef GLW_OBJECT_H
#define GLW_OBJECT_H

#include "./bookkeeping.h"
#include "./noncopyable.h"
#include "./type.h"
#include "./glheaders.h"

namespace glw
{

class Context;

class ObjectArguments
{
	public:

		typedef void            BaseType;
		typedef ObjectArguments ThisType;

		ObjectArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			;
		}
};

class SafeObject
{
	friend class Object;

	public:

		typedef void       BaseType;
		typedef SafeObject ThisType;

		bool isValid(void) const
		{
			return (this->m_name != 0);
		}

		Context * context(void)
		{
			return this->m_context;
		}

		const Context * context(void) const
		{
			return this->m_context;
		}

		GLuint name(void) const
		{
			return this->m_name;
		}

		virtual Type type(void) const = 0;

	protected:

		Context * m_context;
		GLuint    m_name;

		SafeObject(Context * ctx)
			: m_context (ctx)
			, m_name    (0)
		{
			;
		}

		virtual ~SafeObject(void)
		{
			;
		}
};

class Object : public detail::NonCopyable, public virtual SafeObject
{
	friend class Context;

	public:

		typedef void       BaseType;
		typedef SafeObject SafeType;
		typedef Object     ThisType;

	protected:

		GLenum m_target;
		GLint  m_unit;

		Object(Context * ctx)
			: SafeType (ctx)
			, m_target   (GL_NONE)
			, m_unit     (0)
		{
			;
		}

		void destroy(void)
		{
			if (!this->isValid()) return;

			this->doDestroy(this->m_context, this->m_name);

			this->m_context = 0;
			this->m_name    = 0;
		}

		void setBinding(GLenum target, GLint unit) const
		{
			ThisType * that = const_cast<ThisType *>(this);
			that->m_target = target;
			that->m_unit   = unit;
		}

		void bind(void)
		{
			GLW_ASSERT(this->isValid());
			this->doBind();
		}

		void unbind(void)
		{
			GLW_ASSERT(this->isValid());
			this->doUnbind();
		}

		virtual void doDestroy (Context * ctx, GLuint name) = 0;
		virtual void doBind    (void)                       = 0;
		virtual void doUnbind  (void)                       = 0;
};

typedef detail::SafeHandle   <Object> ObjectHandle;
typedef detail::UnsafeHandle <Object> BoundObject;

};

#endif // GLW_OBJECT_H
