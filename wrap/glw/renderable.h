#ifndef GLW_RENDERABLE_H
#define GLW_RENDERABLE_H

#include <string>

#include "./object.h"

namespace glw
{

class RenderableArguments : public ObjectArguments
{
	public:

		typedef ObjectArguments     BaseType;
		typedef RenderableArguments ThisType;

		GLenum format;

		RenderableArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
			this->format = GL_NONE;
		}
};

class SafeRenderable : public virtual SafeObject
{
	public:

		typedef SafeObject     BaseType;
		typedef SafeRenderable ThisType;

		GLenum format(void) const
		{
			return this->m_format;
		}

		virtual int  imageDimensions (void) const = 0;
		virtual bool isArray         (void) const = 0;

	protected:

		GLenum m_format;

		SafeRenderable(Context * ctx)
			: BaseType (ctx)
			, m_format (GL_NONE)
		{
			;
		}
};

class Renderable : public Object, public virtual SafeRenderable
{
	public:

		typedef Object         BaseType;
		typedef SafeRenderable SafeType;
		typedef Renderable     ThisType;

	protected:

		Renderable(Context * ctx)
			: SafeObject (ctx)
			, SafeType   (ctx)
			, BaseType   (ctx)
		{
			;
		}
};

typedef detail::SafeHandle   <Renderable> RenderableHandle;
typedef detail::UnsafeHandle <Renderable> BoundRenderable;

};

#endif // GLW_RENDERABLE_H
