#ifndef GLW_TEXTURE_H
#define GLW_TEXTURE_H

#include "./renderable.h"

namespace glw
{

class TextureSampleMode
{
	public:

		typedef void              BaseType;
		typedef TextureSampleMode ThisType;

		GLenum minFilter;
		GLenum magFilter;
		GLenum wrapS;
		GLenum wrapT;
		GLenum wrapR;

		TextureSampleMode(void)
		{
			this->clear();
		}

		TextureSampleMode(GLenum rMinFilter, GLenum rMagFilter, GLenum rWrapS, GLenum rWrapT, GLenum rWrapR)
			: minFilter (rMinFilter)
			, magFilter (rMagFilter)
			, wrapS     (rWrapS)
			, wrapT     (rWrapT)
			, wrapR     (rWrapR)
		{
			;
		}

		void clear(void)
		{
			this->minFilter = GLW_DONT_CARE;
			this->magFilter = GLW_DONT_CARE;
			this->wrapS     = GLW_DONT_CARE;
			this->wrapT     = GLW_DONT_CARE;
			this->wrapR     = GLW_DONT_CARE;
		}
};

inline TextureSampleMode texSampleMode(GLenum minFilter = GLW_DONT_CARE, GLenum magFilter = GLW_DONT_CARE, GLenum wrapS = GLW_DONT_CARE, GLenum wrapT = GLW_DONT_CARE, GLenum wrapR = GLW_DONT_CARE)
{
	return TextureSampleMode(minFilter, magFilter, wrapS, wrapT, wrapR);
}

class TextureArguments : public RenderableArguments
{
	public:

		typedef RenderableArguments BaseType;
		typedef TextureArguments    ThisType;

		TextureArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
		}
};

class SafeTexture : public virtual SafeRenderable
{
	public:

		typedef SafeRenderable BaseType;
		typedef SafeTexture    ThisType;

	protected:

		SafeTexture(Context * ctx)
			: SafeObject (ctx)
			, BaseType   (ctx)
		{
			;
		}
};

class Texture : public Renderable, public virtual SafeTexture
{
	public:

		typedef Renderable  BaseType;
		typedef SafeTexture SafeType;
		typedef Texture     ThisType;

	protected:

		Texture(Context * ctx)
			: SafeObject     (ctx)
			, SafeRenderable (ctx)
			, SafeType       (ctx)
			, BaseType       (ctx)
		{
			;
		}

		virtual ~Texture(void)
		{
			this->destroy();
		}

		virtual void doDestroy(Context * ctx, GLuint name)
		{
			(void)ctx;
			if (name == 0) return;
			glDeleteTextures(1, &name);
		}

		virtual void doBind(void)
		{
			glActiveTexture(GL_TEXTURE0 + this->m_unit);
			glBindTexture(this->m_target, this->m_name);
		}

		virtual void doUnbind(void)
		{
			glActiveTexture(GL_TEXTURE0 + this->m_unit);
			glBindTexture(this->m_target, 0);
		}
};

typedef detail::SafeHandle   <Texture> TextureHandle;
typedef detail::UnsafeHandle <Texture> BoundTexture;

};

#endif // GLW_TEXTURE_H
