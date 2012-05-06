#ifndef GLW_TEXTURE2D_H
#define GLW_TEXTURE2D_H

#include "./texture.h"

namespace glw
{

class Texture2DArguments : public TextureArguments
{
	public:

		typedef TextureArguments   BaseType;
		typedef Texture2DArguments ThisType;

		GLsizei           width;
		GLsizei           height;
		GLenum            dataFormat;
		GLenum            dataType;
		const void *      data;
		TextureSampleMode sampler;

		Texture2DArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
			this->width      = 0;
			this->height     = 0;
			this->dataFormat = GL_NONE;
			this->dataType   = GL_NONE;
			this->data       = 0;
			this->sampler.clear();
		}
};

class SafeTexture2D : public virtual SafeTexture
{
	public:

		typedef SafeTexture    BaseType;
		typedef SafeTexture2D  ThisType;

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

		SafeTexture2D(Context * ctx)
			: SafeObject     (ctx)
			, SafeRenderable (ctx)
			, BaseType       (ctx)
			, m_width        (0)
			, m_height       (0)
		{
			;
		}
};

class Texture2D : public Texture, public SafeTexture2D
{
	friend class Context;
	friend class detail::SharedObjectBinding<Texture2D>;

	public:

		typedef Texture       BaseType;
		typedef SafeTexture2D SafeType;
		typedef Texture       ThisType;

		virtual Type type(void) const
		{
			return Texture2DType;
		}

		virtual int imageDimensions(void) const
		{
			return 2;
		}

		virtual bool isArray(void) const
		{
			return false;
		}

		void allocateLevel(GLint level, GLsizei width, GLsizei height, GLenum dataFormat, GLenum dataType, const void * data)
		{
			GLW_ASSERT(this->isValid());
			glTexImage2D(this->m_target, level, this->m_format, width, height, 0, dataFormat, dataType, data);
		}

		void setSubImage(GLint level, GLint xoffset, GLint yoffset, GLsizei width, GLsizei height, GLenum dataFormat, GLenum dataType, const void * data)
		{
			GLW_ASSERT(this->isValid());
			glTexSubImage2D(this->m_target, level, xoffset, yoffset, width, height, dataFormat, dataType, data);
		}

		void generateMipmap(void)
		{
			GLW_ASSERT(this->isValid());
			glGenerateMipmap(this->m_target);
		}

		void setSampleMode(const TextureSampleMode & sampler)
		{
			GLW_ASSERT(this->isValid());
			if (GLW_CARE_OF(sampler.minFilter)) glTexParameteri(this->m_target, GL_TEXTURE_MIN_FILTER, sampler.minFilter);
			if (GLW_CARE_OF(sampler.magFilter)) glTexParameteri(this->m_target, GL_TEXTURE_MAG_FILTER, sampler.magFilter);
			if (GLW_CARE_OF(sampler.wrapS    )) glTexParameteri(this->m_target, GL_TEXTURE_WRAP_S,     sampler.wrapS    );
			if (GLW_CARE_OF(sampler.wrapT    )) glTexParameteri(this->m_target, GL_TEXTURE_WRAP_T,     sampler.wrapT    );
		}

	protected:

		Texture2D(Context * ctx)
			: SafeObject     (ctx)
			, SafeRenderable (ctx)
			, SafeTexture    (ctx)
			, BaseType       (ctx)
			, SafeType       (ctx)
		{
			;
		}

		bool create(const Texture2DArguments & args)
		{
			this->destroy();

			GLint activeUnit = 0;
			glGetIntegerv(GL_ACTIVE_TEXTURE, &activeUnit);

			GLint boundName = 0;
			glGetIntegerv(GL_TEXTURE_BINDING_2D, &boundName);

			glGenRenderbuffers(1, &(this->m_name));
			this->setBinding(GL_TEXTURE_2D, 0);
			this->bind();
			this->allocate(args.format, args.width, args.height, args.dataFormat, args.dataType, args.data);
			this->setSampleMode(args.sampler);

			glActiveTexture(activeUnit);
			glBindTexture(GL_TEXTURE_2D, boundName);

			return true;
		}

		void allocate(GLenum format, GLsizei width, GLsizei height, GLenum dataFormat, GLenum dataType, const void * data)
		{
			GLW_ASSERT(this->isValid());
			glTexImage2D(this->m_target, 0, format, width, height, 0, dataFormat, dataType, data);
			this->m_format = format;
			this->m_width  = width;
			this->m_height = height;
		}
};

typedef detail::SafeHandle   <Texture2D> Texture2DHandle;
typedef detail::UnsafeHandle <Texture2D> BoundTexture2D;

} // end namespace glw

#endif // GLW_TEXTURE2D_H
