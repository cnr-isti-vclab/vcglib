#ifndef GLW_BUFFER_H
#define GLW_BUFFER_H

#include "./object.h"

namespace glw
{

class BufferArguments : public ObjectArguments
{
	public:

		typedef ObjectArguments BaseType;
		typedef BufferArguments ThisType;

		GLsizeiptr   size;
		GLenum       usage;
		const void * data;

		BufferArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
			this->size  = 0;
			this->usage = GL_NONE;
			this->data  = 0;
		}
};

class SafeBuffer : public virtual SafeObject
{
	public:

		typedef SafeObject BaseType;
		typedef SafeBuffer ThisType;

		GLsizei size(void) const
		{
			return this->m_size;
		}

		GLenum usage(void) const
		{
			return this->m_usage;
		}

	protected:

		GLsizeiptr m_size;
		GLenum     m_usage;

		SafeBuffer(Context * ctx)
			: BaseType (ctx)
			, m_size   (0)
			, m_usage  (GL_NONE)
		{
			;
		}
};

class Buffer : public Object, public SafeBuffer
{
	friend class Context;
	friend class detail::SharedObjectBinding<Buffer>;

	public:

		typedef Object     BaseType;
		typedef SafeBuffer SafeType;
		typedef Buffer     ThisType;

		virtual Type type(void) const
		{
			return BufferType;
		}

		void setData(const GLsizeiptr size, GLenum usage, const void * data)
		{
			GLW_ASSERT(this->isValid());
			glBufferData(this->m_target, size, data, usage);
			this->m_size  = size;
			this->m_usage = usage;
		}

		void setSubData(GLintptr offset, GLsizeiptr size, const void * data)
		{
			GLW_ASSERT(this->isValid());
			glBufferSubData(this->m_target, offset, size, data);
		}

		void getSubData(GLintptr offset, GLsizeiptr size, void * data)
		{
			GLW_ASSERT(this->isValid());
			glGetBufferSubData(this->m_target, offset, size, data);
		}

		void * map(GLenum access)
		{
			GLW_ASSERT(this->isValid());
			GLW_ASSERT(!this->isMapped());
			void * ptr = glMapBuffer(this->m_target, access);
			if (ptr == 0) return 0;
			this->m_mapAccess = access;
			this->m_mapPtr    = ptr;
			return ptr;
		}

		void unmap(void)
		{
			GLW_ASSERT(this->isValid());
			GLW_ASSERT(this->isMapped());
			glUnmapBuffer(this->m_target);
			this->m_mapAccess = GL_NONE;
			this->m_mapPtr    = 0;
		}

		GLenum mapAccess(void) const
		{
			return this->m_mapAccess;
		}

		bool isMapped(void) const
		{
			return (this->m_mapAccess != GL_NONE);
		}

		void * mapPointer(void) const
		{
			return this->m_mapPtr;
		}

		void vertexAttribPointer(GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const GLvoid * offset)
		{
			GLW_ASSERT(this->isValid());
			GLW_ASSERT(this->m_target == GL_ARRAY_BUFFER);
			glVertexAttribPointer(index, size, type, normalized, stride, offset);
		}

	protected:

		Buffer(Context * ctx)
			: SafeObject  (ctx)
			, SafeType    (ctx)
			, BaseType    (ctx)
			, m_mapAccess (GL_NONE)
			, m_mapPtr    (0)
		{
			;
		}

		virtual ~Buffer(void)
		{
			this->destroy();
		}

		bool create(const BufferArguments & args)
		{
			this->destroy();

			GLint boundName = 0;
			glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &boundName);

			glGenBuffers(1, &(this->m_name));
			this->setBinding(GL_ARRAY_BUFFER, 0);
			this->bind();
			this->setData(args.size, args.usage, args.data);

			glBindBuffer(GL_ARRAY_BUFFER, boundName);

			return true;
		}

		virtual void doDestroy(Context * ctx, GLuint name)
		{
			(void)ctx;
			if (name == 0) return;
			if (this->isMapped()) this->unmap();
			this->m_size      = 0;
			this->m_usage     = GL_NONE;
			this->m_mapAccess = GL_NONE;
			this->m_mapPtr    = 0;
			glDeleteBuffers(1, &name);
		}

		virtual void doBind(void)
		{
			glBindBuffer(this->m_target, this->m_name);
		}

		virtual void doUnbind(void)
		{
			glBindBuffer(this->m_target, 0);
		}

	private:

		GLenum m_mapAccess;
		void * m_mapPtr;
};

typedef detail::SafeHandle   <Buffer> BufferHandle;
typedef detail::UnsafeHandle <Buffer> BoundBuffer;

};

#endif // GLW_BUFFER_H
