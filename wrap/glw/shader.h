#ifndef GLW_SHADER_H
#define GLW_SHADER_H

#include <string>

#include "./object.h"

namespace glw
{

class ShaderArguments : public ObjectArguments
{
	public:

		typedef ObjectArguments BaseType;
		typedef ShaderArguments ThisType;

		std::string source;

		ShaderArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
			this->source.clear();
		}
};

class SafeShader : public virtual SafeObject
{
	public:

		typedef SafeObject BaseType;
		typedef SafeShader ThisType;

		const std::string & source(void) const
		{
			return this->m_source;
		}

		const std::string & log(void) const
		{
			return this->m_log;
		}

		bool isCompiled(void) const
		{
			return this->m_compiled;
		}

	protected:

		std::string  m_source;
		std::string  m_log;
		bool         m_compiled;

		SafeShader(Context * ctx)
			: BaseType   (ctx)
			, m_compiled (false)
		{
			;
		}
};

class Shader : public Object, public virtual SafeShader
{
	public:

		typedef Object     BaseType;
		typedef SafeShader SafeType;
		typedef Shader     ThisType;

	protected:

		Shader(Context * ctx)
			: SafeObject (ctx)
			, SafeType   (ctx)
			, BaseType   (ctx)
		{
			;
		}

		virtual ~Shader(void)
		{
			this->destroy();
		}

		virtual GLenum shaderType(void) const = 0;

		bool create(const ShaderArguments & args)
		{
			this->destroy();

			const GLenum shType = this->shaderType();

			this->m_name = glCreateShader(shType);
			this->setBinding(shType, 0);
			this->bind();
			this->compile(args.source);

			return this->m_compiled;
		}

		virtual void doDestroy(Context * ctx, GLuint name)
		{
			(void)ctx;
			if (name == 0) return;
			this->m_source.clear();
			this->m_log.clear();
			this->m_compiled = false;
			glDeleteShader(name);
		}

		virtual void doBind(void)
		{
			;
		}

		virtual void doUnbind(void)
		{
			;
		}

		void compile(const std::string & source)
		{
			GLW_ASSERT(this->isValid());

			const char * src = source.c_str();
			glShaderSource(this->m_name, 1, &src, 0);
			glCompileShader(this->m_name);

			GLint compileStatus = 0;
			glGetShaderiv(this->m_name, GL_COMPILE_STATUS, &compileStatus);

			this->m_source   = source;
			this->m_log      = ThisType::getInfoLog(this->m_name);
			this->m_compiled = (compileStatus != GL_FALSE);
		}

	private:

		static std::string getInfoLog(GLuint Shader)
		{
			std::string log;

			GLint logLen = 0;
			glGetShaderiv(Shader, GL_INFO_LOG_LENGTH, &logLen);
			if (logLen > 0)
			{
				char * sLog = new char[logLen + 1];
				glGetShaderInfoLog(Shader, logLen, &logLen, sLog);
				if (logLen > 0)
				{
					if (sLog[0] != '\0')
					{
						sLog[logLen - 1] = '\0';
						log = sLog;
					}
				}
				delete [] sLog;
			}
			return log;
		}
};

typedef detail::SafeHandle   <Shader> ShaderHandle;
typedef detail::UnsafeHandle <Shader> BoundShader;

} // end namespace glw

#endif // GLW_SHADER_H
