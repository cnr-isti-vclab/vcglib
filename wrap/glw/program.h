#ifndef GLW_PROGRAM_H
#define GLW_PROGRAM_H

#include <memory.h>

#include <string>
#include <vector>
#include <map>

#include "./vertexshader.h"
#include "./geometryshader.h"
#include "./fragmentshader.h"

namespace glw
{

typedef std::vector<ShaderHandle> ShaderHandleVector;

class VertexAttributeBinding
{
	public:

		typedef void                   BaseType;
		typedef VertexAttributeBinding ThisType;

		typedef std::map<std::string, GLuint> Map;
		typedef Map::const_iterator           ConstIterator;
		typedef Map::iterator                 Iterator;
		typedef Map::value_type               Value;

		Map bindings;

		VertexAttributeBinding(void)
		{
			this->clear();
		}

		void clear(void)
		{
			this->bindings.clear();
		}

		GLuint operator [] (const std::string & attributeName) const
		{
			return this->bindings.find(attributeName)->second;
		}

		GLuint & operator [] (const std::string & attributeName)
		{
			return this->bindings[attributeName];
		}
};

class GeometryStage
{
	public:

		typedef void          BaseType;
		typedef GeometryStage ThisType;

		/*
		GLenum inputPrimitiveType;
		GLenum outputPrimitiveType;
		GLint  maxOutputVertices;
		*/

		GeometryStage(void)
		{
			this->clear();
		}

		void clear(void)
		{
			/*
			this->inputPrimitiveType  = GLW_DONT_CARE;
			this->outputPrimitiveType = GLW_DONT_CARE;
			this->maxOutputVertices   = GLW_DONT_CARE;
			*/
		}
};

class TransformFeedbackStream
{
	public:

		typedef void                    BaseType;
		typedef TransformFeedbackStream ThisType;

		typedef std::vector<std::string> VaryingVector;

		VaryingVector varyings;
		GLenum        bufferMode;

		TransformFeedbackStream(void)
		{
			this->clear();
		}

		void clear(void)
		{
			this->varyings.clear();
			this->bufferMode = GL_INTERLEAVED_ATTRIBS;
		}
};

class RasterizerSettings
{
	public:

		typedef void               BaseType;
		typedef RasterizerSettings ThisType;

		enum RasterizerExecution
		{
			DontCare = 0,
			Autodetect,
			ForceEnabled,
			ForceDisabled
		};

		// TODO
		//RasterizerExecution execution;

		RasterizerSettings(void)
		{
			this->clear();
		}

		void clear(void)
		{
			//this->execution = ThisType::Autodetect;
		}
};

class FragmentOutputBinding
{
	public:

		typedef void                  BaseType;
		typedef FragmentOutputBinding ThisType;

		typedef std::map<std::string, GLuint> Map;
		typedef Map::const_iterator           ConstIterator;
		typedef Map::iterator                 Iterator;
		typedef Map::value_type               Value;

		Map bindings;

		FragmentOutputBinding(void)
		{
			this->clear();
		}

		void clear(void)
		{
			this->bindings.clear();
		}

		GLuint operator [] (const std::string & outName) const
		{
			return this->bindings.find(outName)->second;
		}

		GLuint & operator [] (const std::string & outName)
		{
			return this->bindings[outName];
		}
};

class ProgramArguments : public ObjectArguments
{
	public:

		typedef ObjectArguments    BaseType;
		typedef ProgramArguments   ThisType;

		ShaderHandleVector         shaders;
		VertexAttributeBinding     vertexInputs;
		GeometryStage              geometryStage;
		TransformFeedbackStream    feedbackStream;
		RasterizerSettings         rasterSettings;
		FragmentOutputBinding      fragmentOutputs;

		ProgramArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
			this->shaders         .clear();
			this->vertexInputs    .clear();
			this->geometryStage   .clear();
			this->feedbackStream  .clear();
			this->rasterSettings  .clear();
			this->fragmentOutputs .clear();
		}
};

class SafeProgram : public virtual SafeObject
{
	public:

		typedef SafeObject  BaseType;
		typedef SafeProgram ThisType;

		const ProgramArguments & arguments(void) const
		{
			return this->m_arguments;
		}

		const std::string & log(void) const
		{
			return this->m_log;
		}

		bool isLinked(void) const
		{
			return this->m_linked;
		}

	protected:

		ProgramArguments m_arguments;
		std::string      m_log;
		bool             m_linked;

		SafeProgram(Context * ctx)
			: BaseType (ctx)
			, m_linked (false)
		{
			;
		}
};

class Program : public Object, public SafeProgram
{
	friend class Context;
	friend class detail::SharedObjectBinding<Program>;

	public:

		typedef Object      BaseType;
		typedef SafeProgram SafeType;
		typedef Program     ThisType;

		virtual Type type(void) const
		{
			return ProgramType;
		}

		GLint getUniformLocation(const std::string & name) const
		{
			GLW_ASSERT(this->m_uniforms.count(name) > 0);
			return this->m_uniforms.find(name)->second.location;
		}

#define _GLW_IMPLEMENT_SCALAR_UNIFORM_(TYPE, FUNCION_SUFFIX) \
	void setUniform    (const std::string & name, TYPE x                                                       ) { glUniform1         ## FUNCION_SUFFIX     (this->getUniformLocation(name),                                       x         ); } \
	void setUniform    (const std::string & name, TYPE x, TYPE y                                               ) { glUniform2         ## FUNCION_SUFFIX     (this->getUniformLocation(name),                                       x, y      ); } \
	void setUniform    (const std::string & name, TYPE x, TYPE y, TYPE z                                       ) { glUniform3         ## FUNCION_SUFFIX     (this->getUniformLocation(name),                                       x, y, z   ); } \
	void setUniform    (const std::string & name, TYPE x, TYPE y, TYPE z, TYPE w                               ) { glUniform4         ## FUNCION_SUFFIX     (this->getUniformLocation(name),                                       x, y, z, w); }

#define _GLW_IMPLEMENT_VECTOR_UNIFORM_(TYPE, FUNCION_SUFFIX) \
	void setUniform1   (const std::string & name, const TYPE * v,                                 int count = 1) { glUniform1         ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count),                       v         ); } \
	void setUniform2   (const std::string & name, const TYPE * v,                                 int count = 1) { glUniform2         ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count),                       v         ); } \
	void setUniform3   (const std::string & name, const TYPE * v,                                 int count = 1) { glUniform3         ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count),                       v         ); } \
	void setUniform4   (const std::string & name, const TYPE * v,                                 int count = 1) { glUniform4         ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count),                       v         ); }

#define _GLW_IMPLEMENT_MATRIX_UNIFORM_(TYPE, FUNCION_SUFFIX) \
	void setUniform2x2 (const std::string & name, const TYPE * m,                 bool transpose, int count = 1) { glUniformMatrix2   ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count), GLboolean(transpose), m         ); } \
	void setUniform2x3 (const std::string & name, const TYPE * m,                 bool transpose, int count = 1) { glUniformMatrix2x3 ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count), GLboolean(transpose), m         ); } \
	void setUniform2x4 (const std::string & name, const TYPE * m,                 bool transpose, int count = 1) { glUniformMatrix2x4 ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count), GLboolean(transpose), m         ); } \
	void setUniform3x2 (const std::string & name, const TYPE * m,                 bool transpose, int count = 1) { glUniformMatrix3x2 ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count), GLboolean(transpose), m         ); } \
	void setUniform3x3 (const std::string & name, const TYPE * m,                 bool transpose, int count = 1) { glUniformMatrix3   ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count), GLboolean(transpose), m         ); } \
	void setUniform3x4 (const std::string & name, const TYPE * m,                 bool transpose, int count = 1) { glUniformMatrix3x4 ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count), GLboolean(transpose), m         ); } \
	void setUniform4x2 (const std::string & name, const TYPE * m,                 bool transpose, int count = 1) { glUniformMatrix4x2 ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count), GLboolean(transpose), m         ); } \
	void setUniform4x3 (const std::string & name, const TYPE * m,                 bool transpose, int count = 1) { glUniformMatrix4x3 ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count), GLboolean(transpose), m         ); } \
	void setUniform4x4 (const std::string & name, const TYPE * m,                 bool transpose, int count = 1) { glUniformMatrix4   ## FUNCION_SUFFIX ## v (this->getUniformLocation(name), GLsizei(count), GLboolean(transpose), m         ); }

		_GLW_IMPLEMENT_SCALAR_UNIFORM_(int,          i )
		_GLW_IMPLEMENT_SCALAR_UNIFORM_(unsigned int, ui)
		_GLW_IMPLEMENT_SCALAR_UNIFORM_(float,        f )
		_GLW_IMPLEMENT_VECTOR_UNIFORM_(int,          i )
		_GLW_IMPLEMENT_VECTOR_UNIFORM_(unsigned int, ui)
		_GLW_IMPLEMENT_VECTOR_UNIFORM_(float,        f )
		_GLW_IMPLEMENT_MATRIX_UNIFORM_(float,        f )

		GLW_IMPLEMENT_CUSTOM_UNIFORMS;

#undef _GLW_IMPLEMENT_SCALAR_UNIFORM_
#undef _GLW_IMPLEMENT_VECTOR_UNIFORM_
#undef _GLW_IMPLEMENT_MATRIX_UNIFORM_

	protected:

		Program(Context * ctx)
			: SafeObject (ctx)
			, BaseType   (ctx)
			, SafeType   (ctx)
		{
			;
		}

		virtual ~Program(void)
		{
			this->destroy();
		}

		bool create(const ProgramArguments & args)
		{
			this->destroy();

			this->m_arguments = args;

			GLint boundName = 0;
			glGetIntegerv(GL_CURRENT_PROGRAM, &boundName);

			this->m_name = glCreateProgram();

			// shaders
			{
				for (size_t i=0; i<this->m_arguments.shaders.size(); ++i)
				{
					const ShaderHandle & shader = this->m_arguments.shaders[i];
					if (!shader) continue;
					if (!shader->isCompiled()) continue;
					glAttachShader(this->m_name, shader->name());
				}
			}

			// vertex
			{
				for (VertexAttributeBinding::ConstIterator it=this->m_arguments.vertexInputs.bindings.begin(); it!=this->m_arguments.vertexInputs.bindings.end(); ++it)
				{
					glBindAttribLocation(this->m_name, it->second, it->first.c_str());
				}
			}

			// geometry
			{
				;
			}

			// transform feedback
			{
				const size_t count = this->m_arguments.feedbackStream.varyings.size();
				if (count > 0)
				{
					const char ** varyings = new const char * [count];
					for (size_t i=0; i<count; ++i)
					{
						varyings[i]  = this->m_arguments.feedbackStream.varyings[i].c_str();
					}
					glTransformFeedbackVaryings(this->m_name, GLsizei(count), varyings, this->m_arguments.feedbackStream.bufferMode);
					delete [] varyings;
				}
			}

			// TODO
			// rasterizer
			{
				;
			}

			// fragment
			{
				for (FragmentOutputBinding::ConstIterator it=this->m_arguments.fragmentOutputs.bindings.begin(); it!=this->m_arguments.fragmentOutputs.bindings.end(); ++it)
				{
					glBindFragDataLocation(this->m_name, it->second, it->first.c_str());
				}
			}

			glLinkProgram(this->m_name);

			GLint linkStatus = 0;
			glGetProgramiv(this->m_name, GL_LINK_STATUS, &linkStatus);

			this->m_log    = ThisType::getInfoLog(this->m_name);
			this->m_linked = (linkStatus != GL_FALSE);

			if (this->m_linked)
			{
				this->postLink();
			}

			this->setBinding(GL_CURRENT_PROGRAM, 0);
			this->bind();
			// TODO
			// ... nothing to do ...

			glUseProgram(boundName);

			return this->m_linked;
		}

		virtual void doDestroy(Context * ctx, GLuint name)
		{
			(void)ctx;
			if (name == 0) return;
			glDeleteProgram(name);
		}

		virtual void doBind(void)
		{
			glUseProgram(this->m_name);
		}

		virtual void doUnbind(void)
		{
			glUseProgram(0);
		}

	private:

		class UniformInfo
		{
			public:

				typedef void        BaseType;
				typedef UniformInfo ThisType;

				std::string name;
				GLint       location;
				GLenum      type;
				GLint       size;

				UniformInfo(void)
					: location (-1)
					, type     (GL_NONE)
					, size     (0)
				{
					;
				}
		};

		typedef std::map<std::string, UniformInfo> UniformMap;
		typedef UniformMap::const_iterator         UniformMapConstIterator;
		typedef UniformMap::iterator               UniformMapIterator;
		typedef UniformMap::value_type             UniformMapValue;

		UniformMap m_uniforms;

		static std::string getInfoLog(GLuint Program)
		{
			std::string log;
			GLint logLen = 0;
			glGetProgramiv(Program, GL_INFO_LOG_LENGTH, &logLen);
			if (logLen > 0)
			{
				char * sLog = new char[logLen + 1];
				glGetProgramInfoLog(Program, logLen, &logLen, sLog);
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

		void setupUniforms(void)
		{
			this->m_uniforms.clear();

			GLint ucount = 0;
			glGetProgramiv(this->m_name,  GL_ACTIVE_UNIFORMS, &ucount);
			if (ucount <= 0) return;

			GLint ulen = 0;
			glGetProgramiv(this->m_name, GL_ACTIVE_UNIFORM_MAX_LENGTH, &ulen);
			ulen++; // according to specs, +1 (for null) is already accounted, but some implementations are broken.
			if (ulen <= 0) return;

			UniformInfo info;
			GLchar * uname = new GLchar [ulen + 1];
			for (int i=0; i<int(ucount); ++i)
			{
				GLsizei length = 0;
				glGetActiveUniform(this->m_name, GLuint(i), GLsizei(ulen), &length, &(info.size), &(info.type), uname);
				info.name     = uname;
				info.location = glGetUniformLocation(this->m_name, uname);
				this->m_uniforms.insert(UniformMapValue(info.name, info));
			}
			delete [] uname;
		}

		void postLink(void)
		{
			this->setupUniforms();
		}
};

typedef detail::SafeHandle   <Program> ProgramHandle;
typedef detail::UnsafeHandle <Program> BoundProgram;

} // end namespace glw

#endif // GLW_PROGRAM_H
