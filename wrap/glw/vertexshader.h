#ifndef GLW_VERTEXSHADER_H
#define GLW_VERTEXSHADER_H

#include "./shader.h"

namespace glw
{

class VertexShaderArguments : public ShaderArguments
{
	public:

		typedef ShaderArguments       BaseType;
		typedef VertexShaderArguments ThisType;

		VertexShaderArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
		}
};

class SafeVertexShader : public virtual SafeShader
{
	public:

		typedef SafeShader       BaseType;
		typedef SafeVertexShader ThisType;

	protected:

		SafeVertexShader(Context * ctx)
			: SafeObject (ctx)
			, BaseType   (ctx)
		{
			;
		}
};

class VertexShader : public Shader, public SafeVertexShader
{
	friend class Context;
	friend class detail::SharedObjectBinding<VertexShader>;

	public:

		typedef Shader           BaseType;
		typedef SafeVertexShader SafeType;
		typedef VertexShader     ThisType;

		virtual Type type(void) const
		{
			return VertexShaderType;
		}

	protected:

		VertexShader(Context * ctx)
			: SafeObject (ctx)
			, SafeShader (ctx)
			, SafeType   (ctx)
			, BaseType   (ctx)
		{
			;
		}

		virtual GLenum shaderType(void) const
		{
			return GL_VERTEX_SHADER;
		}

		bool create(const VertexShaderArguments & args)
		{
			return BaseType::create(args);
		}
};

typedef detail::SafeHandle   <VertexShader> VertexShaderHandle;
typedef detail::UnsafeHandle <VertexShader> BoundVertexShader;

};

#endif // GLW_VERTEXSHADER_H
