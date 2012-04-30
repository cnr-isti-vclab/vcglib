#ifndef GLW_FRAGMENTSHADER_H
#define GLW_FRAGMENTSHADER_H

#include "./shader.h"

namespace glw
{

class FragmentShaderArguments : public ShaderArguments
{
	public:

		typedef ShaderArguments         BaseType;
		typedef FragmentShaderArguments ThisType;

		FragmentShaderArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
		}
};

class SafeFragmentShader : public virtual SafeShader
{
	public:

		typedef SafeShader         BaseType;
		typedef SafeFragmentShader ThisType;

	protected:

		SafeFragmentShader(Context * ctx)
			: SafeObject (ctx)
			, BaseType   (ctx)
		{
			;
		}
};

class FragmentShader : public Shader, public SafeFragmentShader
{
	friend class Context;
	friend class detail::SharedObjectBinding<FragmentShader>;

	public:

		typedef Shader             BaseType;
		typedef SafeFragmentShader SafeType;
		typedef FragmentShader     ThisType;

		virtual Type type(void) const
		{
			return FragmentShaderType;
		}

	protected:

		FragmentShader(Context * ctx)
			: SafeObject (ctx)
			, SafeShader (ctx)
			, SafeType   (ctx)
			, BaseType   (ctx)
		{
			;
		}

		virtual GLenum shaderType(void) const
		{
			return GL_FRAGMENT_SHADER;
		}

		bool create(const FragmentShaderArguments & args)
		{
			return BaseType::create(args);
		}
};

typedef detail::SafeHandle   <FragmentShader> FragmentShaderHandle;
typedef detail::UnsafeHandle <FragmentShader> BoundFragmentShader;

};

#endif // GLW_FRAGMENTSHADER_H
