#ifndef GLW_GEOMETRYSHADER_H
#define GLW_GEOMETRYSHADER_H

#include "./shader.h"

namespace glw
{

class GeometryShaderArguments : public ShaderArguments
{
	public:

		typedef ShaderArguments         BaseType;
		typedef GeometryShaderArguments ThisType;

		GeometryShaderArguments(void)
		{
			this->clear();
		}

		void clear(void)
		{
			BaseType::clear();
		}
};

class SafeGeometryShader : public virtual SafeShader
{
	public:

		typedef SafeShader         BaseType;
		typedef SafeGeometryShader ThisType;

	protected:

		SafeGeometryShader(Context * ctx)
			: SafeObject (ctx)
			, BaseType   (ctx)
		{
			;
		}
};

class GeometryShader : public Shader, public SafeGeometryShader
{
	friend class Context;
	friend class detail::SharedObjectBinding<GeometryShader>;

	public:

		typedef Shader             BaseType;
		typedef SafeGeometryShader SafeType;
		typedef GeometryShader     ThisType;

		virtual Type type(void) const
		{
			return GeometryShaderType;
		}

	protected:

		GeometryShader(Context * ctx)
			: SafeObject (ctx)
			, SafeShader (ctx)
			, BaseType   (ctx)
			, SafeType   (ctx)
		{
			;
		}

		virtual GLenum shaderType(void) const
		{
			return GL_GEOMETRY_SHADER;
		}

		bool create(const GeometryShaderArguments & args)
		{
			return BaseType::create(args);
		}
};

typedef detail::SafeHandle   <GeometryShader> GeometryShaderHandle;
typedef detail::UnsafeHandle <GeometryShader> BoundGeometryShader;

}  // end namespace glw

#endif // GLW_GEOMETRYSHADER_H
