#ifndef __GL_OBJECT_H__
#define __GL_OBJECT_H__

#include <GL/gl.h>

class GLObject
{
public:
	GLObject(void)
	{
		this->objectID = 0;
	}

	virtual ~GLObject(void)
	{
	}

	GLuint ObjectID(void) const
	{
		return this->objectID;
	}

	bool ValidObject(void) const
	{
		return (this->objectID != 0);
	}

	virtual void Gen(void) = 0;
	virtual void Del(void) = 0;

protected:
	GLuint objectID;
};

class Bindable
{
public:
	Bindable(void)
	{
		this->bound = false;
	}

	void Bind(void)
	{
		this->bound = true;
		this->DoBind();
	}

	void Unbind(void)
	{
		this->DoUnbind();
		this->bound = false;
	}

	bool IsBound(void) const
	{
		return this->bound;
	}

protected:
	bool bound;

	virtual void DoBind(void) = 0;
	virtual void DoUnbind(void) = 0;
};

#endif __GL_OBJECT_H__
