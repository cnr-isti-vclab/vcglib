#pragma once
#include <QMessageBox>
class checkGLError
{
	public:
	static QString makeString(const char* m)
	{
		QString message(m);
		
		switch(glGetError()) {			
			case GL_NO_ERROR: return QString();
				
			case GL_INVALID_ENUM:					message+=(": invalid enum");				break;				
			case GL_INVALID_VALUE:				message+=(": invalid value");				break;
			case GL_INVALID_OPERATION:		message+=(": invalid operation");		break;
			case GL_STACK_OVERFLOW:				message+=(": stack overflow");			break;
			case GL_STACK_UNDERFLOW:			message+=(": stack underflow");			break;				
			case GL_OUT_OF_MEMORY:				message+=(": out of memory");				break;
		}
		return message;
	}
	
	static void qDebug(const char* m) {
		QString message=makeString(m);
    if(message.isEmpty()) return;
		::qDebug(qPrintable(message));
	}	
	
	static void QMessageBox(const char* m, const char* title) {
		QString message=makeString(m);
		QMessageBox::warning(0, title,message);
	}	
	static void QMessageBox(const char* m) {QMessageBox(m,"GL error");};
};
