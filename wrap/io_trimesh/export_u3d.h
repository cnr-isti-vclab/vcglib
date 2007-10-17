#ifndef __VCGLIB_EXPORTERU3D
#define __VCGLIB_EXPORTERU3D

#include <cstdlib>
#include <string>
#include <QDir.h>
#include<QString.h>

#include "export_idtf.h"

namespace vcg {
namespace tri {
namespace io {

template<typename SaveMeshType>
class ExporterU3D
{
private:

#ifdef WIN32
	static void InvokeConverter(const QString& input_idtf,const QString& output_u3d)
	{
		system(("IDTFConverter.exe -input " + input_idtf + " -output " + output_u3d ).toAscii());
	}

#else
	#ifdef LINUX

		static void InvokeConverter(const QString& input_idtf,const QString& output_u3d)
		{
			system(("IDTFConverter.exe -input " + input_idtf + " -output " + output_u3d ).toAscii()));
		}

	#endif
#endif


public:
	
	static void Save(SaveMeshType& m,const char* outputfile,const int mask)
	{
		QDir dir;
		
		if (!dir.mkdir("__tmp_dir")) return;
		else
		{
			QString curr = QDir::currentPath();
			QString tmp(outputfile);
			tmp = "__tmp_dir/" + tmp + ".idtf";
			vcg::tri::io::ExporterIDTF<SaveMeshType>::Save(m,tmp.toAscii(),mask);
			QDir::setCurrent("../../code/lib/U3D/Bin/Win32/Release");
			InvokeConverter((curr + "/" + tmp).toAscii(),(curr + "/" + outputfile).toAscii());
			QDir::setCurrent(curr);
			dir.remove(tmp.toAscii());
			dir.rmdir("__tmp_dir");
		}

	}
};

}
}
}

#endif