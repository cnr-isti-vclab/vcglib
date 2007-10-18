#ifndef __VCGLIB_EXPORTERU3D
#define __VCGLIB_EXPORTERU3D

#include <cstdlib>
#include <string>
#include <QDir.h>
#include<QString.h>
#include <QProcess.h>

#include "export_idtf.h"

namespace vcg {
namespace tri {
namespace io {

template<typename SaveMeshType>
class ExporterU3D
{
private:

#ifdef WIN32
	static void InvokeConverter(const QString& converter_path,const QString& input_idtf,const QString& output_u3d)
	{
		QProcess p;
		
		p.start(converter_path + "IDTFConverter.exe -input " + input_idtf + " -output " + output_u3d);
		//wait for two minutes
		bool t = p.waitForFinished(120000);
		p.close();
	}

#else
	#ifdef LINUX

		static void InvokeConverter(const char* converter_path,const QString& input_idtf,const QString& output_u3d)
		{
			p.start(converter_path + "IDTFConverter.exe -input " + input_idtf + " -output " + output_u3d);
			//wait for two minutes
			bool t = p.waitForFinished(120000);
			p.close();
		}

	#endif
#endif

	static void SaveLatex(SaveMeshType& m,const QString& file)
	{
		Output_File latex(file.toStdString() + ".tex");
		QString u3df = file + ".u3d";
		latex.write(0,"\\\begin{document}");
		latex.write(0,"\\includemovie[");
		latex.write(1,"poster,");
		latex.write(1,"toolbar, %same as `controls\'");
		latex.write(1,"label=" + file.toStdString() + ",");
		latex.write(1,"text=(" + u3df.toStdString() + "),");
		latex.write(1,"3Daac=60, 3Dcoo=-3.382026195526123 -63.33322525024414 -3.2237343788146973, 3Droo=1.4597717633624103,");
		latex.write(1,"3Dlights=File,");
		latex.write(0,"]{\\linewidth}{\\linewidth}{" + u3df.toStdString() + "}");
		latex.write(0,"\\end{document}");
	}

public:
	
	static void Save(SaveMeshType& m,const char* outputfile,const int mask,const char* converter_path)
	{
		QString curr = QDir::currentPath();
		QString tmp(QDir::tempPath());
		tmp = tmp + "/" + outputfile + ".idtf";
		vcg::tri::io::ExporterIDTF<SaveMeshType>::Save(m,qPrintable(tmp),mask);
		InvokeConverter(converter_path,qPrintable(tmp),qPrintable(curr + "/" + outputfile));
		QDir::setCurrent(curr);
		QString lat (outputfile);
		QStringList l = lat.split(".");
		SaveLatex(m,l[0]);
		QDir dir(QDir::tempPath());
		dir.remove(tmp);
	}
};

}
}
}

#endif