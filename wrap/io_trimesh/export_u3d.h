#ifndef __VCGLIB_EXPORTERU3D
#define __VCGLIB_EXPORTERU3D

#include <cstdlib>
#include <string>
#include <QDir>
#include <QString>
#include <QProcess>

#include "export_idtf.h"

namespace vcg {
namespace tri {
namespace io {

template<typename SaveMeshType>
class ExporterU3D
{

public:
	enum U3DError 
	{
		E_NOERROR,				// 0
		E_ABORTED_CONVERSION //1
	};

	static const char *ErrorMsg(int error)
	{
		static const char * dae_error_msg[] =
		{
			"No errors",
			"Conversion Process From Idtf intermediate file to U3D format aborted"
		};

		if(error>1 || error<0) return "Unknown error";
		else return dae_error_msg[error];
	};

	struct Movie15Parameters
	{
		Movie15Parameters()
		{
			_campar = NULL;
		}

		//WARNING: in movie15 y-axis and z-axis have been inverted!!!
		class CameraParameters
		{ 
		public:
			CameraParameters(const float cam_fov_angle,const float cam_roll_angle,
				const vcg::Point3f& obj_to_cam_dir,const float obj_to_cam_dist,
				const vcg::Point3f& obj_pos = vcg::Point3f(0.0f,0.0f,0.0f))
				:_cam_fov_angle(cam_fov_angle),_cam_roll_angle(cam_roll_angle),_obj_to_cam_dir(obj_to_cam_dir),_obj_to_cam_dist(obj_to_cam_dist),_obj_pos(obj_pos)
			{

			}

			float _cam_fov_angle;
			float _cam_roll_angle;
			vcg::Point3f _obj_to_cam_dir;
			float _obj_to_cam_dist;
			vcg::Point3f _obj_pos;

		};
		CameraParameters* _campar;
	};

private:

	struct IDTFConverterParameters
	{
		const QString& _converter_loc;
		const QString& _input_file;
		const QString& _output_file;

		IDTFConverterParameters(const QString& converter_loc,const QString& input_file,const QString& output_file)
			:_converter_loc(converter_loc),_input_file(input_file),_output_file(output_file)
		{
		}
	};

	static int InvokeConverter(const IDTFConverterParameters& par)
	{
		QProcess p;
		QString convstring = par._converter_loc;
		convstring = convstring + " -input " + par._input_file + " -output " + par._output_file; 
		p.start(convstring);
		//wait until the task has been completed
		bool t = p.waitForFinished(-1);
		p.close();
		return (int) t;
	}

	static void SaveLatex(SaveMeshType& m,const QString& file,const Movie15Parameters& mov_par)
	{
		Output_File latex(file.toStdString() + ".tex");
		QString u3df = file + ".u3d";
		latex.write(0,"\\begin{document}");
		latex.write(0,"\\includemovie[");
		latex.write(1,"poster,");
		latex.write(1,"toolbar, %same as `controls\'");
		latex.write(1,"label=" + file.toStdString() + ",");
		latex.write(1,"text=(" + u3df.toStdString() + "),");
		std::string cam_string;
		Movie15Parameters::CameraParameters* cam = mov_par._campar;
		cam_string = cam_string + "3Daac=" + TextUtility::nmbToStr(cam->_cam_fov_angle) + 
			", 3Droll=" + TextUtility::nmbToStr(cam->_cam_roll_angle) +
			", 3Dc2c=" + TextUtility::nmbToStr(cam->_obj_to_cam_dir.X()) + " " + TextUtility::nmbToStr(cam->_obj_to_cam_dir.Y()) + " " + TextUtility::nmbToStr(cam->_obj_to_cam_dir.Z()) +
			", 3Droo=" + TextUtility::nmbToStr(cam->_obj_to_cam_dist) + 
			", 3Dcoo=" + TextUtility::nmbToStr(cam->_obj_pos.X()) + " " + TextUtility::nmbToStr(cam->_obj_pos.Y()) + " " + TextUtility::nmbToStr(cam->_obj_pos.Z()) + ",";
		latex.write(1,cam_string);
		latex.write(1,"3Dlights=File,");
		latex.write(0,"]{\\linewidth}{\\linewidth}{" + u3df.toStdString() + "}");
		latex.write(0,"\\end{document}");
	}

public:

	static int Save(SaveMeshType& m,const char* output_file,const char* conv_loc,const Movie15Parameters& mov_par,const int mask)
	{
		QString curr = QDir::currentPath();
		QString tmp(QDir::tempPath());
		tmp = tmp + "/" + output_file + ".idtf";
		vcg::tri::io::ExporterIDTF<SaveMeshType>::Save(m,qPrintable(tmp),mask);
		IDTFConverterParameters idtfpar(QString(conv_loc),tmp,QString(output_file));
		int res = InvokeConverter(idtfpar);
		QDir::setCurrent(curr);
		QString lat (output_file);
		QStringList l = lat.split(".");
		SaveLatex(m,l[0],mov_par);
		QDir dir(QDir::tempPath());
		dir.remove(tmp);
		return res;
	}

	static int GetExportMaskCapability()
	{
		int capability = 0;

		//vert
		capability |= MeshModel::IOM_VERTNORMAL;


		////wedg
		capability |= MeshModel::IOM_WEDGTEXCOORD;
		capability |= MeshModel::IOM_WEDGNORMAL;

		return capability;
	}

};

}
}
}

#endif
