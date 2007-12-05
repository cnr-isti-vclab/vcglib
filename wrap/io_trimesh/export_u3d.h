#ifndef __VCGLIB_EXPORTERU3D
#define __VCGLIB_EXPORTERU3D

#include <cstdlib>
#include <string>
#include <QDir>
#include <QString>
#include <QProcess>

#include "export_idtf.h"
#include<vcg/space/point3.h>

namespace vcg {
namespace tri {
namespace io {
namespace u3dparametersclasses
{
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
			CameraParameters()
				:_cam_fov_angle(0.0f),_cam_roll_angle(0.0f),_obj_to_cam_dir(vcg::Point3f(0.0f,0.0f,0.0f)),_obj_to_cam_dist(0.0),_obj_pos(vcg::Point3f(0.0f,0.0f,0.0f))
			{

			}

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
}

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

private:

	static int InvokeConverter(const u3dparametersclasses::IDTFConverterParameters& par)
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

	static void SaveLatex(SaveMeshType& m,const QString& file,const u3dparametersclasses::Movie15Parameters& mov_par)
	{

		Output_File latex(file.toStdString() + ".tex");
		QString u3df = file + ".u3d";
		QStringList file_trim;
		if (u3df.contains("/"))
			file_trim = u3df.split("/");
		else 
			if (u3df.contains("\\"))
				file_trim = u3df.split("/");
		std::string u3d_final =  file_trim.at(file_trim.size() - 1).toStdString();
		latex.write(0,"\\begin{document}");
		latex.write(0,"\\includemovie[");
		latex.write(1,"poster,");
		latex.write(1,"toolbar, %same as `controls\'");
		latex.write(1,"label=" + u3d_final + ",");
		latex.write(1,"text=(" + u3d_final + "),");
		std::string cam_string;
		u3dparametersclasses::Movie15Parameters::CameraParameters* cam = mov_par._campar;
		if (cam != NULL)
		{
			cam_string = cam_string + "3Daac=" + TextUtility::nmbToStr(cam->_cam_fov_angle) + 
				", 3Droll=" + TextUtility::nmbToStr(cam->_cam_roll_angle) +
				", 3Dc2c=" + TextUtility::nmbToStr(cam->_obj_to_cam_dir.X()) + " " + TextUtility::nmbToStr(cam->_obj_to_cam_dir.Y()) + " " + TextUtility::nmbToStr(cam->_obj_to_cam_dir.Z()) +
				", 3Droo=" + TextUtility::nmbToStr(cam->_obj_to_cam_dist) + 
				", 3Dcoo=" + TextUtility::nmbToStr(cam->_obj_pos.X()) + " " + TextUtility::nmbToStr(cam->_obj_pos.Y()) + " " + TextUtility::nmbToStr(cam->_obj_pos.Z()) + ",";
			latex.write(1,cam_string);
		}
		latex.write(1,"3Dlights=File,");
		latex.write(0,"]{\\linewidth}{\\linewidth}{" + u3d_final + "}");
		latex.write(0,"\\end{document}");
	}

public:

	static int Save(SaveMeshType& m,const char* output_file,const char* conv_loc,const u3dparametersclasses::Movie15Parameters& mov_par,const int mask)
	{
		QString curr = QDir::currentPath();
		QString out(output_file);
		QStringList out_trim;
		if (out.contains("/"))
			out_trim = out.split("/");
		else 
			if (out.contains("\\"))
				out_trim = out.split("/");

		QString tmp(QDir::tempPath());
		tmp = tmp + "/" + out_trim.at(out_trim.size() - 1) + ".idtf";

		vcg::tri::io::ExporterIDTF<SaveMeshType>::Save(m,qPrintable(tmp),mask);
		QString conv_loc_st(conv_loc);
		QString output_file_st(output_file);
		u3dparametersclasses::IDTFConverterParameters idtfpar(conv_loc_st,tmp,output_file_st);
		int res = InvokeConverter(idtfpar);
		QDir::setCurrent(curr);
		QString lat (output_file);
		QStringList l = lat.split(".");
		SaveLatex(m,l[0],mov_par);
		QDir dir(QDir::tempPath());
		dir.remove(tmp);
		
		if (res)
			return 0;
		else 
			return 1;
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
