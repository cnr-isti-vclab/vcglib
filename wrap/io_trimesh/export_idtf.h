#ifndef __VCGLIB_EXPORTERIDTF
#define __VCGLIB_EXPORTERIDTF


#include <sstream>
#include <fstream>
#include <ostream>
#include <vcg/complex/trimesh/update/bounding.h>
#include <wrap/io_trimesh/io_mask.h>


class TextUtility
{
public:
	template<typename NUMERICTYPE>
	static std::string nmbToStr(NUMERICTYPE n)
	{
		std::stringstream ss;
		ss << n;
		return ss.str();
	}
};

class Output_File
{
public:
	Output_File(const std::string& file)
		:_file()
	{
		_file.open(file.c_str(),std::ios::out);
	}

	void write(unsigned int tabl,const std::string& st)
	{
		std::string tmp;
		for(unsigned int ii = 0;ii < tabl;++ii)
			tmp += '\t';
		_file << tmp << st << std::endl;
	}

	~Output_File()
	{
		_file.close();
	}

private:
	std::ofstream _file;
	std::string _tab;
};

#include <QString>
#include <QtGlobal>
#include <fstream>
#include <QtGui\QImage>


namespace vcg {
namespace tri {
namespace io {

template<typename SaveMeshType>
class ExporterIDTF
{
public:
typedef typename SaveMeshType::VertexPointer VertexPointer;
typedef typename SaveMeshType::ScalarType ScalarType;
typedef typename SaveMeshType::VertexType VertexType;
typedef typename SaveMeshType::FaceType FaceType;
typedef typename SaveMeshType::ConstVertexIterator ConstVertexIterator;
typedef typename SaveMeshType::VertexIterator VertexIterator;
typedef typename SaveMeshType::FaceIterator FaceIterator;
typedef typename SaveMeshType::ConstFaceIterator ConstFaceIterator;
typedef typename SaveMeshType::CoordType CoordType;

	enum IDTFError 
	{
		E_NOERROR				// 0
	};

	static const char *ErrorMsg(int error)
	{
		static const char * dae_error_msg[] =
		{
			"No errors"
		};

		if(error>0 || error<0) return "Unknown error";
		else return dae_error_msg[error];
	};

	static int Save(SaveMeshType& m,const char* file,const int mask)
	{
		Output_File idtf(file);
		idtf.write(0,"FILE_FORMAT \"IDTF\"");
		idtf.write(0,"FORMAT_VERSION 100\n");

		idtf.write(0,"NODE \"MODEL\" {");
		idtf.write(1,"NODE_NAME \"VcgMesh01\"");
		idtf.write(1,"PARENT_LIST {");
		idtf.write(2,"PARENT_COUNT 1");
		idtf.write(2,"PARENT 0 {");
		idtf.write(3,"PARENT_NAME \"<NULL>\"");
		idtf.write(3,"PARENT_TM {");
		idtf.write(4,"1.000000 0.000000 0.000000 0.000000");
		idtf.write(4,"0.000000 1.000000 0.000000 0.000000");
		idtf.write(4,"0.000000 0.000000 1.000000 0.000000");
		idtf.write(4,"0.000000 0.000000 0.000000 1.000000");
		idtf.write(3,"}");
		idtf.write(2,"}");
		idtf.write(1,"}");
		idtf.write(1,"RESOURCE_NAME \"MyVcgMesh01\"");
		idtf.write(0,"}");


		if (mask & Mask::IOM_WEDGTEXCOORD)
		{
			idtf.write(0,"");
			idtf.write(0,"RESOURCE_LIST \"SHADER\" {");
			idtf.write(1,"RESOURCE_COUNT " + TextUtility::nmbToStr(m.textures.size()));
			for(unsigned int ii = 0; ii < m.textures.size(); ++ii)
			{
				idtf.write(1,"RESOURCE " + TextUtility::nmbToStr(ii) + " {");
				idtf.write(2,"RESOURCE_NAME \"ModelShader" + TextUtility::nmbToStr(ii) +"\"");
				idtf.write(2,"SHADER_MATERIAL_NAME \"Material1\"");
				idtf.write(2,"SHADER_ACTIVE_TEXTURE_COUNT 1");
				idtf.write(2,"SHADER_TEXTURE_LAYER_LIST {");
				idtf.write(3,"TEXTURE_LAYER 0 {");
				idtf.write(4,"TEXTURE_NAME \"Texture" + TextUtility::nmbToStr(ii) +"\"");
				idtf.write(3,"}");
				idtf.write(2,"}");
			}
			idtf.write(1,"}");
			idtf.write(0,"}");
			idtf.write(0,"");
			idtf.write(0,"RESOURCE_LIST \"MATERIAL\" {");
			idtf.write(1,"RESOURCE_COUNT 1");
			idtf.write(1,"RESOURCE 0 {");
			idtf.write(2,"RESOURCE_NAME \"Material1\"");
			idtf.write(2,"MATERIAL_AMBIENT 0.2 0.2 0.2");
			idtf.write(2,"MATERIAL_DIFFUSE 0.8 0.8 0.8");
			idtf.write(2,"MATERIAL_SPECULAR 0.0 0.0 0.0");
			idtf.write(2,"MATERIAL_EMISSIVE 0.0 0.0 0.0");
			idtf.write(2,"MATERIAL_REFLECTIVITY 0.100000");
			idtf.write(2,"MATERIAL_OPACITY 1.000000");
			idtf.write(1,"}");
			idtf.write(0,"}");
			idtf.write(0,"");
			idtf.write(0,"RESOURCE_LIST \"TEXTURE\" {");
			idtf.write(1,"RESOURCE_COUNT " + TextUtility::nmbToStr(m.textures.size()));
			for(unsigned int ii = 0; ii < m.textures.size();++ii)
			{
				idtf.write(1,"RESOURCE " + TextUtility::nmbToStr(ii) + " {");
				idtf.write(2,"RESOURCE_NAME \"Texture" + TextUtility::nmbToStr(ii) + "\"");
				idtf.write(2,"TEXTURE_PATH \"" + m.textures[ii] + "\"");
				idtf.write(1,"}");
			}
			idtf.write(0,"}");
		}
		idtf.write(0,"");
		idtf.write(0,"RESOURCE_LIST \"MODEL\" {");
		idtf.write(1,"RESOURCE_COUNT 1");
		idtf.write(1,"RESOURCE 0 {");
		idtf.write(2,"RESOURCE_NAME \"MyVcgMesh01\"");
		idtf.write(2,"MODEL_TYPE \"MESH\"");
		idtf.write(2,"MESH {");
		idtf.write(3,"FACE_COUNT " + TextUtility::nmbToStr(m.face.size()));
		idtf.write(3,"MODEL_POSITION_COUNT " + TextUtility::nmbToStr(m.vert.size()));
		idtf.write(3,"MODEL_NORMAL_COUNT " + TextUtility::nmbToStr(m.face.size() * 3));
		idtf.write(3,"MODEL_DIFFUSE_COLOR_COUNT 0");
		idtf.write(3,"MODEL_SPECULAR_COLOR_COUNT 0");
		if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD) idtf.write(3,"MODEL_TEXTURE_COORD_COUNT " + TextUtility::nmbToStr(m.face.size() * 3));
		else idtf.write(3,"MODEL_TEXTURE_COORD_COUNT 0");
		idtf.write(3,"MODEL_BONE_COUNT 0");
		unsigned int mod_sha;
		if (m.textures.size() == 0)
			mod_sha = 1;
		else
			mod_sha = m.textures.size();
		idtf.write(3,"MODEL_SHADING_COUNT " + TextUtility::nmbToStr(mod_sha));
		idtf.write(3,"MODEL_SHADING_DESCRIPTION_LIST {");
		unsigned int hh = 0;
		do
		{
			idtf.write(4,"SHADING_DESCRIPTION " + TextUtility::nmbToStr(hh) + " {");
			if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD)
			{
				idtf.write(5,"TEXTURE_LAYER_COUNT 1");
				idtf.write(5,"TEXTURE_COORD_DIMENSION_LIST {");
				idtf.write(6,"TEXTURE_LAYER 0	DIMENSION: 2");
				idtf.write(5,"}");
				idtf.write(5,"SHADER_ID 0");
			}
			else 
			{
				idtf.write(5,"TEXTURE_LAYER_COUNT 0");
				idtf.write(5,"SHADER_ID 0");
			}
			idtf.write(4,"}");
			++hh;
		}
		while(hh < m.textures.size());
		idtf.write(3,"}");
		idtf.write(3,"MESH_FACE_POSITION_LIST {");
		for(ConstFaceIterator fit = m.face.begin();fit != m.face.end();++fit)  
		{
			idtf.write(4,TextUtility::nmbToStr(fit->V(0) - &(*m.vert.begin())) + " " +
				TextUtility::nmbToStr(fit->V(1) - &(*m.vert.begin())) + " " + 
				TextUtility::nmbToStr(fit->V(2) - &(*m.vert.begin())));
		}
		idtf.write(3,"}");

		idtf.write(3,"MESH_FACE_NORMAL_LIST {");
		unsigned int nn = 0;
		for(ConstFaceIterator fit = m.face.begin();fit != m.face.end();++fit)  
		{
			idtf.write(4,TextUtility::nmbToStr(nn) + " " +
				TextUtility::nmbToStr(nn + 1) + " " + 
				TextUtility::nmbToStr(nn + 2));
			nn += 3;
		}
		idtf.write(3,"}");

		idtf.write(3,"MESH_FACE_SHADING_LIST {");
			for(FaceIterator fit = m.face.begin();fit != m.face.end();++fit)  
			{
				unsigned int texind = 0;
				if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD) 
					texind = fit->WT(0).N();
				idtf.write(4,TextUtility::nmbToStr(texind));
			}
		idtf.write(3,"}");
		if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD) 
		{
				idtf.write(3,"MESH_FACE_TEXTURE_COORD_LIST {"); 
				for(unsigned int ii = 0; ii < m.face.size();++ii)
				{
					idtf.write(4,"FACE " + TextUtility::nmbToStr(ii) + " {");
					idtf.write(5,"TEXTURE_LAYER 0 TEX_COORD: " + TextUtility::nmbToStr(ii * 3) + " " + TextUtility::nmbToStr(ii * 3 + 1) + " " + TextUtility::nmbToStr(ii * 3 + 2));
					idtf.write(4,"}");
				}
				idtf.write(3,"}");
		}

		idtf.write(3,"MODEL_POSITION_LIST {");
		vcg::tri::UpdateBounding<SaveMeshType>::Box(m);
		//ScalarType diag = m.bbox.Diag();
		CoordType center = m.bbox.Center();
		for(ConstVertexIterator vit = m.vert.begin();vit != m.vert.end();++vit)  
		{
			CoordType tmp = vit->P();// - center);// /diag;
			idtf.write(4,TextUtility::nmbToStr(-tmp.X()) + " " +
				TextUtility::nmbToStr(tmp.Z()) + " " + 
				TextUtility::nmbToStr(tmp.Y()));
		}
		idtf.write(3,"}");

		idtf.write(3,"MODEL_NORMAL_LIST {");
		for(FaceIterator fitn = m.face.begin();fitn != m.face.end();++fitn)  
		{
			for(unsigned int ii = 0;ii < 3;++ii)
			{
				fitn->N().Normalize();
				idtf.write(4,TextUtility::nmbToStr(fitn->N().X()) + " " +
					TextUtility::nmbToStr(fitn->N().Z()) + " " + 
					TextUtility::nmbToStr(fitn->N().Y()));
			}
		}
		idtf.write(3,"}");

		if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD)
		{
			idtf.write(3,"MODEL_TEXTURE_COORD_LIST {");
			for(FaceIterator fitn = m.face.begin();fitn != m.face.end();++fitn)  
			{
				for(unsigned int ii = 0;ii < 3;++ii)
				{
					idtf.write(4,TextUtility::nmbToStr(fitn->WT(ii).U()) + " " +
						TextUtility::nmbToStr(-fitn->WT(ii).V()) + " " + TextUtility::nmbToStr(0.0f) + " " + TextUtility::nmbToStr(0.0f));
				}
			}
			idtf.write(3,"}");
		}

		idtf.write(2,"}");
		idtf.write(1,"}");
		idtf.write(0,"}");
		if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD)
		{
			idtf.write(0,"");
			idtf.write(0,"MODIFIER \"SHADING\" {");
			idtf.write(1,"MODIFIER_NAME \"VcgMesh01\"");
			idtf.write(1,"PARAMETERS {");
			idtf.write(2,"SHADER_LIST_COUNT " + TextUtility::nmbToStr(m.textures.size()));
			idtf.write(2,"SHADING_GROUP {");
			for(unsigned int ii = 0; ii < m.textures.size();++ii)
			{
				idtf.write(3,"SHADER_LIST " +  TextUtility::nmbToStr(ii) + "{");
				idtf.write(4,"SHADER_COUNT 1");
				idtf.write(4,"SHADER_NAME_LIST {");
				idtf.write(5,"SHADER 0 NAME: \"ModelShader" + TextUtility::nmbToStr(ii) + "\"");
				idtf.write(4,"}");
				idtf.write(3,"}");
			}
			idtf.write(2,"}");
			idtf.write(1,"}");
			idtf.write(0,"}");
		}
		return E_NOERROR;
	}

	static int GetExportMaskCapability()
	{
		int capability = 0;

		//vert
		capability |= Mask::IOM_VERTNORMAL;


		////wedg
		capability |= Mask::IOM_WEDGTEXCOORD;
		capability |= Mask::IOM_WEDGNORMAL;

		return capability;
	}
};
}
}
}

#endif