#ifndef __VCGLIB_EXPORTERIDTF
#define __VCGLIB_EXPORTERIDTF


#include <sstream>
#include <fstream>
#include <ostream>


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
		idtf.write(1,"RESOURCE_NAME \"VcgMesh01\"");
		idtf.write(0,"}");


		if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD)
		{

			idtf.write(0,"RESOURCE_LIST \"TEXTURE\" {");
			idtf.write(1,"RESOURCE_COUNT 1");
			idtf.write(1,"RESOURCE 0 {");
			idtf.write(2,"RESOURCE_NAME \"Texture0\"");
			if (m.textures.size() != 0)
				idtf.write(2,"TEXTURE_PATH \"" + m.textures[0] + "\"");
			idtf.write(1,"}");
			idtf.write(0,"}");
		}

		idtf.write(0,"RESOURCE_LIST \"MODEL\" {");
		idtf.write(1,"RESOURCE_COUNT 1");
		idtf.write(1,"RESOURCE 0 {");
		idtf.write(2,"RESOURCE_NAME \"VcgMesh01\"");
		idtf.write(2,"MODEL_TYPE \"MESH\"");
		idtf.write(2,"MESH {");
		idtf.write(3,"FACE_COUNT " + TextUtility::nmbToStr(m.face.size()));
		idtf.write(3,"MODEL_POSITION_COUNT " + TextUtility::nmbToStr(m.vert.size()));
		idtf.write(3,"MODEL_NORMAL_COUNT " + TextUtility::nmbToStr(m.face.size() * 3));
		idtf.write(3,"MODEL_DIFFUSE_COLOR_COUNT 0");
		idtf.write(3,"MODEL_SPECULAR_COLOR_COUNT 0");
		if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD) idtf.write(3,"MODEL_TEXTURE_COORD_COUNT" + TextUtility::nmbToStr(m.face.size() * 3));
		else idtf.write(3,"MODEL_TEXTURE_COORD_COUNT 0");
		idtf.write(3,"MODEL_BONE_COUNT 0");
		idtf.write(3,"MODEL_SHADING_COUNT 1");
		idtf.write(3,"MODEL_SHADING_DESCRIPTION_LIST {");
		idtf.write(4,"SHADING_DESCRIPTION 0 {");
		if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD)
		{
			idtf.write(5,"TEXTURE_LAYER_COUNT 1");
			idtf.write(5,"TEXTURE_COORD_DIMENSION_LIST {");
			idtf.write(6,"TEXTURE_LAYER 0	DIMENSION: 2");
			idtf.write(5,"}");
		}
		else idtf.write(5,"TEXTURE_LAYER_COUNT 0");
		idtf.write(5,"SHADER_ID 0");
		idtf.write(4,"}");
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

		if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD) 
		{
			idtf.write(3,"MESH_FACE_TEXTURE_COORD_LIST {");
			unsigned int nn = 0;
			for(ConstFaceIterator fit = m.face.begin();fit != m.face.end();++fit)  
			{
				idtf.write(4,"FACE " + TextUtility::nmbToStr(nn) + "{");
				idtf.write(5,"TEXTURE_LAYER 0 TEX_COORD: " + TextUtility::nmbToStr(nn) + " " +
					TextUtility::nmbToStr(nn + 1) + " " + 
					TextUtility::nmbToStr(nn + 2));
				nn += 3;
			}
			idtf.write(3,"}");
		}

		idtf.write(3,"MESH_FACE_SHADING_LIST {");
		for(ConstFaceIterator fit = m.face.begin();fit != m.face.end();++fit)  
		{
			idtf.write(4,TextUtility::nmbToStr(0));
		}
		idtf.write(3,"}");

		idtf.write(3,"MODEL_POSITION_LIST {");
		for(ConstVertexIterator vit = m.vert.begin();vit != m.vert.end();++vit)  
		{
			idtf.write(4,TextUtility::nmbToStr(vit->P().X()) + " " +
				TextUtility::nmbToStr(vit->P().Y()) + " " + 
				TextUtility::nmbToStr(vit->P().Z()));
		}
		idtf.write(3,"}");

		idtf.write(3,"MODEL_NORMAL_LIST {");
		for(FaceIterator fitn = m.face.begin();fitn != m.face.end();++fitn)  
		{
			for(unsigned int ii = 0;ii < 3;++ii)
			{
				idtf.write(4,TextUtility::nmbToStr(fitn->N().X()) + " " +
					TextUtility::nmbToStr(fitn->N().Y()) + " " + 
					TextUtility::nmbToStr(fitn->N().Z()));
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
						TextUtility::nmbToStr(fitn->WT(ii).V()) + TextUtility::nmbToStr(0.0f));
				}
			}
			idtf.write(3,"}");
		}

		idtf.write(2,"}");
		idtf.write(1,"}");
		idtf.write(0,"}");
		return E_NOERROR;
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