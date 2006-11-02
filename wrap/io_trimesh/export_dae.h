#ifndef __VCGLIB_EXPORTERDAE
#define __VCGLIB_EXPORTERDAE

#include<wrap/io_trimesh/util_dae.h>

namespace vcg {
namespace tri {
namespace io {

//FCollada Library assumes that SaveMeshType::ScalarType is always a float

template<typename SaveMeshType>
class ExporterDAE : public UtilDAE
{
public:
	static int Save(SaveMeshType &m, const char * filename,AdditionalInfo*& in)
	{
		/*unsigned int ncomp = sizeof(SaveMeshType::CoordType) / sizeof(SaveMeshType::ScalarType);*/
		assert(in != NULL);

		AdditionalInfoDAE* inf = static_cast<AdditionalInfoDAE*>(in);
		InfoDAE* info = inf->dae;
		QDomNodeList scenelst = info->doc->elementsByTagName("scene"); 
		//removeChildNode(scenelst,"instance_visual_scene");
		for(int vsscn = 0;vsscn < scenelst.size();++vsscn)
		{
			QString url = scenelst.at(vsscn).toElement().attribute("url");
		}

		QDomElement vsnode = info->doc->createElement("instance_visual_scene");
		vsnode.setAttribute("url","#vcg-scene-node");
		scenelst.at(0).appendChild(vsnode);

		QDomNodeList geolib = info->doc->elementsByTagName("library_geometries");
		assert(geolib.size() == 1);
		removeChildNode(geolib.at(0));
		
		
		/*QDomElement mshnode;
		mshnode.setTagName("mesh");*/
		QDomElement mshnode = info->doc->createElement("mesh");
		geolib.at(0).appendChild(mshnode);
		QString st = info->doc->toString();


		QFile file(filename);
		if (!file.open(QIODevice::ReadWrite | QIODevice::Truncate))
			return 1;
		info->doc->setContent(&file);
		file.write(st.toAscii());
		file.close();
		return 0;
	}

	static int GetExportMaskCapability()
	{
		int capability = 0;

		//camera
		//capability |= MeshModel::IOM_CAMERA;

		//vert
		//capability |= MeshModel::IOM_VERTTEXCOORD;

		//face
		capability |= MeshModel::IOM_FACEFLAGS;
		//capability |= MeshModel::IOM_FACECOLOR;
		capability |= MeshModel::IOM_FACENORMAL;

		//wedg
		//capability |= MeshModel::IOM_WEDGTEXCOORD;
		capability |= MeshModel::IOM_WEDGNORMAL;

		return capability;
	}
};
}
}
}

#endif