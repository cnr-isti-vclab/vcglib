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
		removeChildNode(scenelst,"instance_visual_scene");
		for(int vsscn = 0;vsscn < scenelst.size();++vsscn)
		{
			QString url = scenelst.at(vsscn).toElement().attribute("url");
		}

		QDomElement vsnode = info->doc->createElement("instance_visual_scene");
		vsnode.setAttribute("url","#vcg-scene-node");
		scenelst.at(0).appendChild(vsnode);

		QDomNodeList vsscene = info->doc->elementsByTagName("library_visual_scenes"); 
		assert(vsscene.size() == 1);
		QDomElement vslnode = info->doc->createElement("visual_scene");
		vslnode.setAttribute("id","vcg-scene-node");
		vslnode.setAttribute("name","vcg-untitled");

		QDomElement vcgnode = info->doc->createElement("node");
		vcgnode.setAttribute("id","vcg-node");
		vcgnode.setAttribute("name","vcg-untitled");

		QDomElement vcginst = info->doc->createElement("instance_geometry");
		vcginst.setAttribute("url","#vcg-mesh-lib");
		vcgnode.appendChild(vcginst);
		vslnode.appendChild(vcgnode);
		vsscene.at(0).appendChild(vslnode);

		QDomNodeList geolib = info->doc->elementsByTagName("library_geometries");
		assert(geolib.size() == 1);
		//removeChildNode(geolib.at(0));
		
		/*QDomElement mshnode;
		mshnode.setTagName("mesh");*/
		QDomElement geonode = info->doc->createElement("geometry");
		geonode.setAttribute("id","vcg-mesh-lib");
		geonode.setAttribute("name","vcg-mesh");
		
		QDomElement meshnode = info->doc->createElement("mesh");
		
		QDomElement srcposnode = info->doc->createElement("source");
		srcposnode.setAttribute("id","vcg-mesh-positions");
		srcposnode.setAttribute("name","vcg-mesh-positions");

		QDomElement arrayposnode = info->doc->createElement("float_array");
		arrayposnode.setAttribute("id","vcg-mesh-positions-array");
		QString countst;
		arrayposnode.setAttribute("count",countst.number(m.vert.size() * 3));
		QString arrp;
		QString arrn;
		for(SaveMeshType::VertexIterator it = m.vert.begin();it != m.vert.end();++it)
		{
			QString x,y,z;
			arrp.append(x.number(it->P().X()) + " " + y.number(it->P().Y()) + " " + z.number(it->P().Z()));
			arrn.append(x.number(it->N().X()) + " " + y.number(it->N().Y()) + " " + z.number(it->N().Z()));
		}
		QDomText ap = info->doc->createTextNode(arrp);

		QDomElement technode = info->doc->createElement("technique_common");
		QDomElement accnode = info->doc->createElement("accessor");
		accnode.setAttribute("source","#vcg-mesh-positions-array");
		accnode.setAttribute("count",countst.number(m.vert.size()));
		accnode.setAttribute("stride","3");
		
		QDomElement parxnode = info->doc->createElement("param");
		parxnode.setAttribute("name","X");
		parxnode.setAttribute("type","float");
		QDomElement parynode = info->doc->createElement("param");
		parynode.setAttribute("name","Y");
		parynode.setAttribute("type","float");
		QDomElement parznode = info->doc->createElement("param");
		parznode.setAttribute("name","Z");
		parznode.setAttribute("type","float");

		accnode.appendChild(parxnode);
		accnode.appendChild(parynode);
		accnode.appendChild(parznode);
		technode.appendChild(accnode);
		arrayposnode.appendChild(ap);
		srcposnode.appendChild(arrayposnode);
		srcposnode.appendChild(technode);
		
		meshnode.appendChild(srcposnode);
		geonode.appendChild(meshnode);
		geolib.at(0).appendChild(geonode);
		
		QDomElement srcnmnode = info->doc->createElement("source");
		srcnmnode.setAttribute("id","vcg-mesh-normals");
		srcnmnode.setAttribute("name","vcg-mesh-normals");

		QDomElement arraynmnode = info->doc->createElement("float_array");
		arraynmnode.setAttribute("id","vcg-mesh-normals-array");
		arraynmnode.setAttribute("count",countst.number(m.vert.size() * 3));

		QDomElement technode2 = info->doc->createElement("technique_common");
		QDomElement accnode2 = info->doc->createElement("accessor");
		accnode2.setAttribute("source","#vcg-mesh-positions-array");
		accnode2.setAttribute("count",countst.number(m.vert.size()));
		accnode2.setAttribute("stride","3");
		
		QDomElement parxnode2 = info->doc->createElement("param");
		parxnode2.setAttribute("name","X");
		parxnode2.setAttribute("type","float");
		QDomElement parynode2 = info->doc->createElement("param");
		parynode2.setAttribute("name","Y");
		parynode2.setAttribute("type","float");
		QDomElement parznode2 = info->doc->createElement("param");
		parznode2.setAttribute("name","Z");
		parznode2.setAttribute("type","float");

		QDomText an = info->doc->createTextNode(arrn);
		

		accnode2.appendChild(parxnode2);
		accnode2.appendChild(parynode2);
		accnode2.appendChild(parznode2);
		technode2.appendChild(accnode2);
		arraynmnode.appendChild(an);
		srcnmnode.appendChild(arraynmnode);
		srcnmnode.appendChild(technode2);
		
		meshnode.appendChild(srcnmnode);
		
		

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