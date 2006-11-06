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
	static int Save(SaveMeshType &m, const char * filename)
	{
		QDomDocument doc("mydoc");
		QDomElement coll = doc.createElement("COLLADA");
		coll.setAttribute("version","1.4.1");
		coll.setAttribute("xmlns","http://www.collada.org/2005/11/COLLADASchema");
		doc.appendChild(coll);
		QDomElement ass = doc.createElement("asset");
		coll.appendChild(ass);

		QDomElement geolib = doc.createElement("library_geometries");
		
		QDomElement geonode = doc.createElement("geometry");
		geonode.setAttribute("id","vcg-mesh-lib");
		geonode.setAttribute("name","vcg-mesh");
		
		QDomElement meshnode = doc.createElement("mesh");
		
		QDomElement srcposnode = doc.createElement("source");
		srcposnode.setAttribute("id","vcg-mesh-positions");
		srcposnode.setAttribute("name","vcg-mesh-positions");

		QDomElement arrayposnode = doc.createElement("float_array");
		arrayposnode.setAttribute("id","vcg-mesh-positions-array");
		
		QString arrp;
		QString arrn;
		int nvert = 0;
		for(SaveMeshType::VertexIterator it = m.vert.begin();it != m.vert.end();++it)
		{
			if (!(it->IsD()))
			{
				arrp.append(QString::number(it->P().X()) + " " +QString::number(it->P().Y()) + " " + QString::number(it->P().Z()) + " ");
				arrn.append(QString::number(it->N().X()) + " " + QString::number(it->N().Y()) + " " + QString::number(it->N().Z())+ " ");
				++nvert;
			}
		}
		arrayposnode.setAttribute("count",QString::number(nvert * 3));
		QDomText ap = doc.createTextNode(arrp);

		QDomElement technode = doc.createElement("technique_common");
		QDomElement accnode = doc.createElement("accessor");
		accnode.setAttribute("source","#vcg-mesh-positions-array");
		accnode.setAttribute("count",QString::number(nvert));
		accnode.setAttribute("stride","3");
		
		QDomElement parxnode = doc.createElement("param");
		parxnode.setAttribute("name","X");
		parxnode.setAttribute("type","float");
		QDomElement parynode = doc.createElement("param");
		parynode.setAttribute("name","Y");
		parynode.setAttribute("type","float");
		QDomElement parznode = doc.createElement("param");
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
		
		QDomElement srcnmnode = doc.createElement("source");
		srcnmnode.setAttribute("id","vcg-mesh-normals");
		srcnmnode.setAttribute("name","vcg-mesh-normals");

		QDomElement arraynmnode = doc.createElement("float_array");
		arraynmnode.setAttribute("id","vcg-mesh-normals-array");
		arraynmnode.setAttribute("count",QString::number(nvert * 3));

		QDomElement technode2 = doc.createElement("technique_common");
		QDomElement accnode2 = doc.createElement("accessor");
		accnode2.setAttribute("source","#vcg-mesh-normals-array");
		accnode2.setAttribute("count",QString::number(nvert));
		accnode2.setAttribute("stride","3");
		
		QDomElement parxnode2 = doc.createElement("param");
		parxnode2.setAttribute("name","X");
		parxnode2.setAttribute("type","float");
		QDomElement parynode2 = doc.createElement("param");
		parynode2.setAttribute("name","Y");
		parynode2.setAttribute("type","float");
		QDomElement parznode2 = doc.createElement("param");
		parznode2.setAttribute("name","Z");
		parznode2.setAttribute("type","float");

		QDomText an = doc.createTextNode(arrn);
		

		accnode2.appendChild(parxnode2);
		accnode2.appendChild(parynode2);
		accnode2.appendChild(parznode2);
		technode2.appendChild(accnode2);
		arraynmnode.appendChild(an);
		srcnmnode.appendChild(arraynmnode);
		srcnmnode.appendChild(technode2);
		
		meshnode.appendChild(srcnmnode);

		QDomElement vert = doc.createElement("vertices");
		vert.setAttribute("id","vcg-mesh-vertices");
		QDomElement vinp_pos = doc.createElement("input");
		vinp_pos.setAttribute("semantic","POSITION");
		vinp_pos.setAttribute("source","#vcg-mesh-positions");
		QDomElement vinp_nm = doc.createElement("input");
		vinp_nm.setAttribute("semantic","NORMAL");
		vinp_nm.setAttribute("source","#vcg-mesh-normals");

		vert.appendChild(vinp_pos);
		vert.appendChild(vinp_nm);

		meshnode.appendChild(vert);
		
		QDomElement tri = doc.createElement("triangles");
		

		QDomElement tinp_vert = doc.createElement("input");
		tinp_vert.setAttribute("offset","0");
		tinp_vert.setAttribute("semantic","VERTEX");
		tinp_vert.setAttribute("source","#vcg-mesh-vertices");
		QDomElement poly = doc.createElement("p");
		QString triangles_tess;
		int nface = 0;
		for(SaveMeshType::FaceIterator itf = m.face.begin();itf != m.face.end();++itf)
		{
			if (!(itf->IsD()))
			{
				for(unsigned int ii = 0;ii < 3;++ii)
				{
					int ind_v = (*itf).V(ii) - &(m.vert[0]);
					if (triangles_tess == "")
						triangles_tess = QString::number(ind_v);
					else triangles_tess = triangles_tess + " " + QString::number(ind_v);
				}
				++nface;
			}
		}
		tri.setAttribute("count",nface);

		/*QDomElement vslib = doc.createElement("library_visual_scenes");
		doc.appendChild(vslib);*/

		QDomText tri_list = doc.createTextNode(triangles_tess);
		poly.appendChild(tri_list);
		tri.appendChild(tinp_vert);
		tri.appendChild(poly);
		meshnode.appendChild(tri);
		geonode.appendChild(meshnode);
		geolib.appendChild(geonode);
		
		
		coll.appendChild(geolib);	

		QDomElement vslib = doc.createElement("library_visual_scenes");
		QDomElement vslnode = doc.createElement("visual_scene");
		vslnode.setAttribute("id","vcg-scene-node");
		vslnode.setAttribute("name","vcg-untitled");

		QDomElement vcgnode = doc.createElement("node");
		vcgnode.setAttribute("id","vcg-node");
		vcgnode.setAttribute("name","vcg-untitled");

		QDomElement vcginst = doc.createElement("instance_geometry");
		vcginst.setAttribute("url","#vcg-mesh-lib");
		vcgnode.appendChild(vcginst);
		vslnode.appendChild(vcgnode);
		vslib.appendChild(vslnode);

		coll.appendChild(vslib);

		QDomElement vsscn = doc.createElement("scene");
		QDomElement vsnode = doc.createElement("instance_visual_scene");
		vsnode.setAttribute("url","#vcg-scene-node");
		vsscn.appendChild(vsnode);
		coll.appendChild(vsscn);
		
		QString st = doc.toString();
		QFile file(filename);
		if (!file.open(QIODevice::ReadWrite | QIODevice::Truncate))
			return 1;
		doc.setContent(&file);
		file.write(st.toAscii());
		file.close();
		return 0;
		
		
	}

	
	static int Save(SaveMeshType &m, const char * filename,AdditionalInfo*& in)
	{
		/*unsigned int ncomp = sizeof(SaveMeshType::CoordType) / sizeof(SaveMeshType::ScalarType);*/
		assert(in != NULL);

		AdditionalInfoDAE* inf = static_cast<AdditionalInfoDAE*>(in);
		InfoDAE* info = inf->dae;
		QDomNodeList scenelst = info->doc->elementsByTagName("scene"); 
		//removeChildNode(scenelst,"instance_visual_scene");
		assert(scenelst.size() == 1);

		QDomNodeList vsscene = info->doc->elementsByTagName("library_visual_scenes");

		if (info->doc->elementsByTagName("instance_geometry").size() != 0)
		{
			QDomNodeList invisscnlst = info->doc->elementsByTagName("instance_visual_scene");
			assert(invisscnlst.size() == 1);
			QString url;
			referenceToANodeAttribute(invisscnlst.at(0),"url",url); 
			assert(vsscene.size() == 1);
			QDomNode sc = findNodeBySpecificAttributeValue(vsscene.at(0),"visual_scene","id",url);
			//removeChildNode(sc.childNodes(),"instance_geometry");
			for(int no = 0;no < sc.childNodes().size();++no)
			{
				QDomNode oldnode = sc.childNodes().at(no);
				if (sc.childNodes().at(no).toElement().elementsByTagName("instance_geometry").size() != 0)
				{
					sc.removeChild(oldnode);
				}
			}


			QDomElement vcgnode = info->doc->createElement("node");
			vcgnode.setAttribute("id","vcg-node");
			vcgnode.setAttribute("name","vcg-untitled");

			QDomElement vcginst = info->doc->createElement("instance_geometry");
			vcginst.setAttribute("url","#vcg-mesh-lib");
			vcgnode.appendChild(vcginst);
			sc.appendChild(vcgnode);
			vsscene.at(0).appendChild(sc);
			
			QDomNodeList geolib = info->doc->elementsByTagName("library_geometries");
			assert(geolib.size() == 1);
			//removeChildNode(geolib.at(0));
			
			/*QDomElement mshnode;
			mshnode.setTagName("mesh");*/
			
			removeChildNode(geolib.at(0),"geometry","id","vcg-mesh-lib");
			QDomElement geonode = info->doc->createElement("geometry");
			geonode.setAttribute("id","vcg-mesh-lib");
			geonode.setAttribute("name","vcg-mesh");
			
			QDomElement meshnode = info->doc->createElement("mesh");
			
			QDomElement srcposnode = info->doc->createElement("source");
			srcposnode.setAttribute("id","vcg-mesh-positions");
			srcposnode.setAttribute("name","vcg-mesh-positions");

			QDomElement arrayposnode = info->doc->createElement("float_array");
			arrayposnode.setAttribute("id","vcg-mesh-positions-array");
			
			QString arrp;
			arrp.reserve(8 * m.vert.size());
			QString arrn;
			arrn.reserve(8 * m.vert.size());
			int nvert = 0;
			for(SaveMeshType::VertexIterator it = m.vert.begin();it != m.vert.end();++it)
			{
				if (!(it->IsD()))
				{
					arrp.append(QString::number(it->P().X()) + " " +QString::number(it->P().Y()) + " " + QString::number(it->P().Z()) + " ");
					arrn.append(QString::number(it->N().X()) + " " + QString::number(it->N().Y()) + " " + QString::number(it->N().Z())+ " ");
					//arrp.append(QString::number(it->P().X()).append(" ").append(QString::number(it->P().Y())).append(" ").append(QString::number(it->P().Z())).append(" "));
					//arrp.append(QString::number(it->N().X()).append(" ").append(QString::number(it->N().Y())).append(" ").append(QString::number(it->N().Z())).append(" "));
					++nvert;
				}
			}
			arrayposnode.setAttribute("count",QString::number(nvert * 3));
			QDomText ap = info->doc->createTextNode(arrp);

			QDomElement technode = info->doc->createElement("technique_common");
			QDomElement accnode = info->doc->createElement("accessor");
			accnode.setAttribute("source","#vcg-mesh-positions-array");
			accnode.setAttribute("count",QString::number(nvert));
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
			
			QDomElement srcnmnode = info->doc->createElement("source");
			srcnmnode.setAttribute("id","vcg-mesh-normals");
			srcnmnode.setAttribute("name","vcg-mesh-normals");

			QDomElement arraynmnode = info->doc->createElement("float_array");
			arraynmnode.setAttribute("id","vcg-mesh-normals-array");
			arraynmnode.setAttribute("count",QString::number(nvert * 3));

			QDomElement technode2 = info->doc->createElement("technique_common");
			QDomElement accnode2 = info->doc->createElement("accessor");
			accnode2.setAttribute("source","#vcg-mesh-normals-array");
			accnode2.setAttribute("count",QString::number(nvert));
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

			QDomElement vert = info->doc->createElement("vertices");
			vert.setAttribute("id","vcg-mesh-vertices");
			QDomElement vinp_pos = info->doc->createElement("input");
			vinp_pos.setAttribute("semantic","POSITION");
			vinp_pos.setAttribute("source","#vcg-mesh-positions");
			QDomElement vinp_nm = info->doc->createElement("input");
			vinp_nm.setAttribute("semantic","NORMAL");
			vinp_nm.setAttribute("source","#vcg-mesh-normals");

			vert.appendChild(vinp_pos);
			vert.appendChild(vinp_nm);

			meshnode.appendChild(vert);
			
			QDomElement tri = info->doc->createElement("triangles");
			

			QDomElement tinp_vert = info->doc->createElement("input");
			tinp_vert.setAttribute("offset","0");
			tinp_vert.setAttribute("semantic","VERTEX");
			tinp_vert.setAttribute("source","#vcg-mesh-vertices");
			QDomElement poly = info->doc->createElement("p");
			QString triangles_tess;
			int nface = 0;
			triangles_tess.reserve(9*m.face.size());
			for(SaveMeshType::FaceIterator itf = m.face.begin();itf != m.face.end();++itf)
			{
				if (!(itf->IsD()))
				{
					for(unsigned int ii = 0;ii < 3;++ii)
					{
						int ind_v = (*itf).V(ii) - &(m.vert[0]);
						if (triangles_tess == "")
							triangles_tess = QString::number(ind_v);
						else triangles_tess = triangles_tess.append(" ").append(QString::number(ind_v));
					}
					++nface;
				}
			}
			tri.setAttribute("count",nface);

			QDomText tri_list = info->doc->createTextNode(triangles_tess);
			poly.appendChild(tri_list);
			tri.appendChild(tinp_vert);
			tri.appendChild(poly);
			meshnode.appendChild(tri);
			geonode.appendChild(meshnode);
			geolib.at(0).appendChild(geonode);
		}
		else
		{
			removeChildNode(scenelst,"instance_visual_scene");
			for(int vsscn = 0;vsscn < scenelst.size();++vsscn)
			{
				QString url = scenelst.at(vsscn).toElement().attribute("url");
			}

			QDomElement vsnode = info->doc->createElement("instance_visual_scene");
			vsnode.setAttribute("url","#vcg-scene-node");
			scenelst.at(0).appendChild(vsnode);

			
			int vsscene_size = vsscene.size();
			assert(vsscene.size() == 1);
			removeChildNode(vsscene,"visual_scene","id","vcg-scene-node");
			QDomElement vslnode = info->doc->createElement("visual_scene");
			vslnode.setAttribute("id","vcg-scene-node");
			vslnode.setAttribute("name","vcg-untitled");

			QDomElement vcgnode = info->doc->createElement("node");
			vcgnode.setAttribute("id","vcg-node");
			vcgnode.setAttribute("name","vcg-untitled");

			/*QDomNodeList instgeo = info->doc->elementsByTagName("instance_geometry");
			for(int jj = 0;jj < instgeo.size();++jj)
			{
				if (!instgeo.at(jj).isNull())
				{
					QDomNode par = instegeo.at(jj).parent();
					par
			}*/

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
			
			removeChildNode(geolib.at(0),"geometry","id","vcg-mesh-lib");
			QDomElement geonode = info->doc->createElement("geometry");
			geonode.setAttribute("id","vcg-mesh-lib");
			geonode.setAttribute("name","vcg-mesh");
			
			QDomElement meshnode = info->doc->createElement("mesh");
			
			QDomElement srcposnode = info->doc->createElement("source");
			srcposnode.setAttribute("id","vcg-mesh-positions");
			srcposnode.setAttribute("name","vcg-mesh-positions");

			QDomElement arrayposnode = info->doc->createElement("float_array");
			arrayposnode.setAttribute("id","vcg-mesh-positions-array");
			
			QString arrp;
			arrp.reserve(8 * m.vert.size());
			QString arrn;
			arrn.reserve(8 * m.vert.size());
			int nvert = 0;
			for(SaveMeshType::VertexIterator it = m.vert.begin();it != m.vert.end();++it)
			{
				if (!(it->IsD()))
				{
					arrp.append(QString::number(it->P().X()) + " " +QString::number(it->P().Y()) + " " + QString::number(it->P().Z()) + " ");
					arrn.append(QString::number(it->N().X()) + " " + QString::number(it->N().Y()) + " " + QString::number(it->N().Z())+ " ");
					//arrp.append(QString::number(it->P().X()).append(" ").append(QString::number(it->P().Y())).append(" ").append(QString::number(it->P().Z())).append(" "));
					//arrp.append(QString::number(it->N().X()).append(" ").append(QString::number(it->N().Y())).append(" ").append(QString::number(it->N().Z())).append(" "));
					++nvert;
				}
			}
			arrayposnode.setAttribute("count",QString::number(nvert * 3));
			QDomText ap = info->doc->createTextNode(arrp);

			QDomElement technode = info->doc->createElement("technique_common");
			QDomElement accnode = info->doc->createElement("accessor");
			accnode.setAttribute("source","#vcg-mesh-positions-array");
			accnode.setAttribute("count",QString::number(nvert));
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
			
			QDomElement srcnmnode = info->doc->createElement("source");
			srcnmnode.setAttribute("id","vcg-mesh-normals");
			srcnmnode.setAttribute("name","vcg-mesh-normals");

			QDomElement arraynmnode = info->doc->createElement("float_array");
			arraynmnode.setAttribute("id","vcg-mesh-normals-array");
			arraynmnode.setAttribute("count",QString::number(nvert * 3));

			QDomElement technode2 = info->doc->createElement("technique_common");
			QDomElement accnode2 = info->doc->createElement("accessor");
			accnode2.setAttribute("source","#vcg-mesh-normals-array");
			accnode2.setAttribute("count",QString::number(nvert));
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

			QDomElement vert = info->doc->createElement("vertices");
			vert.setAttribute("id","vcg-mesh-vertices");
			QDomElement vinp_pos = info->doc->createElement("input");
			vinp_pos.setAttribute("semantic","POSITION");
			vinp_pos.setAttribute("source","#vcg-mesh-positions");
			QDomElement vinp_nm = info->doc->createElement("input");
			vinp_nm.setAttribute("semantic","NORMAL");
			vinp_nm.setAttribute("source","#vcg-mesh-normals");

			vert.appendChild(vinp_pos);
			vert.appendChild(vinp_nm);

			meshnode.appendChild(vert);
			
			QDomElement tri = info->doc->createElement("triangles");
			

			QDomElement tinp_vert = info->doc->createElement("input");
			tinp_vert.setAttribute("offset","0");
			tinp_vert.setAttribute("semantic","VERTEX");
			tinp_vert.setAttribute("source","#vcg-mesh-vertices");
			QDomElement poly = info->doc->createElement("p");
			QString triangles_tess;
			int nface = 0;
			triangles_tess.reserve(9*m.face.size());
			for(SaveMeshType::FaceIterator itf = m.face.begin();itf != m.face.end();++itf)
			{
				if (!(itf->IsD()))
				{
					for(unsigned int ii = 0;ii < 3;++ii)
					{
						int ind_v = (*itf).V(ii) - &(m.vert[0]);
						if (triangles_tess == "")
							triangles_tess = QString::number(ind_v);
						else triangles_tess = triangles_tess.append(" ").append(QString::number(ind_v));
					}
					++nface;
				}
			}
			tri.setAttribute("count",nface);

			QDomText tri_list = info->doc->createTextNode(triangles_tess);
			poly.appendChild(tri_list);
			tri.appendChild(tinp_vert);
			tri.appendChild(poly);
			meshnode.appendChild(tri);
			geonode.appendChild(meshnode);
			geolib.at(0).appendChild(geonode);
		}

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