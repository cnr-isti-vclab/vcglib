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
private:
	static void CreateVertInput(QDomDocument& doc,QDomNode& vert,const QString& attr,const QString& ref)
	{
		QDomElement vinp_pos = doc.createElement("input");
		vinp_pos.setAttribute("semantic",attr);
		vinp_pos.setAttribute("source",ref);
		vert.appendChild(vinp_pos);
	}

	static void CreateFaceInput(QDomDocument& doc,QDomNode& tri,const QString& attr,const QString& ref,const int offset)
	{
		QDomElement tinp_vert = doc.createElement("input");
		tinp_vert.setAttribute("offset",QString::number(offset));
		tinp_vert.setAttribute("semantic",attr);
		tinp_vert.setAttribute("source",ref);
		tri.appendChild(tinp_vert);
	}

	static void CreateSource(QDomDocument& doc,QDomNode& meshnode,const QString& attr,const QDomText& val,int nvert)
	{
		int nel;
		std::vector<QString> coord;
		if ((attr == "positions") || (attr == "normals") || (attr == "wnornals"))
		{
			nel = 3;
			coord.push_back("X");
			coord.push_back("Y");
			coord.push_back("Z");
		}
		else
		{
			if (attr == "colors")
			{
				nel = 4;
				coord.push_back("R");
				coord.push_back("G");
				coord.push_back("B");
				coord.push_back("A");
			}
			else
			{
				nel = 2;
				coord.push_back("U");
				coord.push_back("V");
			}
		}
		QDomElement srcnmnode = doc.createElement("source");
		srcnmnode.setAttribute("id","vcg-mesh-"+attr);
		srcnmnode.setAttribute("name","vcg-mesh-"+attr);

		QDomElement arraynmnode = doc.createElement("float_array");
		arraynmnode.setAttribute("id","vcg-mesh-"+attr+"-array");
		arraynmnode.setAttribute("count",QString::number(nvert * nel));

		QDomElement technode2 = doc.createElement("technique_common");
		QDomElement accnode2 = doc.createElement("accessor");
		accnode2.setAttribute("source","#vcg-mesh-"+attr+"-array");
		accnode2.setAttribute("count",QString::number(nvert));
		accnode2.setAttribute("stride",QString::number(nel));
		
		for(int jj = 0; jj < nel;++jj)
		{
			QDomElement parxnode2 = doc.createElement("param");
			parxnode2.setAttribute("name",coord[jj]);
			parxnode2.setAttribute("type","float");
			accnode2.appendChild(parxnode2);
		}		

		technode2.appendChild(accnode2);
		arraynmnode.appendChild(val);
		srcnmnode.appendChild(arraynmnode);
		srcnmnode.appendChild(technode2);		
		meshnode.appendChild(srcnmnode);
	}
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
		arrp.reserve(10 * m.vert.size());
		QString arrn;
		arrn.reserve(10 * m.vert.size());
		
		int nvert = 0;
		for(SaveMeshType::VertexIterator it = m.vert.begin();it != m.vert.end();++it)
		{
			if (!(it->IsD()))
			{
				arrp.append(QString::number(it->P().X()).append(" ").append(QString::number(it->P().Y())).append(" ").append(QString::number(it->P().Z())).append(" "));
				arrn.append(QString::number(it->P().X()).append(" ").append(QString::number(it->P().Y())).append(" ").append(QString::number(it->P().Z())).append(" "));
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
		triangles_tess.reserve(10*m.face.size());
		for(SaveMeshType::FaceIterator itf = m.face.begin();itf != m.face.end();++itf)
		{
			if (!(itf->IsD()))
			{
				for(unsigned int ii = 0;ii < 3;++ii)
				{
					int ind_v = (*itf).V(ii) - &(m.vert[0]);
					if (triangles_tess == "")
						triangles_tess = QString::number(ind_v);
					else triangles_tess.append(" ").append(QString::number(ind_v));
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

	
	static int Save(SaveMeshType &m, const char * filename,AdditionalInfo*& in, const int &mask = -1)
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
			
			
			QString arrp;
			arrp.reserve(10 * 3 * m.vert.size());
			QString arrn;
			if(mask & vcg::tri::io::Mask::IOM_VERTNORMAL | mask & vcg::tri::io::Mask::IOM_WEDGNORMAL)
				arrn.reserve(10 * 3 * m.vert.size());
			QString arrt;
			if(mask & vcg::tri::io::Mask::IOM_VERTTEXCOORD)
				arrt.reserve(10 * 2 * m.vert.size());
			QString arrc;
			if(mask & vcg::tri::io::Mask::IOM_VERTCOLOR)
				arrc.reserve(5 * 4 * m.vert.size());
			int nvert = 0;
			for(SaveMeshType::VertexIterator it = m.vert.begin();it != m.vert.end();++it)
			{
				if (!(it->IsD()))
				{
					arrp.append(QString::number(it->P().X()).append(" ").append(QString::number(it->P().Y())).append(" ").append(QString::number(it->P().Z())).append(" "));
					if(mask & vcg::tri::io::Mask::IOM_VERTNORMAL)
						arrn.append(QString::number(it->N().X()).append(" ").append(QString::number(it->N().Y())).append(" ").append(QString::number(it->N().Z())).append(" "));
					if(mask & vcg::tri::io::Mask::IOM_VERTTEXCOORD)
						arrt.append(QString::number(it->T().u()).append(" ").append(QString::number(it->T().v())).append(" "));
					if(mask & vcg::tri::io::Mask::IOM_VERTCOLOR)
						arrc.append(QString::number(it->C().X()).append(" ").append(QString::number(it->C().Y())).append(" ").append(QString::number(it->C().Z())).append(" ").append(QString::number(it->C().W())).append(" "));
					++nvert;
				}
			}

			QDomText ap = info->doc->createTextNode(arrp);
			CreateSource(*(info->doc),meshnode,"positions",ap,nvert);
			
			if(mask & vcg::tri::io::Mask::IOM_VERTNORMAL | mask & vcg::tri::io::Mask::IOM_WEDGNORMAL)
			{
				QDomText an = info->doc->createTextNode(arrn);
				CreateSource(*(info->doc),meshnode,"normals",an,nvert);
			}

			if(mask & vcg::tri::io::Mask::IOM_VERTTEXCOORD)
			{
				QDomText at = info->doc->createTextNode(arrt);
				CreateSource(*(info->doc),meshnode,"textcoords",at,nvert);
			}

			if(mask & vcg::tri::io::Mask::IOM_VERTCOLOR)
			{
				QDomText ac = info->doc->createTextNode(arrc);
				CreateSource(*(info->doc),meshnode,"colors",ac,nvert);
			}

			QDomElement vert = info->doc->createElement("vertices");
			vert.setAttribute("id","vcg-mesh-vertices");
			CreateVertInput(*(info->doc),vert,"POSITION","#vcg-mesh-positions");
			if(mask & vcg::tri::io::Mask::IOM_VERTNORMAL)
				CreateVertInput(*(info->doc),vert,"NORMAL","#vcg-mesh-normals");
			if(mask & vcg::tri::io::Mask::IOM_VERTCOLOR)
				CreateVertInput(*(info->doc),vert,"COLOR","#vcg-mesh-colors");
			if(mask & vcg::tri::io::Mask::IOM_VERTTEXCOORD)
				CreateVertInput(*(info->doc),vert,"TEXCOORD","#vcg-mesh-normals");
			meshnode.appendChild(vert);
			
			QDomElement tri = info->doc->createElement("triangles");
			

			CreateFaceInput(*(info->doc),tri,"VERTEX","#vcg-mesh-vertices",0);
			QDomElement poly = info->doc->createElement("p");
			int nface = 0;
			int nattr = 1;

			QString triangles_wn;
			if (mask & MeshModel::IOM_WEDGNORMAL)
			{
				triangles_wn.reserve(3* 10 * m.face.size());
				CreateFaceInput(*(info->doc),tri,"NORMAL","#vcg-mesh-wnormals",nattr);
				++nattr;
			}
			QString triangles_wt;
			if (mask & MeshModel::IOM_WEDGTEXCOORD)
			{
				triangles_wt.reserve(2 * 10 * m.face.size());
				CreateFaceInput(*(info->doc),tri,"TEXCOORD","#vcg-mesh-wtext",nattr);
				++nattr;
			}
			QString triangles_tess;
			triangles_tess.reserve(nattr * 3 * 10 * m.face.size());
			int wn = 0;
			int wt = 0;
			for(SaveMeshType::FaceIterator itf = m.face.begin();itf != m.face.end();++itf)
			{
				if (!(itf->IsD()))
				{
					for(unsigned int ii = 0;ii < 3;++ii)
					{
						int ind_v = (*itf).V(ii) - &(m.vert[0]);
						if (triangles_tess == "")
							triangles_tess = QString::number(ind_v);
						else triangles_tess.append(" ").append(QString::number(ind_v));
						if (mask & MeshModel::IOM_WEDGNORMAL)
						{
							triangles_tess.append(" ").append(QString::number(wn));
							++wn;
							triangles_wn.append(QString::number((*itf).WN(ii).X()).append(" ").append(QString::number((*itf).WN(ii).Y())).append(" ").append(QString::number((*itf).WN(ii).Z())).append(" "));
						}

						if (mask & MeshModel::IOM_WEDGTEXCOORD)
						{
							triangles_tess.append(" ").append(QString::number(wt));
							++wt;
							triangles_wt.append(QString::number((*itf).WT(ii).u()).append(" ").append(QString::number((*itf).WT(ii).v())).append(" "));
						}
					}
					++nface;
				}
			}

			tri.setAttribute("count",nface);
			if (mask & MeshModel::IOM_WEDGNORMAL)
			{
				QDomText wnt = info->doc->createTextNode(triangles_wn);	
				CreateSource(*(info->doc),meshnode,"wnormals",wnt,nface * 3);
			}

			if (mask & MeshModel::IOM_WEDGTEXCOORD)
			{
				QDomText wtt = info->doc->createTextNode(triangles_wt);	
				CreateSource(*(info->doc),meshnode,"wtext",wtt,nface * 3);
			}

			QDomText tri_list = info->doc->createTextNode(triangles_tess);
			poly.appendChild(tri_list);
			tri.appendChild(poly);
			meshnode.appendChild(tri);
			geonode.appendChild(meshnode);
			geolib.at(0).appendChild(geonode);
		}
		//else
		//{
		//	removeChildNode(scenelst,"instance_visual_scene");
		//	for(int vsscn = 0;vsscn < scenelst.size();++vsscn)
		//	{
		//		QString url = scenelst.at(vsscn).toElement().attribute("url");
		//	}

		//	QDomElement vsnode = info->doc->createElement("instance_visual_scene");
		//	vsnode.setAttribute("url","#vcg-scene-node");
		//	scenelst.at(0).appendChild(vsnode);

		//	
		//	int vsscene_size = vsscene.size();
		//	assert(vsscene.size() == 1);
		//	removeChildNode(vsscene,"visual_scene","id","vcg-scene-node");
		//	QDomElement vslnode = info->doc->createElement("visual_scene");
		//	vslnode.setAttribute("id","vcg-scene-node");
		//	vslnode.setAttribute("name","vcg-untitled");

		//	QDomElement vcgnode = info->doc->createElement("node");
		//	vcgnode.setAttribute("id","vcg-node");
		//	vcgnode.setAttribute("name","vcg-untitled");

		//	/*QDomNodeList instgeo = info->doc->elementsByTagName("instance_geometry");
		//	for(int jj = 0;jj < instgeo.size();++jj)
		//	{
		//		if (!instgeo.at(jj).isNull())
		//		{
		//			QDomNode par = instegeo.at(jj).parent();
		//			par
		//	}*/

		//	QDomElement vcginst = info->doc->createElement("instance_geometry");
		//	vcginst.setAttribute("url","#vcg-mesh-lib");
		//	vcgnode.appendChild(vcginst);
		//	vslnode.appendChild(vcgnode);
		//	vsscene.at(0).appendChild(vslnode);

		//	QDomNodeList geolib = info->doc->elementsByTagName("library_geometries");
		//	assert(geolib.size() == 1);
		//	//removeChildNode(geolib.at(0));
		//	
		//	/*QDomElement mshnode;
		//	mshnode.setTagName("mesh");*/
		//	
		//	removeChildNode(geolib.at(0),"geometry","id","vcg-mesh-lib");
		//	QDomElement geonode = info->doc->createElement("geometry");
		//	geonode.setAttribute("id","vcg-mesh-lib");
		//	geonode.setAttribute("name","vcg-mesh");
		//	
		//	QDomElement meshnode = info->doc->createElement("mesh");
		//	
		//	QDomElement srcposnode = info->doc->createElement("source");
		//	srcposnode.setAttribute("id","vcg-mesh-positions");
		//	srcposnode.setAttribute("name","vcg-mesh-positions");

		//	QDomElement arrayposnode = info->doc->createElement("float_array");
		//	arrayposnode.setAttribute("id","vcg-mesh-positions-array");
		//	
		//	QString arrp;
		//	arrp.reserve(10 * m.vert.size());
		//	QString arrn;
		//	arrn.reserve(10 * m.vert.size());
		//	int nvert = 0;
		//	for(SaveMeshType::VertexIterator it = m.vert.begin();it != m.vert.end();++it)
		//	{
		//		if (!(it->IsD()))
		//		{
		//			//arrp.append(QString::number(it->P().X()) + " " +QString::number(it->P().Y()) + " " + QString::number(it->P().Z()) + " ");
		//			//arrn.append(QString::number(it->N().X()) + " " + QString::number(it->N().Y()) + " " + QString::number(it->N().Z())+ " ");
		//			arrp.append(QString::number(it->P().X()).append(" ").append(QString::number(it->P().Y())).append(" ").append(QString::number(it->P().Z())).append(" "));
		//			arrp.append(QString::number(it->N().X()).append(" ").append(QString::number(it->N().Y())).append(" ").append(QString::number(it->N().Z())).append(" "));
		//			++nvert;
		//		}
		//	}
		//	arrayposnode.setAttribute("count",QString::number(nvert * 3));
		//	QDomText ap = info->doc->createTextNode(arrp);

		//	QDomElement technode = info->doc->createElement("technique_common");
		//	QDomElement accnode = info->doc->createElement("accessor");
		//	accnode.setAttribute("source","#vcg-mesh-positions-array");
		//	accnode.setAttribute("count",QString::number(nvert));
		//	accnode.setAttribute("stride","3");
		//	
		//	QDomElement parxnode = info->doc->createElement("param");
		//	parxnode.setAttribute("name","X");
		//	parxnode.setAttribute("type","float");
		//	QDomElement parynode = info->doc->createElement("param");
		//	parynode.setAttribute("name","Y");
		//	parynode.setAttribute("type","float");
		//	QDomElement parznode = info->doc->createElement("param");
		//	parznode.setAttribute("name","Z");
		//	parznode.setAttribute("type","float");

		//	accnode.appendChild(parxnode);
		//	accnode.appendChild(parynode);
		//	accnode.appendChild(parznode);
		//	technode.appendChild(accnode);
		//	arrayposnode.appendChild(ap);
		//	srcposnode.appendChild(arrayposnode);
		//	srcposnode.appendChild(technode);
		//	
		//	meshnode.appendChild(srcposnode);
		//	
		//	QDomElement srcnmnode = info->doc->createElement("source");
		//	srcnmnode.setAttribute("id","vcg-mesh-normals");
		//	srcnmnode.setAttribute("name","vcg-mesh-normals");

		//	QDomElement arraynmnode = info->doc->createElement("float_array");
		//	arraynmnode.setAttribute("id","vcg-mesh-normals-array");
		//	arraynmnode.setAttribute("count",QString::number(nvert * 3));

		//	QDomElement technode2 = info->doc->createElement("technique_common");
		//	QDomElement accnode2 = info->doc->createElement("accessor");
		//	accnode2.setAttribute("source","#vcg-mesh-normals-array");
		//	accnode2.setAttribute("count",QString::number(nvert));
		//	accnode2.setAttribute("stride","3");
		//	
		//	QDomElement parxnode2 = info->doc->createElement("param");
		//	parxnode2.setAttribute("name","X");
		//	parxnode2.setAttribute("type","float");
		//	QDomElement parynode2 = info->doc->createElement("param");
		//	parynode2.setAttribute("name","Y");
		//	parynode2.setAttribute("type","float");
		//	QDomElement parznode2 = info->doc->createElement("param");
		//	parznode2.setAttribute("name","Z");
		//	parznode2.setAttribute("type","float");

		//	QDomText an = info->doc->createTextNode(arrn);
		//	

		//	accnode2.appendChild(parxnode2);
		//	accnode2.appendChild(parynode2);
		//	accnode2.appendChild(parznode2);
		//	technode2.appendChild(accnode2);
		//	arraynmnode.appendChild(an);
		//	srcnmnode.appendChild(arraynmnode);
		//	srcnmnode.appendChild(technode2);
		//	
		//	meshnode.appendChild(srcnmnode);

		//	QDomElement vert = info->doc->createElement("vertices");
		//	vert.setAttribute("id","vcg-mesh-vertices");
		//	QDomElement vinp_pos = info->doc->createElement("input");
		//	vinp_pos.setAttribute("semantic","POSITION");
		//	vinp_pos.setAttribute("source","#vcg-mesh-positions");
		//	QDomElement vinp_nm = info->doc->createElement("input");
		//	vinp_nm.setAttribute("semantic","NORMAL");
		//	vinp_nm.setAttribute("source","#vcg-mesh-normals");

		//	vert.appendChild(vinp_pos);
		//	vert.appendChild(vinp_nm);

		//	meshnode.appendChild(vert);
		//	
		//	QDomElement tri = info->doc->createElement("triangles");
		//	

		//	QDomElement tinp_vert = info->doc->createElement("input");
		//	tinp_vert.setAttribute("offset","0");
		//	tinp_vert.setAttribute("semantic","VERTEX");
		//	tinp_vert.setAttribute("source","#vcg-mesh-vertices");
		//	QDomElement poly = info->doc->createElement("p");
		//	QString triangles_tess;
		//	int nface = 0;
		//	triangles_tess.reserve(10*m.face.size());
		//	for(SaveMeshType::FaceIterator itf = m.face.begin();itf != m.face.end();++itf)
		//	{
		//		if (!(itf->IsD()))
		//		{
		//			for(unsigned int ii = 0;ii < 3;++ii)
		//			{
		//				int ind_v = (*itf).V(ii) - &(m.vert[0]);
		//				if (triangles_tess == "")
		//					triangles_tess = QString::number(ind_v);
		//				else triangles_tess.append(" ").append(QString::number(ind_v));
		//			}
		//			++nface;
		//		}
		//	}
		//	tri.setAttribute("count",nface);

		//	QDomText tri_list = info->doc->createTextNode(triangles_tess);
		//	poly.appendChild(tri_list);
		//	tri.appendChild(tinp_vert);
		//	tri.appendChild(poly);
		//	meshnode.appendChild(tri);
		//	geonode.appendChild(meshnode);
		//	geolib.at(0).appendChild(geonode);
		//}

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
		capability |= MeshModel::IOM_VERTNORMAL;
		capability |= MeshModel::IOM_VERTTEXCOORD;
		capability |= MeshModel::IOM_VERTCOLOR;
		//capability |= MeshModel::
		////face
		////capability |= MeshModel::IOM_FACEFLAGS;
		////capability |= MeshModel::IOM_FACECOLOR;
		//capability |= MeshModel::IOM_FACENORMAL;

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