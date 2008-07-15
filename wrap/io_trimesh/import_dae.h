/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2008                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef __VCGLIB_IMPORTERDAE
#define __VCGLIB_IMPORTERDAE

//importer for collada's files

#include<wrap/dae/util_dae.h>

// uncomment one of the following line to enable the Verbose Trace for Crease
//#define QDEBUG if(1) ; else printf 
#define QDEBUG qDebug

namespace vcg {
namespace tri {
namespace io {

	template<typename OpenMeshType>
	class ImporterDAE : public UtilDAE
	{

	private:
		static int WedgeNormalAttribute(OpenMeshType& m,const QStringList face,const QStringList wn,const QDomNode wnsrc,const int meshfaceind,const int faceind,const int component)
		{
			int indnm = -1;
			if (!wnsrc.isNull())
			{
				indnm = face.at(faceind).toInt();
				assert(indnm * 3 < wn.size());
				m.face[meshfaceind].WN(component) = vcg::Point3f(wn.at(indnm * 3).toFloat(),wn.at(indnm * 3 + 1).toFloat(),wn.at(indnm * 3 + 2).toFloat());
			}
			return indnm;
		}

		static int WedgeTextureAttribute(OpenMeshType& m,const QStringList face,int ind_txt,const QStringList wt,const QDomNode wtsrc,const int meshfaceind,const int faceind,const int component,const int stride = 2)
		{
			int indtx = -1;
			if (!wtsrc.isNull())
			{
				indtx = face.at(faceind).toInt();
				int num = wt.size(); 
				assert(indtx * stride < wt.size());
				m.face[meshfaceind].WT(component) = vcg::TexCoord2<float>();
				m.face[meshfaceind].WT(component).U() = wt.at(indtx * stride).toFloat();
				m.face[meshfaceind].WT(component).V() = wt.at(indtx * stride + 1).toFloat();
				
				m.face[meshfaceind].WT(component).N() = ind_txt;
				
			}
			return indtx;
		}

		static int WedgeColorAttribute(OpenMeshType& m,const QStringList face,const QStringList wc,const QDomNode wcsrc,const int meshfaceind,const int faceind,const int component)
		{
			int indcl;
			if (!wcsrc.isNull())
			{
				indcl = face.at(faceind).toInt();
				assert(indcl * 4 < wc.size());
				m.face[meshfaceind].WC(component) = vcg::Color4b(wc.at(indcl * 4).toFloat(),wc.at(indcl * 4 + 1).toFloat(),wc.at(indcl * 4 + 2).toFloat(),wc.at(indcl * 4 + 3).toFloat());
			}
			return indcl;
		}

		static void FindStandardWedgeAttributes(WedgeAttribute& wed,const QDomNode nd,const QDomDocument doc)
		{
			wed.wnsrc = findNodeBySpecificAttributeValue(nd,"input","semantic","NORMAL");
			wed.offnm = findStringListAttribute(wed.wn,wed.wnsrc,nd,doc,"NORMAL");

			wed.wtsrc = findNodeBySpecificAttributeValue(nd,"input","semantic","TEXCOORD");
			if (!wed.wtsrc.isNull())
			{
				QDomNode src = attributeSourcePerSimplex(nd,doc,"TEXCOORD");
				if (isThereTag(src,"accessor"))
				{
					QDomNodeList wedatts = src.toElement().elementsByTagName("accessor");
					wed.stride = wedatts.at(0).toElement().attribute("stride").toInt();
				}
				else 
					wed.stride = 2;
			}
			else
				wed.stride = 2;

			wed.offtx = findStringListAttribute(wed.wt,wed.wtsrc,nd,doc,"TEXCOORD"); 

			wed.wcsrc = findNodeBySpecificAttributeValue(nd,"input","semantic","COLOR");
			wed.offcl = findStringListAttribute(wed.wc,wed.wcsrc,nd,doc,"COLOR"); 
		}
	
		static DAEError LoadPolygonalMesh(QDomNodeList& polypatch,OpenMeshType& m,const size_t offset,AdditionalInfoDAE* info)
		{
			return E_NOERROR;
		}

		static DAEError	LoadPolygonalListMesh(QDomNodeList& polylist,OpenMeshType& m,const size_t offset,AdditionalInfoDAE* info)
		{
			typedef PolygonalMesh< MyPolygon<typename OpenMeshType::VertexType> > PolyMesh;
			PolyMesh pm;
			
			//copying vertices 
			for(typename OpenMeshType::VertexIterator itv = m.vert.begin();itv != m.vert.end();++itv)
			{	
				vcg::Point3f p(itv->P().X(),itv->P().Y(),itv->P().Z());
				typename PolyMesh::VertexType v;
				v.P() = p;
				pm.vert.push_back(v);
			}

			int polylist_size = polylist.size();
			for(int pl = 0; pl < polylist_size;++pl)
			{ 
				QString mat =  polylist.at(pl).toElement().attribute(QString("material"));
				QDomNode txt_node = textureFinder(mat,*(info->dae->doc));
				int ind_txt = -1;
				if (!txt_node.isNull())
					ind_txt = indexTextureByImgNode(*(info->dae->doc),txt_node);

				//PolyMesh::PERWEDGEATTRIBUTETYPE att = PolyMesh::NONE;
				WedgeAttribute wa;
				FindStandardWedgeAttributes(wa,polylist.at(pl),*(info->dae->doc));
				QStringList vertcount;
				valueStringList(vertcount,polylist.at(pl),"vcount");
				int indforpol = findOffSetForASingleSimplex(polylist.at(pl));
				int offpols = 0;
				int npolig = vertcount.size();
				QStringList polyind;
				valueStringList(polyind,polylist.at(pl),"p");
				for(int ii = 0;ii < npolig;++ii)
				{
					int nvert = vertcount.at(ii).toInt();
					typename PolyMesh::FaceType p(nvert);
				
					for(int iv = 0;iv < nvert;++iv)
					{
						int index = offset + polyind.at(offpols + iv * indforpol).toInt();
						p._pv[iv] = &(pm.vert[index]);
						int nmindex = -1;

						if (!wa.wnsrc.isNull())
							nmindex = offset + polyind.at(offpols + iv * indforpol + wa.offnm).toInt();

						int txindex = -1;
						if (!wa.wtsrc.isNull())
						{
							txindex = offset + polyind.at(offpols + iv * indforpol + wa.offtx).toInt();
							/*p._txc[iv].U() = wa.wt.at(txindex * 2).toFloat();
							p._txc[iv].V() = wa.wt.at(txindex * 2 + 1).toFloat();
							p._txc[iv].N() = ind_txt;*/
						}
					}
					pm._pols.push_back(p);
					offpols += nvert * indforpol;
				}
			}
			pm.triangulate(m);
			return E_NOERROR;
		}
		
		static DAEError LoadTriangularMesh(QDomNodeList& tripatch,OpenMeshType& m,const size_t offset,AdditionalInfoDAE* info)
		{
			int tripatch_size = tripatch.size();
			for(int tript = 0; tript < tripatch_size;++tript)
			{
				QString mat =  tripatch.at(tript).toElement().attribute(QString("material"));
				QDomNode txt_node = textureFinder(mat,*(info->dae->doc));
				int ind_txt = -1;
				if (!txt_node.isNull())
					ind_txt = indexTextureByImgNode(*(info->dae->doc),txt_node);
				
				int nfcatt = tripatch.at(tript).toElement().elementsByTagName("input").size();

				QStringList face;
				valueStringList(face,tripatch.at(tript),"p");
				int face_size = face.size();
				int offsetface = (int)m.face.size();
				if (face_size != 0) 
				{	
					vcg::tri::Allocator<OpenMeshType>::AddFaces(m,face_size / (nfcatt * 3));
					WedgeAttribute wa;
					FindStandardWedgeAttributes(wa,tripatch.at(tript),*(info->dae->doc));

					int jj = 0;	
					//int dd = m.face.size();
					for(int ff = offsetface;ff < (int) m.face.size();++ff)
					{ 
						
						for(unsigned int tt = 0;tt < 3;++tt)
						{
							int indvt = face.at(jj).toInt();
							assert(indvt + offset < m.vert.size());
							m.face[ff].V(tt) = &(m.vert[indvt + offset]);

							if(tri::HasPerWedgeNormal(m))  
									WedgeNormalAttribute(m,face,wa.wn,wa.wnsrc,ff,jj + wa.offnm,tt);
							if(tri::HasPerWedgeTexCoord(m) && ind_txt != -1)
							{
								WedgeTextureAttribute(m,face,ind_txt,wa.wt,wa.wtsrc,ff,jj + wa.offtx,tt,wa.stride);
							}
							if(tri::HasPerWedgeColor(m))
									WedgeColorAttribute(m,face,wa.wc,wa.wcsrc,ff,jj + wa.offcl,tt);

							jj += nfcatt;
						}
					}
				}
			}
			return E_NOERROR;
		}

		static int LoadControllerMesh(OpenMeshType& m, AdditionalInfoDAE* info, const QDomElement& geo, CallBackPos *cb=0)
		{
			assert(geo.tagName() == "controller");
			QDomNodeList skinList = geo.toElement().elementsByTagName("skin");
			if(skinList.size()!=1) return E_CANTOPEN;
			QDomElement skinNode = skinList.at(0).toElement();
			
			QString geomNode_url;
			referenceToANodeAttribute(skinNode,"source",geomNode_url);
			QDEBUG("Found a controller referencing a skin with url '%s'", qPrintable(geomNode_url));
			QDomNode refNode = findNodeBySpecificAttributeValue(*(info->dae->doc),"geometry","id",geomNode_url);
			LoadMesh(m, info, refNode.toElement());
		}						
		
    /*
		 Basic function that fills a mesh with the coord, norm and tristarting from node that is of kind <geometry>
		 */
		static int LoadMesh(OpenMeshType& m, AdditionalInfoDAE* info, const QDomElement& geo, CallBackPos *cb=0)
		{
			assert(geo.tagName() == "geometry");
			if (isThereTag(geo,"mesh"))
			{
				if ((cb !=NULL) && (((info->numvert + info->numface)%100)==0) && !(*cb)((100*(info->numvert + info->numface))/(info->numvert + info->numface), "Vertex Loading"))
					return E_CANTOPEN;
				/*QDomNodeList geosrc = geo.toElement().elementsByTagName("source");
				int geosrc_size = geosrc.size();
				if (geosrc_size < 1)
				return E_NOVERTEXPOSITION;*/
        QDEBUG("**** Loading a Geometry Mesh ****");
				QDomNodeList vertices = geo.toElement().elementsByTagName("vertices");
				if (vertices.size() != 1) return E_INCOMPATIBLECOLLADA141FORMAT;
				QDomElement vertNode = vertices.at(0).toElement();

				QDomNode srcnode = attributeSourcePerSimplex(vertNode,*(info->dae->doc),"POSITION");
				if (srcnode.isNull()) return E_NOVERTEXPOSITION;

				QStringList geosrcposarr;
				valueStringList(geosrcposarr,srcnode,"float_array");

				int geosrcposarr_size = geosrcposarr.size();
				if ((geosrcposarr_size % 3) != 0)
					return E_CANTOPEN;
				int nvert = geosrcposarr_size / 3;
				size_t offset = m.vert.size();
				if (geosrcposarr_size != 0)
				{
					vcg::tri::Allocator<OpenMeshType>::AddVertices(m,nvert);

					QDomNode srcnodenorm = attributeSourcePerSimplex(vertices.at(0),*(info->dae->doc),"NORMAL");
					QStringList geosrcvertnorm;
					if (!srcnodenorm.isNull())
						valueStringList(geosrcvertnorm,srcnodenorm,"float_array");

					QDomNode srcnodetext = attributeSourcePerSimplex(vertices.at(0),*(info->dae->doc),"TEXCOORD");
					QStringList geosrcverttext;
					if (!srcnodetext.isNull())
						valueStringList(geosrcverttext,srcnodetext,"float_array");

					QDomNode srcnodecolor = attributeSourcePerSimplex(vertices.at(0),*(info->dae->doc),"COLOR");
					QStringList geosrcvertcol;
					if (!srcnodecolor.isNull())
						valueStringList(geosrcvertcol,srcnodecolor,"float_array");

					int ii = 0;
					for(size_t vv = offset;vv < m.vert.size();++vv)
					{						
						Point3f positionCoord(geosrcposarr[ii * 3].toFloat(),geosrcposarr[ii * 3 + 1].toFloat(),geosrcposarr[ii * 3 + 2].toFloat());
						m.vert[vv].P() = positionCoord;

						if (!srcnodenorm.isNull())
						{
							Point3f normalCoord(geosrcvertnorm[ii * 3].toFloat(),
																	geosrcvertnorm[ii * 3 + 1].toFloat(),
																	geosrcvertnorm[ii * 3 + 2].toFloat());
							normalCoord.Normalize();
							m.vert[vv].N() = normalCoord;
						}

						/*if (!srcnodecolor.isNull())
						{
						assert((ii * 4 < geosrcvertcol.size()) && (ii * 4 + 1 < geosrcvertcol.size()) && (ii * 4 + 2 < geosrcvertcol.size()) && (ii * 4 + 1 < geosrcvertcol.size()));
						m.vert[vv].C() = vcg::Color4b(geosrcvertcol[ii * 4].toFloat(),geosrcvertcol[ii * 4 + 1].toFloat(),geosrcvertcol[ii * 4 + 2].toFloat(),geosrcvertcol[ii * 4 + 3].toFloat());
						}*/

						if (!srcnodetext.isNull())
						{

							assert((ii * 2 < geosrcverttext.size()) && (ii * 2 + 1 < geosrcverttext.size()));
							m.vert[vv].T() = vcg::TexCoord2<float>();
							m.vert[vv].T().u() = geosrcverttext[ii * 2].toFloat();
							m.vert[vv].T().v() = geosrcverttext[ii * 2 + 1].toFloat();
						}
						++ii;
					}

					QDomNodeList tripatch = geo.toElement().elementsByTagName("triangles");
					int tripatch_size = tripatch.size();
					QDomNodeList polypatch = geo.toElement().elementsByTagName("polygons");
					int polypatch_size = polypatch.size();
					QDomNodeList polylist = geo.toElement().elementsByTagName("polylist");
					int polylist_size = polylist.size();
					if ((tripatch_size == 0) && (polypatch_size == 0) && (polylist_size == 0))
						return E_NOPOLYGONALMESH;
					
					DAEError err = E_NOERROR;
					if (tripatch_size != 0) 
						err = LoadTriangularMesh(tripatch,m,offset,info);
					else 
						if (polypatch_size != 0) 
							err = LoadPolygonalMesh(polypatch,m,offset,info);
						else
							if (polylist_size != 0) 
								err = LoadPolygonalListMesh(polylist,m,offset,info);
					if (err != E_NOERROR) 
						return err;
				}
				return E_NOERROR;
			}
			else return E_NOMESH;
		}

		static void GetTexCoord(const QDomDocument& doc,AdditionalInfoDAE* inf)
		{
			QDomNodeList txlst = doc.elementsByTagName("library_images");
			for(int img = 0;img < txlst.at(0).childNodes().size();++img)
			{
				QDomNodeList nlst = txlst.at(0).childNodes().at(img).toElement().elementsByTagName("init_from");
				if (nlst.size() > 0)
				{
					inf->texturefile.push_back(nlst.at(0).firstChild().nodeValue());
				}
			}
		}

	// This recursive function add to a mesh the subtree starting from the passed node. 
	// When you start from a visual_scene, you can find nodes. 
  // nodes can be directly instanced or referred from the node library.
		
		static void AddNodeToMesh(QDomElement node, 
															OpenMeshType& m, Matrix44f curTr,
															AdditionalInfoDAE*& info)
		{
				QDEBUG("Starting processing <node> with id %s",qPrintable(node.attribute("id")));
 
				curTr = curTr * getTransfMatrixFromNode(node);
				
				QDomNodeList geomNodeList = node.elementsByTagName("instance_geometry");
				for(int ch = 0;ch < geomNodeList.size();++ch) 
				{
					QDomElement instGeomNode= geomNodeList.at(ch).toElement();
					if(instGeomNode.parentNode()==node) // process only direct child
					{
						QDEBUG("Found a instance_geometry with url %s",qPrintable(instGeomNode.attribute("url")));
						
						QString geomNode_url;
						referenceToANodeAttribute(instGeomNode,"url",geomNode_url);
						QDEBUG("Found a instance_geometry with url '%s'", qPrintable(geomNode_url));
						QDomNode refNode = findNodeBySpecificAttributeValue(*(info->dae->doc),"geometry","id",geomNode_url);
						
						OpenMeshType newMesh;
						LoadMesh(newMesh, info, refNode.toElement());
						tri::UpdatePosition<OpenMeshType>::Matrix(newMesh,curTr);
						tri::Append<OpenMeshType,OpenMeshType>::Mesh(m,newMesh);
					}
				}
				
				QDomNodeList controllerNodeList = node.elementsByTagName("instance_controller");
				for(int ch = 0;ch < controllerNodeList.size();++ch) 
				{
					QDomElement instContrNode= controllerNodeList.at(ch).toElement();
					if(instContrNode.parentNode()==node) // process only direct child
					{
						QDEBUG("Found a instance_controller with url %s",qPrintable(instContrNode.attribute("url")));
					
						QString controllerNode_url;
						referenceToANodeAttribute(instContrNode,"url",controllerNode_url);
						QDEBUG("Found a instance_controller with url '%s'", qPrintable(controllerNode_url));
						QDomNode refNode = findNodeBySpecificAttributeValue(*(info->dae->doc),"controller","id",controllerNode_url);
						
						OpenMeshType newMesh;
						LoadControllerMesh(newMesh, info, refNode.toElement());
						tri::UpdatePosition<OpenMeshType>::Matrix(newMesh,curTr);
						tri::Append<OpenMeshType,OpenMeshType>::Mesh(m,newMesh);
					}
				}
								
				QDomNodeList nodeNodeList = node.elementsByTagName("node");
				for(int ch = 0;ch < nodeNodeList.size();++ch)
				{
					if(nodeNodeList.at(ch).parentNode()==node) // process only direct child
							AddNodeToMesh(nodeNodeList.at(ch).toElement(), m,curTr, info);
				}
				
				QDomNodeList instanceNodeList = node.elementsByTagName("instance_node");
				for(int ch = 0;ch < instanceNodeList.size();++ch)
				{
					if(instanceNodeList.at(ch).parentNode()==node) // process only direct child
					{
						QDomElement instanceNode =  instanceNodeList.at(ch).toElement();
						QString node_url;
						referenceToANodeAttribute(instanceNode,"url",node_url);
						QDEBUG("Found a instance_node with url '%s'", qPrintable(node_url));
						QDomNode refNode = findNodeBySpecificAttributeValue(*(info->dae->doc),"node","id",node_url);
						if(refNode.isNull()) 
							QDEBUG("findNodeBySpecificAttributeValue returned a null node for %s",qPrintable(node_url));
										 
						AddNodeToMesh(refNode.toElement(), m,curTr, info);
					}
				}
		}
		

// Retrieve the transformation matrix that is defined in the childs of a node.
// used during the recursive descent.
static Matrix44f getTransfMatrixFromNode(const QDomElement parentNode)
{
	QDEBUG("getTrans form node with tag %s",qPrintable(parentNode.tagName()));
	assert(parentNode.tagName() == "node");
	
	std::vector<QDomNode> rotationList;
	QDomNode matrixNode;
	QDomNode translationNode;
	for(int ch = 0;ch < parentNode.childNodes().size();++ch)
		{
			if (parentNode.childNodes().at(ch).nodeName() == "rotate")    
				rotationList.push_back(parentNode.childNodes().at(ch));
			if (parentNode.childNodes().at(ch).nodeName() == "translate")	
				translationNode = parentNode.childNodes().at(ch);							
			if (parentNode.childNodes().at(ch).nodeName() == "matrix")	  
				matrixNode = parentNode.childNodes().at(ch);							
		}

		Matrix44f rotM;		   rotM.SetIdentity();
		Matrix44f transM; transM.SetIdentity();

		if (!translationNode.isNull()) ParseTranslation(transM,translationNode);
		if (!rotationList.empty()) ParseRotationMatrix(rotM,rotationList);
		if (!matrixNode.isNull()) 
		{
			ParseMatrixNode(transM,matrixNode);
		  return transM;
		}
	  return transM*rotM;
}

	public:

		//merge all meshes in the collada's file in the templeted mesh m
		//I assume the mesh 
		
		static int Open(OpenMeshType& m,const char* filename,AdditionalInfo*& info, CallBackPos *cb=0)
		{
			QDEBUG("----- Starting the processing of %s ------",filename);
			AdditionalInfoDAE* inf = new AdditionalInfoDAE();
			inf->dae = new InfoDAE(); 
			
			QDomDocument* doc = new QDomDocument(filename);
			QFile file(filename);
			if (!file.open(QIODevice::ReadOnly))
				return E_CANTOPEN;
			if (!doc->setContent(&file)) 
			{
				file.close();
				return E_CANTOPEN;
			}
			file.close();
			
			inf->dae->doc = doc;
			//GetTexture(*(info->doc),inf);
			
			//scene->instance_visual_scene
			QDomNodeList scenes = inf->dae->doc->elementsByTagName("scene");
			int scn_size = scenes.size();
			if (scn_size == 0) 
				return E_NO3DSCENE;
			QDEBUG("File Contains %i Scenes",scenes.size());
			int problem = E_NOERROR;
			bool found_a_mesh = false;
			//Is there geometry in the file? 
			bool geoinst_found = false;
			
			// The main loading loop
			// for each scene in COLLADA FILE
			/*
			Some notes on collada structure.
			the top root is the <scene> that can contains one of more <visual_scene>.
			<visual_scene> can be directly written there (check!) or instanced from their definition in the <library_visual_scene>
			each <visual_scene> contains a hierarchy of <node>	
		  each <node> contains
				transformation
				other node
				instance of geometry 
				instance of controller
			*/
			for(int scn = 0;scn < scn_size;++scn)
			{
				QDomNodeList instscenes = scenes.at(scn).toElement().elementsByTagName("instance_visual_scene");
				int instscn_size = instscenes.size();
				QDEBUG("Scene %i contains %i instance_visual_scene ",scn,instscn_size);
				if (instscn_size == 0)  return E_INCOMPATIBLECOLLADA141FORMAT;

				//for each scene instance in a COLLADA scene
				for(int instscn = 0;instscn < instscn_size; ++instscn)
				{
					QString libscn_url;
					referenceToANodeAttribute(instscenes.at(instscn),"url",libscn_url);	
					QDEBUG("instance_visual_scene %i refers %s ",instscn,qPrintable(libscn_url));
					
					QDomNode nd = QDomNode(*(inf->dae->doc));
					QDomNode visscn = findNodeBySpecificAttributeValue(*(inf->dae->doc),"visual_scene","id",libscn_url);
					if(visscn.isNull()) return E_UNREFERENCEBLEDCOLLADAATTRIBUTE;
					
					//assert (visscn.toElement().Attribute("id") == libscn_url);
					//for each node in the libscn_url visual scene  
					QDomNodeList visscn_child = visscn.childNodes();
					QDEBUG("instance_visual_scene %s has %i children",qPrintable(libscn_url),visscn_child.size());
					
					// for each direct child of a visual scene process it
					for(int chdind = 0; chdind < visscn_child.size();++chdind)
					{
						QDomElement node=visscn_child.at(chdind).toElement();
						if(node.isNull()) continue;
						QDEBUG("Processing Visual Scene child %i - of type '%s'",chdind,qPrintable(node.tagName()));
						Matrix44f baseTr; baseTr.SetIdentity();
						
						if(node.toElement().tagName()=="node")
							AddNodeToMesh(node.toElement(), m, baseTr,inf);
					}	// end for each node of a given scene				
				} // end for each visual scene instance
			} // end for each scene instance 
			return problem;
		}

		static bool LoadMask(const char * filename, AdditionalInfoDAE*& addinfo)
		{
			bool bHasPerWedgeTexCoord = false;
			bool bHasPerWedgeNormal		= false;
			bool bHasPerVertexColor		= false;
			bool bHasPerFaceColor			= false;
			bool bHasPerVertexNormal = false;
			bool bHasPerVertexText = false;
			
			AdditionalInfoDAE* info = new AdditionalInfoDAE();
			info->dae = new InfoDAE(); 
			

			QDomDocument* doc = new QDomDocument(filename);
			QFile file(filename);
			if (!file.open(QIODevice::ReadOnly))
				return false;
			if (!doc->setContent(&file)) 
			{
				file.close();
				return false;
			}
			file.close();
			

			info->dae->doc = doc;
			GetTexCoord(*(info->dae->doc),info);
			QDomNodeList scenes = info->dae->doc->elementsByTagName("scene");
			int scn_size = scenes.size();
			

			//Is there geometry in the file? 
			bool geoinst_found = false;
			//for each scene in COLLADA FILE
			for(int scn = 0;scn < scn_size;++scn)
			{
				QDomNodeList instscenes = scenes.at(scn).toElement().elementsByTagName("instance_visual_scene");
				int instscn_size = instscenes.size();
				if (instscn_size == 0)  return false;

				//for each scene instance in a COLLADA scene
				for(int instscn = 0;instscn < instscn_size; ++instscn)
				{
					QString libscn_url;
					referenceToANodeAttribute(instscenes.at(instscn),"url",libscn_url);	
					QDomNode nd = QDomNode(*(info->dae->doc));
					QDomNode visscn = findNodeBySpecificAttributeValue(*(info->dae->doc),"visual_scene","id",libscn_url);
					if(visscn.isNull()) 	return false;
					
					//for each node in the libscn_url visual scene  
					//QDomNodeList& visscn_child = visscn.childNodes();
					QDomNodeList visscn_child = visscn.childNodes();
					
					//for each direct child of a libscn_url visual scene find if there is some geometry instance
					for(int chdind = 0; chdind < visscn_child.size();++chdind)
					{
						//QDomNodeList& geoinst = visscn_child.at(chdind).toElement().elementsByTagName("instance_geometry");
						QDomNodeList geoinst = visscn_child.at(chdind).toElement().elementsByTagName("instance_geometry");
						int geoinst_size = geoinst.size();
						if (geoinst_size != 0)
						{
							
							geoinst_found |= true;
							QDomNodeList geolib = info->dae->doc->elementsByTagName("library_geometries");
							assert(geolib.size() == 1);
							//!!!!!!!!!!!!!!!!!here will be the code for geometry transformations!!!!!!!!!!!!!!!!!!!!!!
							info->numvert = 0;
							info->numface = 0;
							for(int geoinst_ind = 0;geoinst_ind < geoinst_size;++geoinst_ind)
							{
								QString geo_url;
								referenceToANodeAttribute(geoinst.at(geoinst_ind),"url",geo_url);
								
								QDomNode geo = findNodeBySpecificAttributeValue(geolib.at(0),"geometry","id",geo_url);
								if (geo.isNull())
									return false;
							
								QDomNodeList vertlist = geo.toElement().elementsByTagName("vertices");

								for(int vert = 0;vert < vertlist.size();++vert)
								{
									QDomNode no;
									no = findNodeBySpecificAttributeValue(vertlist.at(vert),"input","semantic","POSITION");
									QString srcurl;
									referenceToANodeAttribute(no,"source",srcurl);
									no = findNodeBySpecificAttributeValue(geo,"source","id",srcurl);
									QDomNodeList fa = no.toElement().elementsByTagName("float_array");
									assert(fa.size() == 1);
									info->numvert += (fa.at(0).toElement().attribute("count").toInt() / 3);
									no = findNodeBySpecificAttributeValue(vertlist.at(vert),"input","semantic","COLOR");									
									if (!no.isNull()) 
										bHasPerVertexColor = true;
									no = findNodeBySpecificAttributeValue(vertlist.at(vert),"input","semantic","NORMAL");									
									if (!no.isNull()) 
										bHasPerVertexNormal = true;
									no = findNodeBySpecificAttributeValue(vertlist.at(vert),"input","semantic","TEXCOORD");									
									if (!no.isNull()) 
										bHasPerVertexText = true;
								}

								const char* arr[] = {"triangles","polylist","polygons"};

								for(unsigned int tt= 0;tt < 3;++tt)
								{
									QDomNodeList facelist = geo.toElement().elementsByTagName(arr[tt]);
									for(int face = 0;face < facelist.size();++face)
									{
										info->numface += facelist.at(face).toElement().attribute("count").toInt() ;
										QDomNode no;
										no = findNodeBySpecificAttributeValue(facelist.at(face),"input","semantic","NORMAL");
										if (!no.isNull()) 
											bHasPerWedgeNormal = true;
										no = findNodeBySpecificAttributeValue(facelist.at(face),"input","semantic","TEXCOORD");
										if (!no.isNull()) 
											bHasPerWedgeTexCoord = true;
									}
								}
							}
						}
					}
				}
			}
			
			if (!geoinst_found)
			{
				QDomNodeList geolib = info->dae->doc->elementsByTagName("library_geometries");
				assert(geolib.size() == 1);
				QDomNodeList geochild = geolib.at(0).toElement().elementsByTagName("geometry");
				//!!!!!!!!!!!!!!!!!here will be the code for geometry transformations!!!!!!!!!!!!!!!!!!!!!!
				info->numvert = 0;
				info->numface = 0;
				for(int geoinst_ind = 0;geoinst_ind < geochild.size();++geoinst_ind)
				{
					QDomNodeList vertlist = geochild.at(geoinst_ind).toElement().elementsByTagName("vertices");

					for(int vert = 0;vert < vertlist.size();++vert)
					{
						QDomNode no;
						no = findNodeBySpecificAttributeValue(vertlist.at(vert),"input","semantic","POSITION");
						QString srcurl;
						referenceToANodeAttribute(no,"source",srcurl);
						no = findNodeBySpecificAttributeValue(geochild.at(geoinst_ind),"source","id",srcurl);
						QDomNodeList fa = no.toElement().elementsByTagName("float_array");
						assert(fa.size() == 1);
						info->numvert += (fa.at(0).toElement().attribute("count").toInt() / 3);
						no = findNodeBySpecificAttributeValue(vertlist.at(vert),"input","semantic","COLOR");									
						if (!no.isNull()) 
							bHasPerVertexColor = true;
						no = findNodeBySpecificAttributeValue(vertlist.at(vert),"input","semantic","NORMAL");									
						if (!no.isNull()) 
							bHasPerVertexNormal = true;
						no = findNodeBySpecificAttributeValue(vertlist.at(vert),"input","semantic","TEXCOORD");									
						if (!no.isNull()) 
							bHasPerVertexText = true;
					}

					QDomNodeList facelist = geochild.at(geoinst_ind).toElement().elementsByTagName("triangles");
					for(int face = 0;face < facelist.size();++face)
					{
						info->numface += facelist.at(face).toElement().attribute("count").toInt() ;
						QDomNode no;
						no = findNodeBySpecificAttributeValue(facelist.at(face),"input","semantic","NORMAL");
						if (!no.isNull()) 
							bHasPerWedgeNormal = true;
						no = findNodeBySpecificAttributeValue(facelist.at(face),"input","semantic","TEXCOORD");
						if (!no.isNull()) 
							bHasPerWedgeTexCoord = true;
					}
				}
			}

			info->mask = 0;
		
			if (bHasPerWedgeTexCoord) 
				info->mask |= vcg::tri::io::Mask::IOM_WEDGTEXCOORD;
			if (bHasPerWedgeNormal) 
				info->mask |= vcg::tri::io::Mask::IOM_WEDGNORMAL;
			if (bHasPerVertexColor)	
				info->mask |= vcg::tri::io::Mask::IOM_VERTCOLOR;
			if (bHasPerFaceColor) 
				info->mask |= vcg::tri::io::Mask::IOM_FACECOLOR;
			if (bHasPerVertexNormal) 
				info->mask |= vcg::tri::io::Mask::IOM_VERTNORMAL;
			if (bHasPerVertexText) 
				info->mask |= vcg::tri::io::Mask::IOM_VERTTEXCOORD;
			
			

			delete (info->dae->doc);
			info->dae->doc = NULL;
			addinfo = info;
			return true;
		}
	};
}
}
}

#endif
