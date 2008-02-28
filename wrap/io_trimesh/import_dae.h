#ifndef __VCGLIB_IMPORTERDAE
#define __VCGLIB_IMPORTERDAE

//importer for collada's files

#include<wrap/dae/util_dae.h>

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


		static int LoadMesh(OpenMeshType& m,AdditionalInfoDAE* info,const QDomNode& geo,const vcg::Matrix44f& t, CallBackPos *cb=0)
		{
			if (isThereTag(geo,"mesh"))
			{
				if ((cb !=NULL) && (((info->numvert + info->numface)%100)==0) && !(*cb)((100*(info->numvert + info->numface))/(info->numvert + info->numface), "Vertex Loading"))
					return E_CANTOPEN;
				/*QDomNodeList geosrc = geo.toElement().elementsByTagName("source");
				int geosrc_size = geosrc.size();
				if (geosrc_size < 1)
				return E_NOVERTEXPOSITION;*/

				QDomNodeList vertices = geo.toElement().elementsByTagName("vertices");
				int vertices_size = vertices.size();
				if (vertices_size != 1)
					return E_INCOMPATIBLECOLLADA141FORMAT;

				QDomNode srcnode = attributeSourcePerSimplex(vertices.at(0),*(info->dae->doc),"POSITION");
				if (srcnode.isNull())
					return E_NOVERTEXPOSITION;

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
						
						assert((ii * 3 < geosrcposarr_size) && (ii * 3 + 1 < geosrcposarr_size) && (ii * 3 + 2 < geosrcposarr_size));
						vcg::Point4f tmp = t * vcg::Point4f(geosrcposarr[ii * 3].toFloat(),geosrcposarr[ii * 3 + 1].toFloat(),geosrcposarr[ii * 3 + 2].toFloat(),1.0f);
						m.vert[vv].P() = vcg::Point3f(tmp.X(),tmp.Y(),tmp.Z());

						if (!srcnodenorm.isNull())
						{
							assert((ii * 3 < geosrcvertnorm.size()) && (ii * 3 + 1 < geosrcvertnorm.size()) && (ii * 3 + 2 < geosrcvertnorm.size()));
							vcg::Matrix44f intr44 = vcg::Inverse(t);
							vcg::Transpose(intr44);
							Matrix33f intr33;
							for(unsigned int rr = 0; rr < 2; ++rr)
							{
								for(unsigned int cc = 0;cc < 2;++cc)
									intr33[rr][cc] = intr44[rr][cc];
							}
							m.vert[vv].N() = (intr33 * vcg::Point3f(geosrcvertnorm[ii * 3].toFloat(),geosrcvertnorm[ii * 3 + 1].toFloat(),geosrcvertnorm[ii * 3 + 2].toFloat())).Normalize();
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

	public:

		//merge all meshes in the collada's file in the templeted mesh m
		//I assume the mesh 
		
		static int Open(OpenMeshType& m,const char* filename,AdditionalInfo*& info, CallBackPos *cb=0)
		{
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

			int problem = E_NOERROR;
			bool found_a_mesh = false;
			//Is there geometry in the file? 
			bool geoinst_found = false;
			//for each scene in COLLADA FILE
			for(int scn = 0;scn < scn_size;++scn)
			{
				QDomNodeList instscenes = scenes.at(scn).toElement().elementsByTagName("instance_visual_scene");
				int instscn_size = instscenes.size();
				if (instscn_size == 0) 
					return E_INCOMPATIBLECOLLADA141FORMAT;

				//for each scene instance in a COLLADA scene
				for(int instscn = 0;instscn < instscn_size; ++instscn)
				{
					QString libscn_url;
					referenceToANodeAttribute(instscenes.at(instscn),"url",libscn_url);	
					QDomNode nd = QDomNode(*(inf->dae->doc));
					QDomNode visscn = findNodeBySpecificAttributeValue(*(inf->dae->doc),"visual_scene","id",libscn_url);
					if(visscn.isNull())
						return E_UNREFERENCEBLEDCOLLADAATTRIBUTE;
					
					//for each node in the libscn_url visual scene  
					QDomNodeList visscn_child = visscn.childNodes();
					
					//for each direct child of a libscn_url visual scene find if there is some geometry instance
					for(int chdind = 0; chdind < visscn_child.size();++chdind)
					{
						QDomNodeList geoinst = visscn_child.at(chdind).toElement().elementsByTagName("instance_geometry");
						int geoinst_size = geoinst.size();
						if (geoinst_size != 0)
						{
							
							geoinst_found |= true;
							QDomNodeList geolib = inf->dae->doc->elementsByTagName("library_geometries");
							assert(geolib.size() == 1);
							//!!!!!!!!!!!!!!!!!here will be the code for geometry transformations!!!!!!!!!!!!!!!!!!!!!!
							for(int geoinst_ind = 0;geoinst_ind < geoinst_size;++geoinst_ind)
							{
								QString geo_url;
								referenceToANodeAttribute(geoinst.at(geoinst_ind),"url",geo_url);
								
								QDomNode geo = findNodeBySpecificAttributeValue(geolib.at(0),"geometry","id",geo_url);
								if (geo.isNull())
									return E_UNREFERENCEBLEDCOLLADAATTRIBUTE;
								vcg::Matrix44f tr;
								tr.SetIdentity();
								TransfMatrix(visscn,geoinst.at(geoinst_ind),tr);
								problem |= LoadMesh(m,inf,geo,tr); 
								if (problem & E_NOMESH)
									found_a_mesh |= false;
								else
									found_a_mesh = true;
							}
						}
					}
					//if there is at least a mesh I clean the problem status variable from E_NOMESH ERROR
				
					if (((problem & E_NOMESH) || (problem & E_NOPOLYGONALMESH)) && (found_a_mesh))
					{	
						if (problem & E_NOMESH) 
							problem = problem & ~E_NOMESH;
						if (problem & E_NOPOLYGONALMESH)
							problem = problem & ~E_NOPOLYGONALMESH;
					}
				}
			}

			if (!geoinst_found)
			{
				QDomNodeList geolib = inf->dae->doc->elementsByTagName("library_geometries");
				assert(geolib.size() == 1);
				QDomNodeList geochild = geolib.at(0).childNodes();
				int geochild_size = geochild.size();
				int problem = 0;
				for(int chd = 0;chd < geochild_size;++chd)
				{
					vcg::Matrix44f tmp;
					tmp.SetIdentity();
					problem |= LoadMesh(m,inf,geochild.at(chd),tmp); 
				}
			}
			//if there is at least a mesh I clean the problem status variable from E_NOMESH or E_NOPOLYGONALMESH ERROR
			if (((problem & E_NOMESH) || (problem & E_NOPOLYGONALMESH)) && (found_a_mesh))
			{	
				if (problem & E_NOMESH) 
					problem = problem & ~E_NOMESH;
				if (problem & E_NOPOLYGONALMESH)
					problem = problem & ~E_NOPOLYGONALMESH;
			}
			info = inf;
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
				if (instscn_size == 0) 
					return false;

				//for each scene instance in a COLLADA scene
				for(int instscn = 0;instscn < instscn_size; ++instscn)
				{
					QString libscn_url;
					referenceToANodeAttribute(instscenes.at(instscn),"url",libscn_url);	
					QDomNode nd = QDomNode(*(info->dae->doc));
					QDomNode visscn = findNodeBySpecificAttributeValue(*(info->dae->doc),"visual_scene","id",libscn_url);
					if(visscn.isNull())
						return false;
					
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
