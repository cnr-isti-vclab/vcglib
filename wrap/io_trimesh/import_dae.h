#ifndef __VCGLIB_IMPORTERDAE
#define __VCGLIB_IMPORTERDAE

//importer for collada's files

#include<wrap/io_trimesh/util_dae.h>

namespace vcg {
namespace tri {
namespace io {

	template<typename OpenMeshType>
	class ImporterDAE : public UtilDAE
	{

	private:

		static int LoadMesh(OpenMeshType& m,InfoDAE* info,const QDomNode& geo,const vcg::Matrix44f& t)
		{
			if (isThereTag(geo,"mesh"))
			{
				/*QDomNodeList geosrc = geo.toElement().elementsByTagName("source");
				int geosrc_size = geosrc.size();
				if (geosrc_size < 1)
				return E_NOVERTEXPOSITION;*/

				QDomNodeList vertices = geo.toElement().elementsByTagName("vertices");
				int vertices_size = vertices.size();
				if (vertices_size != 1)
					return E_INCOMPATIBLECOLLADA141FORMAT;

				QDomNode srcnode = attributeSourcePerSimplex(vertices.at(0),*(info->doc),"POSITION");
				if (srcnode.isNull())
					return E_NOVERTEXPOSITION;

				QStringList geosrcposarr;
				valueStringList(geosrcposarr,srcnode,"float_array");

				int geosrcposarr_size = geosrcposarr.size();
				if ((geosrcposarr_size % 3) != 0)
					return E_CANTOPEN;
				int nvert = geosrcposarr_size / 3;
				size_t offset = m.vert.size();
				vcg::tri::Allocator<OpenMeshType>::AddVertices(m,nvert);

				QDomNode srcnodenorm = attributeSourcePerSimplex(vertices.at(0),*(info->doc),"NORMAL");
				QStringList geosrcvertnorm;
				if (!srcnodenorm.isNull())
					valueStringList(geosrcvertnorm,srcnodenorm,"float_array");

				QDomNode srcnodetext = attributeSourcePerSimplex(vertices.at(0),*(info->doc),"TEXCOORD");
				QStringList geosrcverttext;
				if (!srcnodetext.isNull())
					valueStringList(geosrcverttext,srcnodetext,"float_array");

				QDomNode srcnodecolor = attributeSourcePerSimplex(vertices.at(0),*(info->doc),"COLOR");
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
						vcg::Matrix44f intr44 = vcg::Transpose(vcg::Inverse(t));
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
						m.vert[vv].T() = vcg::TCoord2<float>();
						m.vert[vv].T().u() = geosrcverttext[ii * 2].toFloat();
						m.vert[vv].T().v() = geosrcverttext[ii * 2 + 1].toFloat();
					}
					++ii;
				}

				QDomNodeList tripatch = geo.toElement().elementsByTagName("triangles");
				int tripatch_size = tripatch.size();
				if (tripatch_size == 0)
					return E_NOTRIANGLES;

				for(int tript = 0; tript < tripatch_size;++tript)
				{

					int nfcatt = tripatch.at(tript).toElement().elementsByTagName("input").size();

					QStringList face;
					valueStringList(face,tripatch.at(tript),"p");
					int face_size = face.size();
					int offsetface = (int)m.face.size();
					if (face_size == 0) return E_NOMESH;
					vcg::tri::Allocator<OpenMeshType>::AddFaces(m,face_size / (nfcatt * 3));
					QDomNode wnsrc = QDomNode();
					QStringList wn;
					wnsrc = findNodeBySpecificAttributeValue(tripatch.at(tript),"input","semantic","NORMAL");
					int offnm;
					if (!wnsrc.isNull())
					{
						offnm = wnsrc.toElement().attribute("offset").toInt();
						QDomNode sn = attributeSourcePerSimplex(tripatch.at(tript),*(info->doc),"NORMAL");
						valueStringList(wn,sn,"float_array");
					}

					QDomNode wtsrc = QDomNode();
					QStringList wt;
					wtsrc = findNodeBySpecificAttributeValue(tripatch.at(tript),"input","semantic","TEXCOORD");
					int offtx;
					if (!wtsrc.isNull())
					{
						offtx = wtsrc.toElement().attribute("offset").toInt();
						QDomNode st = attributeSourcePerSimplex(tripatch.at(tript),*(info->doc),"TEXCOORD");
						valueStringList(wt,st,"float_array");
					}

					QDomNode wcsrc = QDomNode();
					QStringList wc;
					wcsrc = findNodeBySpecificAttributeValue(tripatch.at(tript),"input","semantic","COLOR");
					int offcl;
					if (!wcsrc.isNull())
					{
						offcl = wcsrc.toElement().attribute("offset").toInt();
						QDomNode sc = attributeSourcePerSimplex(tripatch.at(tript),*(info->doc),"COLOR");
						valueStringList(wc,sc,"float_array");
					}

					int jj = 0;	
					//int dd = m.face.size();
					for(int ff = offsetface;ff < (int) m.face.size();++ff)
					{ 
						int indvt = face.at(jj).toInt();
						assert(indvt + offset < m.vert.size());
						m.face[ff].V(0) = &(m.vert[indvt + offset]);

						int indnm;
						if (!wnsrc.isNull())
						{
							indnm = face.at(jj + offnm).toInt();
							assert(indnm * 3 < wn.size());
							m.face[ff].WN(0) = vcg::Point3f(wn.at(indnm * 3).toFloat(),wn.at(indnm * 3 + 1).toFloat(),wn.at(indnm * 3 + 2).toFloat());
						}

						int indtx;
						if (!wtsrc.isNull())
						{
							indtx = face.at(jj + offtx).toInt();
							assert(indtx * 2 < wt.size());
							m.face[ff].WT(0) = vcg::TCoord2<float>();
							m.face[ff].WT(0).u() = wt.at(indtx * 2).toFloat();
							m.face[ff].WT(0).v() = wt.at(indtx * 2 + 1).toFloat();	
						}

						/*int indcl;
						if (!wcsrc.isNull())
						{
						indcl = face.at(jj + offcl).toInt();
						assert(indcl * 4 < wc.size());
						m.face[ff].WC(0) = vcg::Color4b(wc.at(indcl * 4).toFloat(),wc.at(indcl * 4 + 1).toFloat(),wc.at(indcl * 4 + 2).toFloat(),wc.at(indcl * 4 + 3).toFloat());
						}*/
						jj += nfcatt;

						indvt = face.at(jj).toInt();
						assert(indvt + offset < m.vert.size());
						m.face[ff].V(1) = &(m.vert[indvt + offset]);
						if (!wnsrc.isNull())
						{
							indnm = face.at(jj + offnm).toInt();
							assert(indnm * 3 < wn.size());
							m.face[ff].WN(1) = vcg::Point3f(wn.at(indnm * 3).toFloat(),wn.at(indnm * 3 + 1).toFloat(),wn.at(indnm * 3 + 2).toFloat());
						}

						if (!wtsrc.isNull())
						{
							indtx = face.at(jj + offtx).toInt();
							assert(indtx * 2 < wt.size());
							m.face[ff].WT(1) = vcg::TCoord2<float>();
							m.face[ff].WT(1).u() = wt.at(indtx * 2).toFloat();
							m.face[ff].WT(1).v() = wt.at(indtx * 2 + 1).toFloat();	
						}

						/*if (!wcsrc.isNull())
						{
						indcl = face.at(jj + offcl).toInt();
						assert(indcl * 4 < wc.size());
						m.face[ff].WC(1) = vcg::Color4b(wc.at(indcl * 4).toFloat(),wc.at(indcl * 4 + 1).toFloat(),wc.at(indcl * 4 + 2).toFloat(),wc.at(indcl * 4 + 3).toFloat());
						}*/
						jj += nfcatt;

						indvt = face.at(jj).toInt();
						assert(indvt + offset < m.vert.size());
						m.face[ff].V(2) = &(m.vert[indvt + offset]);
						if (!wnsrc.isNull())
						{
							indnm = face.at(jj + offnm).toInt();
							assert(indnm * 3 < wn.size());
							m.face[ff].WN(2) = vcg::Point3f(wn.at(indnm * 3).toFloat(),wn.at(indnm * 3 + 1).toFloat(),wn.at(indnm * 3 + 2).toFloat());
						}

						if (!wtsrc.isNull())
						{
							indtx = face.at(jj + offtx).toInt();
							assert(indtx * 2 < wt.size());
							m.face[ff].WT(2) = vcg::TCoord2<float>();
							m.face[ff].WT(2).u() = wt.at(indtx * 2).toFloat();
							m.face[ff].WT(2).v() = wt.at(indtx * 2 + 1).toFloat();	
						}

						/*if (!wcsrc.isNull())
						{
						indcl = face.at(jj + offcl).toInt();
						assert(indcl * 4 < wc.size());
						m.face[ff].WC(2) = vcg::Color4b(wc.at(indcl * 4).toFloat(),wc.at(indcl * 4 + 1).toFloat(),wc.at(indcl * 4 + 2).toFloat(),wc.at(indcl * 4 + 3).toFloat());
						}*/
						jj += nfcatt;

					}
				}
				return E_NOERROR;
			}
			else return E_NOMESH;
		}

		static void GetTexture(const QDomDocument& doc,AdditionalInfoDAE* inf)
		{
			QDomNodeList txlst = doc.elementsByTagName("library_images");
			for(int img = 0;img < txlst.size();++img)
			{
				QDomNodeList nlst = txlst.at(img).toElement().elementsByTagName("init_from");
				if (nlst.size() > 0)
				{
					inf->dae->texturefile.push_back(nlst.at(0).firstChild().nodeValue());
				}
			}
		}
	public:

		//merge all meshes in the collada's file in the templeted mesh m
		//I assume the mesh 
		
		static int Open(OpenMeshType& m,const char* filename,AdditionalInfo*& addinfo)
		{
			AdditionalInfoDAE* inf = new AdditionalInfoDAE();
			inf->dae = new InfoDAE(); 
			InfoDAE* info = inf->dae;

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
			
			info->doc = doc;
			//GetTexture(*(info->doc),inf);

			QDomNodeList& scenes = info->doc->elementsByTagName("scene");
			int scn_size = scenes.size();
			if (scn_size == 0) 
				return E_NO3DSCENE;

			//Is there geometry in the file? 
			bool geoinst_found = false;
			//for each scene in COLLADA FILE
			for(int scn = 0;scn < scn_size;++scn)
			{
				QDomNodeList& instscenes = scenes.at(scn).toElement().elementsByTagName("instance_visual_scene");
				int instscn_size = instscenes.size();
				if (instscn_size == 0) 
					return E_INCOMPATIBLECOLLADA141FORMAT;

				//for each scene instance in a COLLADA scene
				for(int instscn = 0;instscn < instscn_size; ++instscn)
				{
					QString libscn_url;
					referenceToANodeAttribute(instscenes.at(instscn),"url",libscn_url);	
					QDomNode nd = QDomNode(*(info->doc));
					QDomNode visscn = findNodeBySpecificAttributeValue(*(info->doc),"visual_scene","id",libscn_url);
					if(visscn.isNull())
						return E_UNREFERENCEBLEDCOLLADAATTRIBUTE;
					
					//for each node in the libscn_url visual scene  
					QDomNodeList& visscn_child = visscn.childNodes();
					
					//for each direct child of a libscn_url visual scene find if there is some geometry instance
					int problem = 0;
					for(int chdind = 0; chdind < visscn_child.size();++chdind)
					{
						QDomNodeList& geoinst = visscn_child.at(chdind).toElement().elementsByTagName("instance_geometry");
						int geoinst_size = geoinst.size();
						if (geoinst_size != 0)
						{
							
							geoinst_found |= true;
							QDomNodeList& geolib = info->doc->elementsByTagName("library_geometries");
							int geolib_size = geolib.size();
							assert(geolib_size == 1);
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
								problem |= LoadMesh(m,info,geo,tr); 
								if (problem) return problem;
							}
						}
					}
				}
			}

			if (!geoinst_found)
			{
				QDomNodeList& geolib = info->doc->elementsByTagName("library_geometries");
				int geolib_size = geolib.size();
				assert(geolib_size == 1);
				QDomNodeList& geochild = geolib.at(0).childNodes();
				int geochild_size = geochild.size();
				int problem = 0;
				for(int chd = 0;chd < geochild_size;++chd)
				{
					vcg::Matrix44f tmp;
					tmp.SetIdentity();
					problem |= LoadMesh(m,info,geochild.at(chd),tmp); 
					if (problem) return problem;
				}
			}
			addinfo = inf;
			return E_NOERROR;
		}

		static bool LoadMask(const char * filename, AdditionalInfoDAE*& addinfo)
		{
			bool bHasPerWedgeTexCoord = false;
			bool bHasPerWedgeNormal		= false;
			bool bHasPerVertexColor		= false;
			bool bHasPerFaceColor			= false;
			bool bHasPerVertexNormal = false;
			bool bHasPerVertexText = false;
			
			AdditionalInfoDAE* inf = new AdditionalInfoDAE();
			inf->dae = new InfoDAE(); 
			InfoDAE* info = inf->dae;

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
			

			info->doc = doc;
			GetTexture(*(info->doc),inf);
			QDomNodeList& scenes = info->doc->elementsByTagName("scene");
			int scn_size = scenes.size();
			

			//Is there geometry in the file? 
			bool geoinst_found = false;
			//for each scene in COLLADA FILE
			for(int scn = 0;scn < scn_size;++scn)
			{
				QDomNodeList& instscenes = scenes.at(scn).toElement().elementsByTagName("instance_visual_scene");
				int instscn_size = instscenes.size();
				if (instscn_size == 0) 
					return false;

				//for each scene instance in a COLLADA scene
				for(int instscn = 0;instscn < instscn_size; ++instscn)
				{
					QString libscn_url;
					referenceToANodeAttribute(instscenes.at(instscn),"url",libscn_url);	
					QDomNode nd = QDomNode(*(info->doc));
					QDomNode visscn = findNodeBySpecificAttributeValue(*(info->doc),"visual_scene","id",libscn_url);
					if(visscn.isNull())
						return false;
					
					//for each node in the libscn_url visual scene  
					QDomNodeList& visscn_child = visscn.childNodes();
					
					//for each direct child of a libscn_url visual scene find if there is some geometry instance
					int problem = 0;
					for(int chdind = 0; chdind < visscn_child.size();++chdind)
					{
						QDomNodeList& geoinst = visscn_child.at(chdind).toElement().elementsByTagName("instance_geometry");
						int geoinst_size = geoinst.size();
						if (geoinst_size != 0)
						{
							
							geoinst_found |= true;
							QDomNodeList& geolib = info->doc->elementsByTagName("library_geometries");
							int geolib_size = geolib.size();
							assert(geolib_size == 1);
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

								QDomNodeList facelist = geo.toElement().elementsByTagName("triangles");
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
			
			if (!geoinst_found)
			{
				QDomNodeList& geolib = info->doc->elementsByTagName("library_geometries");
				int geolib_size = geolib.size();
				assert(geolib_size == 1);
				QDomNodeList& geochild = geolib.at(0).toElement().elementsByTagName("geometry");
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
			
			

			delete (info->doc);
			addinfo = inf;
			return true;
		}
	};
}
}
}

#endif
