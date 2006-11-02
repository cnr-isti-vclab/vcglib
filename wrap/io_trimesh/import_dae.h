#ifndef __VCGLIB_IMPORTERDAE
#define __VCGLIB_IMPORTERDAE

//importer for collada's files

#include <wrap/io_trimesh/additionalinfo.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/allocate.h>
#include<map>

#include<QtXml/QDomDocument>
#include<QtCore/QFile>
#include <QtCore/QStringList>

#include<vcg/space/point3.h>
#include<vcg/space/tcoord2.h>
#include<vcg/space/color4.h>
//#include <hgrd/hgrd.h>

namespace vcg {
namespace tri {
namespace io {


	class InfoDAE : public AdditionalInfo
	{
		public:

		InfoDAE()
		{
			mask	= 0;
			numvert = 0;
			numface = 0;
			doc = NULL;
		}

		~InfoDAE()
		{
			delete doc;
			texturefile.clear();
		}

		QDomDocument* doc;		
		std::vector<std::string> texturefile; 
	};

	class AdditionalInfoDAE : public AdditionalInfo
	{
	public: 
		vcg::tri::io::InfoDAE* dae;

		AdditionalInfoDAE()
		:AdditionalInfo()
		{
		}

		~AdditionalInfoDAE()
		{
			delete dae;
		}
	};

	template<typename OpenMeshType>
	class ImporterDAE
	{
	public:

		//merge all meshes in the collada's file in the templeted mesh m
		//I assume the mesh 

		enum DAEError 
		{
			E_NOERROR,				// 0
			E_CANTOPEN,				// 1
			E_NOGEOMETRYLIBRARY,     // 2 
			E_NOMESH,      // 3
			E_NOVERTEXPOSITION,            // 4
			E_NO3DVERTEXPOSITION,			// 5
			E_NO3DSCENE, // 6
			E_INCOMPATIBLECOLLADA141FORMAT, //7
			E_UNREFERENCEBLEDCOLLADAATTRIBUTE, // 8
			E_NOTRIANGLES
		};

		static const char *ErrorMsg(int error)
		{
			static const char * dae_error_msg[] =
			{
				"No errors",
				"Can't open file",
				"File without a geometry library",
				"There isn't mesh in file",
				"The meshes in file haven't the vertex position attribute",
				"The importer assumes that the OpenMeshType uses a 3D point for the vertex position",
				"There isn't any scene in Collada file",
				"The input file is not compatible with COLLADA 1.41 standard format",
				"Collada file is trying to referece an attribute that is not in the file",
				"This version of Collada Importer support only triangular mesh file"
			};

			if(error>9 || error<0) return "Unknown error";
			else return dae_error_msg[error];
		};

		
		
	private:
		inline static void referenceToANodeAttribute(const QDomNode& n,const QString& attr,QString& url_st)
		{
			url_st = n.toElement().attribute(attr);
			int sz = url_st.size() - 1;
			url_st = url_st.right(sz);
			assert(url_st.size() != 0);
		}

		inline static QDomNode findNodeBySpecificAttributeValue(const QDomNode& n,const QString& tag,const QString& attrname,const QString& attrvalue)
		{
			QDomNode ndl = n.toElement();
			return findNodeBySpecificAttributeValue((QDomDocument&) ndl,tag,attrname,attrvalue);
		}

		inline static QDomNode findNodeBySpecificAttributeValue(const QDomDocument& n,const QString& tag,const QString& attrname,const QString& attrvalue)
		{
			QDomNodeList ndl = n.elementsByTagName(tag);
			int ndl_size = ndl.size();
			assert(ndl_size != 0);
			int ind = 0;
			while(ind < ndl_size)
			{
				if (ndl.at(ind).toElement().attribute(attrname) == attrvalue)
					return ndl.at(ind);
				++ind;
			}
			return QDomNode();
		}

		inline static bool isThereTag(const QDomNode& n,const QString& tagname)
		{
			QDomNode ndl = n.toElement();
			return isThereTag((QDomDocument&) n,tagname);
		}

		inline static bool isThereTag(const QDomDocument& n,const QString& tagname)
		{
			return ((n.toElement().elementsByTagName(tagname).size() > 0)? true : false);
		}


		inline static QDomNode attributeSourcePerSimplex(const QDomNode& n,const QDomDocument& startpoint,const QString& sem)
		{
			QDomNodeList vertattr = n.toElement().elementsByTagName("input");
			for(int ind = 0;ind < vertattr.size();++ind)
			{
				if (vertattr.at(ind).toElement().attribute("semantic") == sem)
				{
					QString url; 
					referenceToANodeAttribute(vertattr.at(ind),"source",url);
					return findNodeBySpecificAttributeValue(startpoint,"source","id",url);
				}
			}
			return QDomNode();
		}

		inline static void valueStringList(QStringList& res,const QDomNode& srcnode,const QString& tag) 
		{
			QDomNodeList list = srcnode.toElement().elementsByTagName(tag);
			int list_size = list.size();
			assert(list_size == 1);
			QString nd = list.at(0).firstChild().nodeValue();
			res = nd.split(" ");
			if (res.last() == "")
				res.removeLast();
		
		}

	public:
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
										m.vert[vv].P() = vcg::Point3f(geosrcposarr[ii * 3].toFloat(),geosrcposarr[ii * 3 + 1].toFloat(),geosrcposarr[ii * 3 + 2].toFloat());
										
										if (!srcnodenorm.isNull())
										{
											assert((ii * 3 < geosrcvertnorm.size()) && (ii * 3 + 1 < geosrcvertnorm.size()) && (ii * 3 + 2 < geosrcvertnorm.size()));
											m.vert[vv].N() = vcg::Point3f(geosrcvertnorm[ii * 3].toFloat(),geosrcvertnorm[ii * 3 + 1].toFloat(),geosrcvertnorm[ii * 3 + 2].toFloat());
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
								} 
							}
						}
					}
				}
			}

			if (!geoinst_found)
				return E_NOGEOMETRYLIBRARY;
			return E_NOERROR;
		}

		static bool LoadMask(const char * filename, AdditionalInfoDAE &addinfo)
		{
			std::ifstream stream(filename);
			if (stream.fail())
				return false;

			stream.seekg (0, std::ios::end);
			int length = stream.tellg();
			if (length == 0) return false;
			stream.seekg (0, std::ios::beg);

			bool bHasPerWedgeTexCoord = false;
			bool bHasPerWedgeNormal		= false;
			bool bUsingMaterial				= false;
			bool bHasPerVertexColor		= false;
			bool bHasPerFaceColor			= false;

			AdditionalInfoDAE* inf = new AdditionalInfoDAE();
			inf->dae = new InfoDAE(); 
			InfoDAE* info = inf->dae;
			info->doc = new FCDocument();

			unsigned int numvert = 0;
			unsigned int numtriang = 0;

			unsigned int mask = 0;

			FCDGeometryLibrary* geolib = info->doc->GetGeometryLibrary();
			if (geolib->IsEmpty()) return false;
			size_t n = geolib->GetEntityCount();
			std::vector<FCDGeometryMesh*> geomsh(n);


			FUStatus st = info->doc->LoadFromFile(FUStringConversion::ToFString(filename));
			if (st.IsFailure()) 
			{
				delete info->doc;
				info->doc = NULL;
				return false;
			}

			bool amesh = false; 

			//for any mesh in the collada file
			for(unsigned int ii = 0;ii < geomsh.size();++ii)
			{
				if (!geolib->GetEntity(ii)->IsMesh())
				{
					amesh |= false;
				}
				else
				{
					amesh |= true;
					geomsh[ii] = geolib->GetEntity(ii)->GetMesh();

					unsigned int ver;
					if (geomsh[ii]->GetFaceCount() > 0)
					{
						geomsh[ii]->Triangulate(); 

						size_t dim = geomsh[ii]->GetFaceVertexCount() / geomsh[ii]->GetFaceCount();
						assert(dim == 3);
						//MyMesh* msh = new MyMesh();
						//size_t  nattr = geomsh[ii]->GetSourceCount();
						//FCDGeometrySourceList& srclst = geomsh[ii]->GetVertexSources(); 

						FCDGeometrySource* src;
						if ((src = geomsh[ii]->GetPositionSource()) != NULL)
						{
							FloatList& flst = src->GetSourceData();
							unsigned int str = src->GetSourceStride();
							assert(flst.size() % str == 0);
							ver = flst.size() / str;
							numvert += flst.size() / str;  
						}
						else 
						{
							delete info->doc;
							info->doc = NULL;
							return false;
						}

						size_t pol = geomsh[ii]->GetPolygonsCount();

						for(unsigned int pset = 0; pset < pol;++pset)
						{
							FCDGeometryMesh* tmp = geomsh[ii];
							FCDGeometryPolygonsInput* pos = tmp->GetPolygons(pset)->FindInput(FUDaeGeometryInput::POSITION);
							if ((pos == NULL) || (pos->GetSource()->GetSourceStride() != 3)) 
							{
								delete info->doc;
								info->doc = NULL;
								return false;
							}
							//unsigned int hi = pos->indices[1];

							FCDGeometryPolygonsInputList normlist;
							tmp->GetPolygons(pset)->FindInputs(FUDaeGeometryInput::NORMAL,normlist);
							FCDGeometryPolygonsInputList tet; 
							tmp->GetPolygons(pset)->FindInputs(FUDaeGeometryInput::TEXCOORD,tet);

							for(unsigned int kk = 0; kk < tet.size();++kk)
								if ((normlist[0]->GetSource()->GetSourceData().size() == ver * 3 * 3) && (!bHasPerWedgeNormal))
									bHasPerWedgeNormal = true;

							for(unsigned int kk = 0; kk < tet.size();++kk)
								if ((tet[kk]->GetSource()->GetSourceData().size() == ver * 3 * 2) && (!bHasPerTexCoord))
									bHasPerTexCoord = true;
						}
					}
				}
			}
			inf.nvert = 
				addinfo = inf;
			return true;
		}
	};
}
}
}

#endif
