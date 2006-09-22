#ifndef __VCGLIB_IMPORTERDAE
#define __VCGLIB_IMPORTERDAE

//importer for collada's files

#include <FCollada.h>
#include <FUtils/FUStringConversion.h>
#include <FCDocument/FCDocument.h>
#include <FCDocument/FCDLibrary.h>
#include <FCDocument/FCDGeometry.h>
#include <FCDocument/FCDGeometryMesh.h>
#include <FCDocument/FCDGeometrySource.h>
#include <FCDocument/FCDGeometryPolygons.h>
#include <FCDocument/FCDImage.h>


//#include <wrap/gl/trimesh.h>
#include <wrap/additionalinfo.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/allocate.h>

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

		FCDocument* doc;		
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
		E_NO3DVERTEXPOSITION			// 5
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
			"The importer assumes that the OpenMeshType uses a 3D point for the vertex position"
		};

		if(error>5 || error<0) return "Unknown error";
		else return dae_error_msg[error];
	};

	static int Open(OpenMeshType& m,const char* filename)
	{
		assert(filename!=0);
		FCDocument* doc = new FCDocument();
		
		FUStatus st = doc->LoadFromFile(FUStringConversion::ToFString(filename));
		if (st.IsFailure())
		{
			delete doc;
			doc = NULL;
			return E_CANTOPEN;
		}
		FCDGeometryLibrary* geolib = doc->GetGeometryLibrary();
		if (geolib->IsEmpty()) 
		{
			delete doc;
			return E_NOGEOMETRYLIBRARY;
		}
		size_t n = geolib->GetEntityCount();
		std::vector<FCDGeometryMesh*> geomsh(n);
		
		//for any mesh in the collada file

 		for(unsigned int ii = 0;ii < geomsh.size();++ii)
		{
			if (!geolib->GetEntity(ii)->IsMesh())
			{
				delete doc;
				return E_NOMESH;
			}
			else
			{
				geomsh[ii] = geolib->GetEntity(ii)->GetMesh();
				unsigned int offset = m.vert.size();
				if (geomsh[ii]->GetFaceCount() > 0)
				{
					geomsh[ii]->Triangulate(); 
					/*std::vector< std::vector<OpenMeshType::VertexType> > vt(m.face.size());
					for(std::vector< std::vector<OpenMeshType::VertexType> >::iterator
					HGRD<OpenMeshType::VertexType>::Triangulate(*/
					
					//geomsh[ii]->Get
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
						for(unsigned int cont = 0;cont < flst.size();cont += str)
						{
							typename OpenMeshType::VertexIterator vi=vcg::tri::Allocator<OpenMeshType>::AddVertices(m,1);
							vi->P()= vcg::Point3f(flst[cont],flst[cont + 1],flst[cont + 2]);
							vi->N() = vcg::Point3f(0.0,0.0,0.0);
						}
					}

					else 
					{
						delete doc;
						return E_NOVERTEXPOSITION;
					}

					//a single mesh may be composed by a variable numbers of polygons' subsets
					size_t pol = geomsh[ii]->GetPolygonsCount();

					//for any polygons' subset in a single mesh
					for(unsigned int pset = 0; pset < pol;++pset)
					{
						FCDGeometryMesh* tmp = geomsh[ii];
						FCDGeometryPolygonsInput* pos = tmp->GetPolygons(pset)->FindInput(FUDaeGeometryInput::POSITION);
						if ((pos == NULL) || (pos->GetSource()->GetSourceStride() != 3)) 
						{
							delete doc;
							return E_NO3DVERTEXPOSITION;
						}
						//unsigned int hi = pos->indices[1];
						FCDGeometryPolygonsInput* norm = tmp->GetPolygons(pset)->FindInput(FUDaeGeometryInput::NORMAL);
						//unsigned int li = norm->indices[1];
						FCDGeometryPolygonsInput* text = tmp->GetPolygons(pset)->FindInput(FUDaeGeometryInput::TEXCOORD);
						
						bool isvalidwnorm = (m.HasPerWedgeNormal()) && (norm != NULL) && (norm->GetSource()->GetSourceStride() == 3); 
						bool isvalidnorm = (m.HasPerVertexNormal()) && (norm != NULL) && (norm->GetSource()->GetSourceStride() == 3);
						bool isvalidtext = (HasPerWedgeTexture(m)) && (text != NULL) && (text->GetSource()->GetSourceStride() == 2);
						
						FCDGeometryPolygonsInputList normlist;
						tmp->GetPolygons(pset)->FindInputs(FUDaeGeometryInput::NORMAL,normlist);
						FCDGeometryPolygonsInputList tet; 
						tmp->GetPolygons(pset)->FindInputs(FUDaeGeometryInput::TEXCOORD,tet);
						

						for(unsigned int ind = 0;ind < pos->indices.size();++ind)
						{
							typename OpenMeshType::FaceIterator fi=vcg::tri::Allocator<OpenMeshType>::AddFaces(m,1);		
							assert(pos->indices[ind] < m.vert.size());
							fi->V(0) = &m.vert[offset + pos->indices[ind]];
							
							size_t dimn = norm->indices.size();
							if (isvalidnorm) 
							{
								//assert(norm->indices[ind]  * 3  < norm->source->GetSourceData().size());
								UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(normlist[0]);
								//fi->V(0)->N() += vcg::Point3f(norm->source->GetSourceData()[norm->indices[ind] * 3],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								fi->V(0)->N() += vcg::Point3f(normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 2]);
								//++fi->V(0)->incidentfaces;
							}

							if (isvalidtext)
							{
								for(unsigned int hh = 0; hh < tet.size();++hh)
								{
									//NON CAMBIARE!!!!E' L'unico modo in cui restituisce gli indici corretti quando c'e' piu' di un insieme con la stessa semantica!!
									UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(tet[hh]);
									fi->WT(0).t(hh) = vcg::Point2f(tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2],tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1]);
								}
							}

							if (isvalidwnorm)
							{
								fi->WN(0) = vcg::Point3f(norm->GetSource()->GetSourceData()[norm->indices[ind] * 3],norm->GetSource()->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->GetSource()->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								fi->WN(1) = vcg::Point3f(norm->GetSource()->GetSourceData()[norm->indices[ind + 1] * 3],norm->GetSource()->GetSourceData()[norm->indices[ind + 1]  * 3 + 1],norm->GetSource()->GetSourceData()[norm->indices[ind + 1]  * 3 + 2]).Normalize();
								fi->WN(2) = vcg::Point3f(norm->GetSource()->GetSourceData()[norm->indices[ind + 2] * 3],norm->GetSource()->GetSourceData()[norm->indices[ind + 2]  * 3 + 1],norm->GetSource()->GetSourceData()[norm->indices[ind + 2]  * 3 + 2]).Normalize();
							}

							++ind;
							assert(pos->indices[ind] < m.vert.size());
							fi->V(1) = &m.vert[offset + pos->indices[ind]];
							
							if (isvalidnorm) 
							{	
								//assert(norm->indices[ind]  * 3 < norm->source->GetSourceData().size());
								//fi->V(1)->N() += vcg::Point3f(norm->source->GetSourceData()[norm->indices[ind] * 3],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(normlist[0]);
								//fi->V(0)->N() += vcg::Point3f(norm->source->GetSourceData()[norm->indices[ind] * 3],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								fi->V(1)->N() += vcg::Point3f(normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 2]);
								
								//++fi->V(1)->incidentfaces;
							}

							if (isvalidtext)
							{
								for(unsigned int hh = 0; hh < tet.size();++hh)
								{
									UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(tet[hh]);
									fi->WT(1).t(hh) = vcg::Point2f(tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2],tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1]);
								}
							}

							++ind;
							assert(pos->indices[ind] < m.vert.size());
							fi->V(2) = &m.vert[offset + pos->indices[ind]];
							
							if (isvalidnorm)
							{
								//assert(norm->indices[ind]  * 3 < norm->source->GetSourceData().size());
								//fi->V(2)->N() += vcg::Point3f(norm->source->GetSourceData()[norm->indices[ind]  * 3],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(normlist[0]);
								//fi->V(0)->N() += vcg::Point3f(norm->source->GetSourceData()[norm->indices[ind] * 3],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								fi->V(2)->N() += vcg::Point3f(normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 2]);
								//++fi->V(2)->incidentfaces;
							}
							
							if (isvalidtext)
							{
								for(unsigned int hh = 0; hh < tet.size();++hh)
								{
									UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(tet[hh]);
									fi->WT(2).t(hh) = vcg::Point2f(tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2],tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1]);
								}
							}

							if (isvalidnorm) fi->N() = ((fi->V(1)->P() - fi->V(0)->P()) ^ (fi->V(2)->P() - fi->V(0)->P())).Normalize();

							/*FCDGeometryPolygonsInput* posa = tmp->GetPolygons(pset)->FindInput(FUDaeGeometryInput::POSITION);
							FloatList& list = posa->source->GetSourceData();
							int dim = list.size();
							list[0] = -100.0;*/
						}
						//vm.push_back(msh);
						
						if (isvalidnorm)
						{
							vcg::tri::UpdateNormals<OpenMeshType>::PerVertexNormalized(m);
							/*for(MyMesh::VertexIterator vit = msh->vert.begin(); vit != msh->vert.end();++vit)
								vit->N() = (vit->N() / vit->incidentfaces).Normalize();*/
						}
					}
				}
			}
		}

		//doc->WriteToFile("PincoPalla.dae");
		delete doc;
		return E_NOERROR;
	}

	/*this open function should be used when you want to maintain the Collada's XML tree. If the file will correctly opened 
	doc argument in the function's signiture will contain the pointer to XML tree otherwise a NULL pointer*/


	static int Open(OpenMeshType& m,const char* filename,AdditionalInfo*& addinfo)
	{
		AdditionalInfoDAE* inf = new AdditionalInfoDAE();
		inf->dae = new InfoDAE(); 
		InfoDAE* info = inf->dae;
		info->doc = new FCDocument();

		FUStatus st = info->doc->LoadFromFile(FUStringConversion::ToFString(filename));
		if (st.IsFailure()) 
		{
			delete info->doc;
			info->doc = NULL;
			return E_CANTOPEN;
		}

		FCDGeometryLibrary* geolib = info->doc->GetGeometryLibrary();
		if (geolib->IsEmpty()) return E_NOGEOMETRYLIBRARY;
		size_t n = geolib->GetEntityCount();
		std::vector<FCDGeometryMesh*> geomsh(n);
		

		//it tests if there is at least a mesh in the collada file
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
				unsigned int offset = m.vert.size();
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
						for(unsigned int cont = 0;cont < flst.size();cont += str)
						{
							OpenMeshType::VertexIterator vi=vcg::tri::Allocator<OpenMeshType>::AddVertices(m,1);
							vi->P()= vcg::Point3f(flst[cont],flst[cont + 1],flst[cont + 2]);
							vi->N() = vcg::Point3f(0.0,0.0,0.0);
						}
					}
					else 
					{
						delete info->doc;
						info->doc = NULL;
						return E_NOVERTEXPOSITION;
					}

					//a single mesh may be composed by a variable numbers of polygons' subsets
					size_t pol = geomsh[ii]->GetPolygonsCount();

					//unsigned int offset = m.vert.size();

					//for any polygons' subset in a single mesh
					for(unsigned int pset = 0; pset < pol;++pset)
					{
						FCDGeometryMesh* tmp = geomsh[ii];
						FCDGeometryPolygonsInput* pos = tmp->GetPolygons(pset)->FindInput(FUDaeGeometryInput::POSITION);
						if ((pos == NULL) || (pos->GetSource()->GetSourceStride() != 3)) 
						{
							delete info->doc;
							info->doc = NULL;
							return E_NO3DVERTEXPOSITION;
						}
						//unsigned int hi = pos->indices[1];
						FCDGeometryPolygonsInput* norm = tmp->GetPolygons(pset)->FindInput(FUDaeGeometryInput::NORMAL);
						//unsigned int li = norm->indices[1];
						FCDGeometryPolygonsInput* text = tmp->GetPolygons(pset)->FindInput(FUDaeGeometryInput::TEXCOORD);
						
						bool isvalidwnorm = (m.HasPerWedgeNormal()) && (norm != NULL) && (norm->GetSource()->GetSourceStride() == 3); 
						bool isvalidnorm = (m.HasPerVertexNormal()) && (norm != NULL) && (norm->GetSource()->GetSourceStride() == 3);
						bool isvalidtext = (HasPerWedgeTexture(m)) && (text != NULL) && (text->GetSource()->GetSourceStride() == 2);
						
						//m.Enable(

						if (isvalidtext)
						{
							FCDImageLibrary* imlib = NULL;
							imlib = info->doc->GetImageLibrary();
							
							if (imlib != NULL)
							{
								info->texturefile.resize(imlib->GetEntityCount());
								for(unsigned int gg = 0;gg < imlib->GetEntityCount();++gg)
								{
									 FCDImage* img = imlib->GetEntity(gg);
									 info->texturefile[gg] = FUStringConversion::ToString(img->GetFilename());
								}
							}
						}

						FCDGeometryPolygonsInputList normlist;
						tmp->GetPolygons(pset)->FindInputs(FUDaeGeometryInput::NORMAL,normlist);
						FCDGeometryPolygonsInputList tet; 
						tmp->GetPolygons(pset)->FindInputs(FUDaeGeometryInput::TEXCOORD,tet);

						for(unsigned int ind = 0;ind < pos->indices.size();++ind)
						{
							OpenMeshType::FaceIterator fi=vcg::tri::Allocator<OpenMeshType>::AddFaces(m,1);		
							assert(offset + pos->indices[ind]  < m.vert.size());
							fi->V(0) = &m.vert[offset + pos->indices[ind]];
							
							//size_t dimn = norm->indices.size();
							if (isvalidnorm) 
							{
								//assert(norm->indices[ind]  * 3  < norm->source->GetSourceData().size());
								UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(normlist[0]);
								//fi->V(0)->N() += vcg::Point3f(norm->source->GetSourceData()[norm->indices[ind] * 3],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								fi->V(0)->N() += vcg::Point3f(normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 2]);
								//++fi->V(0)->incidentfaces;
							}

							if (isvalidtext)
							{
								for(unsigned int hh = 0; hh < tet.size();++hh)
								{
									//NON CAMBIARE!!!!E' L'unico modo in cui restituisce gli indici corretti quando c'e' piu' di un insieme con la stessa semantica!!
									UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(tet[hh]);
									fi->WT(0).t(hh) = vcg::Point2f(tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2],tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1]);
								}
							}

							if (isvalidwnorm)
							{
								fi->WN(0) = vcg::Point3f(norm->GetSource()->GetSourceData()[norm->indices[ind] * 3],norm->GetSource()->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->GetSource()->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								fi->WN(1) = vcg::Point3f(norm->GetSource()->GetSourceData()[norm->indices[ind + 1] * 3],norm->GetSource()->GetSourceData()[norm->indices[ind + 1]  * 3 + 1],norm->GetSource()->GetSourceData()[norm->indices[ind + 1]  * 3 + 2]).Normalize();
								fi->WN(2) = vcg::Point3f(norm->GetSource()->GetSourceData()[norm->indices[ind + 2] * 3],norm->GetSource()->GetSourceData()[norm->indices[ind + 2]  * 3 + 1],norm->GetSource()->GetSourceData()[norm->indices[ind + 2]  * 3 + 2]).Normalize();
							}

							++ind;
							assert(offset + pos->indices[ind] < m.vert.size());
							fi->V(1) = &m.vert[offset + pos->indices[ind]];
							
							if (isvalidnorm) 
							{	
								//assert(norm->indices[ind]  * 3 < norm->source->GetSourceData().size());
								//fi->V(1)->N() += vcg::Point3f(norm->source->GetSourceData()[norm->indices[ind] * 3],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(normlist[0]);
								//fi->V(0)->N() += vcg::Point3f(norm->source->GetSourceData()[norm->indices[ind] * 3],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								fi->V(1)->N() += vcg::Point3f(normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 2]);
								
								//++fi->V(1)->incidentfaces;
							}

							if (isvalidtext)
							{
								for(unsigned int hh = 0; hh < tet.size();++hh)
								{
									UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(tet[hh]);
									fi->WT(1).t(hh) = vcg::Point2f(tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2],tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1]);
								}
							}

							++ind;
							assert(offset + pos->indices[ind] < m.vert.size());
							fi->V(2) = &m.vert[offset + pos->indices[ind]];
							
							if (isvalidnorm)
							{
								//assert(norm->indices[ind]  * 3 < norm->source->GetSourceData().size());
								//fi->V(2)->N() += vcg::Point3f(norm->source->GetSourceData()[norm->indices[ind]  * 3],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(normlist[0]);
								//fi->V(0)->N() += vcg::Point3f(norm->source->GetSourceData()[norm->indices[ind] * 3],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 1],norm->source->GetSourceData()[norm->indices[ind]  * 3 + 2]).Normalize();
								fi->V(2)->N() += vcg::Point3f(normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1],normlist[0]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 2]);
								//++fi->V(2)->incidentfaces;
							}
							
							if (isvalidtext)
							{
								for(unsigned int hh = 0; hh < tet.size();++hh)
								{
									UInt32List* ls = tmp->GetPolygons(pset)->FindIndices(tet[hh]);
									fi->WT(2).t(hh) = vcg::Point2f(tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2],tet[hh]->GetSource()->GetSourceData()[(*ls)[ind] * 2 + 1]);
								}
							}

							if (isvalidnorm) fi->N() = ((fi->V(1)->P() - fi->V(0)->P()) ^ (fi->V(2)->P() - fi->V(0)->P())).Normalize();

							/*FCDGeometryPolygonsInput* posa = tmp->GetPolygons(pset)->FindInput(FUDaeGeometryInput::POSITION);
							FloatList& list = posa->source->GetSourceData();
							int dim = list.size();
							list[0] = -100.0;*/
						}
						//vm.push_back(msh);
						
						if (isvalidnorm)
						{
							vcg::tri::UpdateNormals<OpenMeshType>::PerVertexNormalized(m);
							/*for(MyMesh::VertexIterator vit = msh->vert.begin(); vit != msh->vert.end();++vit)
								vit->N() = (vit->N() / vit->incidentfaces).Normalize();*/
						}
					}
				}
			}
		}

		addinfo = inf;
		if (!amesh)
		{
			delete info->doc;
			info->doc = NULL;
			return E_NOMESH;
		}
		//doc->WriteToFile("PincoPalla.dae");
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
