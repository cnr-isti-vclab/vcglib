#ifndef __VCGLIB_EXPORTERDAE
#define __VCGLIB_EXPORTERDAE

#include <FCollada.h>
#include <FCDocument/FCDocument.h>
#include <FCDocument/FCDLibrary.h>
#include <FCDocument/FCDGeometry.h>
#include <FCDocument/FCDGeometryMesh.h>
#include <FCDocument/FCDGeometrySource.h>
#include <FCDocument/FCDGeometryPolygons.h>
#include <FCDocument/FCDSceneNode.h>
#include <FCDocument/FCDGeometryInstance.h>

namespace vcg {
namespace tri {
namespace io {

//FCollada Library assumes that SaveMeshType::ScalarType is always a float

template<typename SaveMeshType>
class ExporterDAE
{
public:
	static int Save(SaveMeshType &m, const char * filename)
	{
		unsigned int ncomp = sizeof(SaveMeshType::CoordType) / sizeof(SaveMeshType::ScalarType);
		FCDocument* doc = new FCDocument();
		FCDGeometryLibrary* geolib = doc->GetGeometryLibrary();
		FCDGeometry* geo = geolib->AddEntity();
		geo->SetDaeId("vcgmesh");
		FCDGeometryMesh* mesh = geo->CreateMesh();

		FCDGeometrySource* vsource = mesh->AddVertexSource();
		vsource->SetDaeId("vcg-position");
		vsource->SetSourceType(FUDaeGeometryInput::POSITION);
		FCDGeometryPolygons* polyset = mesh->AddPolygons();
		//for(UInt32List::iterator it = polyinput->indices.begin();it != polyinput->indices.end();++it)
		//	*it = 100;
		
		FloatList plist(m.vert.size() * ncomp);
		
		SaveMeshType::CoordType point;
		unsigned int ii = 0;
		for(FloatList::iterator it = plist.begin();it != plist.end();++it)
		{
			if ((ii % ncomp) == 0)
			{
				point = m.vert[ii / ncomp].P();
			}
			(*it) = point[ii % ncomp];
			++ii;
		}
		
		vsource->SetSourceData(plist,ncomp);

		unsigned int jj = 0;
		FCDGeometryPolygonsInput* ps = mesh->GetPolygons(0)->FindInput(FUDaeGeometryInput::POSITION);
		//ps->indices.resize(m.face.size() * 3);
		for(SaveMeshType::FaceIterator itf = m.face.begin();itf != m.face.end();++itf)
		{
			if( !(*itf).IsD() )
			{
				polyset->AddFace(3);
				ps->indices[jj] = (itf->V(0) - &(m.vert[0]));
				ps->indices[jj + 1] = (itf->V(1) - &(m.vert[0]));
				ps->indices[jj + 2] = (itf->V(2) - &(m.vert[0]));
				jj += 3;
			}
		}

		FloatList nlist;

		FCDGeometrySource* nsource = mesh->AddSource();
		nsource->SetDaeId("vcg-wedgenormal");
		nsource->SetSourceType(FUDaeGeometryInput::NORMAL);
		//in the oldmesh format of vcg library normal is always a Point3
		nlist.resize(m.face.size() * 3 * 3);
		
		unsigned int cc = 0;
		jj = 0;
		
		FCDGeometryPolygonsInput* nm = polyset->AddInput(nsource,1);
		//nm->semantic = FUDaeGeometryInput::NORMAL;
		nm->indices.resize(m.face.size() * 3);
		
		for(SaveMeshType::FaceIterator fit = m.face.begin();fit != m.face.end();++fit)
		{
			if( !(*fit).IsD() )
			{
				vcg::Point3<SaveMeshType::ScalarType> norm = (((fit->V(1))->P() - (fit->V(0))->P()) ^ ((fit->V(2))->P() - (fit->V(0))->P())).Normalize();
				for(unsigned int vv = 0; vv < 3;++vv)
				{
					nm->indices[jj + vv] = cc / 3;
					for(unsigned int hh = 0; hh < 3;++hh)
					{
						nlist[cc] = norm[hh];
						++cc;
					}
				} 
				jj += 3;
			}
		}
		nsource->SetSourceData(nlist,3);
		
		FCDSceneNode* root = doc->AddVisualScene();
		root->SetDaeId("vcg-scenenode");
		FCDSceneNode* scenenod = root->AddChildNode();
		scenenod->AddInstance(geo);
		//root->
		doc->WriteToFile(FUStringConversion::ToFString(filename));
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