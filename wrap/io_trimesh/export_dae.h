#ifndef EXPORTER_DAE_H
#define EXPORTER_DAE_H

#include <wrap/dae/xmldocumentmanaging.h>
#include <wrap/dae/colladaformat.h>
#include <wrap/dae/util_dae.h>

namespace vcg
{
namespace tri
{
namespace io
{
	template<typename MESHMODEL>
	class ExporterDAE
	{
	public:
		static int Save(const MESHMODEL& model,const char* filename,const int mask,const QDomDocument* doc = NULL)
		{
			XMLDocumentWriter stream(filename);
			if (stream.isReliable())
			{
				XMLDocument* document = Collada::DocumentManager::createColladaDocument(model,mask);
				stream.write(*document);
				Collada::DocumentManager::destroyColladaDocument(document);
				return UtilDAE::E_NOERROR;
			}
			else 
				return UtilDAE::E_CANTSAVE;
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
