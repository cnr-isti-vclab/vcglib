/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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

/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.3  2004/03/09 21:26:47  cignoni
cr lf mismatch

Revision 1.2  2004/03/03 15:35:53  cignoni
Yet another cr lf mismatch

Revision 1.3  2004/02/19 15:28:01  ponchio
*** empty log message ***

Revision 1.2  2004/02/13 02:18:57  cignoni
Edited Comments and GPL license


****************************************************************************/

#ifndef __VCGLIB_EXPORT_STL
#define __VCGLIB_EXPORT_STL

#include <stdio.h>

namespace vcg {
namespace tri {
namespace io {

/** 
This class encapsulate a filter for opening stl (sterolitograpy) meshes.
The stl format is quite simple and rather un-flexible. It just stores, in ascii or binary the, unindexed, geometry of the faces.
*/
template <class SaveMeshType>
class ExporterSTL
{
public:
static bool Save(SaveMeshType &m, const char * filename , bool binary =true, const char *objectname=0)
{
  typedef typename SaveMeshType::FaceIterator FaceIterator;
	FILE *fp;

	fp = fopen(filename,"wb");
	if(fp==0)
		return false;

	if(binary)
	{
		// Write Header
		char *header="VCG                                                                                                  ";
		if(objectname)	strncpy(header,objectname,80);
		fwrite(header,80,1,fp);
		// write number of facets
		fwrite(&m.fn,1,sizeof(int),fp); 
		Point3f p;
		unsigned short attributes=0;
    
    FaceIterator fi;		
		for(fi=m.face.begin(); fi!=m.face.end(); ++fi) if( !(*fi).IsD() )
		{
			// For each triangle write the normal, the three coords and a short set to zero
			p.Import(vcg::NormalizedNormal(*fi));
			fwrite(p.V(),3,sizeof(float),fp);
 
			for(int k=0;k<3;++k){
				p.Import((*fi).V(k)->P());
				fwrite(p.V(),3,sizeof(float),fp);
			}
			fwrite(&attributes,1,sizeof(short),fp);
		}
	}
	else
	{
		if(objectname) fprintf(fp,"solid %s\n",objectname);
		else fprintf(fp,"solid vcg\n");

		Point3f p;
		FaceIterator fi;	
		for(fi=m.face.begin(); fi!=m.face.end(); ++fi) if( !(*fi).IsD() )
		{
	  	// For each triangle write the normal, the three coords and a short set to zero
			p.Import(vcg::NormalizedNormal(*fi));
			fprintf(fp,"  facet normal %13e %13e %13e\n",p[0],p[1],p[2]);
			fprintf(fp,"    outer loop\n");
			for(int k=0;k<3;++k){
				p.Import((*fi).V(k)->P());
				fprintf(fp,"      vertex  %13e %13e %13e\n",p[0],p[1],p[2]);			
			}
			fprintf(fp,"    endloop\n");
			fprintf(fp,"  endfacet\n");
		}
		fprintf(fp,"endsolid vcg\n");
	}
	fclose(fp);
	return true;
}
}; // end class

} // end Namespace tri
} // end Namespace io
} // end Namespace vcg



#endif