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

bool Save_STL(const char * filename , bool binary =true, const char *objectname=0)
{
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
		fwrite(&fn,1,sizeof(int),fp); 
		face_iterator fi;
		Point3f p;
		unsigned short attributes=0;

		for(fi=face.begin(); fi!=face.end(); ++fi) if( !(*fi).IsD() )
		{
			// For each triangle write the normal, the three coords and a short set to zero
			p.Import(vcg::NormalizedNormal((*fi).V(0)->P(), (*fi).V(1)->P(), (*fi).V(2)->P()));
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
		face_iterator fi;	
		for(fi=face.begin(); fi!=face.end(); ++fi) if( !(*fi).IsD() )
		{
	  	// For each triangle write the normal, the three coords and a short set to zero
			p.Import(vcg::NormalizedNormal((*fi).V(0)->P(), (*fi).V(1)->P(), (*fi).V(2)->P()));
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

//@}
