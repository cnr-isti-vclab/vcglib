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


/** Function to save to Smf format file.
	@param filename Name of the new Smf file
*/
void Save_Smf(const char * filename)
{
	FILE *fp;
  fp = fopen(filename,"wb");
  fprintf(fp,"#SMF \n" );

	face_iterator fi;
	vertex_iterator vi;
	map<vertex_pointer,int> index;
	int ind;

	for(ind=1,vi=vert.begin(); vi!=vert.end(); ++vi,++ind)
	{
	  fprintf(fp,"v " );
		fprintf(fp,"%f%s",(*vi).P()[0]," " );
    fprintf(fp,"%f%s",(*vi).P()[1]," " );
		fprintf(fp,"%f%s",(*vi).P()[2],"\n");
		index[&*vi] = ind;
	}

  for (fi=face.begin(); fi!=face.end(); ++fi)
  {
		fprintf(fp,"%s","f ");
		for (int j = 0; j < 3; j++)
			fprintf(fp,"%i%s",index[(*fi).V(j)]," ");
		fprintf(fp,"%s","\n");
	}

  fclose(fp);
}
