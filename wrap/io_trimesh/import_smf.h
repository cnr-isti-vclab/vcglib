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
/** Function to load a Smf file
	@param filename Name of Smf file
*/
int Load_Smf( const char * filename )
{
	char buf[1024];
	FILE *fp;
	float x,y,z;
	bool one = true;
	map<int,vertex_pointer> mv;

	fp = fopen(filename,"r");
    
    if(!fp) return -1;
    
    vertex_type v;
	v.Supervisor_Flags() = 0;
	face_type f;
	f.Supervisor_Flags() = 0;
		
	while( fgets(buf,1024,fp) )
	{
		char *vf, *comm_pt;

        if((comm_pt = strstr(buf,"#")) != NULL)
            *comm_pt = '\0';
		if( (vf = strstr(buf,"v")) != NULL )
		{
			sscanf(vf+1,"%f %f %f", &x, &y, &z);
			v.P()[0] = x;
			v.P()[1] = y;
			v.P()[2] = z;
			vert.push_back(v);
		}
		else if( (vf = strstr(buf,"f")) != NULL)
		{
			if(one)
			{
				vertex_iterator vi;
				int ind;
				for(ind=1,vi=vert.begin(); vi!=vert.end(); ++vi,++ind)
					mv[ind]=&*vi;
				one = false;
			}
			int v1,v2,v3;
			sscanf(vf+1,"%d %d %d", &v1, &v2, &v3);
			f.V(0) = mv[v1];
			f.V(1) = mv[v2];
			f.V(2) = mv[v3];
			face.push_back(f);
		}
	}
	vn = vert.size();
	fn = face.size();
	return 0;
}