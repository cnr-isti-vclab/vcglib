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
Revision 1.1  2004/05/27 13:24:08  ganovelli
export_dxf created

****************************************************************************/
#ifndef __VCG_LIB_EXPORTER_DXF
#define __VCG_LIB_EXPORTER_DXF

#include <stdio.h>
/** 
This class encapsulate a filter for saving edge meshes ad polyline in DXF format.
*/
namespace vcg {
namespace edge {
namespace io {



	template <class EdgeMeshType>
class ExporterDXF{
public:
	typedef typename EdgeMeshType::VertexPointer VertexPointer;

	static bool Save( EdgeMeshType  &em, const char * filename){
		FILE * o = fopen(filename,"w");
		if(o==NULL)
			return false;

		// print header
		fprintf(o,"999\nVCGLibraryDXF\n0\nSECTION\n2\nTABLES\n0\nTABLE\n2\nLAYER\n70\n153\n0\nLAYER\n2\nthelayer\n70\n0\n62\n15\n0\nENDTAB\n0\nENDSEC\n0\nSECTION\n2\nENTITIES\n");

		vcg::edge::Pos<typename EdgeMeshType::EdgeType> et;
		typename typename EdgeMeshType::EdgePointer ep = &*em.edges.begin(),start;
		typename typename EdgeMeshType::EdgeIterator ei;

		int i=0,maxc=0,n_=0;
	 
		for(ei = em.edges.begin(); ei != em.edges.end();++ei){(*ei).ClearS();}

		for(ei = em.edges.begin(); ei != em.edges.end();++ei)
			{i=1;
				
				start = &*ei;
				et.Set(&*ei,(*ei).V(0));
				if(!et.e->IsS())
					{// nuovo contorno: trova il bordo se c'e' e posiziona li' l'hal edge
					n_++;
					do{
						ep = et.e;
						if(et.e->EEp(et.Z()) == et.e)
							break;
						et.NextE();
						}while (et.e != start);
					fprintf(o,"0\nPOLYLINE\n10\n0\n70\n0\n8\nthelayer\n");
					start = et.e;
					i=0;
										do{
							if(i++>maxc)
								maxc=i;
								et.e->SetS();
								ep = et.e;
								OutVertex(et.e->V(et.Z()), o);
								et.NextE();
						}while((et.e != ep)&&(et.e !=start)&&(i<em.en));
					fprintf(o,"0\nSEQEND\n");
					}
			}	
				
		fprintf(o,"0\nENDSEC\n0\nEOF\n");
		fclose(o);
		return 0;
		}

private:
	static void OutVertex(const VertexPointer & v, FILE* o){
					fprintf(o,"0\nVERTEX\n100\nAcDb3dPolylineVertex\n70\n32\n");
					int i;
					for(i = 0 ; i < 3; ++i)
					fprintf(o,"%d\n%f\n",(i+1)*10,v->P()[i]);
		}
	};

	};//vcg
	};//edge
	};//io

#endif