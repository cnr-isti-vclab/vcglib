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
****************************************************************************/

#ifndef __VCG_DECIMATION_COLLAPSE
#define __VCG_DECIMATION_COLLAPSE
#include<vcg\complex\tetramesh\decimation\decimation.h>
#include<vcg\complex\tetramesh\modify\edge_collapse.h>


namespace vcg{
namespace tetra{

/** \addtogroup tetramesh */
/*@{*/
/// This Class is specialization for the edge collapse

template<class TETRA_MESH_TYPE>
class Collapse: public vcg::tetra::LocalModification<TETRA_MESH_TYPE>
{
 
  /// The tetrahedral mesh type
  typedef	typename TETRA_MESH_TYPE TetraMeshType;
  /// The tetrahedron type
  typedef	typename TetraMeshType::TetraType TetraType;
	/// The vertex type
	typedef	typename TetraType::VertexType VertexType;
  /// The coordinate type
	typedef	typename TetraType::VertexType::CoordType CoordType;
  /// The scalar type
  typedef	typename TetraMeshType::VertexType::ScalarType ScalarType;
  /////the base type class
  //typedef typename vcg::tetra::LocalModification LocalMod;
  /// The HEdgePos type
  typedef Pos<TetraType> PosType;
  /// The HEdgePos Loop type
  typedef PosLoop<TetraType> PosLType;


private:

///the new point that substitute the edge
Point3<ScalarType> _NewPoint;
///the pointer to edge collapser method
vcg::tetra::EdgeCollapse<TetraMeshType> *_EC;
///mark for up_dating
char _Imark;
///the pos of collapse 
PosType pos;
///pointer to vertex that remain
VertexType *vrem;
public:
/// Default Constructor
	Collapse()
		{
      _EC=NULL;
		}

///Constructor with postype
	Collapse(PosType p)
		{
     /* Imark=mark;*/
      pos=p;
      _EC=NULL;
		}

  ~Collapse()
		{
      if (_EC!=NULL)
        delete(_EC);
		}

private:

///Return the aspect Ratio media of the tetrahedrons
///that share the adge to collapse
ScalarType _AspectRatioMedia(PosType p)
{
  PosLType posl=PosLType(p.T(),p.F(),p.E(),p.V());
  posl.Reset();
  int num=0;
  ScalarType ratio_media=0.f;
  while(!posl.LoopEnd())
  {
    ratio_media+=posl.T()->AspectRatio();
    posl.NextT();
    num++;
  }
  ratio_media=ratio_media/num;
  return (ratio_media);
}


///Modify pos and alfa to obtain the collapse that minimize the error
ScalarType _VolumePreservingError(PosType &pos,CoordType &new_point,int nsteps)
{
  VertexType *ve0=(pos.T()->V(Tetra::VofE(pos.E(),0)));
  VertexType *ve1=(pos.T()->V(Tetra::VofE(pos.E(),1)));
  vrem =ve0;
  bool ext_v0=ve0->IsB();
  bool ext_v1=ve1->IsB();

  ScalarType best_error=0.f;
   if ((ext_v0)&&(!ext_v1))
      new_point=ve0->P();
   else
   if ((!ext_v0)&&(ext_v1))
      new_point=ve1->P();
   else
   if ((!ext_v0)&&(!ext_v1))
     new_point=(ve0->P()+ve1->P())/2.f;
   else
   if ((ext_v0)&&(ext_v1))//both are external vertex
   {
    ScalarType step=1.f/(nsteps-1);
    ScalarType Vol_Original=_EC->VolumeOriginal();
    for (int i=0;i<nsteps;i++)
    {
      best_error=1000000.f;
      ScalarType alfatemp=step*((double)i);
      CoordType newPTemp=(ve0->P()*alfatemp) +(ve1->P()*(1.f-alfatemp));
      //the error is the absolute value of difference of volumes
      ScalarType error=fabs(Vol_Original-_EC->VolumeSimulateCollapse(pos,newPTemp));
      if(error<best_error)
      {
       new_point=newPTemp;
       best_error=error;
      }
    }
   }
   return (best_error);
}



public:


  ScalarType ComputePriority()
  { 
    return (_AspectRatioMedia(this->pos));
  }

  ScalarType ComputeError()
  {
     if (_EC==NULL)
    {
      _EC=new vcg::tetra::EdgeCollapse<TetraMeshType>();
      _EC->FindSets(pos);
    }
      return (_VolumePreservingError(pos,_NewPoint,5));
  }

  void Execute()
  {
    if (_EC==NULL)
    {
      _EC=new vcg::tetra::EdgeCollapse<TetraMeshType>();
      _EC->FindSets(pos);
    }
    _EC->DoCollapse(pos,_NewPoint);
  }
  
  bool PreserveTopology()
  {
    if (_EC==NULL)
    {
      _EC=new vcg::tetra::EdgeCollapse<TetraMeshType>();
      _EC->FindSets(pos);
    }
    return(_EC->CheckPreconditions(pos,_NewPoint));
  }
  
 HeapRetType UpdateHeap()
  {
    HeapRetType h_ret;
    assert(!vrem->IsD());
    VTIterator<TetraType> VTi(vrem->VTb(),vrem->VTi());
    while (!VTi.End())
    {
      for (int j=0;j<6;j++)
      {
        vcg::tetra::Pos<TetraType> p=Pos<TetraType>(VTi.Vt(),Tetra::FofE(j,0),j,Tetra::VofE(j,0));
        h_ret.push_back(HeapRetElem(vcg::tetra::MTEdgeCollapse,p));
      }
      VTi++;
    }
    return (h_ret);
  }

};
}//end namespace tetra
}//end namespace vcg
#endif