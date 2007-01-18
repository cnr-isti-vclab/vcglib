#ifndef __VCG_TRIVIAL_WALKER
#define __VCG_TRIVIAL_WALKER

namespace vcg {
  namespace tri {


// La classe Walker implementa la politica di visita del volume; conoscendo l'ordine di visita del volume
// è conveniente che il Walker stesso si faccia carico del caching dei dati utilizzati durante l'esecuzione 
// degli algoritmi MarchingCubes ed ExtendedMarchingCubes, in particolare il calcolo del volume ai vertici
// delle celle e delle intersezioni della superficie con le celle. In questo esempio il volume da processare
// viene suddiviso in fette; in questo modo se il volume ha dimensione h*l*w (rispettivamente altezza,
// larghezza e profondità), lo spazio richiesto per il caching dei vertici già allocati passa da O(h*l*w)
// a O(h*l). 

template <class MeshType, class VolumeType>
class TrivialWalker
{
private:
	typedef int VertexIndex;
  typedef typename MeshType::ScalarType ScalarType;
  typedef typename MeshType::VertexPointer VertexPointer;
	public:

  // bbox is the portion of the volume to be computed
  // resolution determine the sampling step:
  // should be a divisor of bbox size (e.g. if bbox size is 256^3 resolution could be 128,64, etc)


  void Init(VolumeType &volume)
	{
		_bbox				= Box3i(Point3i(0,0,0),volume.ISize());
		_slice_dimension = _bbox.DimX()*_bbox.DimZ();

		_x_cs = new VertexIndex[ _slice_dimension ];
		_y_cs = new VertexIndex[ _slice_dimension ];
		_z_cs = new VertexIndex[ _slice_dimension ];
		_x_ns = new VertexIndex[ _slice_dimension ];
		_z_ns = new VertexIndex[ _slice_dimension ];
	
	};

	~TrivialWalker()
	{_thr=0;}
  
	template<class EXTRACTOR_TYPE>
	void BuildMesh(MeshType &mesh, VolumeType &volume, EXTRACTOR_TYPE &extractor, const float threshold)
	{
    Init(volume);
		_volume = &volume;
		_mesh		= &mesh;
		_mesh->Clear();
    _thr=threshold;
		vcg::Point3i p1, p2;

		Begin();
		extractor.Initialize();
		for (int j=_bbox.min.Y(); j<(_bbox.max.Y()-1)-1; j+=1)
		{		 if((j%10)==0) 	printf("Marching volume z %i (%i ..%i)\r",j,_bbox.min.Y(),_bbox.max.Y());

			for (int i=_bbox.min.X(); i<(_bbox.max.X()-1)-1; i+=1)
			{
				for (int k=_bbox.min.Z(); k<(_bbox.max.Z()-1)-1; k+=1)
				{
					p1.X()=i;									p1.Y()=j;									p1.Z()=k;
					p2.X()=i+1;	p2.Y()=j+1;	p2.Z()=k+1;
          extractor.ProcessCell(p1, p2);
				}
			}
			NextSlice();
		}
		extractor.Finalize();
		_volume = NULL;
		_mesh		= NULL;
	};

	float V(int pi, int pj, int pk)
	{
    return _volume->Val(pi, pj, pk)-_thr;
	}

	bool Exist(const vcg::Point3i &p0, const vcg::Point3i &p1, VertexPointer &v)
	{ 
		int pos = p0.X()+p0.Z()*_bbox.max.X();
		int vidx;

		if (p0.X()!=p1.X()) // punti allineati lungo l'asse X
			vidx = (p0.Y()==_current_slice) ? _x_cs[pos] : _x_ns[pos];
		else if (p0.Y()!=p1.Y()) // punti allineati lungo l'asse Y
			vidx = _y_cs[pos];
		else if (p0.Z()!=p1.Z()) // punti allineati lungo l'asse Z
			vidx = (p0.Y()==_current_slice)? _z_cs[pos] : _z_ns[pos];
		else
			assert(false);

		v = (vidx!=-1)? &_mesh->vert[vidx] : NULL;
		return v!=NULL;
	}

	void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
	{ 
		int i = p1.X() - _bbox.min.X();
		int z = p1.Z() - _bbox.min.Z();
		VertexIndex index = i+z*_bbox.max.X();
		VertexIndex pos;
		if (p1.Y()==_current_slice)
		{
			if ((pos=_x_cs[index])==-1)
			{
				_x_cs[index] = (VertexIndex) _mesh->vert.size();
				pos = _x_cs[index];
                Allocator<MeshType>::AddVertices( *_mesh, 1 );
				v = &_mesh->vert[pos];
				_volume->GetXIntercept(p1, p2, v, _thr);
				return;
			}
		}
		if (p1.Y()==_current_slice+1)
		{
			if ((pos=_x_ns[index])==-1)
			{
				_x_ns[index] = (VertexIndex) _mesh->vert.size();
				pos = _x_ns[index];
				Allocator<MeshType>::AddVertices( *_mesh, 1 );
				v = &_mesh->vert[pos];
				_volume->GetXIntercept(p1, p2, v,_thr);
				return;
			}
		}
    assert(pos >=0 && pos< _mesh->vert.size());
		v = &_mesh->vert[pos];
	}
	void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
	{
		int i = p1.X() - _bbox.min.X();
		int z = p1.Z() - _bbox.min.Z();
		VertexIndex index = i+z*_bbox.max.X();
		VertexIndex pos;
		if ((pos=_y_cs[index])==-1)
		{
			_y_cs[index] = (VertexIndex) _mesh->vert.size();
			pos = _y_cs[index];
			Allocator<MeshType>::AddVertices( *_mesh, 1);
			v = &_mesh->vert[ pos ];
			_volume->GetYIntercept(p1, p2, v,_thr);
		}
		v = &_mesh->vert[pos];
	}
	void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointer &v) 
	{
		int i = p1.X() - _bbox.min.X();
		int z = p1.Z() - _bbox.min.Z();
		VertexIndex index = i+z*_bbox.max.X();
		VertexIndex pos;
		if (p1.Y()==_current_slice)
		{
			if ((pos=_z_cs[index])==-1)
			{
				_z_cs[index] = (VertexIndex) _mesh->vert.size();
				pos = _z_cs[index];
				Allocator<MeshType>::AddVertices( *_mesh, 1 );
				v = &_mesh->vert[pos];
				_volume->GetZIntercept(p1, p2, v,_thr);
				return;
			}
		}
		if (p1.Y()==_current_slice+1)
		{
			if ((pos=_z_ns[index])==-1)
			{
				_z_ns[index] = (VertexIndex) _mesh->vert.size();
				pos = _z_ns[index];
				Allocator<MeshType>::AddVertices( *_mesh, 1 );
				v = &_mesh->vert[pos];
				_volume->GetZIntercept(p1, p2, v,_thr);
				return;
			}
		}
		v = &_mesh->vert[pos];
	}

protected:
	Box3i		_bbox;

	int _slice_dimension;
	int	_current_slice;
	
	VertexIndex *_x_cs; // indici dell'intersezioni della superficie lungo gli Xedge della fetta corrente
	VertexIndex	*_y_cs; // indici dell'intersezioni della superficie lungo gli Yedge della fetta corrente
	VertexIndex *_z_cs; // indici dell'intersezioni della superficie lungo gli Zedge della fetta corrente
	VertexIndex *_x_ns; // indici dell'intersezioni della superficie lungo gli Xedge della prossima fetta 
	VertexIndex *_z_ns; // indici dell'intersezioni della superficie lungo gli Zedge della prossima fetta 

	MeshType		*_mesh;
	VolumeType	*_volume;

  float _thr;
	void NextSlice() 
	{
		memset(_x_cs, -1, _slice_dimension*sizeof(VertexIndex));
		memset(_y_cs,	-1, _slice_dimension*sizeof(VertexIndex));
		memset(_z_cs, -1, _slice_dimension*sizeof(VertexIndex));

		std::swap(_x_cs, _x_ns);
		std::swap(_z_cs, _z_ns);		
		
		_current_slice += 1;
	}

	void Begin()
	{
		_current_slice = _bbox.min.Y();

		memset(_x_cs, -1, _slice_dimension*sizeof(VertexIndex));
		memset(_y_cs, -1, _slice_dimension*sizeof(VertexIndex));
		memset(_z_cs, -1, _slice_dimension*sizeof(VertexIndex));
		memset(_x_ns, -1, _slice_dimension*sizeof(VertexIndex));
		memset(_z_ns, -1, _slice_dimension*sizeof(VertexIndex));
		
	}
};
} // end namespace 
} // end namespace 
#endif // __VCGTEST_WALKER