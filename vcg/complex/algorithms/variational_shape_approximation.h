#ifndef VCGLIB_REGION_GROWING_VERTEX_
#define VCGLIB_REGION_GROWING_VERTEX_

#include <list>
#include <vector>
#include <deque>

#include <vcg/complex/allocate.h>
#include <vcg/complex/append.h>
#include <vcg/space/fitting3.h>
#include <vcg/space/color4.h>
//#include <vcg/complex/algorithms/update/quality.h>

namespace vcg{
	namespace tri{

/**
Element and associated error
*/
template <class ElemType>
struct ElemError{
	ElemType * f;
	float val;
	ElemError(){};
	ElemError(ElemType * _f,double _val):f(_f),val(_val){}
	const bool operator <(const ElemError & o) const {return val < o.val;}
};

/**
Base class to define a region of a mesh
*/
template<class MeshType, class ElemType>
struct RegionBase{
	RegionBase():isd(false),size(0){}
	typedef typename MeshType MeshType;
	typedef typename ElemType::CoordType CoordType;
    typedef typename ElemType::CoordType::ScalarType ScalarType;
	typedef typename std::list<RegionBase*>::iterator AdjIterator;
	typedef ElemError<ElemType> ElemError;
	typedef typename std::vector<ElemType*>::iterator ElemIte;

	// elem belonging to the region
	std::vector<ElemType*> elems;

	// adjacent regions
	std::list<RegionBase*> adj,nx_adj;

	ElemError  worst;

	int id,size;
	int cleansize;
        ScalarType epsilon;
        ScalarType approx_err;
        ScalarType approx_var;
    int changed;
	bool isd;


	void Clean(){
		adj.sort(); 
		adj.unique();
		nx_adj.clear();
		}

	void Connect( RegionBase * tr){
		adj.push_back(tr);
	}

	 
	void UpdateError( ){
        ElemIte vi;

        float err=0,var=0,max_err=-1.0;
        for(vi=elems.begin(); vi != elems.end(); ++vi)
		{
			float e =Evaluate(**vi);
			err+=e;
			if(e > max_err) {
					worst.val = e;
					worst.f = *vi;
					max_err = e;
			 }
		}
        approx_err =err;
        var = 0.0;
        approx_var = var;
}

	 
 ScalarType Mergeable(  RegionBase & tr){
/*
two regions are mergeable if:
- they have the same planes
- they have an extreme in common
or
- one is surrounded by the other which is much bigger
*/

	//if(!Similar(this->p,tr.p))
	//	return false;

	vcg::Plane3<ScalarType> plane;
	ScalarType fit_error,new_error=0;

	if( (this->adj.size() == 1) && (tr.elems.size()/float(this->elems.size()) > 5) ||
		(tr.adj.size() == 1) && (this->elems.size()/float(tr.elems.size()) > 5)  
		)
		return true;


	std::vector<ElemType*> merged;
	merged.insert(merged.end(),this->elems.begin(),this->elems.end());
	merged.insert(merged.end(),tr.elems.begin(),tr.elems.end());
	Fit(merged,plane,fit_error);
	for(std::vector<ElemType*> ::iterator vi= merged.begin(); vi != merged.end(); ++vi)
		new_error+= Evaluate(**vi);
	new_error/=merged.size();

	return new_error <  this->approx_err * (elems.size()/float(merged.size()))+tr.approx_err * (tr.elems.size()/float(merged.size())) ;

}

 /// pure virtual function to define how to fit the primitive onto this region 
 virtual void  Fit(std::vector<ElemType*>& verts, vcg::Plane3f & plane, float & err) = 0;
 /// pure virtual function to define how to fit the evaluate the actio of adding and element "e" to this region
 virtual ScalarType Evaluate( ElemType & e) = 0;
 /// pure virtual function to init the primitive with the element 
 virtual void Init(ElemType*) = 0;

};

/**
Base class to define a region grower, i.e. a floo algorithm to partion the mesh in regions
*/
template < class RegionType>
struct RegionGrower{
  	typedef typename RegionType RegionType;
  	typedef typename RegionType::MeshType MeshType;

  	typedef typename RegionType::ElemType ElemType;
	typedef typename RegionType::ScalarType ScalarType;

	typedef typename std::list<RegionType> ::iterator TriRegIterator;
    typedef typename std::list<RegionBase<MeshType,ElemType>* >::iterator AdjIterator;
    typedef typename RegionType::ElemError ElemError;

	struct CandiElem{
		ElemType * f;
		float val;
		RegionType * r;
		CandiElem(){};
		CandiElem(ElemType * _f,double _val,RegionType * _r):f(_f),val(_val),r(_r){}
		const bool operator <(const CandiElem & o) const {return val < o.val;}
	};

     struct ErrorEval {

        void Add(const float & v){
            {
            if(v < samples[i_min]) i_min = 0; else {++i_min;
                if(v > samples[i_max])  i_max = 0;else ++i_max;
            }
            if(i_min == samples.size()) {i_min = 0; for(int i= 0; i < samples.size()-1; ++i) if(samples[i]<samples[i_min]) i_min = i; ++i_min; }
            if(i_max == samples.size()) {i_max = 0; for(int i= 0; i < samples.size()-1; ++i) if(samples[i]>samples[i_max]) i_max = i; ++i_max;}
        }
            samples.pop_back();samples.push_front(v);
            boxes.pop_back(); boxes.push_front(vcg::Point2f(samples[i_min],samples[i_max]));
            ++ns ;
        }

        float BoxOverlap(){
            float maxsize =  std::max( boxes.back()[1]-boxes.back()[0], (*boxes.begin())[1]-(*boxes.begin())[0]);
            float overlap =  std::max(0.f, std::min(boxes.back()[1],(*boxes.begin())[1])-std::max(boxes.back()[0],(*boxes.begin())[0]));
            assert(overlap <= maxsize);
            return  (maxsize  > 0.f)?overlap / maxsize:0.0;
        }
        float RelativeDecrease(){
			if(ns<4)  return std::numeric_limits<float>::max();
			return fabs(samples[std::min<unsigned int>(samples.size()-1,ns)]-samples[0]);
        }

        void Init(int n ){
			
			samples.clear();
			boxes.clear();
			for(int i = 0 ; i < n; ++i) samples.push_front(std::numeric_limits<float>::max());
			for(int i = 0 ; i < n; ++i) boxes.push_front(vcg::Point2f(n,n-i));
			i_max = i_min = 0;
			ns = 0;
		}

        private:
            int i_min,i_max;            // index of min and max element in the queue
            std::deque<float> samples;
            std::deque<vcg::Point2f> boxes;
            int ns;

        };


	RegionGrower():lastAdded(NULL),ave_error(-1),ave_var(-1),n_steps(0){}

	std::list<RegionType> regions;
	std::vector<RegionType*> toTunnel;
	std::vector<RegionType*> workingset;

    int n_faces,target_max_regions;
	ElemType * lastAdded;
    ScalarType	 ave_error // average error (over single regions' error)
                ,ave_var// average variance (over single regions' error)
                ,changed	// faces that have changed from previous step ([0..1))
                ,err
                ,target_error; // taget error

	ElemType * worst;
    ScalarType worst_err ;
    MeshType * m;

    ErrorEval erroreval;
	unsigned int n_steps;//  number of steps done
	std::vector<CandiElem> elemheap;
	std::vector<ElemError  > elemerr;

	void Init(MeshType & mesh, int n_seeds,int max_regions, float max_err){
		erroreval.Init(10);
		m = &mesh;
		target_error = mesh.bbox.Diag()*max_err;
		target_max_regions = max_regions;
		regions.clear();
	}


	/// add a region
	void AddRegion(const RegionType  & r){regions.push_back(r);}

	/// remove a region
	void DeleteRegion(const typename std::list<RegionType>::iterator & ri){std::remove(ri);}

				
	void PushHeap(std::vector<ElemType*> & candi, RegionType & r){
					typename std::vector<ElemType*>::iterator ci;
				for(ci = candi.begin(); ci != candi.end(); ++ci)
						{
							elemheap.push_back(CandiElem( *ci,-r.Evaluate(*(*ci)), &r));
							push_heap(elemheap.begin(),elemheap.end());
						}
	}

	/// make two regions  adjacent
	void Connect(RegionType  * r1,RegionType  * r2){			 
		assert(r1!=r2);
		r1->Connect(r2);
		r2->Connect(r1);
	}

	/// initialize a region
	void CreateRegion(ElemType * vi){
		AddRegion(RegionType());
		RegionType & tr =regions.back();
		AddElemToRegion(tr,vi);
		tr.Refit();
		tr.color = vcg::Color4b::Scatter(2000,(int)regions.size());	
	}

	void AddElemToRegion( RegionType & r, ElemType * v){
		r.elems.push_back(v);
		B(v) = (RegionType*) &r;
		if(B(v)!=B_old(v)) ++r.changed;
		B_old(v)  = B(v);
		++r.size;
	}

	/// for each region take the candidates and fill in elemheap
	void Refill( ){
		elemheap.clear();
		typename std::list<RegionType>::iterator ri;
		std::vector<ElemType*>  candi;
		for(ri = regions.begin(); ri != regions.end(); ++ri) if(!(*ri).isd)
			{
				candi.clear();
				Candidates((*ri),candi);
				PushHeap(candi,*ri);
			}
		std::make_heap(elemheap.begin(),elemheap.end());
	}


	/// check if the iteration converged
	   bool IsRelaxed(){
		ScalarType _ave_error=0;
		ScalarType _ave_var= 0;
		ScalarType _changed = 0.0;
		int nr=0;

		typename std::list<RegionType>::iterator ri;
		worst=NULL;;
		worst_err = -1;
		for(ri = regions.begin(); ri != regions.end(); ++ri) if(!(*ri).isd){
				++nr;
				n_faces+=(*ri).elems.size();
				(*ri).UpdateError();
				_ave_error+=(*ri).approx_err;
				_ave_var+=(*ri).approx_var;
				_changed+=(*ri).changed;
				(*ri).changed=0;
				if((*ri).worst.val*(*ri).size > worst_err){
										worst = (*ri).worst.f;
										worst_err = (*ri).worst.val*(*ri).size;
				}
		}


		_ave_error/=nr;
		_ave_var/=nr;
		_changed/=n_faces;
		n_faces = 0;

		erroreval.Add(_ave_error);
		printf("Err: %f ov: %f Dec: %f \n",_ave_error,erroreval.BoxOverlap(),erroreval.RelativeDecrease());


		return  (erroreval.BoxOverlap() > 0.95) || (erroreval.RelativeDecrease() < 0.01*_ave_error);
	}
	

			/// merge two regions
		   void Merge(RegionType & r0,RegionBase<typename RegionType::MeshType,typename RegionType::ElemType> & r1){
					assert(!r1.isd);
				typename RegionType::ElemIte  vi;
				AdjIterator ai;
				for(vi = r1.elems.begin();vi != r1.elems.end(); ++vi)
					AddElemToRegion(r0,(*vi));
				for(ai= r1.adj.begin(); ai != r1.adj.end();++ai) 
					if( !(*ai)->isd && (*ai)!=&r0)
					 r0.nx_adj.push_back(*ai);

				r1.elems.clear();
				r1.isd = true;
				r0.Refit();
			}
		  
		   /// clean degeneracies: look for very small regions in the middle of a single regions (and similar cases)
		   bool IsToTunnel(RegionType * r){
			   return (r->elems.size()<2); 
			   
		   }
		   /// clean degeneracies: look for very small regions in the middle of a single regions (and similar cases)
		   void ComputeToTunnel(){
			   typename std::list<RegionType>::iterator ri;
			   for(ri = regions.begin(); ri != regions.end(); ++ri)
				   if(IsToTunnel(&(*ri)))
					   toTunnel.push_back(&(*ri));
		   }



	   void GrowStepOnce(){
		   if(elemheap.empty()) return;
			CandiElem cf;	
			std::vector<ElemType*> toAdd;
			std::pop_heap(elemheap.begin(),elemheap.end());
			cf = elemheap.back();
			//printf("err:%f\n",cf.val);		
			elemheap.pop_back();
			if (B(cf.f)==NULL){
				AddElemToRegion( *cf.r,cf.f); // ads to region
				lastAdded = &*cf.f;

				for(int i =0; i <  cf.f->SizeNeigh();++i)
					 if((B(cf.f ->Neigh(i)) == NULL) )
						toAdd.push_back(cf.f->Neigh(i));
					 else
						 if(B(cf.f ->Neigh(i))  != B(cf.f))
						Connect((RegionType*)B(cf.f->Neigh(i)),(RegionType*)cf.r);
				PushHeap(toAdd,*cf.r);
			}
			else
			{
				if ( B(cf.f) != (RegionType*)(cf.r) )
				Connect((RegionType*)B(cf.f),(RegionType*)cf.r);
			}
	   }
		/// Execute and iteration
		void GrowStep(){
			CandiElem cf;	
			typename std::list<RegionType>::iterator ri;

			n_steps++;
			for(ri = regions.begin(); ri != regions.end(); ++ri) if(!(*ri).isd)
					(*ri).Clean();

			/*printf("Heap size %d\n", elemheap.size());*/
			while(!elemheap.empty() ){
					std::vector<ElemType*> toAdd;
					std::pop_heap(elemheap.begin(),elemheap.end());
					cf = elemheap.back();
					//printf("err:%f\n",cf.val);		
					elemheap.pop_back();
					if (B(cf.f)==NULL){
						AddElemToRegion( *cf.r,cf.f); // ads to region
						lastAdded = &*cf.f;

						for(int i =0; i <  cf.f->SizeNeigh();++i)
							 if((B(cf.f ->Neigh(i)) == NULL) )
								toAdd.push_back(cf.f->Neigh(i));
							 else
								 if(B(cf.f ->Neigh(i))  != B(cf.f))
								Connect((RegionType*)B(cf.f->Neigh(i)),(RegionType*)cf.r);
						PushHeap(toAdd,*cf.r);
					}
					else
					{
						if ( B(cf.f) != (RegionType*)(cf.r) )
						Connect((RegionType*)B(cf.f),(RegionType*)cf.r);
					}
			}
			int h = (int) elemheap.size();
		 
			
			//printf("----> %d\n",h);
		}

		/// try to merge similare adjacent regions
		unsigned int MergeStep(){
			 TriRegIterator tri;
			 unsigned int merged = 0;
			 typename RegionType::AdjIterator ai;
			 for(tri = regions.begin(); tri != regions.end(); ++tri) if(!(*tri).isd) (*tri).Clean();
			 for(tri = regions.begin(); tri != regions.end(); ++tri)if(!(*tri).isd)
					for(ai = (*tri).adj.begin(); ai != (*tri).adj.end(); ++ai) if(!(*ai)->isd)
						if((*tri).Mergeable(*(*ai))){
							Merge((*tri),*(*ai));
							merged++;
						}


			 for(tri = regions.begin(); tri != regions.end(); ++tri){
				 for(ai = (*tri).nx_adj.begin(); ai != (*tri).nx_adj.end();++ai)
					 if(!(*ai)->isd) 
						 (*tri).adj.push_back(*ai);
				 (*tri).adj.sort();
				 (*tri).adj.unique();
				 }
			return merged;
		}

	/// re-init all for a new itearation
    bool Restart(){
            std::vector<ElemType*>  candi;
            TriRegIterator ri;
            elemheap.clear();

			printf("#regions: %d",this->regions.size());
			ComputeToTunnel();
			if(IsRelaxed()){
				printf("R worst_err:%f, target_error %f \n",worst_err ,target_error);
                if( (worst_err <= target_error) || (regions.size() >= target_max_regions))
                    return false;
                else
                {	printf("Add Region \n");
                    erroreval.Init(10);
                    ElemError   wrs;
                    wrs.f = worst;
                    wrs.val = worst_err;
                    printf("worst triangle error %f\n",worst_err);

                    CreateRegion(wrs.f);// CreateRegion changes wr->r

                    // reset state variables
                     ave_error=-1;
                     ave_var=-1;
                     err=0.0;
                     changed=0;
                }
				printf("\n");
            }
			
			
			//
			printf("tunnelling  %d regions\n",toTunnel.size());
			std::vector<RegionType*>::iterator it;
			for(it = toTunnel.begin(); it != toTunnel.end(); ++it)
			{
				RegionType * r = *it;
				typename RegionType::ElemType *  e;			   
				typename RegionType::ElemIte  vi;
				for(vi = r->elems.begin(); vi != r->elems.end(); ++vi) B(*vi) = NULL;
				r->elems.clear();
				do{e = RandomElem();} while( (B(e)!=NULL) && !IsToTunnel(B(e))); 
				AddElemToRegion(*r,e);
				r->Refit();
				//r->color = vcg::Color4b::Scatter(2000,(int)regions.size());
				//r->Init(e);
				Candidates(*r ,candi);    // take the  candidatees
				PushHeap(candi,(*r ));    // put the faces on to the heap
			}
			toTunnel.clear();
 			//

            for(ri = regions.begin(); ri != regions.end(); )
				if((*ri).isd)
					ri = regions.erase(ri);
				else
					++ri;

             for(ri = regions.begin(); ri != regions.end(); ++ri)
             {
                     candi.clear();
                     (*ri).Refit();            // fit a plane to the region
                     Restart(*ri);             // clear stuff in the region, move the seed to the best fitting to the plabne
                     Candidates(*ri,candi);    // take the (three) faces candidatees
                     PushHeap(candi,(*ri));    // put the faces on to the heap
             }



            return true;
        }


    	/// re-init a specific region
        void Restart(RegionType &r){
                if(!r.elems.empty()){
			 		    r.Refit();
                        r.size=0;
                        //float b_err = r.Evaluate(*(*r.elems.begin())),err;
                         ElemType* b_vert =(*r.elems.begin());
                         typename RegionType::ElemIte  vi;

                        //for( vi = r.elems.begin(); vi != r.elems.end(); ++vi)
                        //{
                        //        err = r.Evaluate(**vi);
                        //        if(err < b_err)
                        //        {
                        //                b_err = err;
                        //                b_vert = *vi;
                        //        }
                        //}

						b_vert = r.elems[r.elems.size()/2];
                        for( vi = r.elems.begin(); vi != r.elems.end(); ++vi) B(*vi) = NULL;
                        r.elems.clear();
                        r.adj.clear();
                        AddElemToRegion(r,b_vert);
						r.Init(b_vert);
                }
            }

		// after init, run all the process
		void MakeCharts(){
			this->Refill();
				while(this->Restart()){
					//do{ 
						this->GrowStep();
						

					//} while(Restart());
				}
				//while(this->MergeStep());
			}

	
		/// pick a random element
		virtual  typename RegionType::ElemType * RandomElem() = 0; 

		/// given an element returns the region to which it is associated (or NULL )
		virtual  RegionType * &B(ElemType *)  = 0;
		/// given an element returns the region to which it was associated (or NULL )
		virtual  RegionType * &B_old(ElemType *)  = 0;
		/// return in c the canidates to be added to the region r
		virtual void   Candidates(RegionType & r, std::vector< typename RegionType::ElemType*> & c) = 0;
};


template < class RegionType>
struct RegionGrowerVertex: public RegionGrower<RegionType> {
  	typedef typename RegionType RegionType;
  	typedef typename RegionType::MeshType MeshType;

  	typedef typename RegionType::ElemType ElemType;
	typedef typename RegionType::ScalarType ScalarType;
	typedef typename std::list<RegionType> ::iterator TriRegIterator;

	typename MeshType:: template  PerVertexAttributeHandle<RegionType*> attr_r;
	typename MeshType:: template  PerVertexAttributeHandle<RegionType*> attr_r_old;

	RegionType * &B(ElemType *e){return (RegionType *)attr_r[e];};
	RegionType * &B_old(ElemType *e){return (RegionType *)attr_r_old[e];}

    void Init(MeshType & mesh, int n_seeds,int max_regions, float max_err){
		RegionGrower<RegionType>::Init(mesh, n_seeds,max_regions,max_err);

		attr_r = vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<RegionType*> (*m,"r");
        if(!vcg::tri::Allocator<MeshType>::IsValidHandle(*m,attr_r))
                attr_r = vcg::tri::Allocator<MeshType>::template AddPerVertexAttribute<RegionType*> (*m,"r");

        attr_r_old = vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<RegionType*> (*m,"r_old");
        if(!vcg::tri::Allocator<MeshType>::IsValidHandle(*m,attr_r_old))
                attr_r_old = vcg::tri::Allocator<MeshType>::template AddPerVertexAttribute<RegionType*> (*m,"r_old");

      

         for(int i = 0; i < m->vert.size(); ++i){
                 attr_r[i] = NULL;
                attr_r_old[i] = NULL;
                if( (i%(m->vn/n_seeds))==0)
                            CreateRegion(&m->vert[i]);

            }
    }

        void   Candidates(RegionType & r, std::vector< typename RegionType::ElemType*> & c){
                typename RegionType::VertexIterator vi;

                for(vi = r.elems.begin(); vi!= r.elems.end(); ++vi)
					 for(int i =0; i <  (*vi)->SizeNeigh();++i)
						 if((attr_r[(*vi)->Neigh(i)] == NULL) )
							c.push_back((*vi)->Neigh(i));

        }

 
		///Export the regions as a collection of meshes
		void Export(std::vector<MeshType*> & submeshes){

			MeshType::PerVertexAttributeHandle<char> sel =  Allocator<MeshType>::AddPerVertexAttribute<char>(*m,"sel");
			assert(Allocator<MeshType>::IsValidHandle<char>(*m,sel));
			for( MeshType::VertexIterator vi = m->vert.begin();vi!= m->vert.end(); ++vi) 
				sel[(*vi)] = (*vi).IsS();

			for(TriRegIterator ti = regions.begin(); ti != regions.end(); ++ti){
				submeshes.push_back(new MeshType());
				for(std::vector<MeshType::VertexPointer>::iterator vi = (*ti).elems.begin();vi!=(*ti).elems.end(); ++vi)
					{(*vi)->SetS(); }
				
				Append<MeshType,MeshType>::Mesh( *submeshes.back(),*m,true);

				for(std::vector<MeshType::VertexPointer>::iterator vi = (*ti).elems.begin();vi!=(*ti).elems.end(); ++vi) 
					{(*vi)->ClearS();}
			}

			for( MeshType::VertexIterator vi = m->vert.begin();vi!= m->vert.end(); ++vi)
				 if(sel[(*vi)]) (*vi).SetS(); else (*vi).ClearS();

		}

		typename MeshType::VertexType * RandomElem(){
			int id = std::max<int>(0,std::min<int>((rand()/float(RAND_MAX))*m->vert.size()-1,m->vert.size()-1));
			return &m->vert[id];
		}

	};



template < class RegionType>
struct RegionGrowerFace: public RegionGrower<RegionType> {
  	typedef typename RegionType RegionType;
  	typedef typename RegionType::MeshType MeshType;

  	typedef typename RegionType::ElemType ElemType;
	typedef typename RegionType::ScalarType ScalarType;
	typedef typename std::list<RegionType> ::iterator TriRegIterator;

	typename MeshType:: template  PerFaceAttributeHandle<RegionType*> attr_r;
	typename MeshType:: template  PerFaceAttributeHandle<RegionType*> attr_r_old;

	RegionType * &B(ElemType *e){return (RegionType *)attr_r[e];};
	RegionType * &B_old(ElemType *e){return (RegionType *)attr_r_old[e];}
	  
    void Init(MeshType & mesh, int n_seeds,int max_regions, float max_err){
		RegionGrower<RegionType>::Init(mesh, n_seeds,max_regions,max_err);

		attr_r = vcg::tri::Allocator<MeshType>::template GetPerFaceAttribute<RegionType*> (*m,"r");
        if(!vcg::tri::Allocator<MeshType>::IsValidHandle(*m,attr_r))
                attr_r = vcg::tri::Allocator<MeshType>::template AddPerFaceAttribute<RegionType*> (*m,"r");

        attr_r_old = vcg::tri::Allocator<MeshType>::template GetPerFaceAttribute<RegionType*> (*m,"r_old");
        if(!vcg::tri::Allocator<MeshType>::IsValidHandle(*m,attr_r_old))
                attr_r_old = vcg::tri::Allocator<MeshType>::template AddPerFaceAttribute<RegionType*> (*m,"r_old");

         for(int i = 0; i < m->face.size(); ++i){
                 attr_r[i] = NULL;
                attr_r_old[i] = NULL;
                if( (i%(m->fn/n_seeds))==0)
                            CreateRegion(&m->face[i]);

            }
    }

        void   Candidates(RegionType & r, std::vector< typename RegionType::ElemType*> & c){
                typename RegionType::FaceIterator fi;

                for(fi = r.elems.begin(); fi!= r.elems.end(); ++fi)
					 for(int i =0; i <  (*fi)->VN();++i)
						 if((attr_r[(*fi)->FFp(i)] == NULL) )
							c.push_back((*fi)->FFp(i));

        }

 
		///Export the regions as a collection of meshes
		void Export(std::vector<MeshType*> & submeshes){
			MeshType::PerFaceAttributeHandle<char> sel =  Allocator<MeshType>::AddPerFaceAttribute<char>(*m,"sel");
			assert(Allocator<MeshType>::IsValidHandle<char>(*m,sel));
			for( MeshType::FaceIterator fi = m->face.begin();fi!= m->face.end(); ++fi) 
				sel[(*fi)] = (*fi).IsS();

			for(TriRegIterator ti = regions.begin(); ti != regions.end(); ++ti){
				submeshes.push_back(new MeshType());
				for(std::vector<MeshType::FacePointer>::iterator fi = (*ti).face.begin();fi!=(*ti).face.end(); ++fi)
					{(*fi)->SetS(); for(unsigned int i = 0; i < 3; ++i) (*fi)->V(i)->SetS();}
				
				Append<MeshType,MeshType>::Mesh( *submeshes.back(),*m,true);

				for(std::vector<MeshType::FacePointer>::iterator fi = (*ti).face.begin();fi!=(*ti).face.end(); ++fi) 
					{(*fi)->ClearS();for(unsigned int i = 0; i < 3; ++i) (*fi)->V(i)->ClearS();}
			}

			for( MeshType::FaceIterator fi = m->face.begin();fi!= m->face.end(); ++fi)
				 if(sel[(*fi)]) (*fi).SetS(); else (*fi).ClearS();
		}

		typename MeshType::FaceType * RandomElem(){
			int id = std::max<int>(0,std::min<int>((rand()/float(RAND_MAX))*m->face.size()-1,m->face.size()-1));
			return &m->face[id];
		}

	};


/** 
	Planar region made of vertices
*/
template<class MeshType>
struct PlanarRegionVertex: public RegionBase<MeshType,typename MeshType::VertexType> {
public:

	typedef PlanarRegionVertex RegionType;
	typedef typename MeshType::VertexType ElemType;
	typedef typename std::vector<ElemType*>::iterator VertexIterator;

	PlanarRegionVertex(): RegionBase<MeshType,typename MeshType::VertexType>(){};


	// planes characterizing the region
    vcg::Plane3<ScalarType,true> p;

	CoordType center;
	vcg::Color4b color;

	ScalarType	PlaneFittingError(std::vector<CoordType> &  samples,vcg::Plane3<float> p);

	// evaluate the gain if the triangle is added to the region
    ScalarType Evaluate( ElemType & f);

	// refit the plane
	void Fit(std::vector<ElemType*>& faces, vcg::Plane3f & plane, float & err);
	void Refit();
	void Init(ElemType * );
};


        template<class MeshType>
typename PlanarRegionVertex<MeshType>::ScalarType
        PlanarRegionVertex< MeshType>:: PlaneFittingError( std::vector<CoordType> &  samples,vcg::Plane3<float> p){

                typename PlanarRegionVertex<MeshType>::ScalarType err =0.0;
                typename std::vector< CoordType>::iterator  si;
                for(si = samples.begin(); si != samples.end(); ++si)
                        err += fabs((*si)*p.Direction()-p.Offset());
                return err/samples.size();
}

// evaluate the gain if the vertex is added to the region
template<class MeshType>
typename PlanarRegionVertex<MeshType>::ScalarType PlanarRegionVertex< MeshType>::Evaluate( ElemType & e){
	//return  fabs(vcg::Distance(p,e.P()));
	return (1-fabs(e.N()*p.Direction()));
}

	// refit the planes
template<class MeshType>
void PlanarRegionVertex< MeshType>::Fit(std::vector<ElemType*>& verts, vcg::Plane3f & plane, float & err){
	VertexIterator vi;
	std::vector<CoordType>  samples;

	if(verts.size()<3){
		samples.push_back( (*verts.begin())->P());
		for(int i  = 0; i < (*verts.begin())->SizeNeigh(); ++i)
			samples.push_back((*verts.begin())->Neigh(i)->P());
		 vcg::PlaneFittingPoints<typename MeshType::ScalarType >(samples,plane);
	}else
	 if(verts.size()==3){
		 plane.Init(verts[0]->P(),verts[1]->P(),verts[2]->P());
	 }else{
			for( vi = verts.begin(); vi != verts.end(); ++vi)
				samples.push_back( (*(*vi)).P());
		 vcg::PlaneFittingPoints<typename MeshType::ScalarType >(samples,plane);
	 }

	 err = PlaneFittingError(samples,plane);
}

template<class MeshType>
void PlanarRegionVertex< MeshType>::Refit(){
	
	Fit(elems,p,this->approx_err);
}

template<class MeshType>
void PlanarRegionVertex< MeshType>::Init(ElemType * b){
	std::vector<CoordType>  samples;
	samples.push_back(b->P());
	for(int i  = 0; i < b->SizeNeigh(); ++i)
		samples.push_back(b->Neigh(i)->P());
	 vcg::PlaneFittingPoints<typename MeshType::ScalarType >(samples,this->p); 
}


/** 
	Planar region made of vertices
*/
template<class MeshType>
struct PlanarRegionFace: public RegionBase<MeshType,typename MeshType::FaceType> {
public:

	typedef PlanarRegionFace RegionType;
	typedef typename MeshType::FaceType ElemType;
	typedef typename std::vector<ElemType*>::iterator FaceIterator;

	PlanarRegionFace(): RegionBase<MeshType,typename MeshType::FaceType>(){};


	// planes characterizing the region
    vcg::Plane3<ScalarType,true> p;

	CoordType center;
	vcg::Color4b color;

	ScalarType	PlaneFittingError(std::vector<CoordType> &  samples,vcg::Plane3<float> p);

	// evaluate the gain if the triangle is added to the region
    ScalarType Evaluate( ElemType & f);

	// refit the plane
	void Fit(std::vector<ElemType*>& faces, vcg::Plane3f & plane, float & err);
	void Refit();
	void Init(ElemType * );
};


template<class MeshType>
typename PlanarRegionFace<MeshType>::ScalarType
        PlanarRegionFace< MeshType>:: PlaneFittingError( std::vector<CoordType> &  samples,vcg::Plane3<float> p){

                typename PlanarRegionFace<MeshType>::ScalarType err =0.0;
                typename std::vector< CoordType>::iterator  si;
                for(si = samples.begin(); si != samples.end(); ++si)
                        err += fabs((*si)*p.Direction()-p.Offset());
                return err/samples.size();
}

// evaluate the gain if the triangle is added to the region
template<class MeshType>
typename PlanarRegionFace<MeshType>::ScalarType PlanarRegionFace< MeshType>::Evaluate( ElemType & f){
//	return (f.N()-p.Direction()).SquaredNorm()*f.Q()*0.5;
	return (vcg::NormalizedNormal(f)-p.Direction()).SquaredNorm()*vcg::DoubleArea(f)*0.5;
//	return vcg::Distance(vcg::Barycenter(f),p)*vcg::DoubleArea(f)*0.5;
}


	// refit the planes
template<class MeshType>
void PlanarRegionFace< MeshType>::Fit(std::vector<ElemType*>& faces, vcg::Plane3f & plane, float & err){
	FaceIterator fi;
	center = CoordType(0.0,0.0,0.0);
	std::vector<CoordType>  samples;

	if(faces.size()<3){
			samples.push_back((*faces.begin())->V(0)->P());
			samples.push_back((*faces.begin())->V(1)->P());
			samples.push_back((*faces.begin())->V(2)->P());
			center+=vcg::Barycenter(*(*faces.begin()))*3;
	}else
	for( fi = faces.begin(); fi != faces.end(); ++fi)
	{
		
		//AddSamples(samples,*(*fi));
		samples.push_back(vcg::Barycenter(*(*fi)));
		center+=vcg::Barycenter(*(*fi));
	}
	center*=1/(typename MeshType::ScalarType)(samples.size());

	if(samples.size()==3){
		plane.SetDirection(vcg::Normal(vcg::Triangle3<typename PlanarRegionFace<MeshType>::ScalarType >(samples[0],samples[1],samples[2])));
		typename MeshType::ScalarType off=samples[0]*plane.Direction();
		plane.SetOffset(off);
	}else
	{ 	vcg::PlaneFittingPoints<typename MeshType::ScalarType >(samples,plane);
	}

	err = PlaneFittingError(samples,plane);
	}

	template<class MeshType>
	void PlanarRegionFace< MeshType>::Refit(){
		Fit(elems,p,this->approx_err);
	}

	template<class MeshType>
	void PlanarRegionFace< MeshType>::Init(ElemType * b){
		p.Init(vcg::Barycenter(*b),b->N());
	}


	} // namespace tri
}// namespace vcg


#endif
