#ifndef __VCG_TETRA_DECIMATION_OPERATION_DEFS
# define __VCG_TETRA_DECIMATION_OPERATION_DEFS

namespace vcg{
  namespace tetra{
  

enum ModifiersType{MTEdgeCollapse,MTEdgeSplit};

  /** \addtogroup tetramesh */
/*@{*/
/// This Class is used to generalize a modifier

template <class TETRA_MESH_TYPE> 
class LocalModification
{

 public:
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
  ///the pos type
  typedef typename Pos<TetraType> PosType;
  typedef typename std::pair<vcg::tetra::ModifiersType,PosType> HeapRetElem;
  //the return type of heap updating 
  typedef typename std::vector<HeapRetElem> HeapRetType;
  /// Default Constructor
	LocalModification()
		{
		};
  
  ~LocalModification()
		{
		};
  
  
  virtual ScalarType ComputePriority()=0;
  virtual ScalarType ComputeError()=0;
  virtual void Execute()=0;
  virtual bool PreserveTopology()=0;
  virtual HeapRetType UpdateHeap()=0;

};//end class local modification

  
  }
}
#endif