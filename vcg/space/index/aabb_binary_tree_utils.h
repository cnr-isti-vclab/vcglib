#ifndef __VCGLIB_AABBBINARYTREEUTILS
#define __VCGLIB_AABBBINARYTREEUTILS

// vcg headers
#include <vcg/space/point3.h>
#include <vcg/space/box3.h>


/***************************************************************************************/

namespace vcg {

/*

Class AABBBinaryTreeUtils

SAMPLE USAGE:

NOTES:

*/

template <class SCALARTYPE>
class AABBBinaryTreeUtils {
	public:
		typedef AABBBinaryTreeUtils<SCALARTYPE> ClassType;
		typedef SCALARTYPE ScalarType;

		template <class S, class T>
		static inline Point3<S> ConvertP3(const Point3<T> & p) {
			return (Point3<S>(S(p[0]), S(p[1]), S(p[2])));
		}

		class EmptyClass {
		};

		class ObjIteratorPtrFunct {
			public:
				template <class T>
				inline T * operator () (T & t) {
					return (&t);
				}

				template <class T *>
				inline T * operator () (T * & t) {
					return (t);
				}
		};

		template <class FACETYPE>
		class FaceBoxFunct {
			public:
				typedef FACETYPE FaceType;

				inline Box3<ScalarType> operator () (const FaceType & f) {
					Box3<ScalarType> box;
					box.SetNull();

					box.Add(AABBBinaryTreeUtils::ConvertP3<ScalarType, FaceType::ScalarType>(f.P(0)));
					box.Add(AABBBinaryTreeUtils::ConvertP3<ScalarType, FaceType::ScalarType>(f.P(1)));
					box.Add(AABBBinaryTreeUtils::ConvertP3<ScalarType, FaceType::ScalarType>(f.P(2)));

					return (box);
				}
		};

		template <class OBJTYPE>
		class ObjBarycenterFunct {
			public:
				typedef OBJTYPE ObjType;

				inline Point3<ScalarType> operator () (const ObjType & obj) {
					return (AABBBinaryTreeUtils::ConvertP3<ScalarType, ObjType::ScalarType>(obj.Barycenter()));
				}
		};

		template <class FACETYPE>
		class FacePointDistanceFunct {
			public:
				typedef FACETYPE FaceType;

				inline bool operator () (const FaceType & f, const Point3<ScalarType> & p, ScalarType & minDist, Point3<ScalarType> & res) {
					FaceType::ScalarType fdist;
					const AFace::CoordType fp = ConvertP3<AFace::ScalarType, ScalarType>(p);
					FaceType::CoordType fres;

					const bool br = face::PointDistance(f, fp, fdist, fres);
					minDist = (ScalarType)(fdist);
					res = ConvertP3<ScalarType, AFace::ScalarType>(fres);

					return (br);
				}
		};

};

/***************************************************************************************/

}	// end namespace vcg

#endif // #ifndef __VCGLIB_AABBBINARYTREEUTILS
