#ifndef __VCGLIB_AABBBINARYTREESHAREDUTILS
#define __VCGLIB_AABBBINARYTREESHAREDUTILS

// vcg headers
#include <vcg/math/base.h>

namespace vcg {

template <class P3TYPE>
inline P3TYPE Abs(const P3TYPE & p) {
	return (P3TYPE(math::Abs(p[0]), math::Abs(p[1]), math::Abs(p[2])));
}

template <class P3TYPE>
inline P3TYPE LowerClamp(const P3TYPE & p, const typename P3TYPE::ScalarType & r) {
	return (P3TYPE(math::Max<typename P3TYPE::ScalarType>(p[0], r), math::Max<typename P3TYPE::ScalarType>(p[1], r), math::Max<typename P3TYPE::ScalarType>(p[2], r)));
}

} // end namespace vcg

#endif // #ifndef __VCGLIB_AABBBINARYTREESHAREDUTILS
