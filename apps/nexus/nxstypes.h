#ifndef NXS_TYPES_H
#define NXS_TYPES_H

#ifdef WIN32
typedef __int64 int64;
#else
typedef unsigned long long int64;
//typedef unsigned long long uint64;
#endif

typedef int int32;
typedef unsigned int uint32;

#endif
