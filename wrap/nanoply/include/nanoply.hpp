#ifndef NANOPLY_HPP
#define NANOPLY_HPP

#include <vector>
#include <tuple>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <cstdio>  // for fwrite()
#include <cmath>   // for fabs(),...
#include <limits>

// Avoid conflicting declaration of min/max macros in windows headers
#if !defined(NOMINMAX) && (defined(_WIN32) || defined(_WIN32_)  || defined(WIN32) || defined(_WIN64))
# define NOMINMAX
# ifdef max
#  undef   max
#  undef   min
# endif
#endif

namespace nanoply
{


typedef enum NNP_ERROR {
  NNP_OK                = 0x0000,
  NNP_UNABLE_TO_OPEN    = 0x0001,
  NNP_MISSING_HEADER    = 0x0002,
  NNP_MISSING_FORMAT    = 0x0004
} ErrorCode;

typedef enum NNP_ENTITY {
  NNP_UNKNOWN_ENTITY = 0x0000,
  NNP_PX		  = 0x0001,
  NNP_PY		  = 0x0002,
  NNP_PZ	      = 0x0004,
  NNP_PXYZ        = 0x0007,

  NNP_NX		  = 0x0010,
  NNP_NY		  = 0x0020,
  NNP_NZ		  = 0x0040,
  NNP_NXYZ        = 0x0070,

  NNP_CR		  = 0x0100,
  NNP_CG		  = 0x0200,
  NNP_CB		  = 0x0400,
  NNP_CA		  = 0x0800,

  NNP_RADIUS	  = 0x1000,
  NNP_SCALE		  = 0x2000
} PlyEntity;

typedef enum NNP_PLYTYPE {
  NNP_UNKNOWN_TYPE= 0x0000,
  NNP_FLOAT32     = 0x0001,
  NNP_FLOAT64     = 0x0002,
  NNP_INT8 	      = 0x0010,
  NNP_INT16       = 0x0020,
  NNP_INT32       = 0x0040,
  NNP_UINT8       = 0x0100,
  NNP_UINT16      = 0x0200,
  NNP_UINT32      = 0x0400,
  NNP_LIST_UINT8_UINT32 = 0x1000,
  NNP_LIST_INT8_UINT32  = 0x2000,
  NNP_LIST_UINT8_INT32  = 0x4000,
  NNP_LIST_INT8_INT32   = 0x8000
} PlyType;

class PlyProperty
{
public:
  PlyProperty(PlyType _t, PlyEntity _e):type(_t),elem(_e){}
  PlyType type;
  PlyEntity elem;
};

// In Ply an element is a collection of properties.
//
class PlyElement
{
public:
  std::string name;
  int cnt;
  std::vector<PlyProperty> propVec;
  bool AddProperty(const char *line);
  bool InitFromHeader(FILE *fp,  char *line);
};
//

bool PlyElement::InitFromHeader(FILE *fp, char *line)
{
  char buf[128];
  for(int i=0;line[i]!=0;++i) buf[i]=tolower(line[i]);
  strtok(buf," \t");    // property
  assert(strstr(buf,"element"));
  char *el = strtok(0," \t\n"); // float
  name = std::string(el);
  sscanf(line,"%*s %*s %i",&cnt);
  printf("Adding Element '%s' (%i) \n",name.c_str(),cnt);
  fgets(line,4096,fp);
  while(strstr(line,"property"))
  {
    AddProperty(line);
    fgets(line,4096,fp);
  }
}
// Skip all the elements in a file
bool PlyElement::SkipAsciiElementsInFile(FILE *fp)
{
  char line[4096];
  for(int i=0;i<this->cnt;++i)
    fgets(line,4096,fp);
}
// Skip all the elements in a file
bool PlyElement::SkipBinaryElementsInFile(FILE *fp)
{
  for(int i=0;i<this->cnt;++i)
  {
    for(j=0;j<this->propVec.size();++j)
    {

    }
  }
}

bool PlyElement::AddProperty(const char *line)
{
  char buf[128];
  for(int i=0;line[i]!=0;++i) buf[i]=tolower(line[i]);
  strtok(buf," \t");    // property
  char *ty = strtok(0," \t\n"); // float
  PlyType type = NNP_UNKNOWN_TYPE;
  PlyEntity prop = NNP_UNKNOWN_ENTITY;

  if(strstr(ty,"float")  || strstr(ty,"float32") ) type = NNP_FLOAT32;
  if(strstr(ty,"double") || strstr(ty,"float64") ) type = NNP_FLOAT64;
  if(strstr(ty,"char")   || strstr(ty,"int8") )    type = NNP_INT8;
  if(strstr(ty,"short")  || strstr(ty,"int16") )   type = NNP_INT16;
  if(strstr(ty,"int")    || strstr(ty,"int32") )   type = NNP_INT32;
  if(strstr(ty,"uchar")  || strstr(ty,"uint8") )    type = NNP_UINT8;
  if(strstr(ty,"ushort") || strstr(ty,"uint16") )   type = NNP_UINT16;
  if(strstr(ty,"uint")   || strstr(ty,"uint32") )   type = NNP_UINT32;
  if(strstr(ty,"list") )  {
    char *ty1 = strtok(0," \t\n");
    char *ty2 = strtok(0," \t\n");
    if( (strstr(ty1,"uchar")||strstr(ty1,"uint8") ) && (strstr(ty2,"uint")|| strstr(ty2,"uint32")) ) type = NNP_LIST_UINT8_UINT32;
    if( (strstr(ty1,"char") || strstr(ty1,"int8") ) && (strstr(ty2,"uint")|| strstr(ty2,"uint32")) ) type = NNP_LIST_INT8_UINT32;
    if( (strstr(ty1,"uchar")||strstr(ty1,"uint8") ) && (strstr(ty2,"int") || strstr(ty2,"int32")) ) type = NNP_LIST_UINT8_INT32;
    if( (strstr(ty1,"char") || strstr(ty1,"int8") ) && (strstr(ty2,"int") || strstr(ty2,"int32")) ) type = NNP_LIST_INT8_INT32;
  }

  assert(type != NNP_UNKNOWN_TYPE);
     char *el = strtok(0," \t\n"); // x

  if(strstr(el,"x")) prop = NNP_PX;
  if(strstr(el,"y")) prop = NNP_PY;
  if(strstr(el,"z")) prop = NNP_PZ;
  if(strstr(el,"nx")) prop = NNP_NX;
  if(strstr(el,"ny")) prop = NNP_NY;
  if(strstr(el,"nz")) prop = NNP_NZ;

//  assert(prop != NNP_UNKNOWN_PROP);
  propVec.push_back(PlyProperty(type,prop));
  printf("Adding Property %i %i\n",type,prop);
}

class Info
{
public:
  ErrorCode errInfo;
  bool binary;
  int mask;
  int vertexNum;
  std::vector<PlyElement> elemVec;
  void Clear() { errInfo=NNP_OK; }
};



inline bool GetInfo(const char *filename, nanoply::Info &info)
{
  FILE *fp=fopen(filename,"r");
  info.Clear();

  if(!fp)
  {
    info.errInfo = NNP_UNABLE_TO_OPEN;
    return false;
  }
  char buf[4096];
  fgets(buf,4096,fp);
  if( (strncmp(buf,"PLY",3)!=0) && (strncmp(buf,"ply",3)!=0) )
  {
    info.errInfo = NNP_MISSING_HEADER;
    return false;
  }
  fgets(buf,4096,fp);
  if(strncmp(buf,"format",strlen("format"))!=0)
    {
      info.errInfo = NNP_MISSING_FORMAT;
      return false;
    }

  if(strstr(buf,"ascii") || strstr(buf,"ASCII"))
    info.binary=false;
  else
  {
    if(strstr(buf,"binary") || strstr(buf,"BINARY"))
    info.binary=true;
    else
    {
      info.errInfo = NNP_MISSING_FORMAT;
      return false;
    }
  }
  fgets(buf,4096,fp);
  while(strncmp(buf,"end_header",strlen("end_header")))
  {
    if(strstr(buf,"comment") || strstr(buf,"COMMENT") )
    {
      fgets(buf,4096,fp);
      continue;
    }
    if(strstr(buf,"element"))
    {
      PlyElement pe;
      pe.InitFromHeader(fp,buf);
      info.elemVec.push_back(pe);
    }
  }
 }

template<class DestinationType>
class DataDescriptor
{
  DestinationType *base;
  int curPos;
  PlyProperty elem;
public:
  DataDescriptor(PlyProperty _e, DestinationType *_b):elem(_e),base(_b),curPos(0)
  {

  }

  DataDescriptor(PlyProperty _e, void *_b ):elem(_e),base((DestinationType *)_b),curPos(0)
  {

  }
  void Restart(){curPos=0;}
  void ReadElem(FILE *fp, int numOfElemToRead)
  {
    fread(base+curPos,sizeof(DestinationType),numOfElemToRead,fp);
  }
};



template <class DataAdaptorTuple>
bool Open(const char *filename, DataAdaptorTuple vertexAdaptor,  nanoply::Info &info)
{
  FILE *fp=fopen(filename,"rb");
  if(!fp)
  {
    return false;
  }
  char buf[4096];
  do
  {
    fgets(buf,4096,fp);
  }
  while(strncmp(buf,"end_header",strlen("end_header")) );
  // Now start the real reading!
  for(int i =0; i< info.elemVec.size();++i)
  {

  }
}

} // end namespace nanoply
#endif // NANOPLY_HPP
