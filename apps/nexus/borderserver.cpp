#include "borderserver.h"
#include <iostream>

using namespace std;
using namespace nxs;

void BorderServer::AddBorder(unsigned short nbord, unsigned int used) {
  BorderEntry entry;
  entry.border_start = Size();
  entry.border_size = nbord;
  entry.border_used = used;
  borders.push_back(entry);
  Resize(entry.border_start + nbord);
}

Border BorderServer::GetBorder(unsigned int border, bool flush) {
  assert(border < borders.size());
  BorderEntry &entry = borders[border];
  Link *start = GetRegion(entry.border_start, entry.border_size, flush);
  return Border(start, entry.border_used, entry.border_size);
}

bool BorderServer::ResizeBorder(unsigned int border, unsigned int nbord) {
  assert(nbord < 65500);
  assert(border < borders.size());
  BorderEntry &entry = borders[border];
  if(nbord > entry.border_size) {
    int capacity = nbord;
    if(capacity < entry.border_size*2) 
      capacity = entry.border_size * 2;
    if(capacity > 65500) 
      capacity = 65500;
    unsigned int newstart = Size(); 
    Resize(newstart + capacity);

    Link *src = GetRegion(entry.border_start, entry.border_size);
    Link *dst = GetRegion(newstart, capacity, false);
    memcpy(dst, src, entry.border_used * sizeof(Link));
    entry.border_start = newstart;
    entry.border_size = capacity;
    entry.border_used = nbord;
    return true;
  }
  entry.border_used = nbord;
  return false;
}

bool BorderServer::ReadEntries(FILE *fp) {
  unsigned int n;
  fread(&n, 1, sizeof(int), fp);
  borders.resize(n);
  fread(&*borders.begin(), n, sizeof(BorderEntry), fp);
  return true;
}

bool BorderServer::WriteEntries(FILE *fp) {
  unsigned int n = borders.size();
  fwrite(&n, 1, sizeof(int), fp);
  fwrite(&*borders.begin(), n, sizeof(BorderEntry), fp);
  return true;
}

  
