#include <iostream>
#include "nexus.h"

using namespace std;
using namespace vcg;
using namespace nxs;

Nexus::Nexus(): index_file(NULL) {}
Nexus::~Nexus() {
  Close();
}

bool Nexus::Create(const string &file) {
  index_file = fopen((file + ".nxs").c_str(), "wb+");
  if(!index_file) {
    cerr << "Could not create file: " << file << ".nxs\n";
    return false;
  }

  totvert = 0;
  totface = 0;
  totchunks = 0;
  totlinks = 0;
  sphere = Sphere3f();

  //Important: chunk_size must be 1 so that i can use Region in VFile.
  if(!patches.Create(file + ".nxp", 1)) {
    cerr << "Could not create file: " << file << ".nxp" << endl;
    return false;
  }
  if(!borders.Create(file + ".nxb")) {
    cerr << "Could not create file: " << file << ".nxp" << endl;
    return false;
  }
  return true;
}

bool Nexus::Load(const string &file) {
  index_file = fopen((file + ".nxs").c_str(), "rb+");
  if(!index_file) return false;
  
  unsigned int readed;
  readed = fread(&totvert, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;
  readed = fread(&totface, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;
  readed = fread(&totchunks, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;
  readed = fread(&totlinks, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;
  readed = fread(&sphere, sizeof(Sphere3f), 1, index_file);
  if(!readed) return false;
  
  unsigned int size; //size of index
  readed = fread(&size, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;

  index.resize(size);
  readed = fread(&index[0], sizeof(Entry), size, index_file);
  if(readed != size) return false;

  if(!patches.Load(file + ".nxp")) return false;
  if(!borders.Load(file + ".nxb")) return false;
  return true;
}

void Nexus::Close() {
  if(!index_file) return;
  rewind(index_file);
  
  fwrite(&totvert, sizeof(unsigned int), 1, index_file);
  fwrite(&totface, sizeof(unsigned int), 1, index_file);
  fwrite(&totchunks, sizeof(unsigned int), 1, index_file);
  fwrite(&totlinks, sizeof(unsigned int), 1, index_file);
  fwrite(&sphere, sizeof(Sphere3f), 1, index_file);
  
  unsigned int size = index.size(); //size of index
  fwrite(&size, sizeof(unsigned int), 1, index_file);
  fwrite(&index[0], sizeof(Entry), size, index_file);
  fclose(index_file);
  patches.Close();
  borders.Close();
}

Patch Nexus::GetPatch(unsigned int patch) {
  Entry &entry = index[patch];
  Chunk *start = patches.GetRegion(entry.patch_offset, entry.patch_size);
  return Patch(start);
}
