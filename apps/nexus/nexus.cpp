#include <iostream>
#include <assert.h>
#include "nexus.h"

using namespace std;
using namespace vcg;
using namespace nxs;

Nexus::Nexus(): index_file(NULL) {}
Nexus::~Nexus() {
  Close();
}

bool Nexus::Create(const string &file, Signature sig, unsigned int c_size) {
  index_file = fopen((file + ".nxs").c_str(), "wb+");
  if(!index_file) {
    cerr << "Could not create file: " << file << ".nxs\n";
    return false;
  }

  signature = sig;
  totvert = 0;
  totface = 0;
  sphere = Sphere3f();
  
  index.clear();

  chunk_size = c_size;

  if(!patches.Create(file + ".nxp", signature, chunk_size)) {
    cerr << "Could not create file: " << file << ".nxp" << endl;
    return false;
  }
  //Important: chunk_size must be 1 so that i can use Region in VFile.
  if(!borders.Create(file + ".nxb", 1)) {
    cerr << "Could not create file: " << file << ".nxb" << endl;
    return false;
  }
  history.clear();
  return true;
}

bool Nexus::Load(const string &file, bool readonly) {
  index_file = fopen((file + ".nxs").c_str(), "rb+");
  if(!index_file) return false;
  
  unsigned int readed;
  readed = fread(&signature, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;
  readed = fread(&totvert, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;
  readed = fread(&totface, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;
  readed = fread(&sphere, sizeof(Sphere3f), 1, index_file);
  if(!readed) return false;
  readed = fread(&chunk_size, sizeof(unsigned int), 1, index_file);
    if(!readed) return false;

  unsigned int size; //size of index
  readed = fread(&size, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;

  index.resize(size);
  readed = fread(&index[0], sizeof(PatchInfo), size, index_file);
  if(readed != size) return false;
  patches.ReadEntries(index_file);
  borders.ReadEntries(index_file);
  
  //history size;
  fread(&size, sizeof(unsigned int), 1, index_file);
  vector<unsigned int> buffer;
  buffer.resize(size);
  fread(&(buffer[0]), sizeof(unsigned int), size, index_file);

  //number of history updates
  size = buffer[0];
  history.resize(size);

  unsigned int pos = 1;
  for(unsigned int i = 0; i < size; i++) {
    unsigned int erased = buffer[pos++];
    unsigned int created = buffer[pos++];
    history[i].erased.resize(erased);
    history[i].created.resize(created);
    for(unsigned int e = 0; e < erased; e++) 
      history[i].erased[e] = buffer[pos++];
    for(unsigned int e = 0; e < created; e++) 
      history[i].created[e] = buffer[pos++];
  }
  
  //TODO support readonly
  if(!patches.Load(file + ".nxp", signature, chunk_size, readonly)) 
    return false;
  if(!borders.Load(file + ".nxb", 1)) return false;
  return true;
}

void Nexus::Close() {  
  patches.Close();
  borders.Close();

  if(!index_file) return;
  rewind(index_file);

  fwrite(&signature, sizeof(unsigned int), 1, index_file);  
  fwrite(&totvert, sizeof(unsigned int), 1, index_file);
  fwrite(&totface, sizeof(unsigned int), 1, index_file);
  fwrite(&sphere, sizeof(Sphere3f), 1, index_file);
  fwrite(&chunk_size, sizeof(unsigned int), 1, index_file);
  
  unsigned int size = index.size(); //size of index
  fwrite(&size, sizeof(unsigned int), 1, index_file);
  fwrite(&(index[0]), sizeof(PatchInfo), size, index_file);

  patches.WriteEntries(index_file);
  borders.WriteEntries(index_file);

  vector<unsigned int> buffer;
  buffer.push_back(history.size());
  for(unsigned int i = 0; i < history.size(); i++) {
    Update &update = history[i];
    buffer.push_back(update.erased.size());
    buffer.push_back(update.created.size());
    for(unsigned int e = 0; e < update.erased.size(); e++)
      buffer.push_back(update.erased[e]);
    for(unsigned int e = 0; e < update.created.size(); e++)
      buffer.push_back(update.created[e]);
  }

  size = buffer.size();
  fwrite(&size, sizeof(unsigned int), 1, index_file);
  fwrite(&(buffer[0]), sizeof(unsigned int), size, index_file);


  fclose(index_file);
  index_file = NULL;
}

Patch &Nexus::GetPatch(unsigned int patch, bool flush) {
  assert(patch < index.size());
  PatchInfo &info = index[patch];
  return patches.GetPatch(patch, info.nvert, info.nface, flush);
}

Border Nexus::GetBorder(unsigned int patch, bool flush) {
  PatchInfo &info = index[patch];
  return borders.GetBorder(patch);
}

void Nexus::AddBorder(unsigned int patch, Link &link) {
  Border border = GetBorder(patch);
	
  unsigned int pos = border.Size();
  if(borders.ResizeBorder(patch, pos+1)) {
    border = GetBorder(patch);
  }
  
  assert(border.Size() < border.Available());
  assert(border.Available() > pos);

  border[pos] = link;  
}

unsigned int Nexus::AddPatch(unsigned int nvert, unsigned int nface,
			     unsigned int nbord) {

  PatchInfo info;
  info.nvert = nvert;
  info.nface = nface;
  
  patches.AddPatch(nvert, nface);
  borders.AddBorder(nbord);

  index.push_back(info);
  totvert += nvert;
  totface += nface;
  return index.size() -1;
}

void Nexus::SetRamBufferSize(unsigned int r_size) {    
  patches.SetRamBufferSize(r_size);
}

/*void Nexus::Join(const std::set<unsigned int> &patches,
		 std::vector<Point3f> &newvert,
		 std::vector<unsigned int> &newface,
		 std::vector<Link> &newbord) {

  map<unsigned int, vector<unsigned int> > remap;
  set<Link> newborders;
  set<unsigned int> erased;
  for(unsigned int u = 0; u < history.size(); u++) 
    for(unsigned int e = 0; e < history[u].erased.size(); e++)
      erased.insert(history[u].erased[e]);

  set<unsigned int>::const_iterator patch_idx;
  for(patch_idx = patches.begin(); patch_idx != patches.end(); patch_idx++) {
    unsigned int patch = *patch_idx;
    PatchInfo &entry = index[patch];
    remap[patch].resize(entry.nvert, 0xffffffff);
  }

  unsigned int vcount = 0;
  unsigned int fcount = 0;
  unsigned int bcount = 0;
  for(patch_idx = patches.begin(); patch_idx != patches.end(); patch_idx++) {
    unsigned int patch = *patch_idx;
    vector<unsigned int> &vmap = remap[patch];

    PatchInfo &entry = index[patch];
    fcount += entry.nface;
    for(unsigned int k = 0; k < entry.nvert; k++) {
      if(vmap[k] == 0xffffffff) { //first time
	      vmap[k] = vcount++;
      }
    }

    Border border = GetBorder(patch);
    for(unsigned int k = 0; k < border.Size(); k++) {
      Link link = border[k];
      if(link.IsNull()) continue;

      assert(link.start_vert < entry.nvert);
      assert(vmap[link.start_vert] != 0xffffffff);

      if(!remap.count(link.end_patch)) { //external
	      //test if erased in history... in wich case we do not add border
	      if(!erased.count(link.end_patch)) {
	        link.start_vert = vmap[link.start_vert];
	        newborders.insert(link);
	      }
	      continue; 
      } 
      //internal
      //TODO unsigned int &rmpv = remap[link.end_patch][link.end_vert];
      if(remap[link.end_patch][link.end_vert] == 0xffffffff) { //first time
	      remap[link.end_patch][link.end_vert] = vmap[link.start_vert];
      }
    }
  }



  newvert.resize(vcount);
  newface.resize(fcount*3);
  newbord.resize(0);

  fcount = 0;
  for(patch_idx = patches.begin(); patch_idx != patches.end(); patch_idx++) {
    Patch patch = GetPatch(*patch_idx);    

    vector<unsigned int> &vmap = remap[*patch_idx];
    assert(vmap.size() == patch.nv);

    for(unsigned int i = 0; i < vmap.size(); i++) {            
      assert(vmap[i] < vcount);
      newvert[vmap[i]] = patch.Vert(i);
    }
        
    for(unsigned int i = 0; i < patch.nf; i++) 
      for(int k = 0; k < 3; k++) {
        //TODO remove this check.
        if(patch.Face(i)[k] >= vmap.size()) {
          cerr << "Face overflow: " << patch.Face(i)[k] << endl;
          exit(0);
        }
        assert(patch.Face(i)[k] < vmap.size());        
	      newface[fcount++] = vmap[patch.Face(i)[k]];
      }
  }  
  
  assert(fcount == newface.size());

  set<Link>::iterator b;
  for(b = newborders.begin(); b != newborders.end(); b++)
    newbord.push_back(*b);

  return;
  }*/


void Nexus::Unify(float threshold) {
  //TODO what if colors or normals or strips?
  //TODO update totvert
  unsigned int duplicated = 0;
  unsigned int degenerate = 0;

  for(unsigned int p = 0; p < index.size(); p++) {
    PatchInfo &entry = index[p];
    Patch patch = GetPatch(p);

    unsigned int vcount = 0;
    map<Point3f, unsigned short> vertices;

    vector<unsigned short> remap;
    remap.resize(patch.nv);

    for(unsigned int i = 0; i < patch.nv; i++) {
      Point3f &point = patch.Vert(i);

      if(!vertices.count(point)) 
        vertices[point] = vcount++;
      else 
        duplicated++;

      remap[i] = vertices[point];
    }
    assert(vertices.size() <= patch.nv);
    if(vertices.size() == patch.nv) //no need to unify
      continue;

    vector<Point3f> newvert;
    newvert.resize(vertices.size());
    map<Point3f, unsigned short>::iterator k;
    for(k = vertices.begin(); k != vertices.end(); k++) {
      newvert[(*k).second] = (*k).first;
    }


    vector<unsigned short> newface;
    //check no degenerate faces get created.
    for(unsigned int f = 0; f < entry.nface; f++) {
      unsigned short *face = patch.Face(f);
      if(face[0] != face[1] && face[1] != face[2] && face[0] != face[2] &&
	      newvert[remap[face[0]]] != newvert[remap[face[1]]] &&
	      newvert[remap[face[0]]] != newvert[remap[face[2]]] &&
	      newvert[remap[face[1]]] != newvert[remap[face[2]]]) {
	      newface.push_back(remap[face[0]]);
	      newface.push_back(remap[face[1]]);
	      newface.push_back(remap[face[2]]);
      } else {
	      degenerate++;
      }
    }

    //rewrite patch now.
    entry.nvert = newvert.size();
    entry.nface = newface.size()/3;
    patch.Init(signature, entry.nvert, entry.nface);

    memcpy(patch.VertBegin(), &(newvert[0]), entry.nvert*sizeof(Point3f));
    memcpy(patch.FaceBegin(), &(newface[0]), entry.nface*3*sizeof(short));

    //testiamo il tutto...  TODO remove this of course
    for(unsigned int i =0; i < patch.nf; i++) {
      for(int k =0 ; k < 3; k++)
        if(patch.Face(i)[k] >= patch.nv) {
          cerr <<" Unify has problems\n";
          exit(0);
        }
    }
    
    //fix patch borders now
    set<unsigned int> close; //bordering pathes
    Border border = GetBorder(p);
    for(unsigned int b = 0; b < border.Size(); b++) {
      if(border[b].IsNull()) continue;
      close.insert(border[b].end_patch);
      border[b].start_vert = remap[border[b].start_vert];
    }

    set<unsigned int>::iterator c;
    for(c = close.begin(); c != close.end(); c++) {
      Border bord = GetBorder(*c);
      for(unsigned int b = 0; b < bord.Size(); b++) {
        if(bord[b].IsNull()) continue;
        if(bord[b].end_patch == p) {
          bord[b].end_vert = remap[bord[b].end_vert];
        }
      }
    }
  }
  //TODO: better to compact directly borders than setting them null.
  //finally: there may be duplicated borders
  for(unsigned int p = 0; p < index.size(); p++) {
    Border border = GetBorder(p);
    set<Link> links;
    for(unsigned int b = 0; b < border.Size(); b++) {
      if(border[b].IsNull()) continue;
      if(links.count(border[b]))
        border[b] = Link();
      else
        links.insert(border[b]);
    }
  }
  
  totvert -= duplicated;
  if(duplicated)
    cerr << "Found " << duplicated << " duplicated vertices" << endl;
  if(degenerate)
    cerr << "Found " << degenerate << " degenerate face while unmifying\n";
}
