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

bool Nexus::Create(const string &file, Signature sig) {
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

  //Important: chunk_size must be 1 so that i can use Region in VFile.
  if(!patches.Create(file + ".nxp", 1)) {
    cerr << "Could not create file: " << file << ".nxp" << endl;
    return false;
  }
  if(!borders.Create(file + ".nxb", 1)) {
    cerr << "Could not create file: " << file << ".nxb" << endl;
    return false;
  }
  history.clear();
  return true;
}

bool Nexus::Load(const string &file) {
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
  
  unsigned int size; //size of index
  readed = fread(&size, sizeof(unsigned int), 1, index_file);
  if(!readed) return false;

  index.resize(size);
  readed = fread(&index[0], sizeof(Entry), size, index_file);
  if(readed != size) return false;

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

  if(!patches.Load(file + ".nxp", 1)) return false;
  if(!borders.Load(file + ".nxb", 1)) return false;
  return true;
}

void Nexus::Close() {
  if(!index_file) return;
  rewind(index_file);

  fwrite(&signature, sizeof(unsigned int), 1, index_file);  
  fwrite(&totvert, sizeof(unsigned int), 1, index_file);
  fwrite(&totface, sizeof(unsigned int), 1, index_file);
  fwrite(&sphere, sizeof(Sphere3f), 1, index_file);
  
  unsigned int size = index.size(); //size of index
  fwrite(&size, sizeof(unsigned int), 1, index_file);
  fwrite(&(index[0]), sizeof(Entry), size, index_file);

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

  patches.Close();
  borders.Close();
}

Patch Nexus::GetPatch(unsigned int patch, bool flush) {
  Entry &entry = index[patch];
  Chunk *start = patches.GetRegion(entry.patch_start, entry.patch_used,flush);
  return Patch(signature, start, entry.nvert, entry.nface);
}

Border Nexus::GetBorder(unsigned int patch, bool flush) {
  Entry &entry = index[patch];
  Link *start = borders.GetRegion(entry.border_start, entry.border_size,flush);
  return Border(start, entry.border_used, entry.border_size);
}


unsigned int Nexus::AddPatch(unsigned int nvert, unsigned int nface,
			     unsigned int nbord) {
  Entry entry;
  entry.patch_start = patches.Size();
  entry.patch_size = Patch::ChunkSize(signature, nvert, nface);
  entry.patch_used = entry.patch_size;
  entry.border_start = borders.Size();
  entry.border_size = nbord;
  entry.border_used = 0;
  entry.nvert = nvert;
  entry.nface = nface;
  
  patches.Resize(patches.Size() + entry.patch_size);
  borders.Resize(borders.Size() + nbord);
  index.push_back(entry);
  totvert += nvert;
  totface += nface;
  return index.size() -1;
}

void Nexus::Join(const std::set<unsigned int> &patches,
		 std::vector<Point3f> &newvert,
		 std::vector<unsigned int> &newface,
		 std::vector<Link> &newbord) {

  map<unsigned int, vector<unsigned int> > remap;
  set<Link> newborders;
  set<unsigned int> erased;
  for(int u = 0; u < history.size(); u++) 
    for(int e = 0; e < history[u].erased.size(); e++)
      erased.insert(history[u].erased[e]);

  set<unsigned int>::const_iterator i;
  for(i = patches.begin(); i != patches.end(); i++) {
    unsigned int patch = *i;
    Nexus::Entry &entry = index[patch];
    remap[*i].resize(entry.nvert, 0xffffffff);
  }

  unsigned int vcount = 0;
  unsigned int fcount = 0;
  unsigned int bcount = 0;
  for(i = patches.begin(); i != patches.end(); i++) {
    unsigned int patch = *i;
    vector<unsigned int> &vmap = remap[*i];

    Nexus::Entry &entry = index[patch];
    fcount += entry.nface;
    //    assert(fcount < 0xffff);
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
      if(remap[link.end_patch][link.end_vert] == 0xffffffff) { //first time
	remap[link.end_patch][link.end_vert] = vmap[link.start_vert];
      }
    }
  }



  newvert.resize(vcount);
  newface.resize(fcount*3);
  newbord.resize(0);

  fcount = 0;
  for(i = patches.begin(); i != patches.end(); i++) {
    Patch patch = GetPatch(*i);
    //    Border border = GetBorder(*i);
    
    vector<unsigned int> &vmap = remap[*i];

    for(unsigned int i = 0; i < patch.nv; i++) {
      assert(vmap[i] < vcount);
      newvert[vmap[i]] = patch.Vert(i);
    }

    for(unsigned int i = 0; i < patch.nf; i++) {
      for(int k = 0; k < 3; k++) {
	newface[3*fcount + k] = vmap[patch.Face(i)[k]];
      }
      assert(patch.Face(i)[0] != patch.Face(i)[1]);
      assert(patch.Face(i)[0] != patch.Face(i)[2]);
      assert(patch.Face(i)[1] != patch.Face(i)[2]);
      assert(newface[3*fcount + 0] != newface[3*fcount + 1]);
      assert(newface[3*fcount + 0] != newface[3*fcount + 2]);
      assert(newface[3*fcount + 1] != newface[3*fcount + 2]);
      
      fcount++;
      assert(fcount *3 <= newface.size());
    }
    /*    for(unsigned int i = 0; i < border.Size(); i++) {
      Link link = border[i];
      if(patches.count(link.end_patch)) continue; 
      link.start_vert = vmap[link.start_vert];
      newbord.push_back(link);*/

      /*      if(remap.count(link.end_patch)) continue;
      Link newlink = link;
      newlink.start_vert = vmap[link.start_vert];
      newbord[bcount++] = newlink;*/
    //    }
  }  
  set<Link>::iterator b;
  for(b = newborders.begin(); b != newborders.end(); b++)
    newbord.push_back(*b);

 /* unsigned int newentry = AddPatch(vcount, fcount, bcount);
  Patch newpatch = GetPatch(newentry);
  Border newborder = GetBorder(newentry);

  memcpy(newpatch.VertBegin(), &(newvert)[0],
	 newvert.size() * sizeof(Point3f));

  memcpy(newpatch.FaceBegin(), &(newface)[0],
	 newface.size() * sizeof(unsigned short));

  memcpy(&(newborder[0]), &(newbord[0]),
	 newbord.size() * sizeof(Link));*/
  return;
}


void Nexus::Unify(float threshold) {
  //TODO what if colors or normals or strips?
  //TODO update totvert
  unsigned int duplicated = 0;
  for(unsigned int p = 0; p < index.size(); p++) {
    Nexus::Entry &entry = index[p];
    Patch patch = GetPatch(p);

    unsigned int vcount = 0;
    map<Point3f, unsigned short> vertices;
    vector<unsigned short> remap;
    remap.resize(patch.nv);
    //    map<unsigned short, unsigned short> remap;
    for(unsigned int i = 0; i < patch.nv; i++) {
      Point3f &point = patch.Vert(i);

      if(!vertices.count(point)) {
        vertices[point] = vcount++;
      } else {
        duplicated++;
      }

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
    newface.resize(patch.nf * 3);
    for(unsigned int f = 0; f < newface.size(); f++)
      newface[f] = remap[patch.FaceBegin()[f]];

    //rewrite patch now.
    entry.nvert = newvert.size();
    patch.Init(signature, entry.nvert, entry.nface);

    memcpy(patch.VertBegin(), &(newvert[0]), entry.nvert*sizeof(Point3f));
    memcpy(patch.FaceBegin(), &(newface[0]), entry.nface*3*sizeof(unsigned short));

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
  //finally: there may be duplicated borders
  for(unsigned int p = 0; p < index.size(); p++) {
    Border border = GetBorder(p);
    //Nexus::Entry &entry = index[p];

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
  //  cout << "Found " << duplicated << " duplicated vertices" << endl;

}
