#include <assert.h>

#include <iostream>
#include <set>

#include "nexus.h"

using namespace std;
using namespace vcg;
using namespace nxs;

Nexus::~Nexus() {
  Close();
}

bool Nexus::Create(const string &file, Signature sig, unsigned int c_size) {
  signature = sig;
  totvert = 0;
  totface = 0;
  sphere = Sphere3f();
  chunk_size = c_size;
  unsigned int header_size = 256;
  if(chunk_size > header_size) header_size = chunk_size;
  
  history.Clear();
  ram_used = 0;
  ram_max = 50 * (1<<20) / chunk_size;

  if(!IndexFile<Entry>::Create(file + ".nxp", header_size)) {
    cerr << "Could not create file: " << file << ".nxp" << endl;
    return false;
  }

  //Important: chunk_size must be 1 so that i can use Region in VFile.
  if(!borders.Create(file + ".nxb")) {
    cerr << "Could not create file: " << file << ".nxb" << endl;
    return false;
  }
  return true;
}


bool Nexus::Load(const string &file, bool rdonly) {
  if(!IndexFile<Entry>::Load(file + ".nxp", rdonly)) return false;
  ram_used = 0;
  ram_max = 50 * (1<<20) / chunk_size;

  history.Clear();
  SetPosition(history_offset);
  unsigned int history_size;
  ReadBuffer(&history_size, sizeof(unsigned int));

  char *buffer = new char[history_size];
  ReadBuffer(buffer, history_size);

  if(!history.Load(history_size, buffer)) {
    cerr << "Error loading history\n";
    return false;
  }

  borders.Load(file + ".nxb", rdonly);
  //TODO on nxsbuilder assure borders are loaded
  return true;
}

void Nexus::Close() { 
  if(!Opened()) return;

  Flush();

  if(!IsReadOnly()) {
    //set history_offset
    if(!size()) history_offset = 0;
    else 
      history_offset = (back().patch_start + back().disk_size);
    history_offset *= chunk_size;

    unsigned int history_size;
    char *mem = history.Save(history_size);
    Redim(history_offset + history_size + sizeof(unsigned int));
    SetPosition(history_offset);
    WriteBuffer(&history_size, sizeof(unsigned int));
    WriteBuffer(mem, history_size);
    delete []mem;
  }
  borders.Close();
  IndexFile<Entry>::Close();
}

void Nexus::SaveHeader() {
  unsigned int magic = 0x3053584e; // nxs0
  WriteBuffer(&magic, sizeof(unsigned int));
  WriteBuffer(&signature, sizeof(unsigned int));
  WriteBuffer(&chunk_size, sizeof(unsigned int));
  WriteBuffer(&offset, sizeof(int64));
  WriteBuffer(&history_offset, sizeof(int64));
  WriteBuffer(&totvert, sizeof(unsigned int));
  WriteBuffer(&totface, sizeof(unsigned int));
  WriteBuffer(&sphere, sizeof(Sphere3f));
}

bool Nexus::LoadHeader() {
  unsigned int magic;
  ReadBuffer(&magic, sizeof(unsigned int));
  if(magic != 0x3053584e) {
    cerr << "Invalid magic. Not a nxs file\n";
    return false;
  }
  ReadBuffer(&signature, sizeof(unsigned int));
  ReadBuffer(&chunk_size, sizeof(unsigned int));
  ReadBuffer(&offset, sizeof(int64));
  ReadBuffer(&history_offset, sizeof(int64));
  ReadBuffer(&totvert, sizeof(unsigned int));
  ReadBuffer(&totface, sizeof(unsigned int));
  ReadBuffer(&sphere, sizeof(Sphere3f));
}

void Nexus::Flush(bool all) {
  if(all) {
    std::map<unsigned int, list<unsigned int>::iterator>::iterator i;
    for(i = index.begin(); i != index.end(); i++) {
      unsigned int patch = (*i).first;
      FlushPatch(patch);
    }
    pqueue.clear();
    index.clear();
  } else {
    while(ram_used > ram_max) {    
      unsigned int to_flush = pqueue.back();
      pqueue.pop_back();
      index.erase(to_flush);        
      FlushPatch(to_flush);
    }
  }
}



Patch &Nexus::GetPatch(unsigned int patch, bool flush) { 
  Entry &entry = operator[](patch);
  if(index.count(patch)) {
    assert(entry.patch);
    list<unsigned int>::iterator &i = index[patch];
    pqueue.erase(i);
    pqueue.push_front(patch);

  } else {
    while(flush && ram_used > ram_max) {    
      unsigned int to_flush = pqueue.back();
      pqueue.pop_back();
      index.erase(to_flush);        
      FlushPatch(to_flush);
    }
    assert(!entry.patch);
    entry.patch = LoadPatch(patch);
    pqueue.push_front(patch);
    list<unsigned int>::iterator i = pqueue.begin();
    index[patch] = i;      
  }                        
  return *(entry.patch);
}

Border Nexus::GetBorder(unsigned int patch, bool flush) {
  return borders.GetBorder(patch);
}

/*void Nexus::AddBorder(unsigned int patch, Link &link) {
  Border border = GetBorder(patch);
	
  unsigned int pos = border.Size();
  if(pos > 65500) {
    cerr << "Exceding border size!!!\n";
    exit(0);
  }
  if(borders.ResizeBorder(patch, pos+1)) {
    border = GetBorder(patch);
  }
  
  assert(border.Available() > pos);

  border[pos] = link;  
  }*/

unsigned int Nexus::AddPatch(unsigned int nvert, unsigned int nface,
			     unsigned int nbord) {

  Entry entry;
  entry.patch_start = 0xffffffff;
  entry.ram_size = Patch::ChunkSize(signature, nvert, nface, chunk_size);
  entry.disk_size = 0xffff;
  entry.nvert = nvert;
  entry.nface = nface;
  entry.error = 0;
  //sphere undefined.
  entry.patch = NULL;
  entry.vbo_array = 0;
  entry.vbo_element = 0;
  
  push_back(entry);
  
  borders.AddBorder(nbord);

  totvert += nvert;
  totface += nface;
  return size() - 1;
}

void Nexus::Unify(float threshold) {
  //TODO what if colors or normals or strips?
  unsigned int duplicated = 0;
  unsigned int degenerate = 0;

  for(unsigned int p = 0; p < size(); p++) {
    Entry &entry = operator[](p);
    Patch &patch = GetPatch(p);

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
  //better to compact directly borders than setting them null.
  //finally: there may be duplicated borders
  for(unsigned int p = 0; p < size(); p++) {
    Border border = GetBorder(p);
    set<Link> links;
    for(unsigned int b = 0; b < border.Size(); b++) {
      Link &link = border[b];
      assert(!link.IsNull());
      //if(border[b].IsNull()) continue;
      links.insert(link);
    }
    int count = 0;
    for(set<Link>::iterator k = links.begin(); k != links.end(); k++)
      border[count++] = *k;      
    
    borders[p].used = links.size();
  }
  
  totvert -= duplicated;
  if(duplicated)
    cerr << "Found " << duplicated << " duplicated vertices" << endl;
  if(degenerate)
    cerr << "Found " << degenerate << " degenerate face while unmifying\n";
}

Patch *Nexus::LoadPatch(unsigned int idx) {
  assert(idx < size());
  Entry &entry = operator[](idx);
  if(entry.patch) return entry.patch;
  
  char *ram = new char[entry.ram_size * chunk_size];
#ifndef NDEBUG
  if(!ram) {
    cerr << "COuld not allocate ram!\n";
    exit(0);
  }
#endif
  
  Patch *patch = new Patch(signature, ram, entry.nvert, entry.nface);
  
  if(entry.patch_start != 0xffffffff) { //was allocated.
    assert(entry.disk_size != 0xffff);
    
    MFile::SetPosition((int64)entry.patch_start * (int64)chunk_size);
    
    if((signature & NXS_COMPRESSED) == 0) { //not compressed
      MFile::ReadBuffer(ram, entry.disk_size * chunk_size);
    } else {
      unsigned char *disk = new unsigned char[entry.disk_size * chunk_size];
      MFile::ReadBuffer(disk, entry.disk_size * chunk_size);
      
      patch->Decompress(entry.ram_size * chunk_size, 
			disk, entry.disk_size * chunk_size);
      delete []disk;
    } 
  }
  ram_used += entry.ram_size;  
  entry.patch = patch;  
  return patch;
}

void Nexus::FlushPatch(unsigned int id) {
  Entry &entry = operator[](id);    
  assert(entry.patch);

  if(!MFile::IsReadOnly()) { //write back patch
    if((signature & NXS_COMPRESSED)) {
      unsigned int compressed_size;
      char *compressed = entry.patch->Compress(entry.ram_size * chunk_size,
					       compressed_size);
      if(entry.disk_size == 0xffff) {//allocate space 
	assert(entry.patch_start == 0xffffffff);
	entry.disk_size = (unsigned int)((compressed_size-1)/chunk_size) + 1;
	entry.patch_start = (unsigned int)(Length()/chunk_size);
	Redim(Length() + entry.disk_size * chunk_size);
      } else {
	//cerr << "OOOOPSPPPS not supported!" << endl;
	exit(-1);
      }
      MFile::SetPosition((int64)entry.patch_start * (int64)chunk_size);
      MFile::WriteBuffer(compressed, entry.disk_size * chunk_size);
      delete []compressed;
    } else {
      if(entry.disk_size == 0xffff) {
	entry.disk_size = entry.ram_size;
	entry.patch_start = (unsigned int)(Length()/chunk_size);
	Redim(Length() + entry.disk_size * chunk_size);
      }
      MFile::SetPosition((int64)entry.patch_start * (int64)chunk_size);
      MFile::WriteBuffer(entry.patch->start, entry.disk_size * chunk_size);
    }
  }

  delete [](entry.patch->start);
  delete entry.patch;  
  entry.patch = NULL;    
  ram_used -= entry.ram_size;      
}
