/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.26  2005/02/22 10:38:14  ponchio
Debug, cleaning and optimization.

Revision 1.25  2005/02/21 17:55:47  ponchio
debug debug debug

Revision 1.24  2005/02/20 19:49:44  ponchio
cleaning (a bit more).

Revision 1.23  2005/02/20 18:07:01  ponchio
cleaning.

Revision 1.22  2005/02/20 00:43:24  ponchio
Less memory x extraction.  (removed frags)

Revision 1.21  2005/02/19 17:14:02  ponchio
History quick by default.

Revision 1.20  2005/02/19 10:45:04  ponchio
Patch generalized and small fixes.

Revision 1.19  2005/02/18 13:04:13  ponchio
Added patch reordering.

Revision 1.18  2005/02/08 12:43:03  ponchio
Added copyright


****************************************************************************/

#ifdef WIN32
#include <wrap/system/getopt.h>
#else
#include <unistd.h>
#endif

#include <assert.h>


#include <iostream>


#include <vcg/simplex/vertex/with/vc.h>
#include <vcg/simplex/face/face.h>
#include <vcg/complex/trimesh/base.h>
//WARNING WARNING this must be included AFTER mesh includes....
#include <wrap/io_trimesh/import_ply.h>
#include <vcg/space/index/grid_static_ptr.h>

#include "nxsalgo.h"
#include "strip.h"
#include "nexus.h"
#include "watch.h"

using namespace nxs;
using namespace vcg;
using namespace std;
using namespace tri;

class CFace;

class CVertex: public VertexVCf<DUMMYEDGETYPE,CFace ,DUMMYTETRATYPE> {};

class CFace: public Face<CVertex, DUMMYEDGETYPE , CFace>{};

class CMesh: public tri::TriMesh<vector<CVertex>, vector<CFace> > {};

string getSuffix(Signature &signature) {
  string suff;
  if(signature.compr)                     suff += "Z";
  if(signature.face == Signature::STRIPS) suff += "S"; 
  if(signature.vcolor)                    suff += "C";
  if(signature.vnorm)                     suff += "N";
  if(signature.vtext)                     suff += "T";
  if(signature.vdata)                     suff += "D";
  return suff;
}

void printInfo(Nexus &nexus, bool verbose, bool dump_history);

int main(int argc, char *argv[]) {
  
  string input;
  string output;
  string plysource;

  bool info = false;
  bool verbose = false;
  bool dump_history = false;

  unsigned int ram_size = 128000000;
  unsigned int chunk_size = 0;

  bool add = false;
  bool add_strips = false;
  bool add_colors = false;
  unsigned char add_normals = 0;
  bool add_textures = false;
  bool add_data = false;

  bool remove = false;
  bool remove_strips = false;
  bool remove_colors = false;
  bool remove_normals = false;
  bool remove_textures = false;
  bool remove_data = false;

  bool compress = false;
  bool uncompress = false;
  bool zsort = false;

  float qvertex = 0;
  float qnormal = 0;
  float qcolor = 0;
  float qtexture = 0;
  float cone_threshold = 0;

  int option;
  while((option = getopt(argc, argv, "ilho:a:r:zxsv:n:k:t:b:c:V:")) != EOF) {
    switch(option) {
    case 'i': info = true; break;
    case 'l': verbose = true; break;
    case 'h': dump_history = true; break;
    case 'o': output = optarg; break;
    case 'p': plysource = optarg; break;

    case 'a': {
      if(strstr(optarg, "strips"))   { add_strips = true;   add = true; }
      if(strstr(optarg, "colors"))   { add_colors = true;   add = true; }
      if(strstr(optarg, "normals"))  { 
	add_normals = Encodings::SHORT4;  add = true; }
      if(strstr(optarg, "normalf"))  { 
	add_normals = Encodings::FLOAT3;  add = true; }
      if(strstr(optarg, "textures")) { add_textures = true; add = true; }
      if(strstr(optarg, "data"))     { add_data = true;     add = true; }
      if(add == false) {
	cerr << "Invalid -a argument: " << optarg << "\n"
	     << "Valid options are: strips, colors, normals, textures, data\n";
	return -1;
      }
      break;
    }
    case 'r': {
      if(strstr(optarg, "strips")) {
	cerr << "Strips removing not supported!\n";
	return -1;
      }
      if(strstr(optarg, "colors"))   { remove_colors = true;   remove = true; }
      if(strstr(optarg, "normals"))  { remove_normals = true;  remove = true; }
      if(strstr(optarg, "textures")) { remove_textures = true; remove = true; }
      if(strstr(optarg, "data"))     { remove_data = true;     remove = true; }
      if(remove == false) {
	cerr << "Invalid -a argument: " << optarg << "\n"
	     << "Valid options are: strip, colors, normals, normalf, "
	     << "textures, data\n";
	return -1;
      }
      break;
    }
      


    case 'z': compress = true; break;
    case 'x': uncompress = true; break;
    case 's': zsort = true; break;

    case 'V': cone_threshold = atof(optarg); break;
    case 'v': qvertex = (float)atof(optarg); 
      if(qvertex == 0) {
	cerr << "Invalid value for quantization: " << optarg << endl;
	return -1;
      }
      break;
    case 'n': qnormal = (float)atof(optarg); 
      if(qnormal == 0) {
	cerr << "Invalid value for quantization: " << optarg << endl;
	return -1;
      }
      break;
    case 'k': qcolor = (float)atof(optarg); 
      if(qcolor == 0) {
	cerr << "Invalid value for quantization: " << optarg << endl;
	return -1;
      }
      break;
    case 't': qtexture = (float)atof(optarg); 
      if(qtexture == 0) {
	cerr << "Invalid value for quantization: " << optarg << endl;
	return -1;
      }
      break;

    case 'b': ram_size = atoi(optarg); 
      if(ram_size == 0) {
	cerr << "Invalid ram_size: " << optarg << endl;
	return -1;
      }
      break;
    case 'c': chunk_size = atoi(optarg); 
      if(chunk_size == 0) {
	cerr << "Invalid chunk_size: " << optarg << endl;
	return -1;
      }
      break;
    default: cerr << "Unknown option: " << (char)option << endl;
      return -1;
    }
  }

  if(optind != argc - 1) {
    cerr << "Usage: " << argv[0] << " <nexus file> [options]\n"
	 << " -i       : display some info about nexus file\n"
	 << " -l       : list nodes\n"
         << " -h       : list history\n"
         << " -o <file>: output filename (default is adding 00 to nexus)\n"
	 << " -a <what>: Add [colors|normals|normalf|strips|textures|data|borders]\n"
         << " -r <what>: As add...\n"
	 << " -p <ply> : Ply source for colors or textures or data\n"
	 << " -z       : compress\n"
	 << " -x       : uncompress\n"
	 << " -s       : sort using zcurve\n"
	 << " -v<float>: Vertex quantization (float is the 0 level amount)\n"
	 << " -n<float>: Normal quantization\n"
         << " -c<float>: Color quantization\n"
	 << " -t<float>: Texture quantization\n"
	 << " -V<float>: Normal cone threshold [0, 1] (0.95 default)\n\n"
	 << "            This option will not create a new nexus file\n";
    return -1;
  }
  input = argv[optind];

  //Sanity test of options...

  if(compress && uncompress) {
    cerr << "x and z are obviously exclusive :P\n";
    return -1;
  }

  if(add_normals && compress) {
    cerr << "Its not possible to add normals and compress in the same step\n"
	 << "Because normals requires 2 passes to be calculated\n\n";
    return -1;
  }
  


  bool compute_cones = false;
  if(!add && !remove && !compress && !uncompress && !zsort &&
     !qvertex && !qcolor && !qnormal && !qtexture && cone_threshold != 0)
    compute_cones = true;
  
  
  Nexus nexus;
  
  if(!nexus.Load(input, true)) {
    cerr << "Could not open nexus file: " << input << "\n";
    return -1;
  }
  nexus.MaxRam() = ram_size / nexus.chunk_size;


  //Sanity tests
  if(remove_strips && !(nexus.signature.face != Signature::STRIPS)) {
    cerr << "Nexus file does not have strips\n";
    return -1;
  }
  if(remove_colors && !nexus.signature.vcolor) {
    cerr << "Nexus file does not have colors\n";
    return -1;
  }
  if(remove_normals && !nexus.signature.vnorm) {
    cerr << "Nexus file does not have normals\n";
    return -1;
  }
  if(remove_textures && !nexus.signature.vtext) {
    cerr << "Nexus file does not have textures\n";
    return -1;
  }
  if(remove_data && !nexus.signature.vdata) {
    cerr << "Nexus file does not have data\n";
    return -1;
  }

  if(add_strips && (nexus.signature.face == Signature::STRIPS)) {
    cerr << "Nexus file already has strips\n";
    return -1;
  }
  if(add_colors && nexus.signature.vcolor) {
    cerr << "Nexus file already has colors\n";
    return -1;
  }
  if(add_normals && nexus.signature.vnorm) {
    cerr << "Nexus file already has normals\n";
    return -1;
  }
  if(add_textures && nexus.signature.vtext) {
    cerr << "Nexus file already has textures\n";
    return -1;
  }
  if(add_data && nexus.signature.vdata) {
    cerr << "Nexus file already has data\n";
    return -1;
  }

  if(nexus.IsCompressed() && compress) {
    cerr << "File already compressed.\n";
    return -1;
  }

  if(!nexus.IsCompressed() && uncompress) {
    cerr << "File not compressed.\n";
    return -1;
  }


  if(info) {
    cout << "Nexus file: " << input << "\n";
    printInfo(nexus, verbose, dump_history);
  }

  
  //determine if we must proceed:
  if(!add && !remove && !compress && !uncompress && !zsort &&
     qvertex == 0 && qnormal == 0 && qcolor == 0 && qtexture == 0 &&
     cone_threshold == 0) {
    nexus.Close();
    return 0;
  }
  
  if(compute_cones) {//just recalculate normal cones
    cerr << "Unimplemented at the moment...\n";

    /*vector<NCone3s> cones;
    //    ComputeCones(Nexus &nexus, float cone_threshold);
    nexus.Close();
    nexus.Load(intput, false);
    for(unsigned int i = 0; i < nexus.size(); i++) {
      nexus[i].cone = cones[i];
      }*/
    nexus.Close();
    return 0;
  }
  
  CMesh mesh;
  GridStaticPtr<CMesh::FaceContainer> grid;
  if(add_colors) {
    if(!plysource.size()) {
      cerr << "No plysource specified when adding color (-p option)\n";
    } else {
      if(!tri::io::ImporterPLY<CMesh>::Open(mesh, plysource.c_str())) {
	cerr << "Could not load ply: " << plysource << endl;
	return -1;
      }
      //calcoliamo il box:
      Box3f box;
      for(unsigned int i = 0; i < mesh.vert.size(); i++)
	box.Add(mesh.vert[i].P());
      grid.SetBBox(box);
      grid.Set(mesh.face);
    }
  }
  
  Signature signature = nexus.signature;
  if(add_strips) signature.face = Signature::STRIPS;
  if(add_normals) signature.vnorm = add_normals;
  if(add_colors) signature.vcolor = Encodings::BYTE4;   

  if(remove_normals) signature.vnorm = 0;
  if(remove_colors) signature.vcolor = 0;

  if(compress) signature.compr = Signature::LZO;
  if(uncompress) signature.compr = 0;

  if(!output.size()) output = input + getSuffix(signature);
  if(output == input) {
    cerr << "Output and input files are the same.\n"
	 << "Use option -o <filename>\n"
	 << "You do not want to overwrite your data. Trust me.\n";
    return -1;
  }

  cout << "Writing to nexus: " << output << endl;

  Nexus out;
  
  if(!chunk_size)
    chunk_size = nexus.chunk_size;

  if(!out.Create(output, signature, chunk_size)) {
    cerr << "Could not open output: " << output << endl;
    return -1;
  }

  //TODO fix this broken interface (you should not care abou chunk_size
  out.MaxRam() = ram_size / out.chunk_size;
  //TODO set rambuffer low (or even direct access!)

  vector<unsigned int> forward;
  vector<unsigned int> backward;
  if(zsort) 
    ZSort(nexus, forward, backward);

  //Fixing history
  assert(nexus.history.IsQuick());

  unsigned int hsize;
  char *buffer = nexus.history.Save(hsize);
  out.history.Load(hsize, buffer);



  if(zsort) {
    assert(0);
    //TODO FIX THIS...
    /*    if(out.history.IsQuick()) {
      for(unsigned int i = 0; i < out.history.n_frags(); i++)
	out.history.frags[i] = backward[out.history.frags[i]];
    } else {
      for(unsigned int i = 0; i < out.history.updates.size(); i++) {
	History::Update &update = out.history.updates[i];
	for(unsigned int k = 0; k < update.created.size(); k++)
	  update.created[k] = backward[update.created[k]];

	for(unsigned int k = 0; k < update.erased.size(); k++)
	  update.erased[k] = backward[update.erased[k]];
      }
      }*/
  }

  Report report(nexus.size());
  cout << "Copying and allocating...\n";
  for(unsigned int p = 0; p < nexus.size(); p++) {
    unsigned int patch = p;
    report.Step(patch);

    if(zsort) patch = forward[patch];

    Entry &src_entry = nexus[patch];
    Patch &src_patch = nexus.GetPatch(patch);
    Border &src_border = nexus.GetBorder(patch);


    vector<unsigned short> strip;
    if(add_strips) {
      ComputeTriStrip(src_patch.nf, src_patch.FaceBegin(), strip);
      assert(strip.size() < 65000);
      out.AddPatch(src_entry.nvert, strip.size(), src_border.Capacity());
      if(verbose) {
	cerr << "tri: " << src_patch.nf << " strip: " << strip.size() 
	     << " ratio: " << (float)strip.size()/(float)src_patch.nf  
	     << endl;
      }
    } else
      out.AddPatch(src_entry.nvert, src_entry.nface, src_border.Capacity());


    Entry &dst_entry = out[p];
    Patch &dst_patch = out.GetPatch(p);

    //copy vertices: 
    assert(out.signature.vert == Signature::POINT3F);
    assert(out.signature.vert = nexus.signature.vert);
    memcpy(dst_patch.Vert3fBegin(), src_patch.Vert3fBegin(), 
	   src_patch.nv * sizeof(Point3f));

    if(qvertex && !add_normals) {
      float *ptr = (float *)dst_patch.Vert3fBegin();
      for(int i = 0; i < dst_patch.nv*3; i++) 
	ptr[i] =  qvertex * (int)(ptr[i]/qvertex);
    } 

    
    //now faces.
    if(add_strips) {
      assert(out.signature.face == Signature::STRIPS);
      memcpy(dst_patch.FaceBegin(), &*strip.begin(), 
	     strip.size() * sizeof(short));
    } else {
      assert(nexus.signature.face == out.signature.face);
      if(nexus.signature.face == Signature::STRIPS) {
	memcpy(dst_patch.FaceBegin(), src_patch.FaceBegin(), 
	       src_patch.nf * sizeof(unsigned short));
      } else if(nexus.signature.face == Signature::TRIANGLES) {
	memcpy(dst_patch.FaceBegin(), src_patch.FaceBegin(), 
	       src_patch.nf * sizeof(unsigned short) * 3);
      } else {
	assert(0);
      }
    }

    if(nexus.signature.vcolor) {
      if(nexus.signature.vcolor == out.signature.vcolor) {
	memcpy(dst_patch.VColorBegin(), src_patch.VColorBegin(), 
	       Patch::encodings[out.signature.vcolor].size(dst_patch.nv));
      } else {
	assert(0);
      }
    }

    if(nexus.signature.vnorm) {
      if(nexus.signature.vnorm == out.signature.vnorm) {
	memcpy(dst_patch.VNormBegin(), src_patch.VNormBegin(), 
	       Patch::encodings[out.signature.vnorm].size(dst_patch.nv));
      } else {
	assert(0);
      }
    }

    //copying entry information;
    dst_entry.sphere = src_entry.sphere;
    dst_entry.error = src_entry.error;
    dst_entry.cone = src_entry.cone;

    out.borders.ResizeBorder(p, src_border.Size());
    Border &dst_border = out.GetBorder(p);
    memcpy(dst_border.Start(), src_border.Start(), 
	   src_border.Size() * sizeof(Link));    
    //TODO test this
    if(zsort)
      for(unsigned i = 0; i < dst_border.Size(); i++)
	dst_border[i].end_patch = backward[dst_border[i].end_patch];

  }
  report.Finish();  

  //TODO this is ok only if we have faces still!
  if(add_normals) {
    cout << "Computing normals" << endl;
    ComputeNormals(out);
  }

  if(add_colors) {
    //    if(!plysource.size()) 
    //source of color:
    //    cerr << "Unsupported color\n";
    //    return -1;
  }

  if(qvertex && add_normals) { 
    report.Init(nexus.size());
    cout << "Quantizing vertices\n";
    for(unsigned int patch = 0; patch < nexus.size(); patch++) {
      report.Step(patch);
      Patch src_patch = nexus.GetPatch(patch);

      float *ptr = (float *)src_patch.Vert3fBegin();
      for(int i = 0; i < src_patch.nv*3; i++) 
	      ptr[i] =  qvertex * (int)(ptr[i]/qvertex);
    }
    report.Finish();
  }

  out.sphere = nexus.sphere;
  
  out.Close();
  nexus.Close();
  return 0;
}


void printInfo(Nexus &nexus, bool verbose, bool dump_history) {
  //perform locality statistics
  double meandist = 0;
  vcg::Sphere3f last = nexus[0].sphere;
  for(unsigned int i = 1; i < nexus.size(); i++) {
    vcg::Sphere3f &sphere = nexus[i].sphere;
    double dist = vcg::Distance(last.Center(), sphere.Center());
    meandist += dist;
    last = sphere;
  }
  meandist /= nexus.size() -1;
  cout << "\n\tCompressed: " << nexus.IsCompressed() 
       << "\n\tStripped  : " 
       << (int)(nexus.signature.face == Signature::STRIPS)
       << "\n\tColor     : " << (int)(nexus.signature.vcolor)
       << "\n\tNormal    : " << (int)(nexus.signature.vnorm)
       << "\n\tTexture   : " << (int)(nexus.signature.vtext)
       << "\n\tData      : " << (int)(nexus.signature.vdata)
       << "\n\n\tVertices: " << nexus.totvert 
       << "\tFaces       : " << nexus.totface
       << "\tPatches     : " << nexus.size() 
       << "\n\tSphere    : " 
       << nexus.sphere.Center()[0] << " "
       << nexus.sphere.Center()[1] << " "
       << nexus.sphere.Center()[2] << " R: "
       << nexus.sphere.Radius()
       << "\n\tAverage distance: " << meandist
       << "\n\tChunk size " << nexus.chunk_size << endl;
   
  if(dump_history) {
    if(nexus.history.IsQuick()) {
      cout << "Quick format\n";
      for(unsigned int i = 0; i < nexus.history.n_nodes(); i++) {
	cout << "Node: " << i << " out: ";
	History::Node node = nexus.history.nodes[i];
	for(History::Link *l = node.out_begin; l != node.out_end; l++) {
	  cout << ".";
	  for(unsigned int p = l->begin; p != l->end; p++) {
	    cout << p << " ";
	  }
	}
	cout << " in: ";
	for(History::Link *j = node.in_begin; j != node.in_end; j++) {
	  cout << ".";
	  for(unsigned int p = j->begin; p != j->end; p++) {
	    cout << p << " ";
	  }
	}
	cout << endl;
      }

    } else {
      cout << "Update format\n";
      for(unsigned int i = 0; i < nexus.history.updates.size(); i++) {
	History::Update &update = nexus.history.updates[i];
	cout << "Created: ";
	for(unsigned int k = 0; k < update.created.size(); k++) {
	  cout << update.created[k] << " ";
	}
	cout << "\nErased: ";
	for(unsigned int k = 0; k < update.erased.size(); k++) {
	  cout << update.erased[k] << " ";
	}
	cout << "\n\n";
      }
    }
  }

  if(verbose) {
    for(unsigned int i = 0; i < nexus.size(); i++) {
      Entry &entry = nexus[i];
      cout << i << " -> nv: " << entry.nvert << " nf: " << entry.nface 
	   << " error: " << entry.error 
	   << " disk_size: " << entry.disk_size 
	   << " start: " << entry.patch_start << endl;
      cout << " Cone: " << entry.cone.n[0] << " "
	   << entry.cone.n[1] << " "
	   << entry.cone.n[2] << " "
	   << entry.cone.n[3] << "\n";
	   
    }
    cout << endl;
  }
}
