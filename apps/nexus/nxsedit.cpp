#include <iostream>
using namespace std;

#include "nxsalgo.h"
#include "nexus.h"

using namespace nxs;
using namespace vcg;

#ifdef WIN32
#include "getopt.h"
#else
#include <unistd.h>
#endif

int main(int argc, char *argv[]) {
  
  string input;
  string output;
  string plysource;
  bool info = false;

  unsigned int add = 0;
  bool add_strip = false;
  bool add_colors = false;
  bool add_normals = false;
  bool add_textures = false;
  bool add_data = false;

  unsigned int remove = 0;
  bool remove_strip = false;
  bool remove_colors = false;
  bool remove_normals = false;
  bool remove_textures = false;
  bool remove_data = false;

  bool compress = false;
  bool uncompress = false;
  float qvertex = 0;
  float qnormal = 0;
  float qcolor = 0;
  float qtexture = 0;

  int option;
  while((option = getopt(argc, argv, "io:a:r:zxv:n:c:t:")) != EOF) {
    switch(option) {
    case 'i': info = true; break;
    case 'o': output = optarg; break;
    case 'a': {
      if(strstr(optarg, "strip")) {
	add_strip = true;
	add |= NXS_STRIP;
	remove |= NXS_FACES;
      }
      if(strstr(optarg, "colors")) {
	add_colors = true;
	add |= NXS_COLORS;
      }
      if(strstr(optarg, "normals")) {
	add_normals = true;
	add |= NXS_NORMALS_SHORT;
      }
      if(strstr(optarg, "textures")) {
	add_textures = true;
	add |= NXS_TEXTURES_SHORT;
      }
      if(strstr(optarg, "data")) {
	add_data = true;
	add |= NXS_DATA32;
      }
      if(add == 0) {
	cerr << "Invalid -a argument: " << optarg << "\n"
	     << "Valid options are: strip, colors, normals, textures, data\n";
	return -1;
      }
      break;
    }
      
    case 'r': {
      if(strstr(optarg, "strip")) {
	cerr << "Strip reming not supported!\n";
	return -1;
	remove_strip = true;
	add |= NXS_FACES;
	remove |= NXS_STRIP;
      }
      if(strstr(optarg, "colors")) {
	remove_colors = true;
	remove |= NXS_COLORS;
      }
      if(strstr(optarg, "normals")) {
	remove_normals = true;
	remove |= NXS_NORMALS_SHORT;
      }
      if(strstr(optarg, "textures")) {
	remove_textures = true;
	remove |= NXS_TEXTURES_SHORT;
      }
      if(strstr(optarg, "data")) {
	remove_data = true;
	remove |= NXS_DATA32;
      }
      if(remove == 0) {
	cerr << "Invalid -a argument: " << optarg << "\n"
	     << "Valid options are: strip, colors, normals, textures, data\n";
	return -1;
      }
      break;
    }
      
    case 'p': plysource = optarg; break;
    case 'z': compress = true; break;
    case 'x': uncompress = true; break;

    case 'v': qvertex = atof(optarg); 
      if(qvertex == 0) {
	cerr << "Invalid value for quantization: " << optarg << endl;
	return -1;
      }
      break;
    case 'n': qnormal = atof(optarg); 
      if(qnormal == 0) {
	cerr << "Invalid value for quantization: " << optarg << endl;
	return -1;
      }
      break;
    case 'c': qcolor = atof(optarg); 
      if(qcolor == 0) {
	cerr << "Invalid value for quantization: " << optarg << endl;
	return -1;
      }
      break;
    case 't': qtexture = atof(optarg); 
      if(qtexture == 0) {
	cerr << "Invalid value for quantization: " << optarg << endl;
	return -1;
      }
      break;

    default: cerr << "Unknown option: " << (char)option << endl;
      return -1;
    }
  }

  if(compress && uncompress) {
    cerr << "x and z are obviously exclusive :P\n";
    return -1;
  }

  if(optind != argc - 1) {
    cerr << "Usage: " << argv[0] << " <nexus file> [options]\n"
	 << " -i       : display some info about nexus file\n"
         << " -o <file>: output filename (default is adding 00 to nexus)\n"
	 << " -a <what>: Add [colors|normals|strip|textures|data|borders]\n"
         << " -r <what>: As add...\n"
	 << " -p <ply> : Ply source for colors or textures or data\n"
	 << " -z       : compress\n"
	 << " -x       : uncompress\n"
	 << " -v<float>: Vertex quantization (float is the 0 level amount)\n"
	 << " -n<float>: Normal quantization\n"
         << " -c<float>: Color quantization\n"
	 << " -t<float>: Texture quantization\n\n";
  }
  input = argv[optind];
  if(!output.size()) output = input + "00";

  Nexus nexus;
  if(!nexus.Load(input)) {
    cerr << "Could not open nexus file: " << input << ".mt\n";
    return -1;
  }

  //Sanity tests
  if(remove_strip && !(nexus.signature & NXS_STRIP)) {
    cerr << "Nexus file does not have strips\n";
    return -1;
  }
  if(remove_colors && !(nexus.signature & NXS_COLORS)) {
    cerr << "Nexus file does not have colors\n";
    return -1;
  }
  if(remove_normals && !(nexus.signature & NXS_NORMALS_SHORT)) {
    cerr << "Nexus file does not have normals\n";
    return -1;
  }
  if(remove_textures && !(nexus.signature & NXS_TEXTURES_SHORT)) {
    cerr << "Nexus file does not have textures\n";
    return -1;
  }
  if(remove_data && !(nexus.signature & NXS_DATA32)) {
    cerr << "Nexus file does not have data\n";
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
    cerr << "Nexus file: " << input << "\n\n"
	 << "\n\tCompressed: " << nexus.IsCompressed() 
	 << "\n\tStripped: " << (int)((nexus.signature&NXS_STRIP) !=0)
	 << "\n\tColor   : " << (int)((nexus.signature&NXS_COLORS) !=0)
	 << "\n\tNormal  : " << (int)((nexus.signature&NXS_NORMALS_SHORT) !=0)
	 << "\n\tTexture : " << (int)((nexus.signature&NXS_TEXTURES_SHORT) !=0)
	 << "\n\tData    : " << (int)((nexus.signature&NXS_DATA32) !=0)
	 << "\n\n\tVertices: " << nexus.totvert 
	 << "\tFaces: " << nexus.totface
	 << "\tPatches: " << nexus.index.size() << "\n\n";
  }
  
  //determine if we must proceed:
  if(add == 0 && remove == 0 && !compress && !uncompress &&
     qvertex == 0 && qnormal == 0 && qcolor == 0 && qtexture == 0) {
    nexus.Close();
    return 0;
  }

  unsigned int signature = nexus.signature;
  signature |= add;
  signature &= ~remove;

  Nexus out;
  if(!out.Create(output, (Signature)signature)) {
    cerr << "Could not open output: " << output << endl;
    return -1;
  }

  for(unsigned int patch = 0; patch < nexus.index.size(); patch++) {
    Nexus::Entry &src_entry = nexus.index[patch];
    Patch src_patch = nexus.GetPatch(patch);
    Border src_border = nexus.GetBorder(patch);


    vector<unsigned short> strip;
    if(add_strip) {
      ComputeTriStrip(src_patch.nf, src_patch.FaceBegin(), strip);
      out.AddPatch(src_entry.nvert, strip.size(), src_entry.border_size);
    } else
      out.AddPatch(src_entry.nvert, src_entry.nface, src_entry.border_size);


    Nexus::Entry &dst_entry = out.index[patch];
    Patch dst_patch = out.GetPatch(patch);
    Border dst_border = out.GetBorder(patch);

    
    //copy vertices: //no clustering
    memcpy(dst_patch.VertBegin(), src_patch.VertBegin(), 
	   src_patch.nv * sizeof(Point3f));

    //now faces.
    if(add_strip) {
      memcpy(dst_patch.FaceBegin(), &*strip.begin(), 
	     strip.size() * sizeof(short));
    } else {
      memcpy(dst_patch.FaceBegin(), src_patch.FaceBegin(), 
	     src_patch.nf * sizeof(unsigned short) *3);
    }

    if((nexus.signature & NXS_COLORS) && (out.signature & NXS_COLORS))
      memcpy(dst_patch.ColorBegin(), src_patch.ColorBegin(), 
	     src_patch.nv * sizeof(unsigned int));

    if((nexus.signature & NXS_NORMALS_SHORT) && 
       (out.signature & NXS_NORMALS_SHORT))
      memcpy(dst_patch.Norm16Begin(), src_patch.Norm16Begin(), 
	     src_patch.nv * sizeof(short)*4);
    

    //copying entry information;
    dst_entry.sphere = src_entry.sphere;
    dst_entry.error = src_entry.error;

    //adding borders.
    dst_entry.border_used = src_entry.border_used;
    memcpy(dst_border.Start(), src_border.Start(), 
	   src_border.Size() * sizeof(Link));
  }

  //TODO this is ok only if we have faces still!
  if(add_normals) {
    cerr << "Computing normals" << endl;
    ComputeNormals(out);
  }

  if(add_colors) {
    //source of color:
    cerr << "Unsupported color\n";
    return -1;
  }

  //fixing sphere.
  out.sphere = nexus.sphere;
  //copying history:
  out.history = nexus.history;

  out.Close();
  nexus.Close();
  return 0;
}

