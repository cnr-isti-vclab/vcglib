#include <iostream>
using namespace std;

#include "nxsalgo.h"
#include "nexus.h"
#include "watch.h"

using namespace nxs;
using namespace vcg;

#ifdef WIN32
#include <wrap/system/getopt.h>
#else
#include <unistd.h>
#endif

#include <assert.h>


string getSuffix(unsigned int signature) {
  string suff;
  if(signature&NXS_COMPRESSED)     suff += "Z";
  if(signature&NXS_STRIP)          suff += "S"; 
  if(signature&NXS_COLORS)         suff += "C";
  if(signature&NXS_NORMALS_SHORT)  suff += "N";
  if(signature&NXS_TEXTURES_SHORT) suff += "T";
  if(signature&NXS_DATA32)         suff += "D";
  return suff;
}

int main(int argc, char *argv[]) {
  
  string input;
  string output;
  string plysource;
  bool info = false;
  unsigned int ram_size = 128000000;
  unsigned int chunk_size = 0;

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
  while((option = getopt(argc, argv, "io:a:r:zxv:n:k:t:b:c:")) != EOF) {
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
    case 'k': qcolor = atof(optarg); 
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
    return -1;
  }
  input = argv[optind];


  Nexus nexus;
  nexus.patches.SetRamBufferSize(ram_size);
  if(!nexus.Load(input, true)) {
    cerr << "Could not open nexus file: " << input << "\n";
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
	 << "\tPatches: " << nexus.index.size() 
	 << "\n\tSphere: " 
	 << nexus.sphere.Center()[0] << " "
	 << nexus.sphere.Center()[1] << " "
	 << nexus.sphere.Center()[2] << " R: "
	 << nexus.sphere.Radius()
	 << "\n\tChunk size " << nexus.chunk_size << endl;
    
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
  if(compress) signature |= NXS_COMPRESSED;
  if(uncompress) signature &= ~NXS_COMPRESSED;

  if(!output.size()) output = input + getSuffix(signature);

  cout << "Writing to nexus: " << output << endl;

  Nexus out;
  out.patches.SetRamBufferSize(ram_size);
  if(!chunk_size)
    chunk_size = nexus.patches.chunk_size;

  if(!out.Create(output, (Signature)signature, chunk_size)) {
    cerr << "Could not open output: " << output << endl;
    return -1;
  }

  //TODO set rambuffer low (or even direct access!)

  Report report(nexus.index.size());
  cout << "Copying and allocating...\n";
  for(unsigned int patch = 0; patch < nexus.index.size(); patch++) {
    report.Step(patch);
    Nexus::PatchInfo &src_entry = nexus.index[patch];
    Patch src_patch = nexus.GetPatch(patch);
    Border src_border = nexus.GetBorder(patch);


    vector<unsigned short> strip;
    if(add_strip) {
      ComputeTriStrip(src_patch.nf, src_patch.FaceBegin(), strip);
      assert(strip.size() < 32767);
      out.AddPatch(src_entry.nvert, strip.size(), src_border.Available());
    } else
      out.AddPatch(src_entry.nvert, src_entry.nface, src_border.Available());


    Nexus::PatchInfo &dst_entry = out.index[patch];
    Patch dst_patch = out.GetPatch(patch);

    //copy vertices: 
    memcpy(dst_patch.VertBegin(), src_patch.VertBegin(), 
	   src_patch.nv * sizeof(Point3f));

    if(qvertex) {
      float *ptr = (float *)dst_patch.VertBegin();
      for(unsigned int i = 0; i < dst_patch.nv*3; i++) {
	ptr[i] =  qvertex * nearbyintf(ptr[i]/qvertex);
	//ptr[i] = 0;
      }
    } 


    //now faces.
    if(add_strip) {
      memcpy(dst_patch.FaceBegin(), &*strip.begin(), 
	     strip.size() * sizeof(short));
    } else {
      if(nexus.signature & NXS_STRIP) {
	memcpy(dst_patch.FaceBegin(), src_patch.FaceBegin(), 
	       src_patch.nf * sizeof(unsigned short));
      } else {
	memcpy(dst_patch.FaceBegin(), src_patch.FaceBegin(), 
	       src_patch.nf * sizeof(unsigned short) * 3);
      }
    }

    if((nexus.signature & NXS_COLORS) && (out.signature & NXS_COLORS))
      memcpy(dst_patch.ColorBegin(), src_patch.ColorBegin(), 
	     src_patch.nv * sizeof(unsigned int));

    if((nexus.signature & NXS_NORMALS_SHORT) && 
       (out.signature & NXS_NORMALS_SHORT))
      memcpy(dst_patch.Norm16Begin(), src_patch.Norm16Begin(), 
	     src_patch.nv * sizeof(short)*4);

    //reordering
    //WATCH OUT BORDERS!
    //    Reorder(out.signature, dst_patch);
    //copying entry information;
    dst_entry.sphere = src_entry.sphere;
    dst_entry.error = src_entry.error;

    //adding borders.
    for(unsigned int i = 0; i < src_border.Size(); i++) {
      Link &link = src_border[i];
      if(link.IsNull()) continue;
      assert(link.end_patch < nexus.index.size());
    }
    Border dst_border = out.GetBorder(patch);
    out.borders.ResizeBorder(patch, src_border.Size());
    memcpy(dst_border.Start(), src_border.Start(), 
	   src_border.Size() * sizeof(Link));
  }
  report.Finish();

  //TODO this is ok only if we have faces still!
  if(add_normals) {
    cout << "Computing normals" << endl;
    ComputeNormals(out);
  }

  if(add_colors) {
    //source of color:
    cerr << "Unsupported color\n";
    return -1;
  }
  if(qvertex) { 
    report.Init(nexus.index.size());
    cout << "Quantizing vertices\n";
    for(unsigned int patch = 0; patch < nexus.index.size(); patch++) {
      report.Step(patch);
      Patch src_patch = nexus.GetPatch(patch);

      float *ptr = (float *)src_patch.VertBegin();
      for(unsigned int i = 0; i < src_patch.nv*3; i++) 
	ptr[i] =  qvertex * nearbyintf(ptr[i]/qvertex);
    }
    report.Finish();
  }

  //fixing sphere.
  out.sphere = nexus.sphere;
  //copying history:
  out.history = nexus.history;

  out.Close();
  nexus.Close();
  return 0;
}

