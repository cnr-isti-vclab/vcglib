#include "vchain.h"
#include <stdio.h>
#include <iostream>

using namespace std;
using namespace vcg;
using namespace nxs;

bool VChain::Save(const string &file) {
  FILE *fp = fopen(file.c_str(), "wb+");
  if(!fp) {
    cerr << "Could not save vchain data.\n";
    return false;
  }

  unsigned int nlevels = size();
  fwrite(&nlevels, sizeof(unsigned int), 1, fp);
  for(int i = 0; i < nlevels; i++) {
    VPartition &level = operator[](i);
    unsigned int npoints = level.size();
    fwrite(&npoints, sizeof(unsigned int), 1, fp);
    fwrite(&(level[0]), sizeof(Point3f), npoints, fp);
  }
  //writing fragments

  unsigned int nfrag = newfragments.size();
  fwrite(&nfrag, sizeof(unsigned int), 1, fp);

  std::map<unsigned int, std::set<unsigned int> >::iterator i;
  for(i = newfragments.begin(); i != newfragments.end(); i++) {
    unsigned int n = (*i).second.size();
    fwrite(&((*i).first), sizeof(unsigned int), 1, fp);
    fwrite(&n, sizeof(unsigned int), 1, fp);
    set<unsigned int>::iterator k;
      for(k = (*i).second.begin(); k != (*i).second.end(); k++)
        fwrite(&*k, sizeof(unsigned int), 1, fp);
  }
  nfrag = oldfragments.size();
  fwrite(&nfrag, sizeof(unsigned int), 1, fp);    

  for(i = oldfragments.begin(); i != oldfragments.end(); i++) {
    unsigned int n = (*i).second.size();
    fwrite(&((*i).first), sizeof(unsigned int), 1, fp);
    fwrite(&n, sizeof(unsigned int), 1, fp);
    set<unsigned int>::iterator k;
    for(k = (*i).second.begin(); k != (*i).second.end(); k++)
      fwrite(&*k, sizeof(unsigned int), 1, fp);
  }
  fclose(fp);
  return true;
}

bool VChain::Load(const string &file) {
  FILE *fp = fopen(file.c_str(), "rb");
  if(!fp) {
    cerr << "Could not load vchain data\n";
    return false;
  }
  unsigned int nlevels;
  fread(&nlevels, sizeof(unsigned int), 1, fp);
  for(int i = 0; i < nlevels; i++) {
    push_back(VPartition());
    VPartition &level = back();

    unsigned int npoints;
    fread(&npoints, sizeof(unsigned int), 1, fp);
    level.resize(npoints);
    fread(&(level[0]), sizeof(Point3f), npoints, fp);
    level.Init();
  }
  //reading fragments
  unsigned int nfrag;
  fread(&nfrag, sizeof(unsigned int), 1, fp);
  for(unsigned int i = 0; i < nfrag; i++) {
    unsigned int p, n;
    fread(&p, sizeof(unsigned int), 1, fp);
    set<unsigned int> &s = newfragments[p];
    fread(&n, sizeof(unsigned int), 1, fp);
    for(unsigned int k = 0; k < n; k++) {
      unsigned int j;
      fread(&j, sizeof(unsigned int), 1, fp);        
      s.insert(j);
    }
  }
  
  fread(&nfrag, sizeof(unsigned int), 1, fp);
  for(unsigned int i = 0; i < nfrag; i++) {
    unsigned int p, n;
    fread(&p, sizeof(unsigned int), 1, fp);
    set<unsigned int> &s = oldfragments[p];
    fread(&n, sizeof(unsigned int), 1, fp);
    for(unsigned int k = 0; k < n; k++) {
      unsigned int j;
      fread(&j, sizeof(unsigned int), 1, fp);        
      s.insert(j);
    }
  }
  fclose(fp);
  return true;
}
