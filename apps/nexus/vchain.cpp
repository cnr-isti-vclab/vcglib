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

****************************************************************************/

#include "vchain.h"
#include <stdio.h>
#include <iostream>

using namespace std;
using namespace vcg;
using namespace nxs;

VChain::~VChain() {
  for(iterator i = begin(); i != end(); i++)
    delete *i;
}

bool VChain::Save(const string &file) {
  FILE *fp = fopen(file.c_str(), "wb+");
  if(!fp) {
    cerr << "Could not save vchain data.\n";
    return false;
  }

  unsigned int nlevels = size();
  fwrite(&nlevels, sizeof(unsigned int), 1, fp);
  for(unsigned int i = 0; i < nlevels; i++) {
    VPartition &level = *operator[](i);
    unsigned int npoints = level.size();
    fwrite(&npoints, sizeof(unsigned int), 1, fp);
    fwrite(&(level[0]), sizeof(Point3f), npoints, fp);
  }
  //writing fragments

  unsigned int nfrag = newfragments.size();
  fwrite(&nfrag, sizeof(unsigned int), 1, fp);

  std::map<unsigned int, std::set<unsigned int> >::iterator j;
  for(j = newfragments.begin(); j != newfragments.end(); j++) {
    unsigned int n = (*j).second.size();
    fwrite(&((*j).first), sizeof(unsigned int), 1, fp);
    fwrite(&n, sizeof(unsigned int), 1, fp);
    set<unsigned int>::iterator k;
      for(k = (*j).second.begin(); k != (*j).second.end(); k++)
        fwrite(&*k, sizeof(unsigned int), 1, fp);
  }
  nfrag = oldfragments.size();
  fwrite(&nfrag, sizeof(unsigned int), 1, fp);    

  for(j = oldfragments.begin(); j != oldfragments.end(); j++) {
    unsigned int n = (*j).second.size();
    fwrite(&((*j).first), sizeof(unsigned int), 1, fp);
    fwrite(&n, sizeof(unsigned int), 1, fp);
    set<unsigned int>::iterator k;
    for(k = (*j).second.begin(); k != (*j).second.end(); k++)
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
  for(unsigned int i = 0; i < nlevels; i++) {
    push_back(new VPartition());
    VPartition &level = *back();

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
