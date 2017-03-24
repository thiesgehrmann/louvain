// File: zahn.h
// -- quality functions (for Zahn-Condorcet criterion) header file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// And based on the article
// Copyright (C) 2013 R. Campigotto, P. Conde Céspedes, J.-L. Guillaume
//
// This file is part of Louvain algorithm.
// 
// Louvain algorithm is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Louvain algorithm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : July 2013
//-----------------------------------------------------------------------------
// see README.txt for more details


#ifndef ZAHN_H
#define ZAHN_H

#include "quality.h"

using namespace std;


class Zahn: public Quality {
 public:
  
  // used to compute the quality participation of each community
  vector<float> in;
  vector<int> w;
  float max; // biggest weight on links

  Zahn(Graph & gr, float max_w);
  ~Zahn();

  inline void remove(int node, int comm, float dnodecomm);

  inline void insert(int node, int comm, float dnodecomm);

  inline float gain(int node, int comm, float dnodecomm, float w_degree);

  float quality();
};


inline void
Zahn::remove(int node, int comm, float dnodecomm) {
  assert(node>=0 && node<size);

  in[comm] -= 2.0*dnodecomm + g.nb_selfloops(node);
  w[comm]  -= g.nodes_w[node];
  
  n2c[node] = -1;
}

inline void
Zahn::insert(int node, int comm, float dnodecomm) {
  assert(node>=0 && node<size);
  
  in[comm] += 2.0*dnodecomm + g.nb_selfloops(node);
  w[comm]  += g.nodes_w[node];
  
  n2c[node] = comm;
}

inline float
Zahn::gain(int node, int comm, float dnc, float degc) {
  assert(node>=0 && node<size);
  
  float wc = (float)w[comm];
  float wu = (float)g.nodes_w[node];
  
  float gain = 2.0*dnc - wu*wc*max;

  return gain;
}


#endif // ZAHN_H
