// File: devind.h
// -- quality functions (for Deviation to Indetermination criterion) header file
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


#ifndef DEVIND_H
#define DEVIND_H

#include "quality.h"

using namespace std;


class DevInd: public Quality {
 public:

  vector<float> in, tot; // used to compute the quality participation of each community
  vector<int> w; // used to store size of communities

  DevInd(Graph & gr);
  ~DevInd();

  inline void remove(int node, int comm, float dnodecomm);

  inline void insert(int node, int comm, float dnodecomm);

  inline float gain(int node, int comm, float dnodecomm, float w_degree);

  float quality();
};


inline void
DevInd::remove(int node, int comm, float dnodecomm) {
  assert(node>=0 && node<size);

  in[comm]  -= 2.0*dnodecomm + g.nb_selfloops(node);
  tot[comm] -= g.weighted_degree(node);
  w[comm]   -= g.nodes_w[node];

  n2c[node] = -1;
}

inline void
DevInd::insert(int node, int comm, float dnodecomm) {
  assert(node>=0 && node<size);

  in[comm]  += 2.0*dnodecomm + g.nb_selfloops(node);
  tot[comm] += g.weighted_degree(node);
  w[comm]   += g.nodes_w[node];

  n2c[node] = comm;
}

inline float
DevInd::gain(int node, int comm, float dnc, float degc) {
  assert(node>=0 && node<size);

  float totc = tot[comm];
  float wc   = (float)w[comm];
  float wu   = (float)g.nodes_w[node];

  float m2   = g.total_weight;
  float n    = (float)g.sum_nodes_w;

  float gain = dnc - (totc*wu + degc*wc)/n + (m2*wu*wc)/(n*n);

  return gain;
}


#endif // DEVIND_H
