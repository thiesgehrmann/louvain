// File: dp.h
// -- quality functions (for Profile Difference criterion) header file
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


#ifndef DP_H
#define DP_H

#include "quality.h"

using namespace std;


class DP: public Quality {
 public:

  vector<float> in; // used to compute the quality participation of each community
  vector<int> w; // used to store size of communities
  int kappa; // number of communities

  float sum_sq; // used to store Sum Âij^2
  float max;

  DP(Graph & gr, float sum, float max_w);
  ~DP();

  // change the weight of each link ij in the graph, from Aij to 2Aij / (d(i)+d(j))
  // return the result of Sum Âij^2
  static float graph_weighting(Graph *g);

  inline void remove(int node, int comm, float dnodecomm);

  inline void insert(int node, int comm, float dnodecomm);

  inline float gain(int node, int comm, float dnodecomm, float w_degree);

  float quality();
};


inline void
DP::remove(int node, int comm, float dnodecomm) {
  assert(node>=0 && node<size);

  in[comm] -= 2.0*dnodecomm + g.nb_selfloops(node);

  w[comm] -= g.nodes_w[node];

  if (w[comm] == 0)
    kappa--;

  n2c[node] = -1;
}

inline void
DP::insert(int node, int comm, float dnodecomm) {
  assert(node>=0 && node<size);

  in[comm] += 2.0*dnodecomm + g.nb_selfloops(node);
  
  if (w[comm] == 0)
    kappa++;
  
  w[comm] += g.nodes_w[node];

  n2c[node] = comm;
}

inline float
DP::gain(int node, int comm, float dnc, float degc) {
  assert(node>=0 && node<size);

  float inc  = in[comm];
  float self = g.nb_selfloops(node);
  float wc   = (float)w[comm];
  float wu   = (float)g.nodes_w[node];
  
  float gain;

  if (wc == 0.0) {
    gain  = (2.0*dnc+self) / wu;
    gain -= 1.0/2.0;
  }
  else {
    gain  = (inc + 2.0*dnc+self) / (wc+wu);
    gain -= inc/wc;
  }

  return gain;
}


#endif // DP_H
