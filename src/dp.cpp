// File: dp.cpp
// -- quality functions (for DP criterion) source file
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


#include "dp.h"

using namespace std;


DP::DP(Graph & gr, float sum, float max_w):Quality(gr,"Profile Difference"),kappa(size),sum_sq(sum),max(max_w) {
  n2c.resize(size);

  in.resize(size);
  w.resize(size);

  // initialization
  for (int i=0 ; i<size ; i++) {
    n2c[i] = i;
    in[i]  = g.nb_selfloops(i);
    w[i]   = g.nodes_w[i];
  }
}

DP::~DP() {
  in.clear();
  w.clear();
}

float
DP::graph_weighting(Graph *g) {
  float sum_sq = 0.0;
  
  vector<float> aux_weights;
  
  // foreach weight, change Aij to 2Aij / (d(i)+d(j))
  for (int u=0 ; u < g->nb_nodes ; u++) {
    pair<vector<int>::iterator, vector<float>::iterator> p = g->neighbors(u);
    int deg = g->nb_neighbors(u);
    for (int i=0 ; i < deg ; i++) {
      int neigh = *(p.first+i);
      float neigh_w = 0.0;
      
      if (g->weights.size() == 0)
	neigh_w = 2.0 / ((float)deg + (float)(g->nb_neighbors(neigh)));
      else {
	float old_neigh = (float)*(p.second+i);
	neigh_w = 2.0*old_neigh / ((float)deg + (float)(g->nb_neighbors(neigh)));
      }
      
      aux_weights.push_back(neigh_w);
      
      sum_sq += neigh_w*neigh_w;
    }
  }
  
  g->weights.clear();
  g->weights = aux_weights;
  
  g->total_weight = 0.0;
  
  // Compute total weight
  for (int i=0 ; i < g->nb_nodes ; i++)
    g->total_weight += (float)(g->weighted_degree(i));

  aux_weights.clear();

  return sum_sq;
}


float
DP::quality() {
  float q = 0.0;
  float n = (float)g.sum_nodes_w;

  for (int i=0 ; i < size ; i++) {
    float wc = (float)w[i];
    if (wc > 0.0)
      q += 2.0*in[i] / wc;
  }
  
  q -= (sum_sq + (float)kappa);

  q /= n*n*max;
  
  return q;
}
