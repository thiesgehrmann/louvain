// File: devuni.cpp
// -- quality functions (for Deviation to Uniformity criterion) source file
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


#include "devuni.h"

using namespace std;


DevUni::DevUni(Graph & gr):Quality(gr,"Deviation to Uniformity") {
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

DevUni::~DevUni() {
  in.clear();
  w.clear();
}

float
DevUni::quality() {
  float q  = 0.0;
  float n  = (float)g.sum_nodes_w;
  float m2 = g.total_weight;

  float sum = 0.0;

  for (int i=0 ; i<size ; i++) {
    float wc = (float)w[i];
    if (wc > 0.0) {
      q += in[i];
      sum += wc*wc;
    }
  }
  
  q -= sum * (m2/(n*n));

  q /= m2;
  
  return q;
}
