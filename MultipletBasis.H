// Copyright (C) 2021 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#ifndef MULTIPLETBASIS_H
#define MULTIPLETBASIS_H

#include "GenBasis.H"
#include "TraceBasis.H"

using namespace std;

namespace SPECTRUM {

//*****************************************************************************
//
// Class SU3Irrep
//
//*****************************************************************************

class SU3Irrep {

 private:
  // Labels.
  int m_p;
  int m_q;
  
  // Dimension.
  int m_dim;

 public:
  // Default constructor by labels.
  SU3Irrep(int p, int q) {
    m_p = p; m_q = q;
    m_dim = 0.5*(m_p+1)*(m_q+1)*(m_p+m_q+2);
    if (m_q > m_p) m_dim = -m_dim;
  }

  // Member functions.
  vector<SU3Irrep> operator*(SU3Irrep rhs) {
    vector<SU3Irrep> newIrreps;
    int p2 = rhs.m_p, q2 = rhs.m_q;
    
    for (int i(0); i <= min(m_p,q2); ++i) {
      for (int j(0); j <= min(p2,m_q); ++j) {
        int n1 = m_p-i;
        int n2 = p2-j;
        int m1 = m_q-j;
        int m2 = q2-i;
        
        if (n1 < 0 or m1 < 0 or n2 < 0 or m2 < 0) break;
        
        // Represent highest rep as irrep.
        newIrreps.push_back(SU3Irrep(n1+n2,m1+m2));
        
        // Save the others if existent.
        for (int k(1); k <= min(n1,n2); ++k) {
          int pNew = n1 + n2 - 2*k;
          int qNew = m1 + m2 + k;
          if (qNew >= 0 and pNew >= 0)
            newIrreps.push_back(SU3Irrep(pNew, qNew));
        }
        for (int k(1); k <= min(m1,m2); ++k) {
          int pNew = n1 + n2 + k;
          int qNew = m1 + m2 - 2*k;
          if (qNew >= 0 and pNew >= 0)
            newIrreps.push_back(SU3Irrep(pNew, qNew));
        }
      }
    }
    
    return newIrreps;
  }
  vector<SU3Irrep> cg_decomposition(vector<SU3Irrep> multiplets) {
    vector<SU3Irrep> decomposition;
    size_t nMultiplets = multiplets.size();

    // Simple cases.
    if (nMultiplets == 1) return multiplets;
    if (nMultiplets == 2) return multiplets.at(0)*multiplets.at(1);

    decomposition = multiplets.at(0)*multiplets.at(1);
    for (int i(2); i < nMultiplets; ++i) {
      vector<SU3Irrep> tempDecomp;
      for (int j(0); j < decomposition.size(); ++j) {
        vector<SU3Irrep> tmp(decomposition.at(j)*multiplets.at(i));
        for (auto irrep : tmp) tempDecomp.push_back(irrep);
      }
      decomposition.clear();
      decomposition = tempDecomp;
    }

    return decomposition;
  }

  // Getters.
  size_t dim() { return m_dim; }
  string get_string() { return "["+to_string(m_dim)+"]"; }
};

//*****************************************************************************
//
// Class MultipletBasis
//
//*****************************************************************************

class MultipletBasis : public GenBasis {
    
  matrix m_bcm;
  size_t m_ng, m_nqp;
  
 public:
  MultipletBasis(size_t n_g, size_t n_qp);
  ~MultipletBasis();

  // Create Multiplet Basis.
  //  void make_basis(int irrepA, int irrepB, vector<int> irrepsFinal) {;}

  matrix bcm();

  friend class CBasis;
};

}

#endif
