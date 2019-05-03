// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#include<fstream>
#include "multiplet_basis.h"
#include "trace_basis.h"
#include "gen_basis.h"

string get_filename(size_t n_g, size_t n_qp);

// member functions of multiplet_basis class
multiplet_basis::multiplet_basis(size_t n_g, size_t n_qp) : gen_basis(get_filename(n_g, n_qp)) {
  // set basis type
  m_btype=1;
    
  // set number of particles
  m_ng=n_g;
  m_nqp=n_qp;
    
  // initialise basis change matrix
  m_bcm=matrix();
}

multiplet_basis::~multiplet_basis() {
    
}

matrix multiplet_basis::bcm() {
  trace_basis tr_basis(m_ng, m_nqp);
    
  this->normalise();
  tr_basis.normalise();
    
  // get permutations and normalisations for trace basis
  m_amp_perms=tr_basis.m_amp_perms;
  m_normalisations=tr_basis.m_normalisations;
  m_confact=tr_basis.m_confact;
    
  size_t tr_dim(tr_basis.dim());
    
  matrix tmp_bcm(m_dim, vector<complex<double>>(tr_dim, 0.));
    
  for (size_t i(0); i<m_dim; i++) {
    for (size_t j(0); j<tr_dim; j++) {
      tmp_bcm.at(i).at(j) = m_ca_basis.at(i).scprod(tr_basis.m_ca_basis.at(j));
    }
  }
    
  m_bcm = tmp_bcm;
    
  return m_bcm;
}

// helper function
string get_filename(size_t n_g, size_t n_qp) {
  string filename="";
  if (n_qp!=0) filename+=to_string(n_qp)+"qqb";
  if (n_g!=0) filename+=to_string(n_g)+"g";
    
  filename+="-multiplet.txt";
  
  ifstream fin("precalc_multiplet_bases/"+filename);
  if (!fin) {
    cerr<<"Error: No precalculated multiplet basis file "<< filename <<" found in /precalc_multiplet_bases/"<<endl;
    exit(EXIT_FAILURE);
  }
  filename="precalc_multiplet_bases/"+filename;
  return filename;
}
