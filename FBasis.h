// Copyright (C) 2021 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#ifndef FBASIS_H
#define FBASIS_H

#include "CBasis.h"

using namespace std;

namespace SPECTRUM {

//*****************************************************************************
//
// Class FType
//
//*****************************************************************************

class FType {
  
  // Vector of adjoint indices
  vector<size_t> m_g;

 public:
  // Default constructor
  FType(vector<size_t> g_inds = {});
  ~FType();

  // Method to add one gluon in all possible ways
  vector<FType> add_one_gluon(size_t g_ind);

  // Getters
  vector<size_t> get_indices();
  size_t no_g();
  bool is_not_empty();
  bool vanishes();

  // Operator definitions
  bool operator>(FType& rhs);
  bool operator==(FType rhs);

  // Printers
  CTerm build_ct(size_t start_ind);
  void print();
};

//*****************************************************************************
//
// Class FVec
//
//*****************************************************************************

class FVec {
  
  // Vector of FTypes - representing one amplitude consisting only of f_ijk's
  vector<FType> m_FVec;
    
 public:
  // Default constructor
  FVec(FType f = FType(vector<size_t>({})));
  ~FVec();
  
  // Method to add one gluon in all possible ways - returns all possible FVecs
  vector<FVec> add_one_gluon(size_t g_ind);

  // Method to order the FTypes
  void order();

  // Setters
  FType& at(size_t i);
  void push_back(FType f);

  // Getters
  vector<size_t> get_indices();
  FType at(size_t i) const;
  size_t no_groups();
  bool has_sg();

  // Operators
  bool operator>(FVec& rhs);

  // Builders
  CAmplitude build_ca();
  void print();
};

//*****************************************************************************
//
// Class FBasis
//
//*****************************************************************************

class FBasis : public CBasis {
  
  // Number of gluons
  size_t m_ng;

  // All gluonic indices
  vector<size_t> m_g_indices;
  
  // Complete FBasis in terms of FVecs
  vector<FVec> m_FBasis;
  
  // Method to make all permutations
  void make_perms();

  // Method to create a colour amplitude basis
  void make_ca_basis();
  
  // Remove single traced gluons
  void remove_sg();
  
  // Method to order the basis
  void normal_order();
  
 public:
  // Constructor by number of gluons
  FBasis(size_t n_g);
  ~FBasis();
};

}

#endif
