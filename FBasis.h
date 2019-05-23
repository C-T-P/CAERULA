// Copyright (C) 2018 Christian T Preuss
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

//*********************************************************************************************************
//
// Class FType
//
//*********************************************************************************************************

class FType {
    
  vector<size_t> m_g;

 public:
  FType(vector<size_t> g_inds = {});
  ~FType();
  vector<FType> add_one_gluon(size_t g_ind);
  vector<size_t> get_indices();
  size_t no_g();
  bool is_not_empty();
  bool vanishes();
  bool operator>(FType& rhs);
  bool operator==(FType rhs);
  CTerm build_ct(size_t start_ind);
  void print();
};

//*********************************************************************************************************
//
// Class FVec
//
//*********************************************************************************************************

class FVec {
  
  vector<FType> m_FVec;
    
 public:
  FVec(FType f = FType(vector<size_t>({})));
  ~FVec();
  void push_back(FType f);
  FType& at(size_t i);
  FType at(size_t i) const;
  vector<FVec> add_one_gluon(size_t g_ind);
  vector<size_t> get_indices();
  size_t no_groups();
  //        bool is_tree_level();
  bool has_sg();
  void order();
  bool operator>(FVec& rhs);
  CAmplitude build_ca();
  void print();
};

//*********************************************************************************************************
//
// Class FBasis
//
//*********************************************************************************************************

class FBasis : public CBasis {
  
  vector<FVec> m_FBasis;
  vector<size_t> m_g_indices;
  size_t m_ng;
  
  void make_perms();
  void make_ca_basis();
  
  void remove_sg();
  void normal_order();
  
 public:
  FBasis(size_t n_g);
  ~FBasis();
};

#endif
