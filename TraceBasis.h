// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#ifndef TRACEBASIS_H
#define TRACEBASIS_H

#include "CBasis.h"
#include "Colourtools.h"

using namespace std;

//*********************************************************************************************************
//
// Class TraceType
//
//********************************************************************************************************* 

class TraceType {
  // The indices of the quark and antiquark
  size_t m_qb, m_q;

  // The indices of the gluons
  vector<size_t> m_g;

 public:
  // Default constructor
  TraceType(vector<size_t> g_inds = {}, size_t q_ind = 0, size_t qb_ind = 0);
  
  // Constructor from quark and antiquark indices
  TraceType(size_t q_ind, size_t qb_ind);

  //Destructor
  ~TraceType();
  
  // Method to add one gluonic index
  vector<TraceType> add_one_gluon(size_t g_ind);

  // Method to get all indices
  vector<size_t> get_indices();

  // Method to compute conjugate
  TraceType conj();

  // Method to get number of gluons
  size_t no_g();

  // Method to get number of quark pairs
  size_t no_qp();

  // Method to check whether trace has indices
  bool is_not_empty();
  
  // Method to check whether the trace vanishes
  bool vanishes();

  // Method to check whether trace type is bigger than a second one
  // TODO: describe what '>' means
  bool operator>(TraceType& rhs);

  // Method to check if two trace types are the same
  bool operator==(TraceType& rhs);

  // Method to build colour amplitude from trace type
  // (A start index is needed)
  CAmplitude build_ca(size_t start_ind);

  // Method to print string expression of trace type to terminal
  void print();
};

//*********************************************************************************************************
//
// Class TraceVec
//
//********************************************************************************************************* 

class TraceVec {
  
  vector<TraceType> m_tr_vec;
    
 public:
  TraceVec(TraceType tr = TraceType({},0,0));
  ~TraceVec();
  void push_back(TraceType tr);
  TraceType& at(size_t i);
  TraceType at(size_t i) const;
  vector<TraceVec> add_one_gluon(size_t g_ind);
  vector<TraceVec> conjugates();
  vector<size_t> get_indices();
  size_t no_groups();
  //        bool is_tree_level();
  bool has_sg();
  void order();
  bool operator>(TraceVec& rhs);
  bool operator==(TraceVec& rhs);
  CAmplitude build_ca();
  void print();
};

class TraceBasis : public CBasis {
    
  vector<TraceVec> m_tr_basis;
  vector<size_t> m_g_indices, m_q_indices, m_qb_indices;
  size_t m_ng, m_nqp;
  
  void make_perms();
  void make_ca_basis();
  
  void remove_sg();
  void remove_conj();
  void normal_order();
  
 public:
  TraceBasis(size_t n_g, size_t n_qp);
  ~TraceBasis();
  
  friend class MultipletBasis;
};

#endif
