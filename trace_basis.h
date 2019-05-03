// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#ifndef TRACE_BASIS_H
#define TRACE_BASIS_H

#include "c_basis.h"
#include "colourtools.h"

using namespace std;

class trace_t {
  // The indices of the quark and antiquark
  size_t m_qb, m_q;

  // The indices of the gluons
  vector<size_t> m_g;

 public:
  // Default constructor
  trace_t(vector<size_t> g_inds = {}, size_t q_ind = 0, size_t qb_ind = 0);
  
  // Constructor from quark and antiquark indices
  trace_t(size_t q_ind, size_t qb_ind);

  //Destructor
  ~trace_t();
  
  // Method to add one gluonic index
  vector<trace_t> add_one_gluon(size_t g_ind);

  // Method to get all indices
  vector<size_t> get_indices();

  // Method to compute conjugate
  trace_t conj();

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
  bool operator>(trace_t& rhs);

  // Method to check if two trace types are the same
  bool operator==(trace_t& rhs);

  // Method to build colour amplitude from trace type
  // (A start index is needed)
  c_amplitude build_ca(size_t start_ind);

  // Method to print string expression of trace type to terminal
  void print();
};

class trace_vec {
  
  vector<trace_t> m_tr_vec;
    
 public:
  trace_vec(trace_t tr = trace_t({},0,0));
  ~trace_vec();
  void push_back(trace_t tr);
  trace_t& at(size_t i);
  trace_t at(size_t i) const;
  vector<trace_vec> add_one_gluon(size_t g_ind);
  vector<trace_vec> conjugates();
  vector<size_t> get_indices();
  size_t no_groups();
  //        bool is_tree_level();
  bool has_sg();
  void order();
  bool operator>(trace_vec& rhs);
  bool operator==(trace_vec& rhs);
  c_amplitude build_ca();
  void print();
};

class trace_basis : public c_basis {
    
  vector<trace_vec> m_tr_basis;
  vector<size_t> m_g_indices, m_q_indices, m_qb_indices;
  size_t m_ng, m_nqp;
  
  void make_perms();
  void make_ca_basis();
  
  void remove_sg();
  void remove_conj();
  void normal_order();
  
 public:
  trace_basis(size_t n_g, size_t n_qp);
  ~trace_basis();
  
  friend class multiplet_basis;
};

#endif
