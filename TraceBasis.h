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
  //   If combineConj = true, the CAmplitude will be combined with its conjugate
  //   Note: a start index is needed
  CAmplitude build_ca(size_t start_ind, bool combineConj);

  // Method to print string expression of trace type to terminal
  void print();
};

//*********************************************************************************************************
//
// Class TraceVec
//
//********************************************************************************************************* 

class TraceVec {

  // Vector of TraceType
  vector<TraceType> m_tr_vec;
    
 public:
  
  // Constructor
  TraceVec(TraceType tr = TraceType({},0,0));

  // Destructor
  ~TraceVec();

  // Push new TraceType to the back
  void push_back(TraceType tr);

  // Write element
  TraceType& at(size_t i);

  // Access element
  TraceType at(size_t i) const;

  // Add one gluonic index in all possible ways
  vector<TraceVec> add_one_gluon(size_t g_ind);

  // Get all conjugate index permutations
  vector<TraceVec> conjugates();

  // Get all indices in order
  vector<size_t> get_indices();

  // Get the number of groups
  size_t no_groups();

  // Check if this index combination corresponds to a tree level amplitude
  //        bool is_tree_level();

  // Check if this index combination has a single gluon attached to a quark ring
  bool has_sg();

  // Order the groups in the index permutation
  void order();
  
  // Comparison >
  bool operator>(TraceVec& rhs);

  // Comparison ==
  bool operator==(TraceVec& rhs);

  // Build CAmplitude for this index permutation 
  CAmplitude build_ca(bool reduceDim);

  // Print index permutation to terminal
  void print();
};

//*********************************************************************************************************
//
// Class TraceBasis
//
//********************************************************************************************************* 

class TraceBasis : public CBasis {
  
  // Basis in terms of TraceVecs
  vector<TraceVec> m_tr_basis;

  // Index permutations of basis vectors
  vector<size_t> m_g_indices, m_q_indices, m_qb_indices;

  // Number of gluons and quark pairs
  size_t m_ng, m_nqp;
  
  // Save all index permutations in m_amp_perms
  void make_perms();

  // Build basis in terms of colour amplitudes
  //  if reduceDim = true, the dimension of the basis will be reduced by a 
  //  linear combination of conjugate basis vectors
  void make_ca_basis();
  
  // Remove quark rings with single gluons attached
  void remove_sg();

  // Remove conjugate basis vectors
  void remove_conj();

  // Normal order the basis
  void normal_order();

  // Do we want to reduce the dimensionality of the basis?
  bool m_reduce_dim = true;
  
 public:

  // Constructor
  TraceBasis(size_t n_g, size_t n_qp, bool reduceDim = true);

  // Destructor
  ~TraceBasis();
  
  friend class MultipletBasis;
};

#endif
