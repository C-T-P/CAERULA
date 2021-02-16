// Copyright (C) 2021 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.                                                                               

#ifndef SPECTRUM_H
#define SPECTRUM_H

// Standard includes.
#include<cstring>
#include<iostream>

// Spectrum includes.
#include "CMatrix.h"
#include "CBasis.h"
#include "Colourtools.h"
#include "TraceBasis.h"
#include "FBasis.h"
#include "GenBasis.h"
#include "MultipletBasis.h"

using namespace std;

namespace SPECTRUM {

//*****************************************************************************
//
// Base class Spectrum.
//
//*****************************************************************************

class Spectrum {

  // Number of gluons and quark pairs
  size_t m_n_g, m_n_qp;

  // Evaluate expressions in large-NC limit?
  bool m_largeNC;

  // Input expression
  string m_expr;

  // Output filename
  string m_out_filename;

  // Colour basis information
  bool m_traceBasis;
  bool m_adjointBasis;
  bool m_multipletBasis;
  bool m_reduceDim;
  
  // Normalise basis?
  bool m_norm_basis;

  // Construct basis change matrix to trace basis?
  bool m_construct_bcm;

  // Calculate determinant of soft matrix?
  bool m_calcDet;

  // Verbosity
  //   0: only banner
  //   1: only errors and option messages
  //   2: debug info
  //   3: full information
  int m_verbose;

  // Pointers
  CBasis* m_basis;  

  //----------------------------------------------------------------------------
  // Private Member Functions
  //----------------------------------------------------------------------------

  // Print banner to terminal
  void print_banner();

 public:
  // Default constructor
  Spectrum() {
    print_banner();
    m_n_g = 0; m_n_qp = 0; m_largeNC = false;
    m_traceBasis = true; m_adjointBasis = false; m_multipletBasis = false;
    m_reduceDim = true; m_norm_basis = true; m_construct_bcm = false; 
    m_calcDet = false; m_verbose = 1;
    m_basis = NULL;}

  // Evaluate a given colour amplitude
  bool evaluate_ca(string expr = "");

  // Simplify a given colour amplitude
  bool simplify_ca(string expr = "");

  // Read basis from file
  bool read_basis(string filename);

  // Construct basis from number of quark pairs and gluons
  bool construct_basis(int nqp, int ng);

  // Calculate the soft matrix
  CMatrix calculate_soft_matrix();
  
  // Calculate colour correlators
  vector<CMatrix> calculate_colour_correlators();

  // Save everything to a file
  void save_to_file(string filename = "");
  
  //----------------------------------------------------------------------------
  // Setters
  //----------------------------------------------------------------------------
  void set_verbose(int verbosity) {
    if (verbosity < 0) m_verbose = 0;
    else if (verbosity > 3) m_verbose = 3;
    else m_verbose = verbosity;};
  void set_ng(int n_g) { m_n_g = n_g; };
  void set_nqp(int n_qp) { m_n_qp = n_qp; };
  void set_largeNC(bool isLC) { m_largeNC = isLC; };
  void set_use_trace_basis(bool use_tr_basis) { m_traceBasis = use_tr_basis; };
  void set_use_adjoint_basis(bool use_adj_basis) { 
    m_adjointBasis = use_adj_basis;};
  void set_use_multiplet_basis(bool use_orth_basis) { 
    m_multipletBasis = use_orth_basis;};
  void set_reduce_trace_dim(bool reduce_dim) { m_reduceDim = reduce_dim; };
  void set_normalise_basis(bool is_normalised) { m_norm_basis = is_normalised; };
  void set_construct_bcm(bool constr_bcm) { m_construct_bcm = constr_bcm; };
  void set_calc_det(bool calcDet) { m_calcDet = calcDet; };

  //---------------------------------------------------------------------------
  // Getters
  //---------------------------------------------------------------------------
  int get_verbose(int verbosity) { return m_verbose; };
  size_t get_ng() { return m_n_g; };
  size_t get_nqp() { return m_n_qp; };
  bool is_largeNC() { return m_largeNC; };
  bool get_use_trace_basis() { return m_traceBasis; };
  bool get_use_adjoint_basis() { return m_adjointBasis; };
  bool get_use_multiplet_basis() { return m_multipletBasis; };
  bool get_reduce_trace_dim() { return m_reduceDim; };
  bool get_normalise_basis() { return m_norm_basis; };
  bool get_construct_bcm() { return m_construct_bcm; };
  bool get_calc_det() { return m_calcDet; };

  //---------------------------------------------------------------------------
  // Printers
  //---------------------------------------------------------------------------
  void print_settings();
  void print_basis();
};

// Some small helper functions.
string bool2str(bool in);
string int2str(int in, int pad=0);

}

#endif
