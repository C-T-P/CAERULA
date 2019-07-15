// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.                                                                               

#ifndef SPECTRUM_H
#define SPECTRUM_H

#include<cstring>
#include<iostream>

using namespace std;

//*********************************************************************************************************
//
// Base class Spectrum
//
//*********************************************************************************************************

class Spectrum {

  // Run option
  int m_runoption;

  // Number of gluons and quark pairs
  size_t m_n_g, m_n_qp;

  // Evaluate expressions in large-NC limit?
  bool m_largeNC;

  // Input expression
  string m_expr;

  // Input filename
  string m_in_filename;

  // Output filename
  string m_out_filename;

  // Save output file?
  bool m_no_outfile;

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

  // Print banner to terminal
  void print_banner();

  // Print help/information message
  void print_help_message();

  // Run error due to contrary input
  void run_error();

  // Basis error as two different bases were specified
  void basis_error();

 public:
  // Default constructor
  Spectrum() {
    print_banner();
    m_runoption = 0; m_n_g = 0; m_n_qp = 0; m_largeNC = false; 
    m_expr = ""; m_in_filename = ""; m_out_filename = ""; m_no_outfile = false;
    m_traceBasis = true; m_adjointBasis = false; m_multipletBasis = false;
    m_reduceDim = true; m_norm_basis = true; m_construct_bcm = false; 
    m_calcDet = false; m_verbose = 1;
  }

  // Initialise from vector of arguments
  void init(int argc, char **argv);

  // Start automated colour computation
  bool start();
  
  // Evaluate a given colour amplitude
  bool evaluate_ca(string expr = "");

  // Simplify a given colour amplitude
  bool simplify_ca(string expr = "");

  //----------------------------------------------------------------------------------
  // Setters
  //----------------------------------------------------------------------------------
  void set_verbose(int verbosity) {
    if (verbosity < 0) m_verbose = 0;
    else if (verbosity > 3) m_verbose = 3;
    else m_verbose = verbosity;
  };
  void set_ng(size_t ng) { m_n_g = ng; };
  void set_nqp(size_t nqp) { m_n_qp = nqp; };
  void set_largeNC(bool isLC) { m_largeNC = isLC; };
  void set_in_filename(string filename) { m_in_filename = filename; };
  void set_out_filename(string filename) { m_out_filename = filename; };
  void print_outfile(bool print_file) { m_no_outfile = !print_file; };
  void use_trace_basis(bool use_tr_basis) { m_traceBasis = use_tr_basis; };
  void use_adjoint_basis(bool use_adj_basis) { m_adjointBasis = use_adj_basis; };
  void use_multiplet_basis(bool use_orth_basis) { m_multipletBasis = use_orth_basis; };
  void reduce_trace_dim(bool reduce_dim) { m_reduceDim = reduce_dim; };
  void normalise_basis(bool is_normalised) { m_norm_basis = is_normalised; };
  void construct_bcm(bool constr_bcm) { m_construct_bcm = constr_bcm; };
  void set_calc_det(bool calcDet) { m_calcDet = calcDet; };

  //----------------------------------------------------------------------------------
  // Getters
  //----------------------------------------------------------------------------------
  int get_verbose(int verbosity) { return m_verbose; };
  size_t get_ng() { return m_n_g; };
  size_t get_nqp() { return m_n_qp; };
  bool is_largeNC() { return m_largeNC; };
  string get_in_filename() { return m_in_filename; };
  string get_out_filename() { return m_out_filename; };
  bool print_outfile() { return !m_no_outfile; };
  bool use_trace_basis() { return m_traceBasis; };
  bool use_adjoint_basis() { return m_adjointBasis; };
  bool use_multiplet_basis() { return m_multipletBasis; };
  bool reduce_trace_dim() { return m_reduceDim; };
  bool normalise_basis() { return m_norm_basis; };
  bool construct_bcm() { return m_construct_bcm; };
  bool get_calc_det() { return m_calcDet; };

};

#endif
