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
  
  // Normalise basis?
  bool m_norm_basis;

  // Construct basis change matrix to trace basis?
  bool m_construct_bcm;

  // Verbosity
  // 0: only banner
  // 1: only errors and option messages
  // 2: debug info
  // 3: full information
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
    m_norm_basis = true; m_construct_bcm = false; m_verbose = 1;
  }

  // Constructor from vector of arguments
  Spectrum(int argc, char **argv);

  // Start colour computation
  bool start();

  // TODO: setters for process, basis, etc.
};

#endif
