// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#ifndef C_BASIS_H
#define C_BASIS_H

#include "colourtools.h"
#include "c_matrix.h"

//*********************************************************************************************************
//
// Class c_basis
//
//********************************************************************************************************* 

class c_basis { 
 protected:
  // Information about process
  process m_process;

  // Vector of colour amplitudes to represent colour basis
  vector<c_amplitude> m_ca_basis;
  
  // Dimension of the basis
  size_t m_dim;

  // Conversion factor to colour flow basis
  double m_confact;

  // All amplitude permutations (representing colour flow in the process)
  // (needed to get colour ordered amplitudes from a ME generator)
  vector<vector<size_t>> m_amp_perms;

  // Normalisation factors of basis vectors
  vector<double> m_normalisations;

  // The scalar product (soft) matrix
  c_matrix m_smat;

  // The colour change matrices (colour soft anomalous dimension matrices)
  vector<c_matrix> m_ccmats;

  // The basis type
  // 0: general basis from file
  // 1: multiplet basis
  // 2: trace basis
  // 3: adjoint basis
  size_t m_btype;    
 
 public:
  // Method to normalise basis vectors
  void normalise(bool to_LC = false);
  
  // Method to get dimension of basis
  size_t dim();

  // Method to print basis in string representation to terminal
  void print();

  // Method to print all information to a text file
  // if to_LC = true, an additional statement will be printed
  // at the beginning of the file
  void print_to_file(string filename = "", bool to_LC = false);
  
  // Method to calculate and store the soft matrix
  // if to_LC = true, a unit matrix will be constructed (no actual computation)
  c_matrix sm(bool to_LC = false);

  // Method to calculate a colour change matrix between two legs
  // if to_LC = true, all scalar products will be evaluated at leading colour order
  c_matrix ccm(size_t lno1, size_t lno2, bool to_LC = false);

  // Method to calculate all colour change matrices for insertions between all legs
  // if to_LC = true, all scalar products will be evaluated at leading colour order
  vector<c_matrix> ccms(bool to_LC = false);

  // Method to check colour conservation for colour change matrices
  // Note: does not make sense for LC computations
  bool check_colourcons();
};

#endif
