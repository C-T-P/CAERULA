// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#ifndef MULTIPLET_BASIS_H
#define MULTIPLET_BASIS_H

#include "gen_basis.h"

using namespace std;

class multiplet_basis : public gen_basis {
    
  matrix m_bcm;
  size_t m_ng, m_nqp;
  
 public:
  multiplet_basis(size_t n_g, size_t n_qp);
  ~multiplet_basis();
    
  matrix bcm();
  
  friend class c_basis;
};

#endif
