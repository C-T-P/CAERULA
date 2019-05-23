// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#ifndef MULTIPLETBASIS_H
#define MULTIPLETBASIS_H

#include "GenBasis.h"

using namespace std;

class MultipletBasis : public GenBasis {
    
  matrix m_bcm;
  size_t m_ng, m_nqp;
  
 public:
  MultipletBasis(size_t n_g, size_t n_qp);
  ~MultipletBasis();
    
  matrix bcm();
  
  friend class CBasis;
};

#endif
