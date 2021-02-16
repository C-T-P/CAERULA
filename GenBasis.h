// Copyright (C) 2021 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#ifndef GENBASIS_H
#define GENBASIS_H

#include "CBasis.h"

using namespace std;

namespace SPECTRUM {

//*****************************************************************************
//
// Class GenBasis.
//
//*****************************************************************************

class GenBasis : public CBasis {
  
 public:
  GenBasis(string filename);
  ~GenBasis();
};

}

#endif
