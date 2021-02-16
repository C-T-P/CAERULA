// Copyright (C) 2021 Christian T Preuss 
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#ifndef INSERT_H
#define INSERT_H

#include "CBasis.H"
#include "Colourtools.H"

namespace SPECTRUM {

CAmplitude construct_insertion_op(Process proc, size_t lno1, size_t lno2);

}

#endif
