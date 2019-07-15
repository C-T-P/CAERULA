// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#include "Spectrum.h"

int main(int argc, char **argv) {
  Spectrum  spectrum;
  spectrum.init(argc,argv);
  return spectrum.start();
}
