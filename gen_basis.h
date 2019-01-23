// Copyright (C) 2018 Christian T Preuss 
// This file is part of Spectrum.                                                                                                                                                                                                                                              
//                                                                                                                                                                                                                                                                             
// Spectrum is free software: you can redistribute it and/or modify                                                                                                                                                                                                            
// it under the terms of the GNU General Public License as published by                                                                                                                                                                                                        
// the Free Software Foundation, either version 3 of the License, or                                                                                                                                                                                                           // any later version.     

#ifndef GEN_BASIS_H
#define GEN_BASIS_H

#include "c_basis.h"

using namespace std;

class gen_basis : public c_basis {
    
    public:
        gen_basis(string filename);
        ~gen_basis();
};

#endif
