// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.                                                                                                                                                                                                                                              
//                                                                                                                                                                                                                                                                             
// Spectrum is free software: you can redistribute it and/or modify                                                                                                                                                                                                            
// it under the terms of the GNU General Public License as published by                                                                                                                                                                                                        
// the Free Software Foundation, either version 3 of the License, or                                                                                                                                                                                                           // any later version.     

#ifndef C_BASIS_H
#define C_BASIS_H

#include "colourtools.h"
#include "c_matrix.h"

class c_basis { 
    protected:
        process m_process;
        vector<c_amplitude> m_ca_basis;
        size_t m_dim;
        double m_confact;
        vector<vector<size_t>> m_amp_perms;
        vector<double> m_normalisations;
        c_matrix m_smat;
        vector<c_matrix> m_ccmats;
        size_t m_btype;
        /*
         basis type:
         0: general basis from file
         1: multiplet basis
         2: trace basis
         3: adjoint basis
         */
    
    public:
        void normalise();
    
        size_t dim();
        void print();
        void print_to_file(string filename = "");
    
        c_matrix sm();
        c_matrix ccm(size_t lno1, size_t lno2, size_t up_to_NC = INT_MAX);
        vector<c_matrix> ccms(size_t up_to_NC = INT_MAX);
        bool check_colourcons();
};

#endif
