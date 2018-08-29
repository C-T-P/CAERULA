#ifndef C_BASIS_H
#define C_BASIS_H

#include "I3NSERT.h"
#include "c_matrix.h"
#include "colourtools.h"

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
        vector<c_matrix> get_ccms(size_t up_to_NC = INT_MAX);
};

#endif
