#ifndef C_BASIS_H
#define C_BASIS_H

#include "I3NSERT.h"
#include "c_matrix.h"
#include "colourtools.h"

class c_basis {
    protected:
        process m_process;
        vector<c_amplitude> m_ca_basis;
        vector<vector<size_t>> m_amp_perms;
        vector<double> m_normalisations;
        size_t m_dim;
    
    public:
        void normalise();
    
        size_t dim();
        void print();
    
        c_matrix sm();
        c_matrix ccm(size_t lno1, size_t lno2);
        vector<c_matrix> get_ccms();
};

#endif
