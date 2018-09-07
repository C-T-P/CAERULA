#ifndef MULTIPLET_BASIS_H
#define MULTIPLET_BASIS_H

#include "c_basis.h"

using namespace std;

class multiplet_basis : public c_basis {
    
    matrix m_bcm;
    size_t m_ng, m_nqp;
    
    public:
        multiplet_basis(size_t n_g, size_t n_qp);
        ~multiplet_basis();
    
        matrix bcm();
    
        friend class c_basis;
};

#endif
