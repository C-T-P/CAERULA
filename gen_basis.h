#ifndef GEN_BASIS_H
#define GEN_BASIS_H

#include "c_basis.h"

using namespace std;

class gen_basis : public c_basis {
    
    public:
        gen_basis(string filename);
        ~gen_basis();
    
        friend class multiplet_basis;
};

#endif
