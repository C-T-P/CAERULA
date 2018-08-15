#include "I3NSERT.h"

// construct insertion operator for gluon insertion between leg lno1 and lno2 as a colour term
c_amplitude construct_insertion_op(process proc, size_t lno1, size_t lno2) {
    c_term ins_op;
    complex<double> prefactor(1.);
    
    if (proc.leg(lno1).second=="g") {
        ins_op.push_back(antisymmetric(lno1,1001,lno1+2000));
        prefactor*=complex<double>(0.,1.);
    }
    else if ((proc.leg(lno1).second=="q" and proc.is_in_leg(lno1)) or (proc.leg(lno1).second=="qb" and !proc.is_in_leg(lno1))) {
        ins_op.push_back(fundamental(1001,lno1+2000,lno1));
        prefactor*=-1;
    }
    else if ((proc.leg(lno1).second=="q" and !proc.is_in_leg(lno1)) or (proc.leg(lno1).second=="qb" and proc.is_in_leg(lno1))) ins_op.push_back(fundamental(1001,lno1,lno1+2000));
    else {
        cerr << "Error constructing the insertion operator between leg " << lno1 << " and leg " << lno2 << ": leg " << lno1 << " is not a quark, anti quark, or gluon." << endl;
        exit(EXIT_FAILURE);
    }
    
    if (proc.leg(lno2).second=="g") {
        ins_op.push_back(antisymmetric(lno2,1001,lno2+2000));
        prefactor*=complex<double>(0.,1.);
    }
    else if ((proc.leg(lno2).second=="q" and proc.is_in_leg(lno2)) or (proc.leg(lno2).second=="qb" and !proc.is_in_leg(lno2))) {
        ins_op.push_back(fundamental(1001,lno2+2000,lno2));
        prefactor*=-1;
    }
    else if ((proc.leg(lno2).second=="q" and !proc.is_in_leg(lno2)) or (proc.leg(lno2).second=="qb" and proc.is_in_leg(lno2))) ins_op.push_back(fundamental(1001,lno2,lno2+2000));
    else { 
        cerr << "Error constructing the insertion operator between leg " << lno1 << " and leg " << lno2 << ": leg " << lno2 << " is not a quark, anti quark, or gluon." << endl;
        exit(EXIT_FAILURE);
    }
    
    for (unsigned int lno(1);lno<=proc.no_of_legs();lno++) {
        if (lno!=lno1 and lno!=lno2) {
            if (proc.leg(lno).second=="g") ins_op.push_back(delta(lno+2000,lno,true));
            else ins_op.push_back(delta(lno+2000,lno,false));
        }
    }
    
    ins_op.cnumber(prefactor);

    return c_amplitude(ins_op);
}

// TODO find algorithm to determine which insertions are equivalent
