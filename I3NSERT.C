#include "colourtools.h"
#include "Main.h"
#include "CONTRACT.h"
#include "I3NSERT.h"

// construct insertion operator for gluon insertion between leg lno1 and lno2 as a colour term
colour_term construct_insertion_op(process m_process, unsigned int lno1, unsigned int lno2) {
    three_ind symmetric;
    three_ind antisymmetric;
    three_ind fundamental;
    two_ind kronecker;
    complex<double> prefactor(1.);
    
    if (m_process.leg(lno1).second=="g") {
        antisymmetric.set_indices(lno1,1001,lno1+2000);
        prefactor*=complex<double>(0.,1.);
    }
    else if ((m_process.leg(lno1).second=="q" and m_process.is_in_leg(lno1)) or (m_process.leg(lno1).second=="qb" and !m_process.is_in_leg(lno1))) { 
        fundamental.set_indices(1001,lno1+2000,lno1);
        prefactor*=-1;
    }
    else if ((m_process.leg(lno1).second=="q" and !m_process.is_in_leg(lno1)) or (m_process.leg(lno1).second=="qb" and m_process.is_in_leg(lno1))) fundamental.set_indices(1001,lno1,lno1+2000);
    else {
        cerr << "Error constructing the insertion operator between leg " << lno1 << " and leg " << lno2 << ": leg " << lno1 << " is not a quark, anti quark, or gluon." << endl;
        exit(EXIT_FAILURE);
    }
    
    if (m_process.leg(lno2).second=="g") { 
        antisymmetric.set_indices(lno2,1001,lno2+2000);
        prefactor*=complex<double>(0.,1.);
    }
    else if ((m_process.leg(lno2).second=="q" and m_process.is_in_leg(lno2)) or (m_process.leg(lno2).second=="qb" and !m_process.is_in_leg(lno2))) {
        fundamental.set_indices(1001,lno2+2000,lno2);
        prefactor*=-1;
    }
    else if ((m_process.leg(lno2).second=="q" and !m_process.is_in_leg(lno2)) or (m_process.leg(lno2).second=="qb" and m_process.is_in_leg(lno2))) fundamental.set_indices(1001,lno2,lno2+2000);
    else { 
        cerr << "Error constructing the insertion operator between leg " << lno1 << " and leg " << lno2 << ": leg " << lno2 << " is not a quark, anti quark, or gluon." << endl;
        exit(EXIT_FAILURE);
    }
     
    bool gluonic(false);
    for (unsigned int lno(1);lno<=m_process.no_of_legs();lno++) {
        if (lno!=lno1 and lno!=lno2) {
            if (m_process.leg(lno).second=="g") gluonic=true;
            else gluonic=false;
            kronecker.set_indices(lno+2000,lno,gluonic);
        }
    }

    colour_term insertion_op;
    insertion_op.sym.push_back(symmetric);
    insertion_op.asym.push_back(antisymmetric);
    insertion_op.fund.push_back(fundamental);
    insertion_op.kron.push_back(kronecker);
    insertion_op.pref.push_back(prefactor);
    insertion_op.NC_ctr.push_back(0);
    
    return insertion_op;
}
