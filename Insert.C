// Copyright (C) 2021 Christian T Preuss 
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#include "Insert.H"

namespace SPECTRUM {

// Construct Insertion operator for gluon Insertion between leg lno1 and lno2 as a colour term.
CAmplitude construct_insertion_op(process proc, size_t lno1, size_t lno2) {
    CTerm ins_op;
    ColourFactor prefactor = ColourFactor(1., 0, 0, 0, 0);

    if (proc.leg(lno1)=="g") {
        ins_op.push_back(Antisymmetric(lno1,1001,lno1+2000));
        prefactor*=complex<double>(0.,1.);
    }
    else if ((proc.leg(lno1)=="q" and proc.is_in_leg(lno1)) or (proc.leg(lno1)=="qb" and !proc.is_in_leg(lno1))) {
        ins_op.push_back(Fundamental(1001,lno1+2000,lno1));
        prefactor*=-1;
    }
    else if ((proc.leg(lno1)=="q" and !proc.is_in_leg(lno1)) or (proc.leg(lno1)=="qb" and proc.is_in_leg(lno1))) ins_op.push_back(Fundamental(1001,lno1,lno1+2000));
    else {
        cerr << "Error constructing the insertion operator between leg " << lno1 << " and leg " << lno2 << ": leg " << lno1 << " is not a quark, anti quark, or gluon." << endl;
        exit(EXIT_FAILURE);
    }
    
    if (proc.leg(lno2)=="g") {
        ins_op.push_back(Antisymmetric(lno2,1001,lno2+2000));
        prefactor*=complex<double>(0.,1.);
    }
    else if ((proc.leg(lno2)=="q" and proc.is_in_leg(lno2)) or (proc.leg(lno2)=="qb" and !proc.is_in_leg(lno2))) {
        ins_op.push_back(Fundamental(1001,lno2+2000,lno2));
        prefactor*=-1;
    }
    else if ((proc.leg(lno2)=="q" and !proc.is_in_leg(lno2)) or (proc.leg(lno2)=="qb" and proc.is_in_leg(lno2))) ins_op.push_back(Fundamental(1001,lno2,lno2+2000));
    else { 
        cerr << "Error constructing the insertion operator between leg " << lno1 << " and leg " << lno2 << ": leg " << lno2 << " is not a quark, anti quark, or gluon." << endl;
        exit(EXIT_FAILURE);
    }
    
    for (unsigned int lno(1);lno<=proc.no_of_legs();lno++) {
        if (lno!=lno1 and lno!=lno2) {
            if (proc.leg(lno)=="g") ins_op.push_back(Delta(lno+2000,lno,true));
            else ins_op.push_back(Delta(lno+2000,lno,false));
        }
    }
    
    ins_op.set_cnumber(prefactor);

    return CAmplitude(ins_op);
}

}
