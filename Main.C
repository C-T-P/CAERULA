// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#include<fstream>
#include<iomanip>
#include<climits>
#include<ctime>
#include "c_matrix.h"
#include "colourtools.h"
#include "trace_basis.h"
#include "f_basis.h"
#include "gen_basis.h"
#include "multiplet_basis.h"

void run_error();
void print_help_message();
void print_startup_message();

int main(int argc, char **argv) {
    print_startup_message();

    c_basis* basis=NULL;
    size_t n_g(0), n_qp(0);
    string expr, filename;
    bool adjoint_basis(false), ortho_basis(false), norm_basis(true), construct_bcm(false), 
      print_to_console(false), no_output(false), is_LC(false);
    
    // read in run parameters
    if (argc==1) print_help_message();
    int runopt(0);
    for(int i(1); i<argc; ++i) {
        if (strcmp(argv[i], "-help")==0 or strcmp(argv[i], "-h")==0) {
            print_help_message();
        }
        else if (strcmp(argv[i], "-evaluate")==0 or strcmp(argv[i], "-e")==0) {
            if (runopt==0) runopt=1;
            else run_error();
            
            expr=argv[i+1];
            
            i++;
        }
        else if (strcmp(argv[i], "-simplify")==0 or strcmp(argv[i], "-s")==0) {
            if (runopt==0) runopt=2;
            else run_error();
            
            expr=argv[i+1];
            
            i++;
        }
        else if (strcmp(argv[i], "-f")==0 or strcmp(argv[i], "-file")==0) {
            if (runopt==0) runopt=3;
            else run_error();
            
            filename=argv[i+1];
            
            i++;
        }
        else if (strcmp(argv[i],"-ng")==0) {
            if (runopt==0 or runopt==4) runopt=4;
            else run_error();
            
            n_g=stoi(argv[i+1]);
            
            i++;
        }
        else if (strcmp(argv[i],"-nqp")==0) {
            if (runopt==0 or runopt==4) runopt=4;
            else run_error();
            
            n_qp=stoi(argv[i+1]);
            
            i++;
        }
        else if (strcmp(argv[i],"-adj")==0) {
            if (!ortho_basis) adjoint_basis=true;
            else run_error();
        }
        else if (strcmp(argv[i],"-multiplet")==0) {
            if (!adjoint_basis) ortho_basis=true;
            else run_error();
        }
        else if (strcmp(argv[i],"-dnorm")==0) {
            norm_basis=false;
        }
        else if (strcmp(argv[i],"-bcm")==0) {
            construct_bcm=true;
        }
        else if (strcmp(argv[i],"-printall")==0) {
            print_to_console=true;
            cout<<"Print to console set to true."<<endl;
        }
        else if (strcmp(argv[i],"-noout")==0) {
            no_output=true;
            cout<<"No output file will be saved."<<endl;
        }
	else if (strcmp(argv[i],"-LC")==0) {
            is_LC=true;
        }
    }
    
    // print warning that only LC will be used
    if (is_LC)
      cout<<"\n\033[1;31mLeading Colour enabled.\033[0m\n" << endl;
    
    // perform colour calculations depending on the given input
    switch (runopt) {
        case 1: {
            c_amplitude ca(expr);
            ca.print();
	    if (is_LC) ca.evaluate_LC();
	    else ca.evaluate();
            cout<<"= "<<ca.result()<<endl;
            
            return 0;
        }
        case 2: {
            c_amplitude ca(expr);
            ca.print();
            ca.simplify();
            cout<<"= ";
            ca.print();
            
            return 0;
        }
        case 3: {
            cout<<"Will construct colour basis from file "<<filename<<"."<<endl;
            cout<<"\nNOTE: amplitude permutations cannot be computed for bases read from files!"<<endl;
            
            basis = new gen_basis(filename);
            
            break;
        }
        case 4: {
            clock_t t(clock());
            if (ortho_basis) {
                cout<<"Will read in multiplet basis for "<<n_g<<" gluons and "<<n_qp<<" quark pairs."<<endl;
                
                multiplet_basis* m_basis = new multiplet_basis(n_g, n_qp);
                
                if (construct_bcm) {
                    m_basis->bcm();
                }
                
                basis = m_basis;
            }
            else if (adjoint_basis) {
                cout<<"Will construct adjoint basis for "<<n_g<<" gluons."<<endl;
                basis = new f_basis(n_g);
            }
            else {
                cout<<"Will construct trace basis for "<<n_g<<" gluons and "<<n_qp<<" quark pairs."<<endl;
                basis = new trace_basis(n_g, n_qp);
            }
            
            t=clock()-t;
            cout<<"Computation time for basis construction: "<<(float)t/CLOCKS_PER_SEC<<"s."<<endl;
            
            break;
        }
        default: {
            cerr<<"Error: Couldn't perform colour calculations! Check -h for usage instructions."<<endl;
            exit(EXIT_FAILURE);
        }
    }

    
    if (norm_basis and !construct_bcm) basis->normalise();
    if (print_to_console) {
        cout<<endl;
        if (construct_bcm or norm_basis) cout<<"Normalised ";
        cout<<"Basis Vectors:"<<endl;
        basis->print();
    }
        
    clock_t t(clock());
    cout<<"\nCalculating the soft matrix..."<<endl;
    c_matrix soft_matrix(basis->sm(is_LC));
    if (print_to_console) {
        cout<<"\nSoft Matrix:"<<endl;
        soft_matrix.print();
    }
    t=clock()-t;
    cout<<"Computation time: "<<(float)t/CLOCKS_PER_SEC<<"s."<<endl;

    t=clock();
    cout<<"\nCalculating the colour change matrices..."<<endl;
    vector<c_matrix> cc_mats(basis->ccms(is_LC));
    if (print_to_console) {
        cout<<"\nColour Change Matrices:"<<endl;
        for (auto& ccm : cc_mats) {
            ccm.print();
            cout<<endl;
        }
    }
    t=clock()-t;
    cout<<"Computation time: "<<(float)t/CLOCKS_PER_SEC<<"s."<<endl;
    
    cout<<"\nIs colour conserved? ";
    if(basis->check_colourcons()) cout<<"Yes!"<<endl;
    else cout<<"NO!"<<endl;
    
    if (!no_output) {
        cout<<"\nPrinting to file..."<<endl;
        basis->print_to_file("", is_LC);
        cout<<"Done!"<<endl;
    }
    
    return 0;
}

void run_error() {
    cerr<<"Please specify EITHER a colour term to be simplified (-s)/evaluated (-e) OR perform a colour space calculation by specifying a basis through a file name (-f) OR the process by the number of quark pairs (-nqp) and number of gluons (-ng). See also the help menu for more information (-h)."<<endl;
    exit(EXIT_FAILURE);
}
void print_help_message() {
    cout<<"\nOptions:"<<endl;
    cout<<"\t-h, -help:\t\tPrint this message and exit."<<endl;
    cout<<endl;
    cout<<"\t-e, -evaluate:\t\tEvaluate colour amplitude in which all indices must be contracted."<<endl;
    cout<<"\t-s, -simplify:\t\tSimplify colour term by replacing internal quark and gluon rings. Returns colour amplitude."<<endl;
    cout<<"\t\t\t\tReturns complex number or NAN if not all indices are contracted."<<endl;
    cout<<endl;
    cout<<"\t-f, -file:\t\tRead process and colour basis from file and compute the soft matrix and all colour change matrices."<<endl;
    cout<<"\t-ng:\t\t\tSpecify number of gluons to construct the trace basis and compute the soft matrix and all colour change matrices."<<endl;
    cout<<"\t-nqp:\t\t\tSpecify number of quark pairs to construct the colour flow/trace basis and compute the soft matrix and all colour change matrices."<<endl;
    cout<<"\t-adj:\t\t\tBuild adjoint basis (f-basis) instead of trace basis.\n\t\t\t\tWorks only for pure gluon processes with ng>=3."<<endl;
    cout<<"\t-multiplet:\t\tBuild multiplet basis (orthogonal basis).\n\t\t\t\tWorks only if a precalculated multiplet basis for this process is provided."<<endl;
    cout<<endl;
    cout<<"\t-bcm:\t\t\tBuild basis change matrix from trace basis to multiplet basis.\n\t\t\t\tWorks only together with the -multiplet option and if a precalculated multiplet basis for the process is provided."<<endl;
    //cout<<"\t-NC:\t\t\tSpecify order in 1/NC to which all colour products shall be evaluated."<<endl;
    cout<<"\t-dnorm:\t\t\tDeactivates normalisation of basis vectors."<<endl;
    cout<<endl;
    exit(EXIT_SUCCESS);
}
void print_startup_message() {
    cout<<endl;
    cout<<"============================================================"<<endl;
    cout<<"============            \e[31mSpe\e[32mctr\e[34mum\e[0m             ==============="<<endl;
    cout<<"============---------------------------------==============="<<endl;
    cout<<"============  Calculations in Colour Space   ==============="<<endl;
    cout<<"============================================================"<<endl;
    cout<<endl;
    cout<<"Spectrum version 1.0, Copyright (C) 2018 Christian T Preuss"<<endl;
    cout<<endl;
    cout<<"=== Author =================================================\n"<<endl;
    cout<<"\tChristian Tobias Preuss\n\tSchool of Physics, Monash University,\n\t3800 Melbourne, Australia\n\temail: christian.preuss@monash.edu\n"<<endl;
    cout<<"============================================================\n"<<endl;
}
