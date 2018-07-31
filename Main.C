#include<fstream>
#include<iomanip> 
#include<ctime>
#include "colourtools.h"
#include "c_matrix.h"
#include "BASIS.h"
#include "I3NSERT.h"
#include "CONTRACT.h"
#include "FCONTRACT.h"
#include "Main.h"

int main(int argc, char **argv) {
    int NC_order(INT_MAX), n_g(0), n_qp(0);
    vector<colour_term> basis;
    colour_term ct;
    process m_process;
    string filename, out_filename;
    
    // read in run parameters
    int runopt(0);
    for(int i(1); i<argc; ++i) {
        if (strcmp(argv[i], "-help")==0 or strcmp(argv[i], "-h")==0) {
            cout<<"\nOptions:"<<endl;
            cout<<"\t-h, -help:\t\tPrint this message and exit."<<endl;
            cout<<"\t-s, -simplify:\t\tSimplify colour term with open indices."<<endl;
            cout<<"\t\t\t\tReturns colour term."<<endl;
            cout<<"\t-e, -evaluate:\t\tEvaluate colour term in which all indices must be contracted."<<endl;
            cout<<"\t\t\t\tReturns complex number or NAN if not all indices are contracted."<<endl;
            cout<<"\t-f, -file:\t\tRead process and colour basis from file and compute the soft matrix and all colour change matrices."<<endl;
            cout<<"\t-ng:\t\t\tSpecify number of gluons to construct the trace basis and compute the soft matrix and all colour change matrices."<<endl;
            cout<<"\t-nqp:\t\t\tSpecify number of quark pairs to construct the colour flow/trace basis and compute the soft matrix and all colour change matrices."<<endl;
            cout<<"\t-NC:\t\t\tSpecify order in 1/NC to which all colour products shall be evaluated."<<endl;
            cout<<endl;
            exit(EXIT_SUCCESS);
        }
        else if (strcmp(argv[i], "-simplify")==0 or strcmp(argv[i], "-s")==0) {
            if (runopt==0) runopt=1;
            else run_error();
            
            string expr=argv[i+1];
            ct=decompose_terms(expr,m_process);
            
            i++;
        }
        else if (strcmp(argv[i], "-evaluate")==0 or strcmp(argv[i], "-e")==0) {
            if (runopt==0) runopt=2;
            else run_error();
            
            string expr=argv[i+1];
            ct=decompose_terms(expr,m_process);
            
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
        else if (strcmp(argv[i],"-NC")==0) {
            NC_order=stoi(argv[i+1]);
            i++;
        }
    }
    
    // print order in 1/NC to which all terms are evaluated
    cout<<"Order of 1/NC set to ";
    if (NC_order!=INT_MAX) cout << NC_order << "." << endl;
    else cout << "Infinity.\n" << endl;
    
    
    // perform colour calculations depending on the given input
    switch (runopt) {
        case 1: {
            cout<<ct.build_string()<<endl;
            simplify_colour_term(ct);
            cout<<"= "<<ct.build_string()<<endl;
            break;
        }
        case 2: {
            cout<<ct.build_string()<<endl;
            complex<double> d=fast_evaluate_colour_term_to_order(ct, NC_order);
            cout<<"= "<<d<<endl;
            break;
        }
        case 3: {
            cout<<"Will construct colour basis from file "<<filename<<"."<<endl;
            for (size_t lno(1);lno<=m_process.no_of_legs();lno++) out_filename+=m_process.leg(lno).second;
            basis=read_basis(filename, m_process);
            colour_calc(basis, m_process, NC_order, out_filename);
            break;
        }
        case 4: {
            if (n_g==0) cout<<"Will construct colour flow basis for "<<n_qp<<" quark pairs."<<endl;
            else cout<<"Will construct trace basis for "<<n_g<<" gluons and "<<n_qp<<" quark pairs."<<endl;
            for (int n(0);n<n_qp;n++) out_filename+="qqb";
            for (int n(0);n<n_g;n++) out_filename+="g";
            basis=construct_basis(n_qp, n_g, m_process);
            colour_calc(basis, m_process, NC_order, out_filename);
            break;
        }
        default: {
            run_error();
            break;
        }
    }
    
    return 0;
}

void colour_calc(vector<colour_term> basis, process m_process, const int NC_order, string out_filename) {
    clock_t t1, t2;
    
    t1=clock();
    // normalise basis
    basis=normalise_basis(basis,NC_order);
    unsigned int DIM(basis.size());
    
    // print leg indices and particles
    cout << "\nProcess:\n" << "leg\t" << "particle\t" << "in/out (1/0)" << endl;
    for (size_t i(1);i<=m_process.no_of_legs();i++) cout << m_process.leg(i).first << "\t" << m_process.leg(i).second << "\t\t" << m_process.is_in_leg(i) <<endl;
    
    // print normalised basis vectors
    cout<<"\nNormalised Basis Vectors:"<<endl;
    for (size_t i(0);i<DIM;i++)
        cout<<"b_"<<i+1<<" = "<<basis.at(i).build_string()<<endl;
    
    // print computation time for basis construction
    t1=clock()-t1;
    float runtime=(float)t1/CLOCKS_PER_SEC;
    cout << "\ncomputation time for basis construction: " << (int)runtime/3600 << " h, " << (int)runtime%60/60 << " m, " << (int)runtime%60+runtime-(int)runtime << " s"<< endl;
    
    t2=clock();
    
    // calculate and print soft matrix (to file)
    c_matrix soft_matrix=calc_soft_matrix(basis,NC_order);
    cout << "\nSoft Matrix:" << endl;
    soft_matrix.print();
    ofstream file;
    file.open(out_filename+"_met.dat");
    file<<"# Dimension of Square Matrix\n"<<DIM<<"\n"<<endl;
    file<<"# Soft Matrix / Colour Metric"<<endl;
    for (size_t i(0);i<DIM;i++)
        for (size_t j(0);j<DIM;j++)
            file<<fixed<<setprecision(17)<<soft_matrix[i][j].real()<<" ";
    file.close();
    
    // calculate and print inverse soft matrix
    c_matrix inv_soft_matrix=calc_inv_soft_matrix(soft_matrix);
    cout << "\nInverse Soft Matrix:" << endl;
    inv_soft_matrix.print();
    
    // calculate and print inverse soft matrix times soft matrix
    c_matrix unit_matrix(DIM);
    cout<<"\nInverse times Soft Matrix"<<endl;
    unit_matrix=soft_matrix*inv_soft_matrix;
    unit_matrix.print();
    
    // calculate and give out colour change matrices for all possible insertions
    bool multiply_with_inv_sm=false;
    vector<c_matrix> colour_change_matrices;
    file.open(out_filename+".dat");
    file<<"# Dimension of Square Matrices\n"<<DIM<<"\n"<<endl;
    cout<<"\nColour Change Matrices (";
    file<<"# Colour Change Matrices (";
    if (!multiply_with_inv_sm) {
        cout<<"not ";
        file<<"not ";
    }
    cout<<"multiplied with inverse colour metric)"<<endl;
    file<<"multiplied with inverse colour metric)"<<endl;
    for (unsigned int lno1(1);lno1<=m_process.no_of_legs();lno1++) {
        for (unsigned int lno2(lno1+1);lno2<=m_process.no_of_legs();lno2++) {
            colour_change_matrices.push_back(calc_colour_change_matrix(basis,soft_matrix,m_process,lno1,lno2,NC_order));
            
            // multiply with inverse soft matrix if wanted
            if (multiply_with_inv_sm) colour_change_matrices.back()=inv_soft_matrix*colour_change_matrices.back();
            
            cout << "C_(" << lno1 << "," << lno2 << ") = " << endl;
            colour_change_matrices.back().print();
            file<<"// C_("<<lno1<<","<<lno2<<") with PIDs "<<m_process.leg(lno1).second<<" and "<<m_process.leg(lno2).second<<endl;
            for (size_t i(0);i<DIM;i++) {
                for (size_t j(0);j<DIM;j++)
                    file<<fixed<<setprecision(17)<<colour_change_matrices.back()[i][j].real()<<" ";
                file<<endl;
            }
            cout<<endl;
        }
    }
    file.close();
    
    // print computation time for colour insertions
    t2=clock()-t2;
    runtime=(float)t2/CLOCKS_PER_SEC;
    cout << "computation time for colour insertions: " << runtime << " s"<< endl;
}

void run_error() {
    cerr<<"Error: competing input given. Please specify EITHER a colour term to be simplified (-s)/evaluated (-e) OR perform a colour space calculation by specifying a basis through a file name (-f) OR the process by the number of quark pairs (-nqp) and number of gluons (-ng). See also the help menu for more information (-h)."<<endl;
    exit(EXIT_FAILURE);
}
