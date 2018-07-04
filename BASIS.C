#include<fstream>
#include<iomanip> 
#include "tensortools.h"
#include "Main.h"
#include "CONTRACT.h"
#include "I3NSERT.h"
#include "BASIS.h"


typedef numeric_limits<double> dbl;

// calculate soft matrix
vector<vector<complex<double>>> calc_soft_matrix(vector<colour_term> basis, int NC_order) {
    unsigned int DIM=basis.size();
    vector<vector<complex<double>>> soft_matrix(DIM, vector<complex<double>>(DIM,0.));
    colour_term ct;
    for (size_t i(0);i<DIM;i++) {
        for (size_t j(0);j<DIM;j++) {
            if(i<=j) {
//                 cout << "\n<b_" << i+1 << ",b_" << j+1 << "> = " << endl;
                ct.delete_all_terms();;
                ct=basis[i].scprod(basis[j]);
//                 cout << ct.build_string() << endl;
                soft_matrix[i][j]=evaluate_colour_term_to_order(ct,NC_order);
//                 cout << soft_matrix[i][j] << endl;
            }
            else soft_matrix[i][j]=conj(soft_matrix[j][i]);
        }
    }
    return soft_matrix;
} 

// calculate inverse soft matrix
vector<vector<complex<double>>> calc_inv_soft_matrix(vector<vector<complex<double>>> soft_matrix) {
    vector<vector<complex<double>>> inv_soft_matrix(soft_matrix.size(), vector<complex<double>>(soft_matrix.size(),0.));
    
    gsl_matrix_complex *S=gsl_matrix_complex_alloc(soft_matrix.size(),soft_matrix.size());
    gsl_matrix_complex *IS=gsl_matrix_complex_alloc(soft_matrix.size(),soft_matrix.size());
    gsl_permutation *p=gsl_permutation_alloc(soft_matrix.size());
    int n;
    for (size_t j(0);j<soft_matrix.size();j++)
        for (size_t k(0);k<soft_matrix.size();k++) gsl_matrix_complex_set(S,j,k,gsl_complex_rect(soft_matrix[j][k].real(),soft_matrix[j][k].imag())); 
    gsl_linalg_complex_LU_decomp(S,p,&n);
    gsl_linalg_complex_LU_invert(S,p,IS);
    
    for (size_t i(0);i<soft_matrix.size();i++)
        for (size_t j(0);j<soft_matrix.size();j++) inv_soft_matrix[i][j]=complex<double>(GSL_REAL(gsl_matrix_complex_get(IS,i,j)),GSL_IMAG(gsl_matrix_complex_get(IS,i,j)));
    
    return inv_soft_matrix;
} 

// calculate colour change matrix for insertion between leg lno1 and lno2
vector<vector<complex<double>>> calc_colour_change_matrix(vector<colour_term> basis, vector<vector<complex<double>>> soft_matrix, diagram process, unsigned int lno1, unsigned int lno2, int NC_order, bool inv_mult) {
    colour_term insertion_op=construct_insertion_op(process,lno1,lno2);
    colour_term ct;
    unsigned int DIM=basis.size();
    vector<vector<complex<double>>> ccm(DIM,vector<complex<double>>(DIM,0.));
    complex<double> t_ccm[DIM][DIM];
    for (size_t i(0);i<DIM;i++) {
        for (size_t j(0);j<DIM;j++) {
            cout << "\rC_(" << lno1 << "," << lno2 << ")[" << i << "][" << j << "] being calculated..." << flush;
            ct=basis[i].scprod(insertion_op.multiply(make_internal(process,basis[j])));
            t_ccm[i][j]=evaluate_colour_term_to_order(ct,NC_order);
            if (isnan(t_ccm[i][j].real())) cerr << "Error: C_(" << lno1 << "," << lno2 << ")[" << i << "][" << j << "] = " <<ct.build_string() << "\n" << endl;
        }
    }
    cout << endl;
    
    if (inv_mult) {
    // multiply with inverse metric
        gsl_vector_complex *b=gsl_vector_complex_alloc(DIM), *x=gsl_vector_complex_alloc(DIM);
        gsl_matrix_complex *S=gsl_matrix_complex_alloc(DIM,DIM);
        gsl_permutation *p=gsl_permutation_alloc(DIM);
        int n;
        for (size_t j(0);j<DIM;j++)
            for (size_t k(0);k<DIM;k++) gsl_matrix_complex_set(S,j,k,gsl_complex_rect(soft_matrix[j][k].real(),soft_matrix[j][k].imag()));    
        for (size_t i(0);i<DIM;i++) {
            for (size_t j(0);j<DIM;j++) gsl_vector_complex_set(b,j,gsl_complex_rect(t_ccm[i][j].real(),t_ccm[i][j].imag()));
                
            gsl_linalg_complex_LU_decomp(S, p, &n);
            gsl_linalg_complex_LU_solve(S, p, b, x);
            for (size_t j(0);j<DIM;j++) ccm[j][i]=complex<double>(GSL_REAL(gsl_vector_complex_get(x,j)),GSL_IMAG(gsl_vector_complex_get(x,j)));
        }
        gsl_permutation_free(p);
        gsl_vector_complex_free(x);
        gsl_vector_complex_free(b);
        gsl_matrix_complex_free(S);
    }
    else {
        for (size_t i(0);i<DIM;i++) 
            for (size_t j(0);j<DIM;j++) ccm[i][j]=t_ccm[i][j];
    }
    
    return ccm;
}

// shift external to internal indices
colour_term make_internal(diagram process, colour_term expr) {
    for (unsigned int lno(1);lno<=process.no_of_legs();lno++) {
        for (size_t t_it(0);t_it<expr.no_of_terms();t_it++) {
            expr.sym[t_it].find_and_rep_indices(lno,lno+2000);
            expr.asym[t_it].find_and_rep_indices(lno,lno+2000);
            expr.fund[t_it].find_and_rep_indices(lno,lno+2000);
            expr.kron[t_it].find_and_rep_indices(lno,lno+2000);
        }
    }
    return expr;
}

// print soft matrix
// void save_colour_to_file(vector<vector<vector<complex<double>>>> colour_change_matrices, vector<vector<complex<double>>> soft_matrix, diagram process) {
//     // TODO unambiguous filename definition
//     string filename;
//     for (size_t lno(1);lno<=process.no_of_legs();lno++) {
//         if (process.leg(lno).second==21) filename+='g';
//         else if ((process.leg(lno).second>=1 and process.leg(lno).second<=6 and process.is_in_leg(lno)) or (process.leg(lno).second>=-6 and process.leg(lno).second<=-1 and !process.is_in_leg(lno))) filename+='q';
//         else if ((process.leg(lno).second>=-6 and process.leg(lno).second<=-1 and process.is_in_leg(lno)) or (process.leg(lno).second>=1 and process.leg(lno).second<=6 and !process.is_in_leg(lno))) filename+="qb";
//         else {
//             cerr<<"Error saving Soft Matrix: Leg number "<<lno<<" should either be a gluon or (anti-)quark, but has PID "<<process.leg(lno).second<<endl;
//             exit(EXIT_FAILURE);
//         }
//     }
//     // TODO print only real parts?
//     
//     file.open(filename+"_met.dat");
//     file<<soft_matrix.size()<<"\n\n";
//     for (size_t i(0);i<soft_matrix.size();i++)
//         for (size_t j(0);j<soft_matrix[i].size();j++)
//             file<<fixed<<setprecision(17)<<soft_matrix[i][j].real()<<" ";
//     file.close();
// }





