#include<fstream>
#include<iomanip> 
#include "tensortools.h"
#include "c_matrix.h"
#include "Main.h"
#include "CONTRACT.h"
#include "I3NSERT.h"
#include "BASIS.h"

// construct basis
vector<colour_term> construct_basis(process m_process, vector<string> basis_strs) {
    vector<colour_term> basis;
    basis.clear();
    for (size_t i(0);i<basis_strs.size();i++)
        basis.push_back(decompose_terms(basis_strs.at(i),m_process));
    return basis;
}

// normalise Basis
vector<colour_term> normalise_basis(vector<colour_term> basis, int NC_order) {
    unsigned int DIM(basis.size());
    vector<colour_term> n_basis;
    colour_term ct;
    complex<double> norm;
    for (size_t i(0);i<DIM;i++) {
        ct.delete_all_terms();
        ct=basis[i].scprod(basis[i]);
        norm=evaluate_colour_term_to_order(ct,NC_order);
        ct=basis[i];
        for (size_t t_it(0);t_it<ct.no_of_terms();t_it++) ct.pref[t_it]*=1./sqrt(norm);
        n_basis.push_back(ct);
    }
    return n_basis;
}

// calculate soft matrix
c_matrix calc_soft_matrix(vector<colour_term> basis, int NC_order) {
    unsigned int DIM(basis.size());
    c_matrix soft_matrix(DIM);
    colour_term ct;
    for (size_t i(0);i<DIM;i++) {
        for (size_t j(0);j<DIM;j++) {
            if(i<=j) {
                ct.delete_all_terms();
                ct=basis[i].scprod(basis[j]);
                soft_matrix[i][j]=evaluate_colour_term_to_order(ct,NC_order);
            }
            else soft_matrix[i][j]=conj(soft_matrix[j][i]);
        }
    }
    return soft_matrix;
} 

// calculate inverse soft matrix
c_matrix calc_inv_soft_matrix(c_matrix soft_matrix) {
    unsigned int DIM(soft_matrix.dim());
    c_matrix inv_soft_matrix(DIM);
    
    gsl_matrix_complex *S=gsl_matrix_complex_alloc(DIM,DIM);
    gsl_matrix_complex *IS=gsl_matrix_complex_alloc(DIM,DIM);
    gsl_permutation *p=gsl_permutation_alloc(DIM);
    int n;
    for (size_t j(0);j<DIM;j++)
        for (size_t k(0);k<DIM;k++) gsl_matrix_complex_set(S,j,k,gsl_complex_rect(soft_matrix[j][k].real(),soft_matrix[j][k].imag())); 
    gsl_linalg_complex_LU_decomp(S,p,&n);
    
    gsl_complex gsl_det(gsl_linalg_complex_LU_det(S,n));
    complex<double> det(complex<double>(GSL_REAL(gsl_det),GSL_IMAG(gsl_det)));
    cout<<"\ndet(S) = "<<det<<endl;
    
    gsl_linalg_complex_LU_invert(S,p,IS);
    
    for (size_t i(0);i<DIM;i++)
        for (size_t j(0);j<DIM;j++) inv_soft_matrix[i][j]=complex<double>(GSL_REAL(gsl_matrix_complex_get(IS,i,j)),GSL_IMAG(gsl_matrix_complex_get(IS,i,j)));
    
    return inv_soft_matrix;
} 

// calculate colour change matrix for insertion between leg lno1 and lno2
c_matrix calc_colour_change_matrix(vector<colour_term> basis, c_matrix soft_matrix, process m_process, unsigned int lno1, unsigned int lno2, int NC_order) {
    colour_term insertion_op=construct_insertion_op(m_process,lno1,lno2);
    colour_term ct;
    unsigned int DIM(basis.size());
    c_matrix ccm(DIM);
    for (size_t i(0);i<DIM;i++) {
        for (size_t j(0);j<DIM;j++) {
            cout << "\rC_(" << lno1 << "," << lno2 << ")[" << i << "][" << j << "] being calculated..." << flush;
            ct=basis[i].scprod(insertion_op.multiply(make_internal(m_process,basis[j])));
            ccm[i][j]=evaluate_colour_term_to_order(ct,NC_order);
            if (isnan(ccm[i][j].real())) cerr << "Error: C_(" << lno1 << "," << lno2 << ")[" << i << "][" << j << "] = " <<ct.build_string() << "\n" << endl;
        }
    }
    cout << endl;
    return ccm;
}

// shift external to internal indices
colour_term make_internal(process m_process, colour_term expr) {
    for (unsigned int lno(1);lno<=m_process.no_of_legs();lno++) {
        for (size_t t_it(0);t_it<expr.no_of_terms();t_it++) {
            expr.sym[t_it].find_and_rep_indices(lno,lno+2000);
            expr.asym[t_it].find_and_rep_indices(lno,lno+2000);
            expr.fund[t_it].find_and_rep_indices(lno,lno+2000);
            expr.kron[t_it].find_and_rep_indices(lno,lno+2000);
        }
    }
    return expr;
}
