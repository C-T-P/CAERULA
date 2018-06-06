#include "tensortools.h"
#include "Main.h"
#include "CONTRACT.h"
#include "I3NSERT.h"
#include "BASIS.h"

// construct trace basis for processes with n gluons and up to 1 quark and 1 antiquark
std::vector<colour_term> construct_trace_basis(diagram process) {
    std::vector<colour_term> trace_basis;
    three_ind symmetric;
    three_ind antisymmetric;
    three_ind fundamental;
    two_ind kronecker;
    std::complex<float> prefactor(1.);
    
    
    
    return trace_basis;
}

// calculate soft matrix
std::vector<std::vector<std::complex<float>>> calc_soft_matrix(std::vector<colour_term> basis) {
    std::vector<std::vector<std::complex<float>>> soft_matrix(basis.size(), std::vector<std::complex<float>>(basis.size(),0.));
    colour_term ct;
    ct.delete_all_terms();
    for (size_t i(0);i<basis.size();i++) {
        for (size_t j(0);j<basis.size();j++) {
            if(i<=j) {
                // transpose of basis[i] ?!
                cout << i << " " << j << endl;
                ct=basis[i].scprod(basis[j]);
                cout << ct.build_string() << endl;
                evaluate(ct);
                cout << " = " << ct.build_string() << "\n" << endl;
                soft_matrix[i][j]=ct.build_complex();
            }
            else soft_matrix[i][j]=conj(soft_matrix[j][i]);
        }
    }
    return soft_matrix;
} 

// calculate colour change matrix for insertion between leg lno1 and lno2
std::vector<std::vector<std::complex<float>>> calc_colour_change_matrix(std::vector<colour_term> basis, std::vector<std::vector<std::complex<float>>> soft_matrix, diagram process, unsigned int lno1, unsigned int lno2) {
    colour_term insertion_op=construct_insertion_op(process,lno1,lno2);
    colour_term ct;
    std::vector<std::vector<std::complex<float>>> ccm(basis.size(),std::vector<std::complex<float>>(basis.size(),0.));
    std::complex<double> t_ccm[basis.size()][basis.size()];
    for (size_t i(0);i<basis.size();i++) {
        for (size_t j(0);j<basis.size();j++) {
            ct=basis[i].scprod(insertion_op.multiply(make_internal(process,basis[j])));
            evaluate(ct);
            t_ccm[i][j]=ct.build_complex();
            if (isnan(t_ccm[i][j].real())) cerr << "Error: C_(" << lno1 << "," << lno2 << ")[" << i << "][" << j << "] = " <<ct.build_string() << "\n" << endl;
        }
    }
    
    // express in terms of basis
    gsl_vector_complex *b=gsl_vector_complex_alloc(basis.size()), *x=gsl_vector_complex_alloc(basis.size());
    gsl_matrix_complex *S=gsl_matrix_complex_alloc(basis.size(),basis.size());
    gsl_permutation *p=gsl_permutation_alloc(basis.size());
    int n;
    for (size_t j(0);j<basis.size();j++)
        for (size_t k(0);k<basis.size();k++) gsl_matrix_complex_set(S,j,k,gsl_complex_rect(soft_matrix[j][k].real(),soft_matrix[j][k].imag()));    
    for (size_t i(0);i<basis.size();i++) {
        for (size_t j(0);j<basis.size();j++) gsl_vector_complex_set(b,j,gsl_complex_rect(t_ccm[i][j].real(),t_ccm[i][j].imag()));
            
        gsl_linalg_complex_LU_decomp(S, p, &n);
        gsl_linalg_complex_LU_solve(S, p, b, x);
        for (size_t j(0);j<basis.size();j++) ccm[j][i]=std::complex<float>(GSL_REAL(gsl_vector_complex_get(x,j)),GSL_IMAG(gsl_vector_complex_get(x,j)));
    }
    
    gsl_permutation_free(p);
    gsl_vector_complex_free(x);
    gsl_vector_complex_free(b);
    gsl_matrix_complex_free(S);

    return ccm; 
}

// shift external to internal indices (by multiplying with Kroneckers)
colour_term make_internal(diagram process, colour_term expr) {
    bool gluonic(false);
    for (unsigned int lno(1);lno<=process.no_of_legs();lno++) {
        for (size_t tno(0);tno<expr.no_of_terms();tno++) {
            if (process.leg(lno).second==21) gluonic=true;
            else gluonic=false;
            expr.kron.at(tno).set_indices(lno+2000,lno,gluonic);
        }
    }
    evaluate(expr);
    return expr;
}
