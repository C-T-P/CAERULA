#ifndef C_MATRIX_H
#define C_MATRIX_H

#include<vector>
#include<complex>
#include<cmath>
using namespace std;

class c_matrix {
    vector<vector<complex<double>>> m_mat;
    public:
    c_matrix(unsigned int dim) {
        for (size_t i(0);i<dim;i++)
            m_mat.push_back(vector<complex<double>>(dim,0.));
    }
    size_t dim() {
        return m_mat.size();
    }
    complex<double> det() {
        int DIM(m_mat.size());
        gsl_matrix_complex *S=gsl_matrix_complex_alloc(DIM,DIM);
        for (size_t j(0);j<DIM;j++)
            for (size_t k(0);k<DIM;k++) gsl_matrix_complex_set(S,j,k,gsl_complex_rect(m_mat[j][k].real(),m_mat[j][k].imag()));
        int n;
        gsl_permutation *p=gsl_permutation_alloc(DIM);
        gsl_linalg_complex_LU_decomp(S, p, &n);
    }
    vector<complex<double>>& operator[](size_t i) {
        return m_mat[i];
    }
    vector<complex<double>> operator[](size_t i) const {
        return m_mat[i];
    }
    c_matrix operator*(c_matrix m_mat2) {
        size_t dim(m_mat2.dim());
        c_matrix m_matr(dim);
        if (dim==m_mat.size()) {
            for (size_t i(0);i<dim;i++)
                for (size_t j(0);j<dim;j++)
                    for (size_t k(0);k<dim;k++)
                        m_matr[i][j]+=m_mat[i][k]*m_mat2[k][j];
        }
        else {
            cout<<"Error mutliplying square matrices "<<&m_mat<<" and "<< &m_mat2<<": mismatch of dimensions: "<<m_mat.size()<<" vs. "<<dim<<endl;
            exit(EXIT_FAILURE);
        }
        return m_matr;
    }
    void print() {
        size_t dim(m_mat.size());
        for (size_t i(0);i<dim;i++) {
            for (size_t j(0);j<dim;j++)
                cout<<m_mat[i][j]<<"\t";
            cout<<endl;
        }
    }
};

#endif
