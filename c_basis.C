#include "c_basis.h"


double eps = 1.e-4;

void c_basis::normalise() {
    for (size_t i(0);i<m_dim;i++) {
        complex<double> z=m_ca_basis.at(i).scprod(m_ca_basis.at(i));
        m_normalisations.at(i)=sqrt(abs(z));
        m_ca_basis.at(i)=m_ca_basis.at(i)*complex<double>(1./m_normalisations.at(i),0.);
    }
}
size_t c_basis::dim() {
    return m_ca_basis.size();
}
void c_basis::print() {
    for (size_t i(0); i<m_dim; i++) {
        cout<< "b_" << i+1 << " = ";
        m_ca_basis.at(i).print();
    }
}
c_matrix c_basis::sm() {
    c_matrix soft_mat(m_dim);
    
    for (size_t i(0); i<m_dim; i++) {
        for (size_t j(0); j<m_dim; j++) {
            cout << "\rCalculating S" << "[" << i << "][" << j << "]..." << flush;
            complex<double> r;
            if (j>=i) r = m_ca_basis.at(i).scprod(m_ca_basis.at(j));
            else r = soft_mat[j][i];
            soft_mat[i][j] = (abs(r) < eps ? 0. : r);
        }
    }
    
    return soft_mat;
}
c_matrix c_basis::ccm(size_t lno1, size_t lno2) {
    c_matrix colour_cm(m_dim);
    
    c_amplitude ins_op=construct_insertion_op(m_process,lno1,lno2);
    
    for (size_t i(0);i<m_dim;i++) {
        for (size_t j(0);j<m_dim;j++) {
            cout << "\rCalculating C_(" << lno1 << "," << lno2 << ")[" << i << "][" << j << "]..." << flush;
            c_amplitude bvj(m_ca_basis.at(j).shift_to_internal(2000));
            bvj.push_back(ins_op);
            complex<double> r(m_ca_basis.at(i).scprod(bvj));
            colour_cm[i][j] = (abs(r) < eps ? 0. : r);
        }
    }
    
    return colour_cm;
}

vector<c_matrix> c_basis::get_ccms() {
    vector<c_matrix> colour_change_mats;
    
    for (unsigned int lno1(1);lno1<=m_process.no_of_legs();lno1++) {
        for (unsigned int lno2(lno1+1);lno2<=m_process.no_of_legs();lno2++) {
            colour_change_mats.push_back(this->ccm(lno1, lno2));
        }
    }
    
    return colour_change_mats;
}
