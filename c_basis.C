#include "c_basis.h"
#include<iomanip>
#include<fstream>

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
void c_basis::print_to_file(string filename) {
    string out_filename="";
    if (filename=="") {
        size_t n_q(0), n_qb(0), n_g(0);
        for (size_t lno(1);lno<=m_process.no_of_legs();lno++) {
            string ptcl=m_process.leg(lno).second;
            
            if (ptcl=="q") n_q++;
            else if (ptcl=="qb") n_qb++;
            else if (ptcl=="g") n_g++;
        }
        if (n_q!=0) out_filename+=to_string(n_q)+"q";
        if (n_qb!=0) out_filename+=to_string(n_qb)+"qb";
        if (n_g!=0) out_filename+=to_string(n_g)+"g";
    }
    else out_filename=filename;
    
    ofstream file;
    file.open(out_filename+".dat");
    
    if (m_amp_perms.size()>0) {
        file<<"% number of basis vectors with non-zero coefficient"<<endl;
        file<<"SIZE_CONNECTED = "<<m_amp_perms.size()<<";\n"<<endl;
        file<<"% prefactors of the partial amplitude in total amplitude"<<endl;
        for (size_t i(0);i<m_amp_perms.size();i++) {
            // TODO: distinguish between adjoint and trace basis here!
            file<<"a_"<<i<<" = "<<fixed<<setprecision(17)<<m_normalisations.at(i)*m_confact<<";"<<endl;
        }
        
        file<<"\n% permutations defining the partial amplitudes"<<endl;
        for (size_t i(0);i<m_amp_perms.size();i++) {
            file<<"A_"<<i<<" =";
            for (const auto& ind : m_amp_perms[i]) file<<" "<<ind-1;
            file<<";"<<endl;
        }
        file<<"\n"<<endl;
    }
    
    file<<"% dimension of the colour space"<<endl;
    file<<"DIM = "<<m_dim<<";"<<endl;
    
    file<<"\n% colour metric aka soft matrix"<<endl;
    file<<"METRIC =";
    for (size_t m(0);m<m_dim;m++)
        for (size_t n(0);n<m_dim;n++) file<<fixed<<setprecision(17)<<" "<<m_smat[m][n].real();
    file<<";"<<endl;
    
    file<<"\n% colour change matrices"<<endl;
    size_t i(0);
    for (unsigned int lno1(1);lno1<=m_process.no_of_legs();lno1++) {
        for (unsigned int lno2(lno1+1);lno2<=m_process.no_of_legs();lno2++) {
            file<<"C_"<<lno1-1<<lno2-1<<" =";
            c_matrix ccm(m_ccmats.at(i));
            for (size_t m(0);m<m_dim;m++)
                for (size_t n(0);n<m_dim;n++)
                    file<<fixed<<setprecision(17)<<" "<<ccm[m][n].real();
            file<<";"<<endl;
            i++;
        }
    }
    
    file.close();
}
c_matrix c_basis::sm() {
    for (size_t i(0); i<m_dim; i++) {
        for (size_t j(0); j<m_dim; j++) {
            cout << "\rCalculating S" << "[" << i << "][" << j << "]..." << flush;
            complex<double> r;
            if (j>=i) r = m_ca_basis.at(i).scprod(m_ca_basis.at(j));
            else r = m_smat[j][i];
            m_smat[i][j] = (abs(r) < eps ? 0. : r);
        }
    }
    
    return m_smat;
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
    for (unsigned int lno1(1);lno1<=m_process.no_of_legs();lno1++) {
        for (unsigned int lno2(lno1+1);lno2<=m_process.no_of_legs();lno2++) {
            m_ccmats.push_back(this->ccm(lno1, lno2));
        }
    }
    return m_ccmats;
}
