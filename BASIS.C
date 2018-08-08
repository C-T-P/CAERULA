#include<fstream>
#include<iomanip> 
#include "colourtools.h"
#include "c_matrix.h"
#include "Main.h"
#include "CONTRACT.h"
#include "FCONTRACT.h"
#include "I3NSERT.h"
#include "BASIS.h"

/*
===========================================================================
========================= Basis from external file ========================
===========================================================================
*/

vector<colour_term> read_basis(string filename, process& m_process) {
    vector<string> basis_strs;
    vector<colour_term> basis;
    
    read_in_process(filename, m_process, basis_strs);
    
    for (size_t i(0);i<basis_strs.size();i++)
        basis.push_back(decompose_terms(basis_strs.at(i),m_process));
    
    return basis;
}

// read in process file, store process information in process variable and basis vectors as strings
void read_in_process(string filename, process& m_process, vector<string>& basis_strs) {
    ifstream fin(filename);
    if (!fin) {
        cerr << "Error reading in process: file " << filename << " could not be opened." << endl;
        exit(EXIT_FAILURE);
    }
    else {
        string line;
        while (getline(fin,line)) {
            if (line.at(0)=='l') {
                string direc("");
                direc=line.substr(2,line.find("\t",2)-1);
                direc.erase(remove_if(direc.begin(),direc.end(),::isspace),direc.end());
                line.erase(0,line.find("\t",2));
                line.erase(remove_if(line.begin(),line.end(),::isspace),line.end());
                if (direc=="in") m_process.add_in_leg(line);
                else if (direc=="out") m_process.add_out_leg(line);
                else { 
                    cerr << "Error reading in process: leg needs either direction \"in\" or \"out\" , but was given direction \"" << direc << "\"." << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else if (line.at(0)=='b') {
                line.erase(remove_if(line.begin(),line.end(),::isspace),line.end());
                basis_strs.push_back(line.substr(1));
            }
        }
    }
}

// decompose input into colour term according to the indices given in the process specification
colour_term decompose_terms(string& input, process m_process) {
    // create mock process in case none was given
    if (m_process.no_of_legs()==0) {
        cout<<"\nNOTE: all Kronecker deltas k_[.,.] will be treated as having quark indices. Use c_[2.,0.]*t_[.,101,102]*t_[.,102,101] for gluonic deltas.\n"<<endl;
        for (int i(0);i<200;i++) m_process.add_out_leg("q");
    }
    
    colour_term expr;
    three_ind symmetric;
    three_ind antisymmetric;
    three_ind fundamental;
    two_ind kronecker;
    int NC_order_c;
    int NC_order_r;
    complex<double> prefactor = 1.;
    for (size_t i(0), mpos(input.find('+'));mpos!=string::npos || input.length()>0;mpos=input.find('+')) {
        ++i;
        // decompose input into summands
        string summand;
        if (mpos==string::npos) {
            summand=input;
            input="";
        }
        else {
            summand=input.substr(0,mpos);
            input=input.substr(mpos+1);
        }
        // decompose summand into factors
        kronecker.clear_indices();
        symmetric.clear_indices();
        antisymmetric.clear_indices();
        fundamental.clear_indices();
        prefactor=1.;
        NC_order_c=0;
        NC_order_r=0;
        for (size_t j(0), mpos(summand.find('*'));mpos!=string::npos || summand.length()>0;mpos=summand.find('*')) {
            ++j;
            string factor;
            if (mpos==string::npos) {
                factor=summand;
                summand="";
            }
            else {
                factor=summand.substr(0,mpos);
                summand=summand.substr(mpos+1);
            }
            if (factor.find("c_[")==0 && factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==string::npos || factor.find(',',cpos+1)!=string::npos) {
                    cerr << "Invalid prefactor." << endl;
                    exit(EXIT_FAILURE);
                }
                complex<double> prfct(0.);
                size_t spos(factor.substr(3,cpos-3).find("1/NC"));
                if (spos==string::npos) prfct+=stod(factor.substr(3,cpos-3));
                else {
                    size_t ppos=factor.substr(3,cpos-3).find('^');
                    double exponent(1);
                    if (ppos!=string::npos) exponent=stod(factor.substr(3+ppos+1,cpos-ppos-2));
                    NC_order_r+=exponent;
                    prfct+=1./pow(NC,exponent);
                }
                spos=factor.substr(cpos+1,factor.length()-cpos-2).find("1/NC");
                if (spos==string::npos) prfct+=stod(factor.substr(cpos+1,factor.length()-cpos-2))*complex<double>(0.,1.);
                else {
                    size_t ppos=factor.substr(cpos+1,factor.length()-cpos-2).find('^');
                    double exponent(1);
                    if (ppos!=string::npos) exponent=stod(factor.substr(cpos+ppos+2,factor.length()-ppos-3));
                    NC_order_c+=exponent;
                    prfct+=(double)1./pow(NC,exponent)*complex<double>(0.,1.);
                }
                // only same orders in 1/NC in real and imaginary part are supported
                if (NC_order_c!=NC_order_r and NC_order_c!=0 and NC_order_r!=0) {
                    cerr << "Error: Expected equal orders in 1/NC for real and imaginary part of prefactor, but real part has order " << NC_order_r << " and imaginary part has order " << NC_order_c << ".\n Specify real and imaginary part seperately to avoid this error." << endl;
                    exit(EXIT_FAILURE);
                }
                prefactor*=prfct;
            }
            else if(factor.find("f_[")==0 && factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==string::npos) {
                    cerr << "Invalid number of indices for f." << endl;
                    exit(EXIT_FAILURE);
                }
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==string::npos || factor.find(',',c2pos+1)!=string::npos) {
                    cerr << "Invalid number of indices for f." << endl;
                    exit(EXIT_FAILURE);
                }
                antisymmetric.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if(factor.find("d_[")==0 && factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==string::npos) {
                    cerr << "Invalid number of indices for d." << endl;
                    exit(EXIT_FAILURE);
                }
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==string::npos || factor.find(',',c2pos+1)!=string::npos) {
                    cerr << "Invalid number of indices for d." << endl;
                    exit(EXIT_FAILURE);
                }
                symmetric.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if(factor.find("t_[")==0 && factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==string::npos) {
                    cerr << "Invalid number of indices for t." << endl;
                    exit(EXIT_FAILURE);
                }
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==string::npos || factor.find(',',c2pos+1)!=string::npos) {
                    cerr << "Invalid number of indices for t." << endl;
                    exit(EXIT_FAILURE);
                }
                fundamental.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if (factor.find("k_[")==0 && factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==string::npos || factor.find(',',cpos+1)!=string::npos) {
                    cerr << "Invalid number of indices for k." << endl;
                    exit(EXIT_FAILURE);
                }
                int ind_i=stoi(factor.substr(3,cpos-3)), ind_j=stoi(factor.substr(cpos+1,factor.length()-cpos-2));
                bool gluonic(false);
                if (m_process.leg((unsigned int)ind_i).second=="g" or m_process.leg((unsigned int)ind_j).second=="g") gluonic=true;
                kronecker.set_indices(ind_i,ind_j,gluonic);
            }
            else cerr << "Invalid input." << endl;
        }
        expr.sym.push_back(symmetric);
        expr.asym.push_back(antisymmetric);
        expr.fund.push_back(fundamental);
        expr.kron.push_back(kronecker);
        expr.pref.push_back(prefactor);
        expr.NC_ctr.push_back(NC_order_r);
    }
    return expr;
}

// normalise Basis
vector<colour_term> normalise_basis(vector<colour_term> basis, int NC_order, vector<complex<double>>& normalisations) {
    unsigned int DIM(basis.size());
    vector<colour_term> n_basis;
    colour_term ct;
    for (size_t i(0);i<DIM;i++) {
        ct.delete_all_terms();
        ct=basis[i].scprod(basis[i]);
        normalisations.at(i)=sqrt(fast_evaluate_colour_term_to_order(ct,NC_order));
        ct=basis[i];
        for (size_t t_it(0);t_it<ct.no_of_terms();t_it++) ct.pref[t_it]*=1./normalisations.at(i);
        n_basis.push_back(ct);
    }
    return n_basis;
}

/*
===========================================================================
==================== Computation of Soft & Hard Matrix ====================
===========================================================================
 */

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
//                soft_matrix[i][j]=evaluate_colour_term_to_order(ct,NC_order);
                soft_matrix[i][j]=fast_evaluate_colour_term_to_order(ct,NC_order);
            }
            else soft_matrix[i][j]=conj(soft_matrix[j][i]);
        }
    }
    return soft_matrix;
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
//            ccm[i][j]=evaluate_colour_term_to_order(ct,NC_order);
            ccm[i][j]=fast_evaluate_colour_term_to_order(ct,NC_order);
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
