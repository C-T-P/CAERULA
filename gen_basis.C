#include "colourtools.h"
#include "trace_basis.h"
#include "f_basis.h"
#include "gen_basis.h"

// member functions of gen_basis class
gen_basis::gen_basis(string filename) {
    vector<string> basis_strs;
    
    // read in process file, store process information and basis vectors as strings
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
                if (direc=="in") m_in_legs.push_back(line);
                else if (direc=="out") m_out_legs.push_back(line);
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
    
    // decompose basis strings to colour objects
    for (auto& input : basis_strs) {
        colour_term expr;
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
            three_ind symmetric;
            three_ind antisymmetric;
            three_ind fundamental;
            two_ind kronecker;
            complex<double> prefactor = 1.;
            int NC_order(0);
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
                    int NC_order_c=0;
                    int NC_order_r=0;
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
                    NC_order=NC_order_r;
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
                    kronecker.set_indices(ind_i,ind_j,false);
                }
                else if (factor.find("K_[")==0 && factor[factor.length()-1]==']') {
                    size_t cpos(factor.find(','));
                    if (cpos==string::npos || factor.find(',',cpos+1)!=string::npos) {
                        cerr << "Invalid number of indices for K." << endl;
                        exit(EXIT_FAILURE);
                    }
                    int ind_i=stoi(factor.substr(3,cpos-3)), ind_j=stoi(factor.substr(cpos+1,factor.length()-cpos-2));
                    kronecker.set_indices(ind_i,ind_j,true);
                }
                else cerr << "Invalid input." << endl;
            }
            expr.sym.push_back(symmetric);
            expr.asym.push_back(antisymmetric);
            expr.fund.push_back(fundamental);
            expr.kron.push_back(kronecker);
            expr.pref.push_back(prefactor);
            expr.NC_ctr.push_back(NC_order);
        }
        m_gen_basis.push_back(expr);
    }
}
gen_basis::~gen_basis() {
    
}
size_t gen_basis::dim() {
    return m_gen_basis.size();
}
process gen_basis::proc() {
    process m_proc;
    
    for (const auto& l : m_in_legs) m_proc.add_in_leg(l);
    for (const auto& l : m_out_legs) m_proc.add_out_leg(l);
    return m_proc;
}
vector<colour_term> gen_basis::ct_basis() {
    return m_gen_basis;
}
void gen_basis::print() {
    for (auto& bv : m_gen_basis) cout<<bv.build_string()<<endl;
}
