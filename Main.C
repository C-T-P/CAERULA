#include<fstream>
#include "tensortools.h"
#include "BASIS.h"
#include "CONTRACT.h"
#include "I3NSERT.h"
#include "Main.h"

int main(int argc, char **argv) {
    // read in process file and basis vecs as strings and store process specs
    std::string filename;
    std::vector<std::string> basis_strs;
    basis_strs.clear();
    diagram process;
    for (int i(1);i<argc;i++) filename+=argv[i];
    read_in_process(filename, process, basis_strs);
    
    // debugging
//     colour_term ct;
//     string str="c_[0.000000,-1.000000]*d_[103,104,5]*f_[101,102,5]*t_[101,1,2]*t_[102,4,3]*t_[103,2,1]*t_[104,3,4]";
//     ct=decompose_terms(str,process);
//     cout << ct.build_string() << endl;
//     evaluate(ct);
//     cout << ct.build_string() << endl;
    
    // give out leg indices and particle ids (particle specs)
    cout << "\nProcess:\n" << "leg\t" << "PID" << endl;
    for (size_t i(1);i<=process.no_of_legs();i++) cout << process.leg(i).first << "\t" << process.leg(i).second << endl;
    
    // save basis vectors as colour terms
    std::vector<colour_term> basis;
    basis.clear();
    cout << "\nBasis vectors:" << endl;
    for (size_t i(0);i<basis_strs.size();i++) {
        basis.push_back(decompose_terms(basis_strs.at(i),process));
        cout << "b_" << i+1 << " = " << basis[i].build_string() << endl;
    }
    
    // calculate and give out soft matrix
    std::vector<std::vector<std::complex<float>>> soft_matrix=calc_soft_matrix(basis);
    cout << "\nSoft Matrix:" << endl;
    for (size_t i(0);i<soft_matrix.size();i++) {
        for (size_t j(0);j<soft_matrix[i].size();j++) 
            cout << soft_matrix[i][j] << "\t";
        cout << endl;
    }
    
    // calculate and give out colour change matrices for all possible insertions
    std::vector<std::vector<std::vector<std::complex<float>>>> colour_change_matrices;
    cout << "\nColour Change Matrices:" << endl;
    for (unsigned int lno1(1);lno1<=process.no_of_legs();lno1++) {
        for (unsigned int lno2(lno1+1);lno2<=process.no_of_legs();lno2++) {
            colour_change_matrices.push_back(calc_colour_change_matrix(basis,soft_matrix,process,lno1,lno2));
            cout << "C_(" << lno1 << "," << lno2 << ") = " << endl;
            for (size_t i(0);i<colour_change_matrices.back().size();i++) {
                for (size_t j(0);j<colour_change_matrices.back()[i].size();j++) 
                    cout << colour_change_matrices.back()[i][j] << "\t";
                cout << endl;
            }
            cout << endl;
        }
    }
    
    return 0;
} 

// read in the process file, store process information in process variable and basis vectors as strings
void read_in_process(std::string filename, diagram& process, std::vector<std::string>& basis_strs) {
    std::ifstream fin(filename);
    if (!fin) {
        cerr << "Error reading in process: file " << filename << " could not be opened." << endl;
        exit(EXIT_FAILURE);
    }
    else {
        std::string line;
        while (getline(fin,line)) {
            if (line.at(0)=='l') {
                std::string direc("");
                direc=line.substr(2,line.find("\t",2)-1);
                direc.erase(remove_if(direc.begin(),direc.end(),::isspace),direc.end());
                line.erase(0,line.find("\t",2));
                line.erase(remove_if(line.begin(),line.end(),::isspace),line.end());
                if (direc=="in") process.add_in_leg(std::stoi(line));
                else if (direc=="out") process.add_out_leg(std::stoi(line));
                else { 
                    cerr << "Error saving process specifications: leg needs either direction \"in\" or \"out\" , but was given direction \"" << direc << "\"." << endl;
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

//
colour_term decompose_terms(std::string& input, diagram process) {
    colour_term expr;
    three_ind symmetric;
    three_ind antisymmetric;
    three_ind fundamental;
    two_ind kronecker;
    std::complex<float> prefactor = 1.;
    for (size_t i(0), mpos(input.find('+'));mpos!=std::string::npos || input.length()>0;mpos=input.find('+')) {
        ++i;
        // decompose input into summands
        std::string summand;
        if (mpos==std::string::npos) {
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
        prefactor = 1.;
        for (size_t j(0), mpos(summand.find('*'));mpos!=std::string::npos || summand.length()>0;mpos=summand.find('*')) {
            ++j;
            std::string factor;
            if (mpos==std::string::npos) {
                factor=summand;
                summand="";
            }
            else {
                factor=summand.substr(0,mpos);
                summand=summand.substr(mpos+1);
            }
            if (factor.find("c_[")==0 && factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==std::string::npos || factor.find(',',cpos+1)!=std::string::npos) {
                    cerr << "Invalid prefactor." << endl;
                    exit(EXIT_FAILURE);
                }
                std::complex<float> prfct=stof(factor.substr(3,cpos-3))+stof(factor.substr(cpos+1,factor.length()-cpos-2))*complex<float>(0.,1.);
                prefactor*=prfct;
            }
            else if(factor.find("f_[")==0 && factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==std::string::npos) {
                    cerr << "Invalid number of indices for f." << endl;
                    exit(EXIT_FAILURE);
                }
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==std::string::npos || factor.find(',',c2pos+1)!=std::string::npos) {
                    cerr << "Invalid number of indices for f." << endl;
                    exit(EXIT_FAILURE);
                }
                antisymmetric.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if(factor.find("d_[")==0 && factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==std::string::npos) {
                    cerr << "Invalid number of indices for d." << endl;
                    exit(EXIT_FAILURE);
                }
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==std::string::npos || factor.find(',',c2pos+1)!=std::string::npos) {
                    cerr << "Invalid number of indices for d." << endl;
                    exit(EXIT_FAILURE);
                }
                symmetric.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if(factor.find("t_[")==0 && factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==std::string::npos) {
                    cerr << "Invalid number of indices for t." << endl;
                    exit(EXIT_FAILURE);
                }
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==std::string::npos || factor.find(',',c2pos+1)!=std::string::npos) {
                    cerr << "Invalid number of indices for t." << endl;
                    exit(EXIT_FAILURE);
                }
                fundamental.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if (factor.find("k_[")==0 && factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==std::string::npos || factor.find(',',cpos+1)!=std::string::npos) {
                    cerr << "Invalid number of indices for k." << endl;
                    exit(EXIT_FAILURE);
                }
                int ind_i=stoi(factor.substr(3,cpos-3)), ind_j=stoi(factor.substr(cpos+1,factor.length()-cpos-2));
                bool gluonic(false);
                if (process.leg((unsigned int)ind_i).second==21 and process.leg((unsigned int)ind_j).second==21) gluonic=true;
                kronecker.set_indices(ind_i,ind_j,gluonic);
            }
            else cerr << "Invalid input." << endl;
        }
        expr.sym.push_back(symmetric);
        expr.asym.push_back(antisymmetric);
        expr.fund.push_back(fundamental);
        expr.kron.push_back(kronecker);
        expr.pref.push_back(prefactor);
    }
    return expr;
}
