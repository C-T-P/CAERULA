#include<fstream>
#include "tensortools.h"
#include "BASIS.h"
#include "CONTRACT.h"
#include "I3NSERT.h"
#include "Main.h"

int main(int argc, char **argv) {
    // define order in 1/NC to which the terms shall be evaluated
    const int NC_order=INT_MAX;
    cout << "Order of 1/NC set to ";
    if (NC_order!=INT_MAX) cout << NC_order << "." << endl;
    else cout << "Infinity." << endl;
    
    // read in process file and basis vecs as strings and store process specs
    string filename;
    vector<string> basis_strs;
    basis_strs.clear();
    diagram process;
    for (int i(1);i<argc;i++) filename+=argv[i];
    read_in_process(filename, process, basis_strs);
    
    // give out leg indices and particle ids
    cout << "\nProcess:\n" << "leg\t" << "PID" << endl;
    for (size_t i(1);i<=process.no_of_legs();i++) cout << process.leg(i).first << "\t" << process.leg(i).second << endl;
    
    
    // save basis vectors as colour terms
    vector<colour_term> basis;
    basis.clear();
    cout << "\nBasis vectors:" << endl;
    for (size_t i(0);i<basis_strs.size();i++) {
        basis.push_back(decompose_terms(basis_strs.at(i),process));
        cout << "b_" << i+1 << " = " << basis[i].build_string() << endl;
    }
    
    // calculate and print out soft matrix
    vector<vector<complex<double>>> soft_matrix=calc_soft_matrix(basis,NC_order);
    cout << "\nSoft Matrix:" << endl;
    for (size_t i(0);i<soft_matrix.size();i++) {
        for (size_t j(0);j<soft_matrix[i].size();j++) cout << soft_matrix[i][j] << "\t";
        cout << endl;  
    }

    // calculate and print out inverse soft matrix
//     vector<vector<complex<double>>> inv_soft_matrix={{{268823., 268822., 268822., -194174., -194174., -194174.}, {268822., 
//   268823., 268822., -194174., -194174., -194174.}, {268822., 268822., 
//   268823., -194174., -194174., -194174.}, {-194174., -194174.,
// -194174., 140256., 140254., 140254.}, {-194174., -194174., -194174., 
//   140254., 140256., 140254.}, {-194174., -194174., -194174., 140254., 
//   140254., 140256.}}};
//     vector<vector<complex<double>>> inv_soft_matrix={{{851853., 851851., 851851., -615306., -615305., -615306.},{851851., 851853., 851851., -615306., -615306., -615305.},{851851., 851851., 851853., -615305., -615306., -615306.,},{-615306., -615306., -615305., 444446., 444444., 444444.},{-615305., -615306., -615306., 444444., 444446., 444444.,},{-615306., -615305., -615306., 444444., 444444., 444446.}}};
    vector<vector<complex<double>>> inv_soft_matrix=calc_inv_soft_matrix(soft_matrix);
    cout << "\nInverse Soft Matrix:" << endl;
    for (size_t i(0);i<inv_soft_matrix.size();i++) {
        for (size_t j(0);j<inv_soft_matrix[i].size();j++) cout << inv_soft_matrix[i][j] << "\t";
        cout << endl;
    }
    
    // inverse soft matrix times soft matrix
    vector<vector<complex<double>>> unit_matrix(soft_matrix.size(), vector<complex<double>>(soft_matrix.size(),0.));
    cout<<"\nInverse times Soft Matrix"<<endl;
    for (size_t i(0);i<soft_matrix.size();i++) {
        for (size_t j(0);j<soft_matrix.size();j++) {
            for (size_t k(0);k<soft_matrix.size();k++) {
                unit_matrix[i][j]+=soft_matrix[i][k]*inv_soft_matrix[k][j];
            }
            cout<<unit_matrix[i][j]<<"\t";
        }
        cout<<endl;
    }
    
    // calculate and give out colour change matrices for all possible insertions
    bool multiply_with_inv_sm=false;
    vector<vector<vector<complex<double>>>> colour_change_matrices;
    cout<<"\nColour Change Matrices (";
    if (!multiply_with_inv_sm) cout<<"not ";
    cout<<"multiplied with inverse colour metric)"<<endl;
    for (unsigned int lno1(1);lno1<=process.no_of_legs();lno1++) {
        for (unsigned int lno2(lno1+1);lno2<=process.no_of_legs();lno2++) {
            colour_change_matrices.push_back(calc_colour_change_matrix(basis,soft_matrix,process,lno1,lno2,NC_order,multiply_with_inv_sm));
            cout << "C_(" << lno1 << "," << lno2 << ") = " << endl;
            for (size_t i(0);i<colour_change_matrices.back().size();i++) {
                for (size_t j(0);j<colour_change_matrices.back()[i].size();j++) cout << colour_change_matrices.back()[i][j] << "\t";
                cout << endl;
            }
            cout << endl;
        }
    }
    
    // Debugging: check accuracy of SoftMatrix*ColourChangeMatrix=TProduct
//     vector<vector<vector<complex<double>>>> t_prods(colour_change_matrices.size(),vector<vector<complex<double>>>(colour_change_matrices[0].size(),vector<complex<double>>(colour_change_matrices[0][0].size(),0.)));
//     for (size_t t_it(0);t_it<colour_change_matrices.size();t_it++) {
//         for (size_t i(0);i<soft_matrix.size();i++) {
//             for (size_t j(0);j<soft_matrix.size();j++) {
//                 for (size_t k(0);k<soft_matrix.size();k++) {
//                     t_prods[t_it][i][j]+=soft_matrix[i][k]*colour_change_matrices[t_it][k][j];
//                 }
//                 cout<<t_prods[t_it][i][j]<<"\t";
//             }
//             cout<<endl;
//         }
//         cout<<endl;
//     }
    
    // calculate and print out sum of all colour change matrices
    vector<vector<complex<double>>> casimir(colour_change_matrices[0].size(),vector<complex<double>>(colour_change_matrices[0][0].size(),0.));
    for (size_t t_it(0);t_it<colour_change_matrices.size();t_it++) {
        for (size_t i(0);i<colour_change_matrices[t_it].size();i++) {
            for (size_t j(0);j<colour_change_matrices[t_it][i].size();j++) casimir[i][j]+=colour_change_matrices[t_it][i][j];
        }
    }
    cout<<"Casimir:"<<endl;
    for (size_t i(0);i<casimir.size();i++) {
        for (size_t j(0);j<casimir[i].size();j++) cout<<casimir[i][j]<<"\t";
        cout<<endl;
    }
    
    // save colour change matrices and soft matrix to files
    cout<<"\nPrinting information to files..."<<endl;
    save_colour_to_file(colour_change_matrices,soft_matrix,process);
    cout<<"Done!"<<endl;
    
    return 0;
} 

// read in the process file, store process information in process variable and basis vectors as strings
void read_in_process(string filename, diagram& process, vector<string>& basis_strs) {
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
                if (direc=="in") process.add_in_leg(stoi(line));
                else if (direc=="out") process.add_out_leg(stoi(line));
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

// decompose input into colour term according to the indices given in the process specification
colour_term decompose_terms(string& input, diagram process) {
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
                if (process.leg((unsigned int)ind_i).second==21 or process.leg((unsigned int)ind_j).second==21) gluonic=true;
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
