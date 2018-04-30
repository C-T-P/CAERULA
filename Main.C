#include<string>
#include<iostream>
#include<vector>
#include<algorithm>
#include<complex>
#include<math.h>
#include "tensortools.h"
#include "BASIS.h"
#include "CONTRACT.h"
#include "Main.h"
using namespace std;

std::string build_string(terms& expr);
void decompose_terms(std::string& input, terms& expr);

int main(int argc, char **argv) {
    std::string input;
    terms expr;
    for (int i(1);i<argc;i++) input+=argv[i];
    //cout << factorise(input) << endl;
    decompose_terms(input, expr);
    evaluate(expr);
    fully_simplify(expr);
    cout << "= " << build_string(expr) << endl;
}
// std::string factorise(std::string string) {
//     if (std::count(string.begin(),string.end(), '(')==std::count(string.begin(),string.end(), ')')) {
//         for (size_t i(0), mpos(string.find('+'));mpos!=std::string::npos || string.length()>0;mpos=string.find('+')) {
//             ++i;
//             if (string.rfind(')',mpos)>=string.rfind('(',mpos) && string.find('(',mpos)<=string.rfind(')',mpos)) {
//                 std::string summand;
//                 if (mpos==std::string::npos) {
//                     summand=string;
//                     string="";
//                 }
//                 else {
//                     summand=string.substr(0,mpos);
//                     string=string.substr(mpos+1);
//                 }
//                 cout << summand << endl;
//             }
//         }
//     }
//     else {
//         cout << "Error: unmatched brackets in expression." << endl;
//         string="";
//     }
//     return string;
// }
void decompose_terms(std::string& input, terms& expr) {
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
            // TODO implement better error handling
            if (factor.find("c_[")==0 && factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==std::string::npos || factor.find(',',cpos+1)!=std::string::npos)
                    cout << "Invalid prefactor." << endl;
                std::complex<float> prfct=stof(factor.substr(3,cpos-3))+stof(factor.substr(cpos+1,factor.length()-cpos-2))*1.i;
                prefactor*=prfct;
            }
            else if(factor.find("f_[")==0 && factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==std::string::npos)
                    cout << "Invalid number of indices for f." << endl;
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==std::string::npos || factor.find(',',c2pos+1)!=std::string::npos)
                    cout << "Invalid number of indices for f." << endl;
                antisymmetric.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if(factor.find("d_[")==0 && factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==std::string::npos)
                    cout << "Invalid number of indices for d." << endl;
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==std::string::npos || factor.find(',',c2pos+1)!=std::string::npos)
                    cout << "Invalid number of indices for d." << endl;
                symmetric.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if(factor.find("t_[")==0 && factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==std::string::npos)
                    cout << "Invalid number of indices for t." << endl;
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==std::string::npos || factor.find(',',c2pos+1)!=std::string::npos)
                    cout << "Invalid number of indices for t." << endl;
                fundamental.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if (factor.find("k_[")==0 && factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==std::string::npos || factor.find(',',cpos+1)!=std::string::npos)
                    cout << "Invalid number of indices for k." << endl;
                kronecker.set_indices(stoi(factor.substr(3,cpos-3)),stoi(factor.substr(cpos+1,factor.length()-cpos-2)));
            }
            else cout << "Invalid input." << endl;
        }
        expr.sym.push_back(symmetric);
        expr.asym.push_back(antisymmetric);
        expr.fund.push_back(fundamental);
        expr.kron.push_back(kronecker);
        expr.pref.push_back(prefactor);
    }
}
void delete_term(int j, terms& expr) {
    expr.sym.erase(expr.sym.begin()+j);
    expr.asym.erase(expr.asym.begin()+j);
    expr.fund.erase(expr.fund.begin()+j);
    expr.kron.erase(expr.kron.begin()+j);
    expr.pref.erase(expr.pref.begin()+j);
}
std::string build_string(terms& expr) {
    std::string return_str;
    for (size_t it(0);it<expr.sym.size();it++) {
        std::string str="";
        if (expr.pref[it].real()!=0. || expr.pref[it].imag()!=0.) str+="c_["+to_string(expr.pref[it].real())+","+to_string(expr.pref[it].imag())+"]";
        if (expr.sym[it].len()>0) 
            for (size_t i(0); i<expr.sym[it].len(); i++) str+="*d_["+to_string(expr.sym[it].index(i,0))+","+to_string(expr.sym[it].index(i,1))+","+to_string(expr.sym[it].index(i,2))+"]";
        if (expr.asym[it].len()>0) 
            for (size_t i(0); i<expr.asym[it].len(); i++) str+="*f_["+to_string(expr.asym[it].index(i,0))+","+to_string(expr.asym[it].index(i,1))+","+to_string(expr.asym[it].index(i,2))+"]"; 
        if (expr.fund[it].len()>0)
            for (size_t i(0); i<expr.fund[it].len(); i++) str+="*t_["+to_string(expr.fund[it].index(i,0))+","+to_string(expr.fund[it].index(i,1))+","+to_string(expr.fund[it].index(i,2))+"]"; 
        if (expr.kron[it].len()>0)
            for (size_t i(0); i<expr.kron[it].len(); i++) str+="*k_["+to_string(expr.kron[it].index(i,0))+","+to_string(expr.kron[it].index(i,1))+"]";
        if (str!="" && return_str!="") return_str+="+";
        if (str!="") return_str+=str;
    }
    return return_str;
}
