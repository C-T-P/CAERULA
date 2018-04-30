#include<string>
#include<iostream>
#include<vector>
#include<algorithm>
#include<complex>
#include<math.h>
#include "tensortools.h"
#include "Main.h"
#include "CONTRACT.h"
using namespace std;

three_ind symmetric;
three_ind antisymmetric;
three_ind fundamental;
two_ind kronecker;
two_ind kronecker_eval;
std::complex<float> prefactor = 1.;
static float NC = 3.; // number of colours

void evaluate(terms& expr) {
    replace_fund(expr);
    for (size_t i(0);i<expr.no_of_terms();i++) {
        reset();
        symmetric=expr.sym[i];
        antisymmetric=expr.asym[i];
        fundamental=expr.fund[i];
        kronecker=expr.kron[i];
        prefactor=expr.pref[i];
        
        // evaluate factor
        fully_contract();
    
        if (nonvanishing_expr()) {
            expr.sym[i]=symmetric;
            expr.asym[i]=antisymmetric;
            expr.fund[i]=fundamental;
            expr.kron[i]=kronecker_eval;
            expr.pref[i]=prefactor;
        }
        else delete_term(i, expr);
    }
}
void fully_contract() {
    replace_k();
    replace_f();
    replace_d();
    replace_2fd();
    replace_2df();
    
    // check for Tr(T_i)
    size_t i=0;
    while (i<fundamental.len()) {
        if (fundamental.index(i,1)==fundamental.index(i,2)) {
            del_all_indices();
            prefactor=0;
        }
        else i++;
    }
}
void replace_fund(terms& expr) {
    for (size_t i(0);i<expr.no_of_terms();i++) {
        size_t it1(0);
        bool evaluated(false);
        while(!evaluated && expr.fund[i].len()>0) {
            std::pair<size_t,size_t> itf1(it1+1,0);
            if (expr.fund[i].index(it1,1)==expr.fund[i].index(it1,2)) delete_term(i,expr);
            else if (expr.fund[i].count_index(expr.fund[i].index(it1,0))>1) {
                itf1=expr.fund[i].find_index(expr.fund[i].index(it1,0),itf1.first);
                if (itf1.second==0) {
                    int a=expr.fund[i].index(it1,1), b=expr.fund[i].index(it1,2), 
                    c=expr.fund[i].index(itf1.first,1), d=expr.fund[i].index(itf1.first,2);
                    expr.fund[i].del_indices(itf1.first);
                    expr.fund[i].del_indices(it1);
                    // second term
                    expr.sym.push_back(expr.sym[i]);
                    expr.asym.push_back(expr.asym[i]);
                    expr.fund.push_back(expr.fund[i]);
                    expr.kron.push_back(expr.kron[i]);
                    expr.pref.push_back(expr.pref[i]);
                    expr.pref[expr.no_of_terms()-1]*=-1./(2.*NC);
                    expr.kron[expr.no_of_terms()-1].set_indices(a,b);
                    expr.kron[expr.no_of_terms()-1].set_indices(c,d);
                    // first term
                    expr.kron[i].set_indices(a,d);
                    expr.kron[i].set_indices(c,b);
                    expr.pref[i]*=0.5;
                }
                else it1++;
            }
            else if (expr.fund[i].count_index(expr.fund[i].index(it1,2))>1) {
                itf1=expr.fund[i].find_index(expr.fund[i].index(it1,2),itf1.first);
                if (itf1.second==1) {
                    int a=expr.fund[i].index(it1,1), c=expr.fund[i].index(itf1.first,2), 
                    j=expr.fund[i].index(it1,0), k=expr.fund[i].index(itf1.first,0), x=expr.fund[i].index(it1,2);
                    expr.fund[i].del_indices(itf1.first);
                    expr.fund[i].del_indices(it1);
                    // third term
                    std::complex<float> z(0.,1./2.);
                    expr.sym.push_back(expr.sym[i]);
                    expr.asym.push_back(expr.asym[i]);
                    expr.fund.push_back(expr.fund[i]);
                    expr.kron.push_back(expr.kron[i]);
                    expr.pref.push_back(expr.pref[i]);
                    expr.pref[expr.no_of_terms()-1]*=z;
                    expr.asym[expr.no_of_terms()-1].set_indices(j,k,x);
                    expr.fund[expr.no_of_terms()-1].set_indices(x,a,c);
                    // second term
                    expr.sym.push_back(expr.sym[i]);
                    expr.asym.push_back(expr.asym[i]);
                    expr.fund.push_back(expr.fund[i]);
                    expr.kron.push_back(expr.kron[i]);
                    expr.pref.push_back(expr.pref[i]);
                    expr.pref[expr.no_of_terms()-1]*=1./2.;
                    expr.sym[expr.no_of_terms()-1].set_indices(j,k,x);
                    expr.fund[expr.no_of_terms()-1].set_indices(x,a,c);
                    // first term
                    expr.pref[i]*=1./(2.*NC);
                    expr.kron[i].set_indices(j,k);
                    expr.kron[i].set_indices(a,c);
                }
                else it1++;
            }
            else if (expr.fund[i].count_index(expr.fund[i].index(it1,1))>1) {
                itf1=expr.fund[i].find_index(expr.fund[i].index(it1,1),itf1.first);
                if (itf1.second==2) {
                    int a=expr.fund[i].index(itf1.first,1), c=expr.fund[i].index(it1,2), 
                    j=expr.fund[i].index(itf1.first,0), k=expr.fund[i].index(it1,0), x=expr.fund[i].index(it1,1);
                    expr.fund[i].del_indices(itf1.first);
                    expr.fund[i].del_indices(it1);
                    // third term
                    std::complex<float> z(0.,1./2.);
                    expr.sym.push_back(expr.sym[i]);
                    expr.asym.push_back(expr.asym[i]);
                    expr.fund.push_back(expr.fund[i]);
                    expr.kron.push_back(expr.kron[i]);
                    expr.pref.push_back(expr.pref[i]);
                    expr.pref[expr.no_of_terms()-1]*=z;
                    expr.asym[expr.no_of_terms()-1].set_indices(j,k,x);
                    expr.fund[expr.no_of_terms()-1].set_indices(x,a,c);
                    // second term
                    expr.sym.push_back(expr.sym[i]);
                    expr.asym.push_back(expr.asym[i]);
                    expr.fund.push_back(expr.fund[i]);
                    expr.kron.push_back(expr.kron[i]);
                    expr.pref.push_back(expr.pref[i]);
                    expr.pref[expr.no_of_terms()-1]*=1./2.;
                    expr.sym[expr.no_of_terms()-1].set_indices(j,k,x);
                    expr.fund[expr.no_of_terms()-1].set_indices(x,a,c);
                    // first term
                    expr.pref[i]*=1./(2.*NC);
                    expr.kron[i].set_indices(j,k);
                    expr.kron[i].set_indices(a,c);
                }
            }
             else it1++;
             if (it1>=expr.fund[i].len()) evaluated=true;
        }
    }
}
void replace_k() {
    while(kronecker.len()>0) {
        int ind_i=kronecker.index(0,0), ind_j=kronecker.index(0,1);
        if (!(is_free_index(ind_i))) {
            if (ind_i==ind_j) prefactor = prefactor*NC;
            else {
                symmetric.find_and_rep_indices(ind_i,ind_j);
                antisymmetric.find_and_rep_indices(ind_i,ind_j);
                fundamental.find_and_rep_indices(ind_i,ind_j);
                kronecker.find_and_rep_indices(ind_i,ind_j);
            }
        }
        else if (!(is_free_index(ind_j))) {
            symmetric.find_and_rep_indices(ind_j,ind_i);
            antisymmetric.find_and_rep_indices(ind_j,ind_i);
            fundamental.find_and_rep_indices(ind_j,ind_i);
            kronecker.find_and_rep_indices(ind_j,ind_i);
        }
        else kronecker_eval.set_indices(ind_i,ind_j);
        kronecker.del_indices(0);
    }
}
void replace_f() {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated && antisymmetric.len()>0) {
        if (antisymmetric.count_index_at(antisymmetric.index(it1,0),it1)>1 || antisymmetric.count_index_at(antisymmetric.index(it1,1),it1)>1) {
            prefactor=0;
            del_all_indices();
        }
        else if (antisymmetric.len()-it1>=2) {
            std::pair<size_t,size_t> itf1(it1+1,0), itf2(it1+1,0);
            bool part_evaluated(false);
            int n(0);
            while (!part_evaluated) {
                if (antisymmetric.count_index(antisymmetric.index(it1,1))>1 && antisymmetric.count_index(antisymmetric.index(it1,2))>1) {
                    itf1=antisymmetric.find_index(antisymmetric.index(it1,1),itf1.first);
                    itf2=antisymmetric.find_index(antisymmetric.index(it1,2),itf2.first);
                    // replace 2 f's by kronecker
                    if (antisymmetric.matching_indices(it1,itf1.first)>=2) {
                        while (antisymmetric.index(itf1.first,1)!=antisymmetric.index(it1,1)) antisymmetric.rotate_indices_at(itf1.first);
                        if (antisymmetric.index(it1,2)==antisymmetric.index(itf1.first,0)) {
                            prefactor*=-1;
                            antisymmetric.swap_indices_at(itf1.first,0,2);
                        }
                        else if (antisymmetric.index(it1,0)==antisymmetric.index(itf1.first,2)) {
                            prefactor*=-1;
                            antisymmetric.swap_indices_at(it1,0,2); 
                        }
                        else if (antisymmetric.index(it1,0)==antisymmetric.index(itf1.first,0)) {
                            antisymmetric.swap_indices_at(itf1.first,0,2);
                            antisymmetric.swap_indices_at(it1,0,2);
                        }
                        prefactor*=NC;
                        kronecker.set_indices(antisymmetric.index(it1,0),antisymmetric.index(itf1.first,0));
                        antisymmetric.del_indices(itf1.first);
                        antisymmetric.del_indices(it1);
                        fully_contract(); 
                        part_evaluated=true;
                    }
                    // replace 3 f's by one f
                    else if (antisymmetric.matching_indices(itf1.first,itf2.first)==1) {
                        if (itf1.second!=1) {
                            antisymmetric.swap_indices_at(itf1.first,itf1.second,1);
                            prefactor*=-1;
                        }
                        if (itf2.second!=2) {
                            antisymmetric.swap_indices_at(itf2.first,itf2.second,2);
                            prefactor*=-1;
                        }
                        while (antisymmetric.count_index_at(antisymmetric.index(itf1.first,2),itf2.first)==0) {
                            antisymmetric.swap_indices_at(itf1.first,0,2);
                            prefactor*=-1;
                        }
                        while (antisymmetric.index(itf1.first,2)!=antisymmetric.index(itf2.first,1)) {
                            antisymmetric.swap_indices_at(itf2.first,0,1);
                            prefactor*=-1;
                        }
                        antisymmetric.set_indices(antisymmetric.index(it1,0),antisymmetric.index(itf1.first,0),antisymmetric.index(itf2.first,0));
                        prefactor*=-NC/2.;
                        if (itf2.first>itf1.first) {
                            antisymmetric.del_indices(itf2.first);
                            antisymmetric.del_indices(itf1.first);
                        }
                        else {
                            antisymmetric.del_indices(itf1.first);
                            antisymmetric.del_indices(itf2.first);
                        }
                        antisymmetric.del_indices(it1);
                        fully_contract();
                        part_evaluated=true;
                    }
                    else {
                        it1++;
                        part_evaluated=true;
                    }
                }
                else if (n<3) {
                    antisymmetric.rotate_indices_at(it1);
                    n++;
                }
                else {
                    it1++;
                    part_evaluated=true;
                }
            }
        }
        else it1++;
        if (it1>=antisymmetric.len()) evaluated=true;
    }
}
void replace_d() {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated && symmetric.len()>0) {
        if (symmetric.count_index_at(symmetric.index(it1,0),it1)>1 || symmetric.count_index_at(symmetric.index(it1,1),it1)>1) {
            prefactor=0;
            del_all_indices();
        }
        else if (symmetric.len()-it1>=2) {
            std::pair<size_t,size_t> itf1(it1+1,0), itf2(it1+1,0);
            bool part_evaluated(false);
            int n(0);
            while (!part_evaluated) {
                if (symmetric.count_index(symmetric.index(it1,1))>1 && symmetric.count_index(symmetric.index(it1,2))>1) {
                    itf1=symmetric.find_index(symmetric.index(it1,1),itf1.first);
                    itf2=symmetric.find_index(symmetric.index(it1,2),itf2.first);                    
                    // replace 2 d's by kronecker
                    if (symmetric.matching_indices(it1,itf1.first)>=2) {
                        while (symmetric.index(itf1.first,1)!=symmetric.index(it1,1)) symmetric.rotate_indices_at(itf1.first);
                        if (symmetric.index(it1,2)==symmetric.index(itf1.first,0)) symmetric.swap_indices_at(itf1.first,0,2);
                        else if (symmetric.index(it1,0)==symmetric.index(itf1.first,2)) symmetric.swap_indices_at(it1,0,2); 
                        else if (symmetric.index(it1,0)==symmetric.index(itf1.first,0)) {
                            symmetric.swap_indices_at(itf1.first,0,2);
                            symmetric.swap_indices_at(it1,0,2);
                        }
                        prefactor*=(pow(NC,2)-4)/NC;
                        kronecker.set_indices(symmetric.index(it1,0),symmetric.index(itf1.first,0));
                        symmetric.del_indices(itf1.first);
                        symmetric.del_indices(it1);
                        fully_contract();
                        part_evaluated=true;
                    }
                    // replace 3 d's by one d
                    else if (symmetric.matching_indices(itf1.first,itf2.first)==1) {
                        if (itf1.second!=1) symmetric.swap_indices_at(itf1.first,itf1.second,1);
                        if (itf2.second!=2) symmetric.swap_indices_at(itf2.first,itf2.second,2);
                        while (symmetric.count_index_at(symmetric.index(itf1.first,2),itf2.first)==0) symmetric.swap_indices_at(itf1.first,0,2);
                        while (symmetric.index(itf1.first,2)!=symmetric.index(itf2.first,1)) symmetric.swap_indices_at(itf2.first,0,1);
                        symmetric.set_indices(symmetric.index(it1,0),symmetric.index(itf1.first,0),symmetric.index(itf2.first,0));
                        prefactor*=-(pow(NC,2)-12.0)/(2.*NC);
                        if (itf2.first>itf1.first) {
                            symmetric.del_indices(itf2.first);
                            symmetric.del_indices(itf1.first);
                        }
                        else {
                            symmetric.del_indices(itf1.first);
                            symmetric.del_indices(itf2.first);
                        }
                        symmetric.del_indices(it1);
                        fully_contract();
                        part_evaluated=true;
                    }
                    else { 
                        it1++;
                        part_evaluated=true;
                    }
                }
                else if (n<3) {
                    symmetric.rotate_indices_at(it1);
                    n++;
                }
                else {
                    it1++;
                    part_evaluated=true;
                }
            }
        }
        else it1++;
        if (it1>=symmetric.len()) evaluated=true;
    }
}
void replace_2fd() {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated && antisymmetric.len()>0 && symmetric.len()>0) {
        std::pair<size_t,size_t> itf0(0,0), itf1(0,0), itf2(0,0);
        bool part_evaluated(false);
        int n(0);
        while (!part_evaluated && antisymmetric.len()>0) {
            if (antisymmetric.count_index(symmetric.index(it1,1))>0 && antisymmetric.count_index(symmetric.index(it1,2))>0) {
                itf0=antisymmetric.find_index(symmetric.index(it1,0),itf0.first);
                itf1=antisymmetric.find_index(symmetric.index(it1,1),itf1.first);
                itf2=antisymmetric.find_index(symmetric.index(it1,2),itf2.first);
                // product of d and f with 2 contracted indices yields 0
                if (itf0.first==itf1.first || itf0.first==itf2.first || itf1.first==itf2.first) {
                    prefactor=0;
                    del_all_indices();
                }
                // replace 1 d and 2 f's by one d
                else if (antisymmetric.matching_indices(itf1.first,itf2.first)>=1) {
                    if (itf1.second!=1) {
                        antisymmetric.swap_indices_at(itf1.first,itf1.second,1);
                        prefactor*=-1;
                    }
                    if (itf2.second!=2) {
                        antisymmetric.swap_indices_at(itf2.first,itf2.second,2);
                        prefactor*=-1;
                    }
                    while (antisymmetric.count_index_at(antisymmetric.index(itf1.first,2),itf2.first)==0) {
                        antisymmetric.swap_indices_at(itf1.first,0,2);
                        prefactor*=-1;
                    }
                    while (antisymmetric.index(itf1.first,2)!=antisymmetric.index(itf2.first,1)) {
                        antisymmetric.swap_indices_at(itf2.first,0,1);
                        prefactor*=-1;
                    }
                    symmetric.set_indices(symmetric.index(it1,0),antisymmetric.index(itf1.first,0),antisymmetric.index(itf2.first,0));
                    prefactor*=-NC/2.;
                    if (itf2.first>itf1.first) {
                        antisymmetric.del_indices(itf2.first);
                        antisymmetric.del_indices(itf1.first);
                    }
                    else {
                        antisymmetric.del_indices(itf1.first);
                        antisymmetric.del_indices(itf2.first);
                    }
                    symmetric.del_indices(it1);
                    fully_contract();
                    part_evaluated=true;
                }
                else {
                    it1++;
                    part_evaluated=true;
                }
            }
            else if (n<3) {
                symmetric.rotate_indices_at(it1);
                n++;
            }
            else {
                it1++;
                part_evaluated=true;
            }
        }
        if (it1>=symmetric.len()) evaluated=true;
    }
}
void replace_2df() {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated && antisymmetric.len()>0 && symmetric.len()>0) {
        std::pair<size_t,size_t> itf0(0,0), itf1(0,0), itf2(0,0);
        bool part_evaluated(false);
        int n(0);
        while (!part_evaluated && symmetric.len()>0) {
            if (symmetric.count_index(antisymmetric.index(it1,1))>0 && symmetric.count_index(antisymmetric.index(it1,2))>0) {
                itf0=symmetric.find_index(antisymmetric.index(it1,0),itf0.first);
                itf1=symmetric.find_index(antisymmetric.index(it1,1),itf1.first);
                itf2=symmetric.find_index(antisymmetric.index(it1,2),itf2.first);
                // product of d and f with 2 contracted indices yields 0
                if (itf0.first==itf1.first || itf0.first == itf2.first || itf1.first==itf2.first) {
                    prefactor=0;
                    del_all_indices();
                }
                // replace 1 d and 2 f's by one d
                else if (symmetric.matching_indices(itf1.first,itf2.first)>=1) {
                    if (itf1.second!=1) symmetric.swap_indices_at(itf1.first,itf1.second,1);
                    if (itf2.second!=2) symmetric.swap_indices_at(itf2.first,itf2.second,2);
                    while (symmetric.count_index_at(symmetric.index(itf1.first,2),itf2.first)==0) symmetric.swap_indices_at(itf1.first,0,2);
                    while (symmetric.index(itf1.first,2)!=symmetric.index(itf2.first,1)) symmetric.swap_indices_at(itf2.first,0,1);
                    antisymmetric.set_indices(antisymmetric.index(it1,0),symmetric.index(itf1.first,0),symmetric.index(itf2.first,0));
                    prefactor*=(pow(NC,2)-4.)/(2.*NC);
                    if (itf2.first>itf1.first) {
                        symmetric.del_indices(itf2.first);
                        symmetric.del_indices(itf1.first);
                    }
                    else {
                        symmetric.del_indices(itf1.first);
                        symmetric.del_indices(itf2.first);
                    }
                    antisymmetric.del_indices(it1);
                    fully_contract(); 
                    part_evaluated=true;
                }
                else {
                    it1++;
                    part_evaluated=true;
                }
            }
            else if (n<3) {
                antisymmetric.rotate_indices_at(it1);
                n++;
            }
            else {
                it1++;
                part_evaluated=true;
            }
        }
        if (it1>=antisymmetric.len()) evaluated=true;
    }
}
bool is_free_index(int index) {
    int cntr(0);
    cntr+=kronecker.count_index(index);
    cntr+=symmetric.count_index(index);
    cntr+=antisymmetric.count_index(index);
    cntr+=fundamental.count_index(index);
    if (cntr<=1) return true;
    else return false;
}
void del_all_indices() {
    kronecker.clear_indices();
    kronecker_eval.clear_indices();
    symmetric.clear_indices();
    antisymmetric.clear_indices();
    fundamental.clear_indices();
}
void reset() {
    prefactor=1.;
    del_all_indices();
}
bool nonvanishing_expr() {
    if (prefactor.real()!=0. && prefactor.imag()!=0. && symmetric.len()!=0 && antisymmetric.len()!=0 && fundamental.len()!=0 && kronecker_eval.len()!=0) return false;
    else return true;
}
