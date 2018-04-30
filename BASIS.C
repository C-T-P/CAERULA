#include<string>
#include<iostream>
#include<vector>
#include<algorithm>
#include<complex>
#include<math.h>
#include "tensortools.h"
#include "Main.h"
#include "BASIS.h"
using namespace std;

void fully_simplify(terms& expr) {
    sort_indices(expr);
    sort_tensors(expr);
    simplify_terms(expr);
}
void sort_indices(terms &expr) {
    for (size_t i(0);i<expr.kron.size();i++) {
        for (size_t j(0);j<expr.kron[i].len();j++) if (expr.kron[i].index(j,0)>expr.kron[i].index(j,1)) expr.kron[i].swap_indices_at(j);
        for (size_t j(0);j<expr.sym[i].len();j++) {
            while (expr.sym[i].index(j,0)>expr.sym[i].index(j,1) || expr.sym[i].index(j,0)>expr.sym[i].index(j,2)) expr.sym[i].rotate_indices_at(j);
            if (expr.sym[i].index(j,1)>expr.sym[i].index(j,2)) expr.sym[i].swap_indices_at(j,1,2);
        }
        for (size_t j(0);j<expr.asym[i].len();j++) {
            while (expr.asym[i].index(j,0)>expr.asym[i].index(j,1) || expr.asym[i].index(j,0)>expr.asym[i].index(j,2)) expr.asym[i].rotate_indices_at(j);
            if (expr.asym[i].index(j,1)>expr.asym[i].index(j,2)) {
                expr.asym[i].swap_indices_at(j,1,2);
                expr.pref[i]*=-1.;
            }
        }
    }
}
void sort_tensors(terms &expr) {
    for (size_t i(0);i<expr.sym.size();i++) {
        expr.sym[i].sort_list();
        expr.asym[i].sort_list();
        expr.fund[i].sort_list();
        expr.kron[i].sort_list();
    }
}
void simplify_terms(terms &expr) {
    size_t i(0);
    float eps=1.e-6;
    while (i<expr.no_of_terms()) {
        size_t j(i+1);
        while (j<expr.no_of_terms()) {
            if (expr.sym[i].get_all_indices()==expr.sym[j].get_all_indices() && expr.asym[i].get_all_indices()==expr.asym[j].get_all_indices() && expr.fund[i].get_all_indices()==expr.fund[j].get_all_indices() && expr.kron[i].get_all_indices()==expr.kron[j].get_all_indices()) {
                expr.pref[i]+=expr.pref[j];
                delete_term(j, expr);
            }
            else j++;
        }
        i++;
    }
    for (size_t i(0);i<expr.pref.size(); i++) if (abs(expr.pref[i].real())<eps && abs(expr.pref[i].imag())<eps) delete_term(i,expr);
}
