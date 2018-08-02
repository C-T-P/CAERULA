#include "colourtools.h"
#include "Main.h"
#include "CONTRACT.h"
#include "FCONTRACT.h"

static double eps(1.e-4); // cutoff value - all doubles smaller than this will be treated as 0

// evaluates colour term to order NC_order in 1/NC
// if NC_order is set to INT_MAX, the expression is evaluated exactly to all orders in 1/NC
// expects expression with completely contracted indices and returns complex double - will return NAN if not all indices are contracted
complex<double> fast_evaluate_colour_term_to_order(colour_term& expr, unsigned int NC_order) {
    complex<double> cf(0);
    for (size_t t_it(0);t_it<expr.no_of_terms();t_it++)
        contract_product(expr.pref[t_it], expr.kron[t_it], expr.sym[t_it], expr.asym[t_it], expr.fund[t_it]);
//    clock_t t;
//    t=clock();
//    cout<<"start replacing"<<endl;
    replace_d_and_f_by_t(expr);
//    t=clock()-t;
//    cout<<"time to replace d and f: "<<(float)t/CLOCKS_PER_SEC<<endl;
//    t=clock();
    replace_fund_by_fierz(expr,NC_order);
//    t=clock()-t;
//    cout<<"time to replace fund by Fierz: "<<(float)t/CLOCKS_PER_SEC<<endl;
//    t=clock();
    size_t no_ot(expr.no_of_terms());
//    cout<<no_ot<<endl;
    for (size_t tno(0);tno<no_ot;tno++) {
        fast_replace_k(expr.pref[tno],expr.kron[tno]);
        cf+=expr.term(tno).build_complex();
    }
//    expr.delete_all_terms();
    if (abs(cf.real())<eps and abs(cf.imag())<eps) cf=0.;
//    t=clock()-t;
//    cout<<"time to evaluate Kroneckers: "<<(float)t/CLOCKS_PER_SEC<<endl;
    return cf;
}

// replace products of fundamental generators by Fierz rule
void replace_fund_by_fierz(colour_term& expr, int NC_order) {
//    cout<<"no of terms: "<<expr.no_of_terms()<<endl;
    for (size_t i(0);i<expr.no_of_terms();i++) {
        size_t it1(0);
        bool evaluated(false);
        while(!evaluated and expr.fund[i].len()>0) {
            pair<size_t,size_t> itf1(it1+1,0);
            if (expr.fund[i].index(it1,1)==expr.fund[i].index(it1,2)) {
                expr.delete_term(i);
                i--;
                if (i>=expr.no_of_terms()) evaluated=true;
            }
            else if (expr.fund[i].count_index(expr.fund[i].index(it1,0))>1) {
                itf1=expr.fund[i].find_index(expr.fund[i].index(it1,0),itf1.first);
                if (itf1.second==0) {
//                     clock_t t0(clock());
                    int a=expr.fund[i].index(it1,1), b=expr.fund[i].index(it1,2),
                    c=expr.fund[i].index(itf1.first,1), d=expr.fund[i].index(itf1.first,2);
                    expr.fund[i].del_indices(itf1.first);
                    expr.fund[i].del_indices(it1);
                    // second term
                    if (NC_order==INT_MAX or expr.NC_ctr[i]<NC_order) {
                        expr.add_term(expr.sym[i],expr.asym[i],expr.fund[i],expr.kron[i],expr.pref[i],expr.NC_ctr[i]);
                        expr.pref.back()*=-1./(2.*NC);
                        expr.NC_ctr.back()++;
                        expr.kron.back().set_indices(a,b,false);
                        expr.kron.back().set_indices(c,d,false);
                    }

                    // first term
                    expr.kron[i].set_indices(a,d,false);
                    expr.kron[i].set_indices(c,b,false);
                    expr.pref[i]*=0.5;
//                     t0=clock()-t0;
//                     if ((float)t0/CLOCKS_PER_SEC>0.1) cout<<(float)t0/CLOCKS_PER_SEC<<endl;
                }
                else it1++;
            }
            else it1++;
            if (it1>=expr.fund[i].len()-1) evaluated=true;
        }
    }
}

// contract indices of Kronecker deltas under the assumption that only kronecker deltas are left in the product and all indices are contracted
void fast_replace_k(complex<double>& prefactor, two_ind& kronecker) {
    while(kronecker.len()>0) {
        int ind_i=kronecker.index(0,0), ind_j=kronecker.index(0,1);
        if (ind_i==ind_j) {
            if (kronecker.is_gluonic(0)) prefactor*=(pow(NC,2)-1.);
            else prefactor*=NC;
        }
        else kronecker.find_and_rep_indices(ind_i,ind_j);
        kronecker.del_indices(0);
    }
}
