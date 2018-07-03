#include "tensortools.h"
#include "Main.h"
#include "CONTRACT.h"

static double eps(1.e-4); // cutoff value - all doubles smaller than this will be treated as 0

// evaluates colour term to order NC_order in 1/NC
// if NC_order is set to INT_MAX, the expression is evaluated exactly to all orders in 1/NC
// expects expression with completely contracted indices and returns complex double - will return NAN if not all indices are contracted
complex<double> evaluate_colour_term_to_order(colour_term& expr, unsigned int NC_order) {
    complex<double> cf(0);
    if (NC_order==INT_MAX) {
        for (size_t t_it(0);t_it<expr.no_of_terms();t_it++) {
            contract_product(expr.pref[t_it], expr.kron[t_it], expr.sym[t_it], expr.asym[t_it], expr.fund[t_it]);
            if (vanishing_expr(expr.pref[t_it], expr.kron[t_it], expr.sym[t_it], expr.asym[t_it], expr.fund[t_it])) {
                expr.delete_term(t_it);
                t_it--;
            }
        }
        replace_d_and_f_by_t(expr);
        sort_colour_term(expr);
        replace_fund(expr,INT_MAX);
    }
    else {
        for (size_t t_it(0);t_it<expr.no_of_terms();t_it++) replace_k(expr.pref[t_it], expr.kron[t_it], expr.sym[t_it], expr.asym[t_it], expr.fund[t_it]);
        replace_d_and_f_by_t(expr);
        replace_fund(expr,NC_order);
    }
    while (expr.no_of_terms()>0) {
        contract_product(expr.pref[0],expr.kron[0],expr.sym[0],expr.asym[0],expr.fund[0]);
        cf+=expr.term(0).build_complex();
        expr.delete_term(0);
    }
    if (abs(cf.real())<eps and abs(cf.imag())<eps) cf=0.;
    return cf;
}

// simplifies colour term by contracting repeated indices - yields colour term with smaller amount of internal indices
// accepts colour terms with uncontracted external indices
void simplify_colour_term(colour_term& expr) {
    for (size_t i(0);i<expr.no_of_terms();i++) {
        contract_product(expr.pref[i], expr.kron[i], expr.sym[i], expr.asym[i], expr.fund[i]);
        if (vanishing_expr(expr.pref[i], expr.kron[i], expr.sym[i], expr.asym[i], expr.fund[i])) {
            expr.delete_term(i);
            i--;
        }
    }
    sort_colour_term(expr);
    if(replace_fund(expr,INT_MAX)) simplify_colour_term(expr);
}

// sort indices and terms and add terms of identical tensors in an expression
void sort_colour_term(colour_term& expr) {
    sort_indices(expr);
    sort_tensors(expr);
    add_terms(expr);
}

// sort indices with increasing value
void sort_indices(colour_term& expr) {
    for (size_t t_it(0);t_it<expr.no_of_terms();t_it++) {
        for (size_t p_it(0);p_it<expr.kron[t_it].len();t_it++) expr.kron[t_it].sort_indices_at(p_it);
        for (size_t p_it(0);p_it<expr.sym[t_it].len();t_it++) expr.sym[t_it].sort_indices_at(p_it);
        for (size_t p_it(0);p_it<expr.asym[t_it].len();t_it++) expr.asym[t_it].sort_indices_at(p_it);
    }
}

// sort tensors in addends with increasing value of their first indices
void sort_tensors(colour_term& expr) {
    for (size_t i(0);i<expr.sym.size();i++) {
        expr.sym[i].sort_list();
        expr.asym[i].sort_list();
        expr.fund[i].sort_list();
        expr.kron[i].sort_list();
    }
}

// add terms of identical tensors
void add_terms(colour_term& expr) {
    size_t i(0);
    while (i<expr.no_of_terms()) {
        size_t j(i+1);
        while (j<expr.no_of_terms()) {
            if (expr.sym[i].get_all_indices()==expr.sym[j].get_all_indices() and expr.asym[i].get_all_indices()==expr.asym[j].get_all_indices() and expr.fund[i].get_all_indices()==expr.fund[j].get_all_indices() and expr.kron[i].get_all_indices()==expr.kron[j].get_all_indices()) {
                expr.pref[i]+=expr.pref[j];
                expr.delete_term(j);
            }
            else j++;
        }
        i++;
    }
    for (size_t i(0);i<expr.pref.size();i++) if (abs(expr.pref[i].real())<eps and abs(expr.pref[i].imag())<eps) expr.delete_term(i);
}


// contract all repeated indices in one product according to identities of delta, f, d, and t
void contract_product(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    replace_k(prefactor, kronecker, symmetric, antisymmetric, fundamental);
    check_fund_trace(prefactor, kronecker, symmetric, antisymmetric, fundamental);
    replace_f(prefactor, kronecker, symmetric, antisymmetric, fundamental);
    replace_d(prefactor, kronecker, symmetric, antisymmetric, fundamental);
    replace_2fd(prefactor, kronecker, symmetric, antisymmetric, fundamental);
    replace_2df(prefactor, kronecker, symmetric, antisymmetric, fundamental);
    check_fund_trace(prefactor, kronecker, symmetric, antisymmetric, fundamental);
}

// replace products of symmetric and antisymmetric structure constants by fundamental generators
void replace_d_and_f_by_t(colour_term& expr){
    // replace symmetric structure constants
    for (size_t t_it(0);t_it<expr.no_of_terms();t_it++) {
        while (expr.sym[t_it].len()>0) {
            int ind_i=expr.sym[t_it].index(0,0), ind_j=expr.sym[t_it].index(0,1), ind_k=expr.sym[t_it].index(0,2);
            int int_ind_a=find_free_internal_ind(101, expr.kron[t_it], expr.sym[t_it], expr.asym[t_it], expr.fund[t_it]), int_ind_b=find_free_internal_ind(int_ind_a+1, expr.kron[t_it], expr.sym[t_it], expr.asym[t_it], expr.fund[t_it]), int_ind_c=find_free_internal_ind(int_ind_b+1, expr.kron[t_it], expr.sym[t_it], expr.asym[t_it], expr.fund[t_it]);
            expr.sym[t_it].del_indices(0);
    
            // second term
            expr.sym.push_back(expr.sym[t_it]);
            expr.asym.push_back(expr.asym[t_it]);
            expr.fund.push_back(expr.fund[t_it]);
            expr.kron.push_back(expr.kron[t_it]);
            expr.pref.push_back(expr.pref[t_it]);
            expr.NC_ctr.push_back(expr.NC_ctr[t_it]);
            expr.pref[expr.no_of_terms()-1]*=complex<double>(2.,0.);
            expr.fund[expr.no_of_terms()-1].set_indices(ind_j, int_ind_a, int_ind_b);
            expr.fund[expr.no_of_terms()-1].set_indices(ind_i, int_ind_b, int_ind_c);
            expr.fund[expr.no_of_terms()-1].set_indices(ind_k, int_ind_c, int_ind_a);
            
            // first term
            expr.pref[t_it]*=complex<double>(2.,0.);
            expr.fund[t_it].set_indices(ind_i, int_ind_a, int_ind_b);
            expr.fund[t_it].set_indices(ind_j, int_ind_b, int_ind_c);
            expr.fund[t_it].set_indices(ind_k, int_ind_c, int_ind_a);
        }
    }
    
    // replace antisymmetric structure constants
    for (size_t t_it(0);t_it<expr.no_of_terms();t_it++) {
        while (expr.asym[t_it].len()>0) {
            int ind_i=expr.asym[t_it].index(0,0), ind_j=expr.asym[t_it].index(0,1), ind_k=expr.asym[t_it].index(0,2);
            int int_ind_a=find_free_internal_ind(101, expr.kron[t_it], expr.sym[t_it], expr.asym[t_it], expr.fund[t_it]), int_ind_b=find_free_internal_ind(int_ind_a+1, expr.kron[t_it], expr.sym[t_it], expr.asym[t_it], expr.fund[t_it]), int_ind_c=find_free_internal_ind(int_ind_b+1, expr.kron[t_it], expr.sym[t_it], expr.asym[t_it], expr.fund[t_it]);
            expr.asym[t_it].del_indices(0);
    
            // second term
            expr.sym.push_back(expr.sym[t_it]);
            expr.asym.push_back(expr.asym[t_it]);
            expr.fund.push_back(expr.fund[t_it]);
            expr.kron.push_back(expr.kron[t_it]);
            expr.pref.push_back(expr.pref[t_it]);
            expr.NC_ctr.push_back(expr.NC_ctr[t_it]);
            expr.pref[expr.no_of_terms()-1]*=complex<double>(0.,2.);
            expr.fund[expr.no_of_terms()-1].set_indices(ind_j, int_ind_a, int_ind_b);
            expr.fund[expr.no_of_terms()-1].set_indices(ind_i, int_ind_b, int_ind_c);
            expr.fund[expr.no_of_terms()-1].set_indices(ind_k, int_ind_c, int_ind_a);
            
            // first term
            expr.pref[t_it]*=complex<double>(0.,-2.);
            expr.fund[t_it].set_indices(ind_i, int_ind_a, int_ind_b);
            expr.fund[t_it].set_indices(ind_j, int_ind_b, int_ind_c);
            expr.fund[t_it].set_indices(ind_k, int_ind_c, int_ind_a);
        }
    }
}

// replace products of fundamental generators
bool replace_fund(colour_term& expr, int NC_order) {
    bool replacements_made(false);
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
                    int a=expr.fund[i].index(it1,1), b=expr.fund[i].index(it1,2), 
                    c=expr.fund[i].index(itf1.first,1), d=expr.fund[i].index(itf1.first,2);
                    expr.fund[i].del_indices(itf1.first);
                    expr.fund[i].del_indices(it1);
                    // second term
                    if (NC_order==INT_MAX or expr.NC_ctr[i]<NC_order) {
                        expr.sym.push_back(expr.sym[i]);
                        expr.asym.push_back(expr.asym[i]);
                        expr.fund.push_back(expr.fund[i]);
                        expr.kron.push_back(expr.kron[i]);
                        expr.pref.push_back(expr.pref[i]);
                        expr.NC_ctr.push_back(expr.NC_ctr[i]);
                        expr.pref[expr.no_of_terms()-1]*=-1./(2.*NC);
                        expr.NC_ctr[expr.no_of_terms()-1]++;
                        expr.kron[expr.no_of_terms()-1].set_indices(a,b,false);
                        expr.kron[expr.no_of_terms()-1].set_indices(c,d,false);
                    }
                    
                    // first term
                    expr.kron[i].set_indices(a,d,false);
                    expr.kron[i].set_indices(c,b,false);
                    expr.pref[i]*=0.5;
                    replacements_made=true;
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
                    expr.sym.push_back(expr.sym[i]);
                    expr.asym.push_back(expr.asym[i]);
                    expr.fund.push_back(expr.fund[i]);
                    expr.kron.push_back(expr.kron[i]);
                    expr.pref.push_back(expr.pref[i]);
                    expr.NC_ctr.push_back(expr.NC_ctr[i]);
                    expr.pref[expr.no_of_terms()-1]*=complex<double>(0.,1./2.);
                    expr.asym[expr.no_of_terms()-1].set_indices(j,k,x);
                    expr.fund[expr.no_of_terms()-1].set_indices(x,a,c);
                    // second term
                    expr.sym.push_back(expr.sym[i]);
                    expr.asym.push_back(expr.asym[i]);
                    expr.fund.push_back(expr.fund[i]);
                    expr.kron.push_back(expr.kron[i]);
                    expr.pref.push_back(expr.pref[i]);
                    expr.NC_ctr.push_back(expr.NC_ctr[i]);
                    expr.pref[expr.no_of_terms()-1]*=1./2.;
                    expr.sym[expr.no_of_terms()-1].set_indices(j,k,x);
                    expr.fund[expr.no_of_terms()-1].set_indices(x,a,c);
                    // first term
                    expr.pref[i]*=1./(2.*NC);
                    expr.kron[i].set_indices(j,k,true);
                    expr.kron[i].set_indices(a,c,false);
                    replacements_made=true;
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
                    expr.sym.push_back(expr.sym[i]);
                    expr.asym.push_back(expr.asym[i]);
                    expr.fund.push_back(expr.fund[i]);
                    expr.kron.push_back(expr.kron[i]);
                    expr.pref.push_back(expr.pref[i]);
                    expr.NC_ctr.push_back(expr.NC_ctr[i]);
                    expr.pref[expr.no_of_terms()-1]*=complex<double>(0.,1./2.);;
                    expr.asym[expr.no_of_terms()-1].set_indices(j,k,x);
                    expr.fund[expr.no_of_terms()-1].set_indices(x,a,c);
                    // second term
                    expr.sym.push_back(expr.sym[i]);
                    expr.asym.push_back(expr.asym[i]);
                    expr.fund.push_back(expr.fund[i]);
                    expr.kron.push_back(expr.kron[i]);
                    expr.pref.push_back(expr.pref[i]);
                    expr.NC_ctr.push_back(expr.NC_ctr[i]);
                    expr.pref[expr.no_of_terms()-1]*=1./2.;
                    expr.sym[expr.no_of_terms()-1].set_indices(j,k,x);
                    expr.fund[expr.no_of_terms()-1].set_indices(x,a,c);
                    // first term
                    expr.pref[i]*=1./(2.*NC);
                    expr.kron[i].set_indices(j,k,true);
                    expr.kron[i].set_indices(a,c,false);
                    replacements_made=true;
                }
                else it1++;
            }
            else it1++;
            if (it1>=expr.fund[i].len()-1) evaluated=true;
        }
    }
    return replacements_made;
}

// check for Tr(T_i)
void check_fund_trace(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    size_t it(0);
    while (it<fundamental.len()) {
        if (fundamental.index(it,1)==fundamental.index(it,2))
            clear_indices_and_prefactor(prefactor, kronecker, symmetric, antisymmetric, fundamental);
        else it++;
    }
}

// contract indices of Kronecker deltas
void replace_k(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated and kronecker.len()>0) {
        int ind_i=kronecker.index(it1,0), ind_j=kronecker.index(it1,1);
        if (!(is_free_index(ind_i, kronecker, symmetric, antisymmetric, fundamental))) {
            if (ind_i==ind_j) {
                if (kronecker.is_gluonic(it1)) prefactor*=(pow(NC,2)-1.); 
                else prefactor*=NC;
            }
            else {
                symmetric.find_and_rep_indices(ind_i,ind_j);
                antisymmetric.find_and_rep_indices(ind_i,ind_j);
                fundamental.find_and_rep_indices(ind_i,ind_j);
                kronecker.find_and_rep_indices(ind_i,ind_j);
            }
            kronecker.del_indices(it1);
        }
        else if (!(is_free_index(ind_j, kronecker, symmetric, antisymmetric, fundamental))) {
            symmetric.find_and_rep_indices(ind_j,ind_i);
            antisymmetric.find_and_rep_indices(ind_j,ind_i);
            fundamental.find_and_rep_indices(ind_j,ind_i);
            kronecker.find_and_rep_indices(ind_j,ind_i);
            kronecker.del_indices(it1);
        }
        else it1++;
        if (it1>=kronecker.len()) evaluated=true;
    }
}

// replace products of 2 or 3 antisymmetric structure constants
void replace_f(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated and antisymmetric.len()>0) {
        if (antisymmetric.count_index_at(antisymmetric.index(it1,0),it1)>1 or antisymmetric.count_index_at(antisymmetric.index(it1,1),it1)>1) {
            prefactor=0;
            clear_indices_and_prefactor(prefactor, kronecker, symmetric, antisymmetric, fundamental);
        }
        else if (antisymmetric.len()-it1>=2) {
            pair<size_t,size_t> itf1(it1+1,0), itf2(it1+1,0);
            bool part_evaluated(false);
            int n(0);
            while (!part_evaluated) {
                if (antisymmetric.count_index(antisymmetric.index(it1,1))>1 and antisymmetric.count_index(antisymmetric.index(it1,2))>1) {
                    itf1=antisymmetric.find_index(antisymmetric.index(it1,1),it1+1);
                    itf2=antisymmetric.find_index(antisymmetric.index(it1,2),it1+1);
                    // replace 3 f's by one f
                    if (antisymmetric.matching_indices(itf1.first,itf2.first)==1) {
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
                        contract_product(prefactor, kronecker, symmetric, antisymmetric, fundamental);
                        part_evaluated=true;
                    }
                    // replace 2 f's by kronecker - f at it1 and f at itf1
                    else if (antisymmetric.matching_indices(it1,itf1.first)>=2) {
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
                        kronecker.set_indices(antisymmetric.index(it1,0),antisymmetric.index(itf1.first,0),true);
                        antisymmetric.del_indices(itf1.first);
                        antisymmetric.del_indices(it1);
                        contract_product(prefactor, kronecker, symmetric, antisymmetric, fundamental); 
                        part_evaluated=true;
                    }
                    // replace 2 f's by kronecker - f at it1 and f at itf2
                    else if (antisymmetric.matching_indices(it1,itf2.first)>=2) {
                        while (antisymmetric.index(itf2.first,2)!=antisymmetric.index(it1,2)) antisymmetric.rotate_indices_at(itf2.first);
                        if (antisymmetric.index(it1,1)==antisymmetric.index(itf2.first,0)) {
                            prefactor*=-1;
                            antisymmetric.swap_indices_at(itf2.first,0,1);
                        }
                        else if (antisymmetric.index(it1,0)==antisymmetric.index(itf2.first,1)) {
                            prefactor*=-1;
                            antisymmetric.swap_indices_at(it1,0,1); 
                        }
                        else if (antisymmetric.index(it1,0)==antisymmetric.index(itf2.first,0)) {
                            antisymmetric.swap_indices_at(itf2.first,0,1);
                            antisymmetric.swap_indices_at(it1,0,1);
                        }
                        prefactor*=NC;
                        kronecker.set_indices(antisymmetric.index(it1,0),antisymmetric.index(itf2.first,0),true);
                        antisymmetric.del_indices(itf2.first);
                        antisymmetric.del_indices(it1);
                        contract_product(prefactor, kronecker, symmetric, antisymmetric, fundamental); 
                        part_evaluated=true;
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

// replace products of 2 or 3 symmetric structure constants
void replace_d(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated and symmetric.len()>0) {
        if (symmetric.count_index_at(symmetric.index(it1,0),it1)>1 or symmetric.count_index_at(symmetric.index(it1,1),it1)>1) {
            prefactor=0;
            clear_indices_and_prefactor(prefactor, kronecker, symmetric, antisymmetric, fundamental);
        }
        else if (symmetric.len()-it1>=2) {
            pair<size_t,size_t> itf1(it1+1,0), itf2(it1+1,0);
            bool part_evaluated(false);
            int n(0);
            while (!part_evaluated) {
                if (symmetric.count_index(symmetric.index(it1,1))>1 and symmetric.count_index(symmetric.index(it1,2))>1) {
                    itf1=symmetric.find_index(symmetric.index(it1,1),it1+1);
                    itf2=symmetric.find_index(symmetric.index(it1,2),it1+1); 
                    // replace 3 d's by one d
                    if (symmetric.matching_indices(itf1.first,itf2.first)==1) {
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
                        contract_product(prefactor, kronecker, symmetric, antisymmetric, fundamental);
                        part_evaluated=true;
                    }
                    // replace 2 d's by kronecker - d at it1 and d at itf1
                    else if (symmetric.matching_indices(it1,itf1.first)>=2) {
                        while (symmetric.index(itf1.first,1)!=symmetric.index(it1,1)) symmetric.rotate_indices_at(itf1.first);
                        if (symmetric.index(it1,2)==symmetric.index(itf1.first,0)) symmetric.swap_indices_at(itf1.first,0,2);
                        else if (symmetric.index(it1,0)==symmetric.index(itf1.first,2)) symmetric.swap_indices_at(it1,0,2); 
                        else if (symmetric.index(it1,0)==symmetric.index(itf1.first,0)) {
                            symmetric.swap_indices_at(itf1.first,0,2);
                            symmetric.swap_indices_at(it1,0,2);
                        }
                        prefactor*=(pow(NC,2)-4)/NC;
                        kronecker.set_indices(symmetric.index(it1,0),symmetric.index(itf1.first,0),true);
                        symmetric.del_indices(itf1.first);
                        symmetric.del_indices(it1);
                        contract_product(prefactor, kronecker, symmetric, antisymmetric, fundamental);
                        part_evaluated=true;
                    }
                    // replace 2 d's by kronecker - d at it1 and d at itf2
                    else if (symmetric.matching_indices(it1,itf2.first)>=2) {
                        while (symmetric.index(itf2.first,2)!=symmetric.index(it1,2)) symmetric.rotate_indices_at(itf2.first);
                        if (symmetric.index(it1,1)==symmetric.index(itf2.first,0)) symmetric.swap_indices_at(itf2.first,0,1);
                        else if (symmetric.index(it1,0)==symmetric.index(itf2.first,1)) symmetric.swap_indices_at(it1,0,2); 
                        else if (symmetric.index(it1,0)==symmetric.index(itf2.first,0)) {
                            symmetric.swap_indices_at(itf2.first,0,1);
                            symmetric.swap_indices_at(it1,0,1);
                        }
                        prefactor*=NC;
                        kronecker.set_indices(symmetric.index(it1,0),symmetric.index(itf2.first,0),true);
                        symmetric.del_indices(itf2.first);
                        symmetric.del_indices(it1);
                        contract_product(prefactor, kronecker, symmetric, antisymmetric, fundamental); 
                        part_evaluated=true;
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

// replace products of 2 antisymmetric and 1 symmetric structure constant
void replace_2fd(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated and antisymmetric.len()>0 and symmetric.len()>0) {
        pair<size_t,size_t> itf0(0,0), itf1(0,0), itf2(0,0);
        bool part_evaluated(false);
        int n(0);
        while (!part_evaluated and antisymmetric.len()>0) {
            if (antisymmetric.count_index(symmetric.index(it1,1))>0 and antisymmetric.count_index(symmetric.index(it1,2))>0) {
                itf0=antisymmetric.find_index(symmetric.index(it1,0),0);
                itf1=antisymmetric.find_index(symmetric.index(it1,1),0);
                itf2=antisymmetric.find_index(symmetric.index(it1,2),0);
                // product of d and f with 2 contracted indices yields 0
                if (itf0.first==itf1.first or itf0.first==itf2.first or itf1.first==itf2.first) {
                    prefactor=0;
                    clear_indices_and_prefactor(prefactor, kronecker, symmetric, antisymmetric, fundamental);
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
                    contract_product(prefactor, kronecker, symmetric, antisymmetric, fundamental);
                    part_evaluated=true;
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

// replace product of 2 symmetric and 1 antisymmetric structure constant
void replace_2df(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated and antisymmetric.len()>0 and symmetric.len()>0) {
        pair<size_t,size_t> itf0(0,0), itf1(0,0), itf2(0,0);
        bool part_evaluated(false);
        int n(0);
        while (!part_evaluated and symmetric.len()>0) {
            if (symmetric.count_index(antisymmetric.index(it1,1))>0 and symmetric.count_index(antisymmetric.index(it1,2))>0) {
                itf0=symmetric.find_index(antisymmetric.index(it1,0),0);
                itf1=symmetric.find_index(antisymmetric.index(it1,1),0);
                itf2=symmetric.find_index(antisymmetric.index(it1,2),0);
                // product of d and f with 2 contracted indices yields 0
                if (itf0.first==itf1.first or itf0.first == itf2.first or itf1.first==itf2.first) {
                    prefactor=0;
                    clear_indices_and_prefactor(prefactor, kronecker, symmetric, antisymmetric, fundamental);
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
                    contract_product(prefactor, kronecker, symmetric, antisymmetric, fundamental); 
                    part_evaluated=true;
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

// check if certain index is repeated in one product
bool is_free_index(int index, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    int cntr(0);
    cntr+=kronecker.count_index(index);
    cntr+=symmetric.count_index(index);
    cntr+=antisymmetric.count_index(index);
    cntr+=fundamental.count_index(index);
    if (cntr<=1) return true;
    else return false;
}

// find a free internal index
int find_free_internal_ind(int internal_ind, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    //int internal_ind(101);
    while (kronecker.count_index(internal_ind)+symmetric.count_index(internal_ind)+antisymmetric.count_index(internal_ind)+fundamental.count_index(internal_ind)!=0) internal_ind++;
    return internal_ind;
}

// delete all quantities of one product
void clear_indices_and_prefactor(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    kronecker.clear_indices();
    symmetric.clear_indices();
    antisymmetric.clear_indices();
    fundamental.clear_indices();
    prefactor=0.;
}

// check if a product vanishes
bool vanishing_expr(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental) {
    if (abs(prefactor.real())>eps and abs(prefactor.imag())>eps and kronecker.len()!=0 and symmetric.len()!=0 and antisymmetric.len()!=0 and fundamental.len()!=0 ) return true;
    else return false;
}
