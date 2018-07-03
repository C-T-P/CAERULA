#ifndef CONTRACT_H
#define CONTRACT_H

complex<double> evaluate_colour_term_to_order(colour_term& expr, unsigned int NC_order);
void simplify_colour_term(colour_term& expr);
void sort_colour_term(colour_term& expr);
void sort_indices(colour_term& expr);
void sort_tensors(colour_term& expr);
void add_terms(colour_term& expr);
void contract_product(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_d_and_f_by_t(colour_term& expr);
bool replace_fund(colour_term& expr, int NC_order);
void check_fund_trace(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_k(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_f(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_d(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_2fd(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_2df(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
bool is_free_index(int index, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void clear_indices_and_prefactor(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
int find_free_internal_ind(int internal_ind, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
bool vanishing_expr(complex<double>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);

#endif
