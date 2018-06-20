#ifndef CONTRACT_H
#define CONTRACT_H

std::complex<float> evaluate_colour_term_to_order(colour_term& expr, unsigned int NC_order);
void simplify_colour_term(colour_term& expr);
void sort_colour_term(colour_term& expr);
void sort_indices(colour_term& expr);
void sort_tensors(colour_term& expr);
void add_terms(colour_term& expr);
void contract_product(std::complex<float>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_d_and_f_by_t(colour_term& expr);
bool replace_fund(colour_term& expr);
void check_fund_trace(std::complex<float>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_k(std::complex<float>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_f(std::complex<float>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_d(std::complex<float>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_2fd(std::complex<float>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void replace_2df(std::complex<float>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
bool is_free_index(int index, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
void clear_indices_and_prefactor(std::complex<float>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
int find_free_internal_ind(int internal_ind, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);
bool vanishing_expr(std::complex<float>& prefactor, two_ind& kronecker, three_ind& symmetric, three_ind& antisymmetric, three_ind& fundamental);

#endif
