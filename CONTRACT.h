#ifndef CONTRACT_H
#define CONTRACT_H

void evaluate(colour_term& expr);
void sort(colour_term& expr);
void sort_indices(colour_term& expr);
void sort_tensors(colour_term& expr);
void add_terms(colour_term& expr);
void fully_contract();
bool replace_fund(colour_term& expr);
void replace_k();
void replace_f();
void replace_d();
void replace_2fd();
void replace_2df();
bool is_free_index(int index);
void del_all_indices();
void reset();
bool nonvanishing_expr();

#endif
