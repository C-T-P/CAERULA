#ifndef CONTRACT_H
#define CONTRACT_H

void evaluate(terms& expr);
void fully_contract();
void replace_fund(terms& expr);
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
