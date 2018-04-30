#ifndef BASIS_H
#define BASIS_H

void fully_simplify(terms& expr);
void sort_indices(terms& expr);
void sort_tensors(terms& expr);
void simplify_terms(terms& expr);

#endif
