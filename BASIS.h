#ifndef BASIS_H
#define BASIS_H

sq_matrix calc_soft_matrix(vector<colour_term> basis, int NC_order);
sq_matrix calc_inv_soft_matrix(sq_matrix soft_matrix);
sq_matrix calc_colour_change_matrix(vector<colour_term> basis, sq_matrix soft_matrix, diagram process, unsigned int lno1, unsigned int lno2, int NC_order);
colour_term make_internal(diagram process, colour_term expr);

#endif
