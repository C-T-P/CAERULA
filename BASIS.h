#ifndef BASIS_H
#define BASIS_H

vector<colour_term> normalise_basis(vector<colour_term> basis, int NC_order, vector<complex<double>>& normalisations);
c_matrix calc_soft_matrix(vector<colour_term> basis, int NC_order);
c_matrix calc_colour_change_matrix(vector<colour_term> basis, c_matrix soft_matrix, process m_process, unsigned int lno1, unsigned int lno2, int NC_order);
colour_term make_internal(process m_process, colour_term expr);

#endif
