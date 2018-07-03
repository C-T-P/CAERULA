#ifndef BASIS_H
#define BASIS_H

vector<vector<complex<double>>> calc_soft_matrix(vector<colour_term> basis, int NC_order);
vector<vector<complex<double>>> calc_inv_soft_matrix(vector<vector<complex<double>>> soft_matrix);
vector<vector<complex<double>>> calc_colour_change_matrix(vector<colour_term> basis, vector<vector<complex<double>>> soft_matrix, diagram process, unsigned int lno1, unsigned int lno2, int NC_order, bool inv_mult);
colour_term make_internal(diagram process, colour_term expr);
void save_colour_to_file(vector<vector<vector<complex<double>>>> colour_change_matrices, vector<vector<complex<double>>> soft_matrix, diagram process);

#endif
