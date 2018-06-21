#ifndef BASIS_H
#define BASIS_H

std::vector<std::vector<std::complex<float>>> calc_soft_matrix(std::vector<colour_term> basis, int NC_order);
std::vector<std::vector<std::complex<float>>> calc_colour_change_matrix(std::vector<colour_term> basis, std::vector<std::vector<std::complex<float>>> soft_matrix, diagram process, unsigned int lno1, unsigned int lno2, int NC_order);
colour_term make_internal(diagram process, colour_term expr);

#endif
