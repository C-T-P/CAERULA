#ifndef FCONTRACT_H
#define FCONTRACT_H

complex<double> fast_evaluate_colour_term_to_order(colour_term& expr, unsigned int NC_order);
void replace_fund_by_fierz(colour_term& expr, int NC_order);
void fast_replace_k(complex<double>& prefactor, two_ind& kronecker);

#endif
