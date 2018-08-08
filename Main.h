#ifndef MAIN_H
#define MAIN_H

string generate_out_filename(process m_process);
void colour_calc(vector<colour_term>& basis, process& m_process, vector<vector<size_t>>& amp_perms, const int& NC_order, string& out_filename, bool& multiply_with_inv_sm, bool& norm_b);
void construct_bcm(vector<colour_term> multiplet_basis, size_t n_qp, size_t n_g, int NC_order);
void run_error();

#endif
