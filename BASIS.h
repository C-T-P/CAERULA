#ifndef BASIS_H
#define BASIS_H

vector<colour_term> read_basis(string filename, process& m_process);
void read_in_process(string filename, process& m_process, vector<string>& basis_strs);
colour_term decompose_terms(string& input, process m_process);

vector<colour_term> construct_basis(int n_qp, int n_g, process& m_process, vector<vector<int>>& amp_perms);
vector<colour_term> normalise_basis(vector<colour_term> basis, int NC_order, vector<complex<double>>& normalisations);

vector<colour_term> build_qqbg_basis(int n_qp, int n_g, vector<vector<int>>& amp_perms);
vector<colour_term> arrange_qngqb_ind (vector<vector<int>>& qqb_ind_combos, vector<int>& g_indices, vector<int>& g_part, vector<vector<int>>& amp_perms);
vector<colour_term> connect_qngqb (vector<int> qqb_inds, vector<int> g_inds, vector<int> g_part, vector<vector<int>>& amp_perms);
colour_term trace_connected_qngqb (int q_ind, vector<int> g_inds, int qb_ind);

vector<colour_term> build_q_basis(int n_qp, vector<vector<int>>& amp_perms);
vector<vector<int>> get_q_ind_combinations(vector<int> q_inds, vector<int> qb_inds);
vector<colour_term> colourflow_q(vector<vector<int>> qqb_ind_combos);

vector<colour_term> build_g_basis(int n_g, vector<vector<int>>& amp_perms);
vector<vector<int>> get_arranged_g_ind_for_part(vector<int> ind, vector<int> g_partition);
vector<vector<int>> get_g_ind_combinations(int N, int K, vector<int> inds);
vector<vector<int>> arrange_g_ind(vector<int> ind, vector<int> g_partition, size_t level);
colour_term trace_connected_g(vector<int> indices);
colour_term write_tr_basis_el(vector<int> ind);
vector<vector<int>> get_g_partitions(int n);
void int_partitions(int n, vector<int>& v, vector<vector<int>>& r_v, int level);

c_matrix calc_soft_matrix(vector<colour_term> basis, int NC_order);
c_matrix calc_inv_soft_matrix(c_matrix soft_matrix);
c_matrix calc_colour_change_matrix(vector<colour_term> basis, c_matrix soft_matrix, process m_process, unsigned int lno1, unsigned int lno2, int NC_order);
colour_term make_internal(process m_process, colour_term expr);

#endif
