#ifndef MAIN_H
#define MAIN_H

void read_in_process(string filename, process& m_process, vector<string>& basis_str);
colour_term decompose_terms(string& input, process m_process);
string generate_out_filename(process m_process);

#endif
