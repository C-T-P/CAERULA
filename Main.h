#ifndef MAIN_H
#define MAIN_H

void read_in_process(string filename, diagram& process, vector<string>& basis_str);
colour_term decompose_terms(string& input, diagram process);
string generate_out_filename(diagram process);

#endif
