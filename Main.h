#ifndef MAIN_H
#define MAIN_H

void read_in_process(std::string filename, diagram& process, std::vector<std::string>& basis_str);
colour_term decompose_terms(std::string& input, diagram process);

#endif
