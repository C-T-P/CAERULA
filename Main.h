#ifndef MAIN_H
#define MAIN_H

void read_in_process(std::string filename, diagram& process, std::vector<std::string>& basis_str);
void decompose_terms(std::string& input, colour_term& expr);

#endif
