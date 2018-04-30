#ifndef MAIN_H
#define MAIN_H

std::string factorise(std::string string);
void decompose_terms(std::string& input, terms& expr);
void delete_term(int j, terms& expr);
std::string build_string(terms& expr);

#endif
