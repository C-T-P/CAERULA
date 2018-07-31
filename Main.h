#ifndef MAIN_H
#define MAIN_H

string generate_out_filename(process m_process);
void colour_calc(vector<colour_term> basis, process m_process, const int NC_order, string out_filename);
void run_error();

#endif
