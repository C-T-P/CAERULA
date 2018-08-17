#ifndef COLOURTOOLS_H
#define COLOURTOOLS_H

#include<cstring>
#include<iostream>
#include<vector>
#include<algorithm>
#include<complex>
#include<cmath>
using namespace std;

static double NC(3.); // number of colours
static double CF((NC*NC-1)/(2.*NC)); // Fundamental Casimir
static double TR(1./2.);

class process {
    /* 
     incoming and outgoing legs: an index is assigned to each leg (first component) and the particle id according to pdg is stored in second component
     */
    vector<pair<size_t,string>> m_in_legs;
    vector<pair<size_t,string>> m_out_legs;
    public:
        process();
        ~process();
        void add_in_leg(string ptcl);
        void add_out_leg(string ptcl);
        void delete_all_legs();
        size_t no_of_legs();
        pair<size_t,string> leg(size_t index);
        bool is_in_leg(size_t index);
};


class delta {
    public:
        size_t m_i, m_j;
        bool m_adj;
        delta(size_t i, size_t j, bool adj);
        ~delta();
        bool is_free(size_t ind);
        string build_string();
};
class fundamental {
    public:
        size_t m_a, m_i, m_j;
        fundamental(size_t a, size_t i, size_t j);
        ~fundamental();
        bool is_free(size_t ind);
        string build_string();
};
class antisymmetric {
    public:
        size_t m_a, m_b, m_c;
        antisymmetric(size_t a, size_t b, size_t c);
        ~antisymmetric();
        bool is_free(size_t ind);
        string build_string();
};
class symmetric {
    public:
        size_t m_a, m_b, m_c;
        symmetric(size_t a, size_t b, size_t c);
        ~symmetric();
        bool is_free(size_t ind);
        string build_string();
};
class c_term {
    complex<double> m_cnum;
    vector<delta> m_k_vec;
    vector<fundamental> m_t_vec;
    vector<antisymmetric> m_f_vec;
    vector<symmetric> m_d_vec;
    int m_NC_order;
    size_t m_fi;
    public:
        c_term();
        c_term(delta k, fundamental t, antisymmetric f, symmetric d, complex<double> c = complex<double>(0.,0.), int NC_order = 0);
        ~c_term();
    
        void push_back(c_term ct);
        void push_back(delta k);
        void push_back(fundamental t);
        void push_back(antisymmetric f);
        void push_back(symmetric d);
        void cnumber(complex<double> c);
        void NC_order(int NC_o);
    
        void evaluate_deltas();
        void shift_inds(size_t by, bool all);
        c_term hconj();
        c_term operator*(c_term ct);
        complex<double> result();
        void clear();
    
        string build_string();
        void print();
    
        friend class c_amplitude;
};
class c_amplitude {
    vector<c_term> m_cterm_vec;
    complex<double> m_result;
    public:
        c_amplitude();
        c_amplitude(c_term ct);
        c_amplitude(string expr);
        ~c_amplitude();
    
        void add(c_term ct);
        void push_back(c_amplitude ca);
        c_amplitude hconj();
        c_amplitude shift_to_internal(size_t by);
        c_amplitude operator*(complex<double> z);
        c_amplitude operator*(c_amplitude ca);
        complex<double> scprod(c_amplitude ca, size_t up_to_NC = INT_MAX);
        void clear();
    
        void evaluate();
        void evaluate(size_t up_to_NC);
        complex<double> result();
    
        size_t no_of_terms();
        string build_string();
        void print();
};

#endif
