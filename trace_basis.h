#include<string>
#include<iostream>
#include<vector>
#include<algorithm>
#include<complex>
#include<cmath>
#include "colourtools.h"


using namespace std;

class trace_t {

    size_t m_qb, m_q;
    vector<size_t> m_g;

    public:
        trace_t(vector<size_t> g_inds = {}, size_t qb_ind = 0, size_t q_ind = 0);
        trace_t(size_t qb_ind, size_t q_ind);
        ~trace_t();
        vector<trace_t> add_one_gluon(size_t g_ind);
        vector<size_t> get_indices();
        trace_t conj();
        size_t no_g();
        size_t no_qp();
        bool is_non_zero();
        bool vanishes();
        bool comp(trace_t& tr_t);
        bool operator==(trace_t& rhs);
        colour_term build_ct();
        void print();
};

class trace_vec {
    
    vector<trace_t> m_tr_vec;
    
    public:
        trace_vec(trace_t tr = trace_t({},0,0));
        ~trace_vec();
        void push_back(trace_t tr);
        trace_t& at(size_t i);
        trace_t at(size_t i) const;
        vector<trace_vec> add_one_gluon(size_t g_ind);
        vector<trace_vec> conjugates();
        vector<size_t> get_indices();
        size_t no_groups();
        bool is_connected();
        bool has_sg();
        void order();
        bool comp(trace_vec& tr_v);
        bool operator==(trace_vec& rhs);
        colour_term build_ct();
        void print();
};

class trace_basis {
    
    vector<trace_vec> m_tr_basis;
    vector<size_t> m_g_indices, m_q_indices, m_qb_indices;
    size_t m_ng, m_nqp;
    
    public:
        trace_basis(size_t n_g, size_t n_qp);
        ~trace_basis();
    
        void remove_sg();
        void remove_conj();
        void normal_order();
    
        size_t dim();
        process proc();
        vector<vector<size_t>> perms();
        vector<colour_term> ct_basis();
        void print();
};
