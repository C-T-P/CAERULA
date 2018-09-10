#ifndef TRACE_BASIS_H
#define TRACE_BASIS_H

#include "c_basis.h"
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
        bool is_not_empty();
        bool vanishes();
        bool operator>(trace_t& rhs);
        bool operator==(trace_t& rhs);
        c_amplitude build_ca(size_t start_ind);
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
        bool is_tree_level();
        bool has_sg();
        void order();
        bool operator>(trace_vec& rhs);
        bool operator==(trace_vec& rhs);
        c_amplitude build_ca();
        void print();
};

class trace_basis : public c_basis {
    
    vector<trace_vec> m_tr_basis;
    vector<size_t> m_g_indices, m_q_indices, m_qb_indices;
    size_t m_ng, m_nqp;
    
    void make_perms();
    void make_ca_basis();
    
    void remove_sg();
    void remove_conj();
    void normal_order();
    
    public:
        trace_basis(size_t n_g, size_t n_qp);
        ~trace_basis();
    
        friend class multiplet_basis;
};

#endif
