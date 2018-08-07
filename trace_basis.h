#include<string>
#include<iostream>
#include<vector>
#include<algorithm>
#include<complex>
#include<cmath>


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
        size_t no_g();
        size_t no_qp();
        bool is_non_zero();
        bool vanishes();
        bool comp(trace_t& tr_t);
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
        vector<size_t> get_indices();
        bool is_connected();
        bool has_sg();
        void order();
        bool comp(trace_vec& tr_v);
        void print();
};

class trace_basis {
    
    vector<trace_vec> m_tr_basis;
    
    public:
        trace_basis(size_t n_g, size_t n_qp);
        ~trace_basis();
//         void add_basis_vec(trace_vec tr_v);
//         trace_vec& at(size_t i);
//         trace_vec at(size_t i) const;
        void remove_sg();
        size_t dim();
        void normal_order();
        vector<vector<size_t>> get_perms();
        void print();
};
