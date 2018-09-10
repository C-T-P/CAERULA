#ifndef F_BASIS_H
#define F_BASIS_H

#include "c_basis.h"

using namespace std;

class f_type {
    
    vector<size_t> m_g;

    public:
        f_type(vector<size_t> g_inds = {});
        ~f_type();
        vector<f_type> add_one_gluon(size_t g_ind);
        vector<size_t> get_indices();
        size_t no_g();
        bool is_not_empty();
        bool vanishes();
        bool operator>(f_type& rhs);
        bool operator==(f_type rhs);
        c_term build_ct(size_t start_ind);
        void print();
};

class f_vec {

    vector<f_type> m_f_vec;

    public:
        f_vec(f_type f = f_type(vector<size_t>({})));
        ~f_vec();
        void push_back(f_type f);
        f_type& at(size_t i);
        f_type at(size_t i) const;
        vector<f_vec> add_one_gluon(size_t g_ind);
        vector<size_t> get_indices();
        size_t no_groups();
        bool is_tree_level();
        bool has_sg();
        void order();
        bool operator>(f_vec& rhs);
        c_amplitude build_ca();
        void print();
};

class f_basis : public c_basis {

    vector<f_vec> m_f_basis;
    vector<size_t> m_g_indices;
    size_t m_ng;
    
    void make_perms();
    void make_ca_basis();
    
    void remove_sg();
    void normal_order();

    public:
        f_basis(size_t n_g);
        ~f_basis();
};

#endif
