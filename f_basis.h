#include<string>
#include<iostream>
#include<vector>
#include<algorithm>
#include<complex>
#include<cmath>
#include "colourtools.h"


using namespace std;

class f_type {
    
    vector<size_t> m_g;

    public:
        f_type(vector<size_t> g_inds = {});
        ~f_type();
        vector<f_type> add_one_gluon(size_t g_ind);
        vector<size_t> get_indices();
        size_t no_g();
        bool is_non_zero();
        bool vanishes();
        bool comp(f_type& f_t);
        bool operator==(f_type rhs);
        colour_term build_ct();
        void print();
};

class f_vec {

    vector<f_type> m_f_vec;

    public:
        f_vec(f_type tr = f_type(vector<size_t>({})));
        ~f_vec();
        void push_back(f_type tr);
        f_type& at(size_t i);
        f_type at(size_t i) const;
        vector<f_vec> add_one_gluon(size_t g_ind);
        vector<size_t> get_indices();
        size_t no_groups();
        bool is_connected();
        bool has_sg();
        void order();
        bool comp(f_vec& tr_v);
        colour_term build_ct();
        void print();
};

class f_basis {

    vector<f_vec> m_f_basis;
    vector<size_t> m_g_indices;
    size_t m_ng;

    public:
        f_basis(size_t n_g);
        ~f_basis();

        void remove_sg();
        void normal_order();

        size_t dim();
        process proc();
        vector<vector<size_t>> perms();
        vector<colour_term> ct_basis();
        void print();
};
