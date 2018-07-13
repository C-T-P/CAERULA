#ifndef TENSORTOOLS_H
#define TENSORTOOLS_H

#include<gsl/gsl_linalg.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<string>
#include<iostream>
#include<vector>
#include<algorithm>
#include<complex>
#include<cmath>
using namespace std;

static double NC(3.); // number of colours

class process {
    /* 
     incoming and outgoing legs: an index is assigned to each leg (first component) and the particle id according to pdg is stored in second component
     */
    vector<pair<unsigned int,int>> in_legs;
    vector<pair<unsigned int,int>> out_legs;
    public:
        process();
        ~process();
        void add_in_leg(int ptcl_id);
        void add_out_leg(int ptcl_id);
        void delete_all_legs();
        unsigned int no_of_legs();
        pair<unsigned int,int> leg(unsigned int index);
        bool is_in_leg(unsigned int index);
};

class three_ind {
    vector<vector<int>> ind;
    public:
        three_ind();
        ~three_ind();
        void set_indices(int i, int j, int k);
        void append_by(vector<vector<int>> ind_v);
        vector<vector<int>> get_all_indices();
        void del_indices(size_t it);
        void clear_indices();
        void find_and_rep_indices(int old_ind, int new_ind);
        int matching_indices(size_t it1, size_t it2);
        void swap_indices_at(size_t pos, size_t it1, size_t it2);
        void rotate_indices_at(size_t it);
        void sort_indices_at(size_t it);
        int count_index(int index);
        int count_index_at(int index, size_t pos);
        int index(size_t it0, size_t it1);
        size_t len();
        void sort_list();
        pair<size_t,size_t> find_index(int index, size_t start);
        bool has_index_at(int index, size_t it);
};
class two_ind {
    struct flagged_indices {
        vector<int> ind;
        bool gluonic_k;
        flagged_indices(vector<int> k_indices, bool is_gluonic) {
            ind=k_indices;
            gluonic_k=is_gluonic;
        }
    };
    vector<flagged_indices> indices;
    public:
        two_ind();
        ~two_ind();
        void set_indices(int i, int j, bool gluonic);
        bool is_gluonic(size_t it);
        void append_by(two_ind tensor);
        vector<vector<int>> get_all_indices();
        vector<bool> get_all_flags();
        void del_indices(int it);
        void clear_indices();
        void find_and_rep_indices(int old_ind, int new_ind);
        int count_index(int index);
        void swap_indices_at(size_t pos);
        void sort_indices_at(size_t it);
        int index(size_t it0, size_t it1);
        size_t len();
        void sort_list();
        pair<size_t,size_t> find_index(int index, size_t start);
};
struct colour_term {
    vector<three_ind> sym;
    vector<three_ind> asym;
    vector<three_ind> fund;
    vector<two_ind> kron;
    vector<complex<double>> pref;
    vector<int> NC_ctr;
    
    // member functions
    colour_term();
    ~colour_term();
    size_t no_of_terms();
    int count_index_in_term(int index, size_t t_no);
    colour_term multiply(colour_term ct);
    colour_term cconj();
    // Scalar Product < A | B > = Tr(A^+ B) (anti-linear in first argument)
    // NOTE: Returns colour term, which has to be evaluated !!!
    colour_term scprod(colour_term ct2);
    colour_term term(size_t termno);
    void delete_term(size_t j);
    void delete_all_terms();
    string build_string();
    complex<double> build_complex();
};

#endif
