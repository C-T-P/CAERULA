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
    size_t no_of_terms() {
        return pref.size();
    };
    int count_index_in_term(int index, size_t t_no) {
        int counter(0);
        if (t_no<sym.size()) counter+=sym[t_no].count_index(index);
        if (t_no<asym.size()) counter+=asym[t_no].count_index(index);
        if (t_no<fund.size()) counter+=fund[t_no].count_index(index);
        if (t_no<kron.size()) counter+=kron[t_no].count_index(index);
        return counter;
    }
    colour_term multiply(colour_term ct) {
        colour_term ct_r;
        three_ind symmmetric;
        three_ind asymmmetric;
        three_ind fundamental;
        two_ind kronecker;
        complex<double> prefactor;
        
        for (size_t t_it1(0);t_it1<no_of_terms();t_it1++) {
            // check if the terms have common internal indices
            vector<vector<vector<int>>> all_ind;
            all_ind.push_back(sym[t_it1].get_all_indices());
            all_ind.push_back(asym[t_it1].get_all_indices());
            all_ind.push_back(fund[t_it1].get_all_indices());
            all_ind.push_back(kron[t_it1].get_all_indices());
            for (size_t ty_it(0);ty_it<4;ty_it++) {
                for (size_t f_it(0);f_it<all_ind[ty_it].size();f_it++) {
                    for (size_t i_it(0);i_it<all_ind[ty_it][f_it].size();i_it++) {
                        for (size_t t_it2(0);t_it2<ct.no_of_terms();t_it2++) {
                            int index=all_ind[ty_it][f_it][i_it];
                            if (ct.count_index_in_term(index,t_it2)>1) {
                                while (ct.count_index_in_term(index,t_it2)>0) index++;
                                ct.sym[t_it2].find_and_rep_indices(all_ind[ty_it][f_it][i_it],index);
                                ct.asym[t_it2].find_and_rep_indices(all_ind[ty_it][f_it][i_it],index);
                                ct.fund[t_it2].find_and_rep_indices(all_ind[ty_it][f_it][i_it],index);
                                ct.kron[t_it2].find_and_rep_indices(all_ind[ty_it][f_it][i_it],index);
                            }
                        }
                    }
                }
            }
            
            for (size_t t_it2(0);t_it2<ct.no_of_terms();t_it2++) {
                symmmetric=sym.at(t_it1);
                asymmmetric=asym.at(t_it1);
                fundamental=fund.at(t_it1);
                kronecker=kron.at(t_it1);
                prefactor=pref.at(t_it1)*ct.pref.at(t_it2);
                symmmetric.append_by(ct.sym.at(t_it2).get_all_indices());
                asymmmetric.append_by(ct.asym.at(t_it2).get_all_indices());
                fundamental.append_by(ct.fund.at(t_it2).get_all_indices());
                kronecker.append_by(ct.kron.at(t_it2));
                ct_r.sym.push_back(symmmetric);
                ct_r.asym.push_back(asymmmetric);
                ct_r.fund.push_back(fundamental);
                ct_r.kron.push_back(kronecker);
                ct_r.pref.push_back(prefactor);
                ct_r.NC_ctr.push_back(NC_ctr.at(t_it1)+ct.NC_ctr.at(t_it2));
            }
        }
        return ct_r;
    }
    colour_term cconj() {
        colour_term ct;
        ct.sym=sym;
        ct.asym=asym;
        ct.fund=fund;
        ct.kron=kron;
        ct.pref=pref;
        ct.NC_ctr=NC_ctr;
        for (size_t i(0);i<pref.size();i++) ct.pref[i]=conj(ct.pref[i]);
        return ct;
    }
    // Scalar Product < A | B > = Tr(A^+ B) (anti-linear in first argument)
    // NOTE: Returns colour term, which has to be evaluated !!!
    colour_term scprod(colour_term ct2) {
        colour_term ct1;
        ct1.sym=sym;
        ct1.asym=asym;
        ct1.fund=fund;
        ct1.kron=kron;
        ct1.pref=pref;
        ct1.NC_ctr=NC_ctr;
        for (size_t t_it(0);t_it<ct1.no_of_terms();t_it++)
            for (size_t f_it(0);f_it<ct1.fund.at(t_it).len();f_it++)
                ct1.fund.at(t_it).swap_indices_at(f_it,1,2);
        return ct1.cconj().multiply(ct2);
    }
    colour_term term(size_t termno) {
        colour_term ct;
        ct.delete_all_terms();
        ct.sym.push_back(sym[termno]);
        ct.asym.push_back(asym[termno]);
        ct.fund.push_back(fund[termno]);
        ct.kron.push_back(kron[termno]);
        ct.pref.push_back(pref[termno]);
        ct.NC_ctr.push_back(NC_ctr[termno]);
        return ct;
    }
    void delete_term(size_t j) {
        sym.erase(sym.begin()+j);
        asym.erase(asym.begin()+j);
        fund.erase(fund.begin()+j);
        kron.erase(kron.begin()+j);
        pref.erase(pref.begin()+j);
        NC_ctr.erase(NC_ctr.begin()+j);
    }
    void delete_all_terms() {
        sym.clear();
        asym.clear();
        fund.clear();
        kron.clear();
        pref.clear();
        NC_ctr.clear();
    }
    string build_string() {
        string return_str;
        for (size_t it(0);it<no_of_terms();it++) {
            string str="";
            if (pref.at(it).real()!=0. || pref.at(it).imag()!=0.) str+="c_["+to_string(pref.at(it).real())+","+to_string(pref.at(it).imag())+"]";
            if (sym.at(it).len()>0) 
                for (size_t i(0); i<sym.at(it).len(); i++) str+="*d_["+to_string(sym.at(it).index(i,0))+","+to_string(sym.at(it).index(i,1))+","+to_string(sym.at(it).index(i,2))+"]";
            if (asym.at(it).len()>0) 
                for (size_t i(0); i<asym.at(it).len(); i++) str+="*f_["+to_string(asym.at(it).index(i,0))+","+to_string(asym.at(it).index(i,1))+","+to_string(asym.at(it).index(i,2))+"]"; 
            if (fund.at(it).len()>0)
                for (size_t i(0); i<fund.at(it).len(); i++) str+="*t_["+to_string(fund.at(it).index(i,0))+","+to_string(fund.at(it).index(i,1))+","+to_string(fund.at(it).index(i,2))+"]"; 
            if (kron.at(it).len()>0)
                for (size_t i(0); i<kron.at(it).len(); i++) str+="*k_["+to_string(kron.at(it).index(i,0))+","+to_string(kron.at(it).index(i,1))+"]";
            if (str!="" && return_str!="") return_str+="+";
            if (str!="") return_str+=str;
        }
        return return_str;
    }
    complex<double> build_complex() {
        if (no_of_terms()==0) {
            return complex<double>(0.,0.);
        }
        else if (no_of_terms()==1) {
            if (sym.at(0).len()==0 and asym.at(0).len()==0 and fund.at(0).len()==0 and kron.at(0).len()==0) return pref.at(0);
            else return complex<double>(NAN,NAN);
        }
        else return complex<double>(NAN,NAN);
    }
};

#endif
