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
#include<math.h>
using namespace std;

class diagram {
    /* 
     incoming and outgoing legs: an index is assigned to each leg (first component) and the particle id according to pdg is stored in secont component
     */
    std::vector<std::pair<unsigned int,int>> in_legs;
    std::vector<std::pair<unsigned int,int>> out_legs;
    public:
        void add_in_leg(int ptcl_id) {
            int index=in_legs.size()+1;
            if (out_legs.size()>0) {
                for (size_t i(0);i<out_legs.size();i++) {
                    out_legs.at(i).first+=1;
                }
            }
            in_legs.push_back(std::pair<int,int>(index,ptcl_id));
        }
        void add_out_leg(int ptcl_id) {
            int index=in_legs.size()+out_legs.size()+1;
            in_legs.push_back(std::pair<int,int>(index,ptcl_id));
        }
        void delete_all_legs() {
            in_legs.clear();
            out_legs.clear();
        }
        unsigned int no_of_legs() {
            return in_legs.size()+out_legs.size()+1;
        }
        // leg numbering starts at 1 !!!
        std::pair<unsigned int,int> leg(unsigned int index) {
            index-=1;
            if (index<in_legs.size()) return in_legs.at(index);
            else if (index<out_legs.size()) return in_legs.at(index);
            else {
                cerr << "Leg " << index << " does not exist in diagram." << endl;
                return std::pair<int,int>(0,0);
            }
        }
        bool is_in_leg(unsigned int index) {
            index-=1;
            if (index<in_legs.size()) return true;
            else return false;
        }
};

class three_ind {
    std::vector<std::vector<int>> ind;
    public:
        void set_indices(int i, int j, int k) {
            std::vector<int> new_ind {i, j, k};
            ind.push_back(new_ind);
        }
        void append_by(std::vector<std::vector<int>> ind_v) {
            ind.insert(ind.end(),ind_v.begin(),ind_v.end());
        }
        std::vector<std::vector<int>> get_all_indices() {
            return ind;
        }
        void del_indices(size_t it) {
            ind.erase(ind.begin()+it);
        }
        void clear_indices() {
            ind.clear();
        }
        void find_and_rep_indices(int old_ind, int new_ind) {
            for (size_t it(0); it<ind.size(); ++it) replace(ind[it].begin(), ind[it].end(), old_ind, new_ind);
        }
        int matching_indices(size_t it1, size_t it2) {
            int cntr(0);
            if (it1<ind.size() && it2<ind.size())
                for (size_t i(0);i<ind[it2].size();i++) cntr+=count(ind[it1].begin(), ind[it1].end()+1, ind[it2][i]);
            return cntr;
        }
        void swap_indices_at(size_t pos, size_t it1, size_t it2) {
            int dummy=ind[pos][it1];
            ind[pos][it1]=ind[pos][it2];
            ind[pos][it2]=dummy;
        }
        void rotate_indices_at(size_t it) {
            std::rotate(ind[it].begin(), ind[it].begin()+1,ind[it].end());
        }
        int count_index(int index) {
            int cntr(0);
            for (size_t it(0); it<ind.size(); ++it) cntr+=count(ind[it].begin(), ind[it].end(), index);
            return cntr;
        }
        int count_index_at(int index, size_t pos) {
            int cntr(0);
            cntr+=count(ind[pos].begin(), ind[pos].end()+1, index);
            return cntr;
        }
        int index(size_t it0, size_t it1) {
            return (ind.at(it0)).at(it1);
        }
        size_t len() {
            return ind.size();
        }
        void sort_list() {            
            std::sort(ind.begin(), ind.end(),[](const std::vector<int>& ind1, const std::vector<int>& ind2) {
                if (ind1[0]>ind2[0]) return false;
                else if (ind1[0]==ind2[0]) {
                    if (ind1[1]>ind2[1]) return false;
                    else if (ind1[1]==ind2[1]) {
                        if (ind1[2]>ind2[2]) return false;
                        else return true;
                    }
                    else return true;
                }
                else return true;
            });
        }
        std::pair<size_t,size_t> find_index(int index, size_t start) {
            if (ind.size()>0) {
                size_t it(start);
                size_t f=find(ind[start].begin(), ind[start].end()+1, index)-ind[start].begin();
                while (it<ind.size() && f>=3) {
                    if ((f=find(ind[it].begin(), ind[it].end()+1, index)-ind[it].begin())>=3) it++;
                }
                return std::pair<size_t,size_t>(it,f);
            }
            else return std::pair<size_t,size_t>(1,3);
        }
        bool has_index_at(int index, size_t it) {
            if (find(ind[it].begin(), ind[it].end()+1, index)<ind[it].end()) return true;
            else return false;
        }
};
class two_ind {
    std::vector<std::vector<int>> ind;
    public:
        void set_indices(int i, int j) {
            std::vector<int> new_ind {i, j};
            ind.push_back(new_ind);
        }
        void append_by(std::vector<std::vector<int>> ind_v) {
            ind.insert(ind.end(),ind_v.begin(),ind_v.end());
        }
        std::vector<std::vector<int>> get_all_indices() {
            return ind;
        }
        void del_indices(int it) {
            ind.erase(ind.begin()+it);
        }
        void clear_indices() {
            ind.clear();
        }
        void find_and_rep_indices(int old_ind, int new_ind) {
            for (size_t it(0); it<ind.size(); ++it) replace(ind[it].begin(), ind[it].end(), old_ind, new_ind);
        }
        int count_index(int index) {
            int cntr(0);
            for (size_t it(0); it<ind.size(); ++it) cntr+=count(ind[it].begin(), ind[it].end(), index);
            return cntr;
        }
        void swap_indices_at(size_t pos) {
            int dummy=ind[pos][1];
            ind[pos][1]=ind[pos][0];
            ind[pos][0]=dummy;
        }
        int index(size_t it0, size_t it1) {
            return (ind.at(it0)).at(it1);
        }
        size_t len() {
            return ind.size();
        }
        void sort_list() {            
            std::sort(ind.begin(), ind.end(),[](const std::vector<int>& ind1, const std::vector<int>& ind2) {
                if (ind1[0]>ind2[0]) return false;
                else if (ind1[0]==ind2[0]) {
                    if (ind1[1]>ind2[1]) return false;
                    else return true;
                }
                else return true;
            });
        }
        std::pair<size_t,size_t> find_index(int index, size_t start) {
            if (ind.size()>0) {
                size_t it(start+1);
                size_t f=find(ind[start].begin(), ind[start].end(), index)-ind[start].begin();
                while (it<ind.size() && f>=2) {
                    f=find(ind[it].begin(), ind[it].end(), index)-ind[it].begin();
                    ++it;
                }
                return std::pair<size_t,size_t>(it-1,f);
            }
            else return std::pair<size_t,size_t>(1,2);
        }
};
struct colour_term {
    std::vector<three_ind> sym;
    std::vector<three_ind> asym;
    std::vector<three_ind> fund;
    std::vector<two_ind> kron;
    std::vector<std::complex<float>> pref;
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
        std::complex<float> prefactor;
        
        for (size_t t_it1(0);t_it1<no_of_terms();t_it1++) {
            // check if the terms have common internal indices
            std::vector<std::vector<std::vector<int>>> all_ind;
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
                prefactor=pref.at(t_it1);
                symmmetric.append_by(ct.sym.at(t_it2).get_all_indices());
                asymmmetric.append_by(ct.asym.at(t_it2).get_all_indices());
                fundamental.append_by(ct.fund.at(t_it2).get_all_indices());
                kronecker.append_by(ct.kron.at(t_it2).get_all_indices());
                prefactor*=ct.pref.at(t_it2);
                ct_r.sym.push_back(symmmetric);
                ct_r.asym.push_back(asymmmetric);
                ct_r.fund.push_back(fundamental);
                ct_r.kron.push_back(kronecker);
                ct_r.pref.push_back(prefactor);
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
        for (size_t i(0);i<pref.size();i++) ct.pref[i]=conj(ct.pref[i]);
        return ct;
    }
    void delete_term(size_t j) {
        sym.erase(sym.begin()+j);
        asym.erase(asym.begin()+j);
        fund.erase(fund.begin()+j);
        kron.erase(kron.begin()+j);
        pref.erase(pref.begin()+j);
    }
    void delete_all_terms() {
        sym.clear();
        asym.clear();
        fund.clear();
        kron.clear();
        pref.clear();
    }
    std::string build_string() {
        std::string return_str;
        for (size_t it(0);it<no_of_terms();it++) {
            std::string str="";
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
    std::complex<float> build_complex() {
        if (no_of_terms()==0) {
            return std::complex<float>(0.,0.);
        }
        else if (no_of_terms()==1) {
            if (sym.at(0).len()==0 and asym.at(0).len()==0 and fund.at(0).len()==0 and kron.at(0).len()==0) return pref.at(0);
            else return std::complex<float>(NAN,NAN);
        }
        else return std::complex<float>(NAN,NAN);
    }
};

#endif