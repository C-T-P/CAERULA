#include<string>
#include<iostream>
#include<vector>
#include<algorithm>
#include<math.h>
using namespace std;

void evaluate(string expr);
void replace_k();
void replace_f();
void replace_d();
void replace_2fd();
void replace_2df();
bool is_free_index(int index);
void del_all_indices();
void print_expr();

class three_ind {
    std::vector<std::vector<int>> ind;
    public:
        void set_indices(int i, int j, int k) {
            std::vector<int> new_ind {i, j, k};
            ind.push_back(new_ind);
        }
        void del_indices(size_t it) {
            ind.erase(ind.begin()+it);
        }
        void clear_indices() {
            ind.clear();
        }
        void find_and_rep_indices(int old_ind, int new_ind) {
            for (int it(0); it<ind.size(); ++it) replace(ind[it].begin(), ind[it].end(), old_ind, new_ind);
        }
        int matching_indices(size_t it1, size_t it2) {
            int cntr(0);
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
        void del_indices(int it) {
            ind.erase(ind.begin()+it);
        }
        void clear_indices() {
            ind.clear();
        }
        void find_and_rep_indices(int old_ind, int new_ind) {
            for (int it(0); it<ind.size(); ++it) replace(ind[it].begin(), ind[it].end(), old_ind, new_ind);
        }
        int count_index(int index) {
            int cntr(0);
            for (size_t it(0); it<ind.size(); ++it) cntr+=count(ind[it].begin(), ind[it].end(), index);
            return cntr;
        }
        int index(size_t it0, size_t it1) {
            return (ind.at(it0)).at(it1);
        }
        size_t len() {
            return ind.size();
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

three_ind symmetric;
three_ind symmetric_eval;
three_ind antisymmetric;
three_ind antisymmetric_eval;
two_ind kronecker;
two_ind kronecker_eval;
float prefactor = 1.0;
static int NC = 3; // number of colours

int main(int argc, char **argv) {
    std::string expr;
    for (int i=1;i<argc;i++) expr+=argv[i];
    evaluate(expr);
    
    /* replace this by overall routine */
    for (int i(0); i<3; i++) {
        replace_k();
        replace_f();
        replace_d();
        replace_2fd();
        replace_2df();
    }
    /***********************************/
    
    cout << expr;
    print_expr();
}
void evaluate(string expr) {
    for (size_t i(0), mpos(expr.find('*'));mpos!=std::string::npos || expr.length()>0;mpos=expr.find('*')) {
        ++i;
        std::string factor;
        if (mpos==std::string::npos) {
            factor=expr;
            expr="";
        }
        else {
            factor=expr.substr(0,mpos);
            expr=expr.substr(mpos+1);
        }
        if(factor.find("f_[")==0 && factor[factor.length()-1]==']') {
            size_t c1pos(factor.find(','));
            if (c1pos==std::string::npos)
                cout << "Invalid number of indices for f." << endl;
            size_t c2pos(factor.find(',',c1pos+1));
            if (c2pos==std::string::npos || factor.find(',',c2pos+1)!=std::string::npos)
                cout << "Invalid number of indices for f." << endl;
            antisymmetric.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
        }
        else if(factor.find("d_[")==0 && factor[factor.length()-1]==']') {
            size_t c1pos(factor.find(','));
            if (c1pos==std::string::npos)
                cout << "Invalid number of indices for d." << endl;
            size_t c2pos(factor.find(',',c1pos+1));
            if (c2pos==std::string::npos || factor.find(',',c2pos+1)!=std::string::npos)
                cout << "Invalid number of indices for d." << endl;
            symmetric.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
        }
        else if (factor.find("k_[")==0 && factor[factor.length()-1]==']') {
            size_t cpos(factor.find(','));
            if (cpos==std::string::npos || factor.find(',',cpos+1)!=std::string::npos)
                cout << "Invalid number of indices for k." << endl;
            kronecker.set_indices(stoi(factor.substr(3,cpos-3)),stoi(factor.substr(cpos+1,factor.length()-cpos-2)));
        }
    }
}
void replace_k() {
    while(kronecker.len()>0) {
        int ind_i=kronecker.index(0,0), ind_j=kronecker.index(0,1);
        if (!(is_free_index(ind_i))) {
            if (ind_i==ind_j) {
                prefactor = prefactor*NC;
            }
            else {
                symmetric.find_and_rep_indices(ind_i,ind_j);
                antisymmetric.find_and_rep_indices(ind_i,ind_j);
                kronecker.find_and_rep_indices(ind_i,ind_j);
            }
        }
        else if (!(is_free_index(ind_j))) {
            symmetric.find_and_rep_indices(ind_j,ind_i);
            antisymmetric.find_and_rep_indices(ind_j,ind_i);
            kronecker.find_and_rep_indices(ind_j,ind_i);
        }
        else {
            kronecker_eval.set_indices(ind_i,ind_j);
        }
        kronecker.del_indices(0);
    }
}
void replace_f() {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated && antisymmetric.len()>0) {
        if (antisymmetric.count_index_at(antisymmetric.index(it1,0),it1)>1 || antisymmetric.count_index_at(antisymmetric.index(it1,1),it1)>1) {
            prefactor=0;
            del_all_indices();
        }
        else if (antisymmetric.len()-it1>=2) {
            std::pair<size_t,size_t> itf1(it1+1,0), itf2(it1+1,0);
            bool part_evaluated(false);
            int n(0);
            while (!part_evaluated) {
                if (antisymmetric.count_index(antisymmetric.index(it1,1))>1 && antisymmetric.count_index(antisymmetric.index(it1,2))>1) {
                    itf1=antisymmetric.find_index(antisymmetric.index(it1,1),itf1.first);
                    itf2=antisymmetric.find_index(antisymmetric.index(it1,2),itf2.first);
                    // replace 2 f's by kronecker
                    if (itf1.first==itf2.first) {
                        if (itf1.second!=1) {
                            antisymmetric.swap_indices_at(itf1.first,itf1.second,1);
                            prefactor*=-1;
                        }
                        if (itf2.second!=2) {
                            antisymmetric.swap_indices_at(itf2.first,itf2.second,2);
                            prefactor*=-1;
                        }
                        prefactor*=NC;
                        kronecker.set_indices(antisymmetric.index(it1,0),antisymmetric.index(itf1.first,0));
                        antisymmetric.del_indices(itf1.first);
                        antisymmetric.del_indices(it1);
                        part_evaluated=true;
                    }
                    // replace 3 f's by one f
                    else if (antisymmetric.matching_indices(itf1.first,itf2.first)>=1) {
                        if (itf1.second!=1) {
                            antisymmetric.swap_indices_at(itf1.first,itf1.second,1);
                            prefactor*=-1;
                        }
                        if (itf2.second!=2) {
                            antisymmetric.swap_indices_at(itf2.first,itf2.second,2);
                            prefactor*=-1;
                        }
                        while (antisymmetric.count_index_at(antisymmetric.index(itf1.first,2),itf2.first)==0) {
                            antisymmetric.swap_indices_at(itf1.first,0,2);
                            prefactor*=-1;
                        }
                        while (antisymmetric.index(itf1.first,2)!=antisymmetric.index(itf2.first,1)) {
                            antisymmetric.swap_indices_at(itf2.first,0,1);
                            prefactor*=-1;
                        }
                        antisymmetric.set_indices(antisymmetric.index(it1,0),antisymmetric.index(itf1.first,0),antisymmetric.index(itf2.first,0));
                        prefactor*=-NC/2.;
                        if (itf2.first>itf1.first) {
                            antisymmetric.del_indices(itf2.first);
                            antisymmetric.del_indices(itf1.first);
                        }
                        else {
                            antisymmetric.del_indices(itf1.first);
                            antisymmetric.del_indices(itf2.first);
                        }
                        antisymmetric.del_indices(it1);
                        part_evaluated=true;
                    }
                    // for replacement of 4 f's insert another else if, find a 4th f and check indices
                    // can also implement replacement of 2 f's with only one contracted index if sums are implemented
                    else it1++;
                }
                else if (n<3) {
                    antisymmetric.rotate_indices_at(it1);
                    n++;
                }
                else {
                    it1++;
                    part_evaluated=true;
                }
            }
        }
        else it1++;
        if (it1>=antisymmetric.len()) evaluated=true;
    }
}
void replace_d() {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated && symmetric.len()>0) {
        if (symmetric.count_index_at(symmetric.index(it1,0),it1)>1 || symmetric.count_index_at(symmetric.index(it1,1),it1)>1) {
            prefactor=0;
            del_all_indices();
        }
        else if (symmetric.len()-it1>=2) {
            std::pair<size_t,size_t> itf1(it1+1,0), itf2(it1+1,0);
            bool part_evaluated(false);
            int n(0);
            while (!part_evaluated) {
                if (symmetric.count_index(symmetric.index(it1,1))>1 && symmetric.count_index(symmetric.index(it1,2))>1) {
                    itf1=symmetric.find_index(symmetric.index(it1,1),itf1.first);
                    itf2=symmetric.find_index(symmetric.index(it1,2),itf2.first);                    
                    // replace 2 d's by kronecker
                    if (itf1.first==itf2.first) {
                        if (itf1.second!=1) symmetric.swap_indices_at(itf1.first,itf1.second,1);
                        if (itf2.second!=2) symmetric.swap_indices_at(itf2.first,itf2.second,2);
                        prefactor*=(pow(NC,2)-4.0)/NC;
                        kronecker.set_indices(symmetric.index(it1,0),symmetric.index(itf1.first,0));
                        symmetric.del_indices(itf1.first);
                        symmetric.del_indices(it1);
                        part_evaluated=true;
                    }
                    // replace 3 d's by one d
                    else if (symmetric.matching_indices(itf1.first,itf2.first)>=1) {
                        if (itf1.second!=1) symmetric.swap_indices_at(itf1.first,itf1.second,1);
                        if (itf2.second!=2) symmetric.swap_indices_at(itf2.first,itf2.second,2);
                        while (symmetric.count_index_at(symmetric.index(itf1.first,2),itf2.first)==0) symmetric.swap_indices_at(itf1.first,0,2);
                        while (symmetric.index(itf1.first,2)!=symmetric.index(itf2.first,1)) symmetric.swap_indices_at(itf2.first,0,1);
                        symmetric.set_indices(symmetric.index(it1,0),symmetric.index(itf1.first,0),symmetric.index(itf2.first,0));
                        prefactor*=-(pow(NC,2)-12.0)/(2.*NC);
                        if (itf2.first>itf1.first) {
                            symmetric.del_indices(itf2.first);
                            symmetric.del_indices(itf1.first);
                        }
                        else {
                            symmetric.del_indices(itf1.first);
                            symmetric.del_indices(itf2.first);
                        }
                        symmetric.del_indices(it1);
                        part_evaluated=true;
                    }
                    // for replacement of 4 d's insert another else if, find a 4th d and check indices
                    // can also implement replacement of 2 d's with only one contracted index if sums are implemented
                    else it1++;
                }
                else if (n<3) {
                    symmetric.rotate_indices_at(it1);
                    n++;
                }
                else {
                    it1++;
                    part_evaluated=true;
                }
            }
        }
        else it1++;
        if (it1>=symmetric.len()) evaluated=true;
    }
}
void replace_2fd() {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated && antisymmetric.len()>0 && symmetric.len()>0) {
        std::pair<size_t,size_t> itf1(0,0), itf2(0,0);
        bool part_evaluated(false);
        int n(0);
        while (!part_evaluated && antisymmetric.len()>0) {
            if (antisymmetric.count_index(symmetric.index(it1,1))>0 && antisymmetric.count_index(symmetric.index(it1,2))>0) {
                itf1=antisymmetric.find_index(symmetric.index(it1,1),itf1.first);
                itf2=antisymmetric.find_index(symmetric.index(it1,2),itf2.first);
                // product of d and f with 2 contracted indices yields 0
                if (itf1.first==itf2.first) {
                    prefactor=0;
                    del_all_indices();
                }
                // replace 1 d and 2 f's by one d
                else if (antisymmetric.matching_indices(itf1.first,itf2.first)>=1) {
                    if (itf1.second!=1) {
                        antisymmetric.swap_indices_at(itf1.first,itf1.second,1);
                        prefactor*=-1;
                    }
                    if (itf2.second!=2) {
                        antisymmetric.swap_indices_at(itf2.first,itf2.second,2);
                        prefactor*=-1;
                    }
                    while (antisymmetric.count_index_at(antisymmetric.index(itf1.first,2),itf2.first)==0) {
                        antisymmetric.swap_indices_at(itf1.first,0,2);
                        prefactor*=-1;
                    }
                    while (antisymmetric.index(itf1.first,2)!=antisymmetric.index(itf2.first,1)) {
                        antisymmetric.swap_indices_at(itf2.first,0,1);
                        prefactor*=-1;
                    }
                    symmetric.set_indices(symmetric.index(it1,0),antisymmetric.index(itf1.first,0),antisymmetric.index(itf2.first,0));
                    prefactor*=-NC/2.;
                    if (itf2.first>itf1.first) {
                        antisymmetric.del_indices(itf2.first);
                        antisymmetric.del_indices(itf1.first);
                    }
                    else {
                        antisymmetric.del_indices(itf1.first);
                        antisymmetric.del_indices(itf2.first);
                    }
                    symmetric.del_indices(it1);
                    part_evaluated=true;
                }
                // for replacement of 2 f's and 2 d's insert another else if, find a 2nd d and check indices
                else it1++;
            }
            else if (n<3) {
                symmetric.rotate_indices_at(it1);
                n++;
            }
            else {
                it1++;
                part_evaluated=true;
            }
        }
        if (it1>=symmetric.len()) evaluated=true;
    }
}
void replace_2df() {
    size_t it1(0);
    bool evaluated(false);
    while(!evaluated && antisymmetric.len()>0 && symmetric.len()>0) {
        std::pair<size_t,size_t> itf1(0,0), itf2(0,0);
        bool part_evaluated(false);
        int n(0);
        while (!part_evaluated && symmetric.len()>0) {
            if (symmetric.count_index(antisymmetric.index(it1,1))>0 && symmetric.count_index(antisymmetric.index(it1,2))>0) {
                itf1=symmetric.find_index(antisymmetric.index(it1,1),itf1.first);
                itf2=symmetric.find_index(antisymmetric.index(it1,2),itf2.first);
                // product of d and f with 2 contracted indices yields 0
                if (itf1.first==itf2.first) {
                    prefactor=0;
                    del_all_indices();
                }
                // replace 1 d and 2 f's by one d
                else if (symmetric.matching_indices(itf1.first,itf2.first)>=1) {
                    if (itf1.second!=1) symmetric.swap_indices_at(itf1.first,itf1.second,1);
                    if (itf2.second!=2) symmetric.swap_indices_at(itf2.first,itf2.second,2);
                    while (symmetric.count_index_at(symmetric.index(itf1.first,2),itf2.first)==0) symmetric.swap_indices_at(itf1.first,0,2);
                    while (symmetric.index(itf1.first,2)!=symmetric.index(itf2.first,1)) symmetric.swap_indices_at(itf2.first,0,1);
                    antisymmetric.set_indices(antisymmetric.index(it1,0),symmetric.index(itf1.first,0),symmetric.index(itf2.first,0));
                    prefactor*=(pow(NC,2)-4.)/(2.*NC);
                    if (itf2.first>itf1.first) {
                        symmetric.del_indices(itf2.first);
                        symmetric.del_indices(itf1.first);
                    }
                    else {
                        symmetric.del_indices(itf1.first);
                        symmetric.del_indices(itf2.first);
                    }
                    antisymmetric.del_indices(it1);
                    part_evaluated=true;
                }
                // for replacement of 2 f's and 2 d's insert another else if, find a 2nd f and check indices
                else it1++;
            }
            else if (n<3) {
                antisymmetric.rotate_indices_at(it1);
                n++;
            }
            else {
                it1++;
                part_evaluated=true;
            }
        }
        if (it1>=antisymmetric.len()) evaluated=true;
    }
}
void del_all_indices() {
    kronecker.clear_indices();
    symmetric.clear_indices();
    antisymmetric.clear_indices();
}
bool is_free_index(int index) {
    int cntr(0);
    cntr+=kronecker.count_index(index);
    cntr+=symmetric.count_index(index);
    cntr+=antisymmetric.count_index(index);
    if (cntr==1) return true;
    else return false;
}
void print_expr() {
    string eval;
    eval+=to_string(prefactor);
    if (symmetric.len() > 0) {
        for (int i(0); i<symmetric.len(); i++) eval+="*d_["+to_string(symmetric.index(i,0))+","+to_string(symmetric.index(i,1))+","+to_string(symmetric.index(i,2))+"]";
    }
    if (antisymmetric.len() > 0) {
        for (int i(0); i<antisymmetric.len(); i++) eval+="*f_["+to_string(antisymmetric.index(i,0))+","+to_string(antisymmetric.index(i,1))+","+to_string(antisymmetric.index(i,2))+"]"; 
    }
    if (kronecker_eval.len() > 0) {
        for (int i(0); i<kronecker_eval.len(); i++) eval+="*k_["+to_string(kronecker_eval.index(i,0))+","+to_string(kronecker_eval.index(i,1))+"]";
    }
    cout << " = " << eval << endl;
}
