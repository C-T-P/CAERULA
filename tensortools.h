#ifndef TENSORTOOLS_H
#define TENSORTOOLS_H

class three_ind {
    std::vector<std::vector<int>> ind;
    public:
        void set_indices(int i, int j, int k) {
            std::vector<int> new_ind {i, j, k};
            ind.push_back(new_ind);
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
struct terms {
    std::vector<three_ind> sym;
    std::vector<three_ind> asym;
    std::vector<three_ind> fund;
    std::vector<two_ind> kron;
    std::vector<std::complex<float>> pref;
    size_t no_of_terms () {
        return pref.size();
    };
};

#endif
