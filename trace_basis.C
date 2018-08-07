#include "trace_basis.h"

vector<vector<size_t>> get_q_ind_combinations(vector<size_t> q_inds, vector<size_t> qb_inds);

int main(int argc, char **argv) {
    size_t n_g, n_qp;
    
    cout<<"Enter number of gluons:"<<endl;
    cin>>n_g;
    
    cout<<"Enter number of quark pairs:"<<endl;
    cin>>n_qp;
    
    trace_basis t_basis(n_g, n_qp);
    t_basis.print();
    
    cout<<"Dimension: "<<t_basis.dim()<<endl;
    
    cout<<"Permutations:"<<endl;
    vector<vector<size_t>> amp_perms(t_basis.get_perms());
    for (const auto& p : amp_perms) {
        for (const auto& i : p) cout<<i<<" ";
        cout<<endl;
    }
}


// member functions of trace_t class
trace_t::trace_t(vector<size_t> g_inds, size_t qb_ind, size_t q_ind) {
    m_q=q_ind;
    m_qb=qb_ind;
    m_g=g_inds;
}
trace_t::trace_t(size_t qb_ind, size_t q_ind) {
    m_q=q_ind;
    m_qb=qb_ind;
    m_g=vector<size_t>();
}
trace_t::~trace_t(void) {

}
vector<trace_t> trace_t::add_one_gluon(size_t g_ind) {
    size_t n_g(m_g.size()), pos(0);
    if (m_qb==0) pos++;
    
    vector<trace_t> new_trs;
    while (pos<=n_g) {
        vector<size_t> new_g_inds(m_g);
        new_g_inds.insert(new_g_inds.begin()+pos,g_ind);
        pos++;
        new_trs.push_back(trace_t(new_g_inds,m_qb,m_q));
    }
    return new_trs;
}
vector<size_t> trace_t::get_indices() {
    vector<size_t> indices;
    if (m_qb!=0) indices.push_back(m_qb);
    for (vector<size_t>::iterator g_it(m_g.begin());g_it!=m_g.end();g_it++) indices.push_back(*g_it);
    if (m_q!=0) indices.push_back(m_q);
    
    return indices;
}
size_t trace_t::no_g() {
    return m_g.size();
}
size_t trace_t::no_qp() {
    if (m_qb!=0 and m_q!=0) return 1;
    return 0;
}
bool trace_t::is_non_zero() {
    if ((m_qb!=0 and m_q!=0) or m_g.size()>0) return true;
    else return false;
}
bool trace_t::vanishes() {
    if (m_qb==0 and m_q==0 and m_g.size()<=1) return true;
    else return false;
}
void trace_t::print() {
    if (m_qb!=0) cout<<"{"<<m_qb;
    else cout<<"(";
    for (vector<size_t>::iterator g_it(m_g.begin());g_it!=m_g.end();g_it++) {
        if (g_it!=m_g.begin() or m_qb!=0) cout<<",";
        cout<<*g_it;
    }
    if (m_q!=0) cout<<","<<m_q<<"}";
    else cout<<")";
}
bool trace_t::comp(trace_t& tr_t) {
    vector<size_t> indices1((*this).get_indices()), indices2(tr_t.get_indices());
    size_t s_1(indices1.size()), s_2(indices2.size());
    
    if (s_1>s_2) return false;
    else if (s_1==s_2) {
        size_t i(0);
        while (i<s_1 and indices1.at(i)==indices2.at(i)) i++;
        if (i<s_1 and indices1.at(i)>indices2.at(i)) return true;
        return false;
    }
    else return true;
}



// member functions of trace_vec class
trace_vec::trace_vec(trace_t tr) {
    if (tr.is_non_zero())
        m_tr_vec.push_back(tr);
}
trace_vec::~trace_vec(void) {
    
}
void trace_vec::push_back(trace_t tr) {
    if (tr.is_non_zero()) 
        m_tr_vec.push_back(tr);
    else {
        cerr<<"Error: can't push empty trace_t to trace_vec."<<endl;
        exit(EXIT_FAILURE);
    }
}
trace_t& trace_vec::at(size_t i) {
    return m_tr_vec.at(i);
}
trace_t trace_vec::at(size_t i) const {
    return m_tr_vec.at(i);
}
vector<trace_vec> trace_vec::add_one_gluon(size_t g_ind) {
    size_t no_of_tr_vecs(m_tr_vec.size());
    
    vector<trace_vec> new_tr_vecs;
    for (size_t i(0);i<no_of_tr_vecs;i++) {
        trace_t tr(m_tr_vec.at(i));
        vector<trace_t> new_tr_ts(tr.add_one_gluon(g_ind));
        for (auto& tr_t : new_tr_ts) {
            trace_vec tmp_tr_vec(*this);
            tmp_tr_vec.at(i)=tr_t;
            new_tr_vecs.push_back(tmp_tr_vec);
        }
    }
    
    // open new trace_t
    trace_vec cpy(*this);
    cpy.push_back(trace_t({g_ind}));
    new_tr_vecs.push_back(cpy);
    
    return new_tr_vecs;
}
vector<size_t> trace_vec::get_indices() {
    vector<size_t> indices;
    
    for (auto& tr_t : m_tr_vec) {
        vector<size_t> tmp(tr_t.get_indices());
        indices.insert(indices.end(),tmp.begin(),tmp.end());
    }
    
    return indices;
}
bool trace_vec::is_connected() {
    size_t n_ql(0), n_con_g(0), n_qlg(0);
    for (auto & tr_t : m_tr_vec) {
        size_t ng(tr_t.no_g()), nqp(tr_t.no_qp());
        if (ng==0 and nqp!=0) n_ql++;
        else if (ng!=0 and nqp==0) n_con_g++;
        else n_qlg++;
    }
    if (n_ql>0 and n_con_g==0 and n_qlg==0) return true;
    if (n_ql==0 and n_con_g==1 and n_qlg==0) return true;
    if (n_ql>0 and n_con_g==0 and n_qlg==1) return true;
    return false;
}
bool trace_vec::has_sg() {
    for (auto& tr_t : m_tr_vec) 
        if (tr_t.vanishes()) return true;
    return false;
}
void trace_vec::order() {
    sort(m_tr_vec.begin(), m_tr_vec.end(), [ ]( trace_t& lhs, trace_t& rhs )
    {
        return !lhs.comp(rhs);
    });
}
bool trace_vec::comp(trace_vec& tr_v) {
    (*this).order();
    tr_v.order();
    
    return (*this).at(0).comp(tr_v.at(0));
}
void trace_vec::print() {
    cout<<"[";
    for (auto& tr : m_tr_vec) tr.print();
    cout<<"]"<<endl;
}



// member functions of trace_basis class
trace_basis::trace_basis(size_t n_g, size_t n_qp) {
    
    // initialise quark, anti quark and gluon indices 
    // the process is ordered as q,...,q,qb,...,qb,g,...,g
    vector<size_t> g_indices, q_indices, qb_indices;
    for (size_t n(1);n<=n_qp;n++) q_indices.push_back(n);
    for (size_t n(n_qp+1);n<=2*n_qp;n++) qb_indices.push_back(n);
    for (size_t n(2*n_qp+1);n<=2*n_qp+n_g;n++) g_indices.push_back(n);
    
//     cout<<"q indices:"<<endl;
//     for (const auto& i : q_indices) cout<<i<<" ";
//     cout<<endl;
//     cout<<"qb indices:"<<endl;
//     for (const auto& i : qb_indices) cout<<i<<" ";
//     cout<<endl;
//     cout<<"g indices:"<<endl;
//     for (const auto& i : g_indices) cout<<i<<" ";
//     cout<<endl;
    
    // build all possible quark lines
    if (n_qp!=0) {
        vector<vector<size_t>> qqb_ind_combos(get_q_ind_combinations(q_indices,qb_indices));
        
        for (const auto& qqb_inds : qqb_ind_combos) {
            trace_vec tmp_trv;
            for (size_t qp_no(0);qp_no<2*n_qp;qp_no+=2)
                tmp_trv.push_back(trace_t(qqb_inds.at(qp_no+1),qqb_inds.at(qp_no)));
            m_tr_basis.push_back(tmp_trv);
        }
    }
    // if there are no quarks
    else {
        if (n_g==1) {
            cerr<<"Error: no trace basis to build out of 0 quark pairs and 1 gluon."<<endl;
            exit(EXIT_FAILURE);
        }
        
        trace_vec tmp_trv(trace_t(vector<size_t>(g_indices.begin(),g_indices.begin()+1)));
        m_tr_basis.push_back(tmp_trv);
        g_indices=vector<size_t>(g_indices.begin()+1,g_indices.end());
    }
    
    
    // successively add all n_g gluons
    for (const auto& g : g_indices) {
        vector<trace_vec> trb_cpy;
        
        for (auto& bv : m_tr_basis) {
            vector<trace_vec> new_bvs(bv.add_one_gluon(g));
            for (auto& v : new_bvs) trb_cpy.push_back(v); 
        }
        m_tr_basis=trb_cpy;
    }
    
    (*this).remove_sg();
    (*this).normal_order();
}
trace_basis::~trace_basis() {
    
}
void trace_basis::remove_sg() {
    for (size_t i(0);i<m_tr_basis.size();i++) {
        if (m_tr_basis.at(i).has_sg()) {
            m_tr_basis.erase(m_tr_basis.begin()+i);
            i--;
        }
    }
}
size_t trace_basis::dim() {
    return m_tr_basis.size();
}
void trace_basis::normal_order() {
    sort(m_tr_basis.begin(), m_tr_basis.end(), [ ]( trace_vec& lhs, trace_vec& rhs )
    {
        return !lhs.comp(rhs);
    });
}
vector<vector<size_t>> trace_basis::get_perms() {
    vector<vector<size_t>> perms;
    
    for (auto& v : m_tr_basis)
        if (v.is_connected()) perms.push_back(v.get_indices());
    
    return perms;
}
void trace_basis::print() {
    for (auto& bv : m_tr_basis) bv.print();
}



// helper functions
vector<vector<size_t>> get_q_ind_combinations(vector<size_t> q_inds, vector<size_t> qb_inds) {
    vector<vector<size_t>> q_ind_perms, qqb_ind_combos;
    vector<size_t> tmp;
    q_ind_perms.push_back(q_inds);
    while (next_permutation(q_inds.begin(),q_inds.end()))
        q_ind_perms.push_back(q_inds);
    
    for (const auto& q_i : q_ind_perms) {
        tmp.clear();
        for (size_t i(0);i<q_inds.size();i++) {
            tmp.push_back(q_i[i]);
            tmp.push_back(qb_inds[i]);
        }
        qqb_ind_combos.push_back(tmp);
    }
    return qqb_ind_combos;
}



















