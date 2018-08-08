#include "trace_basis.h"

vector<vector<size_t>> get_q_ind_combinations(vector<size_t> q_inds, vector<size_t> qb_inds);

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
trace_t trace_t::conj() {
    if (m_qb==0 and m_q==0 and m_g.size()>2) {
        size_t n_g(m_g.size());
        vector<size_t> new_m_g(n_g,m_g.at(0));
        reverse_copy(m_g.begin()+1,m_g.end(),new_m_g.begin()+1);
        return trace_t(new_m_g);
    }
    return trace_t(0,0);
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
bool trace_t::comp(trace_t& tr_t) {
    vector<size_t> indices1((*this).get_indices()), indices2(tr_t.get_indices());
    size_t s_1(indices1.size()), s_2(indices2.size());
    
    if (s_1>s_2) return true;
    else if (s_1==s_2) {
        size_t i(0);
        while (i<s_1 and indices1.at(i)==indices2.at(i)) i++;
        if (i<s_1 and indices1.at(i)>indices2.at(i)) return true;
        return false;
    }
    else return false;
}
bool trace_t::operator==(trace_t& rhs) {
    vector<size_t> indices1((*this).get_indices()), indices2(rhs.get_indices());
    size_t s_1(indices1.size()), s_2(indices2.size());
    if (s_1!=s_2) return false;
    size_t i(0);
    while (i<s_1 and indices1.at(i)==indices2.at(i)) i++;
    if (i==s_1) return true;
    return false;
}
colour_term trace_t::build_ct() {
    colour_term ct;
    three_ind symmetric;
    three_ind antisymmetric;
    three_ind fundamental;
    two_ind kronecker;
    size_t n_g(m_g.size());
    
    if (n_g==0) {
        kronecker.set_indices(m_qb, m_q, false);
        ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
        return ct;
    }
    if (m_qb==0 and m_q==0) {
        if (n_g==2) {
            kronecker.set_indices(m_g.at(0), m_g.at(1), true);
            ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
            return ct;
        }
        
        size_t start_ind(101), incr(0);
        for (size_t i(0);i<n_g;i++) {
            int c_ind;
            if (i==n_g-1) c_ind=start_ind;
            else c_ind=start_ind+incr+1;
            fundamental.set_indices(m_g.at(i),start_ind+incr,c_ind);
            incr++;
        }
        ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
        
        fundamental.clear_indices();
        vector<size_t> refl_ind(n_g,m_g.at(0));
        reverse_copy(m_g.begin()+1,m_g.end(),refl_ind.begin()+1);
        
        start_ind=101, incr=0;
        for (size_t i(0);i<n_g;i++) {
            int c_ind;
            if (i==n_g-1) c_ind=start_ind;
            else c_ind=start_ind+incr+1;
            fundamental.set_indices(refl_ind.at(i),start_ind+incr,c_ind);
            incr++;
        }
        double pref(1.);
        if (n_g%2!=0) pref=-1;
        ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(pref,0.),0);
        
        return ct;
    }
    if (n_g==1) {
        fundamental.set_indices(m_g.at(0),m_qb,m_q);
        ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
        return ct;
    }
    size_t start_ind(101), incr(0);
    fundamental.set_indices(m_g.at(0),m_qb,start_ind);
    while (incr<n_g-2) {
        fundamental.set_indices(m_g.at(incr+1),start_ind+incr,start_ind+incr+1);
        incr++;
    }
    fundamental.set_indices(m_g.back(),start_ind+incr,m_q);
    ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
    return ct;
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
vector<trace_vec> trace_vec::conjugates() {
    vector<trace_vec> conjugate_tr_vecs({*this});
    
    for (size_t i(0);i<m_tr_vec.size();i++) {
        trace_t refl(m_tr_vec.at(i).conj());
        if (refl.is_non_zero()) {
            size_t curr_s(conjugate_tr_vecs.size());
            for (size_t j(0);j<curr_s;j++) {
                trace_vec tmp_tr_vec(conjugate_tr_vecs.at(j));
                tmp_tr_vec.at(i)=refl;
                conjugate_tr_vecs.push_back(tmp_tr_vec);
            }
        }
    }
    
    conjugate_tr_vecs.erase(conjugate_tr_vecs.begin());
    return conjugate_tr_vecs;
}
vector<size_t> trace_vec::get_indices() {
    vector<size_t> indices;
    
    for (auto& tr_t : m_tr_vec) {
        vector<size_t> tmp(tr_t.get_indices());
        indices.insert(indices.end(),tmp.begin(),tmp.end());
    }
    
    return indices;
}
size_t trace_vec::no_groups() {
    return m_tr_vec.size();
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
    
    size_t s_1((*this).no_groups()), s_2(tr_v.no_groups());
    if (s_1>s_2) return !true;
    if (s_1<s_2) return !false;
    
    size_t i(0);
    while (i<s_1 and (*this).at(i)==tr_v.at(i)) i++;
    if (i==s_1) return false;
    return !(*this).at(i).comp(tr_v.at(i));
}
bool trace_vec::operator==(trace_vec& rhs) {
    size_t i(0), s_1(m_tr_vec.size()), s_2(rhs.no_groups());
    if (s_1!=s_2) return false;
    while (i<s_1 and (*this).at(i)==rhs.at(i)) i++;
    if (i==s_1) return true;
    return false;
}
colour_term trace_vec::build_ct() {
    colour_term ct;
    for (auto& tr_t : m_tr_vec) {
        if (ct.no_of_terms()==0) ct=tr_t.build_ct();
        else ct=ct.multiply(tr_t.build_ct());
    }
    return ct;
}
void trace_vec::print() {
    cout<<"[";
    for (auto& tr : m_tr_vec) tr.print();
    cout<<"]"<<endl;
}



// member functions of trace_basis class
trace_basis::trace_basis(size_t n_g, size_t n_qp) {
    m_ng=n_g;
    m_nqp=n_qp;
    
    // initialise quark, anti quark and gluon indices 
    // the process is ordered as q,...,q,qb,...,qb,g,...,g
    for (size_t n(1);n<=m_nqp;n++) m_q_indices.push_back(n);
    for (size_t n(m_nqp+1);n<=2*m_nqp;n++) m_qb_indices.push_back(n);
    for (size_t n(2*m_nqp+1);n<=2*m_nqp+m_ng;n++) m_g_indices.push_back(n);
    size_t g_start(0);
    
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
    if (m_nqp!=0) {
        vector<vector<size_t>> qqb_ind_combos(get_q_ind_combinations(m_q_indices,m_qb_indices));
        
        for (const auto& qqb_inds : qqb_ind_combos) {
            trace_vec tmp_trv;
            for (size_t qp_no(0);qp_no<2*m_nqp;qp_no+=2)
                tmp_trv.push_back(trace_t(qqb_inds.at(qp_no+1),qqb_inds.at(qp_no)));
            m_tr_basis.push_back(tmp_trv);
        }
    }
    // if there are no quarks
    else {
        if (m_ng==1) {
            cerr<<"Error: no trace basis to build out of 0 quark pairs and 1 gluon."<<endl;
            exit(EXIT_FAILURE);
        }
        
        trace_vec tmp_trv(trace_t(vector<size_t>(m_g_indices.begin(),m_g_indices.begin()+1)));
        m_tr_basis.push_back(tmp_trv);
        g_start++;
    }
    
    
    // successively add all n_g gluons
    for (size_t i(g_start);i<m_ng;i++) {
        size_t g(m_g_indices.at(i));
        vector<trace_vec> trb_cpy;
        
        for (auto& bv : m_tr_basis) {
            vector<trace_vec> new_bvs(bv.add_one_gluon(g));
            for (auto& v : new_bvs) trb_cpy.push_back(v); 
        }
        m_tr_basis=trb_cpy;
    }
    
    (*this).remove_sg();
    (*this).normal_order();
    (*this).remove_conj();
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
void trace_basis::remove_conj() {
    for (size_t i(0);i<m_tr_basis.size();i++) {
        vector<trace_vec> conjugate_tvs(m_tr_basis.at(i).conjugates());
        
        for (auto& con : conjugate_tvs) {
            for (size_t j(0);j<m_tr_basis.size();j++) {
                if (m_tr_basis.at(j)==con) {
                    m_tr_basis.erase(m_tr_basis.begin()+j);
                    break;
                }
            }
        }
    }
}
void trace_basis::normal_order() {
    sort(m_tr_basis.begin(), m_tr_basis.end(), [ ]( trace_vec& lhs, trace_vec& rhs )
    {
        return lhs.comp(rhs);
    });
}
size_t trace_basis::dim() {
    return m_tr_basis.size();
}
process trace_basis::proc() {
    process proc;
    for (size_t n(0);n<m_nqp;n++) proc.add_out_leg("q");
    for (size_t n(0);n<m_nqp;n++) proc.add_out_leg("qb");
    for (size_t n(0);n<m_ng;n++) proc.add_out_leg("g");
    return proc;
}
vector<vector<size_t>> trace_basis::perms() {
    vector<vector<size_t>> perms;
    
    for (auto& v : m_tr_basis)
        if (v.is_connected()) perms.push_back(v.get_indices());
    
    return perms;
}
vector<colour_term> trace_basis::ct_basis() {
    vector<colour_term> ct_basis;
    
    for (auto& bv : m_tr_basis)
        ct_basis.push_back(bv.build_ct());
    
    return ct_basis;
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
















