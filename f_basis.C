#include "f_basis.h"

// member functions of f_type class
f_type::f_type(vector<size_t> g_inds) {
    m_g=g_inds;
}
f_type::~f_type(void) {

}
vector<f_type> f_type::add_one_gluon(size_t g_ind) {
    size_t n_g(m_g.size());
    
    vector<f_type> new_fs;
    if (n_g==1) {
        new_fs.push_back(f_type(vector<size_t>({m_g.at(0),g_ind})));
    }
    else {
        for (size_t i(1);i<n_g;i++) {
            vector<size_t> new_g_inds(m_g);
            new_g_inds.insert(new_g_inds.begin()+i,g_ind);
            new_fs.push_back(f_type(new_g_inds));
        }
    }
    return new_fs;
}
vector<size_t> f_type::get_indices() {
    return m_g;
}
size_t f_type::no_g() {
    return m_g.size();
}
bool f_type::is_non_zero() {
    if (m_g.size()>0) return true;
    else return false;
}
bool f_type::vanishes() {
    if (m_g.size()<=1) return true;
    else return false;
}
bool f_type::comp(f_type& tr_t) {
    vector<size_t> indices1((*this).get_indices()), indices2(tr_t.get_indices());
    size_t s_1(indices1.size()), s_2(indices2.size());
    
    if (s_1>s_2) return true;
    else if (s_1==s_2) {
        size_t i(0);
        while (i<s_1 and indices1.at(i)==indices2.at(i)) i++;
        if (i<s_1 and indices1.at(i)>indices2.at(i)) return false;
        return true;
    }
    else return false;
}
bool f_type::operator==(f_type rhs) {
    vector<size_t> indices1((*this).get_indices()), indices2(rhs.get_indices());
    size_t s_1(indices1.size()), s_2(indices2.size());
    if (s_1!=s_2) return false;
    size_t i(0);
    while (i<s_1 and indices1.at(i)==indices2.at(i)) i++;
    if (i==s_1) return true;
    return false;
}
colour_term f_type::build_ct() {
    colour_term ct;
    three_ind symmetric;
    three_ind antisymmetric;
    three_ind fundamental;
    two_ind kronecker;
    size_t n_g(m_g.size());
    
    if (!(*this).vanishes()) {
        if (n_g==2) {
            kronecker.set_indices(m_g.at(0), m_g.at(1), true);
            ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
            return ct;
        }
        if (n_g==3) {
            antisymmetric.set_indices(m_g.at(0),m_g.at(1),m_g.at(2));
            ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
            return ct;
        }
        size_t start_ind(101), incr(0);
        antisymmetric.set_indices(m_g.at(0),m_g.at(1),start_ind);
        for (size_t i(2);i<n_g-1;i++) {
            int c_ind;
            if (i==n_g-2) c_ind=m_g.at(i+1);
            else c_ind=start_ind+incr+1;
            antisymmetric.set_indices(start_ind+incr,m_g.at(i),c_ind);
            incr++;
        }
        ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
        return ct;
    }
    else {
        cerr<<"Error: Cannot build colour_term from empty f_type."<<endl;
        exit(EXIT_FAILURE);
    }
}
void f_type::print() {
    cout<<"(";
    for (vector<size_t>::iterator g_it(m_g.begin());g_it!=m_g.end();g_it++) {
        if (g_it!=m_g.begin()) cout<<",";
        cout<<*g_it;
    }
    cout<<")";
}


// member functions of f_vec class
f_vec::f_vec(f_type f) {
    if (f.is_non_zero())
        m_f_vec.push_back(f);
}
f_vec::~f_vec(void) {

}
void f_vec::push_back(f_type f) {
    if (f.is_non_zero())
        m_f_vec.push_back(f);
    else {
        cerr<<"Error: can't push empty f_type to f_vec."<<endl;
        exit(EXIT_FAILURE);
    }
}
f_type& f_vec::at(size_t i) {
    return m_f_vec.at(i);
}
f_type f_vec::at(size_t i) const {
    return m_f_vec.at(i);
}
vector<f_vec> f_vec::add_one_gluon(size_t g_ind) {
    size_t no_of_f_vecs(m_f_vec.size());

    vector<f_vec> new_f_vecs;
    
    for (size_t i(0);i<no_of_f_vecs;i++) {
        f_type f(m_f_vec.at(i));
        vector<f_type> new_f_ts(f.add_one_gluon(g_ind));
        for (auto& f_t : new_f_ts) {
            f_vec tmp_f_vec(*this);
            tmp_f_vec.at(i)=f_t;
            new_f_vecs.push_back(tmp_f_vec);
        }
    }

    // open new f_type
    f_vec cpy(*this);
    cpy.push_back(f_type({g_ind}));
    new_f_vecs.push_back(cpy);

    return new_f_vecs;
}
vector<size_t> f_vec::get_indices() {
    vector<size_t> indices;

    for (auto& f_t : m_f_vec) {
        vector<size_t> tmp(f_t.get_indices());
        indices.insert(indices.end(),tmp.begin(),tmp.end());
    }

    return indices;
}
size_t f_vec::no_groups() {
    return m_f_vec.size();
}
bool f_vec::is_connected() {
    if (m_f_vec.size()==1) return true;
    else return false;
}
bool f_vec::has_sg() {
    for (auto& f_t : m_f_vec)
        if (f_t.vanishes()) return true;
    return false;
}
void f_vec::order() {
    sort(m_f_vec.begin(), m_f_vec.end(), [ ]( f_type& lhs, f_type& rhs )
    {
        return lhs.comp(rhs);
    });
}
bool f_vec::comp(f_vec& f_v) {
    (*this).order();
    f_v.order();

    size_t s_1((*this).no_groups()), s_2(f_v.no_groups());
    if (s_1>s_2) return !true;
    if (s_1<s_2) return !false;

    size_t i(0);
    while (i<s_1 and (*this).at(i)==f_v.at(i)) i++;
    if (i==s_1) return false;
    return !(*this).at(i).comp(f_v.at(i));
}
colour_term f_vec::build_ct() {
    colour_term ct;
    for (auto& f_t : m_f_vec) {
        if (ct.no_of_terms()==0) ct=f_t.build_ct();
        else ct=ct.multiply(f_t.build_ct());
    }
    return ct;
}
void f_vec::print() {
    cout<<"[";
    for (auto& f : m_f_vec) f.print();
    cout<<"]"<<endl;
}


//// member functions of f_basis class
f_basis::f_basis(size_t n_g) {
    m_ng=n_g;

    // initialise gluon indices
    for (size_t n(1);n<=m_ng;n++) m_g_indices.push_back(n);

    if (m_ng<3) {
        cerr<<"Error: the adjoint basis needs at least 3 gluons."<<endl;
        exit(EXIT_FAILURE);
    }
    
    m_f_basis.push_back(f_vec(f_type(vector<size_t>({m_g_indices.at(0)}))));

    // successively add all n_g gluons
    for (size_t i(1);i<m_ng;i++) {
        size_t g(m_g_indices.at(i));
        vector<f_vec> fb_cpy;

        for (auto& bv : m_f_basis) {
            vector<f_vec> new_bvs(bv.add_one_gluon(g));
            for (auto& v : new_bvs) fb_cpy.push_back(v);
        }
        m_f_basis=fb_cpy;
    }

    (*this).remove_sg();
    (*this).normal_order();
}
f_basis::~f_basis() {

}
void f_basis::remove_sg() {
    for (size_t i(0);i<m_f_basis.size();i++) {
        if (m_f_basis.at(i).has_sg()) {
            m_f_basis.erase(m_f_basis.begin()+i);
            i--;
        }
    }
}
void f_basis::normal_order() {
    sort(m_f_basis.begin(), m_f_basis.end(), [ ]( f_vec& lhs, f_vec& rhs )
    {
        return lhs.comp(rhs);
    });
}
size_t f_basis::dim() {
    return m_f_basis.size();
}
process f_basis::proc() {
    process proc;
    for (size_t n(0);n<m_ng;n++) proc.add_out_leg("g");
    return proc;
}
vector<vector<size_t>> f_basis::perms() {
    vector<vector<size_t>> perms;

    for (auto& v : m_f_basis)
        if (v.is_connected()) perms.push_back(v.get_indices());

    return perms;
}
vector<colour_term> f_basis::ct_basis() {
    vector<colour_term> ct_basis;

    for (auto& bv : m_f_basis)
        ct_basis.push_back(bv.build_ct());

    return ct_basis;
}
void f_basis::print() {
    for (auto& bv : m_f_basis) bv.print();
}
