#include "colourtools.h"

// member functions of class process
process::process(void) {
    
}
process::~process(void) {
    
}
void process::add_in_leg(string ptcl) {
    size_t index=m_in_legs.size()+1;
    if (m_out_legs.size()>0) {
        for (size_t i(0);i<m_out_legs.size();i++) {
            m_out_legs.at(i).first+=1;
        }
    }
    m_in_legs.push_back(pair<size_t,string>(index,ptcl));
}
void process::add_out_leg(string ptcl) {
    size_t index=m_in_legs.size()+m_out_legs.size()+1;
    m_out_legs.push_back(pair<size_t,string>(index,ptcl));
}
void process::delete_all_legs() {
    m_in_legs.clear();
    m_out_legs.clear();
}
size_t process::no_of_legs() {
    return m_in_legs.size()+m_out_legs.size();
}
pair<size_t,string> process::leg(size_t index) {
    // leg numbering starts at 1 !
    if (index<=m_in_legs.size()) return m_in_legs.at(index-1);
    else if (index<=m_out_legs.size()+m_in_legs.size()) return m_out_legs.at(index-m_in_legs.size()-1);
    else {
        cerr << "Leg " << index << " does not exist in diagram." << endl;
        return pair<int,string>(0,"");
    }
}
bool process::is_in_leg(size_t index) {
    if (index<=m_in_legs.size()) return true;
    else return false;
}


// member functions of class delta
delta::delta(size_t i, size_t j, bool adj) {
    m_i=i;
    m_j=j;
    m_adj=adj;
}
delta::~delta() {
    
}
bool delta::is_free(size_t ind) {
    if (ind != m_i and ind != m_j) return true;
    else return false;
}
string delta::build_string() {
    string str;
    if (m_adj) str+="K_[";
    else str+="k_[";
    str+=to_string(m_i)+","+to_string(m_j)+"]";
    return str;
}

// member functions of class fundamental
fundamental::fundamental(size_t a, size_t i, size_t j) {
    m_a=a;
    m_i=i;
    m_j=j;
}
fundamental::~fundamental() {
    
}
bool fundamental::is_free(size_t ind) {
    if (ind != m_a and ind != m_i and ind != m_j) return true;
    else return false;
}
string fundamental::build_string() {
    string str("t_[");
    str+=to_string(m_a)+","+to_string(m_i)+","+to_string(m_j)+"]";
    return str;
}

// member functions of class antisymmetric
antisymmetric::antisymmetric(size_t a, size_t b, size_t c) {
    m_a=a;
    m_b=b;
    m_c=c;
}
antisymmetric::~antisymmetric() {
    
}
bool antisymmetric::is_free(size_t ind) {
    if (ind != m_a and ind != m_b and ind != m_c) return true;
    else return false;
}
string antisymmetric::build_string() {
    string str("f_[");
    str+=to_string(m_a)+","+to_string(m_b)+","+to_string(m_c)+"]";
    return str;
}

// member functions of class symmetric
symmetric::symmetric(size_t a, size_t b, size_t c) {
    m_a=a;
    m_b=b;
    m_c=c;
}
symmetric::~symmetric() {
    
}
bool symmetric::is_free(size_t ind) {
    if (ind != m_a and ind != m_b and ind != m_c) return true;
    else return false;
}
string symmetric::build_string() {
    string str("d_[");
    str+=to_string(m_a)+","+to_string(m_b)+","+to_string(m_c)+"]";
    return str;
}

// member functions of class c_term
c_term::c_term() {
    m_cnum=1.;
    m_NC_order=0;
    m_fi=1001;
}
c_term::c_term(delta k, fundamental t, antisymmetric f, symmetric d, complex<double> c, int NC) {
    m_cnum=c;
    m_NC_order=NC;
    m_k_vec.push_back(k);
    m_t_vec.push_back(t);
    m_f_vec.push_back(f);
    m_d_vec.push_back(d);
    m_fi=1001;
    while (!k.is_free(m_fi) or !t.is_free(m_fi) or !f.is_free(m_fi) or !d.is_free(m_fi)) m_fi++;
}
c_term::~c_term() {
    
}
void c_term::push_back(c_term ct) {
    m_cnum*=ct.m_cnum;
    m_NC_order+=ct.m_NC_order;
    for (auto& k : ct.m_k_vec) m_k_vec.push_back(k);
    for (auto& t : ct.m_t_vec) m_t_vec.push_back(t);
    for (auto& f : ct.m_f_vec) m_f_vec.push_back(f);
    for (auto& d : ct.m_d_vec) m_d_vec.push_back(d);
}
void c_term::push_back(delta k) {
    m_k_vec.push_back(k);
    while (!k.is_free(m_fi)) m_fi++;
}
void c_term::push_back(fundamental t) {
    m_t_vec.push_back(t);
    while (!t.is_free(m_fi)) m_fi++;
}
void c_term::push_back(antisymmetric f) {
    m_f_vec.push_back(f);
    while (!f.is_free(m_fi)) m_fi++;
}
void c_term::push_back(symmetric d) {
    m_d_vec.push_back(d);
    while (!d.is_free(m_fi)) m_fi++;
}
void c_term::cnumber(complex<double> c) {
    m_cnum=c;
}
void c_term::NC_order(int NCo) {
    m_NC_order=NCo;
}
void c_term::evaluate_deltas() {
    for (vector<delta>::iterator k_it(m_k_vec.begin()); k_it!=m_k_vec.end(); ++k_it) {
        bool ev(false);
        
        // replace adjoint indices
        if (k_it->m_adj==true) {
            if (k_it->m_i == k_it->m_j) {
                m_cnum*=(NC*NC-1);
                ev=true;
            }
            
            // replace indices in deltas
            for (vector<delta>::iterator k_it2(k_it+1);!ev and k_it2!=m_k_vec.end() ; ++k_it2) {
                if (k_it->m_i == k_it2->m_i) {
                    if (k_it->m_j == k_it2->m_j) {
                        m_cnum*=(NC*NC-1);
                        m_k_vec.erase(k_it2);
                        ev=true;
                    }
                    else {
                        k_it2->m_i=k_it->m_j;
                        ev=true;
                    }
                }
                else if (k_it->m_j == k_it2->m_i) {
                    if (k_it->m_i == k_it2->m_j) {
                        m_cnum*=(NC*NC-1);
                        m_k_vec.erase(k_it2);
                        ev=true;
                    }
                    else {
                        k_it2->m_i=k_it->m_i;
                        ev=true;
                    }
                }
                else if (k_it->m_i == k_it2->m_j) {
                    if (k_it->m_j == k_it2->m_i) {
                        m_cnum*=(NC*NC-1);
                        m_k_vec.erase(k_it2);
                        ev=true;
                    }
                    else {
                        k_it2->m_j=k_it->m_j;
                        ev=true;
                    }
                }
                else if (k_it->m_j == k_it2->m_j) {
                    if (k_it->m_i == k_it2->m_i) {
                        m_cnum*=(NC*NC-1);
                        m_k_vec.erase(k_it2);
                        ev=true;
                    }
                    else {
                        k_it2->m_j=k_it->m_i;
                        ev=true;
                    }
                }
            }
            // replace indices in fundamentals
            for (vector<fundamental>::iterator t_it(m_t_vec.begin()); !ev and t_it!=m_t_vec.end(); t_it++) {
                if (k_it->m_i == t_it->m_a) {
                    t_it->m_a=k_it->m_j;
                    ev=true;
                }
                else if (k_it->m_j == t_it->m_a) {
                    t_it->m_a=k_it->m_i;
                    ev=true;
                }
            }
            // replace indices in antisymmetrics
            for (vector<antisymmetric>::iterator f_it(m_f_vec.begin()); !ev and f_it!=m_f_vec.end(); f_it++) {
                if (k_it->m_i == f_it->m_a) {
                    f_it->m_a=k_it->m_j;
                    ev=true;
                }
                else if (k_it->m_j == f_it->m_a) {
                    f_it->m_a=k_it->m_i;
                    ev=true;
                }
                else if (k_it->m_i == f_it->m_b) {
                    f_it->m_b=k_it->m_j;
                    ev=true;
                }
                else if (k_it->m_j == f_it->m_b) {
                    f_it->m_b=k_it->m_i;
                    ev=true;
                }
                else if (k_it->m_i == f_it->m_c) {
                    f_it->m_c=k_it->m_j;
                    ev=true;
                }
                else if (k_it->m_j == f_it->m_c) {
                    f_it->m_c=k_it->m_i;
                    ev=true;
                }
            }
            // replace indices in symmetrics
            for (vector<symmetric>::iterator d_it(m_d_vec.begin()); !ev and d_it!=m_d_vec.end(); d_it++) {
                if (k_it->m_i == d_it->m_a) {
                    d_it->m_a=k_it->m_j;
                    ev=true;
                }
                else if (k_it->m_j == d_it->m_a) {
                    d_it->m_a=k_it->m_i;
                    ev=true;
                }
                else if (k_it->m_i == d_it->m_b) {
                    d_it->m_b=k_it->m_j;
                    ev=true;
                }
                else if (k_it->m_j == d_it->m_b) {
                    d_it->m_b=k_it->m_i;
                    ev=true;
                }
                else if (k_it->m_i == d_it->m_c) {
                    d_it->m_c=k_it->m_j;
                    ev=true;
                }
                else if (k_it->m_j == d_it->m_c) {
                    d_it->m_c=k_it->m_i;
                    ev=true;
                }
            }
            
            if (ev) {
                m_k_vec.erase(k_it);
                k_it--;
            }
        }
        
        // replace fundamental indices
        else {
            if (k_it->m_i == k_it->m_j) {
                m_cnum*=NC;
                ev=true;
            }
            
            // replace indices in deltas
            for (vector<delta>::iterator k_it2(k_it+1); !ev and k_it2!=m_k_vec.end(); ++k_it2) {
                if (k_it->m_i == k_it2->m_i) {
                    if (k_it->m_j == k_it2->m_j) {
                        m_cnum*=NC;
                        m_k_vec.erase(k_it2);
                        ev=true;
                    }
                    else {
                        k_it2->m_i=k_it->m_j;
                        ev=true;
                    }
                }
                else if (k_it->m_j == k_it2->m_i) {
                    if (k_it->m_i == k_it2->m_j) {
                        m_cnum*=NC;
                        m_k_vec.erase(k_it2);
                        ev=true;
                    }
                    else {
                        k_it2->m_i=k_it->m_i;
                        ev=true;
                    }
                }
                else if (k_it->m_i == k_it2->m_j) {
                    if (k_it->m_j == k_it2->m_i) {
                        m_cnum*=NC;
                        m_k_vec.erase(k_it2);
                        ev=true;
                    }
                    else {
                        k_it2->m_j=k_it->m_j;
                        ev=true;
                    }
                }
                else if (k_it->m_j == k_it2->m_j) {
                    if (k_it->m_i == k_it2->m_i) {
                        m_cnum*=NC;
                        m_k_vec.erase(k_it2);
                        ev=true;
                    }
                    else {
                        k_it2->m_j=k_it->m_i;
                        ev=true;
                    }
                }
            }
            
            // replace indices in fundamentals
            for (vector<fundamental>::iterator t_it(m_t_vec.begin()); !ev and t_it!=m_t_vec.end(); t_it++) {
                if (k_it->m_i == t_it->m_i) {
                    t_it->m_i=k_it->m_j;
                    ev=true;
                }
                else if (k_it->m_i == t_it->m_j) {
                    t_it->m_j=k_it->m_j;
                    ev=true;
                }
                else if (k_it->m_j == t_it->m_i) {
                    t_it->m_i=k_it->m_i;
                    ev=true;
                }
                else if (k_it->m_j == t_it->m_j) {
                    t_it->m_j=k_it->m_i;
                    ev=true;
                }
            }
            
            if (ev) {
                m_k_vec.erase(k_it);
                k_it--;
            }
        }
    }
}
void c_term::shift_inds(size_t by, bool all) {
    for (vector<delta>::iterator k_it(m_k_vec.begin()); k_it!=m_k_vec.end(); k_it++) {
        if (k_it->m_i>100 or all) k_it->m_i+=by;
        if (k_it->m_j>100 or all) k_it->m_j+=by;
    }
    for (vector<fundamental>::iterator t_it(m_t_vec.begin()); t_it!=m_t_vec.end(); t_it++) {
        if (t_it->m_a>100 or all) t_it->m_a+=by;
        if (t_it->m_i>100 or all) t_it->m_i+=by;
        if (t_it->m_j>100 or all) t_it->m_j+=by;
    }
    for (vector<antisymmetric>::iterator f_it(m_f_vec.begin()); f_it!=m_f_vec.end(); f_it++) {
        if (f_it->m_a>100 or all) f_it->m_a+=by;
        if (f_it->m_b>100 or all) f_it->m_b+=by;
        if (f_it->m_c>100 or all) f_it->m_c+=by;
    }
    for (vector<symmetric>::iterator d_it(m_d_vec.begin()); d_it!=m_d_vec.end(); d_it++) {
        if (d_it->m_a>100 or all) d_it->m_a+=by;
        if (d_it->m_b>100 or all) d_it->m_b+=by;
        if (d_it->m_c>100 or all) d_it->m_c+=by;
    }
    m_fi+=by;
}
c_term c_term::hconj() {
    c_term ct(*this);
    for (auto & t : ct.m_t_vec) {
        swap(t.m_i,t.m_j);
    }
    ct.m_cnum=conj(ct.m_cnum);
    return ct;
}
c_term c_term::operator*(c_term ct) {
    ct.m_cnum*=m_cnum;
    ct.m_NC_order+=m_NC_order;
    
    ct.shift_inds(10000,false);
    for (auto& k : m_k_vec) ct.push_back(k);
    for (auto& t : m_t_vec) ct.push_back(t);
    for (auto& f : m_f_vec) ct.push_back(f);
    for (auto& d : m_d_vec) ct.push_back(d);
    
    return ct;
}
complex<double> c_term::result() {
    if (m_k_vec.size()==0 and m_t_vec.size()==0 and m_f_vec.size()==0 and m_d_vec.size()==0) return m_cnum;
    else return complex<double>(NAN,NAN);
}
void c_term::clear() {
    m_k_vec.clear();
    m_t_vec.clear();
    m_f_vec.clear();
    m_d_vec.clear();
    m_NC_order=0;
    m_cnum=0.;
    m_fi=1001;
}
string c_term::build_string() {
    string str="";
    if (m_cnum!=1.) str+="c_["+to_string(m_cnum.real())+","+to_string(m_cnum.imag())+"]";
    for (auto& k : m_k_vec) {
        if (str!="") str+="*";
        str+=k.build_string();
    }
    for (auto& t : m_t_vec) {
        if (str!="") str+="*";
        str+=t.build_string();
    }
    for (auto& f : m_f_vec) {
        if (str!="") str+="*";
        str+=f.build_string();
    }
    for (auto& d : m_d_vec) {
        if (str!="") str+="*";
        str+=d.build_string();
    }
    return str;
}
void c_term::print() {
    cout<<build_string()<<endl;
}


// member functions of class c_amplitude
c_amplitude::c_amplitude() {
    m_result=0.;
}
c_amplitude::c_amplitude(c_term ct) {
    m_result=0.;
    m_cterm_vec.push_back(ct);
}
c_amplitude::c_amplitude(string expr) {
    m_result=0.;
    
    for (size_t i(0), mpos(expr.find('+'));mpos!=string::npos or expr.length()>0;mpos=expr.find('+')) {
        ++i;
        
        c_term ct;
        
        // decompose expr into summands
        string summand;
        if (mpos==string::npos) {
            summand=expr;
            expr="";
        }
        else {
            summand=expr.substr(0,mpos);
            expr=expr.substr(mpos+1);
        }
        // decompose summand into factors
        for (size_t j(0), mpos(summand.find('*'));mpos!=string::npos or summand.length()>0;mpos=summand.find('*')) {
            ++j;
            string factor;
            if (mpos==string::npos) {
                factor=summand;
                summand="";
            }
            else {
                factor=summand.substr(0,mpos);
                summand=summand.substr(mpos+1);
            }
            if (factor.find("c_[")==0 and factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==string::npos || factor.find(',',cpos+1)!=string::npos) {
                    cerr << "Invalid prefactor." << endl;
                    exit(EXIT_FAILURE);
                }
                complex<double> prfct(0.);
                int NC_order_c=0;
                int NC_order_r=0;
                size_t spos(factor.substr(3,cpos-3).find("1/NC"));
                if (spos==string::npos) prfct+=stod(factor.substr(3,cpos-3));
                else {
                    size_t ppos=factor.substr(3,cpos-3).find('^');
                    double exponent(1);
                    if (ppos!=string::npos) exponent=stod(factor.substr(3+ppos+1,cpos-ppos-2));
                    NC_order_r+=exponent;
                    prfct+=1./pow(NC,exponent);
                }
                spos=factor.substr(cpos+1,factor.length()-cpos-2).find("1/NC");
                if (spos==string::npos) prfct+=stod(factor.substr(cpos+1,factor.length()-cpos-2))*complex<double>(0.,1.);
                else {
                    size_t ppos=factor.substr(cpos+1,factor.length()-cpos-2).find('^');
                    double exponent(1);
                    if (ppos!=string::npos) exponent=stod(factor.substr(cpos+ppos+2,factor.length()-ppos-3));
                    NC_order_c+=exponent;
                    prfct+=(double)1./pow(NC,exponent)*complex<double>(0.,1.);
                }
                // only same orders in 1/NC in real and imaginary part are supported
                if (NC_order_c!=NC_order_r and NC_order_c!=0 and NC_order_r!=0) {
                    cerr << "Error: Expected equal orders in 1/NC for real and imaginary part of prefactor, but real part has order " << NC_order_r << " and imaginary part has order " << NC_order_c << ".\n Specify real and imaginary part seperately to avoid this error." << endl;
                    exit(EXIT_FAILURE);
                }
                ct.m_NC_order+=NC_order_r;
                ct.m_cnum*=prfct;
            }
            else if(factor.find("f_[")==0 and factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==string::npos) {
                    cerr << "Invalid number of indices for f." << endl;
                    exit(EXIT_FAILURE);
                }
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==string::npos || factor.find(',',c2pos+1)!=string::npos) {
                    cerr << "Invalid number of indices for f." << endl;
                    exit(EXIT_FAILURE);
                }
                ct.push_back(antisymmetric(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2))));
            }
            else if(factor.find("d_[")==0 and factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==string::npos) {
                    cerr << "Invalid number of indices for d." << endl;
                    exit(EXIT_FAILURE);
                }
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==string::npos || factor.find(',',c2pos+1)!=string::npos) {
                    cerr << "Invalid number of indices for d." << endl;
                    exit(EXIT_FAILURE);
                }
                ct.push_back(symmetric(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2))));
            }
            else if(factor.find("t_[")==0 and factor[factor.length()-1]==']') {
                size_t c1pos(factor.find(','));
                if (c1pos==string::npos) {
                    cerr << "Invalid number of indices for t." << endl;
                    exit(EXIT_FAILURE);
                }
                size_t c2pos(factor.find(',',c1pos+1));
                if (c2pos==string::npos || factor.find(',',c2pos+1)!=string::npos) {
                    cerr << "Invalid number of indices for t." << endl;
                    exit(EXIT_FAILURE);
                }
                ct.push_back(fundamental(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2))));
            }
            else if (factor.find("k_[")==0 and factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==string::npos || factor.find(',',cpos+1)!=string::npos) {
                    cerr << "Invalid number of indices for k." << endl;
                    exit(EXIT_FAILURE);
                }
                int ind_i=stoi(factor.substr(3,cpos-3)), ind_j=stoi(factor.substr(cpos+1,factor.length()-cpos-2));
                ct.push_back(delta(ind_i,ind_j,false));
            }
            else if (factor.find("K_[")==0 and factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==string::npos || factor.find(',',cpos+1)!=string::npos) {
                    cerr << "Invalid number of indices for K." << endl;
                    exit(EXIT_FAILURE);
                }
                int ind_i=stoi(factor.substr(3,cpos-3)), ind_j=stoi(factor.substr(cpos+1,factor.length()-cpos-2));
                ct.push_back(delta(ind_i,ind_j,true));
            }
            else cerr << "Invalid input." << endl;
        }
        m_cterm_vec.push_back(ct);
    }
}
c_amplitude::~c_amplitude() {
    
}
void c_amplitude::add(c_term ct) {
    m_cterm_vec.push_back(ct);
}
void c_amplitude::push_back(c_amplitude ca) {
    if (m_cterm_vec.size()==0) {
        for (vector<c_term>::iterator c_it(ca.m_cterm_vec.begin()); c_it!=ca.m_cterm_vec.end(); c_it++) {
            m_cterm_vec.push_back(*c_it);
        }
    }
    else {
        vector<c_term> new_c_terms;
        
        for (vector<c_term>::iterator c_it1(m_cterm_vec.begin()); c_it1!=m_cterm_vec.end(); c_it1++) {
            for (vector<c_term>::iterator c_it2(ca.m_cterm_vec.begin()); c_it2!=ca.m_cterm_vec.end(); c_it2++) {
                c_term ct(*c_it1);
                ct.push_back(*c_it2);
                new_c_terms.push_back(ct);
            }
        }
        
        m_cterm_vec=new_c_terms;
    }
}
c_amplitude c_amplitude::hconj() {
    c_amplitude ca(*this);
    
    for (vector<c_term>::iterator c_it(ca.m_cterm_vec.begin()); c_it!=ca.m_cterm_vec.end(); c_it++) (*c_it)=c_it->hconj();
    
    return ca;
}
c_amplitude c_amplitude::shift_to_internal(size_t by) {
    c_amplitude new_ca;
    for (auto ct : m_cterm_vec) {
        ct.shift_inds(by,true);
        new_ca.add(ct);
    }
    return new_ca;
}
c_amplitude c_amplitude::operator*(complex<double> z) {
    c_amplitude new_ca;
    for (vector<c_term>::iterator c_it(m_cterm_vec.begin()); c_it!=m_cterm_vec.end(); c_it++) {
        c_it->m_cnum*=z;
        new_ca.add(*c_it);
    }
    return new_ca;
}
c_amplitude c_amplitude::operator*(c_amplitude ca){
    if (m_cterm_vec.size()==0) return (*this);
    if (ca.m_cterm_vec.size()==0) return ca;
    
    c_amplitude new_ca;
    
    for (vector<c_term>::iterator c_it1(m_cterm_vec.begin()); c_it1!=m_cterm_vec.end(); c_it1++)
        for (vector<c_term>::iterator c_it2(ca.m_cterm_vec.begin()); c_it2!=ca.m_cterm_vec.end(); c_it2++)
            new_ca.add((*c_it1)*(*c_it2));
    
    return new_ca;
}
complex<double> c_amplitude::scprod(c_amplitude ca) {
    c_amplitude scp(ca.hconj());
    scp=(*this)*scp;
    
//    cout<<"\n---------------------------"<<endl;
//    scp.print();
    scp.evaluate();
//    cout<<" = ";
//    scp.print();
//    cout<<" = "<<scp.result()<<endl;
//    cout<<"---------------------------\n"<<endl;
    
    return scp.result();
}
void c_amplitude::clear() {
    m_cterm_vec.clear();
}
void c_amplitude::evaluate() {
    while (m_cterm_vec.size()>0) {
//        vector<c_term> new_c_terms;
        c_amplitude ca;
        vector<c_term>::iterator c_it(m_cterm_vec.begin());
        
//        cout<<"\n-------------------------"<<endl;
//        c_it->print();
        
    startpoint:
        
        // evaluate deltas
        c_it->evaluate_deltas();
        
//        cout<<"= "<<endl;
//        c_it->print();
        
        // evaluate antisymmetrics
        for (vector<antisymmetric>::iterator f_it(c_it->m_f_vec.begin()); c_it->m_cnum!=0. and f_it!=c_it->m_f_vec.end();f_it++) {
            bool ev(false);
            
            // check if term vanishes
            if (f_it->m_a == f_it->m_b or f_it->m_a == f_it->m_c or f_it->m_b == f_it->m_c) {
                c_it->clear();
                ev=true;
            }
            
            // contract with indices in symmetrics
            for (vector<symmetric>::iterator d_it(c_it->m_d_vec.begin()); !ev and d_it!=c_it->m_d_vec.end() ; d_it++) {
                if (f_it->m_a == d_it->m_b) swap<size_t>(d_it->m_a,d_it->m_b);
                else if (f_it->m_a == d_it->m_c) swap<size_t>(d_it->m_a,d_it->m_c);
                else if (f_it->m_b == d_it->m_a) {
                    swap<size_t>(f_it->m_b,f_it->m_a);
                    swap<size_t>(f_it->m_b,f_it->m_c);
                }
                else if (f_it->m_b == d_it->m_b) {
                    swap<size_t>(f_it->m_b,f_it->m_a);
                    swap<size_t>(f_it->m_b,f_it->m_c);
                    swap<size_t>(d_it->m_b,d_it->m_a);
                }
                else if (f_it->m_b == d_it->m_c) {
                    swap<size_t>(f_it->m_b,f_it->m_a);
                    swap<size_t>(f_it->m_b,f_it->m_c);
                    swap<size_t>(d_it->m_c,d_it->m_a);
                }
                else if (f_it->m_c == d_it->m_a) {
                    swap<size_t>(f_it->m_b,f_it->m_a);
                    swap<size_t>(f_it->m_b,f_it->m_c);
                }
                else if (f_it->m_c == d_it->m_b) {
                    swap<size_t>(f_it->m_c,f_it->m_a);
                    swap<size_t>(f_it->m_b,f_it->m_c);
                    swap<size_t>(d_it->m_b,d_it->m_a);
                }
                else if (f_it->m_c == d_it->m_c) {
                    swap<size_t>(f_it->m_c,f_it->m_a);
                    swap<size_t>(f_it->m_b,f_it->m_c);
                    swap<size_t>(d_it->m_c,d_it->m_a);
                }
                
                if (f_it->m_a == d_it->m_a) {
                    if (f_it->m_b == d_it->m_b or f_it->m_b == d_it->m_c) {
                        c_it->clear();
                        c_it->m_cnum=0.;
                        ev=true;
                    }
                    else {
                        // replace f_{abc}d_{ade}
                        
                        size_t b(f_it->m_b), c(f_it->m_c), d(d_it->m_b), e(d_it->m_c);
                        c_it->m_d_vec.erase(d_it);
                        c_it->m_f_vec.erase(f_it);
                        c_term ct(*c_it);
                        f_it--;
                        ev=true;
                        
                        // first term
                        c_it->m_cnum*=complex<double>(0.,-1./TR);
                        c_it->push_back(fundamental(b,c_it->m_fi,++c_it->m_fi));
                        c_it->push_back(fundamental(c,c_it->m_fi-1,c_it->m_fi));
                        c_it->push_back(fundamental(d,c_it->m_fi-1,c_it->m_fi));
                        c_it->push_back(fundamental(e,c_it->m_fi-1,c_it->m_fi-4));
                        
                        
                        // second term
                        c_term ct1(ct);
                        ct1.m_cnum*=complex<double>(0.,1./TR);
                        ct1.push_back(fundamental(c,ct1.m_fi,++ct1.m_fi));
                        ct1.push_back(fundamental(b,ct1.m_fi-1,ct1.m_fi));
                        ct1.push_back(fundamental(d,ct1.m_fi-1,ct1.m_fi));
                        ct1.push_back(fundamental(e,ct1.m_fi-1,ct1.m_fi-4));
                        //                            new_c_terms.push_back(ct1);
                        ca.add(ct1);
                        
                        // second term
                        ct1=ct;
                        ct1.m_cnum*=complex<double>(0.,-1./TR);
                        ct1.push_back(fundamental(d,ct1.m_fi,++ct1.m_fi));
                        ct1.push_back(fundamental(b,ct1.m_fi-1,ct1.m_fi));
                        ct1.push_back(fundamental(c,ct1.m_fi-1,ct1.m_fi));
                        ct1.push_back(fundamental(e,ct1.m_fi-1,ct1.m_fi-4));
                        //                            new_c_terms.push_back(ct1);
                        ca.add(ct1);
                        
                        // third term
                        ct.m_cnum*=complex<double>(0.,1./TR);
                        ct.push_back(fundamental(d,ct.m_fi,++ct.m_fi));
                        ct.push_back(fundamental(c,ct.m_fi-1,ct.m_fi));
                        ct.push_back(fundamental(b,ct.m_fi-1,ct.m_fi));
                        ct.push_back(fundamental(e,ct.m_fi-1,ct.m_fi-4));
                        //                            new_c_terms.push_back(ct);
                        ca.add(ct);
                    }
                }
            }
            
            // contract with indices in antisymmetrics
            for (vector<antisymmetric>::iterator f_it2(c_it->m_f_vec.begin()); !ev and f_it2!=c_it->m_f_vec.end();f_it2++) {
                if (f_it2!=f_it) {
                    if (f_it->m_a == f_it2->m_b) {
                        swap<size_t>(f_it2->m_a,f_it2->m_b);
                        swap<size_t>(f_it2->m_b,f_it2->m_c);
                    }
                    else if (f_it->m_a == f_it2->m_c) {
                        swap<size_t>(f_it2->m_a,f_it2->m_c);
                        swap<size_t>(f_it2->m_c,f_it2->m_b);
                    }
                    else if (f_it->m_b == f_it2->m_a) {
                        swap<size_t>(f_it->m_b,f_it->m_a);
                        swap<size_t>(f_it->m_c,f_it->m_b);
                    }
                    else if (f_it->m_b == f_it2->m_b) {
                        swap<size_t>(f_it->m_b,f_it->m_a);
                        swap<size_t>(f_it2->m_b,f_it2->m_a);
                    }
                    else if (f_it->m_b == f_it2->m_c) {
                        swap<size_t>(f_it->m_b,f_it->m_a);
                        swap<size_t>(f_it2->m_c,f_it2->m_a);
                    }
                    else if (f_it->m_c == f_it2->m_a) {
                        swap<size_t>(f_it->m_c,f_it->m_a);
                        swap<size_t>(f_it->m_c,f_it->m_b);
                    }
                    else if (f_it->m_c == f_it2->m_b) {
                        swap<size_t>(f_it->m_c,f_it->m_a);
                        swap<size_t>(f_it2->m_b,f_it2->m_a);
                    }
                    else if (f_it->m_c == f_it2->m_c) {
                        swap<size_t>(f_it->m_c,f_it->m_a);
                        swap<size_t>(f_it2->m_c,f_it2->m_a);
                    }
                    
                
                    if (f_it->m_a == f_it2->m_a) {
                        if (f_it->m_b == f_it2->m_b) {
                            if (f_it->m_c == f_it2->m_c) c_it->m_cnum*=NC*(NC*NC-1);
                            else {
                                c_it->m_cnum*=NC;
                                c_it->m_k_vec.push_back(delta(f_it->m_c,f_it2->m_c,true));
                            }
                            c_it->m_f_vec.erase(f_it2);
                            c_it->m_f_vec.erase(f_it);
                            f_it--;
                            ev=true;
                        }
                        else if (f_it->m_b == f_it2->m_c) {
                            if (f_it->m_c == f_it2->m_b) c_it->m_cnum*=-1.*NC*(NC*NC-1);
                            else {
                                c_it->m_cnum*=-1.*NC;
                                c_it->m_k_vec.push_back(delta(f_it->m_c,f_it2->m_b,true));
                            }
                            c_it->m_f_vec.erase(f_it2);
                            c_it->m_f_vec.erase(f_it);
                            f_it--;
                            ev=true;
                        }
                        else if (f_it->m_c == f_it2->m_c) {
                            if (f_it->m_b == f_it2->m_b) c_it->m_cnum*=1.*NC*(NC*NC-1);
                            else {
                                c_it->m_cnum*=1.*NC;
                                c_it->m_k_vec.push_back(delta(f_it->m_b,f_it2->m_b,true));
                            }
                            c_it->m_f_vec.erase(f_it2);
                            c_it->m_f_vec.erase(f_it);
                            f_it--;
                            ev=true;
                        }
                        else if (f_it->m_c == f_it2->m_b) {
                            if (f_it->m_b == f_it2->m_c) c_it->m_cnum*=-1.*NC*(NC*NC-1);
                            else {
                                c_it->m_cnum*=-1.*NC;
                                c_it->m_k_vec.push_back(delta(f_it->m_b,f_it2->m_c,true));
                            }
                            c_it->m_f_vec.erase(f_it2);
                            c_it->m_f_vec.erase(f_it);
                            f_it--;
                            ev=true;
                        }
                        else {
                            // replace f_{abc}f_{ade}

                            size_t b(f_it->m_b), c(f_it->m_c), d(f_it2->m_b), e(f_it2->m_c);
                            c_it->m_f_vec.erase(f_it2);
                            c_it->m_f_vec.erase(f_it);
                            c_term ct(*c_it);
                            f_it--;
                            ev=true;

                            // first term
                            c_it->m_cnum*=-1./TR;
                            c_it->push_back(fundamental(d,c_it->m_fi,++c_it->m_fi));
                            c_it->push_back(fundamental(e,c_it->m_fi-1,c_it->m_fi));
                            c_it->push_back(fundamental(b,c_it->m_fi-1,c_it->m_fi));
                            c_it->push_back(fundamental(c,c_it->m_fi-1,c_it->m_fi-4));


                            // second term
                            c_term ct1(ct);
                            ct1.m_cnum*=1./TR;
                            ct1.push_back(fundamental(b,ct1.m_fi,++ct1.m_fi));
                            ct1.push_back(fundamental(d,ct1.m_fi-1,ct1.m_fi));
                            ct1.push_back(fundamental(e,ct1.m_fi-1,ct1.m_fi));
                            ct1.push_back(fundamental(c,ct1.m_fi-1,ct1.m_fi-4));
//                            new_c_terms.push_back(ct1);
                            ca.add(ct1);

                            // second term
                            ct1=ct;
                            ct1.m_cnum*=1./TR;
                            ct1.push_back(fundamental(e,ct1.m_fi,++ct1.m_fi));
                            ct1.push_back(fundamental(d,ct1.m_fi-1,ct1.m_fi));
                            ct1.push_back(fundamental(b,ct1.m_fi-1,ct1.m_fi));
                            ct1.push_back(fundamental(c,ct1.m_fi-1,ct1.m_fi-4));
//                            new_c_terms.push_back(ct1);
                            ca.add(ct1);

                            // third term
                            ct.m_cnum*=-1./TR;
                            ct.push_back(fundamental(b,ct.m_fi,++ct.m_fi));
                            ct.push_back(fundamental(e,ct.m_fi-1,ct.m_fi));
                            ct.push_back(fundamental(d,ct.m_fi-1,ct.m_fi));
                            ct.push_back(fundamental(c,ct.m_fi-1,ct.m_fi-4));
//                            new_c_terms.push_back(ct);
                            ca.add(ct);
                        }
                    }
                }
            }
            
            // contract with indices in fundamentals
            for (vector<fundamental>::iterator t_it(c_it->m_t_vec.begin()); !ev and t_it!=c_it->m_t_vec.end();t_it++) {
                if (f_it->m_b == t_it->m_a) {
                    swap<size_t>(f_it->m_a,f_it->m_b);
                    swap<size_t>(f_it->m_b,f_it->m_c);
                }
                else if (f_it->m_c == t_it->m_a) {
                    swap<size_t>(f_it->m_a,f_it->m_c);
                    swap<size_t>(f_it->m_c,f_it->m_b);
                }
                
                if (f_it->m_a == t_it->m_a) {
                    size_t b(f_it->m_b), c(f_it->m_c), i(t_it->m_i), j(t_it->m_j);
                    c_it->m_f_vec.erase(f_it);
                    c_it->m_t_vec.erase(t_it);
                    c_term ct(*c_it);
                    f_it--;
                    ev=true;
                    
                    // first term
                    c_it->m_cnum*=complex<double>(0.,-1.);
                    c_it->push_back(fundamental(b,i,c_it->m_fi));
                    c_it->push_back(fundamental(c,c_it->m_fi-1,j));
                    
                    // second term
                    ct.m_cnum*=complex<double>(0.,1.);
                    ct.push_back(fundamental(c,i,ct.m_fi));
                    ct.push_back(fundamental(b,ct.m_fi-1,j));
//                    new_c_terms.push_back(ct);
                    ca.add(ct);
                }
            }
            
//            if (ev) c_it->evaluate_deltas();
            if (ev) goto startpoint;
        }
        
//        cout<<"= "<<endl;
//        c_it->print();
    
        // evaluate symmetrics
        for (vector<symmetric>::iterator d_it(c_it->m_d_vec.begin()); d_it!=c_it->m_d_vec.end();d_it++) {
            bool ev(false);
            
            // replace indices in symmetrics
            for (vector<symmetric>::iterator d_it2(d_it+1); !ev and d_it2!=c_it->m_d_vec.end() ; d_it2++) {
                if (d_it->m_a == d_it2->m_b) swap<size_t>(d_it2->m_a,d_it2->m_b);
                else if (d_it->m_a == d_it2->m_c) swap<size_t>(d_it2->m_a,d_it2->m_c);
                
                if (d_it->m_a == d_it2->m_a) {
                    if (d_it->m_b == d_it2->m_b) {
                        if (d_it->m_c == d_it2->m_c) c_it->m_cnum*=(NC*NC-1.)*(NC*NC-4.)/NC;
                        else {
                            c_it->m_cnum*=(NC*NC-4.)/NC;
                            c_it->m_k_vec.push_back(delta(d_it->m_c,d_it2->m_c,true));
                        }
                        ev=true;
                    }
                    else if (d_it->m_b == d_it2->m_c) {
                        if (d_it->m_c and d_it2->m_b) c_it->m_cnum*=(NC*NC-1.)*(NC*NC-4.)/NC;
                        else {
                            c_it->m_cnum*=(NC*NC-4.)/NC;
                            c_it->m_k_vec.push_back(delta(d_it->m_c,d_it2->m_b,true));
                        }
                        ev=true;
                    }
                    if (ev) {
                        c_it->m_d_vec.erase(d_it2);
                        c_it->m_d_vec.erase(d_it);
                        d_it--;
                    }
                }
            }
            
            if (!ev) {
                // replace d_{abc}
//                cout<<"replace "<<d_it->build_string()<<endl;
                size_t a(d_it->m_a), b(d_it->m_b), c(d_it->m_c);
                c_it->m_d_vec.erase(d_it);
                c_term ct(*c_it);
                d_it--;
                ev=true;

                // first term
                c_it->m_cnum*=1./TR;
                c_it->push_back(fundamental(a,c_it->m_fi,++c_it->m_fi));
                c_it->push_back(fundamental(b,c_it->m_fi-1,c_it->m_fi));
                c_it->push_back(fundamental(c,c_it->m_fi-1,c_it->m_fi-3));
//                cout<<"->";
//                c_it->print();

                // second term
                ct.m_cnum*=1./TR;
                ct.push_back(fundamental(b,ct.m_fi,++ct.m_fi));
                ct.push_back(fundamental(a,ct.m_fi-1,ct.m_fi));
                ct.push_back(fundamental(c,ct.m_fi-1,ct.m_fi-3));
//                cout<<"->";
//                ct.print();
//                new_c_terms.push_back(ct);
                ca.add(ct);
            }
//            if (ev) c_it->evaluate_deltas();
            if (ev) goto startpoint;
        }
        
//        cout<<"= "<<endl;
//        c_it->print();
        
        // replace fundamentals
        for (vector<fundamental>::iterator t_it(c_it->m_t_vec.begin()); c_it->m_cnum!=0. and t_it!=c_it->m_t_vec.end(); t_it++) {
            bool ev(false);
            
            // check if fundamental vanishes
            if (t_it->m_i == t_it->m_j) {
                c_it->clear();
                c_it->m_cnum=0.;
                ev=true;
            }
            
            // contract indices with fundamentals
            for (vector<fundamental>::iterator t_it2(c_it->m_t_vec.begin()); !ev and t_it2!=c_it->m_t_vec.end() and c_it->m_t_vec.size()>0; t_it2++) {
                if (t_it2!=t_it) {
                    if (t_it->m_a == t_it2->m_a) {
                        if (t_it->m_j == t_it2->m_i) {
                            if (t_it->m_i == t_it2->m_j) {
                                c_it->m_cnum*=TR*(NC*NC-1);
                            }
                            else {
                                c_it->m_cnum*=CF;
                                c_it->m_k_vec.push_back(delta(t_it->m_i,t_it2->m_j,false));
                            }
                            c_it->m_t_vec.erase(t_it2);
                            c_it->m_t_vec.erase(t_it);
                            t_it--;
                        }
                        else {
                            // replace (T_a)_{ij}(T_a)_{kl}
                            size_t i(t_it->m_i), j(t_it->m_j), k(t_it2->m_i), l(t_it2->m_j);
                            c_it->m_t_vec.erase(t_it2);
                            c_it->m_t_vec.erase(t_it);
                            c_term ct(*c_it);
                            t_it--;

                            // first term
                            c_it->m_cnum*=TR;
                            c_it->m_k_vec.push_back(delta(i,l,false));
                            c_it->m_k_vec.push_back(delta(j,k,false));

                            // second term
                            ct.NC_order(c_it->m_NC_order+1);
                            ct.m_cnum*=-1./(2.*NC);
                            ct.m_k_vec.push_back(delta(i,j,false));
                            ct.m_k_vec.push_back(delta(k,l,false));
//                                new_c_terms.push_back(ct);
                            ca.add(ct);
                        }
                        ev=true;
                    }
                    else if ((t_it->m_i == t_it2->m_j and t_it->m_j == t_it2->m_i) or (t_it->m_j == t_it2->m_i and t_it->m_i == t_it2->m_j)) {
                        c_it->m_cnum*=TR;
                        c_it->m_k_vec.push_back(delta(t_it->m_a,t_it2->m_a,true));
                        c_it->m_t_vec.erase(t_it2);
                        c_it->m_t_vec.erase(t_it);
                        t_it--;
                        ev=true;
                    }
                }
            }
            
            c_it->evaluate_deltas();
        }
        
        c_it->evaluate_deltas();
        m_result+=c_it->result();
//        cout<<"= "<<endl;
//        c_it->print();
//        cout<<"-------------------------\n"<<endl;
        m_cterm_vec.erase(c_it);
        
        ca.evaluate();
        m_result+=ca.result();
        
//        for (const auto& n : new_c_terms) {
////            m_cterm_vec.push_back(n);
//            c_amplitude ca(n);
//            ca.evaluate();
//            m_result+=ca.result();
//        }
    }
}
void c_amplitude::evaluate(size_t up_to_NC) {
    while (m_cterm_vec.size()>0) {
        //        vector<c_term> new_c_terms;
        c_amplitude ca;
        vector<c_term>::iterator c_it(m_cterm_vec.begin());
        
        // evaluate deltas
        c_it->evaluate_deltas();
        
        // evaluate antisymmetrics
        for (vector<antisymmetric>::iterator f_it(c_it->m_f_vec.begin()); c_it->m_cnum!=0. and f_it!=c_it->m_f_vec.end();f_it++) {
            bool ev(false);
            
            // check if term vanishes
            if (f_it->m_a == f_it->m_b or f_it->m_a == f_it->m_c or f_it->m_b == f_it->m_c) {
                c_it->clear();
                ev=true;
            }
            
            // contract with indices in antisymmetrics
            for (vector<antisymmetric>::iterator f_it2(c_it->m_f_vec.begin()); !ev and f_it2!=c_it->m_f_vec.end();f_it2++) {
                if (f_it2!=f_it) {
                    if (f_it->m_a == f_it2->m_b) {
                        swap<size_t>(f_it2->m_a,f_it2->m_b);
                        swap<size_t>(f_it2->m_b,f_it2->m_c);
                    }
                    else if (f_it->m_a == f_it2->m_c) {
                        swap<size_t>(f_it2->m_a,f_it2->m_c);
                        swap<size_t>(f_it2->m_c,f_it2->m_b);
                    }
                    
                    if (f_it->m_a == f_it2->m_a) {
                            // replace f_{abc}f_{ade}
                            
                            size_t b(f_it->m_b), c(f_it->m_c), d(f_it2->m_b), e(f_it2->m_c);
                            c_it->m_f_vec.erase(f_it2);
                            c_it->m_f_vec.erase(f_it);
                            c_term ct(*c_it);
                            f_it--;
                            ev=true;
                            
                            // first term
                            c_it->m_cnum*=-1./TR;
                            c_it->push_back(fundamental(d,c_it->m_fi,++c_it->m_fi));
                            c_it->push_back(fundamental(e,c_it->m_fi-1,c_it->m_fi));
                            c_it->push_back(fundamental(b,c_it->m_fi-1,c_it->m_fi));
                            c_it->push_back(fundamental(c,c_it->m_fi-1,c_it->m_fi-4));
                            
                            
                            // second term
                            c_term ct1(ct);
                            ct1.m_cnum*=1./TR;
                            ct1.push_back(fundamental(b,ct1.m_fi,++ct1.m_fi));
                            ct1.push_back(fundamental(d,ct1.m_fi-1,ct1.m_fi));
                            ct1.push_back(fundamental(e,ct1.m_fi-1,ct1.m_fi));
                            ct1.push_back(fundamental(c,ct1.m_fi-1,ct1.m_fi-4));
                            //                            new_c_terms.push_back(ct1);
                            ca.add(ct1);
                            
                            // second term
                            ct1=ct;
                            ct1.m_cnum*=1./TR;
                            ct1.push_back(fundamental(e,ct1.m_fi,++ct1.m_fi));
                            ct1.push_back(fundamental(d,ct1.m_fi-1,ct1.m_fi));
                            ct1.push_back(fundamental(b,ct1.m_fi-1,ct1.m_fi));
                            ct1.push_back(fundamental(c,ct1.m_fi-1,ct1.m_fi-4));
                            //                            new_c_terms.push_back(ct1);
                            ca.add(ct1);
                            
                            // third term
                            ct.m_cnum*=-1./TR;
                            ct.push_back(fundamental(b,ct.m_fi,++ct.m_fi));
                            ct.push_back(fundamental(e,ct.m_fi-1,ct.m_fi));
                            ct.push_back(fundamental(d,ct.m_fi-1,ct.m_fi));
                            ct.push_back(fundamental(c,ct.m_fi-1,ct.m_fi-4));
                            //                            new_c_terms.push_back(ct);
                            ca.add(ct);
                    }
                }
            }
            
            // contract with indices in fundamentals
            for (vector<fundamental>::iterator t_it(c_it->m_t_vec.begin()); !ev and t_it!=c_it->m_t_vec.end();t_it++) {
                if (f_it->m_b == t_it->m_a) {
                    swap<size_t>(f_it->m_a,f_it->m_b);
                    swap<size_t>(f_it->m_b,f_it->m_c);
                }
                else if (f_it->m_c == t_it->m_a) {
                    swap<size_t>(f_it->m_a,f_it->m_c);
                    swap<size_t>(f_it->m_c,f_it->m_b);
                }
                
                if (f_it->m_a == t_it->m_a) {
                    size_t b(f_it->m_b), c(f_it->m_c), i(t_it->m_i), j(t_it->m_j);
                    c_it->m_f_vec.erase(f_it);
                    c_it->m_t_vec.erase(t_it);
                    c_term ct(*c_it);
                    f_it--;
                    ev=true;
                    
                    // first term
                    c_it->m_cnum*=complex<double>(0.,-1.);
                    c_it->push_back(fundamental(b,i,c_it->m_fi));
                    c_it->push_back(fundamental(c,c_it->m_fi-1,j));
                    
                    // second term
                    ct.m_cnum*=complex<double>(0.,1.);
                    ct.push_back(fundamental(c,i,ct.m_fi));
                    ct.push_back(fundamental(b,ct.m_fi-1,j));
                    //                    new_c_terms.push_back(ct);
                    ca.add(ct);
                }
            }
            
            if (ev) c_it->evaluate_deltas();
        }
        
        // evaluate symmetrics
        for (vector<symmetric>::iterator d_it(c_it->m_d_vec.begin()); d_it!=c_it->m_d_vec.end();d_it++) {
            // replace d_{abc}
            //                cout<<"replace "<<d_it->build_string()<<endl;
            size_t a(d_it->m_a), b(d_it->m_b), c(d_it->m_c);
            c_it->m_d_vec.erase(d_it);
            c_term ct(*c_it);
            d_it--;
            
            // first term
            c_it->m_cnum*=1./TR;
            c_it->push_back(fundamental(a,c_it->m_fi,++c_it->m_fi));
            c_it->push_back(fundamental(b,c_it->m_fi-1,c_it->m_fi));
            c_it->push_back(fundamental(c,c_it->m_fi-1,c_it->m_fi-3));
            //                cout<<"->";
            //                c_it->print();
            
            // second term
            ct.m_cnum*=1./TR;
            ct.push_back(fundamental(b,ct.m_fi,++ct.m_fi));
            ct.push_back(fundamental(a,ct.m_fi-1,ct.m_fi));
            ct.push_back(fundamental(c,ct.m_fi-1,ct.m_fi-3));
            //                cout<<"->";
            //                ct.print();
            //                new_c_terms.push_back(ct);
            ca.add(ct);
            
            c_it->evaluate_deltas();
        }
        
        
        // replace fundamentals
        for (vector<fundamental>::iterator t_it(c_it->m_t_vec.begin()); c_it->m_cnum!=0. and t_it!=c_it->m_t_vec.end(); t_it++) {
            bool ev(false);
            
            // check if fundamental vanishes
            if (t_it->m_i == t_it->m_j) {
                c_it->clear();
                c_it->m_cnum=0.;
                ev=true;
            }
            
            // contract indices with fundamentals
            for (vector<fundamental>::iterator t_it2(c_it->m_t_vec.begin()); !ev and t_it2!=c_it->m_t_vec.end() and c_it->m_t_vec.size()>0; t_it2++) {
                if (t_it2!=t_it) {
                    if (t_it->m_a == t_it2->m_a) {
                        // replace (T_a)_{ij}(T_a)_{kl}
                        size_t i(t_it->m_i), j(t_it->m_j), k(t_it2->m_i), l(t_it2->m_j);
                        c_it->m_t_vec.erase(t_it2);
                        c_it->m_t_vec.erase(t_it);
                        c_term ct(*c_it);
                        t_it--;
                        
                        // first term
                        c_it->m_cnum*=TR;
                        c_it->m_k_vec.push_back(delta(i,l,false));
                        c_it->m_k_vec.push_back(delta(j,k,false));
                        
                        // second term
                        if (c_it->m_NC_order<up_to_NC) {
                            ct.NC_order(c_it->m_NC_order+1);
                            ct.m_cnum*=-1./(2.*NC);
                            ct.m_k_vec.push_back(delta(i,j,false));
                            ct.m_k_vec.push_back(delta(k,l,false));
                            // new_c_terms.push_back(ct);
                            ca.add(ct);
                        }
                    
                        ev=true;
                    }
                }
            }
            
            c_it->evaluate_deltas();
        }
        
        c_it->evaluate_deltas();
        m_result+=c_it->result();
        m_cterm_vec.erase(c_it);
        
        ca.evaluate(up_to_NC);
        m_result+=ca.result();
    }
}
size_t c_amplitude::no_of_terms() {
    return m_cterm_vec.size();
}
string c_amplitude::build_string() {
    string str("");
    
    for (auto& ct : m_cterm_vec) {
        if (str!="") str+="+";
        str+=ct.build_string();
    }
    
    return str;
}
void c_amplitude::print() {
    cout<<this->build_string()<<endl;
}
complex<double> c_amplitude::result() {
    if (m_cterm_vec.size()==0) return m_result;
    else return complex<double>(NAN,NAN);
}



