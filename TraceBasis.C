// Copyright (C) 2021 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#include "TraceBasis.h"

vector<vector<size_t>> get_q_ind_combinations(vector<size_t> q_inds, vector<size_t> qb_inds);

//*****************************************************************************
//
// Member Functions of TraceType
//
//*****************************************************************************

TraceType::TraceType(vector<size_t> g_inds, size_t q_ind, size_t qb_ind) {
  m_q=q_ind;
  m_qb=qb_ind;
  m_g=g_inds;
}
TraceType::TraceType(size_t q_ind, size_t qb_ind) {
  m_q=q_ind;
  m_qb=qb_ind;
  m_g=vector<size_t>();
}
TraceType::~TraceType(void) {

}
vector<TraceType> TraceType::add_one_gluon(size_t g_ind) {
  size_t n_g(m_g.size()), pos(0);
  if (m_qb==0) pos++;
  
  vector<TraceType> new_trs;
  while (pos<=n_g) {
    vector<size_t> new_g_inds(m_g);
    new_g_inds.insert(new_g_inds.begin()+pos,g_ind);
    pos++;
    new_trs.push_back(TraceType(new_g_inds,m_q,m_qb));
  }
  return new_trs;
}
vector<size_t> TraceType::get_indices() {
  vector<size_t> indices;
  if (m_qb!=0) indices.push_back(m_q);
  for (vector<size_t>::iterator g_it(m_g.begin());g_it!=m_g.end();g_it++) indices.push_back(*g_it);
  if (m_q!=0) indices.push_back(m_qb);
  
  return indices;
}
TraceType TraceType::conj() {
  if (m_qb==0 and m_q==0 and m_g.size()>2) {
    size_t n_g(m_g.size());
    vector<size_t> new_m_g(n_g,m_g.at(0));
    reverse_copy(m_g.begin()+1,m_g.end(),new_m_g.begin()+1);
    return TraceType(new_m_g);
  }
  return TraceType(0,0);
}
size_t TraceType::no_g() {
  return m_g.size();
}
size_t TraceType::no_qp() {
  if (m_qb!=0 and m_q!=0) return 1;
  return 0;
}
bool TraceType::is_not_empty() {
  if ((m_qb!=0 and m_q!=0) or m_g.size()>0) return true;
  else return false;
}
bool TraceType::vanishes() {
  if (m_qb==0 and m_q==0 and m_g.size()<=1) return true;
  else return false;
}
bool TraceType::operator>(TraceType& rhs) {
  vector<size_t> indices1(this->get_indices()), indices2(rhs.get_indices());
  size_t s_1(indices1.size()), s_2(indices2.size());
  
  if (s_1>s_2) return true;
  else if (s_1==s_2) {
    size_t i(0);
    while (i<s_1 and indices1.at(i)==indices2.at(i)) ++i;
    if (i<s_1 and indices1.at(i)<indices2.at(i)) return true;
    return false;
  }
  else return false;
}
bool TraceType::operator==(TraceType& rhs) {
  vector<size_t> indices1((*this).get_indices()), indices2(rhs.get_indices());
  size_t s_1(indices1.size()), s_2(indices2.size());
  if (s_1!=s_2) return false;
  size_t i(0);
  while (i<s_1 and indices1.at(i)==indices2.at(i)) ++i;
  if (i==s_1) return true;
  return false;
}
CAmplitude TraceType::build_ca(size_t start_ind, bool combineConj) {
  CTerm ct;
  size_t n_g(m_g.size());
  
  if (n_g==0) {
    ct.push_back(Delta(m_q, m_qb, false));
    return CAmplitude(ct);
  }
  if (m_qb==0 and m_q==0) {
    if (n_g==2) {
      ct.push_back(Delta(m_g.at(0), m_g.at(1), true));
      return CAmplitude(ct);
    }
    
    size_t incr(0);
    for (size_t i(0);i<n_g;i++) {
      int c_ind;
      if (i==n_g-1) c_ind=start_ind;
      else c_ind=start_ind+incr+1;
      ct.push_back(Fundamental(m_g.at(i),start_ind+incr,c_ind));
      incr++;
    }
    
    CAmplitude ca(ct);

    // Add conjugate amplitude
    if (combineConj) {
      ct.clear();
      
      vector<size_t> refl_ind(n_g,m_g.at(0));
      reverse_copy(m_g.begin()+1,m_g.end(),refl_ind.begin()+1);
    
      incr=0;
      for (size_t i(0);i<n_g;i++) {
	int c_ind;
	if (i==n_g-1) c_ind=start_ind;
	else c_ind=start_ind+incr+1;
	ct.push_back(Fundamental(refl_ind.at(i),start_ind+incr,c_ind));
	incr++;
      }
      if (n_g%2!=0) ct.set_cnumber(ColourFactor(-1., 0, 0, 0, 0));
      else ct.set_cnumber(ColourFactor(1., 0, 0, 0, 0));
      
      ca.add(ct);
    }
    
    return ca;
  }
  if (n_g==1) {
    ct.push_back(Fundamental(m_g.at(0),m_q,m_qb));
    return CAmplitude(ct);
  }
  
  size_t incr(0);
  ct.push_back(Fundamental(m_g.at(0),m_q,start_ind));
  while (incr<n_g-2) {
    ct.push_back(Fundamental(m_g.at(incr+1),start_ind+incr,start_ind+incr+1));
    incr++;
  }
  ct.push_back(Fundamental(m_g.back(),start_ind+incr,m_qb));
  return CAmplitude(ct);
}
void TraceType::print() {
  if (m_qb!=0) cout<<"{"<<m_q;
  else cout<<"(";
  for (vector<size_t>::iterator g_it(m_g.begin());g_it!=m_g.end();g_it++) {
    if (g_it!=m_g.begin() or m_qb!=0) cout<<",";
    cout<<*g_it;
  }
  if (m_qb!=0) cout<<","<<m_qb<<"}";
  else cout<<")";
}

//*****************************************************************************
//
// Member Functions of TraceVec
//
//*****************************************************************************

TraceVec::TraceVec(TraceType tr) {
  if (tr.is_not_empty())
    m_tr_vec.push_back(tr);
  else m_tr_vec = {};
}
TraceVec::~TraceVec(void) {
  
}
void TraceVec::push_back(TraceType tr) {
  if (tr.is_not_empty())
    m_tr_vec.push_back(tr);
  else {
    cerr<<"TraceVec::push_back() Error! Can't push empty TraceType to TraceVec."<<endl;
    exit(EXIT_FAILURE);
  }
}
TraceType& TraceVec::at(size_t i) {
  return m_tr_vec.at(i);
}
TraceType TraceVec::at(size_t i) const {
  return m_tr_vec.at(i);
}
vector<TraceVec> TraceVec::add_one_gluon(size_t g_ind) {
  size_t no_of_tr_vecs(m_tr_vec.size());
  
  vector<TraceVec> new_tr_vecs;
  for (size_t i(0);i<no_of_tr_vecs;i++) {
    TraceType tr(m_tr_vec.at(i));
    vector<TraceType> new_tr_ts(tr.add_one_gluon(g_ind));
    for (auto& tr_t : new_tr_ts) {
      TraceVec tmp_tr_vec(*this);
      tmp_tr_vec.at(i)=tr_t;
      new_tr_vecs.push_back(tmp_tr_vec);
    }
  }
  
  // open new TraceType
  TraceVec cpy(*this);
  cpy.push_back(TraceType({g_ind}));
  new_tr_vecs.push_back(cpy);
  
  return new_tr_vecs;
}
vector<TraceVec> TraceVec::conjugates() {
  vector<TraceVec> conjugate_tr_vecs({*this});
  
  for (size_t i(0);i<m_tr_vec.size();i++) {
    TraceType refl(m_tr_vec.at(i).conj());
    if (refl.is_not_empty()) {
      size_t curr_s(conjugate_tr_vecs.size());
      for (size_t j(0);j<curr_s;j++) {
	TraceVec tmp_tr_vec(conjugate_tr_vecs.at(j));
	tmp_tr_vec.at(i)=refl;
	conjugate_tr_vecs.push_back(tmp_tr_vec);
      }
    }
  }
  
  conjugate_tr_vecs.erase(conjugate_tr_vecs.begin());
  return conjugate_tr_vecs;
}
vector<size_t> TraceVec::get_indices() {
  vector<size_t> indices;
  
  for (auto& tr_t : m_tr_vec) {
    vector<size_t> tmp(tr_t.get_indices());
    indices.insert(indices.end(),tmp.begin(),tmp.end());
  }
  
  return indices;
}
size_t TraceVec::no_groups() {
  return m_tr_vec.size();
}
//bool TraceVec::is_tree_level() {
//    return true;
//    size_t n_ql(0), n_con_g(0), n_qlg(0);
//    for (auto & tr_t : m_tr_vec) {
//        size_t ng(tr_t.no_g()), nqp(tr_t.no_qp());
//        if (ng==0 and nqp!=0) n_ql++;
//        else if (ng!=0 and nqp==0) n_con_g++;
//        else n_qlg++;
//    }
//    if (n_ql>0 and n_con_g==0 and n_qlg==0) return true;
//    if (n_ql==0 and n_con_g==1 and n_qlg==0) return true;
//    if (n_con_g==0 and n_qlg>=1) return true;
//    return false;
//}
bool TraceVec::has_sg() {
  for (auto& tr_t : m_tr_vec) 
    if (tr_t.vanishes()) return true;
  return false;
}
void TraceVec::order() {
  sort(m_tr_vec.begin(), m_tr_vec.end(), [ ]( TraceType& lhs, TraceType& rhs ) {
      return lhs>rhs;
    });
}
bool TraceVec::operator>(TraceVec& rhs) {
  this->order();
  rhs.order();
  
  size_t s_1(this->no_groups()), s_2(rhs.no_groups());
  if (s_1>s_2) return false;
  if (s_1<s_2) return true;
  
  if ((this->at(0)).get_indices().size()>rhs.at(0).get_indices().size()) return true;
  if ((this->at(0)).get_indices().size()<rhs.at(0).get_indices().size()) return false;
  return this->at(0)>(rhs.at(0));
}
bool TraceVec::operator==(TraceVec& rhs) {
  size_t i(0), s_1(m_tr_vec.size()), s_2(rhs.no_groups());
  if (s_1!=s_2) return false;
  while (i<s_1 and (*this).at(i)==rhs.at(i)) i++;
  if (i==s_1) return true;
  return false;
}
CAmplitude TraceVec::build_ca(bool reduceDim) {
  CAmplitude ca;
  size_t start_ind(101);
  for (auto& tr_t : m_tr_vec) {
    ca.multiply(tr_t.build_ca(start_ind,reduceDim));
    
    size_t n_qp(tr_t.no_qp()), n_g(tr_t.no_g());
    if (n_qp!=0 and n_g>0) start_ind+=n_g-1;
    else if (n_qp==0 and n_g>2) start_ind+=n_g;
  }
  return ca;
}
void TraceVec::print() {
  cout<<"[";
  for (auto& tr : m_tr_vec) tr.print();
  cout<<"]"<<endl;
}

//*****************************************************************************
//
// Member Functions of TraceBasis
//
//*****************************************************************************

TraceBasis::TraceBasis(size_t n_g, size_t n_qp, bool reduceDim) {
  // Set basis type
  m_btype=2;
  
  // Set number of particles
  m_ng=n_g;
  m_nqp=n_qp;
    
  // Initialise process specifications
  for (size_t n(0);n<m_nqp;n++) m_process.add_out_leg("q");
  for (size_t n(0);n<m_nqp;n++) m_process.add_out_leg("qb");
  for (size_t n(0);n<m_ng;n++) m_process.add_out_leg("g");
    
  // Set conversion factor to colour flow basis
  m_confact=pow(sqrt(1./TR),m_ng);
    
  // Set whether we reduce the dimension of the basis by combining conjugate amplitudes
  m_reduce_dim = reduceDim;

  // Initialise quark, anti quark and gluon indices 
  // Note: the process is ordered as q,...,q,qb,...,qb,g,...,g
  for (size_t n(1);n<=m_nqp;n++) m_q_indices.push_back(n);
  for (size_t n(m_nqp+1);n<=2*m_nqp;n++) m_qb_indices.push_back(n);
  for (size_t n(2*m_nqp+1);n<=2*m_nqp+m_ng;n++) m_g_indices.push_back(n);
  size_t g_start(0);
  
  if (m_verbose >= 3) {
    cout<<"TraceBasis::TraceBasis() Index summary:"<<endl;
    cout<<" q indices:"<<endl;
    for (const auto& i : m_q_indices) cout<<" "<<i<<" ";
    cout<<endl;
    cout<<" qb indices:"<<endl;
    for (const auto& i : m_qb_indices) cout<<" "<<i<<" ";
    cout<<endl;
    cout<<" g indices:"<<endl;
    for (const auto& i : m_g_indices) cout<<" "<<i<<" ";
    cout<<endl;
  }

  // build all possible quark lines
  if (m_nqp!=0) {
    vector<vector<size_t>> qqb_ind_combos(get_q_ind_combinations(m_q_indices,m_qb_indices));
    
    for (const auto& qqb_inds : qqb_ind_combos) {
      TraceVec tmp_trv;
      for (size_t qp_no(0);qp_no<2*m_nqp;qp_no+=2)
	tmp_trv.push_back(TraceType(qqb_inds.at(qp_no),qqb_inds.at(qp_no+1)));
      m_tr_basis.push_back(tmp_trv);
    }
  }
  // if there are no quarks
  else {
    if (m_ng==1) {
      cerr<<"Error: no trace basis to build out of 0 quark pairs and 1 gluon."<<endl;
      exit(EXIT_FAILURE);
    }
    
    TraceVec tmp_trv(TraceType(vector<size_t>(m_g_indices.begin(),m_g_indices.begin()+1)));
    m_tr_basis.push_back(tmp_trv);
    g_start++;
  }
    
  // successively add all n_g gluons
  for (size_t i(g_start);i<m_ng;i++) {
    size_t g(m_g_indices.at(i));
    vector<TraceVec> trb_cpy;
    
    for (auto& bv : m_tr_basis) {
      vector<TraceVec> new_bvs(bv.add_one_gluon(g));
      for (auto& v : new_bvs) trb_cpy.push_back(v); 
    }
    m_tr_basis=trb_cpy;
  }
    
  this->remove_sg();
  this->normal_order();
  if (m_reduce_dim)
    this->remove_conj();
  //  for (auto& v : m_tr_basis) v.print();
  
  m_dim=m_tr_basis.size();
  for (size_t i(0);i<m_dim;i++) {
    m_normalisations.push_back(1.);
    m_norms2.push_back(ColourSum(ColourFactor(0., 0, 0, 0, 0)));
  }
  this->make_perms();
  this->make_ca_basis();
  
  // initialise matrices
  m_smat=CMatrix(m_dim);
  m_ccmats=vector<CMatrix>();

  // Basis is currently not normalised
  m_is_normalised = false;
}
TraceBasis::~TraceBasis() {
    
}

// Private functions.

void TraceBasis::remove_sg() {
  for (size_t i(0);i<m_tr_basis.size();i++) {
    if (m_tr_basis.at(i).has_sg()) {
      m_tr_basis.erase(m_tr_basis.begin()+i);
      --i;
    }
  }
}
void TraceBasis::remove_conj() {
  for (size_t i(0);i<m_tr_basis.size();i++) {
    vector<TraceVec> conjugate_tvs(m_tr_basis.at(i).conjugates());
    
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
void TraceBasis::normal_order() {
  sort(m_tr_basis.begin(), m_tr_basis.end(), [ ]( TraceVec& lhs, TraceVec& rhs ) {
      return lhs>rhs;
    });
}
void TraceBasis::make_perms() {
  for (auto& v : m_tr_basis) m_amp_perms.push_back(v.get_indices());
  //        if (v.is_tree_level()) m_amp_perms.push_back(v.get_indices());
}
void TraceBasis::make_ca_basis() {
  for (auto& bv : m_tr_basis)
    m_ca_basis.push_back(bv.build_ca(m_reduce_dim));
}

//*****************************************************************************
//
// Helper Functions
//
//*****************************************************************************

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
