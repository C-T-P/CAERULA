// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#include "FBasis.h"
#include "Colourtools.h"

// member functions of FType class
FType::FType(vector<size_t> g_inds) {
  m_g=g_inds;
}
FType::~FType(void) {
  
}
vector<FType> FType::add_one_gluon(size_t g_ind) {
  size_t n_g(m_g.size());
  
  vector<FType> new_fs;
  if (n_g==1) {
    new_fs.push_back(FType(vector<size_t>({m_g.at(0),g_ind})));
  }
  else {
    for (size_t i(1);i<n_g;i++) {
      vector<size_t> new_g_inds(m_g);
      new_g_inds.insert(new_g_inds.begin()+i,g_ind);
      new_fs.push_back(FType(new_g_inds));
    }
  }
  return new_fs;
}
vector<size_t> FType::get_indices() {
  return m_g;
}
size_t FType::no_g() {
  return m_g.size();
}
bool FType::is_not_empty() {
  if (m_g.size()>0) return true;
  else return false;
}
bool FType::vanishes() {
  if (m_g.size()<=1) return true;
  else return false;
}
bool FType::operator>(FType& rhs) {
  vector<size_t> indices1((*this).get_indices()), indices2(rhs.get_indices());
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
bool FType::operator==(FType rhs) {
  vector<size_t> indices1((*this).get_indices()), indices2(rhs.get_indices());
  size_t s_1(indices1.size()), s_2(indices2.size());
  if (s_1!=s_2) return false;
  size_t i(0);
  while (i<s_1 and indices1.at(i)==indices2.at(i)) i++;
  if (i==s_1) return true;
  return false;
}
CTerm FType::build_ct(size_t start_ind) {
    CTerm ct;
    size_t n_g(m_g.size());
    
    if (!this->vanishes()) {
      if (n_g==2) {
	ct.push_back(Delta(m_g.at(0), m_g.at(1), true));
	return ct;
      }
      if (n_g==3) {
	ct.push_back(Antisymmetric(m_g.at(0),m_g.at(1),m_g.at(2)));
	return ct;
      }
      size_t incr(0);
      ct.push_back(Antisymmetric(m_g.at(0),m_g.at(1),start_ind));
      for (size_t i(2);i<n_g-1;i++) {
	int c_ind;
	if (i==n_g-2) c_ind=m_g.at(i+1);
	else c_ind=start_ind+incr+1;
	ct.push_back(Antisymmetric(start_ind+incr,m_g.at(i),c_ind));
	incr++;
      }
      return ct;
    }
    else {
      cerr<<"Error: Cannot build colour_term from empty FType."<<endl;
      exit(EXIT_FAILURE);
    }
}
void FType::print() {
  cout<<"(";
  for (vector<size_t>::iterator g_it(m_g.begin());g_it!=m_g.end();g_it++) {
    if (g_it!=m_g.begin()) cout<<",";
    cout<<*g_it;
  }
  cout<<")";
}


// member functions of FVec class
FVec::FVec(FType f) {
  if (f.is_not_empty())
    m_FVec.push_back(f);
}
FVec::~FVec(void) {
  
}
void FVec::push_back(FType f) {
  if (f.is_not_empty())
    m_FVec.push_back(f);
  else {
    cerr<<"Error: can't push empty FType to FVec."<<endl;
    exit(EXIT_FAILURE);
  }
}
FType& FVec::at(size_t i) {
  return m_FVec.at(i);
}
FType FVec::at(size_t i) const {
  return m_FVec.at(i);
}
vector<FVec> FVec::add_one_gluon(size_t g_ind) {
  size_t no_of_FVecs(m_FVec.size());
  
  vector<FVec> new_FVecs;
  
  for (size_t i(0);i<no_of_FVecs;i++) {
    FType f(m_FVec.at(i));
    vector<FType> new_f_ts(f.add_one_gluon(g_ind));
    for (auto& f_t : new_f_ts) {
      FVec tmp_FVec(*this);
      tmp_FVec.at(i)=f_t;
      new_FVecs.push_back(tmp_FVec);
    }
  }
  
  // open new FType
  FVec cpy(*this);
  cpy.push_back(FType({g_ind}));
  new_FVecs.push_back(cpy);
  
  return new_FVecs;
}
vector<size_t> FVec::get_indices() {
  vector<size_t> indices;
  
  for (auto& f_t : m_FVec) {
    vector<size_t> tmp(f_t.get_indices());
    indices.insert(indices.end(),tmp.begin(),tmp.end());
  }
  
  return indices;
}
size_t FVec::no_groups() {
  return m_FVec.size();
}
//bool FVec::is_tree_level() {
//    if (m_FVec.size()==1) return true;
//    else return false;
//}
bool FVec::has_sg() {
  for (auto& f_t : m_FVec)
    if (f_t.vanishes()) return true;
  return false;
}
void FVec::order() {
  sort(m_FVec.begin(), m_FVec.end(), [ ]( FType& lhs, FType& rhs )
       {
	 return lhs>rhs;
       });
}
bool FVec::operator>(FVec& rhs) {
  this->order();
  rhs.order();
  
  size_t s_1(this->no_groups()), s_2(rhs.no_groups());
  if (s_1>s_2) return false;
  if (s_1<s_2) return true;
  
  if ((this->at(0)).get_indices().size()>rhs.at(0).get_indices().size()) return true;
  if ((this->at(0)).get_indices().size()<rhs.at(0).get_indices().size()) return false;
  return this->at(0)>(rhs.at(0));
}
CAmplitude FVec::build_ca() {
  CTerm ct;
  size_t start_ind(101);
  for (auto& f_t : m_FVec) {
    ct.push_back(f_t.build_ct(start_ind));
    if (f_t.no_g()>3) start_ind+=f_t.no_g()-3;
  }
  return CAmplitude(ct);
}
void FVec::print() {
  cout<<"[";
  for (auto& f : m_FVec) f.print();
  cout<<"]"<<endl;
}

// member functions of FBasis class
FBasis::FBasis(size_t n_g) {
  // set basis type
  m_btype=3;
  
  // set number of particles
  m_ng=n_g;
  
  // initialise process specifications
  for (size_t n(0); n<m_ng; n++) m_process.add_out_leg("g");
  
  // set conversion factor to colour flow basis
  m_confact=pow(sqrt(2),m_ng-2);
  
  // initialise gluon indices
  for (size_t n(1);n<=m_ng;n++) m_g_indices.push_back(n);
  
  if (m_ng<3) {
    cerr<<"Error: the adjoint basis needs at least 3 gluons."<<endl;
    exit(EXIT_FAILURE);
  }
  
  m_FBasis.push_back(FVec(FType(vector<size_t>({m_g_indices.at(0)}))));
  
  // successively add all n_g gluons
  for (size_t i(1);i<m_ng;i++) {
    size_t g(m_g_indices.at(i));
    vector<FVec> fb_cpy;
    
    for (auto& bv : m_FBasis) {
      vector<FVec> new_bvs(bv.add_one_gluon(g));
      for (auto& v : new_bvs) fb_cpy.push_back(v);
    }
    m_FBasis=fb_cpy;
  }
  
  this->remove_sg();
  this->normal_order();
  //    for (auto& v : m_FBasis) v.print();
  
  m_dim=m_FBasis.size();
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
FBasis::~FBasis() {
  
}

// private functions
void FBasis::remove_sg() {
  for (size_t i(0);i<m_FBasis.size();i++) {
    if (m_FBasis.at(i).has_sg()) {
      m_FBasis.erase(m_FBasis.begin()+i);
      i--;
    }
  }
}
void FBasis::normal_order() {
  sort(m_FBasis.begin(), m_FBasis.end(), [ ]( FVec& lhs, FVec& rhs )
       {
	 return lhs>rhs;
       });
}
void FBasis::make_perms() {
  for (auto& v : m_FBasis) m_amp_perms.push_back(v.get_indices());
  //        if (v.is_tree_level()) m_amp_perms.push_back(v.get_indices());
}
void FBasis::make_ca_basis() {
  for (auto& bv : m_FBasis)
    m_ca_basis.push_back(bv.build_ca());
}
