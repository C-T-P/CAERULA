// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.     

#include "Colourtools.h"

//*****************************************************************************
//
// Member functions of class process TODO: put in class Spectrum
//
//*****************************************************************************

process::process(void) {
  m_in_legs={};
  m_out_legs={};
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

string process::leg(size_t lno) {
  // leg numbering starts at 1 !
  if (lno<=m_in_legs.size()) return m_in_legs.at(lno-1).second;
  else if (lno<=m_out_legs.size()+m_in_legs.size()) return m_out_legs.at(lno-m_in_legs.size()-1).second;
  else {
    cerr << "Leg " << lno << " does not exist in diagram." << endl;
    return "";
  }
}

bool process::is_in_leg(size_t lno) {
  if (lno<=m_in_legs.size()) return true;
  else return false;
}

//*****************************************************************************
//
// Member functions of class ColourFactor.
//
//*****************************************************************************

ColourFactor::ColourFactor() {
  m_NC = 0;
  m_TR = 0;
  m_CF = 0;
  m_CA = 0;
  m_cmplx = 1.;
}

ColourFactor::ColourFactor(string expr) {
  m_NC = 0;
  m_TR = 0;
  m_CF = 0;
  m_CA = 0;
  m_cmplx = 1.;

  for (size_t mpos(expr.find('*')), dpos(expr.find('/')); dpos!=string::npos or mpos!=string::npos
	 or expr.length()>0; mpos=expr.find('*'), dpos=expr.find('/')) {
    string factor;
    int pow(dpos < mpos ? -1 : +1);
    if (mpos == string::npos and dpos == string::npos) {
      factor=expr;
      expr="";
    }
    else {
      // Choose which symbol came first
      size_t del(mpos < dpos ? mpos : dpos);
      
      size_t dpos2(expr.find('/',del+1)), mpos2(expr.find('*',del+1));
      if (dpos2 == string::npos && mpos2 == string::npos) {
	factor = expr.substr(del+1);
	expr = expr.substr(0,del);
      }
      else {
	size_t del2(dpos2 < mpos2 ? dpos2 : mpos2);
	factor = expr.substr(del+1,del2-del-1);
	expr.erase(del,del2-del);
      }
    }

    if (factor == "NC") m_NC += pow;
    else if (factor == "TR") m_TR += pow;
    else if (factor == "CF") m_CF += pow;
    else if (factor == "CA") m_CA += pow;
    else cout << "ColourFactor::ColourFactor() unknown expression "
	      << factor << " encountered." << endl;
  }
}

ColourFactor::ColourFactor(complex<double> cnum, int pow_NC, int pow_TR, int pow_CF, int pow_CA) {
  m_cmplx = cnum;
  m_NC = pow_NC;
  m_TR = pow_TR;
  m_CF = pow_CF;
  m_CA = pow_CA;
}

string ColourFactor::get_string() {
  string str("");
  if (fabs(m_cmplx.real()) <= TINY && fabs(m_cmplx.imag()) <= TINY) return str;
  if (fabs(m_cmplx.real()) > TINY) {
    str += to_string(m_cmplx.real());
  }
  if (fabs(m_cmplx.imag()) > TINY) {
    if (m_cmplx.imag() < 0.) str += "-";
    else str += "+";
    str += to_string(abs(m_cmplx.imag()))+"I";
  }
  if (m_NC != 0) str +=  " * NC^" + to_string(m_NC);
  if (m_TR != 0) {
    str += " * TR^" + to_string(m_TR);
  }
  if (m_CF != 0) {
    str += " * CF^" + to_string(m_CF);
 }
 if (m_CA != 0) {
    str += " * CA^" + to_string(m_CA);
 }
 return str;
}

void ColourFactor::operator=(string expr) {
  ColourFactor newColFac(expr);
  
  m_cmplx = newColFac.m_cmplx;
  m_NC = newColFac.m_NC;
  m_TR = newColFac.m_TR;
  m_CF = newColFac.m_CF;
  m_CA = newColFac.m_CA;
}

void ColourFactor::operator=(complex<double> z) {
  m_NC = 0;
  m_TR = 0;
  m_CF = 0;
  m_CA = 0;
  m_cmplx = z;
}

ColourFactor ColourFactor::operator*(string expr) {
  ColourFactor result(expr);
  result *= (*this);
  return result;
}

void ColourFactor::operator*=(string expr) {
  ColourFactor rhs(expr);
  m_NC += rhs.m_NC;
  m_TR += rhs.m_TR;
  m_CF += rhs.m_CF;
  m_CA += rhs.m_CA;
  m_cmplx *= rhs.m_cmplx;
}

ColourFactor ColourFactor::operator*(ColourFactor rhs) {
  ColourFactor result;

  result.m_NC = m_NC + rhs.m_NC;
  result.m_TR = m_TR + rhs.m_TR;
  result.m_CF = m_CF + rhs.m_CF;
  result.m_CA = m_CA + rhs.m_CA;
  result.m_cmplx = m_cmplx * rhs.m_cmplx;

  return result;
}

void ColourFactor::operator*=(ColourFactor rhs) {
  m_NC += rhs.m_NC;
  m_TR += rhs.m_TR;
  m_CF += rhs.m_CF;
  m_CA += rhs.m_CA;
  m_cmplx *= rhs.m_cmplx;
}

ColourFactor ColourFactor::operator*(complex<double> rhs) {
  ColourFactor result(*this);
  result.m_cmplx = rhs;
  return result;
}

void ColourFactor::operator*=(complex<double> rhs) {
  m_cmplx *= rhs;
}

ColourFactor ColourFactor::operator/(ColourFactor rhs) {
  ColourFactor result;

  result.m_NC = m_NC - rhs.m_NC;
  result.m_TR = m_TR - rhs.m_TR;
  result.m_CF = m_CF - rhs.m_CF;
  result.m_CA = m_CA - rhs.m_CA;
  result.m_cmplx = m_cmplx / rhs.m_cmplx;

  return result;
}

void ColourFactor::operator/=(ColourFactor rhs) {
  m_NC -= rhs.m_NC;
  m_TR -= rhs.m_TR;
  m_CF -= rhs.m_CF;
  m_CA -= rhs.m_CA;
  m_cmplx = m_cmplx / rhs.m_cmplx;
}

ColourFactor ColourFactor::cconj() {
  ColourFactor new_gf = ColourFactor(m_cmplx, m_NC, m_TR, m_CF, m_CA);
  new_gf.m_cmplx = conj(new_gf.m_cmplx);
  return new_gf;
}

void ColourFactor::del() {
  m_NC = 0;
  m_TR = 0;
  m_CF = 0;
  m_CA = 0;
  m_cmplx = 0.;
}

void ColourFactor::replace_CA() {
  // Replace CA by 2*TR*NC
  m_NC += m_CA;
  m_TR += m_CA;
  while (m_CA > 0) {
    m_cmplx *= 2.;
    --m_CA;
  }
  while (m_CA < 0) {
    m_cmplx *= 1./2.;
    ++m_CA;
  }
}

void ColourFactor::replace_CF() {
  // Replace CF by TR*NC (large-NC approximation)
  m_NC += m_CF;
  m_TR += m_CF;
  m_CF = 0;
}

void ColourFactor::replace_TR() {
  // Replace TR by whatever value it was assigned
  while (m_TR > 0) {
    m_cmplx *= TR;
    --m_TR;
  }
  while (m_TR < 0) {
    m_cmplx *= 1./TR;
    ++m_TR;
  }
}

int ColourFactor::get_order_NC() {
  int order(m_NC);
  order += m_CF;
  order += m_CA;
  return order;
}

complex<double> ColourFactor::get_cnum() {
  complex<double> result(m_cmplx);
  result *= complex<double>(pow(TR, m_TR),0.);
  result *= complex<double>(pow(NC, m_NC),0.);
  result *= complex<double>(pow(CF, m_CF),0.);
  result *= complex<double>(pow(CA, m_CA),0.);
  return result; 
}

complex<double> ColourFactor::get_cnum_LC() {
  complex<double> result(m_cmplx);
  result *= complex<double>(pow(TR, m_TR),0.);
  result *= complex<double>(pow(NC, m_NC),0.);
  result *= complex<double>(pow(TR*NC, m_CF),0.);
  result *= complex<double>(pow(CA, m_CA),0.);
  return result; 
  
}

complex<double> ColourFactor::get_cnum_large_NC() {
  complex<double> result(m_cmplx);

  // Replace CA first
  replace_CA();

  // Use leading colour approximation for CF: CF -> TR*NC 
  int pow_NC = m_NC + m_CF;
  int pow_TR = m_TR + m_CF;

  // Multiply result by TR^pow_TR
  result *= complex<double>(pow(TR,pow_TR),0.);

  // Now check the power of NC
  if (pow_NC < 0) return complex<double>(0.,0.);
  if (pow_NC > 0) return complex<double>(NAN,NAN);
  return result;
}


//*****************************************************************************
//
// Member functions of class ColourSum.
//
//*****************************************************************************

ColourSum::ColourSum() {
  m_cf_sum = {};
}

ColourSum::ColourSum(ColourFactor cf) {
  m_cf_sum.push_back(cf);
}

ColourSum::ColourSum(string expr) {
  for (size_t mpos(expr.find('+')); mpos!=string::npos 
	 or expr.length()>0; mpos=expr.find('+')) {
    string summand;

    if (mpos == string::npos) {
      summand = expr;
      expr = "";
    }
    else {
      summand = expr.substr(0,mpos);
      expr = expr.substr(mpos+1);
    }
    
    m_cf_sum.push_back(ColourFactor(summand));
  }
}

string ColourSum::get_string() {
  string expr("");
  for (auto& term : m_cf_sum) {
    string strNow = term.get_string();
    if (strNow != "") {
      if (expr != "") expr += " + ";
      expr += strNow;
    }
  }
  return expr;
}

void ColourSum::operator=(complex<double> z) {
  m_cf_sum.clear();
  m_cf_sum.push_back(ColourFactor(z, 0, 0, 0, 0));
}

void ColourSum::operator=(string expr) {
  ColourSum new_cs(expr);
  m_cf_sum = new_cs.m_cf_sum;
}

ColourSum ColourSum::operator+(string expr) {
  ColourSum result(expr);
  result += (*this);
  return result;
}

void ColourSum::operator+=(string expr) {
  (*this) += ColourSum(expr);
}

ColourSum ColourSum::operator+(ColourFactor rhs) {
  ColourSum new_cf_sum(*this);
  new_cf_sum+=rhs;
  return new_cf_sum;
}

void ColourSum::operator+=(ColourFactor rhs) {
  m_cf_sum.push_back(rhs);
}

ColourSum ColourSum::operator+(ColourSum rhs) {
  ColourSum new_cf_sum(rhs);
  for (auto& term : m_cf_sum) new_cf_sum+=term;
  return new_cf_sum;
}

void ColourSum::operator+=(ColourSum rhs) {
  for (auto& term : rhs.m_cf_sum) m_cf_sum.push_back(term); 
}

ColourSum ColourSum::operator*(ColourSum rhs) {
  ColourSum result;
  for (auto& term1 : rhs.m_cf_sum) {
    for (auto& term2 : this->m_cf_sum) {
      result += term1 * term2;
    }
  }
  return result;
}

void ColourSum::operator*=(ColourSum rhs) {
  (*this) = (*this) * rhs;
}

ColourSum ColourSum::operator*(string expr) {
  ColourSum result(expr);
  result *= (*this);
  return result;
}

void ColourSum::operator*=(string expr) {
  (*this) *= ColourSum(expr);
}

ColourSum ColourSum::operator*(complex<double> z) {
  ColourSum result(*this);
  for (auto& term : result.m_cf_sum) term *= z;
  return result;
}

void ColourSum::operator*=(complex<double> z) {
  (*this) = (*this) * z;
}

void ColourSum::simplify() {
  for (vector<ColourFactor>::iterator it(m_cf_sum.begin()); it < m_cf_sum.end(); ++it) {
    if (abs(it->m_cmplx.real()) < TINY && abs(it->m_cmplx.imag()) < TINY) {
      m_cf_sum.erase(it);
      --it;
    }

    for (vector<ColourFactor>::iterator it2(it+1); it2 < m_cf_sum.end(); ++it2) {
      if (it->m_NC == it2->m_NC and it->m_TR == it2->m_TR 
	  and it->m_CF == it2->m_CF and it->m_CA == it2->m_CA) {
	it->m_cmplx += it2->m_cmplx;
	m_cf_sum.erase(it2);
	--it2;
      }
    }
  }
}

ColourSum ColourSum::cconj() {
  ColourSum result(*this);
  for (auto& term : result.m_cf_sum) term = term.cconj();
  return result;
}

void ColourSum::del() {
  m_cf_sum.clear();
}

ColourFactor ColourSum::get_leading_NC() {
  ColourFactor leading_term;
  int current_order(0);
  for (auto& term : m_cf_sum) {
    int new_order(term.get_order_NC());
    if (new_order > current_order) {
      current_order = new_order;
      leading_term = term;
      leading_term.replace_CA();
      leading_term.replace_CF();
      leading_term.replace_TR();
    }
    else if (new_order == current_order) {
      ColourFactor new_term(term);
      new_term.replace_CA();
      new_term.replace_CF();
      new_term.replace_TR();
      leading_term.m_cmplx += new_term.m_cmplx;
    }
  }

  return leading_term;
}

complex<double> ColourSum::get_cnum() {
  complex<double> result(0.);
  for (auto& term : m_cf_sum) result += term.get_cnum();
  return result;
}

complex<double> ColourSum::get_cnum_LC() {
  complex<double> result(0.);
  for (auto& term : m_cf_sum) result += term.get_cnum_LC();
  return result;
}

complex<double> ColourSum::get_cnum_large_NC() {
  complex<double> result(0.);
  for (auto& term : m_cf_sum) result += term.get_cnum_large_NC();
  return result;
}

//*****************************************************************************
//
// Member functions of class Delta.
//
//*****************************************************************************

Delta::Delta(size_t i, size_t j, bool adj) {
  m_i=i;
  m_j=j;
  m_adj=adj;
}

Delta::~Delta() {
    
}

bool Delta::is_free(size_t ind) {
  if (ind > m_i and ind > m_j) return true;
  return false;
}

string Delta::build_string() {
  string str;
  if (m_adj) str+="K_[";
  else str+="k_[";
  str+=to_string(m_i)+","+to_string(m_j)+"]";
  return str;
}

//*****************************************************************************
//
// Member functions of class Fundamental.
//
//*****************************************************************************

Fundamental::Fundamental(size_t i, size_t a, size_t b) {
  m_i=i;
  m_a=a;
  m_b=b;
}

Fundamental::~Fundamental() {
  
}

bool Fundamental::is_free(size_t ind) {
  if (ind > m_i and ind > m_a and ind > m_b) return true;
  else return false;
}

string Fundamental::build_string() {
  return "t_["+to_string(m_i)+","+to_string(m_a)+","+to_string(m_b)+"]";
}

//*****************************************************************************
//
// Member functions of class Antisymmetric.
//
//*****************************************************************************

Antisymmetric::Antisymmetric(size_t i, size_t j, size_t k) {
  m_i=i;
  m_j=j;
  m_k=k;
}

Antisymmetric::~Antisymmetric() {
    
}

bool Antisymmetric::is_free(size_t ind) {
  if (ind > m_i and ind > m_j and ind > m_k) return true;
  else return false;
}

string Antisymmetric::build_string() {
  return "f_["+to_string(m_i)+","+to_string(m_j)+","+to_string(m_k)+"]";
}

//*****************************************************************************
//
// Member functions of class Symmetric.
//
//*****************************************************************************

Symmetric::Symmetric(size_t i, size_t j, size_t k) {
  m_i=i;
  m_j=j;
  m_k=k;
}

Symmetric::~Symmetric() {
  
}

bool Symmetric::is_free(size_t ind) {
  if (ind > m_i and ind > m_j and ind > m_k) return true;
  else return false;
}

string Symmetric::build_string() {
  return "d_["+to_string(m_i)+","+to_string(m_j)+","+to_string(m_k)+"]";
}

//*****************************************************************************
//
// Member functions of class CTerm.
//
//*****************************************************************************

CTerm::CTerm() {
  //   m_cnum=1.;
  m_cnum=ColourFactor();
  m_k_vec = {};
  m_t_vec = {};
  m_f_vec = {};
  m_d_vec = {};
  m_fi=1001;
}

CTerm::CTerm(Delta& k, Fundamental& t, Antisymmetric& f, Symmetric& d, ColourFactor c) {
  m_cnum=c;
  m_k_vec = {k};
  m_t_vec = {t};
  m_f_vec = {f};
  m_d_vec = {d};
  m_fi = 1001;
  while (!k.is_free(m_fi) or !t.is_free(m_fi) or !f.is_free(m_fi) or !d.is_free(m_fi)) ++m_fi;
}

CTerm::~CTerm() {
    
}

void CTerm::push_back(CTerm ct) {
  m_cnum*=ct.m_cnum;
  for (auto& k : ct.m_k_vec) this->push_back(k);
  for (auto& t : ct.m_t_vec) this->push_back(t);
  for (auto& f : ct.m_f_vec) this->push_back(f);
  for (auto& d : ct.m_d_vec) this->push_back(d);
}

void CTerm::push_back(Delta k) {
  m_k_vec.push_back(k);
  while (!k.is_free(m_fi)) ++m_fi;
}

void CTerm::push_back(Fundamental t) {
  m_t_vec.push_back(t);
  while (!t.is_free(m_fi)) ++m_fi;
}

void CTerm::push_back(Antisymmetric f) {
  m_f_vec.push_back(f);
  while (!f.is_free(m_fi)) ++m_fi;
}

void CTerm::push_back(Symmetric d) {
  m_d_vec.push_back(d);
  while (!d.is_free(m_fi)) ++m_fi;
}

void CTerm::set_cnumber(ColourFactor c) {
  m_cnum=c;
}

void CTerm::simplify() {
  do {
    this->evaluate_deltas();
    this->replace_zero();
  } while (this->replace_adjoint());
}

void CTerm::replace_zero() {
  // check if Fundamental vanishes
  for (vector<Fundamental>::iterator t_it(m_t_vec.begin()); m_cnum.get_cnum()!=0. and t_it!=m_t_vec.end(); t_it++) {
    if (t_it->m_a == t_it->m_b) {
      clear();
      m_cnum.del();
    }
  }
  
  // check if adjoint generators vanish
  for (vector<Antisymmetric>::iterator f_it(m_f_vec.begin()); m_cnum.get_cnum()!=0. and f_it!=m_f_vec.end(); f_it++) {
    // check if Antisymmetric structure constant vanishes
    if (f_it->m_i == f_it->m_j or f_it->m_i == f_it->m_k or f_it->m_j == f_it->m_k) {
      clear();
      m_cnum.del();
    }
    
    // check if contraction with Symmetric structure constant exists
    for (vector<Symmetric>::iterator d_it(m_d_vec.begin()); m_cnum.get_cnum()!=0. and d_it!=m_d_vec.end(); d_it++) {
      // check if Symmetric structure constant vanishes
      if (d_it->m_i == d_it->m_j or d_it->m_i == d_it->m_k or d_it->m_j == d_it->m_k) {
	clear();
	m_cnum.del();
      }
      else {
	if (f_it->m_i == d_it->m_j) swap<size_t>(d_it->m_i,d_it->m_j);
	else if (f_it->m_i == d_it->m_k) swap<size_t>(d_it->m_i,d_it->m_k);
	else if (f_it->m_j == d_it->m_i) {
	  swap<size_t>(f_it->m_j,f_it->m_i);
	  swap<size_t>(f_it->m_j,f_it->m_k);
	}
	else if (f_it->m_j == d_it->m_j) {
	  swap<size_t>(f_it->m_j,f_it->m_i);
	  swap<size_t>(f_it->m_j,f_it->m_k);
	  swap<size_t>(d_it->m_j,d_it->m_i);
	}
	else if (f_it->m_j == d_it->m_k) {
	  swap<size_t>(f_it->m_j,f_it->m_i);
	  swap<size_t>(f_it->m_j,f_it->m_k);
	  swap<size_t>(d_it->m_k,d_it->m_i);
	}
	else if (f_it->m_k == d_it->m_i) {
	  swap<size_t>(f_it->m_j,f_it->m_i);
	  swap<size_t>(f_it->m_j,f_it->m_k);
	}
	else if (f_it->m_k == d_it->m_j) {
	  swap<size_t>(f_it->m_k,f_it->m_i);
	  swap<size_t>(f_it->m_j,f_it->m_k);
	  swap<size_t>(d_it->m_j,d_it->m_i);
	}
	else if (f_it->m_k == d_it->m_k) {
	  swap<size_t>(f_it->m_k,f_it->m_i);
	  swap<size_t>(f_it->m_j,f_it->m_k);
	  swap<size_t>(d_it->m_k,d_it->m_i);
	}
        
	if (f_it->m_i == d_it->m_i) {
	  if (f_it->m_j == d_it->m_j or f_it->m_j == d_it->m_k) {
	    clear();
	    m_cnum.del();
	  }
	}
      }
    }
  }
}

bool CTerm::replace_adjoint() {
  bool rep_made(false);
    
  // replace contractions of two Antisymmetric structure constants
  for (vector<Antisymmetric>::iterator f_it(m_f_vec.begin()); f_it!=m_f_vec.end(); f_it++) {
    bool ev(false);
    
    for (vector<Antisymmetric>::iterator f_it2(m_f_vec.begin()); !ev and f_it2!=m_f_vec.end();f_it2++) {
      if (f_it2!=f_it) {
	if (f_it->m_i == f_it2->m_j) {
	  swap<size_t>(f_it2->m_i,f_it2->m_j);
	  swap<size_t>(f_it2->m_j,f_it2->m_k);
	}
	else if (f_it->m_i == f_it2->m_k) {
	  swap<size_t>(f_it2->m_i,f_it2->m_k);
	  swap<size_t>(f_it2->m_k,f_it2->m_j);
	}
	else if (f_it->m_j == f_it2->m_i) {
	  swap<size_t>(f_it->m_j,f_it->m_i);
	  swap<size_t>(f_it->m_k,f_it->m_j);
	}
	else if (f_it->m_j == f_it2->m_j) {
	  swap<size_t>(f_it->m_j,f_it->m_i);
	  swap<size_t>(f_it2->m_j,f_it2->m_i);
	}
	else if (f_it->m_j == f_it2->m_k) {
	  swap<size_t>(f_it->m_j,f_it->m_i);
	  swap<size_t>(f_it2->m_k,f_it2->m_i);
	}
	else if (f_it->m_k == f_it2->m_i) {
	  swap<size_t>(f_it->m_k,f_it->m_i);
	  swap<size_t>(f_it->m_k,f_it->m_j);
	}
	else if (f_it->m_k == f_it2->m_j) {
	  swap<size_t>(f_it->m_k,f_it->m_i);
	  swap<size_t>(f_it2->m_j,f_it2->m_i);
	}
	else if (f_it->m_k == f_it2->m_k) {
	  swap<size_t>(f_it->m_k,f_it->m_i);
	  swap<size_t>(f_it2->m_k,f_it2->m_i);
	}
        
	if (f_it->m_i == f_it2->m_i) {
	  if (f_it->m_j == f_it2->m_j) {
	    if (f_it->m_k == f_it2->m_k) m_cnum*=ColourFactor(1., 1, -1, 1, 1); // CA*(NC*NC-1)
	    else {
	      m_cnum *= "CA"; // ColourFactor(1.,0, 0, 0, 1); // CA
	      m_k_vec.push_back(Delta(f_it->m_k,f_it2->m_k,true));
	    }
	    m_f_vec.erase(f_it2);
	    m_f_vec.erase(f_it);
	    f_it--;
	    ev=true;
	  }
	  else if (f_it->m_j == f_it2->m_k) {
	    if (f_it->m_k == f_it2->m_j) m_cnum*=ColourFactor(-1., 1, -1, 1, 1); // -1.*CA*(NC*NC-1)
	    else {
	      m_cnum *= "CA"; //ColourFactor(-1., 0, 0, 0, 1); // -1.*CA
	      m_cnum *= -1.;
	      m_k_vec.push_back(Delta(f_it->m_k,f_it2->m_j,true));
	    }
	    m_f_vec.erase(f_it2);
	    m_f_vec.erase(f_it);
	    f_it--;
	    ev=true;
	  }
	  else if (f_it->m_k == f_it2->m_k) {
	    if (f_it->m_j == f_it2->m_j) m_cnum*=ColourFactor(1., 1, -1, 1, 1); // CA*(NC*NC-1)
	    else {
	      m_cnum *= "CA"; //ColourFactor(1., 0, 0, 0, 1); // CA
	      m_k_vec.push_back(Delta(f_it->m_j,f_it2->m_j,true));
	    }
	    m_f_vec.erase(f_it2);
	    m_f_vec.erase(f_it);
	    f_it--;
	    ev=true;
	  }
	  else if (f_it->m_k == f_it2->m_j) {
	    if (f_it->m_j == f_it2->m_k) m_cnum*=ColourFactor(-1., 1, -1, 1, 1); // -1.*CA*(NC*NC-1)
	    else {
	      m_cnum *= "CA"; //ColourFactor(-1., 0, 0, 0, 1); // -1.*CA
	      m_cnum *= -1.;
	      m_k_vec.push_back(Delta(f_it->m_j,f_it2->m_k,true));
	    }
	    m_f_vec.erase(f_it2);
	    m_f_vec.erase(f_it);
	    f_it--;
	    ev=true;
	  }
          
	  if (ev) rep_made=true;
	}
      }
    }
  }

  // TODO: express as true ColourFactor!!!
  // replace contractions of Symmetric structure constant
  for (vector<Symmetric>::iterator d_it(m_d_vec.begin()); d_it!=m_d_vec.end(); d_it++) {
    bool ev(false);
        
    for (vector<Symmetric>::iterator d_it2(m_d_vec.begin()); !ev and d_it2!=m_d_vec.end();d_it2++) {
      if (d_it2!=d_it) {
	if (d_it->m_i == d_it2->m_j) swap<size_t>(d_it2->m_i,d_it2->m_j);
	else if (d_it->m_i == d_it2->m_k) swap<size_t>(d_it2->m_i,d_it2->m_k);
	else if (d_it->m_j == d_it2->m_i) swap<size_t>(d_it->m_j,d_it->m_i);
	else if (d_it->m_j == d_it2->m_j) {
	  swap<size_t>(d_it->m_j,d_it->m_i);
	  swap<size_t>(d_it2->m_j,d_it2->m_i);
	}
	else if (d_it->m_j == d_it2->m_k) {
	  swap<size_t>(d_it->m_j,d_it->m_i);
	  swap<size_t>(d_it2->m_k,d_it2->m_i);
	}
	else if (d_it->m_k == d_it2->m_i) {
	  swap<size_t>(d_it->m_k,d_it->m_i);
	  swap<size_t>(d_it->m_k,d_it->m_j);
	}
	else if (d_it->m_k == d_it2->m_j) {
	  swap<size_t>(d_it->m_k,d_it->m_i);
	  swap<size_t>(d_it2->m_j,d_it2->m_i);
	}
	else if (d_it->m_k == d_it2->m_k) {
	  swap<size_t>(d_it->m_k,d_it->m_i);
	  swap<size_t>(d_it2->m_k,d_it2->m_i);
	}
        
	if (d_it->m_i == d_it2->m_i) {
	  if (d_it->m_j == d_it2->m_j) {
	    if (d_it->m_k == d_it2->m_k) m_cnum*=(TR*(2*NC*NC-4)-2)/NC*(NC*NC-1);
	    else {
	      m_cnum*=(TR*(2*NC*NC-4)-2)/NC;
	      m_k_vec.push_back(Delta(d_it->m_k,d_it2->m_k,true));
	    }
	    m_d_vec.erase(d_it2);
	    m_d_vec.erase(d_it);
	    d_it--;
	    ev=true;
	  }
	  else if (d_it->m_j == d_it2->m_k) {
	    if (d_it->m_k == d_it2->m_j) m_cnum*=(TR*(2*NC*NC-4)-2)/NC*(NC*NC-1);
	    else {
	      m_cnum*=(TR*(2*NC*NC-4)-2)/NC;
	      m_k_vec.push_back(Delta(d_it->m_k,d_it2->m_j,true));
	    }
	    m_d_vec.erase(d_it2);
	    m_d_vec.erase(d_it);
	    d_it--;
	    ev=true;
	  }
	  else if (d_it->m_k == d_it2->m_k) {
	    if (d_it->m_j == d_it2->m_j) m_cnum*=(TR*(2*NC*NC-4)-2)/NC*(NC*NC-1);
	    else {
	      m_cnum*=(TR*(2*NC*NC-4)-2)/NC;
	      m_k_vec.push_back(Delta(d_it->m_j,d_it2->m_j,true));
	    }
	    m_d_vec.erase(d_it2);
	    m_d_vec.erase(d_it);
	    d_it--;
	    ev=true;
	  }
	  else if (d_it->m_k == d_it2->m_j) {
	    if (d_it->m_j == d_it2->m_k) m_cnum *= (TR*(2*NC*NC-4)-2)/NC*(NC*NC-1);
	    else {
	      m_cnum*=(TR*(2*NC*NC-4)-2)/NC;
	      m_k_vec.push_back(Delta(d_it->m_j,d_it2->m_k,true));
	    }
	    m_d_vec.erase(d_it2);
	    m_d_vec.erase(d_it);
	    d_it--;
	    ev=true;
	  }
          
	  if (ev) rep_made=true;
	}
      }
    }
  }
  
  // remove traces over two generators
  for (vector<Fundamental>::iterator t_it(m_t_vec.begin()); t_it!=m_t_vec.end(); ++t_it) {
    bool ev(false);
    
    for (vector<Fundamental>::iterator t_it2(t_it+1); !ev and t_it2!=m_t_vec.end(); ++t_it2) {
      if (t_it->m_b == t_it2->m_a and t_it->m_a == t_it2->m_b) {
	if (t_it->m_i == t_it2->m_i) {
	  m_cnum *= "CF*NC"; // TR*(NC*NC-1)
	  //	  cout << "CTerm::replace_zero() m_cnum = " << m_cnum.get_string() << endl;
	}
	else {
	  m_cnum *= "TR";
	  m_k_vec.push_back(Delta(t_it->m_i,t_it2->m_i,true));
	}
	ev=true;
	m_t_vec.erase(t_it2);
	m_t_vec.erase(t_it);
	t_it--;
      }
    }
    
    if (ev) rep_made=true;
  }
  
  return rep_made;
}

void CTerm::evaluate_deltas(bool to_LC) {
  for (vector<Delta>::iterator k_it(m_k_vec.begin()); k_it!=m_k_vec.end(); ++k_it) {
    bool ev(false);
    
    // replace adjoint indices
    if (k_it->m_adj==true) {
      if (k_it->m_i == k_it->m_j) {
	if (!to_LC) m_cnum *= "NC/TR*CF"; // (NC*NC-1)
	else m_cnum *= "NC*NC";
	ev=true;
      }
            
      // replace indices in Deltas
      for (vector<Delta>::iterator k_it2(k_it+1);!ev and k_it2!=m_k_vec.end() ; ++k_it2) {
	if (k_it->m_i == k_it2->m_i) {
	  if (k_it->m_j == k_it2->m_j) {
	    if (!to_LC) m_cnum *= "NC/TR*CF"; // (NC*NC-1)
	    else m_cnum *= "NC*NC";
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
	    if (!to_LC) m_cnum *= "NC/TR*CF"; // (NC*NC-1)
	    else m_cnum *= "NC*NC";
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
	    if(!to_LC) m_cnum *= "NC/TR*CF"; // (NC*NC-1)
	    else m_cnum *= "NC*NC";
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
	    if (!to_LC) m_cnum *= "NC/TR*CF"; // (NC*NC-1)
	    else m_cnum *= "NC*NC";
	    m_k_vec.erase(k_it2);
	    ev=true;
	  }
	  else {
	    k_it2->m_j=k_it->m_i;
	    ev=true;
	  }
	}
      }
      // replace indices in Fundamentals
      for (vector<Fundamental>::iterator t_it(m_t_vec.begin()); !ev and t_it!=m_t_vec.end(); t_it++) {
	if (k_it->m_i == t_it->m_i) {
	  t_it->m_i=k_it->m_j;
	  ev=true;
	}
	else if (k_it->m_j == t_it->m_i) {
	  t_it->m_i=k_it->m_i;
	  ev=true;
	}
      }
      // replace indices in Antisymmetrics
      for (vector<Antisymmetric>::iterator f_it(m_f_vec.begin()); !ev and f_it!=m_f_vec.end(); f_it++) {
	if (k_it->m_i == f_it->m_i) {
	  f_it->m_i=k_it->m_j;
	  ev=true;
	}
	else if (k_it->m_j == f_it->m_i) {
	  f_it->m_i=k_it->m_i;
	  ev=true;
	}
	else if (k_it->m_i == f_it->m_j) {
	  f_it->m_j=k_it->m_j;
	  ev=true;
	}
	else if (k_it->m_j == f_it->m_j) {
	  f_it->m_j=k_it->m_i;
	  ev=true;
	}
	else if (k_it->m_i == f_it->m_k) {
	  f_it->m_k=k_it->m_j;
	  ev=true;
	}
	else if (k_it->m_j == f_it->m_k) {
	  f_it->m_k=k_it->m_i;
	  ev=true;
	}
      }
      // replace indices in Symmetrics
      for (vector<Symmetric>::iterator d_it(m_d_vec.begin()); !ev and d_it!=m_d_vec.end(); d_it++) {
	if (k_it->m_i == d_it->m_i) {
	  d_it->m_i=k_it->m_j;
	  ev=true;
	}
	else if (k_it->m_j == d_it->m_i) {
	  d_it->m_i=k_it->m_i;
	  ev=true;
	}
	else if (k_it->m_i == d_it->m_j) {
	  d_it->m_j=k_it->m_j;
	  ev=true;
	}
	else if (k_it->m_j == d_it->m_j) {
	  d_it->m_j=k_it->m_i;
	  ev=true;
	}
	else if (k_it->m_i == d_it->m_k) {
	  d_it->m_k=k_it->m_j;
	  ev=true;
	}
	else if (k_it->m_j == d_it->m_k) {
	  d_it->m_k=k_it->m_i;
	  ev=true;
	}
      }
      
      if (ev) {
	m_k_vec.erase(k_it);
	k_it--;
      }
    }
        
    // replace Fundamental indices
    else {
      if (k_it->m_i == k_it->m_j) {
	m_cnum *= "NC";
	ev=true;
      }
      
      // replace indices in Deltas
      for (vector<Delta>::iterator k_it2(k_it+1); !ev and k_it2!=m_k_vec.end(); ++k_it2) {
	if (k_it->m_i == k_it2->m_i) {
	  if (k_it->m_j == k_it2->m_j) {
	    m_cnum *= "NC"; 
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
	    m_cnum *= "NC";
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
	    m_cnum *= "NC";
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
	    m_cnum *= "NC";
	    m_k_vec.erase(k_it2);
	    ev=true;
	  }
	  else {
	    k_it2->m_j=k_it->m_i;
	    ev=true;
	  }
	}
      }
      
      // replace indices in Fundamentals
      for (vector<Fundamental>::iterator t_it(m_t_vec.begin()); !ev and t_it!=m_t_vec.end(); t_it++) {
	if (k_it->m_i == t_it->m_a) {
	  t_it->m_a=k_it->m_j;
	  ev=true;
	}
	else if (k_it->m_i == t_it->m_b) {
	  t_it->m_b=k_it->m_j;
	  ev=true;
	}
	else if (k_it->m_j == t_it->m_a) {
	  t_it->m_a=k_it->m_i;
	  ev=true;
	}
	else if (k_it->m_j == t_it->m_b) {
	  t_it->m_b=k_it->m_i;
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

void CTerm::shift_inds(size_t by, bool all) {
  for (vector<Delta>::iterator k_it(m_k_vec.begin()); k_it!=m_k_vec.end(); k_it++) {
    if (k_it->m_i>100 or all) k_it->m_i+=by;
    if (k_it->m_j>100 or all) k_it->m_j+=by;
  }
  for (vector<Fundamental>::iterator t_it(m_t_vec.begin()); t_it!=m_t_vec.end(); t_it++) {
    if (t_it->m_i>100 or all) t_it->m_i+=by;
    if (t_it->m_a>100 or all) t_it->m_a+=by;
    if (t_it->m_b>100 or all) t_it->m_b+=by;
  }
  for (vector<Antisymmetric>::iterator f_it(m_f_vec.begin()); f_it!=m_f_vec.end(); f_it++) {
    if (f_it->m_i>100 or all) f_it->m_i+=by;
    if (f_it->m_j>100 or all) f_it->m_j+=by;
    if (f_it->m_k>100 or all) f_it->m_k+=by;
  }
  for (vector<Symmetric>::iterator d_it(m_d_vec.begin()); d_it!=m_d_vec.end(); d_it++) {
    if (d_it->m_i>100 or all) d_it->m_i+=by;
    if (d_it->m_j>100 or all) d_it->m_j+=by;
    if (d_it->m_k>100 or all) d_it->m_k+=by;
  }
  m_fi+=by;
}

CTerm CTerm::hconj() {
  CTerm ct(*this);
  for (auto & t : ct.m_t_vec) {
    swap(t.m_a,t.m_b);
  }
  ct.m_cnum=ct.m_cnum.cconj();
  return ct;
}

CTerm CTerm::operator*(CTerm ct) {
  ct.m_cnum*=m_cnum;
    
  ct.shift_inds(m_fi,false);
  for (auto& k : m_k_vec) ct.push_back(k);
  for (auto& t : m_t_vec) ct.push_back(t);
  for (auto& f : m_f_vec) ct.push_back(f);
  for (auto& d : m_d_vec) ct.push_back(d);
    
  return ct;
}

ColourSum CTerm::result() {
  if (m_k_vec.size()==0 and m_t_vec.size()==0 and m_f_vec.size()==0 and m_d_vec.size()==0) return m_cnum;
  else return ColourFactor(complex<double>(NAN,NAN), 0, 0, 0, 0);
}

void CTerm::clear() {
  m_k_vec.clear();
  m_t_vec.clear();
  m_f_vec.clear();
  m_d_vec.clear();
  m_cnum.del();
  m_fi=1001;
}

string CTerm::build_string() {
  string str="";
  complex<double> cnumber(m_cnum.get_cnum());
  if (cnumber!=1.) str+="c_["+to_string(cnumber.real())+","+to_string(cnumber.imag())+"]";
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

void CTerm::print() {
  cout<<build_string()<<endl;
}

//*****************************************************************************
//
// Member functions of class CAmplitude
//
//*****************************************************************************

CAmplitude::CAmplitude() {
  m_result=0.;
  m_cterm_vec={};
}

CAmplitude::CAmplitude(CTerm ct) {
  m_result=0.;
  m_cterm_vec.push_back(ct);
}

CAmplitude::CAmplitude(string expr) {
  m_result=0.;
    
  for (size_t i(0), mpos(expr.find('+'));mpos!=string::npos or expr.length()>0;mpos=expr.find('+')) {
    ++i;
        
    CTerm ct;
        
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
	double real(stod(factor.substr(3,cpos-2)));
	double imaginary(stod(factor.substr(cpos+1,factor.length()-cpos-2)));
	ct.m_cnum*=complex<double>(real,imaginary);
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
	ct.push_back(Antisymmetric(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), 
				   stoi(factor.substr(c2pos+1,factor.length()-c2pos-2))));
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
	ct.push_back(Symmetric(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), 
			       stoi(factor.substr(c2pos+1,factor.length()-c2pos-2))));
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
	ct.push_back(Fundamental(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), 
				 stoi(factor.substr(c2pos+1,factor.length()-c2pos-2))));
      }
      else if (factor.find("k_[")==0 and factor[factor.length()-1]==']') {
	size_t cpos(factor.find(','));
	if (cpos==string::npos || factor.find(',',cpos+1)!=string::npos) {
	  cerr << "Invalid number of indices for k." << endl;
	  exit(EXIT_FAILURE);
	}
	int ind_i=stoi(factor.substr(3,cpos-3)), ind_j=stoi(factor.substr(cpos+1,factor.length()-cpos-2));
	ct.push_back(Delta(ind_i,ind_j,false));
      }
      else if (factor.find("K_[")==0 and factor[factor.length()-1]==']') {
	size_t cpos(factor.find(','));
	if (cpos==string::npos || factor.find(',',cpos+1)!=string::npos) {
	  cerr << "Invalid number of indices for K." << endl;
	  exit(EXIT_FAILURE);
	}
	int ind_i=stoi(factor.substr(3,cpos-3)), ind_j=stoi(factor.substr(cpos+1,factor.length()-cpos-2));
	ct.push_back(Delta(ind_i,ind_j,true));
      }
      else cerr << "Invalid input." << endl;
    }
    m_cterm_vec.push_back(ct);
  }
}

CAmplitude::~CAmplitude() {
    
}

void CAmplitude::add(CTerm ct) {
  m_cterm_vec.push_back(ct);
}

CAmplitude CAmplitude::hconj() {
  CAmplitude ca(*this);
    
  for (vector<CTerm>::iterator c_it(ca.m_cterm_vec.begin()); c_it!=ca.m_cterm_vec.end(); c_it++) (*c_it)=c_it->hconj();
    
  return ca;
}

CAmplitude CAmplitude::shift_to_internal(size_t by) {
  CAmplitude new_ca;
  for (auto ct : m_cterm_vec) {
    ct.shift_inds(by,true);
    new_ca.add(ct);
  }
  return new_ca;
}

CAmplitude CAmplitude::operator*(complex<double> z) {
  CAmplitude new_ca;
  for (vector<CTerm>::iterator c_it(m_cterm_vec.begin()); c_it!=m_cterm_vec.end(); c_it++) {
    c_it->m_cnum*=z;
    new_ca.add(*c_it);
  }
  return new_ca;
}

CAmplitude CAmplitude::operator*(CAmplitude ca){
  if (m_cterm_vec.size()==0) return (*this);
  if (ca.m_cterm_vec.size()==0) return ca;
  
  CAmplitude new_ca;
  
  for (vector<CTerm>::iterator c_it1(m_cterm_vec.begin()); c_it1!=m_cterm_vec.end(); c_it1++)
    for (vector<CTerm>::iterator c_it2(ca.m_cterm_vec.begin()); c_it2!=ca.m_cterm_vec.end(); c_it2++)
      new_ca.add((*c_it1)*(*c_it2));
  
  return new_ca;
}

void CAmplitude::multiply(CAmplitude ca) {
  if (m_cterm_vec.size()==0) {
    for (vector<CTerm>::iterator c_it(ca.m_cterm_vec.begin()); c_it!=ca.m_cterm_vec.end(); c_it++) {
      m_cterm_vec.push_back(*c_it);
    }
  }
  else {
    vector<CTerm> new_CTerms;
    
    for (vector<CTerm>::iterator c_it1(m_cterm_vec.begin()); c_it1!=m_cterm_vec.end(); c_it1++) {
      for (vector<CTerm>::iterator c_it2(ca.m_cterm_vec.begin()); c_it2!=ca.m_cterm_vec.end(); c_it2++) {
	CTerm ct(*c_it1);
	ct.push_back(*c_it2);
	new_CTerms.push_back(ct);
      }
    }
    
    m_cterm_vec=new_CTerms;
  }
}

ColourSum CAmplitude::scprod(CAmplitude ca, bool to_LC) {
  CAmplitude scp(ca.hconj());
  scp=(*this)*scp;
  //  cout << scp.build_string() << endl;
  scp.evaluate();
  // if (to_LC) scp.evaluate_LC();
  // else scp.evaluate();

  return scp.result();
}

void CAmplitude::clear() {
  m_cterm_vec.clear();
  m_result = 0.;
}

void CAmplitude::evaluate_LC() {
  while (m_cterm_vec.size()>0) {
    CAmplitude ca;

    vector<CTerm>::iterator c_it(m_cterm_vec.begin());

    // replace Kronecker Deltas and check whether term vanishes
    c_it->evaluate_deltas(true);

    // replace Antisymmetric structure constants by commutator
    for (vector<Antisymmetric>::iterator f_it(c_it->m_f_vec.begin()); f_it!=c_it->m_f_vec.end(); f_it++) {
      // replace f_{ijk} = -i/TR * (T_i)_{xy}(T_j)_{yz}(T_k)_{zx} + i/TR * (T_j)_{xy}(T_i)_{yz}(T_k)_{zx}
      size_t i(f_it->m_i), j(f_it->m_j), k(f_it->m_k);
      c_it->m_f_vec.erase(f_it);
      CTerm ct(*c_it);
      f_it--;
      
      // first term
      c_it->m_cnum*=ColourFactor(complex<double>(0.,-1.), 0, -1, 0, 0); //  complex<double>(0.,-1./TR);
      c_it->push_back(Fundamental(i,c_it->m_fi,c_it->m_fi+1));
      c_it->push_back(Fundamental(j,c_it->m_fi-1,c_it->m_fi));
      c_it->push_back(Fundamental(k,c_it->m_fi-1,c_it->m_fi-3));
              
      // second term
      ct.m_cnum*=ColourFactor(complex<double>(0.,1.), 0, -1, 0, 0); // complex<double>(0.,1./TR);
      ct.push_back(Fundamental(j,ct.m_fi,ct.m_fi+1));
      ct.push_back(Fundamental(i,ct.m_fi-1,ct.m_fi));
      ct.push_back(Fundamental(k,ct.m_fi-1,ct.m_fi-3));
      ca.add(ct);
    }

    // replace Symmetric structure constants by commutator
    for (vector<Symmetric>::iterator d_it(c_it->m_d_vec.begin()); d_it!=c_it->m_d_vec.end(); d_it++) {
      // replace d_{ijk} = 1/TR * (T_i)_{xy}(T_j)_{yz}(T_k)_{zx} + 1/TR * (T_j)_{xy}(T_i)_{yz}(T_k)_{zx}
      size_t i(d_it->m_i), j(d_it->m_j), k(d_it->m_k);
      c_it->m_d_vec.erase(d_it);
      CTerm ct(*c_it);
      d_it--;
      
      // first term
      c_it->m_cnum*=ColourFactor(1., 0, -1, 0, 0); // 1/TR;
      c_it->push_back(Fundamental(i,c_it->m_fi,c_it->m_fi+1));
      c_it->push_back(Fundamental(j,c_it->m_fi-1,c_it->m_fi));
      c_it->push_back(Fundamental(k,c_it->m_fi-1,c_it->m_fi-3));
              
      // second term
      ct.m_cnum*=ColourFactor(1., 0, -1, 0, 0); // 1./TR;
      ct.push_back(Fundamental(j,ct.m_fi,ct.m_fi+1));
      ct.push_back(Fundamental(i,ct.m_fi-1,ct.m_fi));
      ct.push_back(Fundamental(k,ct.m_fi-1,ct.m_fi-3));
      ca.add(ct);
    }

    // replace Fundamentals by Fierz identity
    for (vector<Fundamental>::iterator t_it(c_it->m_t_vec.begin()); c_it->m_cnum.get_cnum()!=0. and t_it!=c_it->m_t_vec.end(); t_it++) {
      bool ev(false);
     
      for (vector<Fundamental>::iterator t_it2(t_it+1); !ev and t_it2!=c_it->m_t_vec.end() and c_it->m_t_vec.size()>0; t_it2++) {
	if (t_it->m_i == t_it2->m_i) {
	  // replace (T_i)_{ab}(T_i)_{cd} with LC expression
	  size_t a(t_it->m_a), b(t_it->m_b), c(t_it2->m_a), d(t_it2->m_b);
	  c_it->m_t_vec.erase(t_it2);
	  c_it->m_t_vec.erase(t_it);
	  CTerm ct(*c_it);
	  t_it--;

	  c_it->m_cnum *= "TR";
	  c_it->m_k_vec.push_back(Delta(a,d,false));
	  c_it->m_k_vec.push_back(Delta(c,b,false));

//           if (t_it->m_b == t_it2->m_a) {
//             if (t_it->m_a == t_it2->m_b) {
//               c_it->m_cnum*=TR*NC*NC;
//             }
//             else {
//               c_it->m_cnum*=TR*NC;
//               c_it->m_k_vec.push_back(Delta(t_it->m_a,t_it2->m_b,false));
//             }
//             c_it->m_t_vec.erase(t_it2);
//             c_it->m_t_vec.erase(t_it);
//             t_it--;
//           }
//           else {
//             // replace (T_i)_{ab}(T_i)_{cd} with LC expression
//             size_t a(t_it->m_a), b(t_it->m_b), c(t_it2->m_a), d(t_it2->m_b);
//             c_it->m_t_vec.erase(t_it2);
//             c_it->m_t_vec.erase(t_it);
//             CTerm ct(*c_it);
//             t_it--;

//             // first term
//             c_it->m_cnum*=TR;
//             c_it->m_k_vec.push_back(Delta(a,d,false));
//             c_it->m_k_vec.push_back(Delta(b,c,false));
//           }
          ev=true;
        }
      }
      if (ev) c_it->evaluate_deltas(true);
    }

    // by now only quark lines and rings exist
    c_it->evaluate_deltas(true);
    m_result+=c_it->result();
    m_cterm_vec.erase(c_it);
        
    ca.evaluate_LC();
    m_result+=ca.result();
    
    if (m_result.get_cnum()==complex<double>(NAN,NAN)) {
      cerr<<"Error while computing a scalar product, result is nan:"<<endl;
      c_it->print();
    }
  } 
}

void CAmplitude::evaluate() {
  while (m_cterm_vec.size()>0) {
    CAmplitude ca;
    vector<CTerm>::iterator c_it(m_cterm_vec.begin());
        
    c_it->simplify();
            
    // replace Antisymmetric structure constants
    for (vector<Antisymmetric>::iterator f_it(c_it->m_f_vec.begin()); f_it!=c_it->m_f_vec.end(); f_it++) {
      bool ev(false);
      
      // gluons have to be attached either to a three gluon vertex or a quark-gluon vertex by now
      
      // replace contractions of two Antisymmetric structure constants
      for (vector<Antisymmetric>::iterator f_it2(f_it+1); !ev and f_it2!=c_it->m_f_vec.end(); f_it2++) {
	if (f_it->m_i == f_it2->m_j) {
	  swap<size_t>(f_it2->m_i,f_it2->m_j);
	  swap<size_t>(f_it2->m_j,f_it2->m_k);
	}
	else if (f_it->m_i == f_it2->m_k) {
	  swap<size_t>(f_it2->m_i,f_it2->m_k);
	  swap<size_t>(f_it2->m_j,f_it2->m_k);
	}
        
	// replace f_{abc}f_{ade}
	if (f_it->m_i == f_it2->m_i) {
	  size_t b(f_it->m_j), c(f_it->m_k), d(f_it2->m_j), e(f_it2->m_k);
	  c_it->m_f_vec.erase(f_it2);
	  c_it->m_f_vec.erase(f_it);
	  CTerm ct(*c_it);
	  f_it--;
	  ev=true;
          
	  // first term
	  c_it->m_cnum *= -1.;
	  c_it->m_cnum *= "TR/TR/TR";
	  c_it->push_back(Fundamental(d,c_it->m_fi,c_it->m_fi+1));
	  c_it->push_back(Fundamental(e,c_it->m_fi-1,c_it->m_fi));
	  c_it->push_back(Fundamental(b,c_it->m_fi-1,c_it->m_fi));
	  c_it->push_back(Fundamental(c,c_it->m_fi-1,c_it->m_fi-4));
          
          
	  // second term
	  CTerm ct1(ct);
	  ct1.m_cnum *= "TR/TR/TR";
	  ct1.push_back(Fundamental(b,ct1.m_fi,ct1.m_fi+1));
	  ct1.push_back(Fundamental(d,ct1.m_fi-1,ct1.m_fi));
	  ct1.push_back(Fundamental(e,ct1.m_fi-1,ct1.m_fi));
	  ct1.push_back(Fundamental(c,ct1.m_fi-1,ct1.m_fi-4));
	  ca.add(ct1);
          
	  // third term
	  ct1=ct;
	  ct1.m_cnum *= "TR/TR/TR";
	  ct1.push_back(Fundamental(e,ct1.m_fi,ct1.m_fi+1));
	  ct1.push_back(Fundamental(d,ct1.m_fi-1,ct1.m_fi));
	  ct1.push_back(Fundamental(b,ct1.m_fi-1,ct1.m_fi));
	  ct1.push_back(Fundamental(c,ct1.m_fi-1,ct1.m_fi-4));
	  ca.add(ct1);
                    
	  // fourth term
	  ct.m_cnum *= -1.;
	  ct.m_cnum *= "TR/TR/TR";
	  ct.push_back(Fundamental(b,ct.m_fi,ct.m_fi+1));
	  ct.push_back(Fundamental(e,ct.m_fi-1,ct.m_fi));
	  ct.push_back(Fundamental(d,ct.m_fi-1,ct.m_fi));
	  ct.push_back(Fundamental(c,ct.m_fi-1,ct.m_fi-4));
	  ca.add(ct);
	}
      }
      
      // replace contractions of an Antisymmetric structure constant with a Symmetric one
      for (vector<Symmetric>::iterator d_it(c_it->m_d_vec.begin()); !ev and d_it!=c_it->m_d_vec.end() ; d_it++) {
	if (f_it->m_i == d_it->m_j) swap<size_t>(d_it->m_i,d_it->m_j);
	else if (f_it->m_i == d_it->m_k) swap<size_t>(d_it->m_i,d_it->m_k);
        
	// replace f_{abc}d_{ade}
	if (f_it->m_i == d_it->m_i) {
	  size_t b(f_it->m_j), c(f_it->m_k), d(d_it->m_j), e(d_it->m_k);
	  c_it->m_d_vec.erase(d_it);
	  c_it->m_f_vec.erase(f_it);
	  CTerm ct(*c_it);
	  f_it--;
	  ev=true;
          
	  // first term
	  c_it->m_cnum *= complex<double>(0.,-1.);
	  c_it->m_cnum *= "TR/TR/TR";
	  c_it->push_back(Fundamental(b,c_it->m_fi,c_it->m_fi+1));
	  c_it->push_back(Fundamental(c,c_it->m_fi-1,c_it->m_fi));
	  c_it->push_back(Fundamental(d,c_it->m_fi-1,c_it->m_fi));
	  c_it->push_back(Fundamental(e,c_it->m_fi-1,c_it->m_fi-4));
          
                    
	  // second term
	  CTerm ct1(ct);
	  ct1.m_cnum *= complex<double>(0.,1.);
	  ct1.m_cnum *= "TR/TR/TR";
	  ct1.push_back(Fundamental(c,ct1.m_fi,ct1.m_fi+1));
	  ct1.push_back(Fundamental(b,ct1.m_fi-1,ct1.m_fi));
	  ct1.push_back(Fundamental(d,ct1.m_fi-1,ct1.m_fi));
	  ct1.push_back(Fundamental(e,ct1.m_fi-1,ct1.m_fi-4));
	  ca.add(ct1);
                    
	  // third term
	  ct1=ct;
	  ct1.m_cnum *= "TR/TR/TR";
	  ct1.m_cnum *= complex<double>(0.,-1.);
	  ct1.push_back(Fundamental(d,ct1.m_fi,ct1.m_fi+1));
	  ct1.push_back(Fundamental(b,ct1.m_fi-1,ct1.m_fi));
	  ct1.push_back(Fundamental(c,ct1.m_fi-1,ct1.m_fi));
	  ct1.push_back(Fundamental(e,ct1.m_fi-1,ct1.m_fi-4));
	  ca.add(ct1);
                    
	  // fourth term
	  ct.m_cnum *= "TR/TR/TR";
	  ct.m_cnum *= complex<double>(0.,1.);
	  ct.push_back(Fundamental(d,ct.m_fi,ct.m_fi+1));
	  ct.push_back(Fundamental(c,ct.m_fi-1,ct.m_fi));
	  ct.push_back(Fundamental(b,ct.m_fi-1,ct.m_fi));
	  ct.push_back(Fundamental(e,ct.m_fi-1,ct.m_fi-4));
	  ca.add(ct);
	}
      }
      
      // replace contractions of an Antisymmetric structure constant with a Fundamental generator
      for (vector<Fundamental>::iterator t_it(c_it->m_t_vec.begin()); !ev and t_it!=c_it->m_t_vec.end(); t_it++) {
	// replace f_{abc}(T_a)_{ij}
	if (f_it->m_i == t_it->m_i) {
	  size_t b(f_it->m_j), c(f_it->m_k), i(t_it->m_a), j(t_it->m_b);
	  c_it->m_f_vec.erase(f_it);
	  c_it->m_t_vec.erase(t_it);
	  CTerm ct(*c_it);
	  f_it--;
	  ev=true;
          
	  // first term
	  c_it->m_cnum *= complex<double>(0.,-1.);
	  c_it->push_back(Fundamental(b,i,c_it->m_fi));
	  c_it->push_back(Fundamental(c,c_it->m_fi-1,j));
          
	  // second term
	  ct.m_cnum *= complex<double>(0.,1.);
	  ct.push_back(Fundamental(c,i,ct.m_fi));
	  ct.push_back(Fundamental(b,ct.m_fi-1,j));
	  ca.add(ct);
	}
      }
    }
    
    // replace remaining Symmetric structure constants
    for (vector<Symmetric>::iterator d_it(c_it->m_d_vec.begin()); d_it!=c_it->m_d_vec.end();d_it++) {
      // replace d_{abc}
      size_t a(d_it->m_i), b(d_it->m_j), c(d_it->m_k);
      c_it->m_d_vec.erase(d_it);
      CTerm ct(*c_it);
      d_it--;
      
      // first term
      c_it->m_cnum *= "TR/TR/TR";
      c_it->push_back(Fundamental(a,c_it->m_fi,c_it->m_fi+1));
      c_it->push_back(Fundamental(b,c_it->m_fi-1,c_it->m_fi));
      c_it->push_back(Fundamental(c,c_it->m_fi-1,c_it->m_fi-3));
      
      // second term
      ct.m_cnum *= "TR/TR/TR";
      ct.push_back(Fundamental(b,ct.m_fi,ct.m_fi+1));
      ct.push_back(Fundamental(a,ct.m_fi-1,ct.m_fi));
      ct.push_back(Fundamental(c,ct.m_fi-1,ct.m_fi-3));
      ca.add(ct);
    }
        
    // replace Fundamentals by Fierz identity
    for (vector<Fundamental>::iterator t_it(c_it->m_t_vec.begin()); c_it->m_cnum.get_cnum()!=0. and t_it!=c_it->m_t_vec.end(); t_it++) {
      bool ev(false);
      
      // contract indices with Fundamentals
      for (vector<Fundamental>::iterator t_it2(t_it+1); !ev and t_it2!=c_it->m_t_vec.end() and c_it->m_t_vec.size()>0; t_it2++) {
	if (t_it->m_i == t_it2->m_i) {
	  if (t_it->m_b == t_it2->m_a) {
	    if (t_it->m_a == t_it2->m_b) {
	      c_it->m_cnum *= "CF*NC"; // TR*(NC*NC-1)
	    }
	    else {
	      c_it->m_cnum *= "CF";
	      c_it->m_k_vec.push_back(Delta(t_it->m_a,t_it2->m_b,false));
	    }
	    c_it->m_t_vec.erase(t_it2);
	    c_it->m_t_vec.erase(t_it);
	    t_it--;
	  }
	  else {
	    // replace (T_i)_{ab}(T_i)_{cd}
	    size_t a(t_it->m_a), b(t_it->m_b), c(t_it2->m_a), d(t_it2->m_b);
	    c_it->m_t_vec.erase(t_it2);
	    c_it->m_t_vec.erase(t_it);
	    CTerm ct(*c_it);
	    t_it--;
                        
	    // first term
	    c_it->m_cnum *= "TR";
	    c_it->m_k_vec.push_back(Delta(a,d,false));
	    c_it->m_k_vec.push_back(Delta(b,c,false));
            
	    // second term
	    ct.m_cnum *= "TR/NC";
	    ct.m_cnum *= -1.; 
	    ct.m_k_vec.push_back(Delta(a,b,false));
	    ct.m_k_vec.push_back(Delta(c,d,false));
	    ca.add(ct);
	  }
	  ev=true;
	}
      }
      if (ev) c_it->evaluate_deltas();
    }
    
    // by now only quark lines and rings exist
    c_it->evaluate_deltas();
    m_result += c_it->result();
    m_cterm_vec.erase(c_it);
        
    ca.evaluate();
    m_result+=ca.result();
    
    if (m_result.get_cnum()==complex<double>(NAN,NAN)) {
      cerr<<"Error while computing a scalar product, result is nan:"<<endl;
      c_it->print();
    }
  }
}

void CAmplitude::simplify() {
  for (vector<CTerm>::iterator c_it(m_cterm_vec.begin()); c_it!=m_cterm_vec.end(); ++c_it)
    c_it->simplify();
}

size_t CAmplitude::no_of_terms() {
  return m_cterm_vec.size();
}

string CAmplitude::build_string() {
  string str("");
  
  for (auto& ct : m_cterm_vec) {
    if (str!="") str+="+";
    str+=ct.build_string();
  }
  
  return str;
}

void CAmplitude::print() {
  cout<<this->build_string()<<endl;
}

ColourSum CAmplitude::result() {
  if (m_cterm_vec.size()==0) return m_result;
  else return ColourSum(ColourFactor(complex<double>(NAN,NAN), 0, 0, 0, 0));
}
