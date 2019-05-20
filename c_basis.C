// Copyright (C) 2018 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#include<iomanip>
#include<fstream>
#include "c_basis.h"
#include "insert.h"
#include "multiplet_basis.h"

double eps = 1.e-4;

void c_basis::normalise(bool to_LC) {
  for (size_t i(0);i<m_dim;i++) {
    ColourSum norm2 = m_ca_basis.at(i).scprod(m_ca_basis.at(i), to_LC);
    //    cout << "c_basis::normalise() " << norm2.get_string() << endl;
    m_normalisations.at(i) = sqrt(abs(norm2.get_cnum()));
    //    cout << "c_basis::normalise() " << m_normalisations.at(i) << endl;
    m_norms2.at(i) = norm2;
    m_ca_basis.at(i)=m_ca_basis.at(i)*complex<double>(1./m_normalisations.at(i),0.);
  }
}

size_t c_basis::dim() {
  return m_ca_basis.size();
}

void c_basis::print() {
  for (size_t i(0); i<m_dim; i++) {
    cout<< "b_" << i+1 << " = ";
    m_ca_basis.at(i).print();
  }
}

void c_basis::print_to_file(string filename, bool to_LC) {
  string out_filename="";
  if (filename=="") {
    size_t n_q(0), n_qb(0), n_g(0);
    for (size_t lno(1);lno<=m_process.no_of_legs();lno++) {
      string ptcl=m_process.leg(lno);
      
      if (ptcl=="q") n_q++;
      else if (ptcl=="qb") n_qb++;
      else if (ptcl=="g") n_g++;
    }
    if (n_q!=0) out_filename+=to_string(n_q)+"q";
    if (n_qb!=0) out_filename+=to_string(n_qb)+"qb";
    if (n_g!=0) out_filename+=to_string(n_g)+"g";
    
    if (m_btype==1) out_filename+="_multiplet";
    else if (m_btype==2) out_filename+="_trace";
    else if (m_btype==3) out_filename+="_adjoint";
  }
  else out_filename=filename;
  
  ofstream file;
  file.open(out_filename+".dat");
    
  // Print info if file is leading colour
  if (to_LC) file << "% NOTE: This is a LC file." << endl;

  // print amplitude permutations and normalisations for trace and adjoint basis (and multiplet basis if calculated)
  if (m_amp_perms.size()>0) {
    file<<"% number of basis vectors with non-zero coefficient"<<endl;
    file<<"SIZE_CONNECTED = "<<m_amp_perms.size()<<";\n"<<endl;
    file<<"% prefactors of the partial amplitude in total amplitude"<<endl;
    for (size_t i(0);i<m_amp_perms.size();i++)
      file<<"a_"<<i<<" = "<<fixed<<setprecision(17)<<m_normalisations.at(i)*m_confact<<";"<<endl;
    
    file<<"\n% permutations defining the partial amplitudes"<<endl;
    for (size_t i(0);i<m_amp_perms.size();i++) {
      file<<"A_"<<i<<" =";
      for (const auto& ind : m_amp_perms[i]) file<<" "<<ind-1;
      file<<";"<<endl;
    }
    file<<"\n"<<endl;
  }
  
  // print basis change matrix if basis is a multiplet basis
  if (m_btype==1) {
    multiplet_basis *ortho_basis = static_cast<multiplet_basis*>(this);
    
    if (ortho_basis->m_bcm.size()>0) {
      file<<"% basis change matrix from trace basis to multiplet basis"<<endl;
      file<<"BCM =";
      for (size_t m(0); m<ortho_basis->m_bcm.size(); m++)
	for (size_t n(0); n<ortho_basis->m_bcm.size(); n++)
	  file<<fixed<<setprecision(17)<<" "<<ortho_basis->m_bcm[m][n].real();
      file<<";\n"<<endl;
    }
  }
  
  file<<"% dimension of the colour space"<<endl;
  file<<"DIM = "<<m_dim<<";"<<endl;
  
  file<<"\n% colour metric aka soft matrix"<<endl;
  file<<"METRIC =";
  for (size_t m(0);m<m_dim;m++)
    for (size_t n(0);n<m_dim;n++) file<<fixed<<setprecision(17)<<" "<<m_smat[m][n].real();
  file<<";"<<endl;
  
  file<<"\n% colour change matrices"<<endl;
  size_t i(0);
  for (unsigned int lno1(1);lno1<=m_process.no_of_legs();lno1++) {
    for (unsigned int lno2(lno1+1);lno2<=m_process.no_of_legs();lno2++) {
      file<<"C_"<<lno1-1<<lno2-1<<" =";
      c_matrix ccm(m_ccmats.at(i));
      for (size_t m(0);m<m_dim;m++)
	for (size_t n(0);n<m_dim;n++)
	  file<<fixed<<setprecision(17)<<" "<<ccm[m][n].real();
      file<<";"<<endl;
      i++;
    }
  }
  
  file.close();
}

c_matrix c_basis::sm(bool to_LC) {
  for (size_t i(0); i<m_dim; ++i) {
    for (size_t j(0); j<m_dim; ++j) {
      if (!to_LC) {
	// cout << "\rCalculating S" << "[" << i << "][" << j << "]..." << flush;
	complex<double> r;
	if (j >= i) r = m_ca_basis.at(i).scprod(m_ca_basis.at(j)).get_cnum();
	else r = m_smat[j][i];
	m_smat[i][j] = (abs(r) < eps ? 0. : r);
      }
      else {
	ColourSum r;
	if ( j >= i) { 
	  r = m_ca_basis.at(i).scprod(m_ca_basis.at(j), to_LC);

	  ColourFactor num(r.get_leading_NC());
	  ColourFactor denom(m_norms2.at(i).get_leading_NC() * m_norms2.at(j).get_leading_NC());
	  ColourFactor frac((num * num)/denom);
	  
	  m_smat[i][j] = sqrt(frac.get_cnum_large_NC());
	}
	else m_smat[i][j] = m_smat[j][i];
      }
    }
  }
  
  return m_smat;
}

c_matrix c_basis::ccm(size_t lno1, size_t lno2, bool to_LC) {
  c_matrix colour_cm(m_dim);
  
  c_amplitude ins_op=construct_insertion_op(m_process,lno1,lno2);
  
  for (size_t i(0); i<m_dim; ++i) {
    for (size_t j(0); j<m_dim; ++j) {
      cout << "\rCalculating C_(" << lno1 << "," << lno2 << ")[" << i << "][" << j << "]..." << flush;
      c_amplitude bvj(m_ca_basis.at(j).shift_to_internal(2000));
      bvj.multiply(ins_op);
      if (!to_LC) {
	complex<double> r(0.);
	r = m_ca_basis.at(i).scprod(bvj).get_cnum();
	colour_cm[i][j] = (abs(r) < eps ? 0. : r);
      }
      else {
	ColourSum r = m_ca_basis.at(i).scprod(bvj, to_LC);
	ColourSum v1 = bvj.scprod(bvj, to_LC);

	ColourFactor num(r.get_leading_NC());
	ColourFactor denom(m_norms2.at(i).get_leading_NC() * v1.get_leading_NC());
	ColourFactor frac((num * num)/denom);
	
	//	cout << "c_basis::ccm() " << frac.get_string() << endl;

	colour_cm[i][j] = sqrt(frac.get_cnum_large_NC());
      }
    }
  }
  
  return colour_cm;
}

vector<c_matrix> c_basis::ccms(bool to_LC) {
  for (unsigned int lno1(1);lno1<=m_process.no_of_legs();lno1++) {
    for (unsigned int lno2(lno1+1);lno2<=m_process.no_of_legs();lno2++) {
      m_ccmats.push_back(this->ccm(lno1, lno2, to_LC));
    }
  }
  return m_ccmats;
}

bool c_basis::check_colourcons() {
  bool is_conserved(true);
  size_t no_legs(m_process.no_of_legs());
  complex<double> casimir(0.);
  for (size_t i(1);i<=no_legs;i++) {
    if (m_process.leg(i)=="q" or m_process.leg(i)=="qb") casimir+=CF;
    else if (m_process.leg(i)=="g") casimir+=CA;
  }
  casimir*=0.5;
    
  //    c_matrix casimir_mat(m_smat);
  //    for (size_t i(0); i<m_dim; ++i)
  //        for (size_t j(0); j<m_dim; ++j)
  //            casimir_mat[i][j]*=casimir;
  //
  //    c_matrix sum_ccmats(m_dim);
  //    for (size_t i(0); i<no_ccmats; ++i)
  //        sum_ccmats+=m_ccmats.at(i);
  
  size_t no_ccmats(m_ccmats.size());
  for (size_t i(0); i<m_dim; ++i) {
    for (size_t j(0); j<m_dim; ++j) {
      complex<double> sum_all(casimir*m_smat[i][j]);
      for (size_t k(0); k<no_ccmats; ++k)
	sum_all+=m_ccmats[k][i][j];
      if (abs(sum_all)>eps) is_conserved=false;
    }
  }
  
  //    cout<<"Sum TProds:"<<endl;
  //    sum_ccmats.print();
  //    cout<<"Casimir Matrix:"<<endl;
  //    casimir_mat.print();
  
  return is_conserved;
}




