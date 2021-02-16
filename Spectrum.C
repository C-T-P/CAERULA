// Copyright (C) 2021 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.                                                                               

#include<fstream>
#include<iomanip>
#include<climits>

#include "Spectrum.H"

namespace SPECTRUM {

// Method to print the banner on startup
void Spectrum::print_banner() {
  cout<<endl;
  cout<<" ============================================================"<<endl;
  cout<<" ============            \e[31mSpe\e[32mctr\e[34mum\e[0m             ";
  cout<<" =============="<<endl;
  cout<<" ============---------------------------------==============="<<endl;
  cout<<" ============  Calculations in Colour Space   ==============="<<endl;
  cout<<" ============================================================"<<endl;
  cout<<endl;
  cout<<" Spectrum version 1.1, Copyright (C) 2021 Christian T Preuss"<<endl;
  cout<<endl;
  cout<<" === Author =================================================\n"<<endl;
  cout<<"   Christian Tobias Preuss" << endl
      <<"   School of Physics and Astronomy" << endl
      <<"   Monash University," << endl
      <<"   3800 Melbourne, Australia" << endl
      <<"   email: christian.preuss@monash.edu\n"<<endl;
  cout<<" ============================================================\n"<<endl;
}

// Method to read basis from file
bool Spectrum::read_basis(string filename) {
  bool constructed(false);
  
  if (m_verbose>=3) {
    cout << " Spectrum::read_basis() Constructing colour basis from file "
         << filename << "." << endl;
    cout << "\n NOTE: amplitude permutations cannot be computed for bases"
         << " read from files!" << endl;
  }
  
  m_basis = new GenBasis(filename);
  
  if (m_basis!=NULL) { 
    constructed = true;
    if (m_norm_basis) m_basis->normalise(m_largeNC);
  }
  
  return constructed;
}

bool Spectrum::construct_basis(int nqp, int ng) {
  bool constructed(false);

  // Save number of quark pairs and gluons
  m_n_qp = nqp;
  m_n_g = ng;

  // Construct the basis
  clock_t t(clock());
  if (m_multipletBasis) {
    if (m_verbose>=3) {
      cout<<" Spectrum::construct_basis() Reading in multiplet basis for ";
      cout<<m_n_g<<" gluons and "<<m_n_qp<<" quark pairs."<<endl;
    }
    
    MultipletBasis* ortho_basis = new MultipletBasis(m_n_g, m_n_qp);
    
    if (ortho_basis != NULL) {
      constructed = true;
    
      // TODO: put construct_bcm in class CBasis to make it available to all bases
      if (m_construct_bcm) {
        ortho_basis->bcm();
      }
      
      m_basis = ortho_basis;
    }
  }  else if (m_adjointBasis) {
    if (m_verbose>=3)
      cout<<" Spectrum::construct_basis() Constructing adjoint basis for "<<m_n_g<<" gluons."<<endl;
    m_basis = new FBasis(m_n_g);
  } else if (m_traceBasis) {
    if (m_verbose>=3) {
      cout<<" Spectrum::construct_basis() Constructing trace basis for "<<m_n_g<<" gluons";
      cout<<" and "<<m_n_qp<<" quark pairs."<<endl;
      if (m_reduceDim) {
        cout<<" Spectrum::construct_basis() Dimension of trace basis will be reduced."<<endl;
      }
    }
    
    m_basis = new TraceBasis(m_n_g, m_n_qp, m_reduceDim);
  } else {
    cerr<<" Spectrum::construct_basis() Error! No basis type specified."<<endl;
    exit(EXIT_FAILURE);
  }
  
  t=clock()-t;
  if (m_verbose>=2) {
    cout<<" Spectrum::construct_basis() Computation time for basis construction: "
        <<(float)t/CLOCKS_PER_SEC<<"s."<<endl;
  }

  if (m_basis != NULL) {
    constructed = true;
    if (m_norm_basis) m_basis->normalise(m_largeNC);
  }
  
  return constructed;
}

CMatrix Spectrum::calculate_soft_matrix() {
  if (m_verbose>=1)
    cout << " Starting calculation of the soft matrix." << endl;

  clock_t t(clock());
  CMatrix soft_matrix(m_basis->sm(m_largeNC));
  t=clock()-t;
  
  if (m_verbose>=3) {
    cout << " Spectrum::calculate_soft_matrix()"
         << "  Soft matrix:" << endl;
    soft_matrix.print();
    cout <<"\n (computation time: "<<(float)t/CLOCKS_PER_SEC<<"s.)\n"<<endl;
  }
  else if (m_verbose>=2) {
    cout << " Spectrum::calculate_soft_matrix()"
         <<"  Computation time: "<<(float)t/CLOCKS_PER_SEC<<"s"<<endl;
  }
  
  if (m_calcDet) {
    t=clock();
    double softDet = soft_matrix.det().real();
    t=clock()-t;
    if (m_verbose>=3)
      cout << " Spectrum::calculate_soft_matrix() Determinant = " << softDet 
           <<"  (computation time: "<<(float)t/CLOCKS_PER_SEC<<"s)"<<endl;
  }
  
  return soft_matrix;
}

vector<CMatrix> Spectrum::calculate_colour_correlators() {
  if (m_verbose>=1)
    cout << " Starting calculation of the colour correlators." << endl;

  clock_t t(clock());
  vector<CMatrix> cc_mats(m_basis->ccms(m_largeNC));
  t=clock()-t;

  if (m_verbose>=3) {
    cout << " Spectrum::calculate_colour_correlators()"
         << " Colour Change Matrices:"<<endl;
    for (auto& ccm : cc_mats) {
      ccm.print();
      cout<<endl;
    }
    cout<<" (computation time: "<<(float)t/CLOCKS_PER_SEC<<"s)\n"<<endl;
  }
  else if (m_verbose>=2) {
    cout << " Spectrum::calculate_colour_correlators()"
         << " Computation time: "<<(float)t/CLOCKS_PER_SEC<<"s"<<endl;
  }

  // Check colour conservation
  if(!m_basis->check_colourcons(m_largeNC)) {
    if (m_verbose>=1) {
      cout << " Spectrum::calculate_colour_correlators()" 
           << " ERROR! Colour not conserved." << endl;
    }
  }
  else if (m_verbose>=2) {
    cout << " Spectrum::calculate_colour_correlators()" 
           << " Colour is conserved." << endl;
  }

  return cc_mats;
}

// Save colour computation to file
void Spectrum::save_to_file(string filename) {
  if (m_verbose>=1)
    cout << " Saving results in file."<<endl;
  m_basis->print_to_file(filename, m_largeNC);
}

bool Spectrum::evaluate_ca(string expr) {
  CAmplitude* ca(NULL);
  if (expr != "") {
    ca = new CAmplitude(expr);
  } else {
    if (m_expr != "") ca = new CAmplitude(m_expr);
    else {
      if (m_verbose >= 1) 
        cerr << " Spectrum::evaluate_ca() Error! No colour amplitude given." << endl;
      return false;
    }
  }
  
  if (ca != NULL) {
    cout << " ";
    ca->print();
    ca->evaluate();
    ColourSum result(ca->result());
    result.simplify();
    cout << " = " << result.get_string() << endl;
    cout << " = " << result.get_cnum() << " (NC = 3)" << endl;
    delete ca;
  }
  
  return true;
}

bool Spectrum::simplify_ca(string expr) {
  CAmplitude* ca(NULL);
  if (expr != "") {
    ca = new CAmplitude(expr);
  } else {
    if (m_expr != "") ca = new CAmplitude(m_expr);
    else {
      if (m_verbose >= 1)
	cerr << " Spectrum::evaluate_ca() Error! No colour amplitude given." << endl;
      return false;
    }
  }

  if (ca != NULL) {
    ca->print();
    ca->simplify();
    cout<<"= ";
    ca->print();
    delete ca;
  }

  return true;
}

void Spectrum::print_settings() {
  if (m_verbose>0) {
    cout<<"\n === Settings ===============================================\n"<<endl;
    cout << "  Parameter      " << " Value "                   << "    " << "Default" << endl;
    cout << "  ----------------------------------------------------"    << endl;
    cout << "  verbose         " << int2str(m_verbose,5)       << "    " << int2str(1,5) << endl;
    cout << "  n_g             " << int2str(m_n_g,5)           << "    " << int2str(0,5) << endl;
    cout << "  n_qp            " << int2str(m_n_qp,5)          << "    " << int2str(0,5) << endl;
    cout << "  largeNC         " << bool2str(m_largeNC)        << "    " << "False" << endl;
    cout << "  traceBasis      " << bool2str(m_traceBasis)     << "    " << "True"  << endl;
    cout << "  adjointBasis    " << bool2str(m_adjointBasis)   << "    " << "False" << endl;
    cout << "  multipletBasis  " << bool2str(m_multipletBasis) << "    " << "False" << endl;
    cout << "  reduceDim       " << bool2str(m_reduceDim)      << "    " << "True"  << endl;
    cout << "  norm_basis      " << bool2str(m_norm_basis)     << "    " << "True"  << endl;
    cout << "  construct_bcm   " << bool2str(m_construct_bcm)  << "    " << "False" << endl;
    cout << "  calcDet         " << bool2str(m_calcDet)        << "    " << "False" << endl;
    cout<<"\n ============================================================\n"<<endl;
  }
}

void Spectrum::print_basis() {
  if (m_verbose>0) {
    cout<<"\n === Basis ==================================================\n"<<endl;
    cout << "  Constructed";
    if (m_basis->is_normalised()) cout << " normalised";
    cout << " basis with dimension " << m_basis->dim() << endl;
    if (m_verbose >= 2) {
      m_basis->print();
    }
    cout<<endl;
    cout<<" ============================================================\n"<<endl;
  }
}

// Some small helper functions.
string bool2str(bool in) {
  return in ? "True " : "False";
}
string int2str(int in, int pad) {
  string str=std::to_string(in);
  for (int i(0); i<pad; ++i) str+=" ";
  return str;
}

}
