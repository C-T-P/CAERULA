#include<fstream>
#include<iomanip>
#include<climits>

#include "Spectrum.h"
#include "CMatrix.h"
#include "Colourtools.h"
#include "TraceBasis.h"
#include "FBasis.h"
#include "GenBasis.h"
#include "MultipletBasis.h"

// Method to print the banner on startup
void Spectrum::print_banner() {
  cout<<endl;
  cout<<"============================================================"<<endl;
  cout<<"============            \e[31mSpe\e[32mctr\e[34mum\e[0m             ";
  cout<<"==============="<<endl;
  cout<<"============---------------------------------==============="<<endl;
  cout<<"============  Calculations in Colour Space   ==============="<<endl;
  cout<<"============================================================"<<endl;
  cout<<endl;
  cout<<"Spectrum version 1.0, Copyright (C) 2019 Christian T Preuss"<<endl;
  cout<<endl;
  cout<<"=== Author =================================================\n"<<endl;
  cout<<"\tChristian Tobias Preuss\n\t";
  cout<<"School of Physics and Astronomy,\n\tMonash University,";
  cout<<"\n\t3800 Melbourne, Australia\n\t";
  cout<<"email: christian.preuss@monash.edu\n"<<endl;
  cout<<"============================================================\n"<<endl;
}

// Method to read basis from file
bool Spectrum::read_basis(string filename) {
  bool constructed(false);
  
  if (m_verbose>=3) {
    cout << "Spectrum::read_basis() Constructing colour basis from file "
         << filename << "." << endl;
    cout << "\nNOTE: amplitude permutations cannot be computed for bases"
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
      cout<<"Spectrum::construct_basis() Reading in multiplet basis for ";
      cout<<m_n_g<<" gluons and "<<m_n_qp<<" quark pairs."<<endl;
    }
    
    MultipletBasis* ortho_basis = new MultipletBasis(m_n_g, m_n_qp);
    
    if (ortho_basis != NULL) {
      constructed = true;
    
      // TODO: put construct_bcm in class CBasis 
      // to make it available to all bases
      if (m_construct_bcm) {
        ortho_basis->bcm();
      }
      
      m_basis = ortho_basis;
    }
  }
  else if (m_adjointBasis) {
    if (m_verbose>=3)
      cout<<"Spectrum::construct_basis() Constructing adjoint basis for "<<m_n_g<<" gluons."<<endl;

    m_basis = new FBasis(m_n_g);
  }
  else if (m_traceBasis) {
    if (m_verbose>=3) {
      cout<<"Spectrum::construct_basis() Constructing trace basis for "<<m_n_g<<" gluons";
      cout<<" and "<<m_n_qp<<" quark pairs."<<endl;
      if (m_reduceDim) {
        cout<<"Spectrum::construct_basis() Dimension of trace basis will be reduced."<<endl;
      }
    }
    
    m_basis = new TraceBasis(m_n_g, m_n_qp, m_reduceDim);
  }
  else {
    cerr<<"Spectrum::construct_basis() Error! No basis type specified."<<endl;
    exit(EXIT_FAILURE);
  }
  
  t=clock()-t;
  if (m_verbose>=3) {
    cout<<"Spectrum::construct_basis() Computation time for basis construction: "
        <<(float)t/CLOCKS_PER_SEC<<"s."<<endl;
  }

  if (m_basis != NULL) {
    constructed = true;

    m_basis->normalise(m_largeNC);
  }
  
  return constructed;
}

CMatrix Spectrum::calculate_soft_matrix() {
  if (m_verbose>=3)
    cout << "Spectrum::calculate_soft_matrix()" 
         << " Starting calculation of the soft matrix." << endl;

  clock_t t(clock());
  CMatrix soft_matrix(m_basis->sm(m_largeNC));
  t=clock()-t;
  
  if (m_verbose>=2) {
    cout << "Spectrum::calculate_soft_matrix()"
         << " Soft matrix:" << endl;
    soft_matrix.print();
    cout <<"\n(computation time: "<<(float)t/CLOCKS_PER_SEC<<"s.)\n"<<endl;
  }
  
  if (m_calcDet) {
    t=clock();
    double softDet = soft_matrix.det().real();
    t=clock()-t;
    if (m_verbose>=3)
      cout << "Spectrum::calculate_soft_matrix() Determinant = " << softDet 
           <<" (computation time: "<<(float)t/CLOCKS_PER_SEC<<"s)"<<endl;
  }
  
  return soft_matrix;
}

vector<CMatrix> Spectrum::calculate_colour_correlators() {
  if (m_verbose>=3) {
    cout << "Spectrum::calculate_colour_correlators()"
         << " Starting calculation of the colour correlator matrices." << endl;
  }

  clock_t t(clock());
  vector<CMatrix> cc_mats(m_basis->ccms(m_largeNC));
  t=clock()-t;

  if (m_verbose>=2) {
    cout << "Spectrum::calculate_colour_correlators()"
         << " Colour Change Matrices:"<<endl;
    for (auto& ccm : cc_mats) {
      ccm.print();
      cout<<endl;
    }
    cout<<"(computation time: "<<(float)t/CLOCKS_PER_SEC<<"s)\n"<<endl;
  }

  // Check colour conservation
  if(!m_basis->check_colourcons(m_largeNC)) {
    if (m_verbose>=1) {
      cout << "Spectrum::calculate_colour_correlators()" 
           << " ERROR! Colour not conserved." << endl;
    }
  }
  else if (m_verbose>=2) {
    cout << "Spectrum::calculate_colour_correlators()" 
           << " Colour is conserved." << endl;
  }

  return cc_mats;
}

// Save colour computation to file
void Spectrum::save_to_file(string filename) {
  if (m_verbose>=3)
    cout << "Spectrum::save_to_file() Printing to file."<<endl;
  
  m_basis->print_to_file(filename, m_largeNC);
}

bool Spectrum::evaluate_ca(string expr) {
  CAmplitude* ca(NULL);
  if (expr != "") {
    ca = new CAmplitude(expr);
  }
  else {
    if (m_expr != "") ca = new CAmplitude(m_expr);
    else {
      if (m_verbose >= 1) {
	cerr << "Spectrum::evaluate_ca() Error! No colour amplitude given." << endl;
	return 1;
      }
    }
  }
  
  if (ca != NULL) {
    ca->print();
    ca->evaluate();
    ColourSum result(ca->result());
    result.simplify();
    cout << "= " << result.get_string() << endl;
    cout << "= " << result.get_cnum() << " (NC = 3)" << endl;
  }
  else return 1;

  return 0;
}

bool Spectrum::simplify_ca(string expr) {
  CAmplitude* ca(NULL);
  if (expr != "") {
    ca = new CAmplitude(expr);
  }
  else {
    if (m_expr != "") ca = new CAmplitude(m_expr);
    else {
      if (m_verbose >= 1) {
	cerr << "Spectrum::evaluate_ca() Error! No colour amplitude given." << endl;
	return 1;
      }
    }
  }
  
  if (ca != NULL) {
    CAmplitude ca(m_expr);
    ca.print();
    ca.simplify();
    cout<<"= ";
    ca.print();
  }
  else return 1;

  return 0;
}

void Spectrum::print_settings() {
  cout<<"\n=== Settings ===============================================\n"<<endl;
  cout << "\tParameter\t\tValue\t\tDefault" << endl;
  cout << "\tverbose\t\t\t" << m_verbose << "\t\t" << 1 << endl;
  if (m_verbose >= 1) {
    cout << "\tn_g\t\t\t" << m_n_g << "\t\t" << 0 << endl;
    cout << "\tn_qp\t\t\t" << m_n_qp << "\t\t" << 0 << endl;
    cout << "\tlargeNC\t\t\t" << m_largeNC << "\t\t" << 0 << endl;
    cout << "\ttraceBasis\t\t" << m_traceBasis << "\t\t" << 1 << endl;
    cout << "\tadjointBasis\t\t" << m_adjointBasis << "\t\t" << 0 << endl;
    cout << "\tmultipletBasis\t\t" << m_multipletBasis << "\t\t" << 0 << endl;
    cout << "\treduceDim\t\t" << m_reduceDim << "\t\t" << 1 << endl;
    cout << "\tnorm_basis\t\t" << m_norm_basis << "\t\t" << 1 << endl;
    cout << "\tconstruct_bcm\t\t" << m_construct_bcm << "\t\t" << 0 << endl;
    cout << "\tcalcDet\t\t\t" << m_calcDet << "\t\t" << 0 << endl;
  }
  cout<<"\n============================================================\n"<<endl;
}

void Spectrum::print_basis() {
  cout<<"\n=== Basis ==================================================\n"<<endl;
  if (m_verbose >= 1) {
    cout << "Constructed ";
    if (m_basis->is_normalised()) cout << "normalised ";
    cout << "basis with dimension " << m_basis->dim() << endl;
  }
  if (m_verbose >= 2) {
    m_basis->print();
  }
  cout<<endl;
  cout<<"============================================================\n"<<endl;
}
