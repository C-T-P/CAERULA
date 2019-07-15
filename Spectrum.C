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

/// function to brint the banner on startup
void Spectrum::print_banner() {
  cout<<endl;
  cout<<"============================================================"<<endl;
  cout<<"============            \e[31mSpe\e[32mctr\e[34mum\e[0m             ";
  cout<<"==============="<<endl;
  cout<<"============---------------------------------==============="<<endl;
  cout<<"============  Calculations in Colour Space   ==============="<<endl;
  cout<<"============================================================"<<endl;
  cout<<endl;
  cout<<"Spectrum version 2.0, Copyright (C) 2019 Christian T Preuss"<<endl;
  cout<<endl;
  cout<<"=== Author =================================================\n"<<endl;
  cout<<"\tChristian Tobias Preuss\n\t";
  cout<<"School of Physics and Astronomy,\n\tMonash University,";
  cout<<"\n\t3800 Melbourne, Australia\n\t";
  cout<<"email: christian.preuss@monash.edu\n"<<endl;
  cout<<"============================================================\n"<<endl;
}

// Function to print a help message to cout
void Spectrum::print_help_message() {
  cout<<"\nOptions:"<<endl;

  cout<<"\t-h, -help:\t\tPrint this message and exit."<<endl;

  cout<<endl;

  cout<<"\t-e, -evaluate:";
  cout<<"\t\tEvaluate colour amplitude in which all indices must be contracted."<<endl;

  cout<<"\t-s, -simplify:";
  cout<<"\t\tSimplify colour term by replacing internal quark and gluon rings.";
  cout<<"Returns colour amplitude."<<endl;
  cout<<"\t\t\t\tReturns complex number or NAN if not all indices are contracted."<<endl;

  cout<<endl;

  cout<<"\t-f, -file:";
  cout<<"\t\tRead process and colour basis from input file and compute the soft matrix";
  cout<<"\n\t\t\t\tand all colour change matrices."<<endl;

  cout<<"\t-ng:\t\t\tSpecify number of gluons to construct the trace basis";
  cout<<"\n\t\t\t\t and compute the soft matrix and all colour change matrices."<<endl;

  cout<<"\t-nqp:\t\t\tSpecify number of quark pairs to construct the colour flow/trace basis";
  cout<<"\n\t\t\t\tand compute the soft matrix and all colour change matrices."<<endl;

  cout<<"\t-adj:\t\t\tBuild adjoint basis (f-basis) instead of trace basis.";
  cout<<"\n\t\t\t\tWorks only for pure gluon processes with ng>=3."<<endl;

  cout<<"\t-multiplet:\t\tBuild multiplet basis (orthogonal basis).";
  cout<<"\n\t\t\t\tWorks only if a pre-calculated multiplet basis is provided."<<endl;

  cout<<"\t-largeNC:\t\tUse large-NC limit.";

  cout<<endl;

  cout<<"\t-bcm:\t\t\tBuild basis change matrix from trace basis to multiplet basis.";
  cout<<"\n\t\t\t\tWorks only together with the -multiplet option and if a ";
  cout<<"\n\t\t\t\tprecalculated multiplet basis for the process is provided."<<endl;


  cout<<"\t-nonorm:\t\tDeactivates normalisation of basis vectors."<<endl;

  cout<<endl;

  cout<<"\t-v, -verbose:";
  cout<<"\t\tVerbosity, choose";
  cout<<"\n\t\t\t\t0: silent mode - no output except for errorst";
  cout<<"\n\t\t\t\t1: default mode - additionally settings and diagnostics";
  cout<<"\n\t\t\t\t2: debug mode - additionally  computation infos";
  cout<<"\n\t\t\t\t3: full mode - complete output for each term"<<endl;

  cout<<"\t-nooutfile:";
  cout<<"\t\tDo not save computation to a file."<<endl;
  
  cout<<endl;

  exit(EXIT_SUCCESS);
}

/// Print input error message
void Spectrum::run_error() {
  cerr<<"Please specify EITHER a colour term to be simplified (-s)/evaluated (-e)";
  cerr<<"OR perform a colour space calculation by specifying a basis through a file name (-f)";
  cerr<<"OR the process by the number of quark pairs (-nqp) and number of gluons (-ng).";
  cerr<<"See also the help menu for more information (-h)."<<endl;
  exit(EXIT_FAILURE);
}

/// Print basis error message
void Spectrum::basis_error() {
  cerr<<"Two different basis types were given.";
  cerr<<"Either the adjoint basis (-adj) OR the multiplet basis (-multiplet) can bes used.";
  exit(EXIT_FAILURE);
}

/// Initialise from vector of arguments
void Spectrum::init(int argc, char **argv) {

  // No input given
  if (argc==1) print_help_message();
  
//   // Standard settings
//   m_runoption = 0; m_n_g = 0; m_n_qp = 0;
//   m_expr = ""; m_in_filename = ""; m_out_filename = ""; m_no_outfile = false;
//   m_traceBasis = true; m_adjointBasis = false; m_multipletBasis = false;
//   m_norm_basis = true; m_construct_bcm = false; m_verbose = 1;
  
  for(int i(1); i<argc; ++i) {
    /// Print help message and exit
    if (strcmp(argv[i], "-help")==0 or strcmp(argv[i], "-h")==0) {
      print_help_message();
    }

    /// Set verbosity
    else if (strcmp(argv[i], "-verbose")==0 or strcmp(argv[i], "-v")==0) {
      m_verbose=stoi(argv[i+1]);
      ++i;
    }

    /// Evaluate a given colour amplitude
    else if (strcmp(argv[i], "-evaluate")==0 or strcmp(argv[i], "-e")==0) {
      if (m_runoption==0) m_runoption=1;
      else run_error();
      m_expr=argv[i+1];
      ++i;
    }

    /// Simplify a given colour amplitude
    else if (strcmp(argv[i], "-simplify")==0 or strcmp(argv[i], "-s")==0) {
      if (m_runoption==0) m_runoption=2;
      else run_error();
      m_expr=argv[i+1];
      ++i;
    }

    /// Do full colour computation for process read from file
    else if (strcmp(argv[i], "-f")==0 or strcmp(argv[i], "-file")==0) {
      if (m_runoption==0) m_runoption=3;
      else run_error();
      m_in_filename=argv[i+1];
      ++i;
    }

    /// Do full colour computation for process specified by partons
    else if (strcmp(argv[i],"-ng")==0) {
      if (m_runoption==0 or m_runoption==4) m_runoption=4;
      else run_error();
      m_n_g=stoi(argv[i+1]);
      ++i;
    }
    else if (strcmp(argv[i],"-nqp")==0) {
      if (m_runoption==0 or m_runoption==4) m_runoption=4;
      else run_error();
      m_n_qp=stoi(argv[i+1]);
      ++i;
    }

    /// Reduce the dimension of the trace basis by summing conjugates
    else if (strcmp(argv[i],"-reduceDim")==0) {
      if (m_traceBasis) m_reduceDim=true;
      else run_error();
      if (m_verbose>=1)
	cout<<"Spectrum::Spectrum() Dimension of trace basis will be reduced."<<endl;
    }

    /// Use the adjoint basis
    else if (strcmp(argv[i],"-adj")==0) {
      if (!m_multipletBasis) m_adjointBasis=true;
      else run_error();
      m_traceBasis=false;
      if (m_verbose>=1)
	cout<<"Spectrum::Spectrum() Adjoint basis will be used."<<endl;
    }

    /// Use the multiplet basis
    else if (strcmp(argv[i],"-multiplet")==0) {
      if (!m_adjointBasis) m_multipletBasis=true;
      else run_error();
      m_traceBasis=false;
      if (m_verbose>=1)
	cout<<"Spectrum::Spectrum() Multiplet basis will be used."<<endl;
    }

    /// Use large-NC limit
    else if (strcmp(argv[i],"-largeNC")==0) {
      m_largeNC = true;
      if (m_verbose>=1)
	cout<<"Spectrum::Spectrum() Large-NC limit will be used."<<endl;
    }   

    /// Do not normalise basis
    else if (strcmp(argv[i],"-nonorm")==0) {
      m_norm_basis=false;
      if (m_verbose>=1)
	cout<<"Spectrum::Spectrum() Basis will not be normalised."<<endl;
    }

    /// Construct basis change matrix to trace basis
    else if (strcmp(argv[i],"-bcm")==0) {
      m_construct_bcm=true;
      if (m_verbose>=1)
	cout<<"Spectrum::Spectrum() Basis change matrix will be computed."<<endl;
    }

    /// Calculate determinant of soft matrix
    else if (strcmp(argv[i],"-det")==0) {
      m_calcDet=true;
      if (m_verbose>=1)
	cout<<"Spectrum::Spectrum() Determinant of soft matrix will be calculated."<<endl;
    }

    /// Do not save calculation to file
    else if (strcmp(argv[i],"-nooutfile")==0) {
      m_no_outfile=true;
      if (m_verbose>=1)
	cout<<"Spectrum::Spectrum() No output file will be saved."<<endl;
    }

    else {
      if (m_verbose >= 1)
	cout << "Spectrum::Spectrum() unrecognised option " << argv[i] << ". Will skip." << endl;
    }
  }  
}

/// Perform colour computation
bool Spectrum::start() {
  CBasis* basis=NULL;

  // Perform colour calculations depending on the given input                                     
  switch (m_runoption) {
    case 1: {
      return evaluate_ca();
    }
    case 2: {
      return simplify_ca();
    }
    case 3: {
      if (m_verbose>=1) {
	cout<<"Spectrum::start() Will construct colour basis from file "<<m_in_filename<<"."<<endl;
	cout<<"                \nNOTE: amplitude permutations cannot be computed for bases read from files!"<<endl;
      }
      basis = new GenBasis(m_in_filename);
      break;
    }
    case 4: {
      clock_t t(clock());
      if (m_multipletBasis) {
	if (m_verbose>=1) {
	  cout<<"Spectrum::start() Will read in multiplet basis for ";
	  cout<<m_n_g<<" gluons and "<<m_n_qp<<" quark pairs."<<endl;
	}

	MultipletBasis* ortho_basis = new MultipletBasis(m_n_g, m_n_qp);

	if (m_construct_bcm) {
	  ortho_basis->bcm();
	}
	
	basis = ortho_basis;
      }
      else if (m_adjointBasis) {
	if (m_verbose>=1)
	  cout<<"Spectrum::start() Will construct adjoint basis for "<<m_n_g<<" gluons."<<endl;
	basis = new FBasis(m_n_g);
      }
      else if (m_traceBasis) {
	if (m_verbose>=1) {
	  cout<<"Spectrum::start() Will construct trace basis for "<<m_n_g<<" gluons";
	  cout<<" and "<<m_n_qp<<" quark pairs."<<endl;
          if (m_reduceDim) {
            cout<<"Spectrum::start() Dimension of trace basis will be reduced."<<endl;
          }
	}
	basis = new TraceBasis(m_n_g, m_n_qp, m_reduceDim);
      }
      else {
        cerr<<"Spectrum::start() Error! No basis type specified."<<endl;
	exit(EXIT_FAILURE);
      }
      
      t=clock()-t;
      if (m_verbose>=1)
	cout<<"Spectrum::start() Computation time for basis construction: "<<(float)t/CLOCKS_PER_SEC<<"s."<<endl;

      break;
    }
    default: {
      cerr<<"Spectrum::start() Error! Couldn't perform colour calculations!";
      cerr<<"                  Check -h for usage instructions."<<endl;
      exit(EXIT_FAILURE);
    }
  }

  if (m_norm_basis) basis->normalise(m_largeNC);
  if (m_verbose>=2) {
    cout<<endl;
    if (basis->is_normalised()) cout << "Normalised ";
    cout<<"Basis Vectors:"<<endl;
    basis->print();
  }

  clock_t t(clock());
  if (m_verbose>=1)
    cout<<"\nCalculating the soft matrix..."<<endl;
  CMatrix soft_matrix(basis->sm(m_largeNC));
  if (m_verbose>=2) {
    cout<<"\nSoft Matrix:"<<endl;
    soft_matrix.print();
  }
  t=clock()-t;
  if (m_verbose>=1)
    cout<<"Computation time: "<<(float)t/CLOCKS_PER_SEC<<"s."<<endl;

  if (m_calcDet) {
    t=clock();
    double softDet = soft_matrix.det().real();
    t=clock()-t;
    if (m_verbose>=1)
      cout<<"Determinant = "<<softDet<<" (computation time: "<<(float)t/CLOCKS_PER_SEC<<"s)"<<endl;
  }

  t=clock();
  if (m_verbose>=1)
    cout<<"\nCalculating the colour change matrices..."<<endl;
  vector<CMatrix> cc_mats(basis->ccms(m_largeNC));
  if (m_verbose>=2) {
    cout<<"\nColour Change Matrices:"<<endl;
    for (auto& ccm : cc_mats) {
      ccm.print();
      cout<<endl;
    }
  }
  t=clock()-t;
  if (m_verbose>=1)
    cout<<"Computation time: "<<(float)t/CLOCKS_PER_SEC<<"s."<<endl;

  if (m_verbose>=1) {
    cout<<"\nIs colour conserved? ";
    if (m_largeNC) cout << "(in large-NC limit) ";
    if(basis->check_colourcons(m_largeNC)) cout<<"Yes!"<<endl;
    else cout<<"NO!"<<endl;
  }

  if (!m_no_outfile) {
    if (m_verbose>=1)
      cout<<"\nPrinting to file..."<<endl;
    basis->print_to_file("", m_largeNC);
    if (m_verbose>=1)
      cout<<"Done!"<<endl;
  }

  return 0;
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
