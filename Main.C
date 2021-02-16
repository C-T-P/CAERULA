// Copyright (C) 2021 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#include "Spectrum.H"

void print_help_message();

int main(int argc, char **argv) {

  // Create instance of Spectrum.
  SPECTRUM::Spectrum  spectrum;

  // Store filenames.
  string in_filename("");
  string out_filename("");

  // Store number of quark pairs and gluons.
  int nqp(0), ng(0);

  // Store whether we want to print everything to a file.
  bool printfile(true);

  // If no input given, print help meassage.
  if (argc==1) print_help_message();
  
  // Determine options from command line input.
  for(int i(1); i<argc; ++i) {
    // Print help message and exit.
    if (strcmp(argv[i], "-help")==0 or strcmp(argv[i], "-h")==0) {
      print_help_message();
    }

    // Set verbosity.
    else if (strcmp(argv[i], "-verbose")==0 or strcmp(argv[i], "-v")==0) {
      spectrum.set_verbose(stoi(argv[i+1]));
      ++i;
    }

    // We want to evaluate a given colour amplitude.
    else if (strcmp(argv[i], "-evaluate")==0 or strcmp(argv[i], "-e")==0) {
      spectrum.evaluate_ca(argv[i+1]);
      return 0;
    }

    // We want to simplify a given colour amplitude.
    else if (strcmp(argv[i], "-simplify")==0 or strcmp(argv[i], "-s")==0) {
      spectrum.simplify_ca(argv[i+1]);
      return 0;
    }

    // We want to do a full colour computation for process read from file.
    else if (strcmp(argv[i], "-f")==0 or strcmp(argv[i], "-file")==0) {
      in_filename = argv[i+1];
      ++i;
    }

    // We want to do a full colour computation for process specified by partons.
    else if (strcmp(argv[i],"-ng")==0) {
      ng = stoi(argv[i+1]);
      spectrum.set_ng(ng);
      ++i;
    }
    else if (strcmp(argv[i],"-nqp")==0) {
      nqp = stoi(argv[i+1]);
      spectrum.set_nqp(nqp);
      ++i;
    }

    // We want to reduce the dimension of the trace basis by summing conjugates.
    else if (strcmp(argv[i],"-reduceDim")==0) {
      if (spectrum.get_use_trace_basis()) spectrum.set_reduce_trace_dim(true);
      else {
        cerr << "Error: Can only reduce dimension of trace bases."
             << "Option -reduceDim will be ignored." << endl;
      }
    }

    // Use the adjoint basis.
    else if (strcmp(argv[i],"-adj")==0) {
      if (!spectrum.get_use_multiplet_basis()) {
        spectrum.set_use_adjoint_basis(true);
        spectrum.set_use_trace_basis(false);
      }
      else {
        cerr << "Error: Two competing basis types specified."
             << "Unclear which one to use." << endl;
        exit(EXIT_FAILURE);
      }
    }

    // Use the multiplet basis.
    else if (strcmp(argv[i],"-multiplet")==0) {
      if (!spectrum.get_use_adjoint_basis()) {
        spectrum.set_use_multiplet_basis(true);
        spectrum.set_use_trace_basis(false);
      }
      else {
        cerr << "Error: Two competing basis types specified."
             << "Unclear which one to use." << endl;
        exit(EXIT_FAILURE);
      }
    }

    // Use large-NC limit.
    else if (strcmp(argv[i],"-largeNC")==0) {
      spectrum.set_largeNC(true);
    }   

    // Do not normalise basis.
    else if (strcmp(argv[i],"-nonorm")==0) {
      spectrum.set_normalise_basis(false);
    }

    // Construct basis change matrix to trace basis.
    else if (strcmp(argv[i],"-bcm")==0) {
      spectrum.set_construct_bcm(true);
    }

    // Calculate determinant of soft matrix.
    else if (strcmp(argv[i],"-det")==0) {
      spectrum.set_calc_det(true);
    }

    // Do not save calculation to file.
    else if (strcmp(argv[i],"-nooutfile")==0) {
      printfile = false;
    }

    else {
      cout << "Unrecognised option " << argv[i] << ". Will skip." << endl;
    }
  }

  // Print current settings of parameters
  spectrum.print_settings();

  // If filename is specified read basis from there.
  // Otherwise construct basis from number of quark pairs and gluons
  if (in_filename != "") {
    if(!spectrum.read_basis(in_filename)) {
      cerr << "ERROR: Could not read basis from file " 
           << in_filename << "." << endl;
      exit(EXIT_FAILURE);
    }
  } else {
    if (!spectrum.construct_basis(nqp, ng)) {
      cerr << "ERROR: Could not construct basis for " 
           << nqp << " quark pairs and " << ng << " gluons." << endl;
      exit(EXIT_FAILURE);
    }
  }
  
  // Print basis
  spectrum.print_basis();

  // Calculate soft matrix
  SPECTRUM::CMatrix softmatrix = spectrum.calculate_soft_matrix();

  // Calculate colour correlators
  vector<SPECTRUM::CMatrix> ccmatrices = spectrum.calculate_colour_correlators();

  // Save everything to file
  if (printfile) spectrum.save_to_file(out_filename);

  return 0;
}


// Method to print a help message to cout
void print_help_message() {
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
