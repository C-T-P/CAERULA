// Copyright (C) 2021 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#include<fstream>
#include "GenBasis.H"

namespace SPECTRUM {

//*****************************************************************************
//
// Member functions of class GenBasis.
//
//*****************************************************************************

GenBasis::GenBasis(string filename) {
  // set basis type
  m_btype=0;
  
  vector<string> basis_strs;
    
  // read in process file, store process information and basis vectors as strings
  ifstream fin(filename);
  if (!fin) {
    cerr << "Error reading in process: file " << filename << " could not be opened." << endl;
    exit(EXIT_FAILURE);
  }
  else {
      string line;
      while (getline(fin,line)) {
	if (line.at(0)=='l') {
	  string direc("");
	  direc=line.substr(2,line.find("\t",2)-1);
	  direc.erase(remove_if(direc.begin(),direc.end(),::isspace),direc.end());
	  line.erase(0,line.find("\t",2));
	  line.erase(remove_if(line.begin(),line.end(),::isspace),line.end());
	  if (direc=="in") m_process.add_in_leg(line);
	  else if (direc=="out") m_process.add_out_leg(line);
	  else {
	    cerr << "Error reading in process:"
                 << " leg needs either direction \"in\" or \"out\" , but was given \""
                 << direc << "\"." << endl;
	    exit(EXIT_FAILURE);
	  }
	}
	else if (line.at(0)=='b') {
	  line.erase(remove_if(line.begin(),line.end(),::isspace),line.end());
	  basis_strs.push_back(line.substr(1));
	}
      }
  }
    
  for (auto& str : basis_strs) m_ca_basis.push_back(CAmplitude(str));
  
  m_dim=m_ca_basis.size();
  m_confact=0.;
  m_amp_perms=vector<vector<size_t>>();
  for (size_t i(0); i<m_dim; i++) m_normalisations.push_back(1.);
  
  // initialise matrices
  m_smat=CMatrix(m_dim);
  m_ccmats=vector<CMatrix>();
}

GenBasis::~GenBasis() {;}

}
