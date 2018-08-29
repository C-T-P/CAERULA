#include "c_basis.h"
#include "gen_basis.h"
#include<fstream>

// member functions of gen_basis class
gen_basis::gen_basis(string filename) {
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
                    cerr << "Error reading in process: leg needs either direction \"in\" or \"out\" , but was given direction \"" << direc << "\"." << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else if (line.at(0)=='b') {
                line.erase(remove_if(line.begin(),line.end(),::isspace),line.end());
                basis_strs.push_back(line.substr(1));
            }
        }
    }
    
    for (auto& str : basis_strs) m_ca_basis.push_back(c_amplitude(str));
    
    m_dim=m_ca_basis.size();
    m_confact=0.;
    m_amp_perms=vector<vector<size_t>>();
    for (size_t i(0); i<m_dim; i++) m_normalisations.push_back(1.);
    
    // initialise matrices
    m_smat=c_matrix(m_dim);
    m_ccmats=vector<c_matrix>();
}
gen_basis::~gen_basis() {
    
}
