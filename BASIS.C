#include<fstream>
#include<iomanip> 
#include "colourtools.h"
#include "c_matrix.h"
#include "Main.h"
#include "CONTRACT.h"
#include "FCONTRACT.h"
#include "I3NSERT.h"
#include "BASIS.h"

/*
===========================================================================
========================= Basis from external file ========================
===========================================================================
*/

vector<colour_term> read_basis(string filename, process& m_process) {
    vector<string> basis_strs;
    vector<colour_term> basis;
    
    read_in_process(filename, m_process, basis_strs);
    
    for (size_t i(0);i<basis_strs.size();i++)
        basis.push_back(decompose_terms(basis_strs.at(i),m_process));
    
    return basis;
}

// read in process file, store process information in process variable and basis vectors as strings
void read_in_process(string filename, process& m_process, vector<string>& basis_strs) {
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
}

// decompose input into colour term according to the indices given in the process specification
colour_term decompose_terms(string& input, process m_process) {
    // create mock process in case none was given
    if (m_process.no_of_legs()==0) {
        cout<<"\nNOTE: all Kronecker deltas k_[.,.] will be treated as having quark indices. Use c_[2.,0.]*t_[.,101,102]*t_[.,102,101] for gluonic deltas.\n"<<endl;
        for (int i(0);i<200;i++) m_process.add_out_leg("q");
    }
    
    colour_term expr;
    three_ind symmetric;
    three_ind antisymmetric;
    three_ind fundamental;
    two_ind kronecker;
    int NC_order_c;
    int NC_order_r;
    complex<double> prefactor = 1.;
    for (size_t i(0), mpos(input.find('+'));mpos!=string::npos || input.length()>0;mpos=input.find('+')) {
        ++i;
        // decompose input into summands
        string summand;
        if (mpos==string::npos) {
            summand=input;
            input="";
        }
        else {
            summand=input.substr(0,mpos);
            input=input.substr(mpos+1);
        }
        // decompose summand into factors
        kronecker.clear_indices();
        symmetric.clear_indices();
        antisymmetric.clear_indices();
        fundamental.clear_indices();
        prefactor=1.;
        NC_order_c=0;
        NC_order_r=0;
        for (size_t j(0), mpos(summand.find('*'));mpos!=string::npos || summand.length()>0;mpos=summand.find('*')) {
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
            if (factor.find("c_[")==0 && factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==string::npos || factor.find(',',cpos+1)!=string::npos) {
                    cerr << "Invalid prefactor." << endl;
                    exit(EXIT_FAILURE);
                }
                complex<double> prfct(0.);
                size_t spos(factor.substr(3,cpos-3).find("1/NC"));
                if (spos==string::npos) prfct+=stod(factor.substr(3,cpos-3));
                else {
                    size_t ppos=factor.substr(3,cpos-3).find('^');
                    double exponent(1);
                    if (ppos!=string::npos) exponent=stod(factor.substr(3+ppos+1,cpos-ppos-2));
                    NC_order_r+=exponent;
                    prfct+=1./pow(NC,exponent);
                }
                spos=factor.substr(cpos+1,factor.length()-cpos-2).find("1/NC");
                if (spos==string::npos) prfct+=stod(factor.substr(cpos+1,factor.length()-cpos-2))*complex<double>(0.,1.);
                else {
                    size_t ppos=factor.substr(cpos+1,factor.length()-cpos-2).find('^');
                    double exponent(1);
                    if (ppos!=string::npos) exponent=stod(factor.substr(cpos+ppos+2,factor.length()-ppos-3));
                    NC_order_c+=exponent;
                    prfct+=(double)1./pow(NC,exponent)*complex<double>(0.,1.);
                }
                // only same orders in 1/NC in real and imaginary part are supported
                if (NC_order_c!=NC_order_r and NC_order_c!=0 and NC_order_r!=0) {
                    cerr << "Error: Expected equal orders in 1/NC for real and imaginary part of prefactor, but real part has order " << NC_order_r << " and imaginary part has order " << NC_order_c << ".\n Specify real and imaginary part seperately to avoid this error." << endl;
                    exit(EXIT_FAILURE);
                }
                prefactor*=prfct;
            }
            else if(factor.find("f_[")==0 && factor[factor.length()-1]==']') {
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
                antisymmetric.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if(factor.find("d_[")==0 && factor[factor.length()-1]==']') {
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
                symmetric.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if(factor.find("t_[")==0 && factor[factor.length()-1]==']') {
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
                fundamental.set_indices(stoi(factor.substr(3,c1pos-3)), stoi(factor.substr(c1pos+1,c2pos-c1pos-1)), stoi(factor.substr(c2pos+1,factor.length()-c2pos-2)));
            }
            else if (factor.find("k_[")==0 && factor[factor.length()-1]==']') {
                size_t cpos(factor.find(','));
                if (cpos==string::npos || factor.find(',',cpos+1)!=string::npos) {
                    cerr << "Invalid number of indices for k." << endl;
                    exit(EXIT_FAILURE);
                }
                int ind_i=stoi(factor.substr(3,cpos-3)), ind_j=stoi(factor.substr(cpos+1,factor.length()-cpos-2));
                bool gluonic(false);
                if (m_process.leg((unsigned int)ind_i).second=="g" or m_process.leg((unsigned int)ind_j).second=="g") gluonic=true;
                kronecker.set_indices(ind_i,ind_j,gluonic);
            }
            else cerr << "Invalid input." << endl;
        }
        expr.sym.push_back(symmetric);
        expr.asym.push_back(antisymmetric);
        expr.fund.push_back(fundamental);
        expr.kron.push_back(kronecker);
        expr.pref.push_back(prefactor);
        expr.NC_ctr.push_back(NC_order_r);
    }
    return expr;
}


/*
===========================================================================
================= Construct Basis from no of g and q pairs ================
===========================================================================
*/

vector<colour_term> construct_basis(int n_qp, int n_g, process& m_process, vector<vector<int>>& amp_perms) {
    vector<colour_term> basis;
    
    // setup process with ordering q...qqb...qbg...g and all legs outgoing
    vector<int> q_indices, qb_indices, g_indices;
    for (int n(1);n<=n_qp;n++) {
        q_indices.push_back(n);
        m_process.add_out_leg("q");
    }
    for (int n(n_qp+1);n<=2*n_qp;n++) {
        qb_indices.push_back(n);
        m_process.add_out_leg("qb");
    }
    for (int n(2*n_qp+1);n<=2*n_qp+n_g;n++) {
        qb_indices.push_back(n);
        m_process.add_out_leg("g");
    }
    
    if (n_qp==0 and n_g==0) {
        cerr<<"No basis to build out of 0 gluons and 0 quark pairs."<<endl;
        exit(EXIT_FAILURE);
    }
    else if (n_qp==0) basis=build_g_basis(n_g, amp_perms);
    else if (n_g==0) basis=build_q_basis(n_qp, amp_perms);
    else {
//         basis=build_qqbg_basis(n_qp, n_g);
        cerr<<"Automatic basis construction for mixed processes not yet supported."<<endl;
        exit(EXIT_FAILURE);
    }
    
    return basis;
}

// TODO construct trace basis for mixed processes

// vector<colour_term> build_qqbg_basis(int n_qp, int n_g) {
//     vector<colour_term> basis;
//     vector<int> q_indices;
//     vector<int> qb_indices;
//     vector<int> g_indices;
//     
//     everything defined as outgoing particles 
//     means: qb is an initial state q and vice versa
//     for (int ind(1);ind<=n_qp;ind++) q_indices.push_back(ind);
//     for (int ind(n_qp+1);ind<=2*n_qp;ind++) qb_indices.push_back(ind);
//     for (int ind(2*n_qp+1);ind<=2*n_qp+n_g;ind++) g_indices.push_back(ind);
//     
//     cout<<"q indices"<<endl;
//     for (const auto& i : q_indices) cout<<i<<" ";
//     cout<<"\nqb indices"<<endl;
//     for (const auto& i : qb_indices) cout<<i<<" ";
//     cout<<"\ng indices"<<endl;
//     for (const auto& i : g_indices) cout<<i<<" ";
//     cout<<endl;
//     
//     cout<<"Quark pairs connected to gluons:"<<endl;
//     vector<vector<int>> qqb_ind_combos(get_q_ind_combinations(q_indices, qb_indices));
//     vector<int> qp_ind;
//     for (int con_g(n_g);con_g>0;con_g--) {
//         vector<vector<int>> g_ind_combos(get_g_ind_combinations(g_indices.size(),con_g, g_indices));
//         for (const auto& g_ind : g_ind_combos) {
//             vector<int> g_res_ind;
//             for (const auto& ind : g_indices) if (find(g_ind.begin(),g_ind.end(),ind)==g_ind.end()) g_res_ind.push_back(ind);
//             
//             for (const auto& i : g_ind) cout<<i<<" ";
//             cout<<endl;
//             for (const auto& i : g_res_ind) cout<<i<<" ";
//             cout<<endl;
//         }
//         
//         for (const auto& qqb_i : qqb_ind_combos) {
//             
//             
//             go through quark pairs
//             for (size_t qp_no(0);qp_no<qqb_i.size();qp_no+=2) {
//                 qp_ind.clear();
//                 qp_ind.push_back(qqb_i[qp_no]);
//                 qp_ind.push_back(qqb_i[qp_no+1]);
//                 trace_qngqb(qp_ind, g_indices, con_g);
//             }
//         }
//     }
//     
//     return basis;
// }

// void trace_qngqb(vector<int> qqb_ind, vector<int> g_ind, int con_g) {
//     if (qqb_ind.size()!=2) {
//         cerr<<"Error: Wrong number of quark indices. n gluons can only be connected to exactly 1 quark pair."<<endl;
//         return;
//     }
//     vector<vector<int>> g_ind_combos(get_g_ind_combinations(g_ind.size(),con_g, g_ind)), arr_qgqb_ind;
//     for (const auto& g_i : g_ind_combos) {
//         arr_qgqb_ind.push_back(connect_qngqb_ind(qqb_ind, g_i));
//         
//     }
// }
// 
// vector<int> connect_qngqb_ind(vector<int> qqb_ind, vector<int> g_ind) {
//     if (g_ind.size()==0) return qqb_ind;
//     else if (qqb_ind.size()==0) return g_ind;
//     else {
//         vector<int> con_qngqb_ind(qqb_ind[0]);
//         for (const auto& ind : g_ind) con_qngqb_ind.push_back(ind);
//         con_qngqb_ind.push_back(qqb_ind[0]);
//         return con_qngqb_ind;
//     }
// }

vector<colour_term> build_q_basis(int n_qp, vector<vector<int>>& amp_perms) {
    vector<int> q_indices;
    vector<int> qb_indices;
    
    // everything defined as outgoing particles 
    for (int ind(1);ind<=n_qp;ind++) q_indices.push_back(ind);
    for (int ind(n_qp+1);ind<=2*n_qp;ind++) qb_indices.push_back(ind);
    
    // get arranged quark indices
    vector<vector<int>> qqb_ind_combos(get_q_ind_combinations(q_indices,qb_indices));
    amp_perms=qqb_ind_combos;
    
    return colourflow_q(qqb_ind_combos);
}

vector<vector<int>> get_q_ind_combinations(vector<int> q_inds, vector<int> qb_inds) {
    vector<vector<int>> q_ind_perms, qqb_ind_combos;
    vector<int> tmp;
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

vector<colour_term> colourflow_q(vector<vector<int>> qqb_ind_combos) {
    vector<colour_term> basis;
    colour_term basis_el;
    three_ind symmetric;
    symmetric.clear_indices();
    three_ind antisymmetric;
    antisymmetric.clear_indices();
    three_ind fundamental;
    fundamental.clear_indices();
    two_ind kronecker;
    
    for (auto& qqb_i : qqb_ind_combos) {
        kronecker.clear_indices();
        basis_el.delete_all_terms();
        while (qqb_i.size()>0) {
            kronecker.set_indices(qqb_i[0],qqb_i[1],false);
            qqb_i.erase(qqb_i.begin(),qqb_i.begin()+2);
        }
        basis_el.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
        basis.push_back(basis_el);
    }
    
    return basis;
}

vector<colour_term> build_g_basis(int n_g, vector<vector<int>>& amp_perms) {
    vector<colour_term> basis;
    colour_term basis_el;
    vector<int> indices;
    for (int i(1);i<=n_g;i++) indices.push_back(i);

    // get gluon partitions
    vector<vector<int>> g_partitions(get_g_partitions(n_g));

    vector<vector<int>> arr_g_ind;
    for (size_t i(0);i<g_partitions.size();i++) {
        arr_g_ind=get_arranged_g_ind_for_part(indices,g_partitions[i]);
        
        // save permutations of connected amplitudes to file
        if (i==0) amp_perms=arr_g_ind;
        
        // construct trace basis elements for all arranged index combinations
        for (const auto& inds : arr_g_ind) {
            basis_el.delete_all_terms();
            int place(0);
            for (const auto& p_size : g_partitions[i]) {
                vector<int> sub_ind(inds.begin()+place,inds.begin()+place+p_size);
                colour_term ct_factor(trace_connected_g(sub_ind));
                if (basis_el.no_of_terms()==0) basis_el.add_colour_term(ct_factor);
                else basis_el=basis_el.multiply(ct_factor);
                place+=p_size;
            }
            basis.push_back(basis_el);
        }
        
        // TODO implement adjoint basis construction e.g. via switch
        
        // print groupings
//         cout<<"\n( ";
//         for (const auto& g : g_partitions[i]) cout<<g<<" ";
//         cout<<"):"<<endl;
//         for (const auto& s : arr_g_ind) {
//             for (const auto& el : s) cout<<el<<" ";
//             cout<<endl;
//         }
//         cout<<"DIM = "<<arr_g_ind.size()<<endl;
    }
    return basis;
}

vector<vector<int>> get_arranged_g_ind_for_part(vector<int> ind, vector<int> g_partition) {
    vector<vector<int>> ind_combos, ind_combo_perms, tmp, arr_g_ind;
    ind_combos=arrange_g_ind(ind,g_partition,0);
    
    // TODO: Review this part!!! Should be doable in an easier way
    for (const auto& inds : ind_combos) {
        int place(0);
        ind_combo_perms.clear();
        for (const auto& p_size : g_partition) {
            vector<int> sub_ind(inds.begin()+place,inds.begin()+place+p_size);
            if (p_size>3) {
                vector<vector<int>> permutations;
                permutations.push_back(sub_ind);
                while (next_permutation(sub_ind.begin()+1,sub_ind.end())) permutations.push_back(sub_ind);
                
                // get physical permutations only: remove reflections
                for (size_t p_it(0);p_it<permutations.size();p_it++) {
                    vector<int> p_ind(permutations[p_it]);
                    vector<int> ind_refl({p_ind[0]});
                    for (size_t it(p_ind.size()-1);it>0;it--) ind_refl.push_back(p_ind[it]);
                    
                    vector<vector<int>>::iterator it(find(permutations.begin(),permutations.end(),ind_refl));
                    if (it!=permutations.end()) {
                        permutations.erase(it);
                        p_it--;
                    }
                }

                tmp=ind_combo_perms;
                ind_combo_perms.clear();
                for (const auto& perm_ind : permutations) {
                    if (tmp.size()!=0) {
                        for (auto& tmp_ind : tmp) {
                            vector<int> t_ind=tmp_ind;
                            t_ind.insert(t_ind.end(),perm_ind.begin(),perm_ind.end());
                            ind_combo_perms.push_back(t_ind);
                        }
                    }
                    else ind_combo_perms.push_back(perm_ind);
                }
            }
            else {
                tmp=ind_combo_perms;
                ind_combo_perms.clear();
                if (tmp.size()!=0) {
                    for (auto& tmp_ind : tmp) {
                        tmp_ind.insert(tmp_ind.end(),sub_ind.begin(),sub_ind.end());
                        ind_combo_perms.push_back(tmp_ind);
                    }
                }
                else ind_combo_perms.push_back(sub_ind);
            }
            place+=p_size;
        }
        arr_g_ind.insert(arr_g_ind.end(),ind_combo_perms.begin(),ind_combo_perms.end());
    }
    return arr_g_ind;
}

vector<vector<int>> get_g_ind_combinations(int N, int K, vector<int> inds) {
    vector<vector<int>> ind_combos;
    vector<int> tmp;
    
    string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's

    // permute bitmask
    do {
        tmp.clear();
        for (int i(0);i<N;i++) if (bitmask[i]) tmp.push_back(inds[i]);
        ind_combos.push_back(tmp);
    } while (prev_permutation(bitmask.begin(), bitmask.end()));
    return ind_combos;
}

vector<vector<int>> arrange_g_ind(vector<int> ind, vector<int> g_partition, size_t level) {
    int n_g(ind.size());
    vector<vector<int>> g_comb;
    
    if (n_g==2 or level>=g_partition.size()) g_comb.push_back(ind);
    else {
        g_comb=get_g_ind_combinations(n_g,g_partition[level],ind);
        
        vector<int> s_ind;
        size_t comb_no(g_comb.size());
        for (size_t it(0);it<comb_no;it++) {
            vector<int> s(g_comb[0]);
            g_comb.erase(g_comb.begin());
            s_ind=ind;
            for (const auto& c : s) s_ind.erase(remove(s_ind.begin(), s_ind.end(), c), s_ind.end()); 
            vector<vector<int>> n_g_comb=arrange_g_ind(s_ind, g_partition, level+1);
            vector<int> tmp;
            for (const auto& new_s : n_g_comb) {
                if (level+1>g_partition.size() or g_partition[level]!=g_partition[level+1] or (new_s.size()>0 and s[0]<new_s[0])) {
                    tmp=s;
                    tmp.insert(tmp.end(),new_s.begin(),new_s.end());
                    g_comb.push_back(tmp);
                }
            }
        }
    }
    return g_comb;
}

colour_term trace_connected_g(vector<int> indices) {
    colour_term ct;
    three_ind symmetric;
    three_ind antisymmetric;
    three_ind fundamental;
    two_ind kronecker;
    
    int n_g(indices.size());
    
    if (n_g==2) {
        kronecker.set_indices(indices[0],indices[1],true);
        ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
    }
    else {
        vector<int> ind_refl;
        double sign(1.);
        if (n_g%2!=0) sign=-1.;
        ct.delete_all_terms();
        ind_refl.clear();
        ct.add_colour_term(write_tr_basis_el(indices));
        ind_refl.push_back(indices[0]);
        for (size_t it(1);it<indices.size();it++) ind_refl.push_back(indices[indices.size()-it]);
        ct.add_colour_term(write_tr_basis_el(ind_refl));
        ct.pref[1]*=sign;
    }
    return ct;
}

colour_term write_tr_basis_el(vector<int> ind) {
    colour_term ct;
    three_ind symmetric;
    three_ind antisymmetric;
    three_ind fundamental;
    two_ind kronecker;
    
    int start_ind(101), incr(0);
    for (int i(0);(unsigned)i<ind.size();i++) {
        int c_ind;
        if ((unsigned)i==ind.size()-1) c_ind=start_ind;
        else c_ind=start_ind+incr+1;
        fundamental.set_indices(ind[i],start_ind+incr,c_ind);
        incr++;
    }
    ct.add_term(symmetric,antisymmetric,fundamental,kronecker,complex<double>(1.,0.),0);
    return ct;
}

vector<vector<int>> get_g_partitions(int n) {
    vector<vector<int>> groupings;
    vector<int> trivial_grouping(n);
    int_partitions(n, trivial_grouping, groupings, 0);
    
    return groupings;
}

void int_partitions(int n, vector<int>& v, vector<vector<int>>& r_v, int level) {
    int first;
    vector<int> tmp;
    bool vanishes(false);

    if(n<1) return ;
    v[level]=n;
    for(int i=0;i<=level;i++) {
        tmp.push_back(v[i]);
        
        // groupings with only one gluon vanish due to Tr(T_a)=0
        if (v[i]==1) vanishes=true;
    }
    if (!vanishes) r_v.push_back(tmp);

    first=(level==0) ? 1 : v[level-1];

    for(int i=first;i<=n/2;i++){
        v[level]=i; /* replace last */
        int_partitions(n-i, v, r_v, level+1);
    }
}

// normalise Basis
vector<colour_term> normalise_basis(vector<colour_term> basis, int NC_order, vector<complex<double>>& normalisations) {
    unsigned int DIM(basis.size());
    vector<colour_term> n_basis;
    colour_term ct;
    for (size_t i(0);i<DIM;i++) {
        ct.delete_all_terms();
        ct=basis[i].scprod(basis[i]);
        normalisations.at(i)=sqrt(fast_evaluate_colour_term_to_order(ct,NC_order));
        ct=basis[i];
        for (size_t t_it(0);t_it<ct.no_of_terms();t_it++) ct.pref[t_it]*=1./normalisations.at(i);
        n_basis.push_back(ct);
    }
    return n_basis;
}



/*
===========================================================================
==================== Computation of Soft & Hard Matrix ====================
===========================================================================
 */

// calculate soft matrix
c_matrix calc_soft_matrix(vector<colour_term> basis, int NC_order) {
    unsigned int DIM(basis.size());
    c_matrix soft_matrix(DIM);
    colour_term ct;
    for (size_t i(0);i<DIM;i++) {
        for (size_t j(0);j<DIM;j++) {
            if(i<=j) {
                ct.delete_all_terms();
                ct=basis[i].scprod(basis[j]);
//                soft_matrix[i][j]=evaluate_colour_term_to_order(ct,NC_order);
                soft_matrix[i][j]=fast_evaluate_colour_term_to_order(ct,NC_order);
            }
            else soft_matrix[i][j]=conj(soft_matrix[j][i]);
        }
    }
    return soft_matrix;
} 

// calculate inverse soft matrix
c_matrix calc_inv_soft_matrix(c_matrix soft_matrix) {
    unsigned int DIM(soft_matrix.dim());
    c_matrix inv_soft_matrix(DIM);
    
    gsl_matrix_complex *S=gsl_matrix_complex_alloc(DIM,DIM);
    gsl_matrix_complex *IS=gsl_matrix_complex_alloc(DIM,DIM);
    gsl_permutation *p=gsl_permutation_alloc(DIM);
    int n;
    for (size_t j(0);j<DIM;j++)
        for (size_t k(0);k<DIM;k++) gsl_matrix_complex_set(S,j,k,gsl_complex_rect(soft_matrix[j][k].real(),soft_matrix[j][k].imag())); 
    gsl_linalg_complex_LU_decomp(S,p,&n);
    
    gsl_complex gsl_det(gsl_linalg_complex_LU_det(S,n));
    complex<double> det(complex<double>(GSL_REAL(gsl_det),GSL_IMAG(gsl_det)));
    cout<<"\ndet(S) = "<<det<<endl;
    
    gsl_linalg_complex_LU_invert(S,p,IS);
    
    for (size_t i(0);i<DIM;i++)
        for (size_t j(0);j<DIM;j++) inv_soft_matrix[i][j]=complex<double>(GSL_REAL(gsl_matrix_complex_get(IS,i,j)),GSL_IMAG(gsl_matrix_complex_get(IS,i,j)));
    
    return inv_soft_matrix;
} 

// calculate colour change matrix for insertion between leg lno1 and lno2
c_matrix calc_colour_change_matrix(vector<colour_term> basis, c_matrix soft_matrix, process m_process, unsigned int lno1, unsigned int lno2, int NC_order) {
    colour_term insertion_op=construct_insertion_op(m_process,lno1,lno2);
    colour_term ct;
    unsigned int DIM(basis.size());
    c_matrix ccm(DIM);
    for (size_t i(0);i<DIM;i++) {
        for (size_t j(0);j<DIM;j++) {
            cout << "\rC_(" << lno1 << "," << lno2 << ")[" << i << "][" << j << "] being calculated..." << flush;
            ct=basis[i].scprod(insertion_op.multiply(make_internal(m_process,basis[j])));
//            ccm[i][j]=evaluate_colour_term_to_order(ct,NC_order);
            ccm[i][j]=fast_evaluate_colour_term_to_order(ct,NC_order);
            if (isnan(ccm[i][j].real())) cerr << "Error: C_(" << lno1 << "," << lno2 << ")[" << i << "][" << j << "] = " <<ct.build_string() << "\n" << endl;
        }
    }
    cout << endl;
    return ccm;
}

// shift external to internal indices
colour_term make_internal(process m_process, colour_term expr) {
    for (unsigned int lno(1);lno<=m_process.no_of_legs();lno++) {
        for (size_t t_it(0);t_it<expr.no_of_terms();t_it++) {
            expr.sym[t_it].find_and_rep_indices(lno,lno+2000);
            expr.asym[t_it].find_and_rep_indices(lno,lno+2000);
            expr.fund[t_it].find_and_rep_indices(lno,lno+2000);
            expr.kron[t_it].find_and_rep_indices(lno,lno+2000);
        }
    }
    return expr;
}
