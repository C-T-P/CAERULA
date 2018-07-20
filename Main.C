#include<fstream>
#include<iomanip> 
#include<ctime>
#include "colourtools.h"
#include "c_matrix.h"
#include "BASIS.h"
#include "CONTRACT.h"
#include "I3NSERT.h"
#include "Main.h"

int main(int argc, char **argv) {
    clock_t t1, t2;
    
    t1=clock();
    
    // define order in 1/NC to which all terms shall be evaluated
    const int NC_order=INT_MAX;
    cout<<"\nOrder of 1/NC set to ";
    if (NC_order!=INT_MAX) cout << NC_order << "." << endl;
    else cout << "Infinity.\n" << endl;
    
    // read run parameters
    string runpar;
    vector<colour_term> basis;
    process m_process;
    string out_filename;
    for (int i(1);i<argc;i++) runpar+=argv[i];
    size_t f_pos(runpar.find("-f")), nqp_pos(runpar.find("-nqp")), ng_pos(runpar.find("-ng"));
    if (f_pos!=string::npos and nqp_pos==string::npos and ng_pos==string::npos) {
        string filename(runpar.substr(f_pos+2));
        cout<<"Will construct basis from file "<<filename<<"."<<endl;
        
        basis=read_basis(filename, m_process);
        for (size_t lno(1);lno<=m_process.no_of_legs();lno++) out_filename+=m_process.leg(lno).second;
    }
    else if (nqp_pos!=string::npos or ng_pos!=string::npos) {
        string qp_str, g_str;
        if (nqp_pos>ng_pos) {
            g_str=runpar.substr(ng_pos+3,nqp_pos-ng_pos-3);
            qp_str=runpar.substr(nqp_pos+4);
        }
        else {
            qp_str=runpar.substr(nqp_pos+4,ng_pos-nqp_pos-4);
            g_str=runpar.substr(ng_pos+3);
        }
        
        int n_qp(stoi(qp_str)), n_g(stoi(g_str));
        cout<<"Will construct basis for "<<n_qp<<" quark pairs and "<<n_g<<" gluons."<<endl;
        basis=construct_basis(n_qp, n_g, m_process);
        for (int n(0);n<n_qp;n++) out_filename+="qqb";
        for (int n(0);n<n_g;n++) out_filename+="g";
    }
    else {
        cerr<<"Error: specify basis either by filename ( -f ) OR by number of quark pairs and gluons ( -nqp / -ng )."<<endl;
        exit(EXIT_FAILURE);
    }
    
    // normalise basis
    basis=normalise_basis(basis,NC_order);
    unsigned int DIM(basis.size());
    
    // print leg indices and particles
    cout << "\nProcess:\n" << "leg\t" << "particle\t" << "in/out (1/0)" << endl;
    for (size_t i(1);i<=m_process.no_of_legs();i++) cout << m_process.leg(i).first << "\t" << m_process.leg(i).second << "\t\t" << m_process.is_in_leg(i) <<endl;
    
    // print normalised basis vectors
    cout<<"\nNormalised Basis Vectors:"<<endl;
    for (size_t i(0);i<DIM;i++)
        cout<<"b_"<<i+1<<" = "<<basis.at(i).build_string()<<endl;
    
    // print computation time for basis construction
    t1=clock()-t1;
    float runtime=(float)t1/CLOCKS_PER_SEC;
    cout << "\ncomputation time for basis construction: " << (int)runtime/3600 << " h, " << (int)runtime%60/60 << " m, " << (int)runtime%60+runtime-(int)runtime << " s"<< endl;
    
    t2=clock();
    
    // calculate and print soft matrix (to file)
    c_matrix soft_matrix=calc_soft_matrix(basis,NC_order);
    cout << "\nSoft Matrix:" << endl;
    soft_matrix.print();
    ofstream file;
    file.open(out_filename+"_met.dat");
    file<<"# Dimension of Square Matrix\n"<<DIM<<"\n"<<endl;
    file<<"# Soft Matrix / Colour Metric"<<endl;
    for (size_t i(0);i<DIM;i++)
        for (size_t j(0);j<DIM;j++)
            file<<fixed<<setprecision(17)<<soft_matrix[i][j].real()<<" ";
    file.close();

    // calculate and print inverse soft matrix
    c_matrix inv_soft_matrix=calc_inv_soft_matrix(soft_matrix);
    cout << "\nInverse Soft Matrix:" << endl;
    inv_soft_matrix.print();
    
    // calculate and print inverse soft matrix times soft matrix
    c_matrix unit_matrix(DIM);
    cout<<"\nInverse times Soft Matrix"<<endl;
    unit_matrix=soft_matrix*inv_soft_matrix;
    unit_matrix.print();
    
    // calculate and give out colour change matrices for all possible insertions
    bool multiply_with_inv_sm=false;
    vector<c_matrix> colour_change_matrices;
    file.open(out_filename+".dat");
    file<<"# Dimension of Square Matrices\n"<<DIM<<"\n"<<endl;
    cout<<"\nColour Change Matrices (";
    file<<"# Colour Change Matrices (";
    if (!multiply_with_inv_sm) {
        cout<<"not ";
        file<<"not ";
    }
    cout<<"multiplied with inverse colour metric)"<<endl;
    file<<"multiplied with inverse colour metric)"<<endl;
    for (unsigned int lno1(1);lno1<=m_process.no_of_legs();lno1++) {
        for (unsigned int lno2(lno1+1);lno2<=m_process.no_of_legs();lno2++) {
            colour_change_matrices.push_back(calc_colour_change_matrix(basis,soft_matrix,m_process,lno1,lno2,NC_order));
            
            // multiply with inverse soft matrix if wanted
            if (multiply_with_inv_sm) colour_change_matrices.back()=inv_soft_matrix*colour_change_matrices.back();
            
            cout << "C_(" << lno1 << "," << lno2 << ") = " << endl;
            colour_change_matrices.back().print();
            file<<"// C_("<<lno1<<","<<lno2<<") with PIDs "<<m_process.leg(lno1).second<<" and "<<m_process.leg(lno2).second<<endl;
            for (size_t i(0);i<DIM;i++) {
                for (size_t j(0);j<DIM;j++)
                    file<<fixed<<setprecision(17)<<colour_change_matrices.back()[i][j].real()<<" ";
                file<<endl;
            }
            cout<<endl;
        }
    }
    file.close();
    
    // print computation time for colour insertions
    t2=clock()-t2;
    runtime=(float)t2/CLOCKS_PER_SEC;
    cout << "computation time for colour insertions: " << (int)runtime/3600 << " h, " << (int)runtime%60/60 << " m, " << (int)runtime%60+runtime-(int)runtime << " s"<< endl;
    return 0;
} 
