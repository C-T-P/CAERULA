#include<string>
#include<fstream>
#include<iostream>
#include<vector>
#include<algorithm>
#include<complex>
#include<cmath>


using namespace std;

class gen_basis {
    
    vector<colour_term> m_gen_basis;
    vector<string> m_in_legs;
    vector<string> m_out_legs;
    
    public:
        gen_basis(string filename = "");
        ~gen_basis();
    
        size_t dim();
        process proc();
        vector<colour_term> ct_basis();
        void print();
};
