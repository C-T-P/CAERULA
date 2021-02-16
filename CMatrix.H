// Copyright (C) 2021 Christian T Preuss
// This file is part of Spectrum.
//
// Spectrum is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// any later version.

#ifndef CMATRIX_H
#define CMATRIX_H

#include "Colourtools.H"

using namespace std;

namespace SPECTRUM {

typedef vector<vector<complex<double>>> matrix;

class CMatrix {
  
  vector<vector<complex<double>>> m_mat;

 public:
  CMatrix() {
    m_mat=vector<vector<complex<double>>>();
  }
  CMatrix(size_t dim) {
    for (size_t i(0);i<dim;i++)
      m_mat.push_back(vector<complex<double>>(dim,0.));
  }
  size_t dim() {
    return m_mat.size();
  }
  complex<double> det() {
    size_t n = dim();
    return determinant(*this, n);
  }
  vector<complex<double>>& operator[](size_t i) {
    return m_mat[i];
  }
  vector<complex<double>> operator[](size_t i) const {
    return m_mat[i];
  }
  CMatrix operator*(CMatrix m_mat2) {
    size_t dim(m_mat2.dim());
    CMatrix m_matr(dim);
    if (dim==m_mat.size()) {
      for (size_t i(0);i<dim;++i)
        for (size_t j(0);j<dim;++j)
          for (size_t k(0);k<dim;++k)
            m_matr[i][j]+=m_mat[i][k]*m_mat2[k][j];
    }
    else {
      cout<<"Error mutliplying square matrices "<<&m_mat<<" and "<< &m_mat2
          <<": mismatch of dimensions: "<<m_mat.size()<<" vs. "<<dim<<endl;
      exit(EXIT_FAILURE);
    }
    return m_matr;
  }
  CMatrix operator+(CMatrix m_mat2) {
    size_t dim(m_mat2.dim());
    CMatrix m_matr(dim);
    if (dim==m_mat.size()) {
      for (size_t i(0);i<dim;++i)
        for (size_t j(0);j<dim;++j)
          m_matr[i][j]=m_mat[i][j]+m_mat2[i][j];
    }
    else {
      cout<<"Error adding square matrices "<<&m_mat<<" and "<< &m_mat2
          <<": mismatch of dimensions: "<<m_mat.size()<<" vs. "<<dim<<endl;
      exit(EXIT_FAILURE);
    }
    return m_matr;
  }
  void operator+=(CMatrix m_mat2) {
    size_t dim(m_mat2.dim());
    if (dim==m_mat.size()) {
      for (size_t i(0);i<dim;++i)
        for (size_t j(0);j<dim;j++)
          m_mat[i][j]+=m_mat2[i][j];
    }
    else {
      cout<<"Error adding square matrices "<<&m_mat<<" and "<< &m_mat2
          <<": mismatch of dimensions: "<<m_mat.size()<<" vs. "<<dim<<endl;
      exit(EXIT_FAILURE);
    }
  }
  void print() {
    size_t dim(m_mat.size());
    for (size_t i(0);i<dim;i++) {
      cout<<"  ";
      for (size_t j(0);j<dim;j++)
        cout<<fixed<<setprecision(8)<<m_mat[i][j].real()<<" ";
      cout<<endl;
    }
  }

 private:
  
  // Helper function
  complex<double> determinant(CMatrix mat, size_t n) {
    if (n == 1) {
      return mat[0][0];
    }
    else if (n == 2) {
      return (mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0]);
    }
    else {
      CMatrix temp(n);
      complex<double> det(0.);
      for(size_t p(0); p<n; ++p) {
        size_t h(0), k(0);
        for(size_t i(1); i<n; ++i) {
          for(size_t j(0); j<n; ++j) {
            if(j == p) {
              continue;
            }
            temp[h][k] = mat[i][j];
            ++k;
            if(k == n-1) {
              ++h;
              k = 0;
            }
          }
        }
        det += mat[0][p]*pow(-1,p)*determinant(temp,n-1);
      }
      return det;
    }
  }
};

}

#endif
