//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_SPARSEVEC_H
#define GLPKTEST_SPARSEVEC_H

class SparseVec {
public:
    int n;          // dimension of vector
    int nnz;        // number of non-zero elements
    int *indices;   // one-based array of indices of non-zero elements
    double *values; // one-based array of values of non-zero elements

    SparseVec(int dimension, int capacity=-1);
    SparseVec(int dimension, const std::map<int,double> &);
    ~SparseVec();

//    double operator [](int i) const;
//    double &operator [](int i);

    void toDense(double *denseVec) const;
    void entry(int k, std::pair<int,double> &);
    void add(int i, double v);
    void clear();
    void setCapacity(int size);
//    double dotProd(double *dense) const;

    int maxNonZeroIndex();

protected:
};

std::ostream &operator <<(std::ostream &out, const SparseVec &sVector);

#endif //GLPKTEST_SPARSEVEC_H
