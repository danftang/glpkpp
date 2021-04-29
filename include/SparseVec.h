//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_SPARSEVEC_H
#define GLPKTEST_SPARSEVEC_H

class SparseVec {
public:
        int *indices;   // one-based array of indices of non-zero elements (zero'th element stores number of non-zeros)
        double *values; // one-based array of values of non-zero elements (zero'th element stores dimension of this vector)


    SparseVec(int dimension, int capacity=-1);
    SparseVec(int dimension, const std::map<int,double> &);
    SparseVec(SparseVec &&rvalue) { // move semantics
        indices = rvalue.indices;
        values = rvalue.values;
        rvalue.indices = NULL;
    }
    SparseVec(const SparseVec &rvalue): SparseVec(rvalue.dimension(), rvalue.sparseSize()) { // copy semantics
        std::copy(rvalue.indices, rvalue.indices + rvalue.sparseSize(), indices);
        std::copy(rvalue.values, rvalue.values + rvalue.sparseSize(), values);
    }

    ~SparseVec();

//    double operator [](int i) const;
//    double &operator [](int i);

    int dimension() const { return (int)values[0]; }
    int sparseSize() const { return indices[0]; }
    int &sparseSize() { return indices[0]; }

    void toDense(double *denseVec) const;
//    std::pair<int,double> entry(int k) const;
    void add(int i, double v);
    void clear();
    void setCapacity(int size);
    void setDimension(int dim) { values[0] = (double)dim; }
//    double dotProd(double *dense) const;

    int maxNonZeroIndex();

protected:
};

std::ostream &operator <<(std::ostream &out, const SparseVec &sVector);

#endif //GLPKTEST_SPARSEVEC_H
