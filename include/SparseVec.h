//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_SPARSEVEC_H
#define GLPKTEST_SPARSEVEC_H

class SparseVec {
public:
    int *indices;   // one-based array of indices of non-zero elements (zero'th element stores number of non-zeros)
    double *values; // one-based array of values of non-zero elements (zero'th element stores dimension of this vector)

    struct Header {
        int dimension;
        int capacity;
    };

    SparseVec(int dimension, int capacity=-1);
    SparseVec(int dimension, const std::map<int,double> &);
    SparseVec(SparseVec &&rvalue) { // move semantics
        moveFrom(rvalue);
    }
    SparseVec(const SparseVec &lvalue): SparseVec(lvalue.dimension(), lvalue.sparseSize()) { // copy semantics
        std::copy(lvalue.indices, lvalue.indices + lvalue.sparseSize(), indices);
        std::copy(lvalue.values, lvalue.values + lvalue.sparseSize(), values);
    }

    ~SparseVec();

//    double operator [](int i) const;
//    double &operator [](int i);

    int sparseSize() const { return indices[0]; }
    int &sparseSize() { return indices[0]; }

    void toDense(double *denseVec) const;
//    std::pair<int,double> entry(int k) const;
    void add(int i, double v);
    void clear();
    int capacity() const { return indices==NULL?0:header().capacity; }
    void ensureCapacity(int size);
    int dimension() const { return indices==NULL?0:header().dimension; }
    void setDimension(int dim) {
        if(indices == NULL) ensureCapacity(0);
        header().dimension = dim;
    }
//    double dotProd(double *dense) const;

    int maxNonZeroIndex();

    SparseVec &operator =(const SparseVec &lvalue) {
        ensureCapacity(lvalue.sparseSize());
        std::copy(lvalue.indices, lvalue.indices + lvalue.sparseSize(), indices);
        std::copy(lvalue.values, lvalue.values + lvalue.sparseSize(), values);
        return *this;
    }

    SparseVec &operator =(SparseVec &&rvalue) {
        moveFrom(rvalue);
        return *this;
    }

protected:
    void moveFrom(SparseVec &rvalue) {
        indices = rvalue.indices;
        values = rvalue.values;
        rvalue.indices = NULL;
    }

    Header &header() {
        return *(Header *)values;
    }

    const Header &header() const {
        return *(const Header *)values;
    }



};

std::ostream &operator <<(std::ostream &out, const SparseVec &sVector);

#endif //GLPKTEST_SPARSEVEC_H
