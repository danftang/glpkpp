//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_SPARSEVEC_H
#define GLPKTEST_SPARSEVEC_H

class SparseVec {
public:
    std::vector<int> indices;   // one-based array of indices of non-zero elements (zero'th element stores number of non-zeros)
    std::vector<double> values; // one-based array of values of non-zero elements (zero'th element stores dimension of this vector)

    SparseVec() { }

    SparseVec(int sparseSize): indices(sparseSize), values(sparseSize) { }

    SparseVec(SparseVec &&rvalue) { // move semantics
        swap(rvalue);
    }

    SparseVec(const SparseVec &lvalue): indices(lvalue.indices), values(lvalue.values) { // copy semantics
    }

//    double operator [](int i) const;
//    double &operator [](int i);

    int sparseSize() const { return indices.size(); }

    void toDense(double *denseVec, int size) const;
    void add(int i, double v);
    void clear();
//    int capacity() const { return std::min(indices.capacity(),values.capacity()); }
//    int dimension() const { return (int)values[0]; }
    void resize(int size) { indices.resize(size); values.resize(size); }

    int *glpkIndexArray() { return indices.data() - 1; }
    double *glpkValueArray() { return values.data() - 1; }

    const int *glpkIndexArray() const { return indices.data() - 1; }
    const double *glpkValueArray() const { return values.data() - 1; }

    //    void setDimension(int dim) { values[0] = (double)dim; }
//    double dotProd(double *dense) const;

    int maxNonZeroIndex() const;

    double operator[](int index) {
        for(int i = 0;i<sparseSize(); ++i) {
            if(indices[i] == index) return values[i];
        }
        return 0.0;
    }

    SparseVec &operator +=(const SparseVec &other) {
        for(int i=0; i < other.sparseSize(); ++i) {
            add(other.indices[i], other.values[i]);
        }
        return *this;
    }

    friend SparseVec operator +(SparseVec lhs, const SparseVec &rhs) {
        lhs += rhs;
        return lhs; // should use move constructor (or elision?)
    }



    // copy semantics
    SparseVec &operator =(const SparseVec &lvalue) {
        indices = lvalue.indices;
        values = lvalue.values;
        return *this;
    }

    // move semantics
    SparseVec &operator =(SparseVec &&rvalue) {
        swap(rvalue);
        return *this;
    }

protected:
    void swap(SparseVec &rvalue) {
//        std::cout << "Moving " << rvalue.indices << std::endl;
        indices.swap(rvalue.indices);
        values.swap(rvalue.values);
    }

};

std::ostream &operator <<(std::ostream &out, const SparseVec &sVector);

#endif //GLPKTEST_SPARSEVEC_H
