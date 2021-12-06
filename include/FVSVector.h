//
// Created by daniel on 09/07/2021.
//

#ifndef GLPKPP_FVSVECTOR_H
#define GLPKPP_FVSVECTOR_H


class FVSVector {
public:
    static constexpr double tol = 1e-8;
    std::vector<int>    indices;    // indices of non-zero values
    std::vector<double> vec;        // dense vector of values

    FVSVector() = default;

    FVSVector(FVSVector &&other): indices(std::move(other.indices)), vec(std::move(other.vec)) { }

    explicit FVSVector(int maxIndex): vec(maxIndex+1,0.0) { };


    explicit FVSVector(std::vector<double> &&dense): vec(std::move(dense)) {
        recalculateIndices();
    };


    void clear() {
        for(int i: indices) vec[i] = 0.0;
        indices.clear();
    }

    int sparseSize() const { return indices.size(); }

    double operator [](int i) const {
        return vec[i];
    }

    FVSVector &operator =(std::vector<double> &&dense) {
        vec = std::move(dense);
        recalculateIndices();
        return *this;
    }

    operator FVS() {
        FVS fvs;
        fvs.n = vec.size()-1;
        fvs.nnz = indices.size();
        fvs.ind = indices.data()-1;
        fvs.vec = vec.data();
        return fvs;
    }

    void recalculateIndices() {
        indices.resize(0);
        for(int i=0; i<vec.size(); ++i) {
            if(std::fabs(vec[i]) > tol) indices.push_back(i);
        }
    }

};


#endif //GLPKPP_FVSVECTOR_H
