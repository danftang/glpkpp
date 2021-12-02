//
// Created by daniel on 09/07/2021.
//

#ifndef GLPKPP_FVSVECTOR_H
#define GLPKPP_FVSVECTOR_H

#include <cmath>
#include "glpkpp.h"

class FVSVector: public FVS {
    static constexpr double tol = 1e-8;

    explicit FVSVector(int size) {
        n = size;
        nnz = 0;
        ind = new int(size+1);
        vec = new double(size+1);
        for(int i=0; i<=n; ++i) vec[i] = 0.0;
    };


    explicit FVSVector(const std::vector<double> &dense) {
        ind = new int(dense.size());
        vec = new double(dense.size());
        nnz = 0;
        n = dense.size()-1;
        for(int i=1; i <= n; ++i) {
            double v = dense[i];
            vec[i] = v;
            if(v != 0.0) ind[++nnz] = i;
        }
    };


    ~FVSVector() {
        delete ind;
        delete vec;
    }


    void clear() {
        for(int z=1; z<=nnz; ++z) {
            vec[ind[z]] = 0.0;
        }
        nnz = 0;
    }


    void recalculateNonZeroIndices() {
        nnz = 0;
        for(int i=0; i<n; ++i) {
            if(fabs(vec[i]) > tol) ind[nnz++] = i;
        }
    }

};


#endif //GLPKPP_FVSVECTOR_H
