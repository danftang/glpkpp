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
        ind = new int(size);
        vec = new double(size);
        for(int i=0; i<n; ++i) vec[i] = 0.0;
    };

    ~FVSVector() {
        delete ind;
        delete vec;
    }


    void recalculateNonZeroIndices() {
        nnz = 0;
        for(int i=0; i<n; ++i) {
            if(fabs(vec[i]) > tol) ind[nnz++] = i;
        }
    }

};


#endif //GLPKPP_FVSVECTOR_H
