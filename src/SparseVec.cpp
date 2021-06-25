//
// Created by daniel on 09/04/2021.
//

#include <iomanip>
#include <iostream>
#include "../include/glpkpp.h"

namespace glp {
//    SparseVec::SparseVec(int dim, int capacity):
//    indices(capacity==-1?dim+1:capacity+1),
//    values(capacity==-1?dim+1:capacity+1) {
//        setDimension(dim);
//    }


//    SparseVec::SparseVec(int dimension, const std::map<int, double> &entries): SparseVec(dimension, entries.size()) {
//        for(auto entry: entries) add(entry.first, entry.second);
//    }



    int SparseVec::maxNonZeroIndex() const {
        int max = 0;
        for(int i=0; i<sparseSize(); ++i) {
            if(indices[i] > max) max = indices[i];
        }
        return max;
    }

//    double &SparseVec::operator[](int i) {
//        for (int k = 1; k <= nnz; k++) {
//            if (ind[k] == i + 1) return vec[k];
//        }
//        ind[++nnz] = i + 1;
//        return vec[nnz]; // allows zero values
//    }
//
//
//    double SparseVec::operator[](int i) const {
//        for (int k = 1; k < nnz; k++) {
//            if (ind[k] == i + 1) return vec[k];
//        }
//        return 0.0;
//    }

    // Adds element,
    // doesn't delete if index already exists
    void SparseVec::add(int i, double v) {
        if (v != 0.0) {
            indices.push_back(i);
            values.push_back(v);
        }
    }

    void SparseVec::clear() {
        indices.clear();
        values.clear();
    }

    // to zero-based dense array
    std::vector<double> SparseVec::toDense() const {
        int vecSize = maxNonZeroIndex() + 1;
        std::vector<double> denseVec(vecSize,0.0);
        for (int i = 0; i < sparseSize(); ++i) {
            denseVec[indices[i]] = values[i];
        }
        return denseVec;
    }



//    void SparseVec::entry(int k, std::pair<int, double> &retEntry) {
//        retEntry.first = indices[k];
//        retEntry.second = values[k];
//    }


//    double SparseVec::dotProd(double *dense) const {
//        double *dense1base = dense - 1;
//        double dp = 0.0;
//        int j;
//        for (int i = 1; i <= nnz; ++i) {
//            j = ind[i];
//            dp += dense1base[j] * vec[i];
//        }
//        return dp;
//    }


    std::ostream &operator<<(std::ostream &out, const SparseVec &sVector) {
//        int i;
//        int iMax = sVector.maxNonZeroIndex();
//        double dense[iMax];
//        sVector.toDense(dense, iMax);
//        for (i = 0; i < iMax; ++i) out << std::setw(12) << dense[i] << "\t";
        for(auto entry: sVector) {
            out << "[" << entry.index << "] -> " << entry.value << "  ";
        }
        return out;
    }
};
