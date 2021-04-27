//
// Created by daniel on 09/04/2021.
//

#include <iomanip>
#include <iostream>
#include "../include/glpkpp.h"

namespace glpkpp {
    SparseVec::SparseVec(int dim, int capacity) {
        n = dim;
        int cap = capacity==-1?dim:capacity;
        nnz = 0;
        indices = (int *)new int[cap] - 1; // one based indexing :-(
        values = (double *)new double[cap] - 1;
    }


    SparseVec::SparseVec(int dimension, const std::map<int, double> &entries): SparseVec(dimension, entries.size()) {
        for(auto entry: entries) add(entry.first, entry.second);
    }


    SparseVec::~SparseVec() {
        delete [](indices + 1);
        delete [](values + 1);
    }

    int SparseVec::maxNonZeroIndex() {
        int max = 0;
        for(int i=1; i<=nnz; ++i) {
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

    void SparseVec::add(int i, double v) {
        if (v != 0.0) {
            indices[++nnz] = i;
            values[nnz] = v;
        } // doesn't delete if index already exists
    }

    void SparseVec::clear() {
        nnz = 0;
    }

    // to zero-based dense array
    void SparseVec::toDense(double *dense) const {
        int i;
        for (i = 0; i < n; ++i) { dense[i] = 0.0; }
        for (i = 1; i <= nnz; ++i) {
            if (indices[i] > n || indices[i] < 1)
                std::cout << "Out of range index[" << i << "] = " << indices[i] << " -> " << values[i] << std::endl;
            dense[indices[i] - 1] = values[i];
        }
    }

    void SparseVec::entry(int k, std::pair<int, double> &retEntry) {
        retEntry.first = indices[k];
        retEntry.second = values[k];
    }

    // sets capacity (also invalidates any data)
    void SparseVec::setCapacity(int size) {
        nnz = 0;
        delete [](indices + 1);
        delete [](values + 1);
        indices = ((int *)new int[size]) - 1;
        values = ((double *)new double[size]) - 1;
    }


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
        int i;
        double dense[sVector.n];
        sVector.toDense(dense);
        for (i = 0; i < sVector.n; ++i) out << std::setw(12) << dense[i] << "\t";
        return out;
    }
};
