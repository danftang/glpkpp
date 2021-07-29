//
// Created by daniel on 27/04/2021.
//
#include <iomanip>
#include <iostream>

#include "../include/glpkpp.h"

namespace glp {
    Constraint LinearSum::operator ==(double c) const {
        return Constraint(c, toSparseVec(), c);
    }

    Constraint LinearSum::operator<=(double c) const {
        return Constraint(-std::numeric_limits<double>::infinity(), toSparseVec(), c);
    }

    Constraint LinearSum::operator>=(double c) const {;
        return Constraint(c, toSparseVec(), std::numeric_limits<double>::infinity());
    }

    Constraint operator<=(double c, const LinearSum &linExp) {
        return linExp >= c;
    }


    std::ostream &operator <<(std::ostream &out, const LinearSum &sum) {
        bool first = true;
        for (auto [varId, coeff] : sum.coefficients) {
            if(first) first = false; else out << " +\t";
            out << std::setw(8) << coeff << "X" << varId;
        }
        return out;
    }


};
