//
// Created by daniel on 27/04/2021.
//
#include <iomanip>
#include <iostream>

#include "../include/glpkpp.h"

namespace glp {
    Constraint LinearSum::operator ==(double c) const {
        return Constraint(c, *this, c);
    }

    Constraint LinearSum::operator<=(double c) const {
        return Constraint(-std::numeric_limits<double>::infinity(), *this, c);
    }

    Constraint LinearSum::operator>=(double c) const {;
        return Constraint(c, *this, std::numeric_limits<double>::infinity());
    }


    std::ostream &operator<<(std::ostream &out, const LinearSum &sVector) {
        int i;
        int iMax = sVector.maxNonZeroIndex();
        double dense[iMax];
        sVector.toDense(dense, iMax);
        for (i = 0; i < iMax; ++i) {
            out << std::setw(8) << dense[i] << "X" << i+1;
            if(i == iMax-1) out << "\t"; else out << " +\t";
        }
        return out;
    }


};
