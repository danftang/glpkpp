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
        std::vector<double> dense = sVector.toDense();
        for (int i = 0; i < dense.size(); ++i) {
            out << std::setw(8) << dense[i] << "X" << i;
            if(i == dense.size()-1) out << "\t"; else out << " +\t";
        }
        return out;
    }


};
