//
// Created by daniel on 26/04/2021.
//

#include "../include/glpkpp.h"
#include "../include/Constraint.h"


namespace glp {

    Constraint &Constraint::operator<=(double upperBound) {
        if(upperBound < this->upperBound) this->upperBound = upperBound;
        return *this;
    }


    Constraint::Constraint(double lowerBound, double upperBound):
    upperBound(upperBound),
    lowerBound(lowerBound) { }


    Constraint::Constraint(double lowerBound, SparseVec sum, double upperBound):
    coefficients(std::move(sum)),
    upperBound(upperBound),
    lowerBound(lowerBound) { }

    bool Constraint::isValidSolution(const std::vector<double> &X) const {
        double mx = coefficients * X;
        return (lowerBound - Simplex::zeroTol <= mx) && (mx <= upperBound + Simplex::zeroTol);
    }

    int Constraint::highestVar() const {
        return coefficients.maxNonZeroIndex();
    }

//    Constraint &Constraint::operator+=(std::pair<double, X> &entry) {
//        coefficients[entry.second.id] += entry.first;
//        return *this;
//    }


    std::ostream &operator <<(std::ostream &out, const Constraint &constraint) {
        out << constraint.lowerBound << " <= ";
        bool first = true;
        for(int term=0; term < constraint.coefficients.sparseSize(); ++term) {
            if(first) first = false; else out << " + \t";
            out << constraint.coefficients.values[term] << "X(" << constraint.coefficients.indices[term] << ")";
        }
        out << " <= " << constraint.upperBound;
        return  out;
    }

}