//
// Created by daniel on 26/04/2021.
//

#include "../include/glpkpp.h"
#include "../include/Constraint.h"


namespace glp {
    int Constraint::glpBoundType() const {
        return lowerBound == -std::numeric_limits<double>::infinity()?
               (upperBound == std::numeric_limits<double>::infinity()?GLP_FR:GLP_UP):
               (upperBound == std::numeric_limits<double>::infinity()?
                    GLP_LO:(upperBound == lowerBound?GLP_FX:GLP_DB)
                    );
    }


    Constraint &Constraint::operator<=(double upperBound) {
        this->upperBound = upperBound;
        return *this;
    }


    Constraint::Constraint(double lowerBound, double upperBound):
    upperBound(upperBound),
    lowerBound(lowerBound) { }


    Constraint::Constraint(double lowerBound, const LinearSum &sum, double upperBound):
    coefficients(sum),
    upperBound(upperBound),
    lowerBound(lowerBound) { }

//    Constraint &Constraint::operator+=(std::pair<double, X> &entry) {
//        coefficients[entry.second.id] += entry.first;
//        return *this;
//    }


    std::ostream &operator <<(std::ostream &out, const Constraint &constraint) {
        out << constraint.lowerBound << " <= " << constraint.coefficients << " <= " << constraint.upperBound;
        return  out;
    }

}