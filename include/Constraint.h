
//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKPP_CONSTRAINT_H
#define GLPKPP_CONSTRAINT_H

class Constraint {
public:

    SparseVec coefficients;
    double upperBound;
    double lowerBound;

    Constraint(): Constraint(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()) {}
    Constraint(double lowerBound, double upperBound);
    Constraint(double lowerBound, SparseVec coefficients, double upperBound);
//    Constraint(const LinearSum &sum): Constraint(-std::numeric_limits<double>::infinity(), sum, std::numeric_limits<double>::infinity()) {}

    Constraint & operator <=(double upperBound);

    bool isValidSolution(const std::vector<double> &X) const;

    int highestVar() const;
//    Constraint & operator +=(std::pair<double,X> &);

//    double operator [](int varId) const {
//        auto pValue = coefficients.find(varId);
//        return (pValue == coefficients.end())?0.0:pValue->second;
//    }
//    double & operator [](int varId) { return coefficients[varId]; }


};


std::ostream &operator <<(std::ostream &out, const Constraint &constraint);

#endif //GLPKPP_CONSTRAINT_H
