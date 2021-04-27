//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKPP_CONSTRAINT_H
#define GLPKPP_CONSTRAINT_H

class LinearSum;

class Constraint {
public:

    std::map<int,double> coefficients;
    double upperBound;
    double lowerBound;

    Constraint(double lowerBound, const LinearSum &sumCoefficients, double upperBound);
    Constraint(const LinearSum &sum): Constraint(-std::numeric_limits<double>::infinity(), sum, std::numeric_limits<double>::infinity()) {}

    Constraint & operator <=(double upperBound);

    int glpBoundType() const;

};

std::ostream &operator <<(std::ostream &out, const Constraint &constraint);

#endif //GLPKPP_CONSTRAINT_H
