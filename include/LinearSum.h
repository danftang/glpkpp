//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKPP_LAZYLINEAREXPRESSION_H
#define GLPKPP_LAZYLINEAREXPRESSION_H


// Represents a sum of terms sum_i c_ix_i
class LinearSum: public std::vector<std::pair<double,int>> {
public:

    LinearSum(double coefficient, const X &variable) {
        push_back(std::pair<double,int>(coefficient, variable.id));
    }


    LinearSum &operator +=(const LinearSum &term) {
        for(auto entry: term) {
            push_back(entry);
        }
        return *this;
    }


    friend LinearSum operator +(LinearSum lhs, const LinearSum &rhs) {
        lhs += rhs;
        return lhs; // should use move constructor (or elision?)
    }


    Constraint operator ==(double c) const {
        return Constraint(c, *this, c);
    }

    Constraint operator <=(double c) const {
        return Constraint(-std::numeric_limits<double>::infinity(), *this, c);
    }

    Constraint operator >=(double c) const {;
        return Constraint(c, *this, std::numeric_limits<double>::infinity());
    }
};

inline Constraint operator <=(double c, const LinearSum &linExp) {
    return linExp >= c;
}

inline LinearSum operator *(double coefficient, const X &variable) {
    return LinearSum(coefficient, variable);
}






#endif //GLPKPP_LAZYLINEAREXPRESSION_H
