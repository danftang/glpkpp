//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKPP_LAZYLINEAREXPRESSION_H
#define GLPKPP_LAZYLINEAREXPRESSION_H


// Represents a sum of terms sum_i c_ix_i
class Constraint;

class LinearSum {
public:
    std::map<int,double> coefficients;

    class Monomial {
    public:
        double multiplier;
        X variableId;

        Monomial(double multiplier, const X &variableId) : multiplier(multiplier), variableId(variableId) {}

    };

    LinearSum() {}

    LinearSum(const Monomial &monomial) {
        coefficients[monomial.variableId] = monomial.multiplier;
    }


    LinearSum(LinearSum &&rvalue): coefficients(std::move(rvalue.coefficients)) { }

    LinearSum(const LinearSum &lvalue) : coefficients(lvalue.coefficients) { }

    double &operator [](int index) {
        std::map<int,double>::iterator it = coefficients.find(index);
        if(it == coefficients.end()) return coefficients[index] = 0.0;
        return it->second;
    }

    LinearSum &operator +=(const LinearSum &other) {
        for(auto [varId, coeff]: other.coefficients) {
            (*this)[varId] += coeff;
        }
        return *this;
    }

    LinearSum &operator -=(const LinearSum &other) {
        for(auto [varId, coeff]: other.coefficients) {
            (*this)[varId] -= coeff;
        }
        return *this;
    }


    LinearSum &operator+=(const Monomial &monomial) {
        (*this)[monomial.variableId] += monomial.multiplier;
        return *this;
    }


    LinearSum &operator-=(const Monomial &monomial) {
        (*this)[monomial.variableId] -= monomial.multiplier;
        return *this;
    }


    LinearSum operator +(const LinearSum &rhs) const & {
        LinearSum lhs(*this);
        lhs += rhs;
        return lhs;
    }

    LinearSum operator +(const LinearSum &rhs) && {
        LinearSum lhs(std::move(*this));
        lhs += rhs;
        return lhs;
    }

    LinearSum operator+(const Monomial &rhs) const &{
        LinearSum lhs(*this);
        lhs += rhs;
        return lhs;
    }

    LinearSum operator+(const Monomial &rhs) &&{
        LinearSum lhs(std::move(*this));
        lhs += rhs;
        return lhs;
    }

    LinearSum operator -(const LinearSum &rhs) const & {
        LinearSum lhs(*this);
        lhs -= rhs;
        return lhs;
    }

    LinearSum operator -(const LinearSum &rhs) && {
        LinearSum lhs(std::move(*this));
        lhs -= rhs;
        return lhs;
    }

    LinearSum operator-(const Monomial &rhs) const &{
        LinearSum lhs(*this);
        lhs -= rhs;
        return lhs;
    }

    LinearSum operator-(const Monomial &rhs) &&{
        LinearSum lhs(std::move(*this));
        lhs -= rhs;
        return lhs;
    }

    SparseVec toSparseVec() const {
        SparseVec v(coefficients.size());
        int i=0;
        for(auto [key,val] : coefficients) {
            v.indices[i] = key;
            v.values[i++] = val;
        }
        return v;
    }


    Constraint operator==(double c) const;

    Constraint operator<=(double c) const;

    Constraint operator>=(double c) const;

    friend Constraint operator<=(double c, const LinearSum &linExp);
};


inline LinearSum::Monomial operator *(double coefficient, const X &variable) {
    return LinearSum::Monomial(coefficient, variable);
}

inline LinearSum operator +(const LinearSum::Monomial &lhs, const LinearSum::Monomial &rhs) {
    LinearSum sum(lhs);
    sum += rhs;
    return sum;
}

inline LinearSum operator -(const LinearSum::Monomial &lhs, const LinearSum::Monomial &rhs) {
    LinearSum result(lhs);
    result -= rhs;
    return result;
}

std::ostream &operator<<(std::ostream &out, const LinearSum &sVector);


#endif //GLPKPP_LAZYLINEAREXPRESSION_H
