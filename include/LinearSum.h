//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKPP_LAZYLINEAREXPRESSION_H
#define GLPKPP_LAZYLINEAREXPRESSION_H


// Represents a sum of terms sum_i c_ix_i
class Constraint;

class LinearSum: public SparseVec {
public:
    LinearSum() { }

//    explicit LinearSum(int sparseSize): SparseVec(sparseSize) { }

    LinearSum(LinearSum &&rvalue) {
        swap(rvalue);
    }

    LinearSum(const LinearSum &lvalue): SparseVec(lvalue) { }

    LinearSum &operator +=(const LinearSum &term) {
        SparseVec::operator +=(term);
        return *this;
    }

    // copy semantics
    LinearSum &operator =(const LinearSum &lvalue) {
        SparseVec::operator=(lvalue);
        return *this;
    }

    // move semantics
    LinearSum &operator =(LinearSum &&rvalue) {
        swap(rvalue);
        return *this;
    }

    friend LinearSum operator +(LinearSum lhs, const LinearSum &rhs) {
        lhs += rhs;
        return lhs; // should use move constructor (or elision?)
    }

    Constraint operator ==(double c) const ;

    Constraint operator <=(double c) const;

    Constraint operator >=(double c) const;
};

inline LinearSum operator *(double coefficient, const X &variable) {
    LinearSum term;
    term.add(variable, coefficient);
    return term;
}

std::ostream &operator<<(std::ostream &out, const LinearSum &sVector);


#endif //GLPKPP_LAZYLINEAREXPRESSION_H
