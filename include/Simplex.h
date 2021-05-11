//
// Created by daniel on 15/04/2021.
//

#ifndef GLPKTEST_GLPSIMPLEX_H
#define GLPKTEST_GLPSIMPLEX_H

// n - total number of variables
// m - total number of constraints
// i - row or basic variable identifier (1..m)
// j - non-basic variable identifier (1..n-m)
// k - variable (column) identifier (1..n)
class Simplex: public SPXLP {
public:
    int *map;       // map of variables in this simplex to those of the original problem
    double *pi;     // reduced objective = c_B*B'*N + c_N = pi*N + c_N

    Simplex(const Problem &prob);
    ~Simplex();

    int nRows() { return this->m; }
    int nVars() { return this->n; }

    std::vector<double> tableauRow(int i);
    std::vector<double> tableauCol(int j);
    double reducedObjective(int j);             // value of the j'th (1 <= j <= n-m) element of the reduced objective
    void pivot(int i, int j, bool leavingVarToUpperBound, const std::vector<double> &pivotCol);
    void pivot(int i, int j, bool leavingVarToUpperBound) { pivot(i,j,leavingVarToUpperBound,tableauCol(j)); }
    bool isAtLowerBound(int j) { return flag[j] == 0; }

};

std::ostream &operator <<(std::ostream &out, Simplex &simplex);


#endif //GLPKTEST_GLPSIMPLEX_H
