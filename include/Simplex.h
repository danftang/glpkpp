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
    static constexpr int excludeFixed = 1;  // exclude non-basic fixed variables
    static constexpr int shift = 0;         // don't shift bounds to zero. Don't change this!

    enum BoundType {
        UNBOUNDED = 0,
        LOWERBOUNDED = 1,
        UPPERBOUNDED = 2,
        DOUBLEBOUNDED = 3
    };

    std::vector<int>    kProbTokSim;    // map from variable ids of the original problem to variable ids in this simplex
    std::vector<int>            kSimTokProb;    // map from variable ids in this simplex to variable ids of the original problem
    std::vector<double>         pi;     // pi = c_B*B' where c_B = objective of basic vars so reduced objective = c_B*B'*N + c_N = pi*N + c_N
    std::vector<double>         lpSolution;     // solution in the original variables (excluding auxiliary vars)

    Simplex(const Problem &prob);
    ~Simplex();

    int nRows() { return this->m; }
    int nNCols() { return this->n; }

    std::vector<double> tableauRow(int i);
    std::vector<double> tableauCol(int j);
    double reducedObjective(int j);             // value of the j'th (1 <= j <= n-m) element of the reduced objective
    void pivot(int i, int j, const std::vector<double> &pivotCol);
    void pivot(int i, int j)                    { pivot(i,j,tableauCol(j)); }
    bool isAtUpperBound(int j)                  { return flag[j]; }
    void isAtUpperBound(int j, bool setUpper)   { flag[j] = setUpper; }
    BoundType boundType(int k);
    void syncWith(Problem &prob);
    const std::vector<double> &X(); // current solution in original problem coordinates (excluding auxiliaries)

protected:

    bool lpSolutionIsValid()                { return lpSolution[0] == 0; }
    void lpSolutionIsValid(bool setValid)   { lpSolution[0] = !setValid; }
    void calculateLpSolution();
    void piIsValid(bool setValid);
    bool piIsValid();
};

std::ostream &operator <<(std::ostream &out, Simplex &simplex);


#endif //GLPKTEST_GLPSIMPLEX_H
