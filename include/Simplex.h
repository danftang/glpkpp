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
    static constexpr double zeroTol = 1e-7; // absolute value under which a variable is considered to be zero

    enum BoundType {
        UNBOUNDED = 0,
        LOWERBOUNDED = 1,
        UPPERBOUNDED = 2,
        DOUBLEBOUNDED = 3
    };


    std::vector<int>    kProbTokSim;    // map from variable ids of the original problem to variable ids in this simplex (0 means removed from Sim)
    std::vector<int>    kSimTokProb;    // map from variable ids in this simplex to variable ids of the original problem
    std::vector<double> pi;             // pi = c_B*B' where c_B = objective of basic vars so reduced objective = c_B*B'*N + c_N = pi*N + c_N
    std::vector<double> lpSolution;     // solution in the original variables (excluding auxiliary vars), (non-zero zero'th element signals needs updating)
    Problem &           originalProblem;
    double *            beta;           // current values of the basic vars
//    std::vector<double> rCost;    // reduced cost

    Simplex(Problem &prob);
    ~Simplex();

    int nBasic() const { return this->m; }
    int nNonBasic() const { return this->n - this->m; }
    int nVars() const { return this->n; }

    std::vector<double> tableauRow(int i);
    std::vector<double> tableauCol(int j);

    void recalculatePi();
    void setObjective(const SparseVec &objective);
    std::vector<double> reducedCost();
    double reducedCost(int j);             // value of the j'th (1 <= j <= n-m) element of the reduced objective

    std::vector<double> piTimesMinusN(const std::vector<double> &pi);
    void pivot(int i, int j, const std::vector<double> &pivotCol, bool leavingVarToUpperBound);
    void pivot(int i, int j, const std::vector<double> &pivotCol);
    void pivot(int i, int j)                    { if(i>0) pivot(i,j,tableauCol(j)); else pivot(i,j,std::vector<double>()); }
    bool isAtUpperBound(int j)                  { return flag[j]; }
    void isAtUpperBound(int j, bool setUpper)   { flag[j] = setUpper; }
    BoundType boundType(int k);
    void syncWithLP();

    const std::vector<double> &X(); // current solution in original problem coordinates (excluding auxiliaries)
    double nonBasicValue(int j) { return isAtUpperBound(j)?u[head[nBasic()+j]]:l[head[nBasic()+j]]; }

    void btran(std::vector<double> &rowVec) { if(!valid) spx_factorize(this); bfd_btran(bfd, rowVec.data()); } // in-place btran
    void ftran(std::vector<double> &colVec) { if(!valid) spx_factorize(this); bfd_ftran(bfd, colVec.data()); } // in-place ftran

protected:

    bool lpSolutionIsValid()                { return lpSolution[0] == 0; }
    void lpSolutionIsValid(bool setValid)   { lpSolution[0] = !setValid; }
    void calculateLpSolution();
    void piIsValid(bool setValid);
    bool piIsValid();
};

std::ostream &operator <<(std::ostream &out, Simplex &simplex);


#endif //GLPKTEST_GLPSIMPLEX_H
