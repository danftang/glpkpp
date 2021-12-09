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
    static constexpr double zeroTol = 1e-7; // absolute val under which a variable is considered to be zero

    enum BoundType {
        UNBOUNDED = 0,
        LOWERBOUNDED = 1,
        UPPERBOUNDED = 2,
        DOUBLEBOUNDED = 3
    };


    std::vector<int>    kProbTokSim;    // map from variable ids of the original problem to variable ids in this simplex (0 means removed from Sim)
    std::vector<int>    kSimTokProb;    // map from variable ids in this simplex to variable ids of the original problem
    std::vector<double> pi;             // pi = c_B*B' where c_B = objective of basic vars so reduced objective = c_B*B'*N + c_N = pi*N + c_N
    std::vector<double> lpSolution;     // solution in the original problem variables (excluding auxiliary vars), (non-zero zero'th element signals needs updating)
//    Problem &           originalProblem;
    double *            beta;           // current values of the basic vars
//    std::vector<double> rCost;    // reduced cost

    Simplex(const Problem &prob);
    Simplex(const Simplex &other) = delete; // copying is usually accidental when passing, so better to disable it.
    Simplex(Simplex &&moveFrom);
    ~Simplex();

    int nBasic() const { return m; }
    int nNonBasic() const { return n - m; }
    int nVars() const { return this->n; }
    bool isAuxiliary(int k) const { return kSimTokProb[k] <= nBasic(); } // was the k'th var auxiliary in the original problem?
    bool isStructural(int k) const { return kSimTokProb[k] > nBasic(); } // was the k'th var structural in the original problem?

    std::vector<double> tableauRow(int i);
    std::vector<double> tableauCol(int j);

    void recalculatePi();
    std::vector<double> reducedCost();
    double reducedCost(int j);             // val of the j'th (1 <= j <= n-m) element of the reduced objective

    std::vector<double> piTimesMinusN(const std::vector<double> &pi);

    template<typename COLUMN>
    void pivot(int i, int j, const COLUMN &pivotCol, bool leavingVarToUpperBound);
    void pivot(int i, int j, const std::vector<double> &pivotCol);
    void pivot(int i, int j)                    { if(i>0) pivot(i,j,tableauCol(j)); else pivot(i,j,std::vector<double>()); }
    bool isAtUpperBound(int j) const            { return flag[j]; }
    void isAtUpperBound(int j, bool setUpper)   { flag[j] = setUpper; }
    BoundType boundType(int k);
    void syncWithLP(Problem &);

    const std::vector<double> &X(); // current exactEndState in original problem coordinates (excluding auxiliaries)
    double nonBasicValue(int j) { return isAtUpperBound(j)?u[head[nBasic()+j]]:l[head[nBasic()+j]]; }

    void btran(std::vector<double> &rowVec) { if(!valid) spx_factorize(this); bfd_btran(bfd, rowVec.data()); } // in-place btran
    void ftran(std::vector<double> &colVec) { if(!valid) spx_factorize(this); bfd_ftran(bfd, colVec.data()); } // in-place ftran

    void evalBeta() { spx_eval_beta(this, beta); }

protected:

    bool lpSolutionIsValid() const          { return lpSolution[0] == 0; }
    void lpSolutionIsValid(bool setValid)   { lpSolution[0] = !setValid; }
    void calculateLpSolution();
    void piIsValid(bool setValid);
    bool piIsValid() const;

    template<typename COLUMN>
    void updateBetaAndLPSolution(int i, bool leavingVarToUpperBound, int j, const COLUMN &pivotCol);

    void updateBetaAndLPSolution(int p, double delta_q, const SparseVec &pivotCol);
    void updateBetaAndLPSolution(int p, double delta_q, const std::vector<double> &pivotCol);
    void updateBetaAndLPSolution(int p, double delta_q, const FVSVector &pivotCol);
};

std::ostream &operator <<(std::ostream &out, Simplex &simplex);


// p is leaving row, q is entering col
template<typename COLUMN>
void Simplex::updateBetaAndLPSolution(int p, bool p_flag, int q, const COLUMN &pivotCol) {
//    assert(pivotCol.size() == m+1);
    int i, k;
    double delta_p, delta_q;
    if (p < 0)
    {  /* special case: xN[q] goes to its opposite bound */
        assert(1 <= q && q <= n-m);
        /* xN[q] should be double-bounded variable */
        k = head[m+q]; /* x[k] = xN[q] */
        assert(l[k] != -DBL_MAX && u[k] != +DBL_MAX && l[k] != u[k]);
        /* determine delta xN[q] */
        if (flag[q])
        {  /* xN[q] goes from its upper bound to its lower bound */
            delta_q = l[k] - u[k];
            int kProb = kSimTokProb[k];
            if(kProb>m) lpSolution[kProb-m] = l[k];
        }
        else
        {  /* xN[q] goes from its lower bound to its upper bound */
            delta_q = u[k] - l[k];
            int kProb = kSimTokProb[k];
            if(kProb>m) lpSolution[kProb-m] = u[k];
        }
    }
    else
    {  /* xB[p] leaves the basis, xN[q] enters the basis */
        assert(1 <= p && p <= m);
        assert(1 <= q && q <= n-m);
        /* determine delta xB[p] */
        k = head[p]; /* x[k] = xB[p] */
        if (p_flag)
        {  /* xB[p] goes to its upper bound */
            assert(l[k] != u[k] && u[k] != +DBL_MAX);
            delta_p = u[k] - beta[p];
            int kProb = kSimTokProb[k];
            if(kProb>m) lpSolution[kProb-m] = u[k];
        }
        else if (l[k] == -DBL_MAX)
        {  /* unbounded xB[p] becomes non-basic (unusual case) */
            assert(u[k] == +DBL_MAX);
            delta_p = 0.0 - beta[p];
            int kProb = kSimTokProb[k];
            if(kProb>m) lpSolution[kProb-m] = 0.0;
        }
        else
        {  /* xB[p] goes to its lower bound or becomes fixed */
            delta_p = l[k] - beta[p];
            int kProb = kSimTokProb[k];
            if(kProb>m) lpSolution[kProb-m] = l[k];
        }
        /* determine delta xN[q] */
        delta_q = delta_p / pivotCol[p];
        /* compute new beta[p], which is the val of xN[q] in the
         * adjacent basis */
        k = head[m+q]; /* x[k] = xN[q] */
        if (flag[q])
        {  /* xN[q] has its upper bound active */
            assert(l[k] != u[k] && u[k] != +DBL_MAX);
            beta[p] = u[k] + delta_q;
        }
        else if (l[k] == -DBL_MAX)
        {  /* xN[q] is non-basic unbounded variable */
            assert(u[k] == +DBL_MAX);
            beta[p] = 0.0 + delta_q;
        }
        else
        {  /* xN[q] has its lower bound active or is fixed (latter
             * case is unusual) */
            beta[p] = l[k] + delta_q;
        }
        int kProb = kSimTokProb[k];
        if(kProb>m) lpSolution[kProb-m] = beta[p];
    }
    /* compute new beta[i] for all i != p */
    updateBetaAndLPSolution(p, delta_q, pivotCol);
}


template<typename COLUMN>
void Simplex::pivot(int i, int j, const COLUMN &pivotCol, bool leavingVarToUpperBound) {
    if(i < 1) { // column j goes to opposite bound
        if(isAtUpperBound(j) != leavingVarToUpperBound) {
            updateBetaAndLPSolution(-1, 0, j, pivotCol);
            isAtUpperBound(j, leavingVarToUpperBound);
        }
    } else { // basis update
        int enteringk = head[m+j];
        assert(fabs(pivotCol[i]) > zeroTol);
        assert(u[head[i]] != l[head[i]]);       // TODO: deal with fixed variables leaving the basis
//        assert(u[enteringk] != l[enteringk]);   // fixed variables should never enter the basis
        updateBetaAndLPSolution(i, leavingVarToUpperBound, j, pivotCol);
        int err = spx_update_invb(this, i, enteringk);
        spx_change_basis(this, i, leavingVarToUpperBound, j);
        if(err != 0) {
            if (err == BFD_ELIMIT) {
                // need to do a refactorization
                int facErr = spx_factorize(this);
                assert(facErr == 0);
                spx_eval_beta(this, beta);
                for(int i=1; i<=nBasic(); ++i) {
                    int kProb = kSimTokProb[head[i]];
                    if(kProb>m) lpSolution[kProb-m] = beta[i];
                }
            } else {
                std::cout << "Unhandled error while pivoting: " << err << std::endl;
                assert(false);
            }
        }
        piIsValid(false);
    }
}


#endif //GLPKTEST_GLPSIMPLEX_H
