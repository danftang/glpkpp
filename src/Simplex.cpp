//
// Created by daniel on 15/04/2021.
//

#include <iomanip>
#include "../include/glpkpp.h"

extern "C" {
    #include "spxprob.h"
};

namespace glp {

    Simplex::Simplex(const Problem &prob) {
        glp_prob *P = prob.lp;
        beta = new double[prob.nConstraints()+1];
        spx_init_lp(this, P, excludeFixed);
        spx_alloc_lp(this);
        pi.resize(nBasic()+1);
        kProbTokSim.resize(prob.nVars() + prob.nConstraints() + 1);
        kSimTokProb.resize(n+1);
        lpSolution.resize(prob.nVars()+1);
        spx_build_lp(this, P, excludeFixed, shift, kProbTokSim.data());
        spx_build_basis(this, P, kProbTokSim.data());
        spx_eval_beta(this, beta);
        assert(m == prob.nConstraints()); // TODO: weaken this constraint (needed for conversion back to original space)
        pi[0] = 0.0;    // indicates not yet evaluated
        lpSolution[0] = 0.0;
        for(int kProb=1; kProb <= prob.nVars() + prob.nConstraints(); ++kProb) {
            int kSim = kProbTokSim[kProb];
            if(kSim > 0) {
                kSimTokProb[kSim] = kProb;
            } else if(kSim == 0) {
                // deleted fixed variable
                int j = kProb - prob.nConstraints();
                if(j>0) lpSolution[j] = prob.getColUb(j);
            } else {
                assert(false); // don't support shifted bounds as yet (loss of information on shift).
            }
        }
        calculateLpSolution();
    }

    Simplex::Simplex(Simplex &&moveFrom):
    kProbTokSim(std::move(moveFrom.kProbTokSim)),
    kSimTokProb(std::move(moveFrom.kSimTokProb)),
    pi(std::move(moveFrom.pi)),
    lpSolution(std::move(moveFrom.lpSolution))
    {
        static_cast<SPXLP &>(*this) = static_cast<SPXLP &>(moveFrom);
        moveFrom.n=0;
        moveFrom.m =0;
        moveFrom.nnz=0;
        spx_alloc_lp(&moveFrom);
        beta = moveFrom.beta;
        moveFrom.beta = new double[1];
    }


    Simplex::~Simplex() {
        spx_free_lp(this);
        delete [] beta;
    }

    std::vector<double> Simplex::tableauRow(int i) {
        std::vector<double> trow(n-m+1);
        double rho[m + 1];
        spx_eval_rho(this, i, rho);
        spx_eval_trow(this, rho, trow.data());
        return trow;
    }

    std::vector<double> Simplex::tableauCol(int j) {
        std::vector<double> denseCol(m+1);
        spx_eval_tcol(this, j, denseCol.data());
        return denseCol;
    }

    double Simplex::reducedCost(int j) {
        if (!piIsValid()) recalculatePi();
        return spx_eval_dj(this, pi.data(), j);
    }

// i,j is the pivot element, 1-based, (the j'th col of N replaces the i'th col of B)
// pivotCol - should be set to the current ftran of the incoming column,
// if i<1 then column j goes to its opposite bound and does not enter the basis
// otherwise, the outgoing variable goes to the bound that moves the incoming variable
// in the feasible direction.
// assumes there are no fixed variables in the tableau.
    void Simplex::pivot(int i, int j, const std::vector<double> &pivotCol) {
        if(i < 1) { // column j goes to opposite bound
            pivot(i,j,pivotCol, !isAtUpperBound(j));
        } else {
            const int upperBoundFlag = (pivotCol[i] > 0.0) ^ isAtUpperBound(j); // move incoming variable in feasible direction
            pivot(i, j, pivotCol, upperBoundFlag);
        }
    }


//    void Simplex::pivot(int i, int j, const std::vector<double> &pivotCol, bool leavingVarToUpperBound) {
//        if(i < 1) { // column j goes to opposite bound
//            if(isAtUpperBound(j) != leavingVarToUpperBound) {
//                updateBetaAndLPSolution(-1, 0, j, pivotCol);
//                isAtUpperBound(j, leavingVarToUpperBound);
//            }
//        } else { // basis update
//            int enteringk = head[m+j];
//            assert(fabs(pivotCol[i]) > zeroTol);
//            assert(u[head[i]] != l[head[i]]);       // TODO: deal with fixed variables leaving the basis
//            assert(u[enteringk] != l[enteringk]);   // fixed variables should never enter the basis
//            updateBetaAndLPSolution(i, leavingVarToUpperBound, j, pivotCol);
//            int err = spx_update_invb(this, i, enteringk);
//            spx_change_basis(this, i, leavingVarToUpperBound, j);
//            if(err != 0) {
//                if (err == BFD_ELIMIT) {
//                    // need to do a refactorization
//                    int facErr = spx_factorize(this);
//                    assert(facErr == 0);
//                    spx_eval_beta(this, beta);
//                } else {
//                    std::cout << "Unhandled error while pivoting: " << err << std::endl;
//                    assert(false);
//                }
//            }
//            piIsValid(false);
//        }
//    }
//
//
//    void Simplex::pivot(int i, int j, const SparseVec &spPivotCol, bool leavingVarToUpperBound) {
//        if(i < 1) { // column j goes to opposite bound
//            if(isAtUpperBound(j) != leavingVarToUpperBound) {
//                updateBetaAndLPSolution( -1, 0, j, spPivotCol);
//                isAtUpperBound(j, leavingVarToUpperBound);
//            }
//        } else { // basis update
//            int enteringk = head[m+j];
//            assert(u[head[i]] != l[head[i]]);       // TODO: deal with fixed variables leaving the basis
//            assert(u[enteringk] != l[enteringk]);   // fixed variables should never enter the basis
//            updateBetaAndLPSolution( i, leavingVarToUpperBound, j, spPivotCol);
//            int err = spx_update_invb(this, i, enteringk);
//            spx_change_basis(this, i, leavingVarToUpperBound, j);
//            if(err != 0) {
//                if (err == BFD_ELIMIT) {
//                    // need to do a refactorization
//                    int facErr = spx_factorize(this);
//                    assert(facErr == 0);
//                    spx_eval_beta(this, beta);
//                } else {
//                    std::cout << "Unhandled error while pivoting: " << err << std::endl;
//                    assert(false);
//                }
//            }
//            piIsValid(false);
//        }
//    }

    // p is leaving row, q is entering col
    void Simplex::updateBetaAndLPSolution(int p, double delta_q, const SparseVec &pivotCol) {
        /* compute new beta[i] for all i != p */
        int i, kProb;
        for (int z = 0; z < pivotCol.sparseSize(); z++) {
            i = pivotCol.indices[z];
            if (i != p) {
                beta[i] += pivotCol.values[z] * delta_q;
                kProb = kSimTokProb[head[i]];
                if(kProb > m) lpSolution[kProb-m] = beta[i];
            }
        }
    }

    void Simplex::updateBetaAndLPSolution(int p, double delta_q, const std::vector<double> &pivotCol) {
        /* compute new beta[i] for all i != p */
        int i, kProb;
        for (i = 1; i <= m; i++) {
            if (i != p) {
                beta[i] += pivotCol[i] * delta_q;
                kProb = kSimTokProb[head[i]];
                if (kProb > m) lpSolution[kProb - m] = beta[i];
            }
        }
    }

    void Simplex::updateBetaAndLPSolution(int p, double delta_q, const FVSVector &pivotCol) {
        /* compute new beta[i] for all i != p */
        int i, kProb;
        for (int z = 0; z < pivotCol.sparseSize(); z++) {
            i = pivotCol.indices[z];
            if (i != p) {
                beta[i] += pivotCol.vec[i] * delta_q;
                kProb = kSimTokProb[head[i]];
                if(kProb > m) lpSolution[kProb-m] = beta[i];
            }
        }
    }


    void Simplex::piIsValid(bool setValid) {
        pi[0] = setValid;
    }

    bool Simplex::piIsValid() const {
        return pi[0];
    }

    Simplex::BoundType Simplex::boundType(int k) {
        return BoundType((l[k] == -std::numeric_limits<double>::max()) + 2*(u[k] == std::numeric_limits<double>::max()));
    }

    // Synchronises the state of this simplex with the Problem object originally
    // passed on conctruction
    void Simplex::syncWithLP(Problem &originalProblem) {
        std::vector<int> daeh(n+1); // inverse of head[]
        spx_store_basis(this, originalProblem.lp, kProbTokSim.data(), daeh.data());
        std::vector<double> d(n-m+1); // reduced objective
        for(int j=1; j<=n-m; ++j) { d[j] = reducedCost(j); }
        spx_store_sol(this, originalProblem.lp, shift, kProbTokSim.data(), daeh.data(), beta, pi.data(), d.data());
    }

    const std::vector<double> &Simplex::X() {
        if(!lpSolutionIsValid()) calculateLpSolution();
        return lpSolution;
    }


    void Simplex::calculateLpSolution() {
        int kProb,kSim;
        int nConstraints = nBasic();
        for(int i=1; i<=m; ++i) {
            kProb = kSimTokProb[head[i]];
            if(kProb > nConstraints) {
                lpSolution[kProb-nConstraints] = beta[i];
            }
        }
        for(int j=1; j <= n-m; ++j) {
            kSim = head[m+j];
            kProb = kSimTokProb[kSim];
            if(kProb > nConstraints) {
                if(isAtUpperBound(j)) {
                    lpSolution[kProb-nConstraints] = u[kSim];
                } else {
                    lpSolution[kProb-nConstraints] = l[kSim];
                }
            }
        }
    }

    std::vector<double> Simplex::reducedCost() {
        std::vector<double> reducedObjective(nNonBasic() + 1);
        if (!piIsValid()) recalculatePi();
        for(int j=1; j<=nNonBasic(); ++j) {
            reducedObjective[j] = spx_eval_dj(this, pi.data(), j);
        }
        return reducedObjective;
    }

    void Simplex::recalculatePi() {
        spx_eval_pi(this, pi.data());
        piIsValid(true);
    }


    void Simplex::setObjective(const SparseVec &objective) {
        for(int k=0; k <= nVars(); ++k) {
            c[k] = 0.0;
        }
        for(int l=0; l<objective.sparseSize(); ++l) {
            c[objective.indices[l]] = objective.values[l];
        }
    }

    std::vector<double> Simplex::piTimesMinusN(const std::vector<double> &pi) {
        std::vector<double> d(nNonBasic()+1, 0.0); // = pi*N
        int k, ptr, end;
        for(int j=1; j<=nNonBasic(); ++j) {
            k = head[m + j]; /* x[k] = xN[j] */
            /* dj := c[k] */
            /* dj := dj - A'[k] * pi */
            ptr = A_ptr[k];
            end = A_ptr[k + 1];
            for (; ptr < end; ptr++)
                d[j] -= A_val[ptr] * pi[A_ind[ptr]];
        }
        return d;
    }



    std::ostream &operator<<(std::ostream &out, Simplex &simplex) {
        const int nCols = simplex.n - simplex.m;


        out << std::setprecision(5);

        // column labels
        out << std::setw(24) << "\t";
        for (int j = 1; j <= nCols; ++j) {
            out << std::setw(12) << "x" << simplex.head[j + simplex.m] << "\t";
        }
        out << std::setw(12) << "Beta" << std::endl;

        // col upper limits
        out << std::setw(22) << "Upper bound = \t";
        for (int j = 1; j <= nCols; ++j) {
            out << std::setw(12) << simplex.u[simplex.head[j + simplex.m]] << "\t";
        }
        out << std::endl;

        // cols on upper limits
        out << std::setw(24) << "\t";
        for (int j = 1; j <= nCols; ++j) {
            out << std::setw(12) << (simplex.isAtUpperBound(j) ? "--------":"\t") << "\t";
        }
        out << std::endl;

        // coefficients
        for (int i = 1; i <= simplex.m; ++i) {
            out << std::setw(12) << simplex.l[simplex.head[i]] << " <= x" << simplex.head[i] << " = \t";
            std::vector<double> row = simplex.tableauRow(i);
            for (int j = 1; j <= nCols; ++j) {
                out << std::setw(12) << row[j] << "\t";
            }
            out << std::setw(12) << simplex.beta[i] << " <= " << simplex.u[simplex.head[i]] << std::endl;
        }

        // reduced objective
        out << std::endl << std::setw(22) << "z = \t";
        for (int j = 1; j <= nCols; ++j) {
            out << std::setw(12) << simplex.reducedCost(j) << "\t";
        }
        out << std::endl;

        // col lower limits
        out << std::endl << std::setw(22) << "Lower bound = \t";
        for (int j = 1; j <= nCols; ++j) {
            out << std::setw(12) << simplex.l[simplex.head[j + simplex.m]] << "\t";
        }
        out << std::endl;
        // cols on lower limits
        out << std::setw(24) << "\t";
        for (int j = 1; j <= nCols; ++j) {
            out << std::setw(12) << (simplex.isAtUpperBound(j) ? "\t":"--------") << "\t";
        }
        out << std::endl;

        return out;
    }
};
