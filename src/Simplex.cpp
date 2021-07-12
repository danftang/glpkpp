//
// Created by daniel on 15/04/2021.
//

#include <iomanip>
#include <cassert>
#include "../include/glpkpp.h"
#include "../include/Simplex.h"


extern "C" {
    #include "spxprob.h"
};

namespace glp {

    Simplex::Simplex(Problem &prob): originalProblem(prob) {
        glp_prob *P = prob.lp;
        beta = new double[prob.nConstraints()];
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
        lpSolutionIsValid(false);
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
// if not present this will be calculated
// if i<1 then column j goes to its opposite bound and does not enter the basis
// assumes there are no fixed variables in the tableau.
    void Simplex::pivot(int i, int j, const std::vector<double> &pivotCol) {
        if(i < 1) { // column j goes to opposite bound
            pivot(i,j,pivotCol, !isAtUpperBound(j));
        } else {
            const int upperBoundFlag = (pivotCol[i] > 0.0) ^isAtUpperBound(j); // leaving variable goes to upper bound?
            pivot(i, j, pivotCol, upperBoundFlag);
        }
    }


    void Simplex::pivot(int i, int j, const std::vector<double> &pivotCol, bool leavingVarToUpperBound) {
        if(i < 1) { // column j goes to opposite bound
            if(isAtUpperBound(j) != leavingVarToUpperBound) {
//                std::cout << "Pivoting to opposite bound" << std::endl;
                spx_update_beta(this, beta, -1, 0, j, pivotCol.data());
                isAtUpperBound(j, leavingVarToUpperBound);
            }
        } else { // basis update
            assert(fabs(pivotCol[i]) > zeroTol);
            assert(u[head[i]] != l[head[i]]);       // TODO: deal with fixed variables leaving the basis
            assert(u[head[m + j]] != l[head[m + j]]);   // fixed variables should never enter the basis
            spx_update_beta(this, beta, i, leavingVarToUpperBound, j, pivotCol.data());
            int err = spx_update_invb(this, i, head[j + m]);
            spx_change_basis(this, i, leavingVarToUpperBound, j);
            if(err != 0) {
                if (err == BFD_ELIMIT) {
                    // need to do a refactorization
                    int facErr = spx_factorize(this);
                    assert(facErr == 0);
                    spx_eval_beta(this, beta);
                } else {
                    std::cout << "Unhandled error while pivoting: " << err << std::endl;
                    assert(false);
                }
            }
            piIsValid(false);
        }
        lpSolutionIsValid(false);
    }


    void Simplex::piIsValid(bool setValid) {
        pi[0] = setValid;
    }

    bool Simplex::piIsValid() {
        return pi[0];
    }

    Simplex::BoundType Simplex::boundType(int k) {
        return BoundType((l[k] == -std::numeric_limits<double>::max()) + 2*(u[k] == std::numeric_limits<double>::max()));
    }

    // Synchronises the state of this simplex with the Problem object originally
    // passed on conctruction
    void Simplex::syncWithLP() {
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
        int nConstraints = originalProblem.nConstraints();
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
        lpSolutionIsValid(true);
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


    std::ostream &operator<<(std::ostream &out, Simplex &simplex) {
        const int nCols = simplex.n - simplex.m;

        out << std::setprecision(5);

        // column labels
        out << std::setw(24) << "\t";
        for (int j = 1; j <= nCols; ++j) {
            out << std::setw(12) << "x" << simplex.head[j + simplex.m] << "\t";
        }
        out << std::setw(12) << "B" << std::endl;

        // col upper limits
        out << std::setw(24) << "\t";
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
            out << std::setw(12) << simplex.b[i] << " <= " << simplex.u[simplex.head[i]] << std::endl;
        }

        // reduced objective
        out << std::endl << std::setw(24) << "\t";
        for (int j = 1; j <= nCols; ++j) {
            out << std::setw(12) << simplex.reducedCost(j) << "\t";
        }
        out << std::endl;

        // col lower limits
        out << std::endl << std::setw(24) << "\t";
        for (int j = 1; j <= nCols; ++j) {
            out << std::setw(12) << simplex.l[simplex.head[j + simplex.m]] << "\t";
        }
        out << std::endl;

        return out;
    }
};
