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

    Simplex::Simplex(const Problem &prob) {
        glp_prob *P = prob.lp;
        spx_init_lp(this, P, 1);
        spx_alloc_lp(this);
        pi.resize(n+1);
        kProbTokSim.resize(prob.nVars() + prob.nConstraints() + 1);
        kSimTokProb.resize(n+1);
        lpSolution.resize(prob.nVars()+1);
        spx_build_lp(this, P, 1, shift, kProbTokSim.data());
        spx_build_basis(this, P, kProbTokSim.data());
        spx_eval_beta(this,b);
        pi[0] = 0.0;    // indicates not yet evaluated
        for(int kProb=1; kProb <= prob.nVars() + prob.nConstraints(); ++kProb) {
            int kSim = kProbTokSim[kProb];
            if(kSim > 0) {
                kSimTokProb[kSim] = kProb;
            } else if(kSim == 0) {
                // deleted fixed variable
                int j = kProb - prob.nConstraints();
                lpSolution[j] = prob.getColUb(j);
            } else {
                assert(false); // don't support shifted bounds as yet (loss of information on shift).
            }
        }
        lpSolutionIsValid(false);
    }


    Simplex::~Simplex() {
        spx_free_lp(this);
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

    double Simplex::reducedObjective(int j) {
        if (!piIsValid()) {
            spx_eval_pi(this, pi.data());
            piIsValid(true);
        }
        return spx_eval_dj(this, pi.data(), j);
    }

// i,j is the pivot element, 1-based, (the j'th col of N replaces the i'th col of B)
// pivotCol - should be set to the current ftran of the incoming column,
// if not present this will be calculated
    void Simplex::pivot(int i, int j, const std::vector<double> &pivotCol) {
        assert(pivotCol[i] != 0.0);
        const int upperBoundFlag = (pivotCol[i] > 0.0) ^ isAtUpperBound(j);  // leaving variable goes to upper bound?;
        spx_update_beta(this, b, i, upperBoundFlag, j, pivotCol.data());
        spx_update_invb(this, i, head[j + m]);
        spx_change_basis(this, i, upperBoundFlag, j);
        piIsValid(false);
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
    void Simplex::syncWith(Problem &prob) {
        std::vector<int> daeh(n+1); // inverse of head[]
        spx_store_basis(this, prob.lp, kProbTokSim.data(), daeh.data());
        std::vector<double> d(n-m+1); // reduced objective
        for(int j=1; j<=n-m; ++j) { d[j] = reducedObjective(j); }
        spx_store_sol(this, prob.lp, shift, kProbTokSim.data(), daeh.data(), b, pi.data(), d.data());
    }

    const std::vector<double> &Simplex::X() {
        if(!lpSolutionIsValid()) calculateLpSolution();
        return lpSolution;
    }

    // assumes no free variables and m == original number of constraints
    void Simplex::calculateLpSolution() {
        int kProb,kSim;
        for(int i=1; i<=m; ++i) {
            kProb = kSimTokProb[head[i]];
            if(kProb > m) {
                lpSolution[kProb-m] = b[i];
            }
        }
        for(int j=1; j <= n-m; ++j) {
            kSim = head[m+j];
            kProb = kSimTokProb[kSim];
            if(kProb > m) {
                if(isAtUpperBound(j)) {
                    lpSolution[kProb-m] = u[kSim];
                } else {
                    lpSolution[kProb-m] = l[kSim];
                }
            }
        }
        lpSolutionIsValid(true);
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
            out << std::setw(12) << simplex.reducedObjective(j) << "\t";
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
