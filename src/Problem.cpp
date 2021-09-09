//
// Created by daniel on 09/04/2021.
//

#include <iostream>
#include <iomanip>
#include <cfloat>
#include "../include/glpkpp.h"
#include "../include/Problem.h"


namespace glp {

    Problem::Problem(const std::vector<Constraint> &constraints, const SparseVec &objective): Problem() {
        setConstraints(constraints);
        setObjective(objective);
    }


    SparseVec Problem::getObjective() const {
        SparseVec obj;
        obj.reserve(nVars()+1);
        for (int j = 1; j <= nVars(); ++j) {
            obj.insert(j, glp_get_obj_coef(lp, j));
        }
        return obj;
    }


    SparseVec Problem::getMatRow(int i) const {
        SparseVec rowVec(nVars());
        int newSize = glp_get_mat_row(lp, i, rowVec.glpkIndexArray(), rowVec.glpkValueArray());
        rowVec.resize(newSize);
        return rowVec;
    }

    SparseVec Problem::getMatCol(int j) const {
        SparseVec colVec(nConstraints());
        int newSize = glp_get_mat_col(lp, j, colVec.glpkIndexArray(), colVec.glpkValueArray());
        colVec.resize(newSize);
        return colVec;
    }

    void Problem::addConstraint(const Constraint &constraint) {
        if(constraint.coefficients.sparseSize() == 1) { // monomial
            int varId = constraint.coefficients.indices[0];
            double coeff = constraint.coefficients.values[0];
            ensureNVars(varId);
            double lowerBound = std::max(constraint.lowerBound/coeff, getColLb(varId));
            double upperBound = std::min(constraint.upperBound/coeff, getColUb(varId));
            setColBnds(varId, lowerBound, upperBound);
        } else {
            int newRow = glp_add_rows(lp, 1);
            ensureNVars(constraint.coefficients.maxNonZeroIndex());
            glp_set_mat_row(lp, newRow, constraint.coefficients.sparseSize(),
                            constraint.coefficients.glpkIndexArray(),
                            constraint.coefficients.glpkValueArray());
            setRowBnds(newRow, constraint.lowerBound, constraint.upperBound);
        }
    }

    void Problem::setConstraints(const std::vector<Constraint> &constraints) {
        this->eraseProb();
        int nVars = 0;
        for(const Constraint &constraint: constraints) {
            int nMax = constraint.highestVar();
            if(nMax > nVars) nVars = nMax;
        }
        ensureNVars(nVars);
        for(const Constraint &constraint: constraints) addConstraint(constraint);
    }

    Constraint Problem::getConstraint(int i) const {
        return Constraint(getRowLb(i), getMatRow(i), getRowUb(i));
    }

    std::vector<Constraint> Problem::getConstraints() const {
        std::vector<Constraint> constraints;
        constraints.reserve(nConstraints());
        for(int c=1; c<=nConstraints(); ++c) {
            constraints.emplace_back(getRowLb(c), getMatRow(c), getRowUb(c));
        }
        return constraints;
    }

    void Problem::ensureNVars(int n) {
        if(n > nVars()) {
            glp_add_cols(lp, n - nVars());
            for(int j=nVars()-n+1; j<=nVars(); ++j) setColBnds(j,-DBL_MAX, DBL_MAX);
        }
    }

    int Problem::glpBoundsType(double lowerBound, double upperBound) {
        return lowerBound <= -DBL_MAX?
               (upperBound >= DBL_MAX?GLP_FR:GLP_UP):
               (upperBound >= DBL_MAX?
                GLP_LO:(upperBound == lowerBound?GLP_FX:GLP_DB)
               );
    }

    void Problem::setColBnds(int j, double lowerBound, double upperBound) {
        glp_set_col_bnds(lp, j, glpBoundsType(lowerBound, upperBound), lowerBound, upperBound);
    }

    void Problem::setRowBnds(int i, double lowerBound, double upperBound) {
        glp_set_row_bnds(lp, i, glpBoundsType(lowerBound, upperBound), lowerBound, upperBound);
    }


    void Problem::setObjective(const SparseVec &sum) {
        for(int j=1; j<=nVars(); ++j) {
            glp_set_obj_coef(lp, j, 0.0);
        }
        for(int i=0; i<sum.sparseSize(); ++i) {
            glp_set_obj_coef(lp, sum.indices[i], sum.values[i]);
        }
    }


    std::vector<double> Problem::primalSolution() const {
        std::vector<double> solution(nVars()+1);
        for(int j=1; j<nVars(); ++j) {
            solution[j] = glp_get_col_prim(lp, j);
        }
        return solution;
    }

    SparseVec Problem::mipSolution() const {
        SparseVec rowVec;
        rowVec.reserve(nVars()+1);
        for(int j=1; j<=nVars(); ++j) {
            rowVec.insert(j, glp_mip_col_val(lp, j));
        }
        return rowVec;
    }

    bool Problem::isValidSolution(const std::vector<double> &X) const {
        constexpr double tol = 1e-8;
        for(int j=1; j < X.size(); ++j) {
            double Xj = X[j];
            if(getColLb(j) - Xj > tol || Xj - getColUb(j) > tol) return false;
        }
        for(int i=1; i<=nConstraints(); ++i) {
            double Bi = getMatRow(i) * X;
            if(getRowLb(i) - Bi > tol || Bi - getRowUb(i) > tol) return false;
        }
        return true;
    }


    // find a basis that has the structural variables given in "solution" by setting the objective
    // function to 1.0 for each structural on its lower bound and -1.0 for each structural
    // on its upper bound and minimising.
    void Problem::solutionBasis(const std::vector<double> &solution) {
        SparseVec originalObjective = getObjective();
        ObjectiveDirection originalDirection = getObjDir();
        SparseVec solutionObjective;
        solutionObjective.reserve(nVars());
        for(int j=1; j<solution.size(); ++j) {
            if(solution[j] == getColLb(j)) {
                solutionObjective.insert(j, 1.0);
            } else if(solution[j] == getColUb(j)) {
                solutionObjective.insert(j, -1.0);
            }
        }
        setObjDir(MINIMISE);
        setObjective(solutionObjective);
        advBasis();
        simplex();
        setObjDir(originalDirection);
        setObjective(originalObjective);
    }



    std::ostream &operator<<(std::ostream &out, const Problem &prob) {
        out << "Problem has " << prob.nConstraints() << " constraints and " << prob.nVars() << " variables." << std::endl;
        switch (glp_get_obj_dir(prob.lp)) {
            case GLP_MIN:
                out << "Minimise ";
                break;
            case GLP_MAX:
                out << "Maximise ";
                break;
            default:
                out << "Unknown objective type ";
        }
        out << prob.getObjective() << std::endl;
        out << "Subject to:" << std::endl;

//        std::vector<double> row(prob.nVars()+1);
        for (int i = 1; i <= prob.nConstraints(); ++i) {
            out << std::setw(13) << prob.getRowLb(i) << " <= ";
            std::vector<double> row = prob.getMatRow(i).toDense(prob.nVars()+1);
            for(int j=1; j < row.size(); ++j) {
                out << std::setw(5) << row[j];
            }
            out << " <= " << prob.getRowUb(i) << std::endl;
        }
        for (int j = 1; j <= prob.nVars(); ++j) {
            out << prob.getColLb(j) << " <= x[" << j << "] <= " << prob.getColUb(j) << std::endl;
        }
        return out;
    }

};