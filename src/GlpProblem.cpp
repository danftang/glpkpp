//
// Created by daniel on 09/04/2021.
//

#include <iostream>
#include "../include/glpkpp.h"

namespace glpkpp {

    void GlpProblem::getObjective(SparseVec &obj) {
        obj.clear();
        int nVars = glp_get_num_cols(lp);
        double c;
        for (int j = 1; j <= nVars; ++j) {
            c = glp_get_obj_coef(lp, j);
            if (c != 0.0) obj.add(j, c);
        }
    }


    void GlpProblem::row(int i, SparseVec &rowVec) {
        rowVec.nnz = glp_get_mat_row(lp, i, rowVec.indices, rowVec.values);
    }

    void GlpProblem::col(int j, SparseVec &colVec) {
        colVec.nnz = glp_get_mat_col(lp, j, colVec.indices, colVec.values);
    }

    double GlpProblem::rowLowerBound(int i) {
        return glp_get_row_lb(lp, i);
    }

    double GlpProblem::rowUpperBound(int i) {
        return glp_get_row_ub(lp, i);
    }

    double GlpProblem::colLowerBound(int j) {
        return glp_get_col_lb(lp, j);
    }

    double GlpProblem::colUpperBound(int j) {
        return glp_get_col_ub(lp, j);
    }

    void GlpProblem::addConstraint(const Constraint &constraint) {
        if(constraint.coefficients.size() == 1) { // monomial
            auto entry = *constraint.coefficients.begin();
            ensurenVars(entry.first);
            glp_set_col_bnds(lp, entry.first, constraint.glpBoundType(),
                             constraint.lowerBound/entry.second,
                             constraint.upperBound/entry.second
                             );
        } else {
            int newRow = glp_add_rows(lp, 1);
            SparseVec sparseRow(nVars(), constraint.coefficients);
            ensurenVars(sparseRow.maxNonZeroIndex());
            glp_set_mat_row(lp, newRow, sparseRow.nnz, sparseRow.indices, sparseRow.values);
            glp_set_row_bnds(lp, newRow, constraint.glpBoundType(), constraint.lowerBound, constraint.upperBound);
        }
    }

    void GlpProblem::ensurenVars(int n) {
        if(n > nVars()) {
            glp_add_cols(lp, n - nVars());
        }
    }

    int GlpProblem::glpBoundsType(double lowerBound, double upperBound) {
        return lowerBound == -std::numeric_limits<double>::infinity()?
               (upperBound == std::numeric_limits<double>::infinity()?GLP_FR:GLP_UP):
               (upperBound == std::numeric_limits<double>::infinity()?
                GLP_LO:(upperBound == lowerBound?GLP_FX:GLP_DB)
               );
    }

    void GlpProblem::setObjective(const Constraint &obj) {
        for(auto entry: obj.coefficients) {
            glp_set_obj_coef(lp, entry.first, entry.second);
        }
    }

    std::ostream &operator<<(std::ostream &out, GlpProblem &prob) {
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
        SparseVec rowVec(prob.nVars());
        prob.getObjective(rowVec);
        out << rowVec << std::endl;
        out << "Subject to:" << std::endl;
        for (int i = 1; i <= prob.nConstraints(); ++i) {
            prob.row(i, rowVec);
            out << prob.rowLowerBound(i) << " <= " << rowVec << " <= " << prob.rowUpperBound(i) << std::endl;
        }
        for (int j = 1; j <= prob.nVars(); ++j) {
            out << prob.colLowerBound(j) << " <= x[" << j << "] <= " << prob.colUpperBound(j) << std::endl;
        }
        return out;
    }

};