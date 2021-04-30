//
// Created by daniel on 09/04/2021.
//

#include <iostream>
#include "../include/glpkpp.h"

namespace glp {

    SparseVec Problem::getObjective() {
        SparseVec obj(nVars());
        int nVars = glp_get_num_cols(lp);
        double c;
        for (int j = 1; j <= nVars; ++j) {
            c = glp_get_obj_coef(lp, j);
            if (c != 0.0) obj.add(j, c);
        }
        return obj;
    }


    SparseVec Problem::getMatRow(int i) {
        SparseVec rowVec(nVars());
        rowVec.sparseSize() = glp_get_mat_row(lp, i, rowVec.indices, rowVec.values);
        return rowVec;
    }

    SparseVec Problem::getMatCol(int j) {
        SparseVec colVec(nConstraints());
        colVec.sparseSize() = glp_get_mat_col(lp, j, colVec.indices, colVec.values);
        return colVec;
    }

    void Problem::addConstraint(const Constraint &constraint) {
        if(constraint.coefficients.size() == 1) { // monomial
            auto entry = *constraint.coefficients.begin();
            ensureNVars(entry.first);
            glp_set_col_bnds(lp, entry.first, constraint.glpBoundType(),
                             constraint.lowerBound/entry.second,
                             constraint.upperBound/entry.second
                             );
        } else {
            int newRow = glp_add_rows(lp, 1);
            SparseVec sparseRow(nVars(), constraint.coefficients);
            ensureNVars(sparseRow.maxNonZeroIndex());
            glp_set_mat_row(lp, newRow, sparseRow.sparseSize(), sparseRow.indices, sparseRow.values);
            glp_set_row_bnds(lp, newRow, constraint.glpBoundType(), constraint.lowerBound, constraint.upperBound);
        }
    }

    void Problem::ensureNVars(int n) {
        if(n > nVars()) {
            glp_add_cols(lp, n - nVars());
        }
    }

    int Problem::glpBoundsType(double lowerBound, double upperBound) {
        return lowerBound == -std::numeric_limits<double>::infinity()?
               (upperBound == std::numeric_limits<double>::infinity()?GLP_FR:GLP_UP):
               (upperBound == std::numeric_limits<double>::infinity()?
                GLP_LO:(upperBound == lowerBound?GLP_FX:GLP_DB)
               );
    }

    void Problem::setObjective(const LinearSum &sum) {
        for(auto entry: sum) {
            glp_set_obj_coef(lp, entry.second, entry.first);
        }
    }

    std::ostream &operator<<(std::ostream &out, Problem &prob) {
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
        for (int i = 1; i <= prob.nConstraints(); ++i) {
            out << prob.getRowLb(i) << " <= " << prob.getMatRow(i) << " <= " << prob.getRowUb(i) << std::endl;
        }
        for (int j = 1; j <= prob.nVars(); ++j) {
            out << prob.getColLb(j) << " <= x[" << j << "] <= " << prob.getColUb(j) << std::endl;
        }
        return out;
    }

};