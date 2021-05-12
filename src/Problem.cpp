//
// Created by daniel on 09/04/2021.
//

#include <iostream>
#include <iomanip>
#include "../include/glpkpp.h"

namespace glp {

    SparseVec Problem::getObjective() {
        SparseVec obj;
        obj.reserve(nVars());
        int nVars = glp_get_num_cols(lp);
        double c;
        for (int j = 1; j <= nVars; ++j) {
            c = glp_get_obj_coef(lp, j);
            if (c != 0.0) obj.add(j, c);
        }
        return obj;
    }


    SparseVec Problem::getMatRow(int i) const {
        SparseVec rowVec(nVars());
        rowVec.resize(nVars());
        int newSize = glp_get_mat_row(lp, i, rowVec.glpkIndexArray(), rowVec.glpkValueArray());
        rowVec.resize(newSize);
        return rowVec;
    }

    SparseVec Problem::getMatCol(int j) const {
        SparseVec colVec(nConstraints());
        colVec.resize(nConstraints());
        int newSize = glp_get_mat_col(lp, j, colVec.glpkIndexArray(), colVec.glpkValueArray());
        colVec.resize(newSize);
        return colVec;
    }

    void Problem::addConstraint(const Constraint &constraint) {
        if(constraint.coefficients.sparseSize() == 1) { // monomial
            int varId = constraint.coefficients.indices[0];
            double coeff = constraint.coefficients.values[0];
            ensureNVars(varId);
            glp_set_col_bnds(lp, varId, constraint.glpBoundType(),
                             constraint.lowerBound/coeff,
                             constraint.upperBound/coeff
                             );
        } else {
            int newRow = glp_add_rows(lp, 1);
//            SparseVec sparseRow(nNCols(), constraint.coefficients);
            ensureNVars(constraint.coefficients.maxNonZeroIndex());
            glp_set_mat_row(lp, newRow, constraint.coefficients.sparseSize(),
                            constraint.coefficients.glpkIndexArray(),
                            constraint.coefficients.glpkValueArray());
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
        for(int i=0; i<sum.sparseSize(); ++i) {
            glp_set_obj_coef(lp, sum.indices[i], sum.values[i]);
        }
    }


    SparseVec Problem::primalSolution() const {
        SparseVec rowVec;
        rowVec.reserve(nVars());
        for(int j=1; j<nVars(); ++j) {
            rowVec.add(j, glp_get_col_prim(lp, j));
        }
        return rowVec;
    }

    SparseVec Problem::mipSolution() const {
        SparseVec rowVec;
        rowVec.reserve(nVars());
        for(int j=1; j<nVars(); ++j) {
            rowVec.add(j, glp_mip_col_val(lp, j));
        }
        return rowVec;
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

        std::vector<double> row(prob.nVars());
        for (int i = 1; i <= prob.nConstraints(); ++i) {
            out << std::setw(13) << prob.getRowLb(i) << " <= ";
            prob.getMatRow(i).toDense(row.data(), row.size());
            for(int j=0; j < row.size(); ++j) {
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