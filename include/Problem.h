
//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_GLPPROBLEM_H
#define GLPKTEST_GLPPROBLEM_H

class Problem {
public:
    glp_prob *lp;

public:
    Problem() {
        lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MIN);
    }

    Problem(Problem &&moveFrom) noexcept {
        lp = moveFrom.lp;
        moveFrom.lp = NULL;
    }

    ~Problem() {
        if(lp != NULL) glp_delete_prob(lp);
    }


    int nConstraints()  { return glp_get_num_rows(lp); }
    int nVars()         { return glp_get_num_cols(lp); } // not including artificial (auxiliary) vars

    void ensureNVars(int n); // ensures that there are at least n columns in the problem matrix

    void addConstraint(const Constraint &constraint);
    void setObjective(const LinearSum &);
    SparseVec getObjective();

    // glpk interface
    SparseVec getMatRow(int i);
    SparseVec getMatCol(int j);
    int addRows(int n) { return glp_add_rows(lp, n); }
    double getRowLb(int i) const { return glp_get_row_lb(lp, i); }
    double getRowUb(int i) const { return glp_get_row_ub(lp, i); }
    double getColLb(int j) const { return glp_get_col_lb(lp, j); }
    double getColUb(int j) const { return glp_get_col_ub(lp, j); }
    void stdBasis() { glp_std_basis(lp); }
    void warmUp()   { glp_warm_up(lp); }

protected:

    static int glpBoundsType(double lowerBound, double upperBound);

};

std::ostream &operator <<(std::ostream &out, Problem &prob);

#endif //GLPKTEST_GLPPROBLEM_H
