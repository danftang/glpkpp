
//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_GLPPROBLEM_H
#define GLPKTEST_GLPPROBLEM_H

class Problem {
public:
    glp_prob *lp;

public:

    enum ObjectiveDirection {
        MAXIMISE = GLP_MAX,
        MINIMISE = GLP_MIN
    };

    enum VariableKind {
        CONTINUOUS = GLP_CV,
        INTEGER = GLP_IV,
        BINARY = GLP_BV
    };

    Problem() {
        lp = glp_create_prob();
        setObjDir(MINIMISE);
    }

    Problem(Problem &&moveFrom) noexcept {
        lp = moveFrom.lp;
        moveFrom.lp = NULL;
    }

    ~Problem() {
        if(lp != NULL) glp_delete_prob(lp);
    }


    int nConstraints() const { return glp_get_num_rows(lp); }
    int nVars() const        { return glp_get_num_cols(lp); } // not including artificial (auxiliary) vars

    void ensureNVars(int n); // ensures that there are at least n columns in the problem matrix

    void addConstraint(const Constraint &constraint);

    // Objective stuff
    void setObjective(const LinearSum &);
    SparseVec getObjective();
    void setObjDir(ObjectiveDirection direction) { glp_set_obj_dir(lp, direction); }

    // glpk interface
    SparseVec getMatRow(int i) const;
    SparseVec getMatCol(int j) const;
    int addRows(int n) { return glp_add_rows(lp, n); }
    double getRowLb(int i) const { return glp_get_row_lb(lp, i); }
    double getRowUb(int i) const { return glp_get_row_ub(lp, i); }
    double getColLb(int j) const { return glp_get_col_lb(lp, j); }
    double getColUb(int j) const { return glp_get_col_ub(lp, j); }
    void stdBasis() { glp_std_basis(lp); }
    void advBasis() { glp_adv_basis(lp, 0); }
    void cpxBasis() { glp_cpx_basis(lp); }
    void warmUp()   { glp_warm_up(lp); }

    // LP stuff
    void simplex(const glp_smcp *parm=NULL) { glp_simplex(lp, parm); } // LP-solve
    SparseVec primalSolution() const;

    // MIP stuff
    void intOpt(const glp_iocp *parm=NULL) { glp_intopt(lp,parm); }
    void setColKind(int col, VariableKind varKind) { glp_set_col_kind(lp, col, varKind); }
    SparseVec mipSolution() const;

    bool isValidSolution(const std::vector<double> &X);

protected:

    static int glpBoundsType(double lowerBound, double upperBound);

};

std::ostream &operator <<(std::ostream &out, Problem &prob);

#endif //GLPKTEST_GLPPROBLEM_H
