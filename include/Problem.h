
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

    enum SolutionStatus {
        OPT = GLP_OPT,
        FEAS = GLP_FEAS,
        INFEAS = GLP_INFEAS,
        NOFEAS = GLP_NOFEAS,
        UNBND = GLP_UNBND,
        UNDEF = GLP_UNDEF
    };

    Problem() {
        lp = glp_create_prob();
        setObjDir(MINIMISE);
    }

    Problem(Problem &&moveFrom) noexcept {
        lp = moveFrom.lp;
        moveFrom.lp = NULL;
    }

    Problem(const std::vector<Constraint> &constraints, const SparseVec &objective = SparseVec());

    ~Problem() {
        if(lp != NULL) glp_delete_prob(lp);
    }


    int nConstraints() const { return glp_get_num_rows(lp); }
    int nVars() const        { return glp_get_num_cols(lp); } // not including artificial (auxiliary) vars

    void ensureNVars(int n); // ensures that there are at least n columns in the problem matrix

    // constraint stuff
    void addConstraint(const Constraint &constraint);
    void setConstraints(const std::vector<Constraint> &constraints);
    Constraint getConstraint(int i) const;
    std::vector<Constraint> getConstraints() const;

    // Objective stuff
    void setObjective(const SparseVec &);
    void setObjective(const std::vector<double> &);
    void setObjective(const LinearSum &sum) { setObjective(sum.toSparseVec()); }
    SparseVec getObjective() const;
    void setObjDir(ObjectiveDirection direction) { glp_set_obj_dir(lp, direction); }
    ObjectiveDirection getObjDir() const { return ObjectiveDirection(glp_get_obj_dir(lp)); }

    // glpk interface
    SparseVec getMatRow(int i) const;
    SparseVec getMatCol(int j) const;
    int addRows(int n) { return glp_add_rows(lp, n); }

    double getRowLb(int i) const { return glp_get_row_lb(lp, i); }
    double getRowUb(int i) const { return glp_get_row_ub(lp, i); }
    double getColLb(int j) const { return glp_get_col_lb(lp, j); }
    double getColUb(int j) const { return glp_get_col_ub(lp, j); }
    int getRowType(int i) const { return glp_get_row_type(lp, i); }
    int getColType(int j) const { return glp_get_col_type(lp, j); }
    int getRowStat(int i) const { return glp_get_row_stat(lp, i); }
    int getColStat(int j) const { return glp_get_col_stat(lp, j); }
    void setColBnds(int j, double lowerBound, double UpperBound);
    void setRowBnds(int i, double lowerBound, double UpperBound);
    void setRowStat(int i, int state) { glp_set_row_stat(lp, i, state); }
    void setColStat(int i, int state) { glp_set_col_stat(lp, i, state); }
    void eraseProb() { glp_erase_prob(lp); }

    SolutionStatus getStatus() const { return SolutionStatus(glp_get_status(lp)); }
    void stdBasis() { glp_std_basis(lp); }
    void advBasis() { glp_adv_basis(lp, 0); }
    void cpxBasis() { glp_cpx_basis(lp); }
    void solutionBasis(const std::vector<double> &);
    void warmUp()   { glp_warm_up(lp); }

    // LP stuff
    void simplex(const glp_smcp *parm=NULL) { glp_simplex(lp, parm); } // LP-solve
    std::vector<double> primalSolution() const;

    // MIP stuff
    void intOpt(const glp_iocp *parm=NULL) { glp_intopt(lp,parm); }
    void setColKind(int col, VariableKind varKind) { glp_set_col_kind(lp, col, varKind); }
    SparseVec mipSolution() const;

    // Solution stuff
    bool isValidSolution(const std::vector<double> &X) const; // X in 1..nVars() i.e. only structural vars

protected:

    static int glpBoundsType(double lowerBound, double upperBound);

};

std::ostream &operator <<(std::ostream &out, const Problem &prob);

#endif //GLPKTEST_GLPPROBLEM_H
