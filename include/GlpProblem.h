
//
// Created by daniel on 09/04/2021.
//

#ifndef GLPKTEST_GLPPROBLEM_H
#define GLPKTEST_GLPPROBLEM_H

class GlpProblem {
public:
    glp_prob *lp;

public:
    GlpProblem() {
        lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MIN);
    }

    ~GlpProblem() {
        glp_delete_prob(lp);
    }


    int nConstraints()  { return glp_get_num_rows(lp); }
    int nVars()         { return glp_get_num_cols(lp); } // not including artificial (auxiliary) vars

    void addConstraint(const Constraint &constraint);
    void setObjective(const Constraint &);

    void getObjective(SparseVec &);
    void row(int i, SparseVec &);
    void col(int j, SparseVec &);
    double rowLowerBound(int i);
    double rowUpperBound(int i);
    double colLowerBound(int j);
    double colUpperBound(int j);

    void stdBasis() { glp_std_basis(lp); }
    void warmUp()   { glp_warm_up(lp); }

protected:

    static int glpBoundsType(double lowerBound, double upperBound);

    void ensurenVars(int n);
};

std::ostream &operator <<(std::ostream &out, GlpProblem &prob);

#endif //GLPKTEST_GLPPROBLEM_H
