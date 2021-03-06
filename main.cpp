#include <iostream>
#include "include/glpkpp.h"

using namespace glp;

int main() {
    Problem myProb;

    myProb.addConstraint(1.0*X(1) + 1.0*X(2) + 1.0*X(3) <= 100.0);
    myProb.addConstraint(10.0*X(1) + 4.0*X(2) + 5.0*X(3) <= 600.0);
    myProb.addConstraint(2.0*X(1) + 2.0*X(2) + 6.0*X(3) <= 300.0);
    myProb.addConstraint(0.0 <= 1.0*X(1));
    myProb.addConstraint(0.0 <= 1.0*X(2));
    myProb.addConstraint(0.0 <= 1.0*X(3));
    myProb.setObjective(-10.0*X(1) - 6.0*X(2) - 4.0*X(3));
    myProb.setColKind(1,glp::Problem::INTEGER);
    myProb.setColKind(2,glp::Problem::INTEGER);
    myProb.setColKind(3,glp::Problem::INTEGER);

    std::cout << myProb;

    myProb.simplex();
    std::cout << "Optimality status = " << (myProb.getStatus() == Problem::SolutionStatus::OPT) << std::endl;
    std::cout << "LP relaxation optimal:" << SparseVec(myProb.primalSolution()) << std::endl;

//    myProb.intOpt();
//    std::cout << "Branch-and-cut optimal:" << myProb.mipSolution() << std::endl;

    myProb.stdBasis();
    myProb.warmUp();
    glp::Simplex mySimplex(myProb);
    std::cout << mySimplex << std::endl;

    for(int i=1; i<mySimplex.X().size(); ++i) {
        std::cout << mySimplex.X()[i] << std::endl;
    }
    std::cout << "Valid = " << myProb.isValidSolution(mySimplex.X()) << std::endl;


    mySimplex.pivot(3,3);

    std::cout << mySimplex << std::endl;

    for(int i=1; i<mySimplex.X().size(); ++i) {
        std::cout << mySimplex.X()[i] << std::endl;
    }
    std::cout << "Valid = " << myProb.isValidSolution(mySimplex.X()) << std::endl;

    return 0;
}
