#include <iostream>
#include "include/glpkpp.h"

using namespace glp;

void myFunc(SparseVec::EntryRef &rr) {

}

int main() {

    SparseVec mVec;

    for(int i=0; i<10; ++i) {
        int r = rand()%100;
        mVec.add(r, r + 0.5 );
    }

    auto it = mVec.begin();

    auto v = *it;
    std::cout << v.index << std::endl;
    std::cout << *it << std::endl;
    it++;
    std::cout << v.index << std::endl;
    std::cout << it->index << std::endl;
    std::cout << *it << std::endl;

    (*it).value += 0.01;
    std::cout << *it << std::endl;
//    mVec.add(4,4.2);
//    mVec.add(2,2.2);
//    mVec.add(1,1.2);
//    mVec.add(3,3.2);

    mVec.sort();

//    Problem myProb;
//
//    myProb.addConstraint(1.0*X(1) + 1.0*X(2) + 1.0*X(3) <= 100.0);
//    myProb.addConstraint(10.0*X(1) + 4.0*X(2) + 5.0*X(3) <= 600.0);
//    myProb.addConstraint(2.0*X(1) + 2.0*X(2) + 6.0*X(3) <= 300.0);
//    myProb.addConstraint(0.0 <= 1.0*X(1));
//    myProb.addConstraint(0.0 <= 1.0*X(2));
//    myProb.addConstraint(0.0 <= 1.0*X(3));
//    myProb.setObjective(10.0*X(1) + 6.0*X(2) + 4.0*X(3));
//    myProb.stdBasis();
//    myProb.warmUp();
//
//    glp::Simplex mySimplex(myProb);
//
//    std::cout << myProb;
//    std::cout << mySimplex << std::endl;
//
//    mySimplex.pivot(3,3,true);
//
//    std::cout << mySimplex << std::endl;
//
    return 0;
}
