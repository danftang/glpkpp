//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKPP_GLPKPP_H
#define GLPKPP_GLPKPP_H

#include <ostream>
#include <cstdlib>
#include <iostream>
#include <map>
#include <functional>
#include <limits>
#include <vector>
#include <algorithm>

#include "glpk.h"
extern "C" {
    #include "spxlp.h"
};


namespace glp {
    #include "SparseVec.h"
    #include "X.h"
    #include "LinearSum.h"
    #include "Constraint.h"
    #include "Problem.h"
    #include "Simplex.h"
};

#endif //GLPKPP_GLPKPP_H
