//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKPP_GLPKPP_H
#define GLPKPP_GLPKPP_H

#include <ostream>
#include <cstdlib>
#include <ostream>
#include <map>
#include <functional>
#include <limits>
#include <vector>

#include "glpk.h"
extern "C" {
    #include "spxlp.h"
};


namespace glpkpp {
    #include "SparseVec.h"
    #include "X.h"
    #include "Constraint.h"
    #include "LinearSum.h"
    #include "GlpProblem.h"
    #include "GlpSimplex.h"
};

#endif //GLPKPP_GLPKPP_H
