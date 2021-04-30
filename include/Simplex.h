//
// Created by daniel on 15/04/2021.
//

#ifndef GLPKTEST_GLPSIMPLEX_H
#define GLPKTEST_GLPSIMPLEX_H


class Simplex: public SPXLP {
public:
    int *map;       // map of variables in this simplex to those of the original problem
    double *pi;     // reduced objective = c_B*B'*N + c_N = pi*N + c_N

    Simplex(const Problem &prob);
    ~Simplex();


    void tableauRow(int i, double row[]);
    void tableauCol(int j, double col[]);
    double reducedObjective(int j);             // value of the j'th (1 <= j <= n-m) element of the reduced objective
    void pivot(int i, int j, bool leavingVarToUpperBound, double *pivotCol=NULL);


};

std::ostream &operator <<(std::ostream &out, Simplex &simplex);


#endif //GLPKTEST_GLPSIMPLEX_H
