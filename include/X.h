//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKPP_X_H
#define GLPKPP_X_H

// Represents the id'th variable in a LinearSum
class X {
public:
    int id;
    X(int id): id(id) { }

    operator int() const { return id; }
};


#endif //GLPKPP_X_H
