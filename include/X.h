//
// Created by daniel on 26/04/2021.
//

#ifndef GLPKPP_X_H
#define GLPKPP_X_H

// Represents a variable
class X {
public:
    int id;
    X(int id) { this->id = id; }

    operator int() const { return id; }
};


#endif //GLPKPP_X_H
