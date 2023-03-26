#include "../include/Matrix.hpp"
#include <cmath>
#include <iostream>

void GenGivensRotator(const double v1, const double v2, double *c, double *s) {
    // sqrt(v1^2 + v2^2)
    double sq = std::hypot(v1, v2);
    *c = v1 / sq;
    *s = v2 / sq;
}

void GivensRotate(const double c, const double s, double *x, double *y) {
    double tmp = c * (*x) + s * (*y);
    *y = -s * (*x) + c * (*y);
    *x = tmp;
}