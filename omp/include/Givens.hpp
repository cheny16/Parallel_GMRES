/**
 * @file Givens.hpp
 * @author cheny16
 * @brief Givens header file
 * @version 1.0
 * @date 2023-03-26
 * 
 */
#ifndef GIVENS_HPP
#define GIVENS_HPP

void GenGivensRotator(const double v1, const double v2, double *c, double *s);
void GivensRotate(const double c, const double s, double *x, double *y);

#endif /* GIVENS_HPP */