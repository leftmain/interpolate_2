#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "msr_matrix.h"

// (a, b)
double scalar_product (const double *a, const double *b,
                       unsigned int from, unsigned int to);

// c = a + b * alpha
void linear_combination (const double *a, const double *b, double alpha,
                         double *c, unsigned int from, unsigned int to);

// b = alpha * a
void mult_vector_on_scalar (const double *a, double *b, double alpha,
                            unsigned int from, unsigned int to);

void multiply_msr_matrix_on_vector (const msr_matrix &a, const double *v,
                                    double *r, unsigned int from,
                                    unsigned int to);

int solve_system (const msr_matrix &a, double *r, double *v,
                  double *u, double *x, int maxit, double stop_cond,
                  int p, int t, pthread_barrier_t *barrier, double *buf);

double norm (const double *a, unsigned int n);

double residual (const msr_matrix &a, const double *x,
                 unsigned int from, unsigned int to,
                 int p, int t, double *buf, pthread_barrier_t *barrier);

void print_vector (const double *a, unsigned int n);

double all_reduce_max (double a, int p, int t, double *buf,
                       pthread_barrier_t *barrier);

double all_reduce_sum (double a, int p, int t, double *buf,
                       pthread_barrier_t *barrier);

#endif  // _SOLVER_H_
