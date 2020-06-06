#ifndef _ALG_H_
#define _ALG_H_

#include "header.h"
#include "msr_matrix.h"
#include "solver.h"

struct quad
{
  double x0 = 0.0;
  double y0 = 0.0;
  double x1 = 0.0;
  double y1 = 0.0;
  double x2 = 0.0;
  double y2 = 0.0;
  double x3 = 0.0;
  double y3 = 0.0;
  double dy0 = 0.0;
  double dy1 = 0.0;
};

struct method_data
{
  // canges
  msr_matrix *matrix[2] = {nullptr, nullptr};
  double *c[2] = {nullptr, nullptr}; // result
  std::unique_ptr<double []> *workspace[2] = {nullptr, nullptr};
  unsigned int workspace_size[2] = {0, 0};

  // not changes
  quad q[2];
  unsigned int nx = 0;
  unsigned int ny = 0;
  Functions f;
  int k = 0;
  int p = 1;

  // solver
  double eps = 1e-14;
  double time = 0.0;
  double residual = 0.0;

  // thread
  std::mutex end_count;
  std::mutex new_count;

  bool end_program = false;
};

struct count_thread_args
{
  method_data *data;
  pthread_barrier_t *local_barrier = nullptr;
  int t = 0; // thread num
};

double get_length (double x1, double y1, double x2, double y2);

void *start_count_thread (void *data);

void *build_approximation (void *);

void build_rnz_matrix (msr_matrix &a, unsigned int nx, unsigned int ny);

void fill_matrix (msr_matrix &a, unsigned int nx, unsigned int ny, quad &q,
                  int p, int t, const Func &f, quad &q1, int id);

unsigned int get_k_by_ij (unsigned int nx, unsigned int ,
                          unsigned int i, unsigned int j);

void get_xy_by_ij (const quad &q, double nx, double ny,
                   double *x, double *y,
                   unsigned int i, unsigned int j);

double residual (const double *c, unsigned int nx, unsigned int ny,
                 const quad &q, const Func &f, int p, int t);

#endif // _ALG_H_
