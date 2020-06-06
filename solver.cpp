#include "solver.h"

 
double scalar_product (const double *a, const double *b,
                       unsigned int from, unsigned int to)
{
  double sum = 0.0;

  for (unsigned int i = from; i < to; ++i)
    {
      //printf ("i = %u: %le\n", i, a[i] * b[i]);
      sum += a[i] * b[i];
    }

  return sum;
}

// c = a + b * alpha
void linear_combination (const double *a, const double *b, double alpha,
                         double *c, unsigned int from, unsigned int to)
{
  for (unsigned int i = from; i < to; ++i)
    {
      c[i] = a[i] + b[i] * alpha;
    }
}

// b = alpha * a
void mult_vector_on_scalar (const double *a, double *b, double alpha,
                            unsigned int from, unsigned int to)
{
  for (unsigned int i = from; i < to; ++i)
    {
      b[i] = a[i] * alpha;
    }
}

void multiply_msr_matrix_on_vector (const msr_matrix &a, const double *v,
                                    double *r, unsigned int from,
                                    unsigned int to)
{
  memset (r + from, 0, (to - from) * sizeof (double));

  for (unsigned int i_row = from; i_row < to; ++i_row)
    {
      double d = 0.0;

      // diag element
      d += a.get_diag (i_row) * v[i_row];

      // offdiag element
      auto row_len = a.get_offdiag_len (i_row);
      for (unsigned int j_col = 0; j_col < row_len; ++j_col)
        {
          d += a.get_offdiag (i_row, j_col)
               * v[a.get_col_num (i_row, j_col)];
        }

      r[i_row] = d;
    }
}

// v = M^-1 * r
void
apply_jacobi (const msr_matrix &a, double *r, double *v,
              unsigned int from, unsigned int to)
{
  auto diag = a.get_data ();

  for (unsigned int i = from; i < to; ++i)
    {
      v[i] = r[i] / diag[i];
    }
}

int solve_system (const msr_matrix &a, double *r, double *v,
                  double *u, double *x, int maxit, double stop_cond,
                  int p, int t, pthread_barrier_t *barrier, double *buf)
{
  auto b = a.get_rhs ();
  auto n = a.get_n_rows ();
  int it = 0;

  unsigned int own_rows = (n + p - 1) / p;
  unsigned int row_begin = std::min (t * own_rows, n);
  unsigned int row_end = std::min ((t + 1) * own_rows, n);

  // x = 0
  memset (x + row_begin, 0, (row_end - row_begin) * sizeof (double));

  // r = Ax - b = -b
  mult_vector_on_scalar (b, r, -1.0, row_begin, row_end);
 
  for (it = 1; it < maxit; ++it)
    {
      // v = M^-1 * r
      apply_jacobi (a, r, v, row_begin, row_end);

      // wait to mult A on v
      pthread_barrier_wait (barrier);

      // u = Av
      multiply_msr_matrix_on_vector (a, v, u, row_begin, row_end);

      // c1 = (u, r)
      double c1 = scalar_product (u, r, row_begin, row_end);
      c1 = all_reduce_sum (c1, p, t, buf, barrier);

      // c2 = (u, u)
      double c2 = scalar_product (u, u, row_begin, row_end);
      c2 = all_reduce_sum (c2, p, t, buf, barrier);

      //
      if (fabs (c1) < stop_cond || c2 < stop_cond)
        break;

      // tau = c1 / c2
      double tau = c1 / c2;

      // x = x - v * tau
      linear_combination (x, v, -tau, x, row_begin, row_end);

      // r = r - u * tau
      linear_combination (r, u, -tau, r, row_begin, row_end);
    }

  pthread_barrier_wait (barrier);
  return it;
}

double norm (const double *a, unsigned int n)
{
  double max = 0.0;
  for (unsigned int i = 0; i < n; ++i)
    {
      if (fabs (a[i]) > max) max = fabs (a[i]);
    }
  return max;
}

double residual (const msr_matrix &a, const double *x,
                 unsigned int from, unsigned int to,
                 int p, int t, double *buf, pthread_barrier_t *barrier)
{
  auto b = a.get_rhs ();
  double max = 0.0;

  for (unsigned int i_row = from; i_row < to; ++i_row)
    {
      double sum = 0.0;

      // diag element
      sum += a.get_diag (i_row) * x[i_row];

      // offdiag element
      auto row_len = a.get_offdiag_len (i_row);
      for (unsigned int j_col = 0; j_col < row_len; ++j_col)
        {
          sum += a.get_offdiag (i_row, j_col)
                 * x[a.get_col_num (i_row, j_col)];
        }

      sum -= b[i_row];

      //printf ("i = %u:\t%le\n", i_row, sum);
      if (fabs (sum) > max)
        max = fabs (sum);
    }

  return all_reduce_max (max, p, t, buf, barrier);
}

void print_vector (const double *a, unsigned int n)
{
  for (unsigned int i = 0; i < n; ++i)
    {
      printf ("%.4le\t", a[i]);
    }
  printf ("\n");
}

double all_reduce_max (double a, int p, int t, double *buf,
                       pthread_barrier_t *barrier)
{
  buf[t] = a;
  pthread_barrier_wait (barrier);
  double max = 0.0;
  if (t == 0)
    {
      max = buf[0];
      for (int i = 1; i < p; ++i)
        {
          if (buf[i] > max) max = buf[i];
        }
    }
  pthread_barrier_wait (barrier);
  return max;

}

double all_reduce_sum (double a, int p, int t, double *buf,
                       pthread_barrier_t *barrier)
{
  buf[t] = a;
  pthread_barrier_wait (barrier);
  if (t == 0)
    {
      buf[p] = 0.0;
      for (int i = 0; i < p; ++i)
        {
          buf[p] += buf[i];
        }
    }
  pthread_barrier_wait (barrier);
  return buf[p];
}

