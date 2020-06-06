#include "alg.h"

double get_cpu_time ()
{
  struct rusage time;
  getrusage (RUSAGE_THREAD, &time);
  return (double)time.ru_utime.tv_sec
       + (double)time.ru_utime.tv_usec * 1e-6;
}

double get_global_time ()
{
  struct timeval t;
  gettimeofday (&t, NULL);
  return (double)t.tv_sec + (double)t.tv_usec * 1e-6;
}

double get_length (double x1, double y1, double x2, double y2)
{
  return sqrt ((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

unsigned int get_k_by_ij (unsigned int nx, unsigned int ,
                          unsigned int i, unsigned int j)
{
  return i * (nx + 1) + j;
}

void get_xy_by_ij (const quad &q, double nx, double ,
                   double *x, double *y,
                   unsigned int i, unsigned int j)
{
#if 0
  auto offset = [] (double x0, double y0, double x1, double y1, double d,
                    double *x, double *y, int i)
    {
      double ax = x1 - x0;
      double ay = y1 - y0;
      double al = sqrt (ax * ax + ay * ay);
      *x = x0 + i * d * ax / al;
      *y = y0 + i * d * ay / al;
    };

  *x = 0;
  *y = 0;

  offset (q.x0, q.y0, q.x1, q.y1, q.dy0, x, y, i);

  double xx = 0.;
  double yy = 0.;

  offset (q.x3, q.y3, q.x2, q.y2, q.dy1, &xx, &yy, i);

  double dx = get_length (*x, *y, xx, yy) / nx;

  offset (*x, *y, xx, yy, dx, x, y, j);
#else
  double al = sqrt ((q.x1 - q.x0) * (q.x1 - q.x0) + (q.y1 - q.y0) * (q.y1 - q.y0));
  *x = q.x0 + i * q.dy0 * (q.x1 - q.x0) / al;
  *y = q.y0 + i * q.dy0 * (q.y1 - q.y0) / al;
  al = sqrt ((q.x2 - q.x3) * (q.x2 - q.x3) + (q.y2 - q.y3) * (q.y2 - q.y3));
  double xx = q.x3 + i * q.dy1 * (q.x2 - q.x3) / al;
  double yy = q.y3 + i * q.dy1 * (q.y2 - q.y3) / al;
  al = sqrt ((*x - xx) * (*x - xx) + (*y - yy) * (*y - yy));
  *x += j * (xx - *x) / nx;
  *y += j * (yy - *y) / nx;
#endif
}

// fi = f (xi, yi)
// gi = g (xi, yi)
// returns integral [ fg, dxdy ] (on triangle)
double
integrate (double x1, double x2, double x3,
           double y1, double y2, double y3,
           double f1, double f2, double f3,
           double g1, double g2, double g3)
{
  double J = fabs ((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3));

  double a = f1 - f3;
  double b = f2 - f3;
  double c = f3;

  double d = g1 - g3;
  double e = g2 - g3;
  double f = g3;

  double res = a * (2 * d + e + 4 * f)
             + b * (d + 2 * e + 4 * f)
             + 4 * c * (d + e + 3 * f);

  return res * J / 24.0;
}

// fi = f (xi, yi)
// gi = g (xi, yi)
// returns integral [ fg, dxdy ] (on triangle)
double
integrate (int p1, int p2, int p3,
           double *x, double *y,
           double f1, double f2, double f3,
           double g1, double g2, double g3)
{
  return integrate (x[p1], x[p2], x[p3],
                    y[p1], y[p2], y[p3],
                    f1, f2, f3,
                    g1, g2, g3);
}

void *start_count_thread (void *data_)
{
  method_data *data = (method_data *)data_;

  int p = data->p;

  // count thread stuff
  std::unique_ptr<count_thread_args []> args;
  std::unique_ptr<pthread_t []> threads;
  args = std::make_unique<count_thread_args []> (p);
  threads = std::make_unique<pthread_t []> (p - 1);
  if (threads.get () == nullptr || args.get () == nullptr)
    {
      printf ("mem error\n");
      return nullptr;
    }
  
  pthread_barrier_t local_barrier;
  pthread_barrier_init (&local_barrier, 0, p);

  for (int i = 0; i < p; ++i)
    {
      args[i].local_barrier = &local_barrier;
      args[i].data = data;
      args[i].t = i;
    }

  for (int i = 1; i < p; ++i)
    {
      int ret = pthread_create (threads.get () + i - 1,
                                NULL, build_approximation, args.get () + i);
      if (ret)
        {
          printf ("cannot create thread\n");
          std::abort ();
        }
    }
  build_approximation (args.get ());
  
  for (int i = 1; i < p; ++i)
    {
      pthread_join (threads[i - 1], NULL);
    }

  pthread_barrier_destroy (&local_barrier);

  //printf ("EEEEEEEEEEEEEEND OF COUNT THREADS\n");

  return nullptr;
}

void *build_approximation (void *args_)
{
  count_thread_args *args = (count_thread_args *)args_;
  auto data = args->data;
  auto t = args->t;
  auto local_barrier = args->local_barrier;
  auto eps = data->eps;
  auto p = data->p;
  int maxit = 100;

  while (true)
    {
#ifdef DRAW_GRAPH
      pthread_barrier_wait (local_barrier);
      if (t == 0)
        data->new_count.lock ();
      pthread_barrier_wait (local_barrier);
      if (data->end_program)
        break;
#endif

      auto &f = data->f.get ();

      auto nx = data->nx;
      auto ny = data->ny;
      unsigned int n_rows = (nx + 1) * (ny + 1);
      unsigned int own_rows = (n_rows + p - 1) / p;
      unsigned int row_begin = std::min (t * own_rows, n_rows);
      unsigned int row_end = std::min ((t + 1) * own_rows, n_rows);

      auto count_quad = [&] (msr_matrix &matrix, quad &q,
                             std::unique_ptr<double []> &workspace,
                             double **result, unsigned int *size,
                             quad &q1, int id)
        {
          if (t == 0)
            {
              // for matrix
              error ret = matrix.init (n_rows, 7 * n_rows);
              if (ret != error::all_ok)
                {
                  printf ("init_msr error\n");
                  std::abort ();
                }

#ifdef DEBIG_PRINT
              printf ("n_rows = %u\n", n_rows);
#endif

              // for solver
              unsigned int new_size = 4 * n_rows + 2 * p + 2;
              if (*size < new_size)
                {
                  workspace = std::make_unique<double []> (new_size);
                  if (workspace.get () == nullptr)
                    {
                      printf ("mem error\n");
                      std::abort ();
                      return ;
                    }
                  *size = new_size;
                }

              build_rnz_matrix (matrix, nx, ny);
            }

          pthread_barrier_wait (local_barrier);
          fill_matrix (matrix, nx, ny, q, p, t, f, q1, id);
          pthread_barrier_wait (local_barrier);
          if (t == 0)
            {
              matrix.print (stdout);
            }

          auto b = matrix.get_rhs ();
          double b_norm = norm (b, n_rows);
          double stop_cond = eps * eps * b_norm * b_norm;
          double *r = workspace.get ();
          double *v = r + n_rows;
          double *u = v + n_rows;
          double *x = u + n_rows;
          double *buf = x + n_rows;

          // time start
          double time = 0.0;
          pthread_barrier_wait (local_barrier);
          if (t == 0) time = get_global_time ();
          double solver_time = get_cpu_time ();

          // solver
          auto it = solve_system (matrix, r, v, u, x, maxit, stop_cond,
                                  p, t, local_barrier, buf);

          // time stop
          solver_time = get_cpu_time () - solver_time;
          pthread_barrier_wait (local_barrier);
          if (t == 0) time = get_global_time () - time;

          double resid = residual (matrix, x, row_begin, row_end,
                                   p, t, buf, local_barrier);

          double approx_resid = residual (x, nx, ny, q, f, p, t);
          all_reduce_max (approx_resid, p, t, buf, local_barrier);

          buf[t] = solver_time;
          pthread_barrier_wait (local_barrier);

          if (t == 0)
            {
              printf ("solver time:\t(global: %.2lf)\n", time);
              for (int i = 0; i < p; ++i)
                {
                  printf ("\tthread %d time: %.2lf\n", i, buf[i]);
                }
              printf ("it = %d\tresidual = %le\n", it, resid);
              printf ("approx residual = %le\n", approx_resid);
              data->residual = std::max (approx_resid, data->residual);
              *result = x;
              //printf ("%p\n", x);
            }
        };

      data->residual = 0.0;
      for (unsigned int i = 0; i < 2; ++i)
        {
          count_quad (*data->matrix[i], data->q[i], *data->workspace[i],
                      &data->c[i], &data->workspace_size[i],
                      data->q[(i + 1) % 2], i);
        }
      pthread_barrier_wait (local_barrier);
#if 1
      // count average on the edge
      auto c = data->c[0];
      auto c1 = data->c[1];
      for (unsigned int i = 0; i <= ny; ++i)
        {
          auto ind = get_k_by_ij (nx, ny, i, nx);
          auto ind1 = get_k_by_ij (nx, ny, i, 0);
          auto z = (c[ind] + c1[ind1]) / 2.0;
          c[ind] = c1[ind1] = z;
        }
#endif

      pthread_barrier_wait (local_barrier);
#ifdef DRAW_GRAPH
      if (t == 0)
        data->end_count.unlock ();
      pthread_barrier_wait (local_barrier);
#else
      break;
#endif
    }

  return nullptr;
}

unsigned int get_row_len_and_fill_rnz (unsigned int nx, unsigned int ny,
                                       unsigned int i, unsigned int j,
                                       unsigned int *rnz)
{
  unsigned int row = i * (nx + 1) + j;
  unsigned int row_size = 0;
  auto rnz_offset = rnz[row];

  if (i > 0 && i < ny && j > 0 && j < nx)
    {
      row_size = 7;
      rnz[rnz_offset + 0] = get_k_by_ij (nx, ny, i - 1, j);
      rnz[rnz_offset + 1] = get_k_by_ij (nx, ny, i - 1, j + 1);
      rnz[rnz_offset + 2] = get_k_by_ij (nx, ny, i, j - 1);
      rnz[rnz_offset + 3] = get_k_by_ij (nx, ny, i, j + 1);
      rnz[rnz_offset + 4] = get_k_by_ij (nx, ny, i + 1, j - 1);
      rnz[rnz_offset + 5] = get_k_by_ij (nx, ny, i + 1, j);
    }
  else if (i == 0 && j > 0 && j < nx)
    {
      row_size = 5;
      rnz[rnz_offset + 0] = get_k_by_ij (nx, ny, i, j - 1);
      rnz[rnz_offset + 1] = get_k_by_ij (nx, ny, i, j + 1);
      rnz[rnz_offset + 2] = get_k_by_ij (nx, ny, i + 1, j - 1);
      rnz[rnz_offset + 3] = get_k_by_ij (nx, ny, i + 1, j);
    }
  else if (i == ny && j > 0 && j < nx)
    {
      row_size = 5;
      rnz[rnz_offset + 0] = get_k_by_ij (nx, ny, i - 1, j);
      rnz[rnz_offset + 1] = get_k_by_ij (nx, ny, i - 1, j + 1);
      rnz[rnz_offset + 2] = get_k_by_ij (nx, ny, i, j - 1);
      rnz[rnz_offset + 3] = get_k_by_ij (nx, ny, i, j + 1);
    }
  else if (j == 0 && i > 0 && i < ny)
    {
      row_size = 5;
      rnz[rnz_offset + 0] = get_k_by_ij (nx, ny, i - 1, j);
      rnz[rnz_offset + 1] = get_k_by_ij (nx, ny, i - 1, j + 1);
      rnz[rnz_offset + 2] = get_k_by_ij (nx, ny, i, j + 1);
      rnz[rnz_offset + 3] = get_k_by_ij (nx, ny, i + 1, j);
    }
  else if (j == nx && i > 0 && i < ny)
    {
      row_size = 5;
      rnz[rnz_offset + 0] = get_k_by_ij (nx, ny, i - 1, j);
      rnz[rnz_offset + 1] = get_k_by_ij (nx, ny, i, j - 1);
      rnz[rnz_offset + 2] = get_k_by_ij (nx, ny, i + 1, j - 1);
      rnz[rnz_offset + 3] = get_k_by_ij (nx, ny, i + 1, j);
    }
  else if (i == ny && j == 0)
    {
      row_size = 4;
      rnz[rnz_offset + 0] = get_k_by_ij (nx, ny, i - 1, j);
      rnz[rnz_offset + 1] = get_k_by_ij (nx, ny, i - 1, j + 1);
      rnz[rnz_offset + 2] = get_k_by_ij (nx, ny, i, j + 1);
    }
  else if (i == 0 && j == nx)
    {
      row_size = 4;
      rnz[rnz_offset + 0] = get_k_by_ij (nx, ny, i, j - 1);
      rnz[rnz_offset + 1] = get_k_by_ij (nx, ny, i + 1, j - 1);
      rnz[rnz_offset + 2] = get_k_by_ij (nx, ny, i + 1, j);
    }
  else if (i == 0 && j == 0)
    {
      row_size = 3;
      rnz[rnz_offset + 0] = get_k_by_ij (nx, ny, i, j + 1);
      rnz[rnz_offset + 1] = get_k_by_ij (nx, ny, i + 1, j);
    }
  else if (i == ny && j == nx)
    {
      row_size = 3;
      rnz[rnz_offset + 0] = get_k_by_ij (nx, ny, i - 1, j);
      rnz[rnz_offset + 1] = get_k_by_ij (nx, ny, i, j - 1);
    }
  else
    {
      assert (false);
    }
  assert (row_size != 0);

  auto offdiag_row_size = row_size - 1;

  rnz[row + 1] = rnz[row] + offdiag_row_size;

  return offdiag_row_size;
}

void build_rnz_matrix (msr_matrix &a, unsigned int nx, unsigned int ny)
{
  unsigned int row = 0;
  for (unsigned int i_row = 0; i_row <= ny; ++i_row)
    {
      for (unsigned int j_col = 0; j_col <= nx; ++j_col)
        {
#ifdef DEBUG_PRINT
          //auto row_len = 
#endif
          get_row_len_and_fill_rnz (nx, ny, i_row, j_col,
                                    a.get_internal_rnz ());
#ifdef DEBUG_PRINT
          //printf ("%u\t", row_len + 1);
#endif
          ++row;
        }
#ifdef DEBUG_PRINT
      //printf ("\n");
#endif
    }
}

///////////////////////////////////////////////////////////////

struct triangle
{
  // p0 = 0
  int p1 = 0;
  int p2 = 0;
};

void fill_matrix (msr_matrix &a, unsigned int nx, unsigned int ny, quad &q,
                  int p, int t, const Func &f, quad &, int )
//                  int p, int t, const Func &f, quad &q1, int id)
{
  auto rnz = a.get_rnz ();
  auto data = a.get_internal_data ();
  auto rhs = a.get_internal_rhs ();

  double x[7];
  double y[7];
  double ff[7];
  double integr[7];

  unsigned int own_i = ((ny + 1) + p - 1) / p;
  unsigned int i_from = std::min (t * own_i, ny + 1);
  unsigned int i_to = std::min (i_from + own_i, ny + 1);

  auto integrate_by_triang2 = [] (double *x, double *y, double *integr,
                                  int d1, int d2, int d3)
    {
      integr[d1] = integrate (0, d1, d2, x, y, 1, 0, 0, 0, 1, 0)
                 + integrate (0, d1, d3, x, y, 1, 0, 0, 0, 1, 0);
    };

  auto integrate_by_triang1 = [] (double *x, double *y, double *integr,
                                  int d1, int d2)
    {
      integr[d1] = integrate (0, d1, d2, x, y, 1, 0, 0, 0, 1, 0);
    };

  auto on_triangle = [&] (int p1, int p2, int p3,
                          double *x, double *y, double *ff)
    {
      double x1 = x[p1], y1 = y[p1], f1 = ff[p1];
      double x2 = x[p2], y2 = y[p2], f2 = ff[p2];
      double x3 = x[p3], y3 = y[p3], f3 = ff[p3];

      double x22 = (x1 + x2) / 2.;
      double y22 = (y1 + y2) / 2.;
      double f22 = f (x22, y22);
      double x33 = (x1 + x3) / 2.;
      double y33 = (y1 + y3) / 2.;
      double f33 = f (x33, y33);
      double x23 = (x2 + x3) / 2.;
      double y23 = (y2 + y3) / 2.;
      double f23 = f (x23, y23);

      return integrate (x1, x22, x33,
                        y1, y22, y33,
                        f1, f22, f33, 1.0, 0.5, 0.5)
           + integrate (x22, x2, x23,
                        y22, y2, y23,
                        f22, f2, f23, 0.5, 0.0, 0.0)
           + integrate (x33, x3, x23,
                        y33, y3, y23,
                        f33, f3, f23, 0.5, 0.0, 0.0)
           + integrate (x22, x33, x23,
                        y22, y33, y23,
                        f22, f33, f23, 0.5, 0.5, 0.0);
    };

  auto integral_by_center = [&] (int d1, int d2, double *x, double *y,
                                double *ff, double *rhs_el)
    {
      *rhs_el += on_triangle (0, d1, d2, x, y, ff);
      return integrate (0, d1, d2, x, y, 1, 0, 0, 1, 0, 0);
    };

  for (unsigned int i = i_from; i < i_to; ++i)
    {
      for (unsigned int j = 0; j <= nx; ++j)
        {
          int dots = 0;
          auto i_row = get_k_by_ij (nx, ny, i, j);

          if (i > 0 && i < ny && j > 0 && j < nx)
            {
              dots = 7;

              get_xy_by_ij (q, nx, ny, x + 0, y + 0, i, j);
              get_xy_by_ij (q, nx, ny, x + 1, y + 1, i - 1, j);
              get_xy_by_ij (q, nx, ny, x + 2, y + 2, i - 1, j + 1);
              get_xy_by_ij (q, nx, ny, x + 3, y + 3, i, j - 1);
              get_xy_by_ij (q, nx, ny, x + 4, y + 4, i, j + 1);
              get_xy_by_ij (q, nx, ny, x + 5, y + 5, i + 1, j - 1);
              get_xy_by_ij (q, nx, ny, x + 6, y + 6, i + 1, j);

              for (int k = 0; k < dots; ++k)
                ff[k] = f (x[k], y[k]);

              double rhs_el = 0.0;
              integr[0] = integral_by_center (1, 2, x, y, ff, &rhs_el)
                        + integral_by_center (2, 4, x, y, ff, &rhs_el)
                        + integral_by_center (4, 6, x, y, ff, &rhs_el)
                        + integral_by_center (6, 5, x, y, ff, &rhs_el)
                        + integral_by_center (5, 3, x, y, ff, &rhs_el)
                        + integral_by_center (3, 1, x, y, ff, &rhs_el);

              rhs[i_row] = rhs_el;

              integrate_by_triang2 (x, y, integr, 1, 3, 2);
              integrate_by_triang2 (x, y, integr, 2, 1, 4);
              integrate_by_triang2 (x, y, integr, 3, 5, 1);
              integrate_by_triang2 (x, y, integr, 4, 2, 6);
              integrate_by_triang2 (x, y, integr, 5, 6, 3);
              integrate_by_triang2 (x, y, integr, 6, 4, 5);
            }
          else if (i == 0 && j > 0 && j < nx)
            {
              dots = 5;
              
              get_xy_by_ij (q, nx, ny, x + 0, y + 0, i, j);
              get_xy_by_ij (q, nx, ny, x + 1, y + 1, i, j - 1);
              get_xy_by_ij (q, nx, ny, x + 2, y + 2, i, j + 1);
              get_xy_by_ij (q, nx, ny, x + 3, y + 3, i + 1, j - 1);
              get_xy_by_ij (q, nx, ny, x + 4, y + 4, i + 1, j);

              for (int k = 0; k < dots; ++k)
                ff[k] = f (x[k], y[k]);

              double rhs_el = 0.0;
              integr[0] = integral_by_center (1, 3, x, y, ff, &rhs_el)
                        + integral_by_center (3, 4, x, y, ff, &rhs_el)
                        + integral_by_center (4, 2, x, y, ff, &rhs_el);

              rhs[i_row] = rhs_el;

              integrate_by_triang1 (x, y, integr, 1, 3);
              integrate_by_triang2 (x, y, integr, 3, 1, 4);
              integrate_by_triang2 (x, y, integr, 4, 3, 2);
              integrate_by_triang1 (x, y, integr, 2, 4);
            }
          else if (i == ny && j > 0 && j < nx)
            {
              dots = 5;
              
              get_xy_by_ij (q, nx, ny, x + 0, y + 0, i, j);
              get_xy_by_ij (q, nx, ny, x + 1, y + 1, i - 1, j);
              get_xy_by_ij (q, nx, ny, x + 2, y + 2, i - 1, j + 1);
              get_xy_by_ij (q, nx, ny, x + 3, y + 3, i, j - 1);
              get_xy_by_ij (q, nx, ny, x + 4, y + 4, i, j + 1);

              for (int k = 0; k < dots; ++k)
                ff[k] = f (x[k], y[k]);

              double rhs_el = 0.0;
              integr[0] = integral_by_center (3, 1, x, y, ff, &rhs_el)
                        + integral_by_center (1, 2, x, y, ff, &rhs_el)
                        + integral_by_center (2, 4, x, y, ff, &rhs_el);

              rhs[i_row] = rhs_el;

              integrate_by_triang1 (x, y, integr, 3, 1);
              integrate_by_triang2 (x, y, integr, 1, 3, 2);
              integrate_by_triang2 (x, y, integr, 2, 1, 4);
              integrate_by_triang1 (x, y, integr, 4, 2);
            }
          else if (j == 0 && i > 0 && i < ny)
            {
              dots = 5;
              
              get_xy_by_ij (q, nx, ny, x + 0, y + 0, i, j);
              get_xy_by_ij (q, nx, ny, x + 1, y + 1, i - 1, j);
              get_xy_by_ij (q, nx, ny, x + 2, y + 2, i - 1, j + 1);
              get_xy_by_ij (q, nx, ny, x + 3, y + 3, i, j + 1);
              get_xy_by_ij (q, nx, ny, x + 4, y + 4, i + 1, j);

              for (int k = 0; k < dots; ++k)
                ff[k] = f (x[k], y[k]);

              double rhs_el = 0.0;
              integr[0] = integral_by_center (1, 2, x, y, ff, &rhs_el)
                        + integral_by_center (2, 3, x, y, ff, &rhs_el)
                        + integral_by_center (3, 4, x, y, ff, &rhs_el);

              rhs[i_row] = rhs_el;

              integrate_by_triang1 (x, y, integr, 1, 2);
              integrate_by_triang2 (x, y, integr, 2, 1, 3);
              integrate_by_triang2 (x, y, integr, 3, 2, 4);
              integrate_by_triang1 (x, y, integr, 4, 3);
            }
          else if (j == nx && i > 0 && i < ny)
            {
              dots = 5;
              
              get_xy_by_ij (q, nx, ny, x + 0, y + 0, i, j);
              get_xy_by_ij (q, nx, ny, x + 1, y + 1, i - 1, j);
              get_xy_by_ij (q, nx, ny, x + 2, y + 2, i, j - 1);
              get_xy_by_ij (q, nx, ny, x + 3, y + 3, i + 1, j - 1);
              get_xy_by_ij (q, nx, ny, x + 4, y + 4, i + 1, j);

              for (int k = 0; k < dots; ++k)
                ff[k] = f (x[k], y[k]);

              double rhs_el = 0.0;
              integr[0] = integral_by_center (1, 2, x, y, ff, &rhs_el)
                        + integral_by_center (2, 3, x, y, ff, &rhs_el)
                        + integral_by_center (3, 4, x, y, ff, &rhs_el);

              rhs[i_row] = rhs_el;

              integrate_by_triang1 (x, y, integr, 1, 2);
              integrate_by_triang2 (x, y, integr, 2, 1, 3);
              integrate_by_triang2 (x, y, integr, 3, 2, 4);
              integrate_by_triang1 (x, y, integr, 4, 3);
            }
          else if (i == ny && j == 0)
            {
              dots = 4;
              
              get_xy_by_ij (q, nx, ny, x + 0, y + 0, i, j);
              get_xy_by_ij (q, nx, ny, x + 1, y + 1, i - 1, j);
              get_xy_by_ij (q, nx, ny, x + 2, y + 2, i - 1, j + 1);
              get_xy_by_ij (q, nx, ny, x + 3, y + 3, i, j + 1);

              for (int k = 0; k < dots; ++k)
                ff[k] = f (x[k], y[k]);

              double rhs_el = 0.0;
              integr[0] = integral_by_center (1, 2, x, y, ff, &rhs_el)
                        + integral_by_center (2, 3, x, y, ff, &rhs_el);

              rhs[i_row] = rhs_el;

              integrate_by_triang1 (x, y, integr, 1, 2);
              integrate_by_triang2 (x, y, integr, 2, 1, 3);
              integrate_by_triang1 (x, y, integr, 3, 2);
            }
          else if (i == 0 && j == nx)
            {
              dots = 4;
              
              get_xy_by_ij (q, nx, ny, x + 0, y + 0, i, j);
              get_xy_by_ij (q, nx, ny, x + 1, y + 1, i, j - 1);
              get_xy_by_ij (q, nx, ny, x + 2, y + 2, i + 1, j - 1);
              get_xy_by_ij (q, nx, ny, x + 3, y + 3, i + 1, j);

              for (int k = 0; k < dots; ++k)
                ff[k] = f (x[k], y[k]);

              double rhs_el = 0.0;
              integr[0] = integral_by_center (1, 2, x, y, ff, &rhs_el)
                        + integral_by_center (2, 3, x, y, ff, &rhs_el);

              rhs[i_row] = rhs_el;

              integrate_by_triang1 (x, y, integr, 1, 2);
              integrate_by_triang2 (x, y, integr, 2, 1, 3);
              integrate_by_triang1 (x, y, integr, 3, 2);
            }
          else if (i == 0 && j == 0)
            {
              dots = 3;
              
              get_xy_by_ij (q, nx, ny, x + 0, y + 0, i, j);
              get_xy_by_ij (q, nx, ny, x + 1, y + 1, i, j + 1);
              get_xy_by_ij (q, nx, ny, x + 2, y + 2, i + 1, j);

              for (int k = 0; k < dots; ++k)
                ff[k] = f (x[k], y[k]);

              double rhs_el = 0.0;
              integr[0] = integral_by_center (1, 2, x, y, ff, &rhs_el);
              rhs[i_row] = rhs_el;

              integrate_by_triang1 (x, y, integr, 1, 2);
              integrate_by_triang1 (x, y, integr, 2, 1);
            }
          else if (i == ny && j == nx)
            {
              dots = 3;
              
              get_xy_by_ij (q, nx, ny, x + 0, y + 0, i, j);
              get_xy_by_ij (q, nx, ny, x + 1, y + 1, i - 1, j);
              get_xy_by_ij (q, nx, ny, x + 2, y + 2, i, j - 1);

              for (int k = 0; k < dots; ++k)
                ff[k] = f (x[k], y[k]);

              double rhs_el = 0.0;
              integr[0] = integral_by_center (1, 2, x, y, ff, &rhs_el);
              rhs[i_row] = rhs_el;

              integrate_by_triang1 (x, y, integr, 1, 2);
              integrate_by_triang1 (x, y, integr, 2, 1);
            }
          else
            {
              assert (false);
            }

          data[i_row] = integr[0];
          auto offset = rnz[i_row];
          for (int j_col = 0; j_col < dots - 1; j_col++)
            {
              data[offset + j_col] = integr[j_col + 1];
            }
        }
    }
}

double residual (const double *c, unsigned int nx, unsigned int ny,
                 const quad &q, const Func &f, int p, int t)
{
  double residual = 0.0;
  double x = 0.0;
  double y = 0.0;

  unsigned int own_i = ((ny + 1) + p - 1) / p;
  unsigned int i_from = std::min (t * own_i, ny + 1);
  unsigned int i_to = std::min (i_from + own_i, ny + 1);

  for (unsigned int i = i_from; i < i_to; ++i)
    {
      for (unsigned int j = 0; j <= nx; ++j)
        {
          auto i_row = get_k_by_ij (nx, ny, i, j);
          auto pf = c[i_row];

          get_xy_by_ij (q, nx, ny, &x, &y, i, j);
          residual = std::max (residual, fabs (f (x, y) - pf));
        }
    }

  return residual;
}

