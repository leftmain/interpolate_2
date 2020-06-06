#include "header.h"
#include "msr_matrix.h"
#include "alg.h"

#ifdef DRAW_GRAPH
#include <qapplication.h>
#include <qgl.h>
#include "glwidget.h"
#endif

error parse_args (int argc, char **argv, method_data &data);
double get_length (double x1, double y1, double x2, double y2);

error parse_file (const char *file_name,
                  double *x1, double *y1,
                  double *x2, double *y2,
                  double *x3, double *y3, double *l);

int main (int argc, char **argv)
{
  // set data
  method_data data;
  error ret = parse_args (argc, argv, data);
  if (ret != error::all_ok)
    return -1;
  const Function funcs[] = {
    {"f(x,y) = 1", [] (double, double){ return 1.; }},
    {"f(x,y) = x", [] (double x, double ) { return x; }},
    {"f(x,y) = y", [] (double , double y) { return y; }},
    {"f(x,y) = x + y", [] (double x, double y) { return x + y; }},
    {"f(x,y) = sqrt(x + y)", [] (double x, double y) { return sqrt (x + y); }},
    {"f(x,y) = x^2 + y^2", [] (double x, double y) { return x * x + y * y; }},
    {"f(x,y) = exp(x^2 - y^2)", [] (double x, double y)
                                  { return exp (x * x - y * y); }},
    {"f(x,y) = 1/(25(x^2 + y^2) + 1)", [] (double x, double y)
                                  { return 1. / (25 * (x * x + y * y) + 1); }}
  };
  int number_of_funcs = sizeof (funcs) / sizeof (Function);
  if (data.k > number_of_funcs - 1)
    {
      printf ("k must be less than %d\n", number_of_funcs);
      return -1;
    }
  data.f.init (funcs, data.k, number_of_funcs);
  msr_matrix matrix_0;
  msr_matrix matrix_1;
  std::unique_ptr<double []> workspace_0;
  std::unique_ptr<double []> workspace_1;
  data.matrix[0] = &matrix_0;
  data.matrix[1] = &matrix_1;
  data.workspace[0] = &workspace_0;
  data.workspace[1] = &workspace_1;

#ifdef DRAW_GRAPH
  data.end_count.lock ();
  data.new_count.lock ();
#endif

  // start count threads
  pthread_t count_thread;
  int res = pthread_create (&count_thread, NULL, start_count_thread, &data);
  if (res)
    {
      printf ("cannot create thread\n");
    }

#ifdef DRAW_GRAPH

  // window stuff
  QApplication app (argc, argv);
  GLWidget glwidget (&data);

  app.setActiveWindow (&glwidget);
  glwidget.show ();

  //pthread_join (count_thread, NULL);
  return app.exec ();

#else

  pthread_join (count_thread, NULL);
  return 0;

#endif
}

error parse_args (int argc, char **argv, method_data &data)
{
  int max_k = 9;
  int max_p = 100;

  if (argc != 7
      || sscanf (argv[2], "%u", &data.nx) != 1
      || sscanf (argv[3], "%u", &data.ny) != 1
      || sscanf (argv[4], "%d", &data.k) != 1
      || sscanf (argv[5], "%lf", &data.eps) != 1
      || sscanf (argv[6], "%d", &data.p) != 1
      || data.k > max_k || data.p > max_p)
    {
      printf ("usage: %s [file] nx ny k eps p\n", argv[0]);
      return error::bad_data;
    }
  if (data.k > max_k || data.k < 0)
    {
      printf ("k must be from 0 to %u\n", max_k);
      return error::bad_data;
    }
  if (data.p > max_p || data.p < 0)
    {
      printf ("p must be from 0 to %u\n", max_p);
      return error::bad_data;
    }
  const char *file_name = argv[1];

  double x1, y1, x2, y2, x3, y3, l;
  error ret = parse_file (file_name, &x1, &y1, &x2, &y2, &x3, &y3, &l);
  if (ret != error::all_ok)
    {
      printf ("problem with file\n");
      return ret;
    }

  double ll = get_length (0., 0., x1, y1);
  if (l > ll || l < 1e-16)
    {
      printf ("l is too big or too small\n");
      return error::bad_data;
    }
  data.q[0].x1 = x1;
  data.q[0].y1 = y1;
  data.q[0].x2 = x2;
  data.q[0].y2 = y2;
  data.q[0].x0 = x1 * l / ll;
  data.q[0].y0 = y1 * l / ll;
  auto root_of_lines = [] (double x1, double y1, double x2, double y2,
                           double x3, double y3, double x4, double y4,
                           double *x, double *y)
    {
      auto get_line_coeff = [] (double x1, double y1,
                                double x2, double y2,
                                double *a, double *b, double *c)
        {
          *a = y2 - y1;
          *b = x1 - x2;
          *c = y1 * (x2 - x1) - x1 * (y2 - y1);
        };

      double a, b, c;
      get_line_coeff (x1, y1, x2, y2, &a, &b, &c);

      double d, e, f;
      get_line_coeff (x3, y3, x4, y4, &d, &e, &f);

      if (fabs (b * d - e * a) < 1e-16)
        {
          printf ("bad dots\n");
          std::abort ();
        }

      *x = - (b * f - c * e) / (b * d - e * a);
      *y =   (a * f - c * d) / (b * d - e * a);
    };

  root_of_lines (0., 0., x2, y2,
                 data.q[0].x0, data.q[0].y0,
                 data.q[0].x0 + (x2 - x1), data.q[0].y0 + (y2 - y1),
                 &data.q[0].x3, &data.q[0].y3);

  data.q[1].x0 = data.q[0].x3;
  data.q[1].y0 = data.q[0].y3;
  data.q[1].x1 = x2;
  data.q[1].y1 = y2;
  data.q[1].x2 = x3;
  data.q[1].y2 = y3;
  root_of_lines (0., 0., x3, y3,
                 data.q[1].x0, data.q[1].y0,
                 data.q[1].x0 + (x3 - x2), data.q[1].y0 + (y3 - y2),
                 &data.q[1].x3, &data.q[1].y3);

#ifdef DEBIG_PRINT
  auto print_q = [] (const quad& q)
    {
      printf ("x0 = %.2lf, y0 = %.2lf\n"
              "x1 = %.2lf, y1 = %.2lf\n"
              "x2 = %.2lf, y2 = %.2lf\n"
              "x3 = %.2lf, y3 = %.2lf\n",
              q.x0, q.y0,
              q.x1, q.y1,
              q.x2, q.y2,
              q.x3, q.y3);
    };
  print_q (data.q[0]);
  print_q (data.q[1]);
#endif

  data.q[0].dy0 = get_length (data.q[0].x0, data.q[0].y0,
                              data.q[0].x1, data.q[0].y1) / data.ny;
  data.q[0].dy1 = get_length (data.q[0].x2, data.q[0].y2,
                              data.q[0].x3, data.q[0].y3) / data.ny;
  data.q[1].dy0 = get_length (data.q[1].x0, data.q[1].y0,
                              data.q[1].x1, data.q[1].y1) / data.ny;
  data.q[1].dy1 = get_length (data.q[1].x2, data.q[1].y2,
                              data.q[1].x3, data.q[1].y3) / data.ny;
  return error::all_ok;
}

error parse_file (const char *file_name,
                  double *x1, double *y1,
                  double *x2, double *y2,
                  double *x3, double *y3, double *l)
{
  file_pointer fp (file_name, "r");
  if (fp.opened ())
    {
      constexpr int buf_size = 100;
      char buf[buf_size];

      while (fgets (buf, buf_size, fp) && (strchr (buf, '#') || strlen (buf) == 0));

      if (sscanf (buf, "%lf%lf", x1, y1) != 2)
        return error::cannot_read;

      while (fgets (buf, buf_size, fp) && (strchr (buf, '#') || strlen (buf) == 0));

      if (sscanf (buf, "%lf%lf", x2, y2) != 2)
        return error::cannot_read;

      while (fgets (buf, buf_size, fp) && (strchr (buf, '#') || strlen (buf) == 0));
      
      if (sscanf (buf, "%lf%lf", x3, y3) != 2)
        return error::cannot_read;

      while (fgets (buf, buf_size, fp) && (strchr (buf, '#') || strlen (buf) == 0));
      
      if (sscanf (buf, "%lf", l) != 1)
        return error::cannot_read;
    }
  else return error::cannot_open;

  return error::all_ok;
}
