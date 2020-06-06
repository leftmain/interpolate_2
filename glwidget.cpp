#include "glwidget.h"
#include "alg.h"

#include <stdio.h>
#include <math.h>

void GLWidget::initializeGL()
{
  glClearColor(1.0, 1.0, 1.0, 1.0);

  setDefaultCamera();
}

void GLWidget::draw_text (QPainter &painter, const char *str, int ih)
{
  const int hh = 20;
  const int ww = 300;
  ih++;
  painter.drawText(5, hh * (ih - 1), ww, hh * ih, Qt::AlignLeft, str);
}

// scale z from [min,max] to [1e-4,1]
double GLWidget::scale (double z)
{
  double eps = (method == 2) ? 1e-16 : 1e-13;
  if (fabs (max - min) < eps)
    {
      if (z < eps)
        return 0.0;
      return 1.0;
    }
  return 1e-4 + (z - min) * (1 - 1e-4) / (max - min);
}

double GLWidget::count (unsigned int i, unsigned int j, double *c, const quad &q,
                        double *xx, double *yy)
{
  double x, y, z;
  auto nx = data->nx;
  auto ny = data->ny;

  get_xy_by_ij (q, nx, ny, &x, &y, i, j);

  if (xx != nullptr && yy != nullptr)
    {
      *xx = x;
      *yy = y;
    }

  if (method % 3 == 0)
    z = f (x, y);
  else if (method % 3 == 1)
    z = c[get_k_by_ij (nx, ny, i, j)];
  else
    z = fabs (f (x, y) - c[get_k_by_ij (nx, ny, i, j)]);

  return z;
}

double GLWidget::count_for_min_max (unsigned int i, unsigned int j,
                                    double *c, const quad &q)
{
  double x, y, z;
  auto nx = data->nx;
  auto ny = data->ny;

  get_xy_by_ij (q, nx, ny, &x, &y, i, j);

  if (method % 3 == 2)
    z = fabs (f (x, y) - c[get_k_by_ij (nx, ny, i, j)]);
  else
    z = f (x, y);

  return z;
}

double GLWidget::count_and_scale (unsigned int i, unsigned int j,
                                  double *c, const quad &q,
                                  double *xx, double *yy)
{
#if 0
  double z = scale (count (i, j, c, q, xx, yy));
  if (z < 0)
    z = 0.0;
  return z;
#else
  return scale (count (i, j, c, q, xx, yy));
#endif
}

void GLWidget::count_min_max ()
{
  auto nx = data->nx;
  auto ny = data->ny;

  auto count_on_quad = [&] (const quad &q, double *c)
    {
      double z = 0.0;
      //for (unsigned int i = 1; i < nx; ++i)
      for (unsigned int i = 0; i <= nx; ++i)
        {
          //for (unsigned int j = 1; j < ny; ++j)
          for (unsigned int j = 0; j <= ny; ++j)
            {
              z = fabs (count_for_min_max (i, j, c, q));
              if (z > max) max = z;
              if (z < min) min = z;
            }
        }
    };

  max = min = fabs (count_for_min_max (0, 0, data->c[0], data->q[0]));
  max = 0.0;
  min = 1.0;
  for (int k = 0; k < 2; ++k)
    count_on_quad (data->q[k], data->c[k]);

  //printf ("min = %le, max = %le\n", min, max);
}

void GLWidget::paintGL()
{
  setProjectionMatrix ();

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnable (GL_DEPTH_TEST);

  glColor3d (0.0, 0.0, 0.0);

  // draw text
  char buf[300];
  QPainter painter(this);
  painter.setPen(Qt::black);
  painter.setFont(QFont("Arial", 12));

  draw_text (painter, f.get_name (), 0);

  if (method % 3 == 0)
    sprintf (buf, "mode: function");
  else if (method % 3 == 1)
    sprintf (buf, "mode: approx");
  else if (method % 3 == 2)
    sprintf (buf, "mode: residual");
  draw_text (painter, buf, 1);

  sprintf (buf, "%le", data->residual);
  draw_text (painter, buf, 2);

  sprintf (buf, "nx = %u, ny = %u", data->nx, data->ny);
  draw_text (painter, buf, 3);

#if 1
  sprintf (buf, "min = %.2le, max = %.2le", min, max);
  draw_text (painter, buf, 4);
#endif

  painter.end();
  // end drawing text

  // draw axis
  auto draw_axis = [&] (int i, int size = 3.0)
    {
      double x[] = {0.0, 0.0, 0.0};
      x[i] = size;

      glBegin (GL_LINES);
        {
          double xx = 0.0;
          glVertex3d (xx, xx, xx);
          glVertex3d (x[0], x[1], x[2]);
        }
      glEnd ();
    };
  glColor3d (0.0, 0.0, 0.0);
  draw_axis (0);
  draw_axis (1);
  draw_axis (2);

  auto draw_quad = [&] (const quad &q, double *c)
    {
      auto draw_quad_by_dots = [&] (const quad &q, int type)
        {
          auto nx = data->nx;
          auto ny = data->ny;

          auto add_point = [&] (unsigned int i, unsigned int j)
            {
              double x, y, z;
              z = count_and_scale (i, j, c, q, &x, &y);
              glVertex3d (x, y, z);
            };
          unsigned int offset = 1;
          if (nx > 70 || ny > 70)
            offset = nx / 80;

          for (unsigned int i = 1; i <= ny; i += offset)
            {
              glBegin (type);
              for (unsigned int j = 0; j < nx; j += offset)
                {
                  add_point (i, j);
                  add_point (i - 1, j);
                  add_point (i - 1, j + 1);
                  add_point (i, j);
                  add_point (i, j + 1);
                  add_point (i - 1, j + 1);
                }
              glEnd ();
            }
        };

      glColor3d (0.0, 0.0, 1.0);
      draw_quad_by_dots (q, GL_LINE_STRIP);
    };

  for (int k = 0; k < 2; ++k)
    draw_quad (data->q[k], data->c[k]);

  glDisable (GL_DEPTH_TEST);
}

void GLWidget::resizeGL (int nWidth, int nHeight)
{
  glViewport (0, 0, nWidth, nHeight);
  aspect = ((double) nWidth) / nHeight;
  updateGL ();
}

void GLWidget::rebuild_approximation ()
{
  data->end_count.try_lock ();
  data->new_count.unlock ();
  data->end_count.lock ();
  count_min_max ();
}

void GLWidget::keyPressEvent (QKeyEvent* event)
{
  auto recount_params = [&] ()
    {
      data->q[0].dy0 = get_length (data->q[0].x0, data->q[0].y0,
                                data->q[0].x1, data->q[0].y1) / data->ny;
      data->q[0].dy1 = get_length (data->q[0].x2, data->q[0].y2,
                                data->q[0].x3, data->q[0].y3) / data->ny;
      data->q[1].dy0 = get_length (data->q[1].x0, data->q[1].y0,
                                 data->q[1].x1, data->q[1].y1) / data->ny;
      data->q[1].dy1 = get_length (data->q[1].x2, data->q[1].y2,
                                 data->q[1].x3, data->q[1].y3) / data->ny;
    };

  switch (event->key ()) {
    case Qt::Key_Q:
      close ();
      break;
    case Qt::Key_C:
      setDefaultCamera ();
      break;
    case Qt::Key_Up:
      angle_v = std::min (angle_v + 5.0, 80.0);
      break;
    case Qt::Key_Down:
      angle_v = std::max (angle_v - 5.0, -80.0);
      break;
    case Qt::Key_Left:
      angle_h -= 5.0;
      break;
    case Qt::Key_Right:
      angle_h += 5.0;
      break;
    case Qt::Key_Plus:
    case Qt::Key_2:
      camera_p = std::max (camera_p - 0.1, 6.0);
      break;
    case Qt::Key_Minus:
    case Qt::Key_3:
      camera_p += 0.1;
      break;
    case Qt::Key_0:
      f.next ();
      rebuild_approximation ();
      break;
    case Qt::Key_1:
      method = (method + 1) % 3;
      count_min_max ();
      break;
    case Qt::Key_4:
      if (data->nx * 2 <= max_nx || data->ny * 2 <= max_ny)
        {
          data->nx *= 2;
          data->ny *= 2;
          recount_params ();
          rebuild_approximation ();
        }
      break;
    case Qt::Key_5:
      if (data->nx / 2 >= min_nx || data->ny / 2 >= min_ny)
        {
          data->nx /= 2;
          data->ny /= 2;
          recount_params ();
          rebuild_approximation ();
        }
      break;
  }

  updateGL ();
}

void GLWidget::setDefaultCamera ()
{
  camera_p = 7;
  angle_h = 40;//45;
  angle_v = 20;//20;
  aspect = 1.0 * width () / height ();
}

void GLWidget::setProjectionMatrix ()
{
  GLdouble view[16] = {0}, projection[16] = {0}, tmp[16] = {0};

  static GLdouble near = 5, top = 2, bottom, left, right;

  bottom = -top;
  right = top * aspect;
  left = -right;

  projection[0] = 2 * near / (right - left);
  projection[2] = (right + left) / (right - left);
  projection[5] = 2 * near / (top - bottom);
  projection[6] = (top + bottom) / (top - bottom);
  projection[10] = - 1;
  projection[11] = - 2 * near;
  projection[14] = -1;

  GLdouble cam_x, cam_y, cam_z;
  cam_x = 0;
  cam_y = 0;
  cam_z = camera_p;

  view[0] = 1;
  view[6] = -1;
  view[9] = 1;
  view[15] = 1;

  view[12] = -cam_x;
  view[13] = -cam_y;
  view[14] = -cam_z;

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  glRotated (angle_h, 0, 0, 1);
  glGetDoublev (GL_PROJECTION_MATRIX, tmp);

  glLoadTransposeMatrixd (projection);
  glMultMatrixd (view);

  glRotated (angle_h, 0, 0, 1);
  glRotated (angle_v, tmp[0], tmp[4], tmp[8]);
}

