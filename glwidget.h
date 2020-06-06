#ifndef _my_widget
#define _my_widget

#include <qgl.h>
#include <QKeyEvent>
#include <qnamespace.h>

#include "header.h"
#include "alg.h"

class GLWidget : public QGLWidget
{
Q_OBJECT

private:
  method_data *data;
  Functions &f;
  int method = 0;
  double min = 0.0;
  double max = 0.0;
  unsigned int max_nx = 2000;
  unsigned int max_ny = 2000;
  unsigned int min_nx = 2;
  unsigned int min_ny = 2;

  void draw_text (QPainter &painter, const char *str, int i = 0);
  void rebuild_approximation ();
  void count_min_max ();
  double count (unsigned int i, unsigned int j, double *c, const quad &q,
                double *xx = nullptr, double *yy = nullptr);
  double count_for_min_max (unsigned int i, unsigned int j,
                            double *c, const quad &q);
  double count_and_scale (unsigned int i, unsigned int j,
                          double *c, const quad &q,
                          double *xx = nullptr, double *yy = nullptr);
  double scale (double z);

public:
  GLWidget (method_data *d)
    :
      QGLWidget (QGLFormat (QGL::SampleBuffers)),
      data (d),
      f (data->f)
  {
    data->new_count.unlock ();
    data->end_count.lock ();
    count_min_max ();
  }
  ~GLWidget ()
  {
    data->end_program = true;
    data->new_count.unlock ();
  }

protected:
  virtual void paintGL ();
  virtual void initializeGL ();
  virtual void resizeGL (int nWidth, int nHeight);
  virtual void keyPressEvent (QKeyEvent *e);

  void setProjectionMatrix ();
  void setDefaultCamera ();

  double angle_h, angle_v;
  double camera_p;
  double aspect;
};

#endif
