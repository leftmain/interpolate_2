#ifndef _HEADER_H_
#define _HEADER_H_

#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <algorithm>
#include <mutex>
#include <condition_variable>

#define DRAW_GRAPH
#define DEBUG_PRINT

enum class error : int
{
  all_ok = 0,
  cannot_open = 1,
  cannot_read = 2,
  cannot_alloc = 3,
  bad_data = 4,
};

class file_pointer
{
private:
  FILE *fp = nullptr;

public:
  file_pointer (const char *file_name, const char *rw)
    {
      fp = fopen (file_name, rw);
    }
  ~file_pointer ()
    {
      if (fp)
        fclose (fp);
    }
  operator const FILE *() const { return fp; }
  operator FILE *() const { return fp; }
  bool opened () const
    {
      return (fp == nullptr) ? false : true;
    }
};

typedef double (*func_type) (double, double);

struct Func
{
  func_type f;

  double operator () (double x, double y) const
    {
      return f (x, y);
    }
  void operator= (const Func &ff)
    {
      f = ff.f;
    }
};

struct Function
{
  const char *name;
  Func f;

  double operator () (double x, double y) const
    {
      return f (x, y);
    }
  const Func &get () const
    {
      return f;
    }
};

class Functions
{
private:
  const Function *funcs = nullptr;
  int id = 0;
  int n = 0;

public:
  Functions (const Function *f, int id_, int n_) : funcs (f), id (id_), n (n_) {}
  Functions () {}

  void next ()
    {
      id = (id + 1) % n;
    }

  void init (const Function *f, int i, int nn)
    {
      funcs = f;
      id = i;
      n = nn;
    }

  const Func &get () const
    {
      return funcs[id].f;
    }

  const char * get_name () const
    {
      return funcs[id].name;
    }

  int get_id () const
    {
      return id;
    }

  double operator () (double x, double y) const
    {
      return funcs[id].f (x, y);
    }
};


#endif // _HEADER_H_

