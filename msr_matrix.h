#ifndef _MSR_MATRIX_H_
#define _MSR_MATRIX_H_

#include <stdlib.h>
#include <memory>
#include <vector>
#include "header.h"

class msr_matrix
{
private:
  std::unique_ptr<unsigned int []> rnz;
  std::unique_ptr<double []> data;
  std::unique_ptr<double []> rhs;
  unsigned int n_rows = 0;
  unsigned int rnz_size = 0;
  unsigned int max_print = 10;

public:
  msr_matrix () {}
  msr_matrix (unsigned int n_rows, unsigned int offdiag_size)
    {
      init (n_rows, offdiag_size);
    }

  error init (unsigned int n, unsigned int offdiag_size)
    {
      if (n < 1)
        return error::bad_data;

      if (n_rows < n || rnz_size < offdiag_size + n_rows + 1)
        {
          n_rows = n;
          rnz_size = n_rows + 1 + offdiag_size;
          rnz = std::make_unique<unsigned int []> (rnz_size);
          data = std::make_unique<double []> (rnz_size);
          rhs = std::make_unique<double []> (rnz_size);

          if (   rnz.get () == nullptr
              || data.get () == nullptr
              || rhs.get () == nullptr)
            {
              return error::cannot_alloc;
            }
        }
      else
        {
          n_rows = n;
          rnz_size = n_rows + 1 + offdiag_size;
        }

      rnz[0] = n_rows + 1;

      return error::all_ok;
    }

  unsigned int get_n_rows () const { return n_rows; }
  unsigned int *get_internal_rnz () const { return rnz.get (); }
  double *get_internal_data () { return data.get (); }
  double *get_internal_rhs () { return rhs.get (); }
  const unsigned int *get_rnz () const { return rnz.get (); }
  const double *get_data () const { return data.get (); }
  const double *get_rhs () const { return rhs.get (); }
  
  void set_row_len (unsigned int i_row, unsigned int len)
    {
      rnz[i_row + 1] = rnz[i_row] + len;
    }
  void set_col (unsigned int i_row, unsigned int offset, unsigned int j_col)
    {
      rnz[rnz[i_row] + offset] = j_col;
    }

  unsigned int get_offdiag_len (unsigned int i_row_begin,
                                unsigned int i_row_end) const
    {
      return rnz[i_row_end] - rnz[i_row_begin];
    }
  unsigned int get_offdiag_len (unsigned int i_row) const
    {
      return get_offdiag_len (i_row, i_row + 1);
    }
  unsigned int get_offdiag_len () const
    {
      return get_offdiag_len (0, n_rows);
    }
  unsigned int get_col_num (unsigned int i_row, unsigned int j_col) const
    {
      return rnz[rnz[i_row] + j_col];
    }
  int find_col_num (unsigned int i_row, unsigned int col) const
    {
      auto ptr_begin = rnz.get () + rnz[i_row];
      auto ptr_end = rnz.get () + rnz[i_row + 1];
      auto len = ptr_end - ptr_begin;
      auto offdiag_i = std::lower_bound (ptr_begin, ptr_end, col) - ptr_begin;
      return offdiag_i < len && ptr_begin[offdiag_i] == col
              ? offdiag_i : -1;
    }

  double get_diag (unsigned int i_row) const
    {
      return data[i_row];
    }
  void set_diag (unsigned int i_row, double v)
    {
      data[i_row] = v;
    }
  double get_offdiag (unsigned int i_row, unsigned int j_col) const
    {
      return data[rnz[i_row] + j_col];
    }
  void set_offdiag (unsigned int i_row, unsigned int j_col, double v)
    {
      data[rnz[i_row] + j_col] = v;
    }

  void print (FILE *fp, const char *prefix = nullptr);
  error print (const char *file_name, const char *prefix = nullptr);
};

#endif // _MSR_MATRIX_H_
