#include "msr_matrix.h"

void msr_matrix::print (FILE *fp, const char *prefix)
{
  if (prefix)
    {
      fprintf (fp, "%s:\n", prefix);
    }
  auto n_rows = std::min (msr_matrix::max_print, msr_matrix::n_rows);
  for (unsigned int i_row = 0; i_row < n_rows; ++i_row)
    {
      for (unsigned int j_col = 0; j_col < n_rows; ++j_col)
        {
          if (j_col == i_row)
            fprintf (fp, "%.4le\t", get_diag (i_row));
          else
            {
              auto j = find_col_num (i_row, j_col);
              if (j >= 0)
                {
                  fprintf (fp, "%.4le\t", get_offdiag (i_row, j));
                }
              else
                {
                  fprintf (fp, "%.4le\t", 0.0);
                }
            }
        }
      if (rhs.get () != nullptr)
        {
          printf ("|\t%.3le", rhs[i_row]);
        }
      printf ("\n");
    }
  printf ("---------------------------------------------------------\n");
}

error msr_matrix::print (const char *file_name, const char *prefix)
{
  file_pointer fp (file_name, "w");
  if (fp.opened ())
    {
      print (fp, prefix);
    }
  else
    {
      return error::cannot_open;
    }
  return error::all_ok;
}

