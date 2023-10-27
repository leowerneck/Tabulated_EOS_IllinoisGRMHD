#ifndef GRHayL_test_H_
#define GRHayL_test_H_

static inline int indexf(const int gridmax, const int i, const int j, const int k) {
  return i + j*gridmax + k*gridmax*gridmax;
}

#define check_file_was_successfully_open(fp, filename) \
  if(!fp) CCTK_VERROR("Could not open file %s.\n", filename);

#endif // GRHayL_test_H
