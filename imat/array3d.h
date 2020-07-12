///////////////////////////////////////////////////////////////////////////////
//   array3d.h: for general 3d array handle
//      author: hjbai@mdl.ipc.pku.edu.cn
//     created: 2008-07-09
//    finished: ...
//     version: 0.0
///////////////////////////////////////////////////////////////////////////////
#ifndef IMAT_ARRAY3D_H
#define IMAT_ARRAY3D_H
#include <vector>
#include "base.h"

namespace imat{

template <class T>
struct Array3D {
  Array3D() : xsize(0), ysize(0), zsize(0), data(NULL) {}
  Array3D(const T *source,
          const isize_t xsize, const isize_t ysize, const isize_t zsize);
  // ~Array3D() {}
  //
  T &       operator() (isize_t i, isize_t j, isize_t k);
  const T & operator() (isize_t i, isize_t j, isize_t k) const;
  void      operator*= (const Array3D<T> &other);
  //
  void Allocate(const isize_t xsize, const isize_t ysize, const isize_t zsize);
  void SetZero();
  void SetAll(const T &x);
  //
  isize_t xsize;
  isize_t ysize;
  isize_t zsize;
  std::vector<T> data;
};

#include "array3d_inl.h"

} // namespace imat
#endif // IMAT_ARRAY3D_H
