#ifndef IMAT_ARRAY3D_INL_H
#define IMAT_ARRAY3D_INL_H

template <class T>
Array3D<T>::Array3D(const T *source,
                    const isize_t rows, const isize_t cols,
                    const isize_t layers) {
  xsize = rows;
  ysize = cols;
  zsize = layers;
  MFOR(i, xsize*ysize*zsize) data.push_back(*(source + i));
}

template <class T>
void Array3D<T>::Allocate(const isize_t xsize0,
                          const isize_t ysize0,
                          const isize_t zsize0) {
  assert(xsize == 0 && ysize == 0 && zsize == 0);
  xsize = xsize0;
  ysize = ysize0;
  zsize = zsize0;
  data.resize(xsize * ysize * zsize);
}

template <class T>
inline T & Array3D<T>::operator() (const isize_t i,
                                   const isize_t j,
                                   const isize_t k) {
  assert(i>=0 and i<xsize);
  assert(j>=0 and j<ysize);
  assert(k>=0 and k<zsize);
  // return data[i * (ysize * zsize) + j * zsize + k];
  return data[(i * ysize + j) * zsize + k];
}

template <class T>
inline const T & Array3D<T>::operator() (const isize_t i,
                                         const isize_t j,
                                         const isize_t k) const {
  assert(i>=0 and i<xsize);
  assert(j>=0 and j<ysize);
  assert(k>=0 and k<zsize);
  // return data[i * (ysize * zsize) + j * zsize + k];
  return data[(i * ysize + j) * zsize + k];
}

// Make sure to use this when T == numerical type
template <class T>
void Array3D<T>::SetZero() {
  MFOR(i, xsize*ysize*zsize) { data[i] = 0; }
}

template <class T>
void Array3D<T>::SetAll(const T &x) {
  MFOR(i, xsize*ysize*zsize) { data[i] = x; }
}

template <class T>
void Array3D<T>::operator*= (const Array3D<T> &other) {
  assert(this->xsize == other.xsize);
  assert(this->ysize == other.ysize);
  assert(this->zsize == other.zsize);
  MFOR(i, xsize*ysize*zsize) { this->data[i] *= other.data[i]; }
}

#endif // IMAT_ARRAY3D_INL_H
