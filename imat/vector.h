#ifndef IMAT_VECTOR_H
#define IMAT_VECTOR_H

#include <cstdarg>
#include <cstring>
#include <cmath>
#include <cassert>

#include <iostream>

#include "base.h"

namespace imat {

// Vector iterator
template <class T, bool isconst>
struct VecIterator {
  typedef typename choose<isconst, const T&, T&>::type reference;
  typedef typename choose<isconst, const T*, T*>::type pointer;

  VecIterator()                        : p(NULL), stride(0)      {}
  VecIterator(T* ptr, isize_t xstride) : p(ptr), stride(xstride) {}
  VecIterator(const VecIterator<T, false> &src);
  VecIterator& operator=(const VecIterator<T, false> &src);

  reference operator*()  const { return *(this->p); }
  pointer   operator->() const { return this->p;    }
  bool operator==(const VecIterator<T, isconst> &other) const;
  bool operator!=(const VecIterator<T, isconst> &other) const;
  VecIterator<T, isconst>& operator++();
  VecIterator<T, isconst>& operator+=(int steps);
  VecIterator<T, isconst>& operator-=(int steps);

  // data
  T *p;
  isize_t stride;
};

// Vector Class:
// Provide a C++ Vector class, whoes data can be exposed to some C numeric
// libraries (eg. GSL lib)
template <class T>
class Vector  {
  // Following functions are not as efficient as operate_and_assign (+=, -= and
  // *=). Try to use operate_and_assign instead of these functions whenever
  // possible.
  template <class TT> 
  friend Vector<TT> operator- (const Vector<TT> &vect1, const Vector<TT> &vect2);
  template <class TT> 
  friend Vector<TT> operator+ (const Vector<TT> &vect1, const Vector<TT> &vect2);
  template <class TT> 
  friend Vector<TT> operator* (const Vector<TT> &vect1, const Vector<TT> &vect2);
  template <class TT> 
  friend Vector<TT> operator* (const Vector<TT> &vec, const TT &x);
  template <class TT> 
  friend Vector<TT> operator* (const TT &x, const Vector<TT> &vec);
  //
  template <class TT> 
  friend TT DotProduct(const Vector<TT> &vect1, const Vector<TT> &vect2);

 public:
  Vector() : size_(0), stride_(0), data_(NULL), owner_(false) {}
  ~Vector() { if (owner_) delete []data_; }
  explicit Vector(isize_t len);
  Vector(const Vector<T> &source);
  Vector(isize_t length, T first, ...);
  Vector<T> & operator= (const Vector<T> &source);
  ListInitializer<T> operator= (const T &value);
  void Allocate(isize_t len);
  void View(T *source, isize_t size, isize_t stride=1);
  void SubCopy(const Vector<T> &src,
               isize_t i0 = kDefault,
               isize_t num = kDefault);
  //
  T &       operator() (isize_t i);
  const T & operator() (isize_t i) const;
  Vector<T> operator-() const;
  bool operator== (const Vector<T> &other) const;
  void operator+= (const Vector<T> &other);
  void operator-= (const Vector<T> &other);
  void operator*= (const Vector<T> &other);
  void operator+= (const T &x);
  void operator-= (const T &x);
  void operator*= (const T &x);
  //
  void Print() const;
  void SetZero();
  void SetAll(const T &x);
  void SetBasis(isize_t i);
  T Sum() const;
  T Min() const;
  T Max() const;
  double Mean() const;
  void Norm();
  T Product() const;
  //
  bool     empty()  const { return (size_ == 0); }
  isize_t  size()   const { return size_;        }
  isize_t  stride() const { return stride_;      }
  T*       data()         { return data_;        }
  const T* data()   const { return data_;        }
  void     clear();
  //
  typedef VecIterator<T, false> iterator;
  typedef VecIterator<T, true> const_iterator;
  iterator begin() { return iterator(data_, stride_);                 }
  iterator end()   { return iterator(data_ + size_*stride_, stride_); }
  const_iterator begin() const { return const_iterator(data_, stride_); }
  const_iterator end()   const ;

 protected:
  bool is_view_() const { return (data_ != NULL and owner_ == false); }
  void copy_(const Vector<T> &source);
  void allocate_(isize_t len);

  isize_t size_;
  isize_t stride_;
  T *data_;
  bool owner_;
};

typedef Vector<int> VectorI;
typedef Vector<double> VectorD;

// Specialized functions
#include "vector_inl.h"

} // namespace imat 
#endif // IMAT_VECTOR_H
