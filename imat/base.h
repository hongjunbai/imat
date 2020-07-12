#ifndef IMAT_BASE_H
#define IMAT_BASE_H

#include <cassert>
#include <cstddef>
#include "macros.h"

namespace imat {

typedef long isize_t;

static const int kDefault = 0;        // Default value of non-negtive int

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TName) \
  TName(const TName&);               \
  void operator=(const TName&);

// Analog for ListInitializer in FLENs
// Initialize objects sequentially by operator,
template <class T>
class ListInitializer {
 public:
  ListInitializer(T *begin, int length, T value);
  ListInitializer<T> operator, (T value);
 private:
  ListInitializer(T *begin, T *end, T value);
  T *it_;
  T *end_;
};

// struct choose used to facilitate the definition of iterator
template <bool flag, class IsTrue, class IsFalse>
struct choose;

template <class IsTrue, class IsFalse>
struct choose <true, IsTrue, IsFalse> {
  typedef IsTrue type;
};

template <class IsTrue, class IsFalse>
struct choose <false, IsTrue, IsFalse> {
  typedef IsFalse type;
};

// Not used by any part of array:
// template <class T> void AddAss (T *x, const T &y);
// template <class T> void SubAss (T *x, const T &y);
// template <class T> void MulAss (T *x, const T &y);
// template <class T> void DivAss (T *x, const T &y);
// template <class T> T Add(const T &x, const T &y);
// template <class T> T Sub(const T &x, const T &y);
// template <class T> T Mul(const T &x, const T &y);
// template <class T> T Div(const T &x, const T &y);

#include "base_inl.h"

} // namespace imat
#endif // IMAT_BASE_H
