#ifndef IMAT_BASE_INL_H
#define IMAT_BASE_INL_H

// ListInitializer
template <class T>
ListInitializer<T>::ListInitializer(T *begin, int length, T value)
    : it_(begin), end_(begin+length) {
  assert(it_ != end_);
  *it_ = value;
}

template <class T>
ListInitializer<T>::ListInitializer(T *begin, T *end, T value)
    : it_(begin), end_(end) {
  assert(it_ != end_);
  *it_ = value;
}

template <class T>
ListInitializer<T> ListInitializer<T>::operator, (T value) {
  ++it_;
  assert(it_ < end_);
  return ListInitializer<T> (it_, end_, value);
}

// template <class T>
// inline void AddAss(T *x, const T &y) { (*x) += y; }
// 
// template <class T>
// inline void SubAss(T *x, const T &y) { (*x) -= y; }
// 
// template <class T>
// inline void MulAss(T *x, const T &y) { (*x) *= y; }
// 
// template <class T>
// inline void DivAss(T *x, const T &y) { (*x) /= y; }
// 
// template <class T>
// inline T Add(const T &x, const T &y) { return x + y; }
// 
// template <class T>
// inline T Sub(const T &x, const T &y) { return x - y; }
// 
// template <class T>
// inline T Mul(const T &x, const T &y) { return x * y; }
// 
// template <class T>
// inline T Div(const T &x, const T &y) { return x / y; }

#endif // IMAT_BASE_INL_H
