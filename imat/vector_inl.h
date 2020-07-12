#ifndef IMAT_VECTOR_INL_H
#define IMAT_VECTOR_INL_H

// Vector iterator
// Member functions
template <class T, bool isconst>
VecIterator<T, isconst>::VecIterator(const VecIterator<T, false> &src) 
    : p(src.p), stride(src.stride) { };

template <class T, bool isconst>
inline VecIterator<T, isconst>& 
VecIterator<T, isconst>::operator=(const VecIterator<T, false> &src) {
  p = src.p;
  stride = src.stride;
  return *this;
}

template <class T, bool isconst>
inline bool VecIterator<T, isconst>::operator==(
    const VecIterator<T, isconst> &other) const {
  return this->p == other.p and this->stride == other.stride;
}

template <class T, bool isconst>
inline bool VecIterator<T, isconst>::operator!=(
    const VecIterator<T, isconst> &other) const {
  return this->p != other.p or this->stride != other.stride;
}

template <class T, bool isconst>
inline VecIterator<T, isconst>& VecIterator<T, isconst>::operator++() { 
  p += stride;
  return *this; 
}

template <class T, bool isconst>
inline VecIterator<T, isconst>& 
VecIterator<T, isconst>::operator+=(int steps) {
  p += stride*steps;
  return *this; 
}

template <class T, bool isconst>
inline VecIterator<T, isconst>& VecIterator<T, isconst>::operator-=(
    int steps) {
  p -= stride*steps;
  return *this; 
}

// Vector 
// Friend functions
template <class T> 
Vector<T> operator- (const Vector<T> &vect1, const Vector<T> &vect2) {
  assert(vect1.size() == vect2.size());
  Vector<T> result(vect1.size());
  MFOR(i, vect1.size()) result(i) = vect1(i) - vect2(i);
  return result;
}

template <class T> 
Vector<T> operator+ (const Vector<T> &vect1, const Vector<T> &vect2) {
  assert(vect1.size() == vect2.size());
  Vector<T> result(vect1.size());
  MFOR(i, vect1.size()) result(i) = vect1(i) + vect2(i);
  return result;
}

template <class T> 
Vector<T> operator* (const Vector<T> &vect1, const Vector<T> &vect2) {
  assert(vect1.size() == vect2.size());
  Vector<T> result(vect1.size());
  MFOR(i, vect1.size()) result(i) = vect1(i) * vect2(i);
  return result;
}

template <class T> 
Vector<T> operator* (const Vector<T> &vect, const T &x) {
  Vector<T> result(vect.size());
  MFOR(i, vect.size()) result(i) = vect(i) * x;
  return result;
}

template <class T> 
Vector<T> operator* (const T &x, const Vector<T> &vect) {
  Vector<T> result(vect.size());
  MFOR(i, vect.size()) {
    result(i) = vect(i) * x;
  }
  return result;
}

template <class T> 
T DotProduct(const Vector<T> &vect1, const Vector<T> &vect2) {
  assert(vect1.size() == vect2.size());
  T accumlator = 0;
  MFOR(i, vect1.size()) accumlator += vect1(i) * vect2(i);
  return accumlator;
}

// Vector
template <class T>
Vector<T>::Vector(isize_t len) {
  allocate_(len);
}

template <class T>
Vector<T>::Vector(const Vector<T> &source) {
  copy_(source);
}

template <class T>
Vector<T>::Vector(isize_t length, T first, ...) {
  assert(length > 0);
  allocate_(length);
  va_list args;
  T val;
  va_start(args, first);
  this->data_[0] = first;
  MFOR(i, length) {
    val = va_arg(args, T);
    this->data_[i] = val;
  }
  va_end(args);
}

template <class T>
Vector<T> & Vector<T>::operator= (const Vector<T> &source) {
  if (this == &source) return *this;
  if (owner_) delete []data_;
  copy_(source);
  return *this;
}

template <class T>
ListInitializer<T> Vector<T>::operator= (const T &value) {
  assert(owner_ or stride_ == 1); // Vector
  return ListInitializer<T>(data_, size_, value);
}

template <class T>
void Vector<T>::Allocate(isize_t len) {
  assert(empty());
  allocate_(len);
}

template <class T>
inline void Vector<T>::View(T *source, isize_t size, isize_t stride) {
  assert(empty());
  data_ = source;
  size_ = size;
  stride_ = stride;
  owner_ = false;
}

template <class T>
void Vector<T>::SubCopy(const Vector<T> &src, isize_t i0, isize_t num) {
  if (num == kDefault) num = src.size_;
  assert(src.size_ >= num);
  assert(size_ >= num);
  MFOR(i, num) (*this)(i) = src(i0 + i);
}

template <class T>
inline T & Vector<T>::operator() (isize_t i) {
  assert(i>=0 and i<size_);
  return *(data_ + i * stride_);
}

template <class T>
inline const T & Vector<T>::operator() (isize_t i) const {
  assert(i>=0 and i<size_);
  return *(data_ + i * stride_);
}

template <class T>
Vector<T> Vector<T>::operator-() const {
  Vector<T> result(size_);
  MFOR(i, size_) result(i) = -(*this)(i);
  return result;
}

template <class T>
bool Vector<T>::operator== (const Vector<T> &other) const  {
  bool is_equal = true;
  if (size_ != other.size_) {
    is_equal = false;
  } else {
    MFOR(i, size_) {
      if ((*this)(i) != other(i)) {
        is_equal = false;
        break;
      }
    }
  }
  return is_equal;
}

template <class T>
void Vector<T>::operator+= (const Vector<T> &other) {
  assert(size_ == other.size_);
  MFOR(i, size_) (*this)(i) += other(i);
}

template <class T>
void Vector<T>::operator-= (const Vector<T> &other) {
  assert(size_ == other.size_);
  MFOR(i, size_) (*this)(i) -= other(i);
}

template <class T>
void Vector<T>::operator*= (const Vector<T> &other) {
  assert(size_ == other.size_);
  MFOR(i, size_) (*this)(i) *= other(i);
}

template <class T>
void Vector<T>::operator+= (const T &x) {
  MFOR(i, size_) (*this)(i) += x;
}

template <class T>
void Vector<T>::operator-= (const T &x) {
  MFOR(i, size_) (*this)(i) -= x;
}

template <class T>
void Vector<T>::operator*= (const T &x) {
  MFOR(i, size_) (*this)(i) *= x;
}

template <class T>
void Vector<T>::Print() const {
  std::cout << "[";
  MFOR(i, size_-1) std::cout << (*this)(i) << ", ";
  std::cout << (*this)(size_-1) << "]" << std::endl;
}

template <class T>
void Vector<T>::SetZero() {
  MFOR(i, size_) (*this)(i) = 0;
}

template <class T>
void Vector<T>::SetAll(const T &x) {
  MFOR(i, size_) (*this)(i) = x;
}

template <class T>
void Vector<T>::SetBasis(isize_t i) {
  SetZero();
  (*this)(i) = 1;
}

template <class T>
T Vector<T>::Sum() const {
  T accumlator = (*this)(0);
  MRANGE (i, 1, size_) accumlator += (*this)(i);
  return accumlator;
}

template <class T>
T Vector<T>::Min() const {
  T min= (*this)(0);
  MRANGE (i, 1, size_) min = min < (*this)(i) ? min : (*this)(i);
  return min;
}

template <class T>
T Vector<T>::Max() const {
  T max= (*this)(0);
  MRANGE (i, 1, size_) max = max > (*this)(i) ? max : (*this)(i);
  return max;
}

template <class T>
double Vector<T>::Mean() const {
  return double(Sum()) / size_;
}

template <class T>
void Vector<T>::Norm() {
  MSSHOW("Only Vector<double> could be normalized.");
}

template <>
inline void Vector<double>::Norm() {
  double accumlator = 0;
  MFOR(i, size_) accumlator += (*this)(i) * (*this)(i);
  accumlator = sqrt(accumlator);
  MFOR(i, size_) (*this)(i) /= accumlator;
}

template <class T>
T Vector<T>::Product() const {
  T accumlator = (*this)(0);
  MRANGE(i, 1, size_) accumlator *= (*this)(i);
  return accumlator;
}

template <class T>
void Vector<T>::clear() {
  if (owner_) delete []data_;
  data_ = NULL;
  size_ = stride_ = 0;
  owner_ = false;
}

template <class T>
inline typename Vector<T>::const_iterator Vector<T>::end() const {
  return const_iterator(data_ + size_ * stride_, stride_);
}

template <class T>
void Vector<T>::copy_(const Vector<T> &source) {
 allocate_(source.size_);
  if (! source.is_view_()) { // Vector
    memcpy(data_, source.data_, sizeof(T) * size_);
  } else {                   // Souce is VectorView or EmptyVecotr
    MFOR(i, size_) (*this)(i) = source(i);
  }
}

template <class T>
inline void Vector<T>::allocate_(isize_t len) {
  assert(len> 0);
  size_ = len;
  stride_ = 1;
  data_ = new T [size_];
  assert(data_);
  owner_ = true;
}

#endif // IMAT_VECTOR_INL_H
