#ifndef IMAT_MATRIX_H
#define IMAT_MATRIX_H

#include <cstdarg>
#include <cstdlib>

#include <fstream>
#include <sstream>

#ifdef GSL
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#endif // GSL

//#include "split_str.h"
#include "vector.h"
#include "base.h"

namespace imat {

enum Dim { kByRow = 0, kByCol = 1 };  // Choose dimention for Min, Max, Mean
                                      // Sum, Covariance and Sort

// Matrix Class:
// Provide a C++ Matrix class, whoes data can be exposed to C numeric libraries
// (eg. GSL lib).
template <class T>
class Matrix {
  // Following functions are not as efficient as operate_and_assign (+=, -= and
  // *=). Try to use operate_and_assign instead of these functions whenever
  // possible.
  template <class TT> 
  friend Matrix<TT> operator- (const Matrix<TT> &lfs, const Matrix<TT> &rhs);
  template <class TT> 
  friend Matrix<TT> operator+ (const Matrix<TT> &lfs, const Matrix<TT> &rhs);
  template <class TT> 
  friend Matrix<TT> operator* (const Matrix<TT> &lfs, const Matrix<TT> &rhs);
  template <class TT> 
  friend Matrix<TT> operator* (const Matrix<TT> &mat, const TT &x);
  template <class TT> 
  friend Matrix<TT> operator* (const TT &x, const Matrix<TT> &mat);
  // 
  // Cross Product
  template <class TT>
  friend void Cross(Matrix<TT> *res, const Matrix<TT> &a, const Matrix<TT> &b);

 public:
  Matrix() : rows_(0), cols_(0), tda_(0), data_(NULL), owner_(false) {}
  ~Matrix()  { if (owner_) delete []data_; }
  Matrix(isize_t rows, isize_t cols);
  Matrix(const Matrix <T> &source);
  Matrix<T> & operator= (const Matrix<T> &source);
  ListInitializer<T> operator= (const T &value);
  void Allocate(isize_t rows, isize_t cols); // for empty matrix only
  // Matrix views: can only be created by mat.View(...)
  // View raw data
  void View(T *src, isize_t rows, isize_t cols, isize_t tda=kDefault);
  // View submatrix: const and non-const version
  void View(Matrix<T> *src,
            isize_t row_0,         isize_t col_0,
            isize_t rows=kDefault, isize_t cols=kDefault);
  void View(const Matrix<T> &src,
            isize_t row_0,         isize_t col_0,
            isize_t rows=kDefault, isize_t cols=kDefault);
  // Copy elements from a submatrix of src
  void SubCopy(const Matrix<T> &src,
               isize_t row_0=kDefault, isize_t col_0=kDefault,
               isize_t rows=kDefault,  isize_t cols=kDefault);

  T &       operator() (isize_t i, isize_t j);
  const T & operator() (isize_t i, isize_t j) const;
  Matrix<T> operator-() const;
  bool operator== (const Matrix<T> &other) const;
  void operator+= (const Matrix<T> &other);
  void operator-= (const Matrix<T> &other);
  void operator*= (const Matrix<T> &other);
  void operator+= (const T &x);
  void operator-= (const T &x);
  void operator*= (const T &x);

  void SetZero();           // Set all elements as 0
  void SetAll(const T &x);  // Set all elements as x
  void SetUnit();           // Set square matrix as unit matrix
  T Sum() const;            // Sumation of all elements
  T Min() const;            // Minimum of all elements
  T Max() const;            // Maximum of all elements
  double Mean() const;      // Mean value of all elements
  void GetRow(Vector<T> *row_view, isize_t i);
  void GetCol(Vector<T> *col_view, isize_t i);
  void Trans();

  // Following functions are designed to operate on rows or cols:
  // eg. A = [[ 1,  2,  3],
  //          [ 4,  5,  6],
  //          [ 7,  8,  9],
  //          [10, 11, 12]]
  // A.Sum(kByRow) ==> [6, 15, 24, 33]
  // A.Sum(kByCol) ==> [22, 26, 30]
  //
  Vector<T> Sum(Dim dimention) const;
  Vector<T> Min(Dim dimention) const;
  Vector<T> Max(Dim dimention) const;
  Vector<double> Mean(Dim dimention) const;
  void Norm(Dim dimention=kByRow);            // Only for T==double

#ifdef GSL
  // Get covariance matrix of data seris stored in cols or rows
  Matrix<T> Cov(Dim dimention=kByCol) const;  // Only for T==double
  void Sort(Dim dimention=kByRow, int key=0); // Only for T==double
  void Eigen(Vector<T> *eigen_value, Matrix<T> *eigen_vector); // Only for T==double
#endif // GSL

  void Print() const;
  void Readf(const char *filename);
  void Printf(const char *filename) const;

  bool     empty() const { return (rows_ == 0 or cols_ == 0); }
  isize_t  rows()  const { return rows_;                      }
  isize_t  cols()  const { return cols_;                      }
  isize_t  tda()   const { return tda_;                       }
  T*       data()        { return data_;                      }
  const T* data()  const { return data_;                      }
  void     clear();

 protected:
  bool is_view_() const { return (data_ != NULL and owner_ == false); }
  void copy_(const Matrix<T> &source);
  void allocate_(isize_t rows, isize_t cols);

  isize_t rows_;
  isize_t cols_;
  isize_t tda_;
  T *data_;
  bool owner_;
};

typedef Matrix<int> MatrixI;
typedef Matrix<double> MatrixD;

#include "matrix_inl.h"

} // namespace imat
#endif // IMAT_MATRIX_H
