#ifndef IMAT_MATRIX_INL_H
#define IMAT_MATRIX_INL_H

#define MFOR_ALL(i, j, c) MFOR((i), (c).rows_) MFOR ((j), (c).cols_)

// Friend functions
template <class T> 
Matrix<T> operator- (const Matrix<T> &lfs, const Matrix<T> &rhs) {
  assert(lfs.rows() == rhs.rows());
  assert(lfs.cols() == rhs.cols());
  Matrix<T> result(lfs.rows(), lfs.cols());
  MFOR_ALL(i, j, lfs) result(i, j) = lfs(i, j) - rhs(i, j);
  return result;
}

template <class T> 
Matrix<T> operator+ (const Matrix<T> &lfs, const Matrix<T> &rhs) {
  assert(lfs.rows() == rhs.rows());
  assert(lfs.cols() == rhs.cols());
  Matrix<T> result(lfs.rows(), lfs.cols());
  MFOR_ALL(i, j, lfs) result(i, j) = lfs(i, j) + rhs(i, j);
  return result;
}

template <class T> 
Matrix<T> operator* (const Matrix<T> &lfs, const Matrix<T> &rhs) {
  assert(lfs.cols() == rhs.rows());
  Matrix<T> result(lfs.rows(), rhs.cols());
  T accumlator;
  MFOR(i, lfs.rows()) {
    MFOR(j, rhs.cols()) {
      accumlator = 0;
      MFOR(k, lfs.cols()) {
        accumlator += lfs(i, k) * rhs(k, j);
      }
      result(i, j) = accumlator;
    }
  }
  return result;
}

template <class T> 
Matrix<T> operator* (const Matrix<T> &mat, const T &x) {
  Matrix<T> result(mat.rows(), mat.cols());
  MFOR_ALL(i, j, mat) result(i, j) = mat(i, j) * x;
  return result;
}

template <class T> 
Matrix<T> operator* (const T &x, const Matrix<T> &mat) {
  Matrix<T> result(mat.rows(), mat.cols());
  MFOR_ALL(i, j, mat) result(i, j) = mat(i, j) * x;
  return result;
}

template <class T> 
void Cross(Matrix<T> *result, const Matrix<T> &lfs, const Matrix<T> &rhs) {
  assert(lfs.cols() == rhs.rows());
  assert(result->rows() >= lfs.rows());
  assert(result->cols() >= rhs.cols());
  T accumlator;
  MFOR(i, lfs.rows()) {
    MFOR(j, rhs.cols()) {
      accumlator = 0;
      MFOR(k, lfs.cols()) { accumlator += lfs(i, k) * rhs(k, j); }
      (*result)(i, j) = accumlator;
    }
  }
}

// Matrix
template <class T>
Matrix<T>::Matrix(isize_t rows, isize_t cols) {
  allocate_(rows, cols);
}

template <class T>
Matrix<T>::Matrix(const Matrix<T> &source) {
  copy_(source);
}

template <class T>
Matrix<T> & Matrix<T>::operator= (const Matrix<T> &source) {
  if (this == &source) return *this;
  if (owner_) delete []data_;
  copy_(source);
  return *this;
}

template <class T>
ListInitializer<T> Matrix<T>::operator= (const T &value) {
  assert(owner_); // Matrix
  return ListInitializer<T>(data_, rows_*cols_, value);
}

template <class T>
void Matrix<T>::Allocate(isize_t rows, isize_t cols) {
  assert(empty());
  allocate_(rows, cols);
}

template <class T>
inline void Matrix<T>::View(T *src, isize_t rows, isize_t cols, isize_t tda) {
  assert(empty());
  if (tda == kDefault) tda = cols;
  data_ = src;
  rows_ = rows;
  cols_ = cols;
  tda_ = tda;
  owner_ = false;
}

template <class T>
inline void Matrix<T>::View(Matrix<T> *src,
                            isize_t row_0, isize_t col_0,
                            isize_t rows, isize_t cols) {
  assert(empty());
  assert(! src->empty());
  if (rows == kDefault and cols == kDefault) {
    rows = src->rows_ - row_0;
    cols = src->cols_ - col_0;
  } else {
    isize_t src_end = (src->rows_ - 1)*src->tda_ + (src->cols_ - 1);
    isize_t view_end = row_0*src->tda_ + col_0 + (rows - 1)*tda_ + (cols -1);
    assert(src_end >= view_end);
  }
  data_ = src->data_ + row_0 * src->tda_ + col_0;
  rows_ = rows;
  cols_ = cols;
  tda_ = src->tda_;
  owner_ = false;
}

template <class T>
inline void Matrix<T>::View(const Matrix<T> &src,
                           isize_t row_0, isize_t col_0,
                           isize_t rows, isize_t cols) {
  assert(empty());
  assert(! src.empty());
  if (rows == kDefault and cols == kDefault) {
    rows = src.rows_ - row_0;
    cols = src.cols_ - col_0;
  } else {
    isize_t src_end = (src.rows_ - 1)*src.tda_ + (src.cols_ - 1);
    isize_t view_end = row_0*src.tda_ + col_0 + (rows - 1)*tda_ + (cols -1);
    assert(src_end >= view_end);
  }
  data_ = src.data_ + row_0 * src.tda_ + col_0;
  rows_ = rows;
  cols_ = cols;
  tda_ = src.tda_;
  owner_ = false;
}

template <class T>
void Matrix<T>::SubCopy(const Matrix<T> &src,
                        isize_t row_0, isize_t col_0,
                        isize_t rows, isize_t cols) {
  assert(! src.empty());
  if (rows == kDefault and cols == kDefault) {
    rows = src.rows_ - row_0;
    cols = src.cols_ - col_0;
  } else {
    isize_t src_end = (src.rows_ - 1)*src.tda_ + (src.cols_ - 1);
    isize_t view_end = row_0*src.tda_ + col_0 + (rows - 1)*tda_ + (cols -1);
    assert(src_end >= view_end);
  }
  assert(rows_ >= rows and cols_ >= cols);

  MFOR(i, rows)
    MFOR(j, cols)
      (*this)(i, j) = src(row_0 + i, col_0 +j);
}

// 
template <class T>
inline T & Matrix<T>::operator() (isize_t i, isize_t j) {
  assert(i>=0 and i<rows_);
  assert(j>=0 and j<cols_);
  return data_[i*tda_ + j];
}

template <class T>
inline const T & Matrix<T>::operator() (isize_t i, isize_t j) const {
  assert(i>=0 and i<rows_);
  assert(j>=0 and j<cols_);
  return data_[i*tda_ + j];
}

template <class T>
Matrix<T> Matrix<T>::operator-() const {
  Matrix<T> result(rows_, cols_);
  MFOR_ALL(i, j, *this) result(i, j) = -(*this)(i, j);
  return result;
}

template <class T>
bool Matrix<T>::operator== (const Matrix<T> &other) const {
  bool is_equal = true;
  if (rows_ != other.rows_ or cols_ != other.cols_) {
    is_equal = false;
  } else {
    MFOR(i, rows_) {
      MFOR(j, cols_) {
        if ((*this)(i, j) != other(i, j)) {
          is_equal = false;
          break;
        }
      }
      if (!is_equal) break;
    }
  }
  return is_equal;
}

template <class T>
void Matrix<T>::operator+= (const Matrix<T> &other) {
  assert(rows_ == other.rows_);
  assert(cols_ == other.cols_);
  MFOR_ALL(i, j, *this) (*this)(i, j) += other(i, j);
}

template <class T>
void Matrix<T>::operator-= (const Matrix<T> &other) {
  assert(rows_ == other.rows_);
  assert(cols_ == other.cols_);
  MFOR_ALL(i, j, *this) (*this)(i, j) -= other(i, j);
}

template <class T>
void Matrix<T>::operator*= (const Matrix<T> &other) {
  assert(rows_ == other.rows_);
  assert(cols_ == other.cols_);
  MFOR_ALL(i, j, *this) (*this)(i, j) *= other(i, j);
}

template <class T>
void Matrix<T>::operator+= (const T &x) {
  MFOR_ALL(i, j, *this) (*this)(i, j) += x;
}

template <class T>
void Matrix<T>::operator-= (const T &x) {
  MFOR_ALL(i, j, *this) (*this)(i, j) -= x;
}

template <class T>
void Matrix<T>::operator*= (const T &x) {
  MFOR_ALL(i, j, *this) (*this)(i, j) *= x;
}

// 
template <class T>
void Matrix<T>::SetZero() {
  MFOR_ALL(i, j, *this) (*this)(i, j) = 0;
}

template <class T>
void Matrix<T>::SetAll(const T &x) {
  MFOR_ALL(i, j, *this) (*this)(i, j) = x;
}

template <class T>
void Matrix<T>::SetUnit() {
  MFOR(i, rows_) {
    (*this)(i, i) = 1;
    MRANGE (j, i+1, cols_) { 
      (*this)(i, j) = 0; (*this)(j, i) = 0;
    }
  }
}

template <class T>
T Matrix<T>::Sum() const {
  T accumlator = 0;
  MFOR_ALL(i, j, *this) accumlator += (*this)(i, j);
  return accumlator;
}

template <class T>
T Matrix<T>::Min() const {
  T min= (*this)(0, 0);
  MFOR_ALL(i, j, *this) min = min < (*this)(i, j) ? min : (*this)(i, j);
  return min;
}

template <class T>
T Matrix<T>::Max() const {
  T max= (*this)(0, 0);
  MFOR_ALL(i, j, *this) max = max > (*this)(i, j) ? max : (*this)(i, j);
  return max;
}

template <class T>
double Matrix<T>::Mean() const {
  return ((double) Sum()) / (rows_ * cols_);
}

template <class T>
inline void Matrix<T>::GetRow(Vector<T> *row_view, isize_t i) {
  row_view->View(data() + i*tda(), cols(), 1);
}

template <class T>
inline void Matrix<T>::GetCol(Vector<T> *col_view, isize_t i) {
  col_view->View(data() + i, rows(), tda());
}


inline void ChooseDim(isize_t *num_base, isize_t *num_cycles,
                      isize_t *step_base, isize_t *step_cycles,
                      isize_t rows, isize_t cols, isize_t tda,
                      Dim dimention) {
  if (dimention == kByRow) {
    *num_base = rows;
    *num_cycles = cols;
    *step_base = tda;
    *step_cycles = 1;
  } else if (dimention == kByCol) {
    *num_base = cols;
    *num_cycles = rows;
    *step_base = 1;
    *step_cycles = tda;
  }
}

template <class T>
Vector<T> Matrix<T>::Sum(Dim dimention) const {
  assert(dimention == 0 or dimention == 1);
  isize_t num_base, num_cycles, step_base, step_cycles;
  ChooseDim(&num_base, &num_cycles, &step_base, &step_cycles, 
            rows_, cols_, tda_,
            dimention);
  Vector<T> result(num_base);
  result.SetZero();
  MFOR(i, num_base) {
    MFOR(j, num_cycles) {
      result(i) += *(data_ + i*step_base + j*step_cycles);
    }
  }
  return result;
}

template <class T>
Vector<T> Matrix<T>::Min(Dim dimention) const {
  assert(dimention == 0 or dimention == 1);
  isize_t num_base, num_cycles, step_base, step_cycles;
  ChooseDim(&num_base, &num_cycles, &step_base, &step_cycles, 
            rows_, cols_, tda_,
            dimention);
  Vector<T> result(num_base);
  result.SetZero();
  MFOR(i, num_base) {result(i) = *(data_ + i * step_base);}
  MFOR(i, num_base) {
    MRANGE (j, 1, num_cycles) {
      result(i) = result(i) < *(data_ + i*step_base + j*step_cycles) ?
                      result(i) : *(data_ + i*step_base + j*step_cycles);
    }
  }
  return result;
}

template <class T>
Vector<T> Matrix<T>::Max(Dim dimention) const {
  assert(dimention == 0 or dimention == 1);
  isize_t num_base, num_cycles, step_base, step_cycles;
  ChooseDim(&num_base, &num_cycles, &step_base, &step_cycles, 
            rows_, cols_, tda_,
            dimention);
  Vector<T> result(num_base);
  result.SetZero();
  MFOR(i, num_base) {result(i) = *(data_ + i * step_base);}
  MFOR(i, num_base) {
    MRANGE (j, 1, num_cycles) {
      result(i) = result(i) > *(data_ + i*step_base + j*step_cycles) ?
                      result(i) : *(data_ + i*step_base + j*step_cycles);
    }
  }
  return result;
}

template <class T>
Vector<double> Matrix<T>::Mean(Dim dimention) const {
  assert(dimention == 0 or dimention == 1);
  isize_t num_base, num_cycles, step_base, step_cycles;
  ChooseDim(&num_base, &num_cycles, &step_base, &step_cycles, 
            rows_, cols_, tda_,
            dimention);
  Vector<double> result(num_base);
  result.SetZero();
  MFOR(i, num_base) {
    MFOR(j, num_cycles) {
      result(i) += *(data_ + i*step_base + j*step_cycles);
    }
  }
  MFOR(i, num_base) {
      result(i) /= num_cycles;
  }
  return result;
}

template <class T>
void Matrix<T>::Norm(Dim dimention) {
  MESHOW("Norm() is not defined for Matrix other than MatrixD.");
}

#ifdef GSL
template <class T>
Matrix<T> Matrix<T>::Cov(Dim dimention) const {
  MESHOW("Covariance() is not defined for Matrix other than MatrixD.");
  return Matrix<T>();
}

template <class T>
void Matrix<T>::Sort(Dim dimention, int key) {
  MESHOW("Sort() is not defined for Matrix other than MatrixD.");
}

template <class T>
void Matrix<T>::Eigen(Vector<T> *eigen_value, Matrix<T> *eigen_vector) {
  MESHOW("Eigen() is not defined for Matrix other than MatrixD.");
}
#endif // GSL

template <>
inline void Matrix<double>::Norm(Dim dimention) {
  assert(dimention == 0 or dimention == 1);
  isize_t num_base, num_cycles, step_base, step_cycles;
  ChooseDim(&num_base, &num_cycles, &step_base, &step_cycles, 
            rows_, cols_, tda_,
            dimention);
  Vector<double> summed(num_base);
  summed.SetZero();
  MFOR(i, num_base) {
    MFOR(j, num_cycles) {
      double element = *(data_ + i*step_base + j*step_cycles);
      summed(i) += element * element;
    }
  }
  MFOR(i, num_base) {summed(i) = sqrt(summed(i));}
  MFOR(i, num_base) {
    MFOR(j, num_cycles) {
      *(data_ + i*step_base + j*step_cycles) /= summed(i);
    }
  }
}

#ifdef GSL
template <>
inline Matrix<double> Matrix<double>::Cov(Dim dimention) const {
  assert(dimention == 0 or dimention == 1);
  isize_t num_base, num_cycles, step_base, step_cycles;
  ChooseDim(&num_base, &num_cycles, &step_base, &step_cycles, 
            rows_, cols_, tda_,
            dimention);
  MatrixD cov_mat(num_base, num_base);
  MFOR(i, num_base) {
    MRANGE (j, i, num_base) {
      double *data1, *data2;
      data1 = data_ + i*step_base;
      data2 = data_ + j*step_base;
      cov_mat(i, j) = gsl_stats_covariance(data1, step_cycles,
                                           data2, step_cycles,
                                           num_cycles);
      if ( i != j) cov_mat(j, i) = cov_mat(i, j);
    }
  }
  return cov_mat;
}

template <>
inline void Matrix<double>::Sort(Dim dimention, int key) {
  assert(dimention == 0 or dimention == 1);
  isize_t num_base, num_cycles, step_base, step_cycles;
  ChooseDim(&num_base, &num_cycles, &step_base, &step_cycles, 
            rows_, cols_, tda_,
            dimention);
  size_t *p = new size_t[num_base];
  gsl_sort_index(p, data_ + key*step_cycles, step_base, num_base);
  Matrix<double> swap(*this);
  isize_t s_num_base, s_num_cycles, s_step_base, s_step_cycles;
  ChooseDim(&s_num_base, &s_num_cycles, &s_step_base, &s_step_cycles, 
            swap.rows_, swap.cols_, swap.tda_,
            dimention);
  MFOR(i, num_base) {
    MFOR(j, num_cycles) {
      *(data_ + i*step_base + j*step_cycles) = 
            *(swap.data_ + p[i]*s_step_base + j*s_step_cycles);
    }
  }
  delete []p;
}

template <>
inline void Matrix<double>::Eigen(Vector<double> *eigen_value, 
                           Matrix<double> *eigen_vector) {
  assert(rows_ == cols_);
  assert(eigen_value->size() == rows_);
  assert(eigen_vector->rows() == rows_);
  assert(eigen_vector->cols() == cols_);
  /*
  MSSHOW("==============================================================");
  MSSHOW(" Solving the eigen problem by MatrixD.Eigen() ...");
  MSSHOW("");
  MSSHOW(" BEWARE!");
  MSSHOW(" 1, Eigen() is only defined for Real Symmetric Matrices; ");
  MSSHOW(" 2, The lower left pair of original matrix will be destroyed");
  MSSHOW("    during the solving process.");
  MSSHOW("==============================================================");
  MSSHOW("");
  */
  gsl_matrix_view m = gsl_matrix_view_array_with_tda(data_,
                                                     rows_,
                                                     cols_,
                                                     tda_);
  gsl_matrix_view evec = gsl_matrix_view_array_with_tda(eigen_vector->data(),
                                                        eigen_vector->rows(),
                                                        eigen_vector->cols(),
                                                        eigen_vector->tda());
  gsl_vector_view eval = gsl_vector_view_array_with_stride(
      eigen_value->data(), 
      eigen_value->stride(),
      eigen_value->size());
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(rows_);
  gsl_eigen_symmv(&m.matrix, &eval.vector, &evec.matrix, w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(&eval.vector, &evec.matrix, GSL_EIGEN_SORT_ABS_ASC);
  /*
  MSSHOW("Eigen vectors (in cols):");
  eigen_vector->Print();
  MSSHOW("");
  MSSHOW("Eigen values:");
  eigen_value->Print();
  MSSHOW("");
  MSSHOW("Array be diagonalized:");
  this->Print();
  MSSHOW("");
  */
}
#endif // GSL

template <class T>
void Matrix<T>::Trans() {
  assert(!is_view_());
  Matrix<T> old_data(*this);
  std::swap(rows_, cols_);
  tda_ = cols_;
  MFOR_ALL(i, j, *this) (*this)(i, j) = old_data(j, i);
}

template <class T>
void Matrix<T>::clear() {
  if (owner_) delete []data_;
  data_ = NULL;
  rows_ = cols_ = tda_ = 0;
  owner_ = false;
}

template <class T>
void Matrix<T>::copy_(const Matrix<T> &source) {
  allocate_(source.rows_, source.cols_);
  if (! source.is_view_()) {    // source is Matrix
    memcpy(data_, source.data_, sizeof(T) * rows_ * cols_);
  } else {                      // source is Matrix view
    MFOR_ALL(i, j, source) (*this)(i, j) = source(i, j);
  }
}

template <class T>
inline void Matrix<T>::allocate_(isize_t rows, isize_t cols) {
  assert(rows > 0 and cols > 0);
  rows_ = rows;
  cols_ = cols;
  tda_ = cols;
  data_ = new T [rows_ * cols_];
  assert(data_);
  owner_ = true;
}

// I/O functions
template <class T>
void Matrix<T>::Print() const {
  // std::cout << "rows: " << rows_ << " cols: " << cols_ << std::endl;
  /*
  if (!owner_) {
    MSSHOW("(Not Owner of Data)");
  } else {
    MSSHOW("(Owner of Data)");
  }
  if (data_ != NULL) {
    std::cout << "data_ (address): " << data_ << std::endl;
  } else {
    MSSHOW("data_ is NULL");
  }
  */
  MFOR(i, rows_) {
    if (i==0) { 
      std::cout << "[[";
    } else {
      std::cout << " [";
    }
    MFOR(j, cols_-1) {
      std::cout.width(8);
      std::cout << (*this)(i, j) << ", ";
    }
    std::cout.width(8);
    if (i==rows_-1) {
      std::cout << (*this)(i, cols_-1) << "]]" << std::endl;
    }
    else {
      std::cout << (*this)(i, cols_-1) << "]," << std::endl;
    }
  }
}

template <class T>
void Matrix<T>::Readf(const char *filename) {
  std::ifstream instream(filename);
  // Head line
  std::string headline;
  getline(instream, headline);
  long rows, cols;
  if (headline.size() >=1 && headline[0] == '%') {
    sscanf(headline.c_str(), "%% %ld\t%ld\n", &rows, &cols);
  } else {
    std::cerr << "A file contains a matrix is needed!  Please check the "
              << "formate of input file: " << filename  << std::endl;
  }
  // Data
  if (owner_) delete []data_;
  allocate_(rows, cols);
  isize_t i = 0;
  while (instream >> data_[i]) {++i;}
  instream.close();
}

template <class T>
void Matrix<T>::Printf(const char *filename) const {
  std::ofstream os(filename);
  assert(os);
  os << "% " << rows_ << '\t' << cols_ << std::endl;
  MFOR(i, rows_) {
    MFOR(j, cols_-1) { os << (*this)(i, j) << " "; }
    os << (*this)(i, cols_-1) << std::endl;
  }
  os.close();
}

#endif // IMAT_MATRIX_INL_H
