#include "imat/matrix.h"

using namespace std;
using namespace imat;

void Show(const char *string) {
  cout << "==============================================================================="<< endl;
  cout << string << " ..." << endl;
}

int main (void) {
  // Example: imat::Matrix
  Show("Example: Matrix(isize_t, isize_t)");
  Matrix<double> X(3, 3);
  Show("Example: Matrix.Print()");
  MSSHOW("After X(3, 3):");
  X.Print();
  Show("Example: Matrix.SetZero()");
  MSSHOW("X.SetZero()");
  X.SetZero();
  Show("Example: Matrix.SetUnit()");
  X.SetUnit();
  X.Print();
  Show("Example: operator+= ");
  MSSHOW("X += 1; X += X");
  X += 1;
  X += X;
  X.Print();
  Show("Example: operator*= ");
  MSSHOW("X *= 2; X *= X");
  X *= 2;
  X *= X;
  X.Print();
  Show("Example: Matrix(const Matrix &)");
  Matrix<double> Y(X);
  MSSHOW("Y(X)");
  Y.Print();
  Show("Example: Matrix() and operator=");
  Matrix<double> Z;
  MSSHOW("Z=Y");
  Z = Y;
  Z.Print();
  MatrixD ZZ;
  MSSHOW("ZZ.SubCopy(Y)");
  ZZ.Allocate(3, 3);
  ZZ.SubCopy(Y);
  ZZ.Print();
  //View(), Min(), Max(), Sum(), Norm(), Mean()
  Show("Example: Matrix.View(const Type *, isize_t, isize_t, isize_t tda)");
  double data[12] = {1.1, 2, 3, 4.0, 5, 6, 7, 8, 9, 10, 11, 12};
  Matrix<double> A;
  A.View(data, 3, 4);
  MSSHOW("Matrix view: A");
  A.Print();
  Show("Example: Matrix.Min(kByRow), Matrix.Min(kByCol), Matrix.Min()");
  Vector<double> B;
  MSSHOW("A.Min(kByRow)");
  B = A.Min(imat::kByRow);
  B.Print();
  MSSHOW("A.Min(kByCol)");
  B = A.Min(imat::kByCol);
  B.Print();
  double Min = A.Min();
  cout << "A.Min(): " << Min << endl;
  Show("Example: Matrix.Max(kByRow), Matrix.Max(kByCol), Matrix.Max()");
  MSSHOW("A.Max(kByRow)");
  B = A.Max(kByRow);
  B.Print();
  MSSHOW("A.Max(kByCol)");
  B = A.Max(kByCol);
  B.Print();
  double Max = A.Max();
  cout << "A.Max(): " << Max << endl;
  Show("Example: Matrix.Sum(kByRow), Matrix.Sum(kByCol), Matrix.Sum()");
  MSSHOW("A.Sum(kByRow)");
  B = A.Sum(kByRow);
  B.Print();
  MSSHOW("A.Sum(kByCol)");
  B = A.Sum(kByCol);
  B.Print();
  double Sum = B.Sum();
  cout << "A.Sum(): " << Sum << endl;
  Show("Example: Matrix.Norm(kByRow), Matrix.Norm(kByCol)");
  MSSHOW("A.Norm(kByRow)");
  A.Norm(kByRow);
  A.Print();
  MSSHOW("A.Norm(kByCol)");
  A.Norm(kByCol);
  A.Print();
  Show("Example: Matrix.Mean(kByRow), Matrix.Mean(kByCol), Matrix.Mean()");
  MSSHOW("A.Mean(kByRow)");
  B = A.Mean(kByRow);
  B.Print();
  MSSHOW("A.Mean(kByCol)");
  B = A.Mean(kByCol);
  B.Print();
  double Mean = A.Mean();
  cout << "A.Mean(): " << Mean << endl;
  //File I/O
  Show("Example: file I/O");
  double a[] = {0.11, 0.12, 0.13,
                0.21, 0.22, 0.23 };
  Matrix<double> F, E;
  F.View(a, 2, 3);
  MSSHOW("Matrix to be printed (to matrix_demo.dat):");
  F.Print();
  F.Printf("matrix_demo.dat");
  E.Readf("matrix_demo.dat");
  MSSHOW("Matrix readed (from matrix_demo.dat):");
  E.Print();
  MSSHOW("Finished Matrix.printf() and Matrix.readf().");
  // Cross product
  Show("Example: Cross(Matrix<T> *res, const Matrix<T> &m1, const Matrix<T> &m2)");
  F = E;
  F.Print();
  F.Trans();
  Matrix<double> D(3, 3);
  D.SetZero();
  MSSHOW("A: ");
  F.Print();
  MSSHOW("B: ");
  E.Print();
  Cross(&D, F, E);
  MSSHOW("D = AxB :");
  D.Print();
  //Tanspose and multiply
  Show("Example: Transpose and multiply");
  Matrix<double> G;
  Matrix<double> H(E);
  E.Trans();
  G = H * E;
  G.Print();
  MSSHOW("Finished Tanspose and multiply.");
  //SubMatrix
  Show("Example: SubMatrix");
  MSSHOW("Original matrix (D):");
  D.Print();
  MSSHOW("Submatrix(D(1:3, 1:3)):");
  Matrix<double> sub_matrix;
  sub_matrix.View(&D, 1, 1, 2, 2);
  sub_matrix.Print();
  sub_matrix.clear();
  sub_matrix.View(&D, 1, 1);
  sub_matrix.Print();
  double aa[] = { 1.11, 1.12, 
                  1.21, 1.22};
  Matrix<double> CC;
  CC.View(aa, 2, 2);
  MSSHOW("Matrix be substracted by submatrix(CC):");
  CC.Print();
  sub_matrix -= CC;
  MSSHOW("D(1:3, 1:3) -= CC:");
  sub_matrix.Print();
  MSSHOW("Original matrix (D):");
  D.Print();
  MSSHOW("Finished SubMatrix.");

  Show("Example: SubCopy(src, row_start, col_start, rows, cols)");
  MSSHOW("Original matrix (D):");
  D.Print();
  Matrix<double> sub_copy(2, 2);
  sub_copy.SubCopy(D, 1, 1, 2, 2);
  sub_copy.Print();
  sub_copy.SubCopy(D, 1, 1);
  sub_copy.Print();
  MSSHOW("Matrix be substracted by submatrix(CC):");
  CC.Print();
  MSSHOW("sub_copy -= CC:");
  sub_copy -= CC;
  sub_copy.Print();
  MSSHOW("Original matrix (D):");
  D.Print();
#ifdef GSL
  Show("Example: Eigen(eigen_value, eigen_vector)");
  double eig_data[9] = {14, 10, 12,
                        10,  9, 11, 
                        12, 11, 19};
  MatrixD my_data(3,3);
  Matrix<double> to_eig;
  to_eig.View(eig_data, 3, 3);
  MSSHOW("Matrix to be used for eigen problem");
  to_eig.Print();
  my_data.SubCopy(to_eig, 0, 0, 2, 3);
  Vector<double> eval(3);
  MatrixD evec(3, 3);
  MSSHOW("Matrix to be diagnolized:");
  to_eig.Print();
  cout << endl;
  to_eig.Eigen(&eval, &evec);
  MSSHOW("Eigen values:");
  eval.Print();
  cout << endl;
  MSSHOW("Eigen vectors (in cols):");
  evec.Print();
  cout << endl;
  MSSHOW("Matrix after diagnolization:");
  to_eig.Print();
  cout << endl;
  Show("Example: Covariance(kByCol)");
  my_data.Print();
  MatrixD my_cov;
  my_cov = my_data.Cov();
  my_cov.Print();
  Show("Example: Matrix.Sort(kByRow), Matrix.Sort(kByCol), Matrix.Sort()");
  A.SubCopy(to_eig);
  A.Sort(kByRow, 0);
  MSSHOW("A.Sort(kByRow, 0): ");
  A.Print();
  A.Sort(kByCol, 0);
  MSSHOW("A.Sort(kByCol, 0)");
  A.Print();
#endif // GSL
  Show("Example: Matrx.GetRow(vec_view, index), Matrx.GetCol(vew_view, index)");
  MSSHOW("A:");
  A.Print();
  Vector<double> vec_view;
  MSSHOW("A.GetRow(&vec_view, 0)");
  A.GetRow(&vec_view, 0);
  vec_view.Print();
  MSSHOW("A.GetCol(&vec_view, 1)");
  vec_view.clear();
  A.GetCol(&vec_view, 1);
  vec_view.Print();

}
