#include "imat/vector.h"

#include <cstdlib>
#include <cstdio>
#include <ctime>

using namespace imat;
using namespace std;

void Show(const char *string) {
  cout << "==============================================================================="<< endl;
  cout << string << " ..." << endl;
}

int main() {
  int SIZE = 10;
  int source[SIZE];
  for (int i=0; i<SIZE; i++) { source[i] = (i + 1) * 2; }
  Show("Example: Vector(int) and operator()");
  VectorI vec1(SIZE);
  for (int i=0; i<SIZE; i++) vec1(i) = source[i];
  vec1.Print();
  printf("Vector.Sum() is %d; Vector.Product() is %u\n", vec1.Sum(), vec1.Product());
  Show("Example: copy-construction function");
  VectorI vec2(vec1);
  vec2.Print();
  Show("Example: operator=");
  VectorI vec3;
  vec3 = vec1;
  vec3.Print();
  Show("Example: Vector.View(int *source, isize_t size)");
  VectorI vec4;
  vec4.View(source, SIZE);
  vec4.Print();
  Show("Example: operator- :");
  VectorI vec5 = vec3 - vec2;
  vec5.Print();
  Show("Example: operator+=, Vector.min(), Vector.max():");
  vec1.Print();
  MSSHOW("vec1 += 1");
  vec1 += 1;
  vec1.Print();
  vec3.Print();
  MSSHOW("vec1 += vec3");
  vec1 += vec3;
  vec1.Print();
  printf("vec1.Min(): %u, vec1.Max(): %u\n", vec1.Min(), vec1.Max());
  Show("Example: ListInitializer<T> operator=(const T &value)");
  VectorD vec6(SIZE);
  vec6 = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9;
  vec6.Print();
  Show("Example: Vector.Norm()");
  vec6.Print();
  vec6.Norm();
  vec6.Print();
  vec1.Norm();
  Show("Example: Vector.SubCopy(src, start, num):");
  VectorI vec7(5);
  vec1.Print();
  MSSHOW(" vec7.SubCopy(vec1, 0, 5)");
  vec7.SubCopy(vec1, 0, 5);
  vec7.Print();
  Show("Example: VectorI::const_iterator:");
  VectorI::const_iterator cit;
  for (cit = vec7.begin(); cit != vec7.end(); ++cit) {
    if (cit == vec7.begin()) MSSHOW(" Has operator== ");
    cout << *cit << " ";
  }
  MSSHOW("");
  if (cit != vec7.begin()) MSSHOW(" Has operator!= ");
  Show("Example: VectorI::iterator:");
  VectorI::iterator it;
  for (it = vec7.begin(); it != vec7.end(); ++it) {
    if (it == vec7.begin()) MSSHOW(" Has operator== ");
    cout << *it << " ";
  }
  MSSHOW("");
  if (it != vec7.begin()) MSSHOW(" Has operator!= ");
  Show("Example: Vector<T>::SubCopy(const Vector<T>&)");
  VectorI vec8(vec7.size());
  vec8.SubCopy(vec7);
  vec8.Print();
  return EXIT_SUCCESS;
}
