#ifndef _GENERAL_H_
#define _GENERAL_H_

#include <set>
#include <string>
#include <vector>

// #define LOG_DEBUG
// #define LOG_VERBOSE

namespace apm {

typedef std::vector<int> vInt;
typedef std::vector<double> vDouble;
typedef std::vector<std::vector<int>> vvInt;
typedef std::vector<std::vector<double>> vvDouble;
typedef std::set<int> sInt;
typedef std::vector<std::string> vString;

typedef vvDouble matrixDense;

/* A sparse vector with basic vector operations.

**Description**
  The vector is stored as two lists. One list of indices and a list of values.
  Allow operations like scalar-multiplication and getting and setting values.
*/
class VectorSparse {
 private:
  vInt index;
  size_t len;

 public:
  vDouble value;

  /* Create a sparse vector.

  **Parameter**
       - indices: vector of indices
       - values: vector of values at the given indices
  */
  VectorSparse(vInt indices, vDouble values);

  /* Create a sparse vector with constant values.

 **Parameter**
      - indices: vector of indices
      - values: value for all entries at the indices
  */
  VectorSparse(vInt indices, double value);

  size_t length() {
    return this->len;
  }

  int getIndex(size_t i) {
    return this->index.at(i);
  }

  double getValue(size_t i) {
    return this->value.at(i);
  }
  /* Get the value at index i.

  **Parameter**
      - i: index of the value to retrieve

  **Return**
      If the index is in the list of vector then return the value. Otherwise
      return 0.

  **Description**
      Loop through the index list to check if the there is a non zero value to
      return. If not match is found return zero.
  */
  double get(int i);

  /* Print the sparse vector.

  **Description**
      Print the non-zero entries as index-value pairs.
  */
  void print();

  /* Scalar multiply the vector with a number.

  **Parameter**
      - num: value to multiply with

  **Return**
      Reference to `this`.

  **Description**
      Multiply all values with `num`.
  */
  VectorSparse &multiply(double num);

  /* Scalar divide the vector by a number.

  **Parameter**
      - num: value to divide by

  **Return**
      Reference to `this`.

  **Description**
      Divide all values by `num`.
  */
  VectorSparse &divide(double num);

  /* Compute the dot-product with the vector. */
  double dot(vDouble v);
};

// A list of sparse vectors (rows of the matrix)
typedef std::vector<VectorSparse> matrixSparse;

// Split the string into substrings delimite dy delim.
vString split(const std::string &str, char delimitor);

}  // namespace apm

#endif  // _GENERAL_H_