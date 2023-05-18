#include "general.h"

#include <iostream>
#include <sstream>

namespace apm {

VectorSparse::VectorSparse(vInt indices, vDouble values) {
  if (indices.size() != values.size()) {
    std::cout << "Error the length of indices and values does not match!"
              << std::endl;
  } else {
    this->index = indices;
    this->value = values;
    this->len = indices.size();
  }
}

VectorSparse::VectorSparse(vInt indices, double value) {
  this->len = indices.size();
  std::vector<double> ones(this->len);
  for (size_t i = 0; i < this->len; i++) {
    ones.at(i) = 1;
  }

  this->index = indices;
  this->value = ones;
}

void VectorSparse::print() {
  for (size_t i = 0; i < this->len; i++) {
    std::cout << '(' << this->index.at(i)
              << "," << this->value.at(i)
              << "), ";
  }
  std::cout << std::endl;
}

VectorSparse &VectorSparse::multiply(double num) {
  for (size_t i = 0; i < this->len; i++) {
    this->value.at(i) *= num;
  }
  return *this;
}

VectorSparse &VectorSparse::divide(double num) {
  for (size_t i = 0; i < this->len; i++) {
    this->value.at(i) /= num;
  }
  return *this;
}

double VectorSparse::get(int i) {
  double x = 0;
  size_t j = 0;

  while ((j < this->len) && (this->index.at(j) != i)) {
    ++j;
  }
  if (j < this->len) {
    x = this->value.at(j);
  }

  return x;
}

double VectorSparse::dot(vDouble v) {
  double res = 0;
  for (size_t i = 0; i < this->len; ++i) {
    res += this->value.at(i) * v.at(this->index.at(i));
  }
  return res;
}

vString split(const std::string &str, char delimitor) {
  std::istringstream iss(str);
  std::string item;
  vString result;
  while (std::getline(iss, item, delimitor)) {
    result.push_back(item);
  }

  return result;
}

}  // namespace apm