#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() noexcept {
  rows_ = 3;
  cols_ = 3;
  matrix_ = new double[cols_ * rows_]{};
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) (*this)(i, j) = i * 3 + j + 1;
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows < 1 || cols < 1) throw std::length_error("Bad size");
  rows_ = rows;
  cols_ = cols;
  matrix_ = new double[cols * rows]{};
}

S21Matrix::S21Matrix(const S21Matrix& other) noexcept {
  cols_ = other.cols_;
  rows_ = other.rows_;
  matrix_ = new double[rows_ * cols_]{};

  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) (*this)(i, j) = other(i, j);
}

S21Matrix::~S21Matrix() noexcept {
  if (matrix_) delete[] matrix_;
  cols_ = 0;
  rows_ = 0;
  matrix_ = nullptr;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  bool res = true;
  if (this->cols_ != other.cols_ && this->rows_ != other.rows_)
    res = false;
  else
    for (int i = 0; i < rows_ && res; i++)
      for (int j = 0; j < cols_ && res; j++)
        if (std::abs((*this)(i, j) - other(i, j) > eps)) res = false;

  return res;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (cols_ != other.cols_ || rows_ != other.rows_)
    throw std::logic_error("Matrices shouls heve the same size");
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) (*this)(i, j) += other(i, j);
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (cols_ != other.cols_ || rows_ != other.rows_)
    throw std::logic_error("Matrices shouls heve the same size");
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) (*this)(i, j) -= other(i, j);
}
void S21Matrix::MulNumber(const double num) noexcept {
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) (*this)(i, j) *= num;
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_)
    throw std::logic_error("Matrices shouls heve the same size");
  S21Matrix res(rows_, other.cols_);
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < other.cols_; j++)
      for (int k = 0; k < cols_; k++) res(i, j) += (*this)(i, k) * other(k, j);
  *this = res;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix res(cols_, rows_);
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) res(j, i) = (*this)(i, j);
  return res;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_)
    throw std::logic_error("Matrices shouls heve the same size");
  S21Matrix calc(rows_, rows_);
  if (rows_ == 1)
    calc(0, 0) = (*this)(0, 0);
  else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < rows_; j++) {
        calc(i, j) = (GetMinor(i, j)).Determinant() * pow(-1, i + j);
      }
    }
  }
  return calc;
}

double S21Matrix::Determinant() {
  double det = 0;
  if (cols_ != rows_)
    throw std::logic_error("Matrices shouls heve the same size");

  if (cols_ == 1)
    det = (*this)(0, 0);
  else {
    int sign = 1;
    for (int i = 0; i < cols_; i++) {
      det += (*this)(i, 0) * sign * (this->GetMinor(i, 0)).Determinant();
      sign *= -1;
    }
  }
  return det;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = this->Determinant();
  if (det == 0) throw std::logic_error("Det = 0");

  return *this;
}
S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(other);
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (cols_ != other.cols_ || rows_ != other.rows_) {
    if (matrix_) delete[] matrix_;
    matrix_ = new double[other.rows_ * other.cols_];
    rows_ = other.rows_;
    cols_ = other.cols_;
  }
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++) (*this)(i, j) = other(i, j);
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
  MulNumber(other);
  return *this;
}

double& S21Matrix::operator()(int i, int j) {
  if (i < 0 || j < 0 || i >= cols_ || j >= rows_)
    throw std::out_of_range("Out of range");
  return matrix_[i * cols_ + j];
}
double& S21Matrix::operator()(int i, int j) const {
  if (i < 0 || j < 0 || i >= cols_ || j >= rows_)
    throw std::out_of_range("Out of range");
  return matrix_[i * cols_ + j];
}

S21Matrix S21Matrix::GetMinor(int row, int col) {
  S21Matrix minor(rows_ - 1, cols_ - 1);
  int i = 0, j = 0;
  for (int ki = 0; ki < rows_; ki++)
    for (int kj = 0; kj < cols_; kj++) {
      if (ki != row && col != kj) {
        minor(i, j) = (*this)(ki, kj);
        j++;
        if (j == cols_ - 1) {
          j = 0;
          i++;
        }
      }
    }
  return minor;
}
void S21Matrix::print() {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      std::cout << matrix_[i * cols_ + j] << ' ';
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

int main() {
  S21Matrix a;
  S21Matrix b;
  a -= b;
  a.print();
}
