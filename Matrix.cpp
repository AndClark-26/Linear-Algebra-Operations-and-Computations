#include <iostream>
#include <iomanip>
#include <utility>
#include <cmath>
#include "Matrix.h"

namespace linear {

	//Constants
	const double EPSILON = 1e-6;

	//Matrix constructors and destructors
	Matrix::Matrix(const int i, const int j) {
		m = new double [i * j];

		rows = new int (i);
		columns = new int (j);
	}

	Matrix::Matrix() {
		rows = new int (0);
		columns = new int (0);
	}

	Matrix::Matrix(const Matrix& A) {
		if (m) {
			delete[] m;
			delete rows;
			delete columns;
		}
		
		rows = new int(*A.rows);
		columns = new int(*A.columns);
		
		m = new double[(*rows) * (*columns)];
		for (int i = 0; i < (*rows) * (*columns); i++) {
			*(m + i) = *(A.m + i);
		}
	}

	Matrix::~Matrix() {
		if (m) {
			delete[] m;
		}

		delete rows;
		delete columns;

		std::cout << "Matrix deleted\n";
	}

	//Matrix public member functions
	void Matrix::fill() {
		for (int i = 0; i < *rows; i++) {
			std::cout << "Enter " << *columns << " values\n";
			for (int j = 0; j < *columns; j++) {
				std::cin >> m[i * (*columns) + j];
			}
		}
	}

	void Matrix::fill(const double arr[], const int size) {
		if (size < *rows * *columns) {
			std::cout << "Insufficient elements\n";
		}
		else {
			int counter = 0;
			for (int i = 0; i < *rows; i++) {
				for (int j = 0; j < *columns; j++) {
					m[i * (*columns) + j] = arr[counter];
					counter++;
				}
			}
		}
	}

	void Matrix::print() const {
		std::cout << std::fixed << std::setprecision(1);

		for (int i = 0; i < *rows; i++) {
			for (int j = 0; j < *columns; j++) {
				std::cout << std::setw(6) << std::setfill(' ') << m[i * (*columns) + j];
			}
			std::cout << "\n";
		}

		std::cout << "\n";
	}

	void Matrix::dimensions() const {
		std::cout << *rows << " x " << *columns << " matrix \n";
	}

	void Matrix::transpose() {
		double* A = new double[(*rows) * (*columns)];

		for (int i = 0; i < *rows; i++) {
			for (int j = 0; j < *columns; j++) {
				A[j * (*rows) + i] = m[i * (*columns) + j];
			}
		}

		delete[] m;
		m = A;

		int* temp = new int (*rows);
		delete rows;
		rows = new int (*columns);
		delete columns;
		columns = new int (*temp);
		delete temp;
	}

	void Matrix::rref() {
		int pivot_row = 0;
		int pivot_column = 0;

		while (pivot_column < *columns && pivot_row < *rows) {
			if (nonzero_row(pivot_row, pivot_column)) {
				leading_one(pivot_row, pivot_column);
				zero_above(pivot_row, pivot_column);
				zero_below(pivot_row, pivot_column);
				pivot_row++;
			}
			pivot_column++;
		}
	}

	void Matrix::rref_inverse() {
		int pivot_row = 0;
		int pivot_column = 0;

		while (pivot_column < *columns / 2 && pivot_row < *rows) {
			if (nonzero_row(pivot_row, pivot_column)) {
				leading_one(pivot_row, pivot_column);
				zero_above(pivot_row, pivot_column);
				zero_below(pivot_row, pivot_column);
				pivot_row++;
			}
			pivot_column++;
		}
	}

	//Matrix protected member functions
	void Matrix::swap(const int R1, const int R2) {
		for(int i = 0; i < *columns; i++) {
			std::swap(m[R1 * (*columns) + i], m[R2 * (*columns) + i]);
		}
	}

	void Matrix::multiply(const int R, const double val) {
		for (int i = 0; i < *columns; i++) {
			m[R * (*columns) + i] *= val;
		}
	}

	void Matrix::divide(const int R, const double val) {
		for (int i = 0; i < *columns; i++) {
			m[R * (*columns) + i] /= val;
		}
	}

	void Matrix::add_row(const int R1, const int R2, const double val) {
		for (int i = 0; i < *columns; i++) {
			m[R1 * (*columns) + i] += m[R2 * (*columns) + i] * val;
		}
	}

	void Matrix::sub_row(const int R1, const int R2, const double val) {
		for (int i = 0; i < *columns; i++) {
			m[R1 * (*columns) + i] -= m[R2 * (*columns) + i] * val;
		}
	}

	bool Matrix::nonzero_row(const int& pivot_row, const int& pivot_column) {
		for (int i = pivot_row; i < *rows; i++) {
			if (fabs(m[i * (*columns) + pivot_column]) > EPSILON) {
				if (i == pivot_row) {
					return true;
				}
				else {
					swap(pivot_row, i);
					return true;
				}
			}
		}

		return false;
	}

	void Matrix::leading_one(const int& pivot_row, const int& pivot_column) {
		if (m[pivot_row * (*columns) + pivot_column] != 1) {
			divide(pivot_row, m[pivot_row * (*columns) + pivot_column]);
		}
	}

	void Matrix::zero_below(const int& pivot_row, const int& pivot_column) {
		for (int i = pivot_row + 1; i < *rows; i++) {
			if (fabs(m[i * (*columns) + pivot_column]) > EPSILON) {
				sub_row(i, pivot_row, m[i * (*columns) + pivot_column]);
				m[i * (*columns) + pivot_column] = 0;
			}
		}
	}

	void Matrix::zero_above(const int& pivot_row, const int& pivot_column) {
		for (int i = pivot_row - 1; i >= 0; i--) {
			if (fabs(m[i * (*columns) + pivot_column]) > EPSILON) {
				sub_row(i, pivot_row, m[i * (*columns) + pivot_column]);
				m[i * (*columns) + pivot_column] = 0;
			}
		}
	}
	
	//Matrix operators
	Matrix operator * (const Matrix& A, const Matrix& B) {
		if (*A.columns != *B.rows) {
			std::cout << "Undefined operation";
			return Matrix(0, 0);
		}
		else {
			Matrix C(*A.rows, *B.columns);

			for (int i = 0; i < *A.rows; i++) {
				for (int j = 0; j < *B.columns; j++) {
					C.m[i * (*C.columns) + j] = 0;

					for (int k = 0; k < *A.columns; k++) {
						C.m[i * (*C.columns) + j] = C.m[i * (*C.columns) + j] + A.m[i * (*A.columns) + k] * B.m[k * (*B.columns) + j];
					}
				}
			}

			return C;
		}
	}
	Matrix operator * (Matrix& A, const int scalar) {
		for (int i = 0; i < *A.rows; i++) {
			for (int j = 0; j < *A.columns; j++) {
				A.m[i * (*A.columns) + j] *= scalar;
			}
		}

		return A;
	}

	Matrix operator + (const Matrix& A, const Matrix& B) {
		if (*A.rows != *B.rows && *A.columns != *B.columns) {
			std::cout << "Undefined operation\n";
			return Matrix(0, 0);
		}
		else {
			Matrix C(*A.rows, *A.columns);

			for (int i = 0; i < *C.rows; i++) {
				for (int j = 0; j < *C.columns; j++) {
					C.m[i * (*C.columns) + j] = 0;
					C.m[i * (*C.columns) + j] = A.m[i * (*A.columns) + j] + B.m[i * (*B.columns) + j];
				}
			}

			return C;
		}
	}

	Matrix operator - (const Matrix& A, const Matrix& B) {
		if (*A.rows != *B.rows && *A.columns != *B.columns) {
			std::cout << "Undefined operation\n";
			return Matrix(0, 0);
		}
		else {
			Matrix C(*A.rows, *A.columns);

			for (int i = 0; i < *A.rows; i++) {
				for (int j = 0; j < *A.columns; j++) {
					C.m[i * (*A.columns) + j] = 0;
					C.m[i * (*A.columns) + j] = A.m[i * (*A.columns) + j] - B.m[i * (*B.columns) + j];
				}
			}

			return C;
		}
	}

    void Matrix::operator=(const Matrix &other)
    {
        if (m) {
			delete[] m;
		}
		
		delete rows;
		delete columns;
			
		rows = new int(*other.rows);
		columns = new int(*other.columns);
			
		m = new double[(*rows) * (*columns)];
		for (int i = 0; i < (*rows) * (*columns); i++) {
			*(m + i) = *(other.m + i);
		}
    }

    //Matrix functions
	Matrix rref(const Matrix& A) {
		Matrix B(A);
		B.rref();
		return B;
	}

	Matrix transpose(const Matrix& A) {
		Matrix B(A);
		B.transpose();
		return B;
	}

	double determinant(const Matrix& A) {
		if (*A.rows != *A.columns) {
			std::cout << "Undefined determinant\n";
			return 0;
		}

		if (*A.rows == 2) {
			return A.m[0] * A.m[3] - A.m[1] * A.m[2];
		}

		double sum = 0;

		for (int i = 0; i < *A.columns; i++) {
			sum += A.m[i] * pow(-1, i + 2) * determinant(minor(A, 0, i));
		}

		return sum;
	}

	Matrix minor(const Matrix& A, const int row, const int column) {
		Matrix M(*A.rows - 1, *A.rows - 1);
		int index = 0;

		for (int i = 0; i < *A.rows; i++) {
			if (i == row) {
				continue;
			}
			for (int j = 0; j < *A.columns; j++) {
				if (j == column) {
					continue;
				}
				M.m[index] = A.m[i * (*A.columns) + j];
				index++;
			}
		}

		return M;
	}

	Matrix inverse(const Matrix& A) {
		if (*A.rows != *A.columns) {
			std::cout << "Undefined inverse matrix\n";
			return Matrix(0, 0);
		}
		
		Matrix B(*A.rows, 2 * (*A.columns));

		for (int i = 0; i < *B.rows; i++) {
			for (int j = 0; j < *B.columns; j++) {
				if (i < *A.rows && j < *A.columns) {
					B.m[i * (*B.columns) + j] = A.m[i * (*A.columns) + j];
				}
				else if (j - (*A.columns) == i) {
					B.m[i * (*B.columns) + j] = 1;
				}
				else {
					B.m[i * (*B.columns) + j] = 0;
				}
			}
		}

		B.rref_inverse();
		Matrix INV = submatrix(B, 0, *A.columns);
		return INV;
	}
	Matrix submatrix(const Matrix& A, const int a, const int b) {
		Matrix SUB(*A.rows - a, *A.columns - b);

		int index = 0;
		for (int i = a; i < *A.rows; i++) {
			for (int j = b; j < *A.columns; j++) {
				SUB.m[index] = A.m[i * (*A.columns) + j];
				index++;
			}
		}

		return SUB;
	}

}
