#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>

namespace linear {

	class Matrix {
	public:
		Matrix(const int i, const int j);
		Matrix();
		Matrix(const Matrix& other);
		~Matrix();
		
		void fill();
		void fill(const double arr[], const int size);
		void print() const;
		void dimensions() const;
		void transpose();
		void rref();
		void rref_inverse();

		friend Matrix operator * (const Matrix& A, const Matrix& B);
		friend Matrix operator * (Matrix& A, const int scalar);
		friend Matrix operator + (const Matrix& A, const Matrix& B);
		friend Matrix operator - (const Matrix& A, const Matrix& B);
		friend Matrix rref(const Matrix& A);
		friend Matrix transpose(const Matrix& A);
		friend double determinant(const Matrix& A);
		friend Matrix minor(const Matrix& A, const int row, const int column);
		friend void identity(Matrix& A, const int n);
		friend Matrix inverse(const Matrix& A);
		friend Matrix submatrix(const Matrix& A, const int a, const int b);

		void operator = (const Matrix& other);

	protected:
		double* m;
		int* rows;
		int* columns;
		
		void swap(const int r1, const int r2);
		void multiply(const int r, const double val);
		void divide(const int r, const double val);
		void add_row(const int r1, const int r2, const double val);
		void sub_row(const int r1, const int r2, const double val);

		bool nonzero_row(const int& pivot_row, const int& pivot_column);
		void leading_one(const int& pivot_row, const int& pivot_column);
		void zero_below(const int& pivot_row, const int& pivot_column);
		void zero_above(const int& pivot_row, const int& pivot_column);

	};

}

#endif
