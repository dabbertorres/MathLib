#ifndef SWM_MATRIX_HPP
#define SWM_MATRIX_HPP

#include <vector>
#include <algorithm>

namespace swm
{
	template<typename T, std::size_t Cols, std::size_t Rows>
	class Matrix
	{
		// Matrix types are friends
		friend class Matrix;

		public:
			Matrix();
			Matrix(std::initializer_list<T> il);

			template<std::size_t Size>
			static Matrix<T, Size, Size> identity();

			static Matrix<T, Rows, Cols> transpose(const Matrix<T, Cols, Rows>& m);

			// element accessor
			T& operator()(int x, int y);
			const T& operator()(int x, int y) const;

			constexpr int cols() const;
			constexpr int rows() const;

			template<std::size_t NewCols, std::size_t NewRows>
			Matrix<T, NewCols, NewRows> slice() const;

			// convert element type
			template<typename U>
			operator Matrix<U, Cols, Rows>() const;

			// matrix math
			template<typename U>
			auto operator+ (const Matrix<U, Cols, Rows>& other);

			template<typename U>
			auto operator- (const Matrix<U, Cols, Rows>& other);

			template<typename U, std::size_t OtherCols>
			auto operator* (const Matrix<U, OtherCols, Cols>& other);

			template<typename U, typename = std::enable_if<std::is_fundamental<U>::value>>
			friend auto operator* (const Matrix<T, Cols, Rows>& other, const U& scalar);

			template<typename U, typename = std::enable_if<std::is_fundamental<U>::value>>
			friend auto operator* (const U& scalar, const Matrix<T, Cols, Rows>& other);

		private:
			std::vector<T> data;
	};

	template<typename T, std::size_t Cols, std::size_t Rows>
	Matrix<T, Cols, Rows>::Matrix()
	:	data(Cols * Rows, 0)
	{}

	template<typename T, std::size_t Cols, std::size_t Rows>
	Matrix<T, Cols, Rows>::Matrix(std::initializer_list<T> il)
	:	data(il)
	{
		data.resize(Cols * Rows);
	}

	template<typename T, std::size_t Cols, std::size_t Rows>
	template<std::size_t Size>
	Matrix<T, Size, Size> Matrix<T, Cols, Rows>::identity()
	{
		Matrix<T, Size, Size> id;

		for(int i = 0; i < Size; ++i)
		{
			id(i, i) = 1;
		}

		return id;
	}

	template<typename T, std::size_t Cols, std::size_t Rows>
	Matrix<T, Rows, Cols> Matrix<T, Cols, Rows>::transpose(const Matrix<T, Cols, Rows>& m)
	{
		Matrix<T, Rows, Cols> t;

		for(int i = 0; i < Rows; ++i)
		{
			for(int j = 0; j < Cols; ++j)
			{
				t(i, j) = m(j, i);
			}
		}

		return t;
	}

	template<typename T, std::size_t Cols, std::size_t Rows>
	T& Matrix<T, Cols, Rows>::operator()(int col, int row)
	{
		return data[col + row * Cols];
	}

	template<typename T, std::size_t Cols, std::size_t Rows>
	const T& Matrix<T, Cols, Rows>::operator()(int col, int row) const
	{
		return data[col + row * Cols];
	}

	template<typename T, std::size_t Cols, std::size_t Rows>
	constexpr int Matrix<T, Cols, Rows>::cols() const
	{
		return Cols;
	}

	template<typename T, std::size_t Cols, std::size_t Rows>
	constexpr int Matrix<T, Cols, Rows>::rows() const
	{
		return Rows;
	}

	template<typename T, std::size_t Cols, std::size_t Rows>
	template<std::size_t NewCols, std::size_t NewRows>
	Matrix<T, NewCols, NewRows> Matrix<T, Cols, Rows>::slice() const
	{
		Matrix<T, NewCols, NewRows> sliced;

		constexpr std::size_t lesserRows = Rows < NewRows ? Rows : NewRows;

		// if the number of columns is changing, we can't copy in order
		if(Cols != NewCols)
		{
			constexpr std::size_t lesserCols = Cols < NewCols ? Cols : NewCols;

			const auto slicedBegin = sliced.data.begin();
			const auto thisBegin = data.begin();

			for(int i = 0; i < lesserRows; ++i)
			{
				const auto thisRowBegin = thisBegin + i * Cols;
				const auto thisRowEnd = thisRowBegin + lesserCols;
				const auto slicedRowBegin = slicedBegin + i * NewCols;

				auto thisRow = thisRowBegin;
				auto slicedRow = slicedRowBegin;
				for(; thisRow != thisRowEnd; ++thisRow, ++slicedRow)
				{
					*slicedRow = *thisRow;
				}
			}
		}
		else
		{
			const auto thisBegin = data.begin();
			const auto end = thisBegin + lesserRows * Cols;

			sliced.data.insert(sliced.data.begin(), thisBegin, end);
		}

		return sliced;
	}

	template<typename T, std::size_t Cols, std::size_t Rows>
	template<typename U>
	Matrix<T, Cols, Rows>::operator Matrix<U, Cols, Rows>() const
	{
		Matrix<U, Cols, Rows> newMatrix;

		std::transform(data.begin(), data.end(), newMatrix.data.begin(), [](const T& in)
		{
			return static_cast<U>(in);
		});

		return newMatrix;
	}

	// matrix math
	template<typename T, std::size_t Cols, std::size_t Rows>
	template<typename U>
	auto Matrix<T, Cols, Rows>::operator+ (const Matrix<U, Cols, Rows>& other)
	{
		using Return = decltype(std::declval<T>() + std::declval<U>());

		Matrix<Return, Cols, Rows> newMatrix;

		std::transform(data.begin(), data.end(), other.data.begin(), newMatrix.data.begin(), [](const T& lhs, const U& rhs)
		{
			return lhs + rhs;
		});

		return newMatrix;
	}

	template<typename T, std::size_t Cols, std::size_t Rows>
	template<typename U>
	auto Matrix<T, Cols, Rows>::operator- (const Matrix<U, Cols, Rows>& other)
	{
		using Return = decltype(std::declval<T>() - std::declval<U>());

		Matrix<Return, Cols, Rows> newMatrix;

		std::transform(data.begin(), data.end(), other.data.begin(), newMatrix.data.begin(), [](const T& lhs, const U& rhs)
		{
			return lhs - rhs;
		});

		return newMatrix;
	}

	template<typename T, std::size_t Cols, std::size_t Rows>
	template<typename U, std::size_t OtherCols>
	auto Matrix<T, Cols, Rows>::operator* (const Matrix<U, OtherCols, Cols>& other)
	{
		using Return = decltype(std::declval<T>() * std::declval<U>());

		Matrix<Return, OtherCols, Rows> newMatrix;

		for(int i = 0; i < newMatrix.rows(); ++i)
		{
			for(int j = 0; j < newMatrix.cols(); ++j)
			{
				for(int k = 0; k < Cols; ++k)
				{
					newMatrix(j, i) += this->operator()(k, i) * other(j, k);
				}
			}
		}

		return newMatrix;
	}

	template<typename T, std::size_t Cols, std::size_t Rows, typename U>
	auto operator* (const Matrix<T, Cols, Rows>& matrix, const U& scalar)
	{
		using Return = decltype(std::declval<T>() * std::declval<U>());

		Matrix<Return, Cols, Rows> newMatrix;

		std::transform(matrix.data.begin(), matrix.data.end(), newMatrix.data.begin(), [](const T& in)
		{
			return in * scalar;
		});

		return newMatrix;
	}

	template<typename T, std::size_t Cols, std::size_t Rows, typename U>
	auto operator* (const U& scalar, const Matrix<T, Cols, Rows>& matrix)
	{
		return matrix * scalar;
	}
}

#ifdef MATRIX_PRINT

#include <string>
#include <iostream>

namespace swm
{
	template<typename T, std::size_t Cols, std::size_t Rows>
	void printMatrix(const Matrix<T, Cols, Rows>& mat, const std::string& rowPrefix = "")
	{
		for(int i = 0; i < Rows; ++i)
		{
			std::cout << rowPrefix;

			for(int j = 0; j < Cols; ++j)
			{
				std::cout << mat(j, i) << ' ';
			}

			std::cout << '\n';
		}

		std::cout << '\n';
	}
}

#endif

#endif
