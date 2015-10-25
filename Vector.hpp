#ifndef SWM_VECTOR_HPP
#define SWM_VECTOR_HPP

#include "Matrix.hpp"

namespace swm
{
	template<typename T>
	using Vector2 = Matrix<T, 1, 2>;

	template<typename T>
	using Vector3 = Matrix<T, 1, 3>;

	template<typename T>
	using Vector4 = Matrix<T, 1, 4>;
}

#endif
